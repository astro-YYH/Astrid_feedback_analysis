import glob
import numpy as np
import os
import argparse
from bigfile import BigFile

# get black hole details from one 'BlackholeDetails'
def get_blackholedetails(file, fields_target=['BHID', 'BHMass', 'KineticFdbkEnergy', 'KEflag', 'z', 'SwallowID'], oldformat=False): # get only BHID, KineticFdbkEnergy, KEflag, z (scale factor)
    fields = ('size','BHID','BHMass','Mdot','Density','timebin','Encounter','MinPos',\
        'MinPot','Entropy','GasVel','acMom','acMass','acBHMass',\
        'Fdbk','SPHID','SwallowID','CountProgs','Swallowed',\
        'BHpos','srDensity','srParticles','srVel','srDisp',\
        'DFAccel','DragAccel','GravAccel','BHvel','Mtrack','Mdyn',\
        'KineticFdbkEnergy','NumDM','V1sumDM','V2sumDM','MgasEnc','KEflag','z','size2')
    dtype = ('i','Q','d','d','d','i','i','3d',\
        'd','d','3d','3d','d','d',\
        'd','Q','Q','i','i',\
        '3d','d','d','3d','d',\
        '3d','3d','3d','3d','d','d',\
        'd','d','3d','d','d','i','d','i') if not oldformat else ('i','Q','d','d','d','i','i','3d',\
        'd','d','3d','3d','d','d',\
        'd','Q','Q','i','i',\
        '3d','d','d','3d','d',\
        '3d','3d','3d','3d','d','d',\
        'd','d','3d','d','d','d','d','i')
    dt = {'names':fields, 'formats':dtype}

    dt_target = {'names':fields_target, 'formats':[dt['formats'][dt['names'].index(f)] for f in fields_target]}

    datadir = file
    bhfile_zoo = glob.glob(os.path.join(datadir,"*"))


    results = []
    for filename in bhfile_zoo:
        data = np.fromfile(filename, dtype=dt, count=-1)
        # results.append(data[fields_target])

        tmp = np.empty(data.shape[0], dtype=dt_target)
        for name in dt_target['names']:
            tmp[name] = data[name]

        results.append(tmp)
        del data
    results = np.concatenate(results)

    bhinfo = results
        
    return bhinfo

# get black hole details from a list of 'BlackholeDetails' files
def get_bhinfo(file_list, fields_target=['BHID', 'BHMass', 'KineticFdbkEnergy', 'KEflag', 'z', 'SwallowID'], oldformat=False):
    all_bhinfo = []
    for file in file_list:
        print(f"Reading BH details from {file}...")
        bhinfo = get_blackholedetails(file, fields_target, oldformat=oldformat)
        all_bhinfo.append(bhinfo)
    all_bhinfo = np.concatenate(all_bhinfo)
    return all_bhinfo

# get activation time, the BH mass at activation and kinetic feedback energy released before a given scale factor for a given BHID
def get_bh_kinetic_info(bhinfo, bhid, scale_factor=1.):
    # filter by BHID and scale factor
    mask = (bhinfo['BHID'] == bhid) & (bhinfo['z'] <= scale_factor) & (bhinfo['KEflag'] != 0)

    # warn if no entries found
    if not np.any(mask):
        print(f"Warning: No entries found for BHID {bhid} at or before scale factor {scale_factor}.")
        return None, None, None
    # filtered_bhinfo = bhinfo[mask]
    # we need a copy for resorting later
    bhinfo_target = np.copy(bhinfo[mask])
    # sort by scale factor ascending
    sorted_indices = np.argsort(bhinfo_target['z'])
    bhinfo_target = bhinfo_target[sorted_indices]

    atime_activation = bhinfo_target['z'][0]
    mass_at_activation = bhinfo_target['BHMass'][0]
    # kinetic feedback energy released
    # find the indices where KEflag == 2 (release event)
    idx_release = np.where(bhinfo_target['KEflag'] == 2)[0]
    if len(idx_release) == 0:
        kenergy = 0.
    else:
        kenergy = np.sum(bhinfo_target['KineticFdbkEnergy'][idx_release-1])  # sum the energy released at each event (the energy is stored in the previous entry)

    return atime_activation, mass_at_activation, kenergy

def get_swallowid(swallowid_orig, part_snap, bhid, bhid_max, end_time, cut_time):
    swallowid = swallowid_orig
    # find SwallowID from PART if SwallowID at EndTime > max(BHID) (this happens when a BH is swallowed when it is inactive and thus not recorded) and the end time is prior to the cut-off time
    if swallowid_orig > bhid_max and end_time < cut_time:
        print(f"Warning: SwallowID {swallowid_orig} of BHID {bhid} at EndTime {end_time:.6f} exceeds max BHID {bhid_max}. Looking up in PART...")
        i_bh = np.where(part_snap['5/ID'][:] == bhid)[0]
        if len(i_bh) == 0:
            # print(f"Error: BHID {unique_bhids[i]} not found in PART snapshot!")
            raise ValueError(f"BHID {bhid} not found in PART snapshot!")
        else:
            swallowid = part_snap['5/BlackholeSwallowID'][i_bh[0]]
            print(f"Found SwallowID {swallowid} for BHID {bhid} from PART.")
    return swallowid

def get_swallower_series(swallowid, part_snap, cut_time):
    current_id = swallowid
    swallower_ids = []
    while True:
        # find SwallowID from PART
        i_bh = np.where(part_snap['5/ID'][:] == current_id)[0]
        if len(i_bh) == 0:
            raise ValueError(f"BHID {current_id} not found in PART snapshot!")
        # check time
        end_time = part_snap['5/BlackholeSwallowTime'][i_bh[0]]
        if end_time >= cut_time:
            break # stop if the end time exceeds cut-off

        # append current_id to swallower_ids
        swallower_ids.append(current_id)

        # get next swallower id
        current_id = part_snap['5/BlackholeSwallowID'][i_bh[0]]
        
    return swallower_ids
    
def get_kinetic_bhinfo(file_list, part_snap, scale_factor=1., oldformat=False):
    data = get_bhinfo(file_list, fields_target=['BHID', 'BHMass', 'KineticFdbkEnergy', 'KEflag', 'z', 'SwallowID'], oldformat=oldformat)

    # find the 'z' value closest to scale_factor as the cut-off
    unique_z = np.unique(data['z'])
    z_cut = unique_z[np.argmin(np.abs(unique_z - scale_factor))]

    print(f"Using scale factor={z_cut:.6f} (redshift={1/z_cut-1:.6f}) as cut-off scale factor.")

    # Filter once
    print("Filtering BH info...")
    # mask = (bhinfo['KEflag'] != 0) & (bhinfo['z'] <= scale_factor)
    mask = (data['z'] <= z_cut) # should not filter out KEflag==0 here, as some BHs may turn on kinetic feedback and then off again before scale_factor, while they still contribute to the total feedback energy
    data = data[mask]

    # Sort by BHID, then z
    order = np.lexsort((data['z'], data['BHID']))
    data = data[order]

    bhid = data['BHID']
    z = data['z']
    ke = data['KineticFdbkEnergy']
    keflag = data['KEflag']
    mass = data['BHMass']

    # Find BHID boundaries
    unique_bhids, start_idx = np.unique(bhid, return_index=True)

    out = {
        'BHID': [],
        'BHMass': [],
        # 'ActivationTime': [], # the original activation time is not that useful (early BHs are always in the low-accretion mode initially but no feedback energy is released), so instead, we focus on the time a BH releases feedback energy for the first time
        # 'ActivationMass': [],
        'FirstReleaseTime': [],
        'FirstReleaseMass': [],
        'ReleasedKE': [],
        'EndTime': [],  # last recorded time (when it is swallowed or last snapshot)
        'SwallowID': []
    }

    print("Processing BH entries...")
    for i, start in enumerate(start_idx):
        end = start_idx[i + 1] if i + 1 < len(start_idx) else len(bhid)

        z_bh = z[start:end]
        # skip if the BH does not reach the cut-off time
        # if z_bh[-1] < z_cut:
        #     continue
        ke_bh = ke[start:end]
        keflag_bh = keflag[start:end]

        # skip if the BH never releases any kinetic feedback energy and does not reach the cut-off time
        if z_bh[-1] < z_cut and not np.any(keflag_bh == 2):
            continue

        mass_bh = mass[start:end]
        swallowid_bh = data['SwallowID'][start:end]

        swallowid = get_swallowid(swallowid_bh[-1], part_snap, unique_bhids[i], unique_bhids[-1], z_bh[-1], z_cut)

        # activation
        out['BHID'].append(unique_bhids[i])
        # out['ActivationTime'].append(z_bh[0])
        # out['ActivationMass'].append(mass_bh[0])
        out['BHMass'].append(mass_bh[-1])
        out['EndTime'].append(z_bh[-1])
        out['SwallowID'].append(swallowid)
        # first release time and mass: the first time KEflag == 2
        idx_first_release = np.where(keflag_bh == 2)[0]
        if len(idx_first_release):
            first_idx = idx_first_release[0]
            out['FirstReleaseTime'].append(z_bh[first_idx])
            out['FirstReleaseMass'].append(mass_bh[first_idx])
        else:
            out['FirstReleaseTime'].append(-1.0)  # indicate no release
            out['FirstReleaseMass'].append(-1.0)

        # released KE (flag==2 → energy in previous entry)
        idx = np.where(keflag_bh == 2)[0]
        if len(idx):
            out['ReleasedKE'].append(np.sum(ke_bh[idx - 1]))
        else:
            out['ReleasedKE'].append(0.0)

    # loop over BH entries in out to check whether the SwallowID exists in out, and if not, find it in data and add it 
    # list of swallower ids (that need to be added later)
    swallower_ids = []
    for i in range(len(out['BHID'])):
        # skip if EndTime == cut_time
        if out['EndTime'][i] == z_cut:
            continue
            
        swallowid = out['SwallowID'][i]
        if swallowid not in out['BHID']:
            # find all swallower ids for this swallowid (the BH swallowed the swallower, and so on)
            swallower_id_series = get_swallower_series(swallowid, part_snap, z_cut)
            swallower_ids.extend(swallower_id_series)

    # remove swallower ids that are already in out['BHID'] or duplicated
    swallower_ids = list(set([sid for sid in swallower_ids if sid not in out['BHID']]))

    print(f"Adding {len(swallower_ids)} swallower BH entries (not found in initial round of extraction)...")

    # find these swallower ids in data (unique BHIDs) and add them to out
    idx_swallowers = np.where(np.isin(unique_bhids, swallower_ids))[0]
    assert len(idx_swallowers) == len(swallower_ids), "Mismatch in number of swallower IDs found (there could be swallowers that cannot be found in BlackholeDetails)."

    for i in idx_swallowers:
        start = start_idx[i]
        end = start_idx[i + 1] if i + 1 < len(start_idx) else len(bhid)

        z_bh = z[start:end]
        ke_bh = ke[start:end]
        keflag_bh = keflag[start:end]
        mass_bh = mass[start:end]
        swallowid_bh = data['SwallowID'][start:end]
        swallowid = get_swallowid(swallowid_bh[-1], part_snap, unique_bhids[i], unique_bhids[-1], z_bh[-1], z_cut)  # can be optimized, since we have already found the swallower series above

        # activation
        out['BHID'].append(unique_bhids[i])
        # out['ActivationTime'].append(z_bh[0])
        # out['ActivationMass'].append(mass_bh[0])
        out['BHMass'].append(mass_bh[-1])
        out['EndTime'].append(z_bh[-1])
        out['SwallowID'].append(swallowid)
        # first release time and mass: the first time KEflag == 2
        idx_first_release = np.where(keflag_bh == 2)[0]
        if len(idx_first_release):
            first_idx = idx_first_release[0]
            out['FirstReleaseTime'].append(z_bh[first_idx])
            out['FirstReleaseMass'].append(mass_bh[first_idx])
        else:
            out['FirstReleaseTime'].append(-1.0)  # indicate no release
            out['FirstReleaseMass'].append(-1.0)
        # released KE (flag==2 → energy in previous entry)
        idx = np.where(keflag_bh == 2)[0]
        if len(idx):
            out['ReleasedKE'].append(np.sum(ke_bh[idx - 1]))
        else:
            out['ReleasedKE'].append(0.0)

        # ReleasedKE should always be zero for swallowers here (otherwise they would have been included in the first loop)
        assert out['ReleasedKE'][-1] == 0.0, f"ReleasedKE for swallower BH {unique_bhids[i]} should be zero."

    # convert lists to arrays
    for k in out:
        out[k] = np.asarray(out[k])

    return out

# example usage
# python extract_kbhinfo.py --inputdir /path/to/BlackholeDetails --basename BlackholeDetails --output kbhinfo.npz 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract KBH info from BlackholeDetails.")
    parser.add_argument("--inputdir", type=str, help="Input directory containing BlackholeDetails.")
    parser.add_argument("--part", type=str, help="PART snapshot file (ensure the time is not before the target cut-off).")
    parser.add_argument("--basename", default="BlackholeDetails", type=str, help="Base name of the BlackholeDetails files.")
    parser.add_argument("--output", type=str, help="Output file to save extracted KBH info.")
    parser.add_argument("--time", type=np.float64, default=1.0, help="Scale factor up to which to extract KBH info.")
    parser.add_argument("--oldformat", action='store_true', help="Use old format extraction method, i.e., KEflag is a double instead of int.")
    args = parser.parse_args()

    # print all arguments
    print("Arguments:")
    for arg in vars(args):
        print(f"  {arg}: {getattr(args, arg)}")

    bd_list = glob.glob(os.path.join(args.inputdir, args.basename + "*"))
    part_snap = BigFile(args.part)

    kBHinfo = get_kinetic_bhinfo(bd_list, part_snap, scale_factor=args.time, oldformat=args.oldformat)

    # save to output file as npz
    print(f"Saving KBH info to {args.output}...")
    np.savez(args.output, **kBHinfo)
    print(f"Saved KBH info to {args.output}.")