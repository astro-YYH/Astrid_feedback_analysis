import numpy as np
import os

def convert_units(data_dict):
    f_mass = 1e10  # 1e10 Msun/h to Msun/h
    f_energy = 1e10 * (1e5)**2 * 1.989e33  # 1e10 Msun/h *(km/s)^2 to erg/h
    c_kmps = 3e5  # speed of light in km/s
    f_energy_in_mass = 1e10 / (c_kmps)**2  # 1e10 Msun/h *(km/s)^2 to Msun/h * c^2  (E=mc^2)
    data_dict['FirstReleaseMass'] *= f_mass  # to Msun/h
    data_dict['BHMass'] *= f_mass
    data_dict['ReleasedKE'] *= f_energy_in_mass  # to Msun/h * c^2

    return data_dict

def plot_with_percentiles(ax, x, y, xbins=30, color='blue', xlower=None, xupper=None, label=None):
    # define bins in logarithmic space
    if xlower is None:
        xlower = np.min(x)
    if xupper is None:
        xupper = np.max(x)
    x_edges = np.logspace(np.log10(xlower), np.log10(xupper), xbins+1)
    x_centers = np.sqrt(x_edges[:-1] * x_edges[1:])

    y_median = []
    y_p16 = []  # 16th percentile in log space
    y_p84 = []

    for i in range(xbins):
        mask = (x >= x_edges[i]) & (x < x_edges[i+1])
        if np.sum(mask) > 0:
            y_bin = y[mask]
            y_median.append(np.median(y_bin))
            y_p16.append(np.percentile(y_bin, 16))
            y_p84.append(np.percentile(y_bin, 84))
        else:
            y_median.append(np.nan)
            y_p16.append(np.nan)
            y_p84.append(np.nan)

    y_median = np.array(y_median)
    y_p16 = np.array(y_p16)
    y_p84 = np.array(y_p84)
    ax.plot(x_centers, y_median, linestyle='-', color=color, label=label)
    # print(f"Median y values: {y_median}")
    ax.fill_between(x_centers, y_p16, y_p84, color=color, alpha=0.3)
    
def get_kBHinfo_all_feedback(kBHinfo):
    """
    Get the total feedback energy released by the primary BH and all its swallowed BHs by the target end time (the max).
    """
    # sort by EndTime to ensure primary BHs come after their swallowed BHs
    sort_indices = np.argsort(kBHinfo['EndTime'])
    kBHinfo_sorted = {key: val[sort_indices] for key, val in kBHinfo.items()}

    # for debugging
    energy_missing = 0.0
    energy_total = np.sum(kBHinfo_sorted['ReleasedKE'])

    # add the feedback energy from swallowed BHs to their primary BHs, until the target end time
    for i in range(len(kBHinfo_sorted['BHID'])):
        # stop when reaching the target end time
        if kBHinfo_sorted['EndTime'][i] == kBHinfo_sorted['EndTime'][-1]:
            break
        # skip if the BH has not released any feedback energy
        if kBHinfo_sorted['ReleasedKE'][i] <= 0:
            continue

        swallower_id = kBHinfo_sorted['SwallowID'][i]

        # add the feedback energy to the swallower BH if it exists
        idx_swallower = np.where(kBHinfo_sorted['BHID'] == swallower_id)[0]

        # ensure len(idx_swallower) <= 1
        assert len(idx_swallower) <= 1, f"Error: Multiple swallower BHs found for BHID {kBHinfo_sorted['BHID'][i]}"

        if len(idx_swallower) > 0:
            idx_swallower = idx_swallower[0]
            # print message if the swallower BH has not released any feedback energy yet
            if kBHinfo_sorted['ReleasedKE'][idx_swallower] <= 0:
                print(f"Info: Swallower BHID {swallower_id} has not released feedback energy yet") # this is possible but rare
            kBHinfo_sorted['ReleasedKE'][idx_swallower] += kBHinfo_sorted['ReleasedKE'][i]
            # update first release time and mass if needed
            if kBHinfo_sorted['FirstReleaseTime'][idx_swallower] > kBHinfo_sorted['FirstReleaseTime'][i] or kBHinfo_sorted['FirstReleaseTime'][idx_swallower] < 0:
                # print the release time
                print(f"Updating FirstReleaseTime of BHID {swallower_id} from {kBHinfo_sorted['FirstReleaseTime'][idx_swallower]} to {kBHinfo_sorted['FirstReleaseTime'][i]} (from swallowed BHID {kBHinfo_sorted['BHID'][i]})")
                kBHinfo_sorted['FirstReleaseTime'][idx_swallower] = kBHinfo_sorted['FirstReleaseTime'][i]
                kBHinfo_sorted['FirstReleaseMass'][idx_swallower] = kBHinfo_sorted['FirstReleaseMass'][i]
        else:
            print(f"Warning: Swallower BHID {swallower_id} not found for BHID {kBHinfo_sorted['BHID'][i]}") # we don't want this to happen, but it can happen if inactive (thus not recorded) BHs are swallowed
            # print end time
            energy_missing += kBHinfo_sorted['ReleasedKE'][i]
            print(f"EndTime of this BH: {kBHinfo_sorted['EndTime'][i]}")

    # print energy summary
    print(f"Total feedback energy: {energy_total:.3e}")
    print(f"Missing feedback energy (from swallowed BHs whose swallower BHs are not found): {energy_missing:.3e}")

    # we only need the primary BHs at the target end time
    mask_primary = kBHinfo_sorted['EndTime'] == kBHinfo_sorted['EndTime'][-1]
    kBHinfo_allfdbk = {key: val[mask_primary] for key, val in kBHinfo_sorted.items()}
    return kBHinfo_allfdbk

def get_kBH_end(filename, all_feedback=True, released_only=True):
    kBHinfo = dict(np.load(filename))
    if all_feedback:
        kBHinfo = get_kBHinfo_all_feedback(kBHinfo)
    # retain only those at the target end time (the max) and not swallowed
    mask = (kBHinfo['EndTime'] == kBHinfo['EndTime'].max()) & (kBHinfo['SwallowID'] > kBHinfo['BHID'].max())
    
    kBHinfo = {key: val[mask] for key, val in kBHinfo.items()}
    if released_only:
        # retain only those that released feedback energy
        mask_feedback = kBHinfo['ReleasedKE'] > 0
        kBHinfo = {key: val[mask_feedback] for key, val in kBHinfo.items()}
    # convert units
    kBHinfo = convert_units(kBHinfo)
    return kBHinfo

# function to plot energy density distribution
def plot_energy_density_distribution(ax, kBHinfo, nbins=25, xlower=1e7, xupper=1e10, volume=50**3, color='blue', label=None, linewidth=2, masstype='BHMass', ls='-'):
    # bins in log space
    mbins = np.logspace(np.log10(xlower), np.log10(xupper), nbins+1)
    de, _ = np.histogram(kBHinfo[masstype], weights=kBHinfo['ReleasedKE'], bins=mbins)
    # energy density distribution
    du_dlogm = de / volume / np.diff(np.log10(mbins))  # [c^2 Msun/h / (Mpc/h)^3]

    mbin_centers = np.sqrt(mbins[:-1] * mbins[1:])
    # cut the leading and trailing zeros for better visualization
    nonzero_mask = du_dlogm > 0
    mbin_centers = mbin_centers[nonzero_mask]
    du_dlogm = du_dlogm[nonzero_mask]
    ax.plot(mbin_centers, du_dlogm, linestyle=ls, color=color, label=label, linewidth=linewidth)
    # print total energy density (sum over all bins)
    total_energy_density = np.sum(de) / volume  # [c^2 Msun/h / (Mpc/h)^3]
    print(f"Total energy density (from {xlower:.1e} to {xupper:.1e} Msun/h): {total_energy_density:.3e} c^2 Msun/h / (Mpc/h)^3")

    # return the total energy density for reference
    return total_energy_density

def get_bhhalo(filename):
    bhhalos = dict(np.load(filename))
    # convert units if needed (Mass in Msun/h)
    f_mass = 1e10  # 1e10 Msun/h to Msun/h
    bhhalos['GroupMass'] *= f_mass
    bhhalos['M200c'] *= f_mass
    bhhalos['M200m'] *= f_mass
    bhhalos['M500c'] *= f_mass
    return bhhalos

def get_kBH_end_zs(sim_name, zs, all_feedback=True, kbh_base='kbhinfo', bhhalo_base='bhhaloinfo', halo_info=False, released_only=True, swallow_count=False, swallow_base='bhswallowcount'):
    kbh_data = {}
    for z in zs:
        if z == 0.0:
            filename = f'{kbh_base}_{sim_name}.npz'
        else:
            filename = f'{kbh_base}_{sim_name}_z{z:.1f}.npz'
        # check if file exists
        if not os.path.isfile(filename):
            print(f'File {filename} does not exist. Skipping z={z}.')
            continue
        kbh_data[z] = get_kBH_end(filename, all_feedback=all_feedback, released_only=released_only)
        # sort by BHID
        sort_indices = np.argsort(kbh_data[z]['BHID'])
        kbh_data[z] = {key: val[sort_indices] for key, val in kbh_data[z].items()}

        if halo_info:
            if z == 0.0:
                halo_filename = f'{bhhalo_base}_{sim_name}.npz'
            else:
                halo_filename = f'{bhhalo_base}_{sim_name}_z{z:.1f}.npz'
            if not os.path.isfile(halo_filename):
                print(f'File {halo_filename} does not exist. Skipping halo info for z={z}.')
                continue
            bhhalos = get_bhhalo(halo_filename)
            # if there are BHs in kbh_data[z] that are not in bhhalos, print a warning, and add a dummy entry with zero mass in bhhalos
            mask = np.isin(kbh_data[z]['BHID'], bhhalos['BHID'])
            if len(kbh_data[z]['BHID']) != np.sum(mask):
                print("Warning: {} BH(s) at z={} not found in bhhalo data. Adding dummy entries with zero mass.".format(len(kbh_data[z]['BHID']) - np.sum(mask), z))
                missing_bhids = kbh_data[z]['BHID'][~mask]
                for bhid in missing_bhids:
                    for key in bhhalos:
                        if key == 'BHID':
                            bhhalos[key] = np.append(bhhalos[key], bhid)
                        elif key == 'Time':
                            continue  # skip Time
                        else:
                            # zero in appropriate dtype
                            bhhalos[key] = np.append(bhhalos[key], 0)  # GroupID==0 means not in any halo

            # extract only those halos corresponding to the kBHs (if release_only is True, these are only those that released feedback)
            mask = np.isin(bhhalos['BHID'], kbh_data[z]['BHID'])
            # print a message about the number of matched BHs and the total number of BHs in bhhalos
            print(f"At z={z}, matched {np.sum(mask)} BHs out of {len(bhhalos['BHID'])} total BHs in bhhalo data.")
            bhhalos = {key: val[mask] for key, val in bhhalos.items() if key != 'Time'}  # skip Time

            # reorder bhhalos and kbh_data[z] to have the same order of BHID
            sort_indices = np.argsort(bhhalos['BHID'])
            bhhalos = {key: val[sort_indices] for key, val in bhhalos.items() if key != 'Time'}  # skip Time

            # ensure that the BHIDs match
            assert np.array_equal(kbh_data[z]['BHID'], bhhalos['BHID']), "Error: BHIDs do not match between kBH data and bhhalo data."

            for key in bhhalos:
                if key == 'Time' or key == 'BHID':
                    continue
                kbh_data[z][key] = bhhalos[key]

        if swallow_count:
            swallow_filename = f'{swallow_base}_{sim_name}_z{z:.1f}.npz'
            if not os.path.isfile(swallow_filename):
                print(f'File {swallow_filename} does not exist. Skipping swallow count for z={z}.')
                continue
            print(f'Loading swallow count from {swallow_filename} for z={z}.')
            swallow_data = dict(np.load(swallow_filename))
            # ensure time matches
            # print('the input swallow file is', swallow_filename) # for debugging
            assert np.isclose(1/swallow_data['Time'] - 1, z), f"Error: Time in swallow data does not match z={z}."

            # confirm all BHIDs in kbh_data[z] are in swallow_data (this should always be true)
            if not np.all(np.isin(kbh_data[z]['BHID'], swallow_data['BHID'])):
                missing_bhids = kbh_data[z]['BHID'][~np.isin(kbh_data[z]['BHID'], swallow_data['BHID'])]
                raise ValueError(f"Error: The following BHIDs in kBH data at z={z} are not found in swallow data: {missing_bhids}")

            # extract only those BHs that are in kbh_data[z]
            mask = np.isin(swallow_data['BHID'], kbh_data[z]['BHID'])
            swallow_data = {key: val[mask] for key, val in swallow_data.items() if key != 'Time'}  # skip Time
            # reorder swallow_data and kbh_data[z] to have the same order of BHID
            sort_indices = np.argsort(swallow_data['BHID'])
            swallow_data = {key: val[sort_indices] for key, val in swallow_data.items() if key != 'Time'}  # skip Time

            # ensure that the BHIDs match
            assert np.array_equal(kbh_data[z]['BHID'], swallow_data['BHID']), "Error: BHIDs do not match between kBH data and swallow data."

            kbh_data[z]['SwallowCount'] = swallow_data['SwallowCount']

    return kbh_data