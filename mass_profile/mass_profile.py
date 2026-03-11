import argparse
from mpi4py import MPI
from bigfile import BigFile, FileMPI
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from scipy.interpolate import interp1d

def print_message(message, rank_print=0):
    """Print message from rank 0 and flush."""
    if rank == rank_print:
        print(message, flush=True)

def calc_halo_mass_profile(pig, halo_idx, nbin, abs_redge_min, rel_redge_max, rho_crit, rho_m_cb, atime, ptypes = [0,1,4,5]):
    """
    Calculate mass profile for a given halo.

    Parameters:
    pig : BigFile
        The BigFile object containing simulation data.
    halo_idx : int
        The index of the halo to process.
    rel_rbin_redge : array_like
        The relative radial bin right edges (in units of R_200m).
    rho_crit : float
        The critical density at the simulation redshift (in 1e10 Msun/h / (kpc/h)^3).
    rho_m_cb : float
        The mean matter density (CDM + baryon) at the simulation redshift (in 1e10 Msun/h / (kpc/h)^3).
    ptypes : list of int
        The particle types to include in the mass profile calculation.
    Returns:
    M_200m : float
        The mass within R_200m (in 1e10 Msun/h).
    R_200m : float
        The radius R_200m (in kpc/h).
    M_200c : float
        The mass within R_200c (in 1e10 Msun/h).
    R_200c : float
        The radius R_200c (in kpc/h).
    M_500c : float
        The mass within R_500c (in 1e10 Msun/h).
    R_500c : float
        The radius R_500c (in kpc/h).
    mass_profile_dict : dict
        Dictionary of cumulative mass profiles for different particle types.
    """
    # R_200m > R_200c > R_500c
    # search between 5 knp/h to 10000 knp/h for R_200m
    # 1000 bins logarithmically spaced (later interpolate to get exact R_200m)
    
    # check if the FOF group has any subhalo
    if pig['FOFGroups/GroupNsubs'][halo_idx] == 0:
        print(f"Warning: Halo {halo_idx} has no subhalo, skipping mass profile calculation.", flush=True)
        mr_profile_dict_dummy = {}
        for ptype in ptypes:
            mr_profile_dict_dummy[ptype] = np.zeros(nbin)
        mr_profile_dict_dummy['r'] = np.zeros(nbin)
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, mr_profile_dict_dummy

    # locate particle start and end indices for this halo
    i_subhalo = pig['FOFGroups/GroupFirstSub'][halo_idx]
    pstart = pig['FOFGroups/OffsetByType'][halo_idx]
    pend = pstart + pig['FOFGroups/LengthByType'][halo_idx]
    Lbox = pig['Header'].attrs['BoxSize'][0]  # in kpc/h
    pos_halo = pig['SubGroups/SubhaloPos'][i_subhalo]  # in kpc/h

    # r bin right edges in kpc/h, sphere radius
    r_redge = np.logspace(np.log10(abs_redge_min), np.log10(10000), 1000)  # comoving radius in kpc/h

    # particle mass, central halo of this group
    # position and mass of particle
    mr = {} # M(<r)
    for ptype in ptypes:
        pos_part = pig[f"{ptype}/Position"][pstart[ptype]:pend[ptype]]
        # apply periodic BC
        del_pos = pos_part - pos_halo
        del_pos[del_pos < -Lbox / 2] += Lbox 
        del_pos[del_pos > Lbox / 2] -= Lbox

        mass_part = pig[f"{ptype}/Mass"][pstart[ptype]:pend[ptype]]

        # distribute particles into radial bins
        del_r = np.linalg.norm(del_pos, axis=-1)
        r_edge = np.concatenate(([0], r_redge))
        m_r, _ = np.histogram(del_r, bins=r_edge, weights=mass_part)
        mr[ptype] = np.cumsum(m_r)

    mr_total  = np.zeros_like(r_redge)
    for ptype in mr:
        mr_total += mr[ptype]

    vr = 4/3 * np.pi * (atime*r_redge)**3 # physical volume in (kpc/h)^3
    rhor = mr_total / vr  # in 1e10 Msun/h / (kpc/h)^3

    # ensure that R_200m is within the searched range: 200*rho_m_cb > rhor[-1]
    if not 200 * rho_m_cb > rhor[-1]:
        print(f"Error: Halo {halo_idx} R_200m is beyond 10000 kpc/h, consider increasing search range.", flush=True)
        raise ValueError("R_200m out of bounds.")

    # compute M200m, R200m, M200c, R200c, M500c, R500c
    # interpolate and find the radii: log-log interpolation
    # check if any rhor <= 0:
    rhor_positive = rhor[rhor > 0]
    if len(rhor_positive) < len(rhor):
        print(f"Warning: Halo {halo_idx} has {len(rhor) - len(rhor_positive)} bins with non-positive density, excluding them from interpolation.", flush=True)

    lg_rrho = interp1d(np.log(rhor), np.log(atime*r_redge), kind='linear') # physical radius instead of comoving
    # exit if out of bounds
    if (np.log(500 * rho_crit) > np.log(rhor).max()):
        # warning message and return zeroes
        print(f"Warning: Halo {halo_idx} R_500c is smaller than the minimum searched radius {atime*r_redge[0]} kpc/h, setting quantity to zero.", flush=True)
        mr_profile_dict_dummy = {}
        for ptype in ptypes:
            mr_profile_dict_dummy[ptype] = np.zeros(nbin)
        mr_profile_dict_dummy['r'] = np.zeros(nbin)
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, mr_profile_dict_dummy
    
    R_200m = np.exp(lg_rrho(np.log(200 * rho_m_cb)))  # in kpc/h
    R_200c = np.exp(lg_rrho(np.log(200 * rho_crit)))  # in kpc/h
    R_500c = np.exp(lg_rrho(np.log(500 * rho_crit)))  # in kpc/h

    # interpolate and find the masses: log-log interpolation
    lg_rmr = interp1d(np.log(atime*r_redge), np.log(mr_total), kind='linear')
    M_200m = np.exp(lg_rmr(np.log(R_200m)))  # in 1e10 Msun/h
    M_200c = np.exp(lg_rmr(np.log(R_200c)))  # in 1e10 Msun/h
    M_500c = np.exp(lg_rmr(np.log(R_500c)))  # in 1e10 Msun/h   

    # determine relative radial bins for mass profile
    rel_redge_min = abs_redge_min / (R_200m / atime)
    rel_rbin_redge = np.logspace(np.log10(rel_redge_min), np.log10(rel_redge_max), nbin)

    # cumulative mass profile for different particle types
    # define radial bins relative to R_200m, here comoving for convenience: Position is comoving
    r_redge_profile = (R_200m/atime) * rel_rbin_redge
    r_bin_profile = np.concatenate(([0], r_redge_profile))

    mr_profile_dict = {}
    for ptype in ptypes:
        pos_part = pig[f"{ptype}/Position"][pstart[ptype]:pend[ptype]]
        # apply periodic BC
        del_pos = pos_part - pos_halo
        del_pos[del_pos < -Lbox / 2] += Lbox 
        del_pos[del_pos > Lbox / 2] -= Lbox

        mass_part = pig[f"{ptype}/Mass"][pstart[ptype]:pend[ptype]]

        # distribute particles into radial bins
        del_r = np.linalg.norm(del_pos, axis=-1)
        m_r, _ = np.histogram(del_r, bins=r_bin_profile, weights=mass_part)
        mr_profile_dict[ptype] = np.cumsum(m_r)

    # r right edges for profile in physical kpc/h
    mr_profile_dict['r'] = r_redge_profile * atime  # in physical kpc/h
        
    return M_200m, R_200m, M_200c, R_200c, M_500c, R_500c, mr_profile_dict

# example usage
# mpirun -n 4 python mass_profile.py --pigdir /path/to/pigfile --savedir /path/to/save/dir --mindmpart 100 --maxrad 3.0 --minrad 0.01 --nrbin 100 --mnu 0.06

if __name__ == "__main__":

    #----------- MPI Init ----------------
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    parser = argparse.ArgumentParser(description='Calculate halo mass profile (cumulative).')
    parser.add_argument('--pigdir', required=True, type=str, help='the path to the pig file')
    parser.add_argument('--savedir', default=None, type=str, help='the path to  save the output')
    parser.add_argument('--mindmpart', default=100, type=int, help='minimum number of dm particles in halo to process')
    parser.add_argument('--maxradrel', default=3.0, type=float, help='maximum radius (in unit of R200m) to calculate mass profile')
    # parser.add_argument('--minrad', default=0.01, type=float, help='minimum radius (in unit of R200m) to calculate mass profile')
    parser.add_argument('--minradabs', default=5, type=float, help='minimum comoving radius (in kpc/h) to calculate mass profile')
    parser.add_argument('--nrbin', default=100, type=int, help='number of radial bins (logarithmic) to calculate mass profile')
    parser.add_argument('--mnu', default=0.06, type=float, help='total neutrino mass in eV (for cosmology)')
    parser.add_argument('--verbose', action='store_true', help='print more messages')

    args = parser.parse_args()

    # set savedir: the same as pigdir if not specified
    if args.savedir is None:
        args.savedir = args.pigdir

    comm.barrier()
    
    # open: read and write
    pig = BigFile(args.pigdir)
    dest_w  = FileMPI(comm, args.savedir, create=True)

    # how many halos to process
    # find the first FoF group with less than mindmpart DM particles
    len_dm = pig["FOFGroups/LengthByType"][:][:,1]
    # where len_dm starts to be < mindmpart
    n_halo_target = (len_dm < args.mindmpart).nonzero()[0]
    if n_halo_target.size > 0:
        n_halo_target = n_halo_target[0]
        # print and flush from rank 0
        print_message(f"The first {n_halo_target} halos contain >= {args.mindmpart} DM particles.")
    else:
        n_halo_target = len(len_dm)
        print_message(f"All {n_halo_target} halos have >= {args.mindmpart} DM particles.")

    comm.barrier()

    # print some info from rank 0
    print_message(f"args {args}")
    print_message(f'saving at: {args.savedir}')
    print_message(f"mpi rank size {size}")

    # prepare output blocks
    # get non-trivial particle types
    part_types = pig["Header"].attrs['NumPartInGroupTotal'].nonzero()[0]
    blockname_list = ['M200m', 'R200m', 'M200c', 'R200c', 'M500c', 'R500c', 'RCumProfile'] + [('MassCumProfileType' + str(ptype)) for ptype in part_types]
    # under 'FOFGroups'
    blockname_list = ['FOFGroups/' + name for name in blockname_list]

    blocktype_list = ['<f4', '<f4', '<f4', '<f4', '<f4', '<f4', ('<f4', (args.nrbin,))] + [('<f4', (args.nrbin,)) for _ in part_types]
    n_file = pig["FOFGroups/LengthByType"].Nfile

    # create or open blocks
    blocklist = []
    for blockidx, name in enumerate(blockname_list):
        # create
        dest_w.create(blockname_list[blockidx], blocktype_list[blockidx], n_halo_target, n_file)
        blocklist.append(dest_w[blockname_list[blockidx]])

    comm.barrier()

    # calculate critical density and mean matter density at this redshift
    hh = pig["Header"].attrs['HubbleParam'][0]
    Om0 = pig["Header"].attrs['Omega0'][0]
    Ob0 = pig["Header"].attrs['OmegaBaryon'][0]
    Tcmb0 = pig["Header"].attrs['CMBTemperature'][0]
    atime = pig["Header"].attrs['Time'][0]
    z = 1.0/atime - 1.0
    cosmo = FlatLambdaCDM(
    H0 = 100 * hh * u.km/u.s/u.Mpc,
    Om0 = Om0,
    Ob0 = Ob0,
    m_nu = [args.mnu, 0, 0] * u.eV,  # total neutrino mass
    Tcmb0 = Tcmb0 * u.K
    )
    rho_crit = cosmo.critical_density(z).to(u.Msun / u.kpc**3).value * 1e-10 / hh**2  # in MP-Gadget units: (1e10 Msun/h) / (kpc/h)^3
    rho_crit0 = cosmo.critical_density0.to(u.Msun / u.kpc**3).value * 1e-10 / hh**2
    Ocb = cosmo.Om(z) - cosmo.Onu(z)
    rho_m_cb = cosmo.critical_density(z).to(u.Msun / u.kpc**3).value * Ocb * 1e-10 / hh**2  # mean matter density (CDM + baryon), in MP-Gadget units; do not include neutrinos in this context

    # split halos among ranks
    i_halo = 0
    local_processed = 0
    # process halos
    while(i_halo < n_halo_target):
        # only process halos assigned to this rank
        if (i_halo % size) == rank:
            # placeholder for mass profile calculation
            M_200m, R_200m, M_200c, R_200c, M_500c, R_500c, mass_profile_dict = \
                calc_halo_mass_profile(pig, i_halo, args.nrbin, args.minradabs, args.maxradrel, rho_crit, rho_m_cb, atime, ptypes=part_types)
            # write to blocks
            blocklist[0].write(i_halo, np.array([M_200m]))
            blocklist[1].write(i_halo, np.array([R_200m]))
            blocklist[2].write(i_halo, np.array([M_200c]))
            blocklist[3].write(i_halo, np.array([R_200c]))
            blocklist[4].write(i_halo, np.array([M_500c]))
            blocklist[5].write(i_halo, np.array([R_500c]))
            blocklist[6].write(i_halo, mass_profile_dict['r'])
            for ptype_idx, ptype in enumerate(part_types):
                mass_profile = mass_profile_dict.get(ptype)
                blocklist[7 + ptype_idx].write(i_halo, mass_profile)
            local_processed += 1
            if args.verbose and (local_processed % 10 == 0) and rank == 1:
                print(f"Rank {rank}: processed {local_processed} halos.", flush=True)

        i_halo += 1

    # print rank finished message
    print(f"Rank {rank} finished.", flush=True)
    comm.barrier()
    
    # add attribute for radial bins
    # if rank == 0:
    #     for ptype_idx, ptype in enumerate(part_types):
    #         blocklist[6 + ptype_idx].attrs['RelRBinREdgeR200m'] = rel_rbin_redge

    comm.barrier()
    print_message("Finished processing halos.")