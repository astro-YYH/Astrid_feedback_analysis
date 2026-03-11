from scipy.interpolate import interp1d
import numpy as np


def loglog_interp(r_src, y_src, r_tgt):
    # exclude the first point where r=0 and then attach it back later
    if r_src[0] == 0:
        r_src_nonzero = r_src[1:]
        y_src_nonzero = y_src[1:]
        attach_zero = True
    else:
        r_src_nonzero = r_src
        y_src_nonzero = y_src
        attach_zero = False

    f = interp1d(np.log(r_src_nonzero), np.log(y_src_nonzero), kind='linear', fill_value='extrapolate')

    if r_tgt[0] == 0 and attach_zero:
        y_tgt_nonzero = np.exp(f(np.log(r_tgt[1:])))
        y_tgt = np.concatenate(([y_src[0]], y_tgt_nonzero))
    else:
        y_tgt = np.exp(f(np.log(r_tgt)))

    return y_tgt

def get_rho_profile_type(pig, ptype):
    r_redge = r_redge = pig['FOFGroups/RCumProfile'][:]

    hubble = pig['Header'].attrs['HubbleParam'][0]

    mr = pig[f'FOFGroups/MassCumProfileType{ptype}'][:]

    r_bin = np.concatenate((np.zeros((r_redge.shape[0], 1)), r_redge), axis=1)

    vol_shell = (4.0/3.0) * np.pi * (r_bin[:, 1:]**3 - r_bin[:, :-1]**3)

    mass_shell = np.zeros_like(mr)
    mass_shell[:, 0] = mr[:, 0]
    mass_shell[:, 1:] = mr[:, 1:] - mr[:, :-1]

    density_prof = mass_shell / vol_shell 

    r_bin_mid = np.sqrt(r_bin[:, 1:] * r_bin[:, :-1]) 

    density_prof *= 1e10 * hubble**2
    r_bin_mid /= hubble
    # print units
    # print(f'Units of rho: Msun/(kpc)^3')
    # print(f'Units of r: kpc')
    # # print Hubble
    # print(f'Hubble parameter: {hubble}')

    return r_bin_mid, density_prof

def get_rho_profile(pig):  # total matter density profile

    # extract available particle types
    # find FOFGroups/MassCumProfileType* in pig.keys()
    ptypes = []
    for key in pig.keys():
        if key.startswith('FOFGroups/MassCumProfileType'):
            ptype = int(key.split('Type')[-1])
            ptypes.append(ptype)

    # get density profile for each particle type
    density_profiles = []
    for ptype in ptypes:
        r_bin_mid, density_prof = get_rho_profile_type(pig, ptype)
        density_profiles.append(density_prof)

    # sum to get total density profile
    density_prof = np.sum(density_profiles, axis=0)

    return r_bin_mid, density_prof

# stacked density profile, halo mass from m_lower to m_upper
# def get_stacked_rho_profile(pig, m_lower, m_upper, mass_def='M500c'):
    # mass_def can be 'M200m', 'M200c' or 'M500c'
    # m_lower and m_upper in units of Msun/h
    halo_masses = pig[f'FOFGroups/{mass_def}'][:]
    # index mask for halos within the mass range
    mask = (halo_masses >= m_lower/1e10) & (halo_masses < m_upper/1e10)
    # warn if no halos found
    if np.sum(mask) == 0:
        print(f'Warning: No halos found in mass range {m_lower:.2e}--{m_upper:.2e} Msun/h')
        return None, None
    # get density profiles for selected halos
    r_bin_mid, density_prof = get_rho_profile(pig)
    selected_density_prof = density_prof[mask, :]
    r_bin_mid = r_bin_mid[mask, :]
    # interpolate to common radial bins: use the smallest range of r_bin_mid across selected halos
    i_common = np.argmin(r_bin_mid[:, -1])  # find where r_bin_mid[:, -1] is minimum
    r_common = r_bin_mid[i_common, :]
    # interpolate in log-log space
    density_prof_interp = np.zeros_like(selected_density_prof)
    for i in range(selected_density_prof.shape[0]):
        density_prof_interp[i, :] = loglog_interp(r_bin_mid[i, :], selected_density_prof[i, :], r_common)
    # average over halos
    stacked_density_prof = np.mean(density_prof_interp, axis=0)

    return r_common, stacked_density_prof

# stacked density profile, halo mass from m_lower to m_upper
def get_stacked_rho_profile(pig, m_lower=1e13, m_upper=1e14, mass_def='M500c', groupids=None):
    """
    Return stacked density profile and 16--84 percentile range.

    Parameters
    ----------
    pig : BigFile
    m_lower, m_upper : float
        Halo mass range in Msun/h
    mass_def : str
        'M200m', 'M200c', or 'M500c'

    Returns
    -------
    r_common : array
        Common radial bins
    rho_med : array
        Median stacked density profile
    rho_p16 : array
        16th percentile density profile
    rho_p84 : array
        84th percentile density profile
    """

    halo_masses = pig[f'FOFGroups/{mass_def}'][:]

    # if groupids is provided, use it to select halos, and ignore m_lower and m_upper
    if groupids is not None:
        mask_idxs = groupids - 1
        mask = np.zeros_like(halo_masses, dtype=bool)
        mask[mask_idxs] = True
    else:
        # mask halos in mass range (stored in 1e10 Msun/h)
        mask = (halo_masses >= m_lower / 1e10) & (halo_masses < m_upper / 1e10)

        if np.sum(mask) == 0:
            print(f'Warning: No halos found in mass range {m_lower:.2e}--{m_upper:.2e} Msun/h')
            return None, None, None, None

    # get individual halo profiles
    r_bin_mid, density_prof = get_rho_profile(pig)

    selected_density_prof = density_prof[mask, :]
    r_bin_mid = r_bin_mid[mask, :]

    # define common radial range: 0, largest of min r_min, smallest of max r_max
    r_nonzero_lower = np.max(r_bin_mid[:, 1])  # exclude r=0
    r_nonzero_upper = np.min(r_bin_mid[:, -1])
    # common r grid: 0, logspace between r_nonzero_lower and r_nonzero_upper
    n_bin = r_bin_mid.shape[1]
    r_common_nonzero = np.logspace(np.log10(r_nonzero_lower), np.log10(r_nonzero_upper), n_bin - 1)
    r_common = np.concatenate(([0.0], r_common_nonzero))

    # interpolate all profiles to common r grid (log-log)
    n_halo = selected_density_prof.shape[0]
    n_bin = len(r_common)
    density_prof_interp = np.zeros((n_halo, n_bin))

    for i in range(n_halo):
        density_prof_interp[i, :] = loglog_interp(
            r_bin_mid[i, :],
            selected_density_prof[i, :],
            r_common
        )
        # for debugging
        # print something if density_prof_interp[i, 4] is nan
        if np.isnan(density_prof_interp[i, 4]):
            print(f'Warning: NaN in interpolated profile for halo {i}')
            print(f'r_bin_mid: {r_bin_mid[i, :]}')
            print(f'selected_density_prof: {selected_density_prof[i, :]}')

    # stack in log space
    log_rho = np.log10(density_prof_interp)

    log_rho_med = np.median(log_rho, axis=0)
    # find the index of the first non-nan value backwards: log_rho_med
    i_nonan = np.where(~np.isnan(log_rho_med))[0][-1]
    # cut off at i_nonan
    log_rho = log_rho[:, :i_nonan+1]
    r_common = r_common[:i_nonan+1]
    
    # log_rho_p16 = np.percentile(log_rho, 16, axis=0)
    # log_rho_p84 = np.percentile(log_rho, 84, axis=0)

    log_rho_med = np.nanmedian(log_rho, axis=0)
    log_rho_p16 = np.nanpercentile(log_rho, 16, axis=0)
    log_rho_p84 = np.nanpercentile(log_rho, 84, axis=0)

    # back to linear space
    rho_med = 10**log_rho_med
    rho_p16 = 10**log_rho_p16
    rho_p84 = 10**log_rho_p84

    return r_common, rho_med, rho_p16, rho_p84

# get ratio
def get_profile_ratio(r1, rho1, r2, rho2):
    # if None is passed in, return None
    if r1 is None or rho1 is None or r2 is None or rho2 is None:
        return None, None
    # find the first index where rho == 0 and truncate
    i_cut1 = np.where(rho1 <= 0)[0][0] if np.any(rho1 <= 0) else len(rho1)
    i_cut2 = np.where(rho2 <= 0)[0][0] if np.any(rho2 <= 0) else len(rho2)

    r1 = r1[:i_cut1]
    rho1 = rho1[:i_cut1]
    r2 = r2[:i_cut2]
    rho2 = rho2[:i_cut2]

    # interpolate to the smaller r grid
    r_common = r1 if r1.max() < r2.max() else r2

    # interpolate in log-log space
    rho1_interp = loglog_interp(r1, rho1, r_common)
    rho2_interp = loglog_interp(r2, rho2, r_common)

    ratio = rho2_interp / rho1_interp
    return r_common, ratio

def loglin_interp(r_src, y_src, r_tgt):
        # exclude the first point where r=0 and then attach it back later
        if r_src[0] == 0:
            r_src_nonzero = r_src[1:]
            y_src_nonzero = y_src[1:]
            attach_zero = True
        else:
            r_src_nonzero = r_src
            y_src_nonzero = y_src
            attach_zero = False

        f = interp1d(np.log(r_src_nonzero), y_src_nonzero, kind='linear', fill_value='extrapolate')

        if r_tgt[0] == 0 and attach_zero:
            y_tgt_nonzero = f(np.log(r_tgt[1:]))
            y_tgt = np.concatenate(([y_src[0]], y_tgt_nonzero))
        else:
            y_tgt = f(np.log(r_tgt))

        return y_tgt

def get_matched_profile_ratio(pig1, pig2, gids1, gids2, m_lower=1e13, m_upper=1e14, mass_def='M500c'):
    # gids1 and gids2 are the matched group ids between pig1 and pig2, in units of GroupID (starting from 1)
    # get density profiles for matched groups and calculate ratio
    gidx1 = gids1 - 1  # convert to index
    gidx2 = gids2 - 1

    r1, rho1 = get_rho_profile(pig1)
    r2, rho2 = get_rho_profile(pig2)

    r1, rho1 = r1[gidx1, :], rho1[gidx1, :]
    r2, rho2 = r2[gidx2, :], rho2[gidx2, :]

    halo_masses1 = pig1[f'FOFGroups/{mass_def}'][:][gidx1]  # consider only pig 1 for mass selection
    mask = (halo_masses1 >= m_lower / 1e10) & (halo_masses1 < m_upper / 1e10)

    # ratio for each matched group
    rs = []
    ratios = []
    for i in range(len(gidx1)):
        if not mask[i]:
            continue
        r_common, ratio = get_profile_ratio(r1[i, :], rho1[i, :], r2[i, :], rho2[i, :])
        rs.append(r_common)
        ratios.append(ratio)
    
    # return r_common and the median of ratios
    if len(ratios) == 0:
        return None, None, None, None
    # stack ratios in log space

    # interpolate all ratios to the common r grid defined by r_min_max (except 0) and r_max_min
    r_min_max = max([r[1] for r in rs])
    r_max_min = min([r[-1] for r in rs])
    r_common = np.logspace(np.log10(r_min_max), np.log10(r_max_min), len(rs[0]) - 1)
    r_common = np.concatenate(([0.0], r_common))

    # interpolate in log-linear space
    ratios_interp = []
    for i, ratio in enumerate(ratios):
        ratios_interp.append(loglin_interp(rs[i], ratio, r_common))

    ratios_interp = np.array(ratios_interp)
    ratio_med = np.median(ratios_interp, axis=0)
    ratio_p16 = np.percentile(ratios_interp, 16, axis=0)
    ratio_p84 = np.percentile(ratios_interp, 84, axis=0)

    return r_common, ratio_med, ratio_p16, ratio_p84

# plot stacked density profile
def plot_stacked(ax, r, rho_median, rho_p16, rho_p84, label_sim, color, ls='-', alpha=0.3):
    ax.fill_between(r, rho_p16, rho_p84, alpha=alpha, color=color)
    ax.plot(r, rho_median, label=label_sim, color=color, linestyle=ls)