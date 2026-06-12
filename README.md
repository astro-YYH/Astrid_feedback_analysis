This repository contains analysis code and notebooks for studying black hole feedback, host halos, and halo density profiles in cosmological simulations run with MP-Gadget (https://github.com/MP-Gadget/MP-Gadget). The code was used in Yang et al. 2026.

The code is organized into two workflows:

1. `mass_profile/`: compute cumulative halo mass profiles from `PIG_*_subfind` catalogs and compare stacked density profiles between matched halos.
2. `black_hole/`: extract black hole evolution information (feedback histories in particular), attach halo properties, and make summary figures using halo quantities produced by the former, such as `M500c`.

## Requirements

The scripts use standard scientific Python packages together with simulation-specific I/O:

- `numpy`
- `scipy`
- `matplotlib`
- `astropy`
- `mpi4py`
- `bigfile`
- Jupyter for the notebooks

Several scripts assume access to MP-Gadget `PART`, `PIG_*`, and `BlackholeDetails*` outputs stored outside this repository.

## Repository Layout

### `mass_profile/`

Code for halo mass-profile extraction and density-profile comparison.

#### Python scripts

- `mass_profile.py`
	MPI-enabled production script that reads a `PIG_*_subfind` catalog and computes cumulative mass profiles for each resolved halo.

	For each halo above a DM particle threshold, it computes:
	- `M200m`, `R200m`
	- `M200c`, `R200c`
	- `M500c`, `R500c`
	- cumulative mass profiles for each available particle type

	The output is written back as BigFile blocks under `FOFGroups/`, so the resulting catalog can be analyzed without recomputing profiles.

	Example:

	```bash
	mpirun -n 4 python mass_profile.py \
		--pigdir /path/to/PIG_021_subfind \
		--savedir /path/to/output_bigfile \
		--mindmpart 300 \
		--maxradrel 3.0 \
		--minradabs 5 \
		--nrbin 100 \
		--mnu 0.06
	```

- `match_group.py`
	Matches halos between two simulations. Candidate matches are filtered by halo mass and central-subhalo position, then ranked using a dark-matter particle-overlap merit function:

	$$
	M = \frac{N_{ij}^2}{N_i N_j}
	$$

	The output is a structured `.npy` array with `GroupID1`, `GroupID2`, `Merit`, and `Distance`.

	This is used when comparing matched systems across two runs, for example between DF and Repos variants.

	Example:

	```bash
	python match_group.py \
		--pig1 /path/to/run1/PIG_020_subfind \
		--pig2 /path/to/run2/PIG_020_subfind \
		--gid1 /path/to/group_ids_run1.npy \
		--gid2 /path/to/group_ids_run2.npy \
		--masstol 0.15 \
		--postol 100 \
		--output matched_gids.npy
	```

- `profile_funcs.py`
	Analysis helpers for the profile workflow. It provides functions to:
	- convert cumulative mass profiles into shell density profiles,
	- stack profiles in a halo-mass bin,
	- interpolate profiles onto a common radial grid,
	- compute hydro-to-DMO profile ratios,
	- combine matched halo lists into stacked comparison curves.

#### Jupyter notebooks

- `density_mathced_profile_ratio.ipynb`
	Loads precomputed profile catalogs and matched halo lists, then compares stacked density-profile ratios between hydro and DMO runs.

	The notebook currently focuses on `z = 0` and plots three halo-mass bins, showing
	the ratio of hydro to DMO density profiles as a function of radius for four simulation variants, with percentile bands and a reference line marking the median DMO `R500c` in each mass bin.

### `black_hole/`

Code for extracting kinetic-feedback black hole catalogs and relating them to halo properties.

#### Python scripts

- `extract_kbhinfo.py`
	Reads one or more `BlackholeDetails*` directories, extracts a reduced black hole history catalog, and writes an `.npz` file containing:
	`BHID`, final `BHMass`, `FirstReleaseTime`, `FirstReleaseMass`, total `ReleasedKE`, `EndTime`, and `SwallowID`.

	Main logic:
    - truncates the catalog at a target scale factor `--time`,
	- follows swallow chains using a matching `PART` snapshot,
	- keeps enough swallowed-BH entries to reconstruct total feedback inherited by surviving BHs later.

	Example:

	```bash
	python extract_kbhinfo.py \
		--inputdir /path/to/output \
		--basename BlackholeDetails \
		--part /path/to/PART_021 \
		--output kbhinfo_kn1-DF.npz \
		--time 1.0
	```

- `extract_bhhalo.py`
	Reads a `PIG_*_subfind` halo catalog and stores host-halo information for every BH particle in an `.npz` file. The output includes `BHID`, `GroupID`, `GroupMass`, `M200c`, `M200m`, `M500c`, and the snapshot scale factor `Time`.

	This is the bridge between BH-centric catalogs and halo-mass-based plots.

	Example:

	```bash
	python extract_bhhalo.py \
		--pig /path/to/PIG_021_subfind \
		--output bhhaloinfo_kn1-DF.npz
	```

- `utils.py`
	Shared utilities for the black-hole analysis and notebooks. In particular it:
	- converts masses and feedback energy into plotting units,
	- aggregates feedback from swallowed BHs onto the surviving primary BH,
	- loads `kbhinfo_*` and `bhhaloinfo_*` products across multiple redshifts,
	- plots kinetic-feedback energy density distributions as a function of BH or halo mass.

#### Jupyter notebooks

- `BH_halo_relation.ipynb`
	Loads precomputed `kbhinfo_*` and `bhhaloinfo_*` products for four simulations and makes a 2x4 comparison plot at `z = 0` and `z = 1`.

	The notebook visualizes BH mass versus halo mass, with point color showing released kinetic feedback energy. Non-kinetic BHs are plotted separately in gray.

- `KE_halo_times.ipynb`
	Uses the processed BH catalogs to plot the kinetic-feedback energy density distribution
	$\mathcal{U}_\mathrm{kf}^\mathrm{halo}(M_{500c}, z)$
	for multiple simulation variants.

	The current notebook configuration compares `z = 0` and `z = 1`, annotating total energy density for each model on the figure.

#### Data products included here

The `.npz` files in this directory are precomputed black hole and BH-host-halo catalogs for several simulation variants and redshifts. They are used directly by the notebooks.

## Typical Workflow

1. Run `mass_profile.py` on halo catalogs to create cumulative profile blocks.
2. Use `mass_profile/match_group.py` and `mass_profile/profile_funcs.py` to compare matched halos across simulations.
3. Run `extract_kbhinfo.py` on `BlackholeDetails*` plus a matching `PART` snapshot.
4. Run `extract_bhhalo.py` on the corresponding `PIG_*_subfind` halo catalog.
5. Use `black_hole/utils.py` and the notebooks to build BH-halo and energy-density figures.

## Notes

- The repository contains both code and precomputed intermediate products.
- Several notebooks assume specific filenames already exist in the working directory.

## License

This software is distributed under the MIT License. See `LICENSE` for details.
