import glob
import numpy as np
import os
import argparse
from bigfile import BigFile

def match_groups(pig1, pig2, gid1, gid2, mass_tol=0.15, pos_tol=100):
    """
    Match groups in two simulations based on mass and positions.

    Parameters:
    - pig1, pig2: BigFile objects for the two simulations.
    - gid1, gid2: Arrays of group IDs to match from pig1 and pig2.
    - mass_tol: Tolerance for mass matching (fractional).
    - pos_tol: Tolerance for position matching (in kpc/h).

    Returns:
    - matched_gids: Array of matched group IDs (gid1, gid2), -1 if no match found.
    """

    # if gid1 or gid2 is None, use all groups in pig1 or pig2
    if gid1 is None:
        gid1 = np.arange(1, pig1['FOFGroups/M200c'].size + 1)  # GroupIDs are 1-based
    if gid2 is None:
        gid2 = np.arange(1, pig2['FOFGroups/M200c'].size + 1)  # GroupIDs are 1-based
    gidx1 = gid1 - 1  # convert to zero-based index
    gidx2 = gid2 - 1
    # keep only gidx >= 0 (gidx == -1 means no group)
    gidx1 = gidx1[gidx1 >= 0]
    gidx2 = gidx2[gidx2 >= 0]

    len_avail_pig1 = pig1['FOFGroups/M200c'].size # only consider groups with available M200c
    len_avail_pig2 = pig2['FOFGroups/M200c'].size
    gidx1 = gidx1[gidx1 < len_avail_pig1]
    gidx2 = gidx2[gidx2 < len_avail_pig2]
    # also exclude those with zero mass
    gidx1 = gidx1[pig1['FOFGroups/M200c'][:][gidx1] > 0]
    gidx2 = gidx2[pig2['FOFGroups/M200c'][:][gidx2] > 0]

    # print gid sizes
    print(f"Number of groups in gid1: {len(gidx1)}")
    print(f"Number of groups in gid2: {len(gidx2)}")

    # loop over gid1 and find matches in gid2
    gidx1_matched = []
    gidx2_matched = []
    merit_matched = []
    dist_matched = []

    for i in range(len(gidx1)):
        # mass: keep only those with similar mass
        print(f"Matching group {i+1}/{len(gidx1)}: Group idx {gidx1[i]}")
        # print('M200c size', pig1['FOFGroups/M200c'].size, pig2['FOFGroups/M200c'].size)
        mass1 = pig1['FOFGroups/M200c'][gidx1[i]]
        gidx2_available = gidx2
        mass2_all = pig2['FOFGroups/M200c'][:][gidx2_available]

        mass_mask = (mass2_all >= mass1 * (1 - mass_tol)) & (mass2_all <= mass1 * (1 + mass_tol))

        if np.sum(mass_mask) == 0:
            print("  No mass match found.")
            continue  # no match found
        
        subidx1 = pig1['FOFGroups/GroupFirstSub'][gidx1[i]]
        
        # filter by distance before computing merit, to save time
        pos1 = pig1['SubGroups/SubhaloPos'][subidx1]
        subidx2_all = pig2['FOFGroups/GroupFirstSub'][:][gidx2_available[mass_mask]]
        # if any subidx2_all is out of bounds, raise error
        if np.any(subidx2_all >= pig2['SubGroups/SubhaloPos'].size):
            # find the first out of bounds index
            first_out_of_bounds = np.where(subidx2_all >= pig2['SubGroups/SubhaloPos'].size)[0][0]
            raise ValueError(f"Subhalo index out of bounds in pig2. the gidx2 is {gidx2_available[mass_mask][first_out_of_bounds]}")
        pos2_all = pig2['SubGroups/SubhaloPos'][:][subidx2_all]

        dists = np.linalg.norm(pos2_all - pos1, axis=1)

        mask_pos = dists <= pos_tol

        if np.sum(mask_pos) == 0:
            print("  No position match found within tolerance.")
            continue  # no match found

        # reserve only those that pass both mass and position mask
        gidx2_available = gidx2_available[mass_mask][mask_pos]

        # merit function: find the halo with largest fractional particle overlap
        # M = Nij / (Ni * Nj)
        # where Nij is the number of shared particles, Ni and Nj are the total number of particles in each halo
        # consider only DM particles (type 1), the first subhalo in each group

        pstart1 = pig1['SubGroups/SubhaloOffsetType'][subidx1][1]
        pend1 = pstart1 + pig1['SubGroups/SubhaloLenType'][subidx1][1]
        part_ids1 = pig1['1/ID'][pstart1:pend1]

        max_merit = -1
        best_j = -1

        for j in range(len(gidx2_available)):
            subidx2 = pig2['FOFGroups/GroupFirstSub'][gidx2_available[j]]

            pstart2 = pig2['SubGroups/SubhaloOffsetType'][subidx2][1]
            pend2 = pstart2 + pig2['SubGroups/SubhaloLenType'][subidx2][1]
            part_ids2 = pig2['1/ID'][pstart2:pend2]

            # compute intersection
            shared_ids = np.intersect1d(part_ids1, part_ids2)
            Nij = len(shared_ids)
            Ni = len(part_ids1)
            Nj = len(part_ids2)

            merit = Nij**2 / (Ni * Nj)

            if merit > max_merit:
                max_merit = merit
                best_j = j

         # print merit
        print(f"  Best merit: {max_merit:.4f}.")
        # skip if max_merit is < .5
        if max_merit < 0.5:
            print("  No good merit match found, skipping.")
            continue

        best_gidx2_merit = gidx2_available[best_j]

        # position: find the closest in position and check within tolerance
        
        # pos1 = pig1['SubGroups/SubhaloPos'][subidx1]

        # subidx2_all = pig2['FOFGroups/GroupFirstSub'][:][gidx2_available[mass_mask]]
        # # if any subidx2_all is out of bounds, raise error
        # if np.any(subidx2_all >= pig2['SubGroups/SubhaloPos'].size):
        #     # find the first out of bounds index
        #     first_out_of_bounds = np.where(subidx2_all >= pig2['SubGroups/SubhaloPos'].size)[0][0]
        #     raise ValueError(f"Subhalo index out of bounds in pig2. the gidx2 is {gidx2_available[mass_mask][first_out_of_bounds]}")
        # pos2_all = pig2['SubGroups/SubhaloPos'][:][subidx2_all]

        # dists = np.linalg.norm(pos2_all - pos1, axis=1)

        # min_dist_idx = np.argmin(dists)
        
        # best_gidx2_pos = gidx2_available[mass_mask][min_dist_idx]

        # # check if it is corresponding to the best_j found from merit
        # if best_gidx2_pos != best_gidx2_merit:
        #     print("  Best position match does not correspond to best merit match, skipping.")
        #     continue  # skip this match for consistency

        # calculate distance between pos1 and best merit match
        subidx2_merit = pig2['FOFGroups/GroupFirstSub'][best_gidx2_merit]
        pos2_merit = pig2['SubGroups/SubhaloPos'][subidx2_merit]
        dist_merit = np.linalg.norm(pos2_merit - pos1)

        gidx1_matched.append(gidx1[i])
        gidx2_matched.append(best_gidx2_merit)
        merit_matched.append(max_merit)
        dist_matched.append(dist_merit)
        print(f"  Matched with group idx {best_gidx2_merit} at distance {dist_merit:.2f} kpc/h.")

        assert dist_merit <= pos_tol, f"Matched group is beyond position tolerance: {dist_merit:.2f} kpc/h."

    gidx1_matched = np.array(gidx1_matched)
    gidx2_matched = np.array(gidx2_matched)
    merit_matched = np.array(merit_matched)
    dist_matched = np.array(dist_matched)
    print(f"Total matched groups: {len(gidx1_matched)}")
    # build output array
    gid1_matched = gidx1_matched + 1  # convert back to GroupID
    gid2_matched = gidx2_matched + 1  # convert back to GroupID
    
    # combine into structured array
    matched_gids = np.zeros(len(gid1_matched), dtype=[('GroupID1', np.uint64), ('GroupID2', np.uint64), ('Merit', np.float64), ('Distance', np.float64)])
    matched_gids['GroupID1'] = gid1_matched
    matched_gids['GroupID2'] = gid2_matched
    matched_gids['Merit'] = merit_matched
    matched_gids['Distance'] = dist_matched

    return matched_gids

# example usage
# python match_group.py --pig1 pig1 --pig2 pig2 --gid1 groupid_primaryswallowcount0_kn1S-Repos1.npy --gid2 groupid_primaryswallowcount0_kn1S-DF.npy

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Match groups in two simulations based on mass and positions.")
    parser.add_argument("--pig1", type=str, help="PIG file 1 containing the halo catalog.")
    parser.add_argument("--pig2", type=str, help="PIG file 2 containing the halo catalog.")
    parser.add_argument("--gid1", type=str, default=None, help="Group ID file (e.g., groupid_primaryswallowcount0_kn1S-Repos1.npy).")
    parser.add_argument("--gid2", type=str, default=None, help="Group ID file (e.g., groupid_primaryswallowcount0_kn1S-DF.npy).")
    # mass_tol and pos_tol can be adjusted if needed
    parser.add_argument("--masstol", default=0.15, type=float, help="Mass tolerance for matching groups.")
    parser.add_argument("--postol", default=100, type=float, help="Position tolerance (kpc/h) for matching groups.")
    parser.add_argument("--output", default="matched_gids.npy", type=str, help="Output file to save matched group IDs.")
    args = parser.parse_args()

    # print all arguments
    print("Arguments:")
    for arg in vars(args):
        print(f"  {arg}: {getattr(args, arg)}")

    pig1 = BigFile(args.pig1)
    pig2 = BigFile(args.pig2)

    gid1 = np.load(args.gid1) if args.gid1 is not None else None
    gid2 = np.load(args.gid2) if args.gid2 is not None else None

    print("Matching groups...")
    # match group IDs
    matched_gids = match_groups(pig1, pig2, gid1, gid2, mass_tol=args.masstol, pos_tol=args.postol)

    # save to output file as npy
    print(f"Saving matched group IDs to {args.output}...")
    np.save(args.output, matched_gids)
    print(f"Saved matched group IDs to {args.output}.")