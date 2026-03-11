# import glob
import numpy as np
# import os
import argparse
from bigfile import BigFile

def get_swallow_count(part_snap,):
    # get the current time
    # current_time = part_snap['Header'].attrs['Time'][0]

    # get BHIDs and their SwallowIDs from PART
    bhids = part_snap['5/ID'][:]
    swallowids = part_snap['5/BlackholeSwallowID'][:]
    # swallow time
    swallowtimes = part_snap['5/BlackholeSwallowTime'][:]
    # sort by swallow time ascending
    order = np.argsort(swallowtimes)
    bhids = bhids[order]
    swallowids = swallowids[order]
    swallowtimes = swallowtimes[order]

    # BHs existing at the current time
    mask_current = swallowtimes == np.inf

    # ensure mask_current is exactly at the end
    assert np.all(swallowtimes[:-np.sum(mask_current)] < np.inf) and np.all(swallowtimes[-np.sum(mask_current):] == np.inf), "swallowtimes array is not properly sorted."

    # ensure the BH count matches those not swallowed
    assert np.sum(mask_current) == np.sum(part_snap['5/Swallowed'][:] == 0), "Mismatch in number of BHs not swallowed."

    bhids_current = bhids[mask_current]

    swallowcount = np.zeros_like(bhids, dtype=np.uint32)

    # add swallowed BHs to their swallowers until reaching a BH that is not swallowed
    for i in range(len(bhids)):
        if swallowtimes[i] == np.inf:
            break  # stop if reaching a BH that is not swallowed

        swallower_id = swallowids[i]
        # find the index of the swallower in bhids (not necessarily in bhids_current)
        idx_swallower = np.where(bhids == swallower_id)[0]
        if len(idx_swallower) == 0:
            raise ValueError(f"Swallower ID {swallower_id} not found in bhids.")
        # increment the swallow count of the swallower
        swallowcount[idx_swallower[0]] += 1 + swallowcount[i]  # add the count of the swallowed BH as well

    swallowcount_current = swallowcount[mask_current]

    # only return current BHs
    out = {
        'BHID': bhids_current,
        'SwallowCount': swallowcount_current,
        'Time': part_snap['Header'].attrs['Time'][0],
    }

    return out

# example usage
# python extract_swallowed.py --part /path/to/PART_snapshot --output bhswallowcount.npz

# count number of swallowed BHs for each BH existing at the current time

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract KBH info from BlackholeDetails.")
    parser.add_argument("--part", type=str, help="PART snapshot file (ensure the time is not before the target cut-off).")
    parser.add_argument("--output", type=str, help="Output file to save extracted swallow info.")  # bhswallowcount.npz
    args = parser.parse_args()

    # print all arguments
    print("Arguments:")
    for arg in vars(args):
        print(f"  {arg}: {getattr(args, arg)}")

    part_snap = BigFile(args.part)

    swallowinfo = get_swallow_count(part_snap,)

    # save to output file as npz
    print(f"Saving swallow info to {args.output}...")
    np.savez(args.output, **swallowinfo)
    print(f"Saved swallow info to {args.output}.")