import glob
import numpy as np
import os
import argparse
from bigfile import BigFile

def get_bhhaloinfo(pig):
    t_pig = pig['Header'].attrs['Time'][0]

    # get halo properties for these BHs: ID, Mass, M200c, M200m, M500c
    ids_bh = pig['5/ID'][:] 
    ids_bhgroup = pig['5/GroupID'][:]

    idces_bhgroup = ids_bhgroup - 1  # group index in pig (assuming GroupID starts from 1)

    # ensure the properties exist
    properties_needed = ['FOFGroups/Mass', 'FOFGroups/M200c', 'FOFGroups/M200m', 'FOFGroups/M500c']
    for prop in properties_needed:
        if prop not in pig.keys():
            raise KeyError(f"Property {prop} not found in PIG file!")

    haloinfo = {
        'BHID': ids_bh,
        'GroupID': ids_bhgroup,
        'GroupMass': pig['FOFGroups/Mass'][:][idces_bhgroup],
        # 'M200c': pig['FOFGroups/M200c'][idx_kgroup],
        # 'M200m': pig['FOFGroups/M200m'][idx_kgroup],
        # 'M500c': pig['FOFGroups/M500c'][idx_kgroup],
    }
    # M200c, M200m, M500c don't have the same length as Mass, so create a zero array first
    for prop_name in properties_needed:
        if prop_name == 'FOFGroups/Mass':
            continue
        prop = np.zeros(pig['FOFGroups/Mass'].size)
        prop[:pig[prop_name].size] = pig[prop_name][:]
        haloinfo[prop_name.split('/')[-1]] = prop[idces_bhgroup]  # such that halo mass is set to 0 if not available

    # also save the time
    haloinfo['Time'] = t_pig

    return haloinfo

# example usage
# python extract_bhhalo.py --pig /path/to/pigfile --output bhhalo.npz

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract BH halo properties from PIG_subfind.")
    parser.add_argument("--pig", type=str, help="PIG file containing the halo catalog.")
    # parser.add_argument("--kbhfile", type=str, help="Input KBH info file (e.g., kbhinfo.npz).")
    parser.add_argument("--output", default="bhhalo.npz", type=str, help="Output file to save extracted BH halo properties.")
    args = parser.parse_args()

    # print all arguments
    print("Arguments:")
    for arg in vars(args):
        print(f"  {arg}: {getattr(args, arg)}")

    pig = BigFile(args.pig)

    haloinfo = get_bhhaloinfo(pig)

    # save to output file as npz
    print(f"Saving BH halo properties to {args.output}...")
    np.savez(args.output, **haloinfo)
    print(f"Saved BH halo properties to {args.output}.")