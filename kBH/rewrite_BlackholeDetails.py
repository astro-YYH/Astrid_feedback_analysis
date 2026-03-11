import glob
import numpy as np
import os
import argparse

def write_new_bd(bd_in, bd_out, scale_factor=1.0, oldformat=False):
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

    datadir = bd_in
    bhfile_zoo = glob.glob(os.path.join(datadir,"*"))

    for filename in bhfile_zoo:
        data = np.fromfile(filename, dtype=dt, count=-1)
        # results.append(data[fields_target])

        # extract entries with z <= scale_factor
        mask = data['z'] <= scale_factor
        data = data[mask]

        # write to new file in bd_out with the same filename
        outfile = os.path.join(bd_out, os.path.basename(filename))
        print(f"Writing truncated BlackholeDetails to {outfile}...")
        data.tofile(outfile)
        
        del data

    print("Done.")

# example usage
# python rewrite_BlackholeDetails.py --inputbd /path/to/BlackholeDetails --outputbd /path/to/output --time .4

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Truncate BlackholeDetails to a given scale factor.")
    parser.add_argument("--inputbd", type=str, help="Input BlackholeDetails.")
    parser.add_argument("--outputbd", type=str, help="Output BlackholeDetails.")
    parser.add_argument("--time", type=np.float64, default=1.0, help="Scale factor up to which to truncate BlackholeDetails.")
    parser.add_argument("--oldformat", action='store_true', help="Use old format extraction method, i.e., KEflag is a double instead of int.")
    args = parser.parse_args()

    bd_in = args.inputbd
    bd_out = args.outputbd

    # mkdir -p output directory if not exists
    if not os.path.exists(bd_out):
        os.makedirs(bd_out)

    write_new_bd(bd_in, bd_out, scale_factor=args.time, oldformat=args.oldformat)