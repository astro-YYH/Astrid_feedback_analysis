#!/bin/bash 
#SBATCH -p small
#SBATCH -A AST21005
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=kbhinfo
#SBATCH --output=%x-%j.out

date

# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-DF/output --output kbhinfo_kn1-DF.npz --part /scratch3/01317/yyang440/small_Astrid/kn1-DF/output/PART_021
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output --output kbhinfo_kn1-Repos1.npz --oldformat --part /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output/PART_021
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output --output kbhinfo_kn1S-Repos1.npz --part /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output/PART_020
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1S-DF/output --output kbhinfo_kn1S-DF.npz --part /scratch3/01317/yyang440/small_Astrid/kn1S-DF/output/PART_020

# z=0.2  for time, keep 6 decimal places to avoid precision issues
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-DF/output --output kbhinfo_kn1-DF_z0.2.npz --time 0.833333 --part /scratch3/01317/yyang440/small_Astrid/kn1-DF/output/PART_021
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output --output kbhinfo_kn1-Repos1_z0.2.npz --time 0.833333 --oldformat --part /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output/PART_021
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output --output kbhinfo_kn1S-Repos1_z0.2.npz --time 0.833333 --part /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output/PART_020
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1S-DF/output --output kbhinfo_kn1S-DF_z0.2.npz --time 0.833333 --part /scratch3/01317/yyang440/small_Astrid/kn1S-DF/output/PART_020

# z=0.5
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-DF/output --output kbhinfo_kn1-DF_z0.5.npz --time 0.666667 --part /scratch3/01317/yyang440/small_Astrid/kn1-DF/output/PART_021
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output --output kbhinfo_kn1-Repos1_z0.5.npz --time 0.666667 --oldformat --part /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output/PART_021
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output --output kbhinfo_kn1S-Repos1_z0.5.npz --time 0.666667 --part /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output/PART_020
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1S-DF/output --output kbhinfo_kn1S-DF_z0.5.npz --time 0.666667 --part /scratch3/01317/yyang440/small_Astrid/kn1S-DF/output/PART_020

# z=1.0
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-DF/output --output kbhinfo_kn1-DF_z1.0.npz --time 0.5 --part /scratch3/01317/yyang440/small_Astrid/kn1-DF/output/PART_021
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output --output kbhinfo_kn1-Repos1_z1.0.npz --time 0.5 --oldformat --part /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output/PART_021
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output --output kbhinfo_kn1S-Repos1_z1.0.npz --time 0.5 --part /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output/PART_020
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1S-DF/output --output kbhinfo_kn1S-DF_z1.0.npz --time 0.5 --part /scratch3/01317/yyang440/small_Astrid/kn1S-DF/output/PART_020

# z=1.5 (a=0.4)
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1S-DF/output --output kbhinfo_kn1S-DF_z1.5.npz --time 0.4 --part /scratch3/01317/yyang440/small_Astrid/kn1S-DF/output/PART_020
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-DF/output --output kbhinfo_kn1-DF_z1.5.npz --time 0.4 --part /scratch3/01317/yyang440/small_Astrid/kn1-DF/output/PART_021
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output --output kbhinfo_kn1-Repos1_z1.5.npz --time 0.4 --oldformat --part /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output/PART_021
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output --output kbhinfo_kn1S-Repos1_z1.5.npz --time 0.4 --part /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output/PART_020

# z=2.3
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-DF/output --output kbhinfo_kn1-DF_z2.3.npz --time 0.3030


# z=2.0
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output --output kbhinfo_kn1-Repos1_z2.0.npz --time 0.333333 --oldformat --part /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output/PART_021
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-DF/output --output kbhinfo_kn1-DF_z2.0.npz --time 0.333333 --part /scratch3/01317/yyang440/small_Astrid/kn1-DF/output/PART_021
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output --output kbhinfo_kn1S-Repos1_z2.0.npz --time 0.333333 --part /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output/PART_020
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1S-DF/output --output kbhinfo_kn1S-DF_z2.0.npz --time 0.333333 --part /scratch3/01317/yyang440/small_Astrid/kn1S-DF/output/PART_020

# extract halo info
# z=2.0
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1-DF/009/PIG_009_subfind --output bhhaloinfo_kn1-DF_z2.0.npz
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1-Repos1/008/PIG_008_subfind --output bhhaloinfo_kn1-Repos1_z2.0.npz

# swallowcount0 matching
# z = 0
python match_group.py --pig1 /scratch3/01317/yyang440/small_subfind/kn1S-Repos1/020/PIG_020_subfind --pig2 /scratch3/01317/yyang440/small_subfind/kn1S-DF/020/PIG_020_subfind --gid1 groupid_primaryswallowcount0_kn1S-Repos1.npy --gid2 groupid_primaryswallowcount0_kn1S-DF.npy

date