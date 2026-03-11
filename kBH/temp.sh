#!/bin/bash

# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-DF/output --output kbhinfo_kn1-DF.npz 
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output --output kbhinfo_kn1-Repos1.npz --oldformat
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn0-Repos1/output --output kbhinfo_kn0-Repos1.npz --oldformat
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output --output kbhinfo_kn1S-Repos1.npz

# z=0.2
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-DF/output --output kbhinfo_kn1-DF_z0.2.npz --time 0.83333
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output --output kbhinfo_kn1-Repos1_z0.2.npz --time 0.83333 --oldformat
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output --output kbhinfo_kn1S-Repos1_z0.2.npz --time 0.83333

# z=0.5
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-DF/output --output kbhinfo_kn1-DF_z0.5.npz --time 0.66667
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output --output kbhinfo_kn1-Repos1_z0.5.npz --time 0.66667 --oldformat
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output --output kbhinfo_kn1S-Repos1_z0.5.npz --time 0.66667

# z=1.0
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-DF/output --output kbhinfo_kn1-DF_z1.0.npz --time 0.5
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output --output kbhinfo_kn1-Repos1_z1.0.npz --time 0.5 --oldformat
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output --output kbhinfo_kn1S-Repos1_z1.0.npz --time 0.5

# z=2.3
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-DF/output --output kbhinfo_kn1-DF_z2.3.npz --time 0.3030
# python extract_kbhinfo.py --inputdir /scratch3/01317/yyang440/small_Astrid/kn1-DF/output --output kbhinfo_kn1-DF_z2.0.npz --time 0.33333

# python rewrite_BlackholeDetails.py --inputbd /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output/BlackholeDetails --outputbd /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output/BlackholeDetails_new --oldformat --time 0.4

# extract halo info
# kn1-DF
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1-DF/009/PIG_009_subfind --output bhhaloinfo_kn1-DF_z2.0.npz
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1-DF/013/PIG_013_subfind --output bhhaloinfo_kn1-DF_z1.0.npz
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1-DF/015/PIG_015_subfind --output bhhaloinfo_kn1-DF_z0.5.npz
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1-DF/018/PIG_018_subfind --output bhhaloinfo_kn1-DF_z0.2.npz
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1-DF/021/PIG_021_subfind --output bhhaloinfo_kn1-DF.npz

# kn1-Repos1
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1-Repos1/008/PIG_008_subfind --output bhhaloinfo_kn1-Repos1_z2.0.npz
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1-Repos1/013/PIG_013_subfind --output bhhaloinfo_kn1-Repos1_z1.0.npz
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1-Repos1/015/PIG_015_subfind --output bhhaloinfo_kn1-Repos1_z0.5.npz
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1-Repos1/017/PIG_017_subfind --output bhhaloinfo_kn1-Repos1_z0.2.npz
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1-Repos1/021/PIG_021_subfind --output bhhaloinfo_kn1-Repos1.npz

# kn1S-DF
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1S-DF/008/PIG_008_subfind --output bhhaloinfo_kn1S-DF_z2.0.npz
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1S-DF/012/PIG_012_subfind --output bhhaloinfo_kn1S-DF_z1.0.npz
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1S-DF/014/PIG_014_subfind --output bhhaloinfo_kn1S-DF_z0.5.npz
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1S-DF/017/PIG_017_subfind --output bhhaloinfo_kn1S-DF_z0.2.npz
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1S-DF/020/PIG_020_subfind --output bhhaloinfo_kn1S-DF.npz

# kn1S-Repos1
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1S-Repos1/009/PIG_009_subfind --output bhhaloinfo_kn1S-Repos1_z2.0.npz
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1S-Repos1/013/PIG_013_subfind --output bhhaloinfo_kn1S-Repos1_z1.0.npz
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1S-Repos1/015/PIG_015_subfind --output bhhaloinfo_kn1S-Repos1_z0.5.npz
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1S-Repos1/017/PIG_017_subfind --output bhhaloinfo_kn1S-Repos1_z0.2.npz
# python extract_bhhalo.py --pig /scratch3/01317/yyang440/small_subfind/kn1S-Repos1/020/PIG_020_subfind --output bhhaloinfo_kn1S-Repos1.npz

# extract swallow count
# z =0
python extract_swallowed.py --part /scratch3/01317/yyang440/small_Astrid/kn1-DF/output/PART_021 --output bhswallowcount_kn1-DF_z0.0.npz
python extract_swallowed.py --part /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output/PART_021 --output bhswallowcount_kn1-Repos1_z0.0.npz
python extract_swallowed.py --part /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output/PART_020 --output bhswallowcount_kn1S-Repos1_z0.0.npz
python extract_swallowed.py --part /scratch3/01317/yyang440/small_Astrid/kn1S-DF/output/PART_020 --output bhswallowcount_kn1S-DF_z0.0.npz

# z = 0.5
python extract_swallowed.py --part /scratch3/01317/yyang440/small_Astrid/kn1-DF/output/PART_015 --output bhswallowcount_kn1-DF_z0.5.npz
python extract_swallowed.py --part /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output/PART_015 --output bhswallowcount_kn1-Repos1_z0.5.npz
python extract_swallowed.py --part /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output/PART_015 --output bhswallowcount_kn1S-Repos1_z0.5.npz
python extract_swallowed.py --part /scratch3/01317/yyang440/small_Astrid/kn1S-DF/output/PART_014 --output bhswallowcount_kn1S-DF_z0.5.npz

# z= 1
python extract_swallowed.py --part /scratch3/01317/yyang440/small_Astrid/kn1-DF/output/PART_013 --output bhswallowcount_kn1-DF_z1.0.npz
python extract_swallowed.py --part /scratch3/01317/yyang440/small_Astrid/kn1-Repos1/output/PART_013 --output bhswallowcount_kn1-Repos1_z1.0.npz
python extract_swallowed.py --part /scratch3/01317/yyang440/small_Astrid/kn1S-Repos1/output/PART_013 --output bhswallowcount_kn1S-Repos1_z1.0.npz
python extract_swallowed.py --part /scratch3/01317/yyang440/small_Astrid/kn1S-DF/output/PART_012 --output bhswallowcount_kn1S-DF_z1.0.npz


# z= 2
# python extract_swallowed.py --part /scratch3/01317/yyang440/small_Astrid/kn1-DF/output/PART_009 --output bhswallowcount_kn1-DF_z2.0.npz

# swallowcount0 matching
python match_group.py --pig1 /scratch3/01317/yyang440/small_subfind/kn1S-Repos1/020/PIG_020_subfind --pig2 /scratch3/01317/yyang440/small_subfind/kn1S-DF/020/PIG_020_subfind --gid1 groupid_primaryswallowcount0_kn1S-Repos1.npy --gid2 groupid_primaryswallowcount0_kn1S-DF.npy

