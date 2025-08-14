#!/bin/bash
# Motif_size=[5,8]; glycine_pos = [3,5]
# If you set the motif_size to 8 then the glycine_pos should be 5
# If you set the motif_size to 5 then the glycine_pos should be 3
motif_size=5
Asp_pos=1
Ser_pos=3
target_motif_pdb_fn=input_template/MOTIFxVAV1xB_${motif_size}res.pdb
output_dir=out_csv_${motif_size}res/
temp_dir=temp_output_${motif_size}res/
mkdir -p $output_dir
mkdir -p $temp_dir
python3 source/mining.py \
        --Asp_pos $Asp_pos \
        --Ser_pos $Ser_pos \
        --target_motif_pdb_fn $target_motif_pdb_fn \
        --pdb_dir test_data/ \
        --af2_dir af2_domains/ \
        --temp_output_dir $temp_dir \
        --output_dir $output_dir test.list 

