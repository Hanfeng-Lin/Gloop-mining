#!/bin/bash
num_processes=20
# Motif_size=[5,8]; glycine_pos = [3,5]
# If you set the motif_size to 8 then the glycine_pos should be 5
# If you set the motif_size to 5 then the glycine_pos should be 3
motif_size=5
Asp_pos=1
Ser_pos=3
target_motif_pdb_fn=input_template/MOTIFxVAV1xB_${motif_size}res.pdb
alignment_dir=alignment_output_${motif_size}res/
output_dir=out_csv_${motif_size}res_test/
pdb_dir=./pdb_models
af2_dir=./af2_models
mkdir -p $output_dir
mkdir -p $alignment_dir
cat split_job_lists/all_split.list | xargs -P $num_processes -n 1 python3 source/mining.py \
        --Asp_pos $Asp_pos \
        --Ser_pos $Ser_pos \
        --target_motif_pdb_fn $target_motif_pdb_fn \
        --pdb_dir $pdb_dir \
        --af2_dir $af2_dir \
        --temp_output_dir $alignment_dir \
        --output_dir $output_dir >> output_${motif_size}res.log 2>> error_${motif_size}res.log

