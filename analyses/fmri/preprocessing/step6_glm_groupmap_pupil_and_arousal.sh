#!/bin/bash


# modified by Ryan from Leili's script

# run from the scripts directory

SCRIPTS_DIR=$(pwd)

cd ..
MAIN_DIR=$(pwd)

IN_DIR=$MAIN_DIR/glm_results/individuals_MIDaffemo
OUT_DIR=$MAIN_DIR/glm_results/group_MIDaffemo

N=24

SUBJECTS=$(cat $MAIN_DIR/data/subjects-pupil-fmri.txt)

MASK=$MAIN_DIR/ROIs/glm_mask_gm.nii.gz
RESAMP_DIM=2.9

OUT_FILES_model55=('param_arous_early' 'pupil_scaled_convolved')
STATS_model55=('27' '31')
# if wm and csf are included as non-nuisance regressors
# STATS_model55=('11' '15' '19' '23' '27' '31' '35')

REG_FILE=glm55_z_b4_mni.nii.gz

cd $OUT_DIR

for i in ${!OUT_FILES_model55[@]}
do

        out_file=${OUT_FILES_model55[$i]}
        stat=${STATS_model55[$i]}
        echo making $out_file
        echo ----- stat: $stat

        # Create the glm list
        glm_sub_list=()
        for sub in $SUBJECTS; do
            glm_sub_list+=("$IN_DIR/${sub}_${REG_FILE}[$stat]")
        done

        3dttest++ -overwrite -mask $MASK -prefix model55_${out_file} -setA "${glm_sub_list[@]}"

        3dmerge -doall -1zscore -prefix model55_${out_file}_${N}z model55_${out_file}+tlrc
        rm -f model55_${out_file}+tlrc.*

        3dAFNItoNIFTI -prefix model55_${out_file}_${N}z.nii.gz model55_${out_file}_${N}z+tlrc
        rm -f model55_${out_file}_${N}z+tlrc.*

done
