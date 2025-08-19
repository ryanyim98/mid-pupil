#!/bin/bash


# modified by Ryan from Leili's script

# run from the scripts directory

SCRIPTS_DIR=$(pwd)

cd ..
MAIN_DIR=$(pwd)

IN_DIR=$MAIN_DIR/glm_results/individuals_MIDaffemo
OUT_DIR=$MAIN_DIR/glm_results/group_MIDaffemo

N=30

SUBJECTS=$(cat $MAIN_DIR/data/subjects-midaffemo-included.txt)

MASK=$MAIN_DIR/ROIs/glm_mask_gm.nii.gz
RESAMP_DIM=2.9

OUT_FILES_model27=('param_arous')
STATS_model27=('27')
# if wm and csf are included as non-nuisance regressors
# STATS_model27=('11' '15' '19' '23' '27' '31' '35')

REG_FILE=glm27_z_b4_mni.nii.gz

cd $OUT_DIR

for i in ${!OUT_FILES_model27[@]}
do

        out_file=${OUT_FILES_model27[$i]}
        stat=${STATS_model27[$i]}
        echo making $out_file
        echo ----- stat: $stat

        # Create the glm list
        glm_sub_list=()
        for sub in $SUBJECTS; do
            glm_sub_list+=("$IN_DIR/${sub}_${REG_FILE}[$stat]")
        done

        3dttest++ -overwrite -mask $MASK -prefix model27_${out_file} -setA "${glm_sub_list[@]}"

        3dmerge -doall -1zscore -prefix model27_${out_file}_${N}z model27_${out_file}+tlrc
        rm -f model27_${out_file}+tlrc.*

        3dAFNItoNIFTI -prefix model27_${out_file}_${N}z.nii.gz model27_${out_file}_${N}z+tlrc
        rm -f model27_${out_file}_${N}z+tlrc.*

done
