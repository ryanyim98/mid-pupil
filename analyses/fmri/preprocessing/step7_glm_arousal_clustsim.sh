#!/bin/bash

SCRIPTS_DIR=$(pwd)
cd ..
MAIN_DIR=$(pwd)

IN_DIR=$MAIN_DIR/glm_results/individuals_MIDaffemo
OUT_DIR=$MAIN_DIR/glm_results/group_MIDaffemo_clustersim

N=30
SUBJECTS=$(cat $MAIN_DIR/data/subjects-midaffemo-included.txt)
MASK=$MAIN_DIR/ROIs/glm_mask_gm.nii.gz
# MASK=$MAIN_DIR/ROIs/anteriorinsula8mmkg_func.nii.gz
# DILATED_MASK="${MAIN_DIR}/ROIs/anteriorinsula8mmkg_func_dilated.nii.gz"

RESAMP_DIM=2.9

OUT_FILES_model27=('param_arous_early')
STATS_model27=('26')
REG_FILE=glm27c_b4_mni.nii.gz

cd $OUT_DIR

for i in ${!OUT_FILES_model27[@]}
do
    out_file=${OUT_FILES_model27[$i]}
    stat=${STATS_model27[$i]}
    echo "Making $out_file (stat: $stat)"

    # Create the glm list
    glm_sub_list=()
    for sub in $SUBJECTS; do
        glm_sub_list+=("$IN_DIR/${sub}_${REG_FILE}[$stat]")
    done

    # Run ClustSim and ETAC together
    3dttest++ -overwrite -mask $MASK -prefix model27c_${out_file} \
              -Clustsim \
              -ETAC \
              -ETAC_opt NN=2:fpr=0.05:pthr=0.001,0.005,0.01:onesided \
              -setA "${glm_sub_list[@]}"

    # Process ClustSim output
    3dmerge -doall -1zscore -prefix model27c_${out_file}_${N}z model27c_${out_file}+tlrc
    rm -f model27c_${out_file}+tlrc.*

    3dAFNItoNIFTI -prefix model27c_${out_file}_${N}z.nii.gz model27c_${out_file}_${N}z+tlrc
    rm -f model27c_${out_file}_${N}z+tlrc.*
done
