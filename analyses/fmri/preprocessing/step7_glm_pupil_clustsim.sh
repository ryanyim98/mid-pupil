#!/bin/bash


# modified by Ryan from Leili's script

# run from the scripts directory

SCRIPTS_DIR=$(pwd)

cd ..
MAIN_DIR=$(pwd)

IN_DIR=$MAIN_DIR/glm_results/individuals_MIDaffemo
OUT_DIR=$MAIN_DIR/glm_results/group_MIDaffemo_clustersim

N=24


SUBJECTS=$(cat $MAIN_DIR/data/subjects-pupil-fmri.txt)

MASK=$MAIN_DIR/ROIs/glm_mask_gm.nii.gz
RESAMP_DIM=2.9

OUT_FILES_model33b=('pupil_scaled')
STATS_model33b=('26')

REG_FILE=glm33b_b4_mni.nii.gz

cd $OUT_DIR

for i in ${!OUT_FILES_model33b[@]}
do

        out_file=${OUT_FILES_model33b[$i]}
        stat=${STATS_model33b[$i]}
        echo making $out_file
        echo ----- stat: $stat

        # Create the glm list
        glm_sub_list=()
        for sub in $SUBJECTS; do
            glm_sub_list+=("$IN_DIR/${sub}_${REG_FILE}[$stat]")
        done

        # Run ClustSim and ETAC together
        3dttest++ -overwrite -mask $MASK -prefix model33b_${out_file} \
                  -Clustsim \
                  -ETAC \
                  -ETAC_opt NN=2:fpr=0.05:pthr=0.001,0.005,0.01:twosided \
                  -setA "${glm_sub_list[@]}"

        # Process ClustSim output
        3dmerge -doall -1zscore -prefix model33b_${out_file}_${N}z model33b_${out_file}+tlrc
        rm -f model33b_${out_file}+tlrc.*

        3dAFNItoNIFTI -prefix model33b_${out_file}_${N}z.nii.gz model33b_${out_file}_${N}z+tlrc
        rm -f model33b_${out_file}_${N}z+tlrc.*
done
