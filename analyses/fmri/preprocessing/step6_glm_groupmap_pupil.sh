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

OUT_FILES_model33b=('pupil_scaled_convolved')
STATS_model33b=('27')

REG_FILE=glm33b_z_b4_mni.nii.gz

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

        3dttest++ -overwrite -mask $MASK -prefix model33b_${out_file} -setA "${glm_sub_list[@]}"

        3dmerge -doall -1zscore -prefix model33b_${out_file}_${N}z model33b_${out_file}+tlrc
        rm -f model33b_${out_file}+tlrc.*

        3dAFNItoNIFTI -prefix model33b_${out_file}_${N}z.nii.gz model33b_${out_file}_${N}z+tlrc
        rm -f model33b_${out_file}_${N}z+tlrc.*

done

# copy the MNI template into the output directory
cp $MAIN_DIR/ROIs/mni_ns.nii.gz $OUT_DIR/.
