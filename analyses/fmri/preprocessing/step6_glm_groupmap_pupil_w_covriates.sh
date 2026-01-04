#!/bin/bash


# modified by Ryan from Leili's script

# run from the scripts directory

SCRIPTS_DIR=$(pwd)

cd ..
MAIN_DIR=$(pwd)

IN_DIR=$MAIN_DIR/glm_results/individuals_MIDaffemo
OUT_DIR=$MAIN_DIR/glm_results/group_MIDaffemo


SUBJECTS=$(cat $MAIN_DIR/data/subjects-pupil-fmri.txt)

N=24

MASK=$MAIN_DIR/ROIs/glm_mask_gm.nii.gz
RESAMP_DIM=2.9

OUT_FILES_model33c=('pupil_scaled_convolved' 'luminance' 'pupilx' 'blink' 'saccade')
STATS_model33c=('27' '31' '35' '39' '43')

REG_FILE=glm33c_z_b4_mni.nii.gz

cd $OUT_DIR

for i in ${!OUT_FILES_model33c[@]}
do

        out_file=${OUT_FILES_model33c[$i]}
        stat=${STATS_model33c[$i]}
        echo making $out_file
        echo ----- stat: $stat

        # Create the glm list
        glm_sub_list=()
        for sub in $SUBJECTS; do
            glm_sub_list+=("$IN_DIR/${sub}_${REG_FILE}[$stat]")
        done

        3dttest++ -overwrite -mask $MASK -prefix model33c_${out_file} -setA "${glm_sub_list[@]}"

        3dmerge -doall -1zscore -prefix model33c_${out_file}_${N}z model33c_${out_file}+tlrc
        rm -f model33c_${out_file}+tlrc.*

        3dAFNItoNIFTI -prefix model33c_${out_file}_${N}z.nii.gz model33c_${out_file}_${N}z+tlrc
        rm -f model33c_${out_file}_${N}z+tlrc.*

done

# copy the MNI template into the output directory
cp $MAIN_DIR/ROIs/mni_ns.nii.gz $OUT_DIR/.
