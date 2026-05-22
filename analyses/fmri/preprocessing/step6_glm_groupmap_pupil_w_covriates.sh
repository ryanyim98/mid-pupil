#!/bin/bash


# modified by Ryan from Leili's script
#
# Group-level inference uses 3dMEMA (mixed-effects meta-analysis; Chen et al., 2012)
# so within-subject noise (per-voxel beta SE via t) is combined with between-subject variance.
#
# Input: per-subject MNI buckets from s7_glm33c_apply (glm33c_b4_mni.nii.gz).
# STATS_model33c = beta sub-brick for each contrast; paired t from 3dDeconvolve -tout is beta + 1.

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
STATS_model33c=('26' '30' '34' '38' '42')

REG_FILE=glm33c_b4_mni.nii.gz

cd $OUT_DIR
export AFNI_DECONFLICT=OVERWRITE

for i in ${!OUT_FILES_model33c[@]}
do

        out_file=${OUT_FILES_model33c[$i]}
        beta_brick=${STATS_model33c[$i]}
        t_brick=$((beta_brick + 1))
        echo making $out_file
        echo ----- beta sub-brick: $beta_brick, t sub-brick: $t_brick

        mema_argv=()
        for sub in $SUBJECTS; do
            mema_argv+=("$sub" "$IN_DIR/${sub}_${REG_FILE}[${beta_brick}]" "$IN_DIR/${sub}_${REG_FILE}[${t_brick}]")
        done

        rm -f model33c_${out_file}+tlrc.* model33c_${out_file}_${N}z+tlrc.* model33c_${out_file}_${N}z.nii.gz

        3dMEMA -overwrite -prefix model33c_${out_file} \
                -jobs 4 \
                -mask $MASK \
                -missing_data 0 \
                -HKtest \
                -model_outliers \
                -residual_Z \
                -set g "${mema_argv[@]}"

        mema_dset="model33c_${out_file}+tlrc"
        tstat_i=""
        for cand in 'g:t' 'g#0_Tstat' 'G1#0_Tstat' 'c1#0_Tstat'; do
                tstat_i=$(3dinfo -label2index "$cand" "$mema_dset" 2>/dev/null | awk 'NF {print $1; exit}')
                if [[ "$tstat_i" =~ ^[0-9]+$ ]]; then
                        echo "MEMA group T-stat sub-brick: [$tstat_i] (label $cand)"
                        break
                fi
                tstat_i=""
        done
        if [[ ! "$tstat_i" =~ ^[0-9]+$ ]]; then
                echo "ERROR: could not resolve MEMA group T-stat sub-brick for $mema_dset" >&2
                3dinfo -label "$mema_dset" >&2
                exit 1
        fi

        3dmerge -overwrite -doall -1zscore -prefix model33c_${out_file}_${N}z "${mema_dset}[${tstat_i}]"
        rm -f model33c_${out_file}+tlrc.*

        3dAFNItoNIFTI -overwrite -prefix model33c_${out_file}_${N}z.nii.gz model33c_${out_file}_${N}z+tlrc
        rm -f model33c_${out_file}_${N}z+tlrc.*

done

# copy the MNI template into the output directory
cp $MAIN_DIR/ROIs/mni_ns.nii.gz $OUT_DIR/.
