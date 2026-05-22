#!/bin/bash


#
# Group-level: 3dMEMA (Chen et al., 2012); per-subject beta + t from 3dDeconvolve -tout.
#
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
STATS_model55=('26' '30')
# if wm and csf are included as non-nuisance regressors
# STATS_model55=('11' '15' '19' '23' '27' '31' '35')

REG_FILE=glm55_b4_mni.nii.gz

cd $OUT_DIR

export AFNI_DECONFLICT=OVERWRITE
for i in ${!OUT_FILES_model55[@]}
do

        out_file=${OUT_FILES_model55[$i]}
        stat=${STATS_model55[$i]}
        echo making $out_file
        echo ----- stat: $stat

        # Create the glm list
        mema_argv=()
        for sub in $SUBJECTS; do
            mema_argv+=("$sub" "$IN_DIR/${sub}_${REG_FILE}[${stat}]" "$IN_DIR/${sub}_${REG_FILE}[$((stat + 1))]")
        done
        3dMEMA -overwrite -prefix model55_${out_file} \
                -jobs 4 \
                -mask $MASK \
                -missing_data 0 \
                -HKtest \
                -model_outliers \
                -residual_Z \
                -set g "${mema_argv[@]}"
        mema_dset="model55_${out_file}+tlrc"
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
        3dmerge -overwrite -doall -1zscore -prefix model55_${out_file}_${N}z "${mema_dset}[${tstat_i}]"

        rm -f model55_${out_file}+tlrc.*

        3dAFNItoNIFTI -overwrite -prefix model55_${out_file}_${N}z.nii.gz model55_${out_file}_${N}z+tlrc
        rm -f model55_${out_file}_${N}z+tlrc.*

done
s
