#!/bin/bash


# modified by Ryan from Leili's script
#
# Group-level inference uses 3dMEMA (mixed-effects meta-analysis; Chen et al., 2012)
# so within-subject noise (per-voxel beta SE via t) is combined with between-subject variance.
#
# Input per subject: same MNI bucket as s7_glm27c_apply.sh (glm27c_b4_mni.nii.gz).
# STATS_model27c sub-brick = beta for param_arous_early; paired t is the next sub-brick
# (3dDeconvolve -tout order: Full_F, R^2, then [beta, tstat] pairs per regressor).

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

OUT_FILES_model27c=('param_arous_early')
STATS_model27c=('26')
# t-stat sub-brick paired with each beta under 3dDeconvolve -tout (beta index + 1)
# if wm and csf are included as non-nuisance regressors
# STATS_model27c=('11' '15' '19' '23' '27' '31' '35')

REG_FILE=glm27c_b4_mni.nii.gz

cd $OUT_DIR
export AFNI_DECONFLICT=OVERWRITE

for i in ${!OUT_FILES_model27c[@]}
do

        out_file=${OUT_FILES_model27c[$i]}
        beta_brick=${STATS_model27c[$i]}
        t_brick=$((beta_brick + 1))
        echo making $out_file
        echo ----- beta sub-brick: $beta_brick, t sub-brick: $t_brick

        mema_argv=()
        for sub in $SUBJECTS; do
            mema_argv+=("$sub" "$IN_DIR/${sub}_${REG_FILE}[${beta_brick}]" "$IN_DIR/${sub}_${REG_FILE}[${t_brick}]")
        done

        # Clean previous AFNI outputs for this contrast
        rm -f model27c_${out_file}+tlrc.* model27c_${out_file}_${N}z+tlrc.* model27c_${out_file}_${N}z.nii.gz

        # MEMA sub-brick names depend on AFNI build; common forms are g:b / g:t (your 3dinfo -label)
        # or g#0_Coef / g#0_Tstat. Use 3dinfo -label2index and pass a numeric selector to 3dmerge.
        3dMEMA -overwrite -prefix model27c_${out_file} \
                -jobs 4 \
                -mask $MASK \
                -missing_data 0 \
                -HKtest \
                -model_outliers \
                -residual_Z \
                -set g "${mema_argv[@]}"

        # AFNI_21.x: 3dmerge does not accept some label selectors (e.g. with '#'); use brick index.
        mema_dset="model27c_${out_file}+tlrc"
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
                echo "Sub-brick labels from: 3dinfo -label $mema_dset" >&2
                3dinfo -label "$mema_dset" >&2
                exit 1
        fi

        # Match prior pipeline: voxelwise z-score of the group map (for display / legacy scripts)
        3dmerge -overwrite -doall -1zscore -prefix model27c_${out_file}_${N}z "${mema_dset}[${tstat_i}]"
        rm -f model27c_${out_file}+tlrc.*

        3dAFNItoNIFTI -overwrite -prefix model27c_${out_file}_${N}z.nii.gz model27c_${out_file}_${N}z+tlrc
        rm -f model27c_${out_file}_${N}z+tlrc.*

done
