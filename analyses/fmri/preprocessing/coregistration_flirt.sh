#!/bin/bash

#run from the scripts directory
cd ..
MAIN_DIR=$(pwd)

cd data/subjects
SUB_DIR=$(pwd)

if [ ! -d $MAIN_DIR/qa ]; then
    mkdir $MAIN_DIR/qa
fi

QA_DIR=$MAIN_DIR/qa/coreg

if [ ! -d $QA_DIR ]; then
    mkdir $QA_DIR
fi

ROI_DIR=$MAIN_DIR/ROIs

# SUBJECTS=$(cat $MAIN_DIR/data/subjects.txt)

3dresample -overwrite -input $ROI_DIR/mni_ns.nii.gz -dxyz 2.9 2.9 2.9 -prefix $ROI_DIR/mni_ns_func.nii.gz
3dresample -overwrite -input $ROI_DIR/ventricle_mask.nii.gz -dxyz 2.9 2.9 2.9 -prefix $ROI_DIR/ventricle_mask_func.nii.gz

ANAT_TEMPLATE=$MAIN_DIR/ROIs/mni_ns.nii.gz
ANAT_TEMPLATE_FUNC=$MAIN_DIR/ROIs/mni_ns_func.nii.gz

#SUBJECTS =$(cat $MAIN_DIR/data/subjects.txt)

#for SUBJECT in $SUBJECTS_mid1st
for SUBJECT in ag240330 nq240330 yl240330
do
    echo; echo;
    echo -----------------------------------------------
    echo ------------ processing $SUBJECT --------------
    echo -----------------------------------------------

    ################### MID ###########################
    if [ ! -d $SUB_DIR/$SUBJECT/func_proc_MIDaffemo ]; then
        mkdir $SUB_DIR/$SUBJECT/func_proc_MIDaffemo
    fi

    cd $SUB_DIR/$SUBJECT/func_proc_MIDaffemo

    if [ ! -d xfs ]; then
        mkdir xfs
    fi

    ### skull-stip T1
    3dSkullStrip -overwrite -prefix t1_ns.nii.gz -input ../raw/T1.nii.gz -push_to_edge

    ### take a reference volume & skull strip it
    3dTcat -overwrite -prefix ref_vol.nii.gz ../raw/MIDaffemo_run-01.nii.gz[0]
    3dSkullStrip -overwrite -input ref_vol.nii.gz -prefix ref_vol_ns.nii.gz

    ### # register t1 to mni:
    flirt -ref $ANAT_TEMPLATE -in t1_ns.nii.gz -out t1_mni -omat xfs/t12mni.mat \
        -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear

    ### coregister ref_vol to mni
    flirt -in ref_vol_ns.nii.gz -ref $ANAT_TEMPLATE_FUNC -out ref_vol_mni \
      -omat xfs/func2mni.mat -bins 256 -cost corratio \
      -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear

    convert_xfm -inverse xfs/func2mni.mat -omat xfs/mni2func.mat

    ### copy the outputs t1_mni to qa folder
    cp t1_mni.nii.gz $QA_DIR/${SUBJECT}_t1_mni.nii.gz
    cp ref_vol_mni.nii.gz $QA_DIR/${SUBJECT}_ref_vol_mni.nii.gz

done # subject loop
