#!/bin/bash

#run from the scripts directory
cd ..
MAIN_DIR=$(pwd)
QA_DIR=$MAIN_DIR/qa/qa_ants
cd data/subjects
SUB_DIR=$(pwd)

SUBJECTS=$(cat $MAIN_DIR/data/subjects-midaffemo.txt)

#for SUBJECT in $SUBJECTS
for SUBJECT in sj240311
do
    echo; echo;
    echo -----------------------------------------------
    echo ------------ processing $SUBJECT --------------
    echo -----------------------------------------------

    cd $SUB_DIR/$SUBJECT/anat_proc

    if [ ! -d xfs ]; then
        mkdir xfs
    fi

    3dresample -overwrite -input t1_ns.nii.gz -dxyz 2.9 2.9 2.9 -prefix t1_ns_2.9.nii.gz
    3dresample -overwrite -input t1_mni.nii.gz -dxyz 2.9 2.9 2.9 -prefix t1_mni_2.9.nii.gz

    c3d_affine_tool -ref t1_mni_2.9.nii.gz -src t1_ns_2.9.nii.gz -itk t12mni_xform_Affine.txt -ras2fsl -o xfs/t12mni_ants.mat
    #c3d_affine_tool -itk t12mni_xform_Affine.txt -o fsl_transformation_file2.mat

    ## test the matrix transformation
    # flirt -in t1_ns.nii.gz -applyxfm -init xfs/t12mni_ants.mat \
    #   -out t1_mni_test.nii.gz -paddingsize 0.0 -interp trilinear \
    #   -ref t1_mni.nii.gz


  ################### MID ###########################

    ### coregister ref_vol_whole to T1
    flirt -ref t1_ns_2.9.nii -in ../func_proc_MIDaffemo/ref_vol_ns.nii.gz -out ref_vol_func2t1 -omat xfs/func2t1.mat\
            -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6  -interp trilinear

    ### concatenate the xfs for whole2t1 and t12mni to get whole2mni
    convert_xfm -concat xfs/t12mni_ants.mat xfs/func2t1.mat -omat xfs/func2mni_ants.mat

    ### apply the part2t1 transformation
    flirt -in ../func_proc_MIDaffemo/ref_vol_ns.nii.gz -applyxfm -init xfs/func2mni_ants.mat \
        -out ref_vol_mni_ants.nii.gz -paddingsize 0.0 -interp trilinear \
        -ref t1_mni_2.9.nii.gz

    ## coreg in another way
    # antsApplyTransforms -d 3 -i anat_proc/t1_acpc_ns.nii.gz -r anat_proc/t1_mni.nii.gz -o anat_proc/t1_mni_test.nii.gz -t anat_proc/t12mni_xform_Affine.txt

    ### copy the outputs t1_mni to qa folder
    cp t1_mni.nii.gz $QA_DIR/${SUBJECT}_t1_mni.nii.gz
    cp ref_vol_mni_ants.nii.gz $QA_DIR/${SUBJECT}_ref_vol_mni_ants.nii.gz


done # subject loop for mid1st
