#!/bin/bash

# converted to bash by Leili Mortazavi, August 17, 2024
# based on Matlab script by Kelly MacNiven

# usage: this script takes segmentation files made using freesurfer in
# .mgz format and converts them to nifti files. Also changes some of the
# segmentation indices for better functionality with mr vista software.

# the script takes in files created from Freesurfer's "recon-all" call, and
# converts them to nifti files.


# t1.class.nii.gz - segmentation file with the following indices:
# white matter L: 3
# white matter R: 4
# gray matter L: 5
# gray matter R: 6
# unlabeled: either 1 or 0

# aparc_aseg - segmentation of cortex and sub-cortical structures.
# aparc.a2009s+aseg - Destrieux atlas segmentation of cortex and
# sub-cortical structures (for vlpfc and anterior insula VOIs)
# See here for look up table of segmentation indices:
#     https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT


# Important Note: make sure $FREESURFER_HOME path exists in .bash_profile or .bashrc

scriptsDir=$(pwd)
cd ..
mainDir=$(pwd)
dataDir=$mainDir/data

SUBJECTS='ab240324 ag240330 ak240314 am240301 an240314 as240303 av240312 ax240308 by240324 ch240307 dm240321 ek240305 el240312 ey240320 gy240322 hc240324 ig240319'
# SUBJECTS=$(cat $mainDir/data/subjects_dti2mm_all.txt)

resample_type=weighted

for SUB in $SUBJECTS
do

	echo
	echo -------------
	echo -------------working on subject: $SUB ------------
	echo -------------
	echo

	cd $dataDir/subjects/$SUB

	inDir=$FREESURFER_HOME/subjects/$SUB/mri

	if [ ! -d "${inDir}" ]; then
		echo SKIPPING subject $SUB!!!
		echo because NO FreeSurfer folder exists for subject $SUB
	else
		mri_convert --out_orientation RAS -i $inDir/T1.mgz -o anat_proc/t1_fs.nii.gz
		mri_convert --out_orientation RAS -i $inDir/ribbon.mgz -o anat_proc/t1_class.nii.gz -rt $resample_type
		mri_convert --out_orientation RAS -i $inDir/aparc+aseg.mgz -o anat_proc/aparc+aseg.nii.gz -rt $resample_type
		mri_convert --out_orientation RAS -i $inDir/aparc.a2009s+aseg.mgz -o anat_proc/aparc.a2009s+aseg.nii.gz -rt $resample_type


		# # Convert freesurfer label values to itkGray label values
        # We want to convert
        #  - Left white:   2 => 3
        #  - Left gray:    3 => 5
        #  - Right white: 41 => 4
        #  - Right gray:  42 => 6

        3dcalc -a anat_proc/t1_class.nii.gz -expr 'equals(a, 2)*3' -prefix anat_proc/t1_class_wmL.nii.gz
        3dcalc -a anat_proc/t1_class.nii.gz -expr 'equals(a, 3)*5' -prefix anat_proc/t1_class_gmL.nii.gz
        3dcalc -a anat_proc/t1_class.nii.gz -expr 'equals(a, 41)*4' -prefix anat_proc/t1_class_wmR.nii.gz
        3dcalc -a anat_proc/t1_class.nii.gz -expr 'equals(a, 42)*6' -prefix anat_proc/t1_class_gmR.nii.gz
        3dcalc -overwrite \
        	   -a anat_proc/t1_class_wmL.nii.gz \
        	   -b anat_proc/t1_class_wmR.nii.gz \
        	   -c anat_proc/t1_class_gmL.nii.gz \
        	   -d anat_proc/t1_class_gmR.nii.gz \
        	   -expr 'a+b+c+d' -prefix anat_proc/t1_class_labeled.nii.gz

        # # make a white-matter mask (cerebral only)
        3dcalc -a anat_proc/t1_class_wmL.nii.gz -b anat_proc/t1_class_wmR.nii.gz -expr 'step(a+b)' -prefix anat_proc/wm_mask.nii.gz

        # make a more comprehensive white-matter + striatal regions
        ## FreeSurfer regions LEFT: 16 2 10 11 12 13 26 28 11107 11113 11154 11116 11117 11118 11124 11139 11148 11149 11155 11163 11164 11165
        3dcalc -a anat_proc/aparc.a2009s+aseg.nii.gz -expr 'or(equals(a, 16), equals(a, 2), equals(a, 10), equals(a, 11), equals(a, 12),equals(a, 13), equals(a, 26), equals(a, 28), equals(a, 11113), equals(a, 11154), equals(a, 11107), equals(a, 11116), equals(a, 11117), equals(a, 11118), equals(a, 11124), equals(a, 11139), equals(a, 11148), equals(a, 11149), equals(a, 11155), equals(a, 11163), equals(a, 11164), equals(a, 11165))' \
        	   -prefix anat_proc/wm_mask_left.nii.gz
        ## FreeSurfer regions RIGHT: 16 41 49 50 51 52 58 60 12106 12113 12115 12116 12117 12118 12124 12139 12148 12149 12155 12163 12164 12165
        3dcalc -a anat_proc/aparc.a2009s+aseg.nii.gz -expr 'or(equals(a, 16), equals(a, 41), equals(a, 49), equals(a, 50), equals(a, 51),equals(a, 52), equals(a, 58), equals(a, 60), equals(a, 12106), equals(a, 12113), equals(a, 12115), equals(a, 12116), equals(a, 12117), equals(a, 12118), equals(a, 12124), equals(a, 12139), equals(a, 12148), equals(a, 12149), equals(a, 12155), equals(a, 12163), equals(a, 12164), equals(a, 12165))' \
        	   -prefix anat_proc/wm_mask_right.nii.gz
	fi

done
