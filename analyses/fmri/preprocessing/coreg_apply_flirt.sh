#!/bin/bash


# assumes directory structure:

# MAIN_DIR ----
# 		   ---- scripts
# 		   ---- qa
# 		   ---- data
# 					---- subjects
# 							--- SUB
# 									---- raw
# 									---- func_proc_mid
# 									---- func_proc_risk


# run from the scripts directory
cd ..
MAIN_DIR=$(pwd)

cd data/subjects
SUB_DIR=$(pwd)

# SUBJECTS=$(cat $MAIN_DIR/data/subjects.txt)

ANAT_TEMPLATE_FUNC=$MAIN_DIR/ROIs/mni_ns_func.nii.gz


###############################################################


#for SUBJECT in $SUBJECTS
for SUBJECT in ag240330 nq240330 yl240330
do
	echo; echo;
	echo -----------------------------------------------
	echo ------------ processing $SUBJECT --------------
	echo -----------------------------------------------

	cd $SUB_DIR/$SUBJECT/func_proc_MIDaffemo

	echo
	echo ------------ MIDaffemo --------------
	echo

	for KERNEL in 0 2 4;
	do
	echo -----------------------------------------------
	echo ------------ processing kernel $KERNEL --------------
	echo -----------------------------------------------
		### apply the func2mni transformation
		flirt -in pp_MIDaffemo_b${KERNEL}_orig.nii.gz -applyxfm -init xfs/func2mni.mat -out pp_MIDaffemo_b${KERNEL}_mni.nii.gz -paddingsize 0.0 -interp trilinear	-ref $ANAT_TEMPLATE_FUNC

		### get rid of the ventricles by applying the mask
		ventricle_mask_func.nii.gz
	done

done # subjects loop
