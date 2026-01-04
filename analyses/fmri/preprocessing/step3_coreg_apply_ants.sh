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
SUBJECTS=$(cat $MAIN_DIR/data/subjects-midaffemo.txt)

ANAT_TEMPLATE_FUNC=$MAIN_DIR/ROIs/mni_ns_func.nii.gz


###############################################################


#for SUBJECT in $SUBJECTS
for SUBJECT in dm240321

do
	echo; echo;
	echo -----------------------------------------------
	echo ------------ processing $SUBJECT --------------
	echo -----------------------------------------------

	cd $SUB_DIR/$SUBJECT/func_proc_MIDaffemo


	for KERNEL in 0 2 4
	do
		### apply the part2mni transformation
		#mv pp_mid_b${KERNEL}_mni.nii.gz pp_mid_b${KERNEL}_mni_ants.nii.gz
		flirt -in pp_MIDaffemo_b${KERNEL}_orig.nii.gz -applyxfm -init ../anat_proc/xfs/func2mni_ants.mat -out pp_MIDaffemo_b${KERNEL}_mni_ants.nii.gz -paddingsize 0.0 -interp trilinear	-ref $ANAT_TEMPLATE_FUNC
	done

done # subjects loop
