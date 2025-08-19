#!/bin/bash
#glm27 is the model with parametric regressors for arousal ratings REGARDLESS of probe type and probe timing (ant-ant, out-out)
# (neural activity 1 TR prior predicting subsequent affect)

SCRIPTS_DIR=$(pwd)

cd ..

MAIN_DIR=$(pwd)

cd data
DATA_DIR=$(pwd)

SUBJECTS=$MAIN_DIR/data/subjects-midaffemo.txt

OUT_DIR=$MAIN_DIR/glm_results/individuals_MIDaffemo

ANAT_TEMPLATE_FUNC=$MAIN_DIR/ROIs/mni_ns_func.nii.gz

#mc240228 cannot get glm because of missing values
#for SUBJECT in sp240228
for SUBJECT in $(cat $SUBJECTS)
do
	echo
	echo
	echo processing $SUBJECT
	cd $DATA_DIR/subjects/$SUBJECT/func_proc_MIDaffemo


		for KERNEL in 0 2 4
		do

		    3dDeconvolve -overwrite -float -input pp_MIDaffemo_b${KERNEL}_orig.nii.gz -nfirst 0 -num_stimts 15 -xjpeg Xmat -polort 2 -concat '1D: 0 245 490 735' \
						-censor motion_censor.1D -mask bmask.nii.gz \
						-stim_file 1 3dmotion.1D'[1]' -stim_base 1 -stim_label 1 roll \
		        -stim_file 2 3dmotion.1D'[2]' -stim_base 2 -stim_label 2 pitch \
		        -stim_file 3 3dmotion.1D'[3]' -stim_base 3 -stim_label 3 yaw \
		        -stim_file 4 3dmotion.1D'[4]' -stim_base 4 -stim_label 4 dS \
		        -stim_file 5 3dmotion.1D'[5]' -stim_base 5 -stim_label 5 dL \
		        -stim_file 6 3dmotion.1D'[6]' -stim_base 6 -stim_label 6 dP \
		        -stim_file 7 roi_ts/MIDaffemo_b${KERNEL}_wm_ants.1D -stim_base 7 -stim_label 7 white_matter \
		        -stim_file 8 roi_ts/MIDaffemo_b${KERNEL}_csf_ants.1D -stim_base 8 -stim_label 8 csf \
						-stim_file 9 regs/ant1st4s_c.1D -stim_label 9 ant4s \
						-stim_file 10 regs/ant2nd2s_c.1D -stim_label 10 ant2nd \
						-stim_file 11 regs/out_c.1D -stim_label 11 out \
						-stim_file 12 regs/button_c.1D -stim_label 12 button \
						-stim_file 13 regs/reportaff_rtmod_c.1D -stim_label 13 reportaff_rtmod \
						-stim_file 14 regs/reportemo_rtmod_c.1D -stim_label 14 reportemo_rtmod \
						-stim_file 15 regs/param_arous_c.1D -stim_label 15 param_arous \
		        -tout -fout -rout -bucket $OUT_DIR/${SUBJECT}_glm27_b${KERNEL}

		    3dmerge -doall -1zscore -prefix $OUT_DIR/${SUBJECT}_glm27_z_b${KERNEL} $OUT_DIR/${SUBJECT}_glm27_b${KERNEL}+orig

				3dAFNItoNIFTI -overwrite -prefix $OUT_DIR/${SUBJECT}_glm27_b${KERNEL}_orig.nii.gz $OUT_DIR/${SUBJECT}_glm27_b${KERNEL}+orig
				3dAFNItoNIFTI -overwrite -prefix $OUT_DIR/${SUBJECT}_glm27_z_b${KERNEL}_orig.nii.gz $OUT_DIR/${SUBJECT}_glm27_z_b${KERNEL}+orig

				#replace the z-scored R2 with actual R2
				temp_vol="temp_vol.nii.gz"
				temp_output="temp_output.nii.gz"
				# Extract the first volume from the input dataset
				3dTcat -overwrite -prefix $temp_vol $OUT_DIR/${SUBJECT}_glm27_b${KERNEL}_orig.nii.gz[0]
				# Replace the first volume of the output dataset with the extracted volume
				3dTcat -overwrite -prefix $temp_output $temp_vol $OUT_DIR/${SUBJECT}_glm27_z_b${KERNEL}_orig.nii.gz[1..$]
				# do some cleaning
				mv $temp_output $OUT_DIR/${SUBJECT}_glm27_z_b${KERNEL}_orig.nii.gz
				rm $temp_vol

				rm $OUT_DIR/${SUBJECT}_glm27_b${KERNEL}+orig.*
				rm $OUT_DIR/${SUBJECT}_glm27_z_b${KERNEL}+orig.*

				# for some reason when we aply flirt to the non-indexed data, it only transforms the 1st volume, but indexing solves the issue (!)
				3dTcat -overwrite -prefix $OUT_DIR/${SUBJECT}_glm27_z_b${KERNEL}_IND_orig.nii.gz $OUT_DIR/${SUBJECT}_glm27_z_b${KERNEL}_orig.nii.gz[0..$]

				### apply the part2mni transformation
				flirt -in $OUT_DIR/${SUBJECT}_glm27_z_b${KERNEL}_IND_orig.nii.gz -init ../anat_proc/xfs/func2mni_ants.mat -applyxfm \
						-out $OUT_DIR/${SUBJECT}_glm27_z_b${KERNEL}_mni.nii.gz -paddingsize 0.0 -interp trilinear \
						-ref $ANAT_TEMPLATE_FUNC

			done # kernel loop
done # subject loop

# copy the MNI template into the output directory
cp $MAIN_DIR/ROIs/mni_ns.nii.gz $OUT_DIR/.
