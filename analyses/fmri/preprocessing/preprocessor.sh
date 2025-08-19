#!/bin/bash


# assumes directory structure:

# MAIN_DIR ----
# 		   ---- scripts
# 		   ---- qa
# 		   ---- data
# 					---- subjects
# 							--- SUB
# 									---- raw
# 									---- func_proc_MIDaffemo
# 									---- func_proc_risk


# run from the scripts directory
cd ..
MAIN_DIR=$(pwd)

cd data/subjects
SUB_DIR=$(pwd)

#SUBJECTS=$(cat $MAIN_DIR/data/subjects.txt)

### qa directory
if [ ! -d  $MAIN_DIR/qa/qa_motion_MIDaffemo ]; then
    mkdir $MAIN_DIR/qa/qa_motion_MIDaffemo
fi
QA_DIR=$MAIN_DIR/qa/qa_motion_MIDaffemo


###############################################################


#for SUBJECT in $SUBJECTS
for SUBJECT in ag240330 nq240330 yl240330
do
	echo; echo;
	echo -----------------------------------------------
	echo ------------ processing $SUBJECT --------------
	echo -----------------------------------------------

	cd $SUB_DIR/$SUBJECT
  # if [ "$SUBJECT" = "lm240209" ]; then
  #   RUNS=("run-01" "run-02" "run-03" "run-04")
  # else
  #   RUNS=("run-01" "run-02")
  # fi

  RUNS=("run-01" "run-02" "run-03" "run-04")

	### subject input & output directories
    inDir=$SUB_DIR/$SUBJECT/raw
    outDir=$SUB_DIR/$SUBJECT/func_proc_MIDaffemo

	### make outDir & cd to it:
    if [ ! -d "$outDir" ]; then
        mkdir $outDir
    fi

    cd $outDir

    ### also make a "xfs" directory to house all xform files
    if [ ! -d xfs ]; then
        mkdir xfs
    fi


  for RUN in ${RUNS[@]}
	do
	    echo
	    echo
	    echo ------WOKING ON $RUN------
	    echo

	    ### trim the first 5 volumes
	    3dTcat -overwrite -prefix MIDaffemo_${RUN}_trim.nii.gz $inDir/MIDaffemo_${RUN}.nii.gz[5..$]

	    ### time slice correction
	    3dTshift -overwrite -slice 0 -tpattern altplus -prefix MIDaffemo_${RUN}_ts.nii.gz MIDaffemo_${RUN}_trim.nii.gz

	    ### motion correction
	    3dvolreg -overwrite -Fourier -twopass -prefix MIDaffemo_${RUN}_m.nii.gz -base ref_vol.nii.gz -dfile 3dmotion_${RUN}.1D MIDaffemo_${RUN}_ts.nii.gz
        1d_tool.py -infile 3dmotion_${RUN}.1D[1..6] -show_censor_count -censor_prev_TR -censor_motion 0.5 MIDaffemo_${RUN}_trim.nii.gz
		  mv MIDaffemo_${RUN}_trim.nii.gz_enorm.1D MIDaffemo_${RUN}_enorm.1D
    	mv MIDaffemo_${RUN}_trim.nii.gz_censor.1D motion_censor_${RUN}.1D
    	mv MIDaffemo_${RUN}_trim.nii.gz_CENSORTR.txt MIDaffemo_${RUN}_CENSORTR.txt

    	### apply different smoothing kernels

    	for KERNEL in 0 2 4
    	do
		    ### spatial smoothing
		    3dmerge -overwrite -prefix MIDaffemo_${RUN}_b${KERNEL}.nii.gz -1blur_fwhm $KERNEL -doall MIDaffemo_${RUN}_m.nii.gz

		    ### convert to percent signal change
	    	3dTstat -overwrite -prefix MIDaffemo_${RUN}_avg_b${KERNEL}.nii.gz MIDaffemo_${RUN}_b${KERNEL}.nii.gz[0..$]
	    	3drefit -overwrite -abuc MIDaffemo_${RUN}_avg_b${KERNEL}.nii.gz
	    	3dcalc -overwrite -datum float -a MIDaffemo_${RUN}_b${KERNEL}.nii.gz[0..$] -b MIDaffemo_${RUN}_avg_b${KERNEL}.nii.gz -expr "((a-b)/b)*100" -prefix MIDaffemo_${RUN}_psc_b${KERNEL}.nii.gz

	    	### high-pass filtering
	    	# 3dFourier -prefix pp_MIDaffemo_${RUN}_b${KERNEL}.nii.gz -highpass .011 MIDaffemo_${RUN}_psc_b${KERNEL}.nii.gz
	    	3dBandpass -overwrite -prefix pp_MIDaffemo_${RUN}_b${KERNEL}.nii.gz .011 9999 MIDaffemo_${RUN}_psc_b${KERNEL}.nii.gz

	    done # kernel loop

    done # run loop

    ### concatenate runs

    for KERNEL in 0 2 4
    do
      3dTcat -overwrite -prefix pp_MIDaffemo_b${KERNEL}_orig.nii.gz pp_MIDaffemo_run-01_b${KERNEL}.nii.gz pp_MIDaffemo_run-02_b${KERNEL}.nii.gz pp_MIDaffemo_run-03_b${KERNEL}.nii.gz pp_MIDaffemo_run-04_b${KERNEL}.nii.gz
      # if [ "$SUBJECT" = "lr240201" ]; then #Luis had 2 runs
      #   3dTcat -overwrite -prefix pp_MIDaffemo_b${KERNEL}_orig.nii.gz pp_MIDaffemo_run-01_b${KERNEL}.nii.gz pp_MIDaffemo_run-02_b${KERNEL}.nii.gz
      # else
      #   if [ "$SUBJECT" = "cw240110" ]; then #Cynthia had 1 run
      #     3dTcat -overwrite -prefix pp_MIDaffemo_b${KERNEL}_orig.nii.gz pp_MIDaffemo_run-01_b${KERNEL}.nii.gz
      # fi
	done # kernel loop

	# cat motion_censor_run-01.1D motion_censor_run-02.1D >> motion_censor.1D
	# cat 3dmotion_run-01.1D 3dmotion_run-02.1D >> 3dmotion.1D
	# 1dplot -png $QA_DIR/${SUBJECT}_motionplot 3dmotion.1D[4..6]
  rm motion_censor.1D
  rm 3dmotion.1D

  cat motion_censor_run-01.1D motion_censor_run-02.1D motion_censor_run-03.1D motion_censor_run-04.1D >> motion_censor.1D
  cat 3dmotion_run-01.1D 3dmotion_run-02.1D 3dmotion_run-03.1D 3dmotion_run-04.1D >> 3dmotion.1D
  1dplot -png $QA_DIR/${SUBJECT}_motionplot 3dmotion.1D[4..6]

  # if [ "$SUBJECT" = "cw240110" ]; then
  #   cat motion_censor_run-01.1D >> motion_censor.1D
  #   cat 3dmotion_run-01.1D >> 3dmotion.1D
  #   1dplot -png $QA_DIR/${SUBJECT}_motionplot 3dmotion_run-01.1D[4..6]
  # else
  #   if [ "$SUBJECT" = "lr240201" ]; then
  #     cat motion_censor_run-01.1D motion_censor_run-02.1D >> motion_censor.1D
  #     cat 3dmotion_run-01.1D 3dmotion_run-02.1D >> 3dmotion.1D
  #     1dplot -png $QA_DIR/${SUBJECT}_motionplot 3dmotion.1D[4..6]
  # fi
done # subject loop
