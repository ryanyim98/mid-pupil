#!/bin/bash


# written by Leili Mortazavi 6/4/2023
# Adapted by Ryan

#######################################################################################

# usage: download raw data using flywheel CLI
#        this script downloads all scans of a single subject at a time
#        after everything is downloaded, it will organize the data with proper names
#        ---(see 2nd for-loop & define as needed)

#######################################################################################

# run from scripts directory of the project
cd ../
MAIN_DIR=$(pwd)

# using a subject list to download the data from flywheel assumes that you have properly input subject IDs on flywheel.
# SUBJECTS=$(cat $MAIN_DIR/data/subjects.txt)

###################### DOWNLOAD ALL DATA, one subject at a time #######################

# # login to flywheel with Leili's API key
fw login cni.flywheel.io:djEL8O-lgvgk12KlKsKl_TJ1kaFAkEtUVCK-UW3ZE_pvHx1Ef3zT6n5SQ

FW_DIR=knutson/fmrimatch/

# for SUB in $SUBJECTS
 for SUB in ag240330 nq240330 yl240330
do
	echo -----------------------------------------
	echo ------------ downaloading $SUB ----------
	echo -----------------------------------------

	cd $MAIN_DIR/data

	if [ ! -d subjects ]; then
		mkdir subjects
	fi

	cd subjects

	# fw download $FW_DIR/$SUB -i nifti -i bval -i bvec
	fw download $FW_DIR/$SUB -i nifti

	tar -xvf ${SUB}.tar

done # subject loop



###################### Organize the directories as needed ######################

cd $MAIN_DIR/data/subjects

## move all subject folders out of the scitran/knutson/mvpa folder structure
mv scitran/knutson/fmrimatch/* .

for SUB in ag240330 nq240330 yl240330
do
	echo -----------------------------------------
	echo ------------ working on $SUB ------------
	echo -----------------------------------------

	cd $MAIN_DIR/data/subjects/$SUB

	if [ ! -d raw ]; then
		mkdir raw
	fi

	# move the properly named files into the raw folder
	## note: the initial folder name (i.e., CNI session ID) is a 5-digit number starting with 2 in this project

	mv 2*/T1*/*.nii.gz raw/T1.nii.gz
	mv 2*/*MIDaffemo_run1*/*.nii.gz raw/MIDaffemo_run-01.nii.gz
	mv 2*/*MIDaffemo_run2*/*.nii.gz raw/MIDaffemo_run-02.nii.gz
  mv 2*/*MIDaffemo_run3*/*.nii.gz raw/MIDaffemo_run-03.nii.gz
  mv 2*/*MIDaffemo_run4*/*.nii.gz raw/MIDaffemo_run-04.nii.gz
  # mv 2*/movie_run1*/*.nii.gz raw/movie_run-01.nii.gz
	# mv 2*/movie_run2*/*.nii.gz raw/movie_run-02.nii.gz

	# delete the initial 5-digit folder
	#### careful before deleting all the downloaded data ####
	rm -r 2*

done

rm -r $MAIN_DIR/data/subjects/scitran
rm $MAIN_DIR/data/subjects/*.tar

# mv rl240312 el240312
      # #mc240228 has 253 TR for run-01 and 247 TR for run-02. Need to trim/pad it to 250
# cd $MAIN_DIR/data/subjects/mc240228
# ### run 1
# scp MIDaffemo_run-01.nii.gz MIDaffemo_run-01-orig.nii.gz
# 3dTcat -overwrite -prefix MIDaffemo_run-01.nii.gz MIDaffemo_run-01.nii.gz[3..$]
# ### run 2
# scp MIDaffemo_run-02.nii.gz MIDaffemo_run-02-orig.nii.gz
# Create a temporary dataset filled with zeros
# 3dcalc -a MIDaffemo_run-02-orig.nii.gz[0] -expr '0' -datum float -prefix temp_zero.nii.gz
#
# # Concatenate the original dataset with the zero-filled datasets
# 3dTcat -overwrite -prefix MIDaffemo_run-02.nii.gz MIDaffemo_run-02-orig.nii.gz temp_zero.nii.gz temp_zero.nii.gz temp_zero.nii.gz
#
# # Clean up temporary file
# rm temp_zero.nii.gz
