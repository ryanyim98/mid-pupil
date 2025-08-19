#!/usr/bin/python

# filename: t12mni_script.py
# script to do the following:

# skull-strip t1 in acpc space
# estimate xform to mni space & save out xforms
# apply xform


import os,sys,glob


##################### define global variables #################################
# EDIT AS NEEDED:

main_dir = os.path.join(os.path.expanduser('~'),'dopa', 'midaffemo')

data_dir = os.path.join(os.path.expanduser('~'),'dopa', 'midaffemo','data', 'subjects')

t1_dir = os.path.join(data_dir,'%s','anat_proc')  # first %s is data_dir & 2nd is subject id

mni_file = os.path.join(main_dir,'ROIs','mni_ns.nii.gz') # %s is data_dir

print('execute commands?')
xc = bool(input('enter 1 for yes, or 0 to only print: '))


## add ants directory to path
# os.system('export PATH=$PATH:'+os.path.join(os.path.expanduser('~')+'/repos/antsbin/bin'))


###############################################################################
############################### DO IT #########################################
###############################################################################


#########  print commands & execute if xc is True, otherwise just print them
def doCommand(cmd):

	print(cmd+'\n')
	if xc is True:
		os.system(cmd)


#########  get main data directory and subjects to process
# def whichSubs():


# 	from getCueSubjects import getsubs
# 	subjects,gi = getsubs()

# 	print ' '.join(subjects)

# 	input_subs = raw_input('subject id(s) (hit enter to process all subs): ')
# 	print '\nyou entered: '+input_subs+'\n'

# 	if input_subs:
# 		subjects=input_subs.split(' ')

# 	return subjects




if __name__ == '__main__':

	# subjects = whichSubs()

	# Open the file and read lines
	subjects_list = os.path.join(main_dir, 'data', 'subjects-midaffemo.txt')

	with open(subjects_list, 'r') as file:
		# Read all lines into a list and strip the newline characters
		subjects = [line.strip() for line in file.readlines()]

	subjects=['ab240324','ag240330','ak240314','am240301','an240314','as240303','av240312','ax240308','by240324','ch240307','dm240321','ek240305',
'el240312','ey240320','gy240322','hc240324','ig240319']

	for subject in subjects:

		print('WORKING ON SUBJECT '+subject+'\n')


		# define subject's raw & pre-processed directories
		t1_fs = os.path.join(t1_dir,'t1_fs.nii.gz') % (subject)# acpc aligned freesurfer t1
		t1_ns = os.path.join(t1_dir,'t1_ns.nii.gz') % (subject) # name for skullstripped t1
		t1_mni =  os.path.join(t1_dir,'t1_mni.nii.gz') % (subject) # name for t1 in mni space
		xform = os.path.join(t1_dir,'t12mni_xform_') % (subject)  # name for t1 > mni xform

		####### skull strip t1
		if os.path.isfile(t1_ns):
			print('\n skull stripped t1 file '+t1_ns+' already exists...\n')
		else:
			cmd = '3dSkullStrip -prefix '+t1_ns+' -input '+t1_fs
			doCommand(cmd)


		######## estimate xform from t1 native space to mni space
		cmd = 'ANTS 3 -m CC['+mni_file+','+t1_ns+',1,4] -r Gauss[3,0] -o '+xform+' -i 100x50x30x10 -t SyN[.25]'
		doCommand(cmd)


		########## apply xform on t1
		cmd = 'WarpImageMultiTransform 3 '+t1_ns+' '+t1_mni+' '+xform+'Warp.nii.gz '+xform+'Affine.txt'
		doCommand(cmd)


		# change header to play nice with afni
		cmd = '3drefit -view tlrc -space MNI '+t1_mni
		doCommand(cmd)


		print('FINISHED SUBJECT '+subject)
