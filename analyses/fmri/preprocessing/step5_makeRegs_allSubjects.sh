#!/bin/bash


SCRIPTS_DIR=$(pwd)

cd ..

MAIN_DIR=$(pwd)

DATA_DIR=$MAIN_DIR/data

REGS='ant1st2s ant1st4s  ant2nd2s out button button_rtmod report
         reportaff  reportemo reportaff_rtmod reportemo_rtmod
         gvnant4s1  lvnant4s1 gvnant2s2 lvnant2s2
         gantparametric1 lantparametric1 gantparametric2 lantparametric2
         all_ant_param_posarous all_ant_param_negarous all_ant_param_arous all_ant_param_valence
         ant_param_arous ant_param_valence ant_param_posarous ant_param_negarous
         report_param_arous report_param_valence
         out_ant_param_arous out_ant_param_valence out_ant_param_posarous out_ant_param_negarous
         out_out_param_arous out_out_param_valence out_out_param_posarous out_out_param_negarous
         gvnout  lvnout goutparametric loutparametric outparametric
         param_arous param_valence param_posarous param_negarous param_arous_early
         param_posarous_1trprior param_negarous_1trprior param_posarous_2trprior param_negarous_2trprior
         param_posarous_1trlater param_negarous_1trlater param_posarous_2trlater param_negarous_2trlater
         LvS LvS2s LvS2sGint LvS2sLint
         pa1 pa2 pa3 na1 na2 na3
         arous1 arous2 arous3 valence1 valence2 valence3
         lparametric2s gparametric2s antparametric_lin antparametric_abs outparametric_lin outparametric_abs
         param_posarous_4s param_negarous_4s report_param_posarous report_param_negarous Anxious Angry Sad Calm Happy Excited
         ant_ev ant_ev_abs out_pe out_pe_pos out_pe_neg ant_p_gain ant_recent_avg_gain ant_recent_avg_winpercent
         ant_ev_pos ant_ev_neg out_pe_abs'


cd $SCRIPTS_DIR
./move_pupil_regs.sh

SUBJECTS=$DATA_DIR/subjects-midaffemo.txt

for SUBJECT in $(cat $SUBJECTS)
#for SUBJECT in ag240330
do
	echo; echo ---------------------------------------
	echo Making vectors for $SUBJECT
	echo ---------------------------------------
	echo

	cd $DATA_DIR/subjects/$SUBJECT/func_proc_MIDaffemo

    if [ ! -d regs ]; then
        mkdir regs
    fi

    echo

    rm regs/*_c.1D
    rm regs/*_tr.1D

    Rscript $SCRIPTS_DIR/makeRegs.R

    cd regs

	############### convolve task regressors ######################
    # this will take the hemodynamic delay into account

    for REG in $REGS
    do
        waver -dt 2.0 -GAM -input ${REG}_tr.1D > ${REG}_c.1D
    done # regs loop

    waver -dt 2.0 -GAM -input pupil_scaled.1D > pupil_scaled_c.1D
    waver -dt 2.0 -GAM -input luminance.1D > luminance_c.1D
    waver -dt 2.0 -GAM -input pupilx.1D > pupilx_c.1D
    waver -dt 2.0 -GAM -input pupily.1D > pupily_c.1D
    waver -dt 2.0 -GAM -input blink.1D > blink_c.1D
    waver -dt 2.0 -GAM -input saccade.1D > saccade_c.1D

done # subject loop
