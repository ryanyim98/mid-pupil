#!/usr/bin/env Rscript

# written by Leili (initially July 2020 -- updated January 2023)

# This script is run from within the func_proc directory of each subject.

library(readr)
library(dplyr)
library(stringr)
library(zoo)

PATH_BEH = "../raw/"
# PATH_BEH = "../data/subjects/cw230819/raw/"
PATH_REG = "regs/REGNAME_tr.1D"

midmatrix_4runs <- NULL
my.runs = c("run-01", "run-02","run-03","run-04")

for (r in 1:4) {
  runnum = my.runs[r]
  PATH_FILE = paste0(PATH_BEH,list.files(PATH_BEH)[grepl(runnum,list.files(PATH_BEH)) & grepl(".csv",list.files(PATH_BEH))])
  
  midmatrix_1run <- read_csv(file = PATH_FILE)  
  
  if (runnum == "run-01" ){
    # create the 5 TRs at the tail of run 1
    end_run1 <- midmatrix_1run %>%
      subset(trial == 24) %>%
      .[1:5,] %>%
      mutate(TR = TR + 10)
  } else if (runnum == "run-02") {
    # create the 5 TRs at the tail of run 2
    end_run2 <- midmatrix_1run %>%
      subset(trial == 24) %>%
      .[1:5,] %>%
      mutate(TR = TR + 9)
  } else if (runnum == "run-03") {
    # create the 5 TRs at the tail of run 3
    end_run3 <- midmatrix_1run %>%
      subset(trial == 24) %>%
      .[1:5,] %>%
      mutate(TR = TR + 10)
  } else {
    # create the 5 TRs at the tail of run 4
    end_run4 <- midmatrix_1run %>%
      subset(trial == 24) %>%
      .[1:5,] %>%
      mutate(TR = TR + 11)
  }
  
  if (runnum == "run-01"){
    midmatrix_1run <- midmatrix_1run %>% 
      # insert 5 rows at the end of each run
      add_row(end_run1, .after = 240)
  } else if (runnum == "run-02"){
    midmatrix_1run <- midmatrix_1run %>% 
      # insert 5 rows at the end of each run
      add_row(end_run2, .after = 240)
  } else if (runnum == "run-03"){
    midmatrix_1run <- midmatrix_1run %>% 
      # insert 5 rows at the end of each run
      add_row(end_run3, .after = 240)
  } else {
    midmatrix_1run <- midmatrix_1run %>% 
      # insert 5 rows at the end of each run
      add_row(end_run4, .after = 240)
  }
  
  
  midmatrix_1run <- midmatrix_1run%>%  
    mutate(cue_value = ifelse(trialtype == 1, "-$0",ifelse(trialtype == 2, "-$1", ifelse(trialtype == 3, "-$5",
                                                                                         ifelse(trialtype == 4, "+$0",ifelse(trialtype == 5, "+$1","+$5"))))),
           trial_gain = if_else(trial_gain == "($1)", "-$1", 
                                if_else(trial_gain == "($5)","-$5", 
                                        if_else(trial_gain %in% c("$0","$1","$5"),paste0("+",trial_gain),trial_gain))),
           run = runnum) %>% 
    #ratings can't be scaled here because 4 runs are only concatenated later after this
    mutate(ant1st2s = ifelse(TR == 2, yes = 1, no = NA),
           ant1st4s = ifelse(TR >= 1 & TR <= 2, yes = 1, no = NA),
           ant2nd2s = ifelse(probe == 1 & TR == 7, yes = 1, no = NA),
           ant = ifelse(ant1st4s>0 | ant2nd2s>0,1,NA),
           out = ifelse(probe == 1 & TR == 9, yes = 1,no = ifelse(probe == 2 & TR == 4, yes = 1, no = 0)),
           report = ifelse(probe == 1 & TR >= 3 & TR <= 6, yes = 1,
                           no = ifelse(probe == 2 & TR >= 5 & TR <= 8, yes = 1, no = 0)), #could be improved by recording RT
           reportaff = ifelse(report == 1 & rating_type == 1, yes = 1,
                             no = 0),
           reportemo = ifelse(report == 1 & rating_type == 2, yes = 1,
                             no = 0),
           tr_of_interest = ifelse(probe == 1 & TR == 2, yes = 1,
                                   no = ifelse(probe == 2 & TR == 4, yes = 1, no = 0)),
           tr_of_interest_early = ifelse(probe == 1 & TR %in% c(1,2), yes = 1,
                                   no = ifelse(probe == 2 & TR %in% c(3,4), yes = 1, no = 0)),
           reportaff_rtmod = ifelse(probe == 1 & rating_type == 1 & report == 1 & 
                                      arating_t > 0 & TR <= 3 + 2 + arating_t/2, 
                                    yes = 1,
                                    no = ifelse(probe == 2 & rating_type == 1 & report == 1 & 
                                                  arating_t > 0 & TR <= 5 + 2 + arating_t/2, 
                                                yes = 1, no =ifelse(arating_t < 0 & rating_type == 1 & report == 1, yes = 1, no =0))),
           reportemo_rtmod = ifelse(probe == 1 & rating_type == 2 & erating_t > 0 & 
                                      report == 1 & TR <= 3 + erating_t/2, yes = 1,
                              no = ifelse(probe == 2 & rating_type == 2 & erating_t > 0 & report == 1 & 
                                            TR <= 5 + erating_t/2, 
                                          yes = 1, no =ifelse(erating_t < 0 & rating_type == 2 & report == 1,
                                                              yes = 1, no = 0))),
           gantparametric1 = ifelse(ant1st4s == 1 & trialtype >= 4, 
                          yes = as.numeric(gsub("\\$","",as.character(cue_value))), no = NA),
           lantparametric1 = ifelse(ant1st4s == 1 & trialtype < 4, 
                              yes = abs(as.numeric(gsub("\\$","",as.character(cue_value)))), no = NA),
           gantparametric2 = ifelse(ant2nd2s == 1 & trialtype >= 4, 
                                    yes = as.numeric(gsub("\\$","",as.character(cue_value))), no = NA),
           lantparametric2 = ifelse(ant2nd2s == 1 & trialtype < 4, 
                                    yes = abs(as.numeric(gsub("\\$","",as.character(cue_value)))), no = NA),
           lparametric2s = ifelse(rating_type == 1 & probe == 1 & TR == 2 & trialtype < 4, 
                                    yes =abs(as.numeric(gsub("\\$","",as.character(cue_value)))),
                                    ifelse(rating_type == 1 & probe == 2 & TR == 4 & trialtype < 4,
                                           yes=abs(as.numeric(gsub("\\$","",as.character(cue_value)))),no=NA)),
           gparametric2s = ifelse(rating_type == 1 & probe == 1 & TR == 2 & trialtype >= 4, 
                                  yes =abs(as.numeric(gsub("\\$","",as.character(cue_value)))),
                                  ifelse(rating_type == 1 & probe == 2 & TR == 4 & trialtype>= 4,
                                         yes=abs(as.numeric(gsub("\\$","",as.character(cue_value)))),no=NA)),
           antparametric_lin = ifelse(ant1st4s == 1, as.numeric(gsub("\\$","",as.character(cue_value))),
                                   NA),
           antparametric_abs = ifelse(ant1st4s == 1, abs(as.numeric(gsub("\\$","",as.character(cue_value)))),
                                       NA),
           outparametric_lin = ifelse(out == 1, as.numeric(gsub("\\$","",as.character(cue_value))),
                                      NA),
           outparametric_abs = ifelse(out == 1, abs(as.numeric(gsub("\\$","",as.character(cue_value)))),
                                      NA),
           gvnant4s1 = ifelse(ant1st4s == 1 & trialtype == 6,
                           yes = 1, no = ifelse(ant1st4s == 1 & trialtype == 4,
                                                yes = -1, no = NA)),
           lvnant4s1 = ifelse(ant1st4s == 1 & trialtype == 3,
                           yes = 1, ifelse(ant1st4s == 1 & trialtype == 1,
                                           yes = -1, no = NA)),
           gvnant2s2 = ifelse(ant2nd2s == 1 & trialtype == 6,
                              yes = 1, no = ifelse(ant2nd2s == 1 & trialtype == 4,
                                                   yes = -1, no = NA)),
           lvnant2s2 = ifelse(ant2nd2s == 1 & trialtype == 3,
                              yes = 1, ifelse(ant2nd2s == 1 & trialtype == 1,
                                              yes = -1, no = NA)),
           LvS = ifelse(ant1st4s == 1 & cue_size == 2, yes = 1, 
                        no = ifelse(ant1st4s == 1 & cue_size == 1, yes = -1, no = NA)),
           LvS2s = ifelse(TR == 1 & cue_size == 2, yes = 1, 
                          no = ifelse(TR == 1 & cue_size == 1, yes = -1, no = NA)),
           LvS2sGint = ifelse(trialtype == 6 & LvS2s == 1, yes = 1,  #large small  * +$5 vs. +$0
                                 no = ifelse(trialtype == 6 & LvS2s == -1, yes = -1,
                                             no = ifelse(trialtype == 4 & LvS2s == 1, yes = -1,
                                                         no = ifelse(trialtype == 4 & LvS2s == -1, yes = 1, no = 0)))),
           LvS2sLint = ifelse(trialtype == 3 & LvS2s == 1, yes = 1,  #large small  * -$5 vs. -$0
                              no = ifelse(trialtype == 3 & LvS2s == -1, yes = -1,
                                          no = ifelse(trialtype == 1 & LvS2s == 1, yes = -1,
                                                      no = ifelse(trialtype == 1 & LvS2s == -1, yes = 1, no = 0)))),
           # press (the third or 7th TR)
           button = ifelse(probe == 1 & TR == 8, yes = 1, no = ifelse(probe == 2 & TR == 3, yes = 1, no = 0)),
           button_rtmod = ifelse(button == 1 & hit == 1, yes = rt, 
                                 ifelse(button == 1 & hit == 0, yes = target_ms + 0.1,
                                        no = 0)), #shorter, the better
           ant1st4s_rtmod = ifelse(ant1st4s == 1 & hit == 1, yes = rt, 
                                 ifelse(ant1st4s == 1 & hit == 0, yes = target_ms + 0.1,
                                        no = 0)), #shorter, the better
           # outcome (not differentiate two probe types)
           gvnout = ifelse(out == 1 & trialtype == 6 & hit == 1, 
                           yes = 1, no = ifelse(out == 1 & trialtype == 6 & hit == 0, 
                                                yes = -1, no = NA)),
           lvnout = ifelse(out == 1 & trialtype == 3 & hit == 0, 
                           yes = 1, no = ifelse(out == 1 & trialtype == 3 & hit == 1, 
                                                yes = -1, no = NA)),
           outparametric = ifelse(out == 1, yes = as.numeric(gsub("\\$","",as.character(trial_gain))), no = NA),
           goutparametric = ifelse(out == 1 & trialtype >= 4, yes = as.numeric(gsub("\\$","",as.character(trial_gain))), no = NA),
           loutparametric = ifelse(out == 1 & trialtype <= 3, 
                                   yes = abs(as.numeric(gsub("\\$","",as.character(trial_gain)))), no = NA)
    )
  
  midmatrix_4runs <- bind_rows(midmatrix_4runs,
                               midmatrix_1run %>% mutate(run = r,
                                                         trial_all = trial + (run-1)*24))
} # run loop


# after combining 4 runs, scale the ratings
midmatrix_trials <- unique(midmatrix_4runs %>% 
  filter(TR == 1)) %>% 
  ungroup() %>% 
  mutate(aratings = arating - mean(arating, na.rm = T),
         vratings = vrating - mean(vrating, na.rm = T),
         posarous = (aratings + vratings)/2,
         negarous = (aratings - vratings)/2) %>% 
  mutate(n_Anxious = sum(erating == 1,na.rm = T),
         n_Angry = sum(erating == 2,na.rm = T),
         n_Sad = sum(erating == 3,na.rm = T),
         n_Calm = sum(erating == 4,na.rm = T),
         n_Happy = sum(erating == 5,na.rm = T),
         n_Excited = sum(erating == 6,na.rm = T),
         n_Emotions = n_Anxious+n_Angry+n_Sad+n_Calm+n_Happy+n_Excited)

if (length(unique(midmatrix_trials$posarous)) < 4) { #3 plus NA
  qpa = NA
  return('PA does not have enough levels')
  break
} else if (length(unique(midmatrix_trials$posarous)) == 4) { #3 plus NA
  qpa<-unique(midmatrix_trials$posarous)%>%sort()
} else {
  qpa<-quantile(unique(midmatrix_trials$posarous),probs = c(1/3,2/3),na.rm = T) #quantile(midmatrix_trials$posarous,probs = c(1/3,2/3),na.rm = T)
}


if (length(unique(midmatrix_trials$negarous)) < 4) { #3 plus NA
  qna = NA
  return('NA does not have enough levels')
  break
} else if (length(unique(midmatrix_trials$negarous)) == 4) {
  qna<-unique(midmatrix_trials$negarous)%>%sort()
} else {
  qna<-quantile(unique(midmatrix_trials$negarous),probs = c(1/3,2/3),na.rm = T)#quantile(midmatrix_trials$negarous,probs = c(1/3,2/3),na.rm = T)
}

if (length(unique(midmatrix_trials$aratings)) < 4 & length(unique(midmatrix_trials$aratings)) >= 2) { #3 plus NA
  # return('Arousal does not have enough levels')
  # break
  qa<-unique(midmatrix_trials$aratings)%>%sort() #only 2 levels
} else if (length(unique(midmatrix_trials$aratings)) == 4) {
  qa<-unique(midmatrix_trials$aratings)%>%sort()
} else {
  qa<-quantile(unique(midmatrix_trials$aratings),probs = c(1/3,2/3),na.rm = T)#quantile(midmatrix_trials$aratings,probs = c(1/3,2/3),na.rm = T)
}

if (length(unique(midmatrix_trials$vratings)) < 4 & length(unique(midmatrix_trials$vratings)) >= 2) { #3 plus NA
  # return('Valence does not have enough levels')
  # break
  qv<-unique(midmatrix_trials$vratings)%>%sort() #only 2 levels
} else if (length(unique(midmatrix_trials$vratings)) == 4) {
  qv<-unique(midmatrix_trials$vratings)%>%sort()
} else {
  qv<-quantile(unique(midmatrix_trials$vratings),probs = c(1/3,2/3),na.rm = T)#quantile(midmatrix_trials$aratings,probs = c(1/3,2/3),na.rm = T)
}


#calculate RPE, EV etc.
midmatrix_trials <- midmatrix_trials %>% 
  mutate(cue_value_num = as.numeric(paste0(substring(cue_value,1,1),substring(cue_value,3,3))),
         trial_gain_num = if_else(hit == 1 & cue_value_num >= 0, cue_value_num, 
                                  if_else(hit == 1 & cue_value_num < 0, 0,
                                          if_else(hit == 0 & cue_value_num >= 0, 0,
                                                  if_else(hit == 0 & cue_value_num < 0, cue_value_num, NA)))),
         last_hit = lag(hit,1),
         last_cue = lag(cue_value,1),
         last_gain = lag(trial_gain,1),
         last_gain_num = lag(trial_gain_num,1),
         last_cue_num = lag(cue_value_num,1),
         last_total_winpercent = lag(total_winpercent,1),
         last_binned_winpercent = lag(binned_winpercent,1),
         recent_avg_gain = lag(rollmean(trial_gain_num, 5,fill = c(0,0,0), align = "right"),1),
         gain_contrast = trial_gain_num - last_gain_num,
         recent_avg_winpercent = lag(rollmean(hit, 5, fill = 0, align = "right"),1),
         rt = ifelse(rt == -1, target_ms+0.1,rt),#get proxy for rt since we did not record missing rt
         last_rt = lag(rt,1))
  

p_gain = NULL
ev = NULL
pe = NULL
r = NULL #actual reward
c = NULL #cue value
hit = NULL
ratio = NULL
for (t in 1:nrow(midmatrix_trials)){
  if (t %in% c(1,25,49,73)){ #fit each run separately
    p_gain[t] = 0.5
  }
  c[t] = midmatrix_trials$cue_value_num[t]
  r[t] = midmatrix_trials$trial_gain_num[t]
  hit[t] = midmatrix_trials$hit[t]
  ev[t] =  ifelse(c[t] > 0,c[t] * p_gain[t],
                  ifelse(c[t] < 0, c[t] * (1-p_gain[t]), 0))
  pe[t] =  r[t] - ev[t]
  ratio[t] = ifelse(c[t]<=0,0,pe[t]/abs(c[t]))
  if (t < nrow(midmatrix_trials)){
    p_gain[t+1]=p_gain[t]+ 0.2 * ratio[t]
  }
}

rpe_calc = data.frame(ev,pe,p_gain) %>% mutate(trial_all = row_number())


midmatrix_4runs <- midmatrix_4runs %>% 
  left_join(midmatrix_trials %>%
              select(trial_all,aratings,vratings,posarous,negarous,n_Anxious:n_Emotions,
                     cue_value_num,
                     last_gain_num,recent_avg_gain,recent_avg_winpercent), by = "trial_all") %>% 
  left_join(rpe_calc %>% 
              select(trial_all,ev,pe,p_gain), by = "trial_all") %>% 
  ungroup() %>% 
  group_by(run) %>% 
  mutate(onetrprior = lead(tr_of_interest,1),
         twotrprior = lead(tr_of_interest,2),
         threetrprior = lead(tr_of_interest,3),
         onetrlater = lag(tr_of_interest,1),
         twotrlater = lag(tr_of_interest,2),
         threetrlater = lag(tr_of_interest,3)) %>% 
  ungroup() %>% 
  mutate(all_ant_param_posarous = ifelse(rating_type == 1 & ant1st4s == 1,  #these are during ANTICIPATION
                                         yes = posarous, no = NA),
         all_ant_param_negarous = ifelse(rating_type == 1 & ant1st4s == 1,  #these are during ANTICIPATION
                                         yes = negarous, no = NA),
         all_ant_param_arous = ifelse(rating_type == 1 & ant1st4s == 1,  #these are during ANTICIPATION
                                         yes = aratings, no = NA),
         all_ant_param_valence = ifelse(rating_type == 1 & ant1st4s == 1,  #these are during ANTICIPATION
                                         yes = vratings, no = NA),
         ant_param_arous = ifelse(rating_type == 1 & probe == 1 & ant1st4s == 1,  #these are during ANTICIPATION
                           yes = aratings, no = NA), #try the scaled version first
         ant_param_valence = ifelse(rating_type == 1 & probe == 1 & ant1st4s == 1, 
                                     yes = vratings, no = NA),
        report_param_arous = ifelse(rating_type == 1 & probe == 1 & TR >= 5 & TR <= 6,  #these are during RATINGS
                                  yes = aratings, 
                                  no = ifelse(rating_type == 1 & probe == 2 & TR >= 7 & TR <= 8,
                                              yes = aratings, 
                                              no = NA)), #try the scaled version first
        report_param_valence = ifelse(rating_type == 1 & probe == 1 & TR >= 3 & TR <= 4, #these are during RATINGS
                                    yes = vratings, 
                                    no = ifelse(rating_type == 1 & probe == 2 & TR >= 5 & TR <= 6,
                                                yes = vratings, 
                                                no = NA)),
        report_param_posarous = ifelse(rating_type == 1 & probe == 1 & TR ==4,  #these are during RATINGS
                                    yes = posarous, 
                                    no = ifelse(rating_type == 1 & probe == 2 & TR ==6,
                                                yes = posarous, 
                                                no = NA)), #try the scaled version first
        report_param_negarous = ifelse(rating_type == 1 & probe == 1 & TR == 4, #these are during RATINGS
                                      yes = negarous, 
                                      no = ifelse(rating_type == 1 & probe == 2 & TR == 6,
                                                  yes = negarous, 
                                                  no = NA)),
          ant_param_posarous = ifelse(rating_type == 1 & probe == 1 & ant1st4s == 1,  #these are during ANTICIPATION
                                      yes = posarous, no = NA), #try the scaled version first
          ant_param_negarous = ifelse(rating_type == 1 & probe == 1 & ant1st4s == 1, 
                                      yes = negarous, no = NA),
          out_ant_param_arous = ifelse(rating_type == 1 & probe == 2 & ant1st4s == 1, 
                                   yes = aratings, no = NA), #try the scaled version first
          out_ant_param_valence = ifelse(rating_type == 1 & probe == 2 & ant1st4s == 1, 
                                     yes = vratings, no = NA),
          out_ant_param_posarous = ifelse(rating_type == 1 & probe == 2 & ant1st4s == 1,  #these are during ANTICIPATION
                                      yes = posarous, no = NA), #try the scaled version first
          out_ant_param_negarous = ifelse(rating_type == 1 & probe == 2 & ant1st4s == 1, 
                                      yes = negarous, no = NA),
         out_out_param_arous = ifelse(rating_type == 1 & probe == 2 & out == 1, 
                                      yes = aratings, no = NA), #try the scaled version first
         out_out_param_valence = ifelse(rating_type == 1 & probe == 2 & out == 1, 
                                        yes = vratings, no = NA),
         out_out_param_posarous = ifelse(rating_type == 1 & probe == 2 & out == 1,  #these are during OUTCOME
                                         yes = posarous, no = NA), #try the scaled version first
         out_out_param_negarous = ifelse(rating_type == 1 & probe == 2 & out == 1, 
                                         yes = negarous, no = NA),
         param_arous = ifelse(rating_type == 1 & tr_of_interest == 1, yes = aratings, no = 0),
         param_arous_early = ifelse(rating_type == 1 & tr_of_interest_early == 1, yes = aratings, no = 0),
         param_valence = ifelse(rating_type == 1 & tr_of_interest == 1, yes = vratings, no = 0),
         pa1 = ifelse(tr_of_interest == 1 & posarous <= qpa[1], yes = 1, no = 0),
         pa2 = ifelse(tr_of_interest == 1 & posarous > qpa[1] & posarous <= qpa[2], yes = 1, no = 0),
         pa3 = ifelse(tr_of_interest == 1 & posarous > qpa[2], yes = 1, no = 0),
         na1 = ifelse(tr_of_interest == 1 & negarous <= qna[1], yes = 1, no= 0),
         na2 = ifelse(tr_of_interest == 1 & negarous > qna[1] & negarous <= qna[2], yes = 1,no = 0),
         na3 = ifelse(tr_of_interest == 1 & negarous > qna[2], yes = 1,no = 0),
         
         arous1 = ifelse(tr_of_interest == 1 & aratings <= qa[1], yes = 1, no = 0),
         arous2 = ifelse(tr_of_interest == 1 & aratings > qa[1] & aratings <= qa[2], yes = 1,no = 0),
         arous3 = ifelse(tr_of_interest == 1 & aratings > qa[2], yes = 1,no = 0),
         
         valence1 = ifelse(tr_of_interest == 1 & vratings <= qv[1], yes = 1,no = 0),
         valence2 = ifelse(tr_of_interest == 1 & vratings > qv[1] & vratings <= qv[2], yes = 1,no = 0),
         valence3 = ifelse(tr_of_interest == 1 & vratings > qv[2], yes = 1,no = 0),
         
         param_posarous = ifelse(tr_of_interest == 1, yes = posarous,  no = 0),
         param_negarous = ifelse(tr_of_interest == 1, yes = negarous,  no = 0),
        
        Anxious = ifelse(tr_of_interest == 1 & erating == 1, yes = 1, 
                         no = ifelse(tr_of_interest == 1 & erating %in% c(2:6), -1*n_Anxious/(n_Emotions-n_Anxious), 0)),
        Angry = ifelse(tr_of_interest == 1 & erating == 2, yes = 1, 
                       no = ifelse(tr_of_interest == 1 & erating %in% c(1,3:6), -1*n_Angry/(n_Emotions-n_Angry), 0)),
        Sad = ifelse(tr_of_interest == 1 & erating == 3, yes = 1, 
                     no = ifelse(tr_of_interest == 1 & erating %in% c(1:2,4:6),  -1*n_Sad/(n_Emotions-n_Sad), 0)),
        Calm = ifelse(tr_of_interest == 1 & erating == 4, yes = 1, 
                      no = ifelse(tr_of_interest == 1 & erating %in% c(1:3,5:6),  -1*n_Calm/(n_Emotions-n_Calm), 0)),
        Happy = ifelse(tr_of_interest == 1 & erating == 5, yes = 1, 
                       no = ifelse(tr_of_interest == 1 & erating %in% c(1:4,6),  -1*n_Happy/(n_Emotions-n_Happy), 0)),
        Excited = ifelse(tr_of_interest == 1 & erating == 6, yes = 1, 
                         no = ifelse(tr_of_interest == 1 & erating %in% c(1:5), -1*n_Excited/(n_Emotions-n_Excited), 0)),
        param_posarous_1trprior = ifelse(onetrprior == 1, yes = posarous, no = 0),
        param_negarous_1trprior = ifelse(onetrprior == 1, yes = negarous, no = 0),
        param_posarous_2trprior= ifelse(twotrprior == 1, yes = posarous, no = 0),
        param_negarous_2trprior = ifelse(twotrprior == 1, yes = negarous, no = 0),
        param_posarous_1trlater = ifelse(onetrlater == 1, yes = posarous, no = 0),
        param_negarous_1trlater = ifelse(onetrlater == 1, yes = negarous, no = 0),
        param_posarous_2trlater= ifelse(twotrlater == 1, yes = posarous, no = 0),
        param_negarous_2trlater = ifelse(twotrlater == 1, yes = negarous, no = 0),
        param_posarous_4s = ifelse(probe == 1 & TR >= 1 & TR <= 2, yes = posarous, 
                                no = ifelse(probe == 2 & TR >= 3 & TR <= 4, yes = posarous, no = 0)),
        param_negarous_4s = ifelse(probe == 1 & TR >= 1 & TR <= 2, yes = negarous, 
                                no = ifelse(probe == 2 & TR >= 3 & TR <= 4, yes = negarous, no = 0)),
         run1 = ifelse(run == 1, 1, 0),
         run2 = ifelse(run == 2, 1, 0),
         run3 = ifelse(run == 3, 1, 0),
         run4 = ifelse(run == 4, 1, 0),
        ant_ev = ifelse(ant1st4s == 1, yes = ev, no = 0),
        ant_ev_abs = abs(ant_ev),
        ant_ev_pos = ifelse(ev > 0, ev, 0),
        ant_ev_neg = ifelse(ev < 0, abs(ev), 0),
        out_pe = ifelse(out == 1, yes = pe, no = 0),
        out_pe_abs = abs(out_pe),
        out_pe_pos = ifelse(out == 1 & out_pe > 0, yes = pe, no = 0),
        out_pe_neg = ifelse(out == 1 & out_pe < 0, yes = abs(pe), no = 0),
        ant_p_gain = ifelse(ant1st4s == 1, yes = p_gain, no = 0),
        ant_recent_avg_gain = ifelse(ant1st4s == 1, yes = recent_avg_gain, no = 0),
        ant_recent_avg_winpercent = ifelse(ant1st4s == 1, yes = recent_avg_winpercent, no = 0),
        ant_abs =ifelse(ant1st4s == 1, abs(as.numeric(gsub("\\$","",as.character(cue_value)))),0),
        ant_valence =ifelse(ant1st4s == 1, as.numeric(gsub("\\$","",as.character(cue_value))),0))

a<-midmatrix_4runs %>%select(onetrprior,twotrprior,posarous,param_posarous_1trprior,param_posarous_2trprior)

regs <- midmatrix_4runs %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  select(ant1st2s, ant1st4s, ant2nd2s,out,button,button_rtmod,report,
         reportaff, reportemo,reportaff_rtmod,reportemo_rtmod,
         gvnant4s1, lvnant4s1,gvnant2s2, lvnant2s2,
         gantparametric1,lantparametric1, gantparametric2,lantparametric2,lparametric2s,gparametric2s,
         all_ant_param_posarous,all_ant_param_negarous,all_ant_param_arous,all_ant_param_valence,
         ant_param_arous,ant_param_valence,ant_param_posarous,ant_param_negarous,
         report_param_arous,report_param_valence,
         out_ant_param_arous,out_ant_param_valence,out_ant_param_posarous,out_ant_param_negarous,
         out_out_param_arous,out_out_param_valence,out_out_param_posarous,out_out_param_negarous,
         gvnout, lvnout,goutparametric,loutparametric,outparametric,
         param_arous,param_valence,param_posarous,param_negarous,
         param_arous_early,
         LvS,LvS2s,LvS2sGint,LvS2sLint,run1,run2,run3,run4,
         ant1st4s_rtmod,
         pa1,pa2,pa3,na1,na2,na3,
         arous1,arous2,arous3,valence1,valence2,valence3,
         antparametric_lin,antparametric_abs,outparametric_lin,outparametric_abs,
         param_posarous_4s,param_negarous_4s,report_param_posarous,report_param_negarous,
         param_posarous_1trprior,param_negarous_1trprior,param_posarous_2trprior,param_negarous_2trprior,
         param_posarous_1trlater,param_negarous_1trlater,param_posarous_2trlater,param_negarous_2trlater,
         Anxious,Angry,Sad,Calm,Happy,Excited,
         ant_ev,ant_ev_abs,out_pe, ant_p_gain,ant_recent_avg_gain,ant_recent_avg_winpercent,
         out_pe_pos,out_pe_neg,ant_ev_pos,ant_ev_neg,out_pe_abs)
    
n_regs = length(colnames(regs))

# write each reg column into a separate .1D file
for(i in 1:n_regs) {
  write.table(regs[i],
              file = str_replace(PATH_REG, "REGNAME", colnames(regs[i])),
              sep = "\n", row.names = FALSE, col.names = FALSE)
  } # reg loop


# Summarise data availability
# Define the list of emotions
emotions <- c("Anxious", "Angry", "Sad", "Calm", "Happy", "Excited")

# Loop through each emotion
for (emotion in emotions) {
  # Check if the file for the current emotion exists
  file_name <- paste0("~/dopa/midaffemo/glm_results/",emotion, ".txt")
  if (!file.exists(file_name)) {
    # If the file does not exist, create it
    file.create(file_name)
  }
  n_emotion = midmatrix_trials[paste0("n_",emotion)][1,1]
  print(n_emotion)
  # Check if the sum of the current emotion column is 0
  if (n_emotion < 7) {
    # If sum is 0, append "SUB" to the file without overwriting
    sub=basename(dirname(getwd()))
    write(sub, file = file_name, append = TRUE)
  }
}


