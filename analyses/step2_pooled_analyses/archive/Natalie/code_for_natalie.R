a <- fmri3_beh_data %>% filter(subject == "ab240324") %>% 
  ungroup() %>% 
  mutate(valence = as.numeric(scale(vrating)),
         arousal = as.numeric(scale(arating)),
         PosArous = (valence + arousal),
         NegArous = (valence-arousal))

b <- df.1sub %>% 
  rowwise() %>% 
  mutate(run_number = as.numeric(strsplit(run,"0")[[1]][2])) %>% 
  mutate(trial = trial + 24 * (run_number-1)) %>% 
  ungroup() %>% 
  mutate(tr_all = row_number())

reg_data <- left_join(a,b,by = "trial") %>% 
  relocate(tr)

reg_data_ant <- reg_data %>% filter(tr == 5, probe == 1, rating_type == 1)

reg_data_out <- reg_data %>% filter(tr == 7, probe == 2, rating_type == 1)

