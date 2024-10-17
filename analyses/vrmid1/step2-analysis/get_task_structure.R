source(here::here("load_libraries.R"))

my.df <-read.csv("../data/reduced_data.csv")%>%
  mutate(probe = ifelse(probe == "anti","out", 
                        ifelse(probe == "out","anti",NA)))

task_struct <- my.df%>%
  filter(Subject == "ay220916")%>%
  select(trial,trial.type,probe,ITI_duration)%>%
  unique()
  
rownames(task_struct) <- NULL

write_csv(task_struct,"../data/task_struct.csv")
