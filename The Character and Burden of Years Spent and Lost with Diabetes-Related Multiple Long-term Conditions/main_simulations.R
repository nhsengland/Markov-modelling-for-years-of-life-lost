# Preamble
library(dplyr)
library(tidyr)
library(expm)
library(purrr)
source("markov_helper_functions.R")
library(ggplot2)
library(scales)
library(readr)

con <- DBI::dbConnect(odbc::odbc(),
                      Driver = "ODBC Driver 18 for SQL Server",
                      Server = ,
                      Trusted_connection = "yes")


# Data load
  #' A table of subsegment IDs and (1/NA) corresponding condition flags for all 
  #' conditions
  raw_ss <- readRDS("markov_subsegment_v22data_v21pop_MF.RDS") 
  
  #' A table of count data for unique age at transition, from state, do state, 
  #' and sex of individual
  raw_markov <- readRDS("markov_transitions_v22data_v21pop_MF.RDS")
  
  #' The outputs of an alive_death markov model (otherwise known as survival 
  #' model). Each line is from state (Alive/Dead), age from, age to, probability 
  #' of being in state "Alive", probability of being in state "Death", and
  #' years betwen from and to, labelled 'age'. Model is the same as that 
  #' calculated in the two vars function below.
  raw_alive_dead <- readRDS("alive_death_all_v22data_v21pop.RDS") 



  
  raw_markov <- raw_markov %>%
    filter(1 == 1) #For whole population
  #filter(Sex == 1) #For Male simulation
  #filter(Sex == 2) #For Male simulation
  
  
 

### Continuation
# Analytics pipeline
twovars <- function(metric1, metric2){
  
  # Define markov states in order
  health_levels <-  c("Neither morbidity",#Better described as 'not both'
                      "Both morbidities", #Better described as 'at least' both
                      "Death")
  
  # Categorise each subsegment into one of the health states
  df_state <- raw_ss %>% 
    mutate(state = case_when(.data[[metric1]] == 1 & 
                               .data[[metric2]] ==1 ~ "Both morbidities",
                             TRUE ~ "Neither morbidity")) %>% 
    select(Subsegment_Combination_ID, state) %>% 
    rbind(c(NA_integer_, "Death")) %>% 
    mutate(Subsegment_Combination_ID = as.integer(Subsegment_Combination_ID))
  
  # Map markov subsegment transitions to health state transitions
  df_pre_processed <- raw_markov %>% 
    left_join(df_state, by = c("xS" = "Subsegment_Combination_ID"), 
              suffix= c("x", "y")) %>% 
    left_join(df_state, by = c("yS" = "Subsegment_Combination_ID"),  
              suffix= c("x", "y")) %>% 
    group_by(statex, statey, xA) %>% 
    summarise(n = sum(n)) %>% 
    arrange(xA, statex, statey)
  
  # Cleanse and supplement data for Markov processing
  dat_transitions_long <- df_pre_processed %>% 
    rename(age = xA,
           from = statex, 
           to = statey) %>% 
    filter(age < 101) %>% 
    filter(!is.na(n)) %>% 
    filter(!is.na(to)) %>%   
    arrange(age, from, to) %>% 
    as.data.frame() %>% 
    add_row(age = 9999, from = "Death", to ="Death") %>% 
    complete(age, from, to, fill = list(n = 0)) %>% 
    filter(age != 9999) %>% 
    as.data.frame() %>% 
    mutate(n = if_else(from == "Death" & to == "Death", 1L, n))
  
  # Ensure multiple rows for the same transitions states are amalgamated
  dat_transitions_long_fix <- dat_transitions_long %>%  
    group_by(age, from) %>%  
    mutate(tot_age_from = sum(n)) %>% 
    ungroup() %>%  
    mutate(n = if_else(tot_age_from == 0 & from == to, 1L, n)) %>% 
    select(-tot_age_from)
  
  dat_trans_stage2 <- dat_transitions_long_fix %>% 
    calculate_transition_probabilities()

  generate_all_table_internal <- function(df){
    compute_at_age <- function(i, df){
      print(paste0("Doing my best to compute ", i))
      dat_r <- lapply(i:100, get_age_transition_matrix, df = df)  %>% 
        purrr::accumulate(matrixproduct) %>% 
        map_df(revert_df_trans, .id = "age") %>% 
        select(age, from, health_levels) %>% 
        mutate(age = as.numeric(age)) %>% 
        mutate(age_from = i) %>% 
        mutate(age_to = age_from + age)
    }
    out <- map_df(0:100, compute_at_age, df = df) %>% 
      pivot_longer(cols = `Neither morbidity`:`Death`, 
                   names_to = "to", values_to = "prop")
  }
  
  outlist <- lapply(0:100, get_age_transition_matrix, df = dat_trans_stage2)
  ww <- generate_all_table_internal(dat_trans_stage2)  %>% 
    mutate(type = paste0(metric1, "-", metric2))
  
  readr::write_csv(ww, paste0("inter_mf_", dat_LTC_counts[i,1], "-", 
                              dat_LTC_counts[i,2], ".csv"))

  
  fun_yll <- function(agei) {
    measure_mm <- ww %>% 
      filter(from == "Both morbidities") %>% 
      filter(age_from == agei) %>% 
      filter(to == "Both morbidities") %>% 
      pull(prop) %>% 
      sum()
    
    measure_alive <- raw_alive_dead %>% 
      pivot_longer(Alive:Death, names_to = "to", values_to = "prop") %>% 
      filter(from == "Alive") %>% 
      filter(age_from == agei) %>% 
      filter(to == "Alive") %>% 
      pull(prop) %>% sum()
    
    return(measure_alive - measure_mm)
  }
  
  fun_ymm <- function(agei){
    measure_mm <- ww %>% 
      filter(from == "Both morbidities") %>% 
      filter(age_from == agei) %>% 
      filter(to == "Both morbidities") %>% 
      pull(prop) %>% 
      sum()
    
    return(measure_mm)
  }
  
  fun_le_mean <- function(agei){
    measure_alive <- raw_alive_dead %>% 
      pivot_longer(Alive:Death, names_to = "to", values_to = "prop") %>% 
      filter(from == "Alive") %>% 
      filter(age_from == agei) %>% 
      filter(to == "Alive") %>% 
      pull(prop) %>% sum()
    
    return(measure_alive)
  }
  
  fun_inner_le_mean <- function(agei){
    ww %>% 
      filter(from == "Neither morbidity") %>% 
      filter(age_from == agei) %>% 
      filter(to != "Death") %>%  
      group_by(age_to) %>% 
      summarise(prop = sum(prop)) %>%  
      pull(prop) %>% 
      sum()
  }
  
  fun_le_median <- function(agei){
    rr <- raw_alive_dead %>%  
      filter(from == "Alive") %>% 
      filter(age_from == agei) 
    
    out <- ifelse(nrow(rr) > 1, approx(rr$Alive, rr$age_to, xout = 0.5)$y, -100)
    return(out)
  }
  
  prop_acquiring <- (filter(dat_trans_stage2, 
                            from == "Neither morbidity" & 
                              to == "Both morbidities") %>% 
                       pull(prop)) * (ww %>% filter(age_from == 0) %>% 
       filter(from == "Neither morbidity" & to == "Neither morbidity") %>% 
       pull(prop))
  
  prop_acquiring_wrong <- dat_transitions_long_fix %>% 
    filter(from == "Neither morbidity" & to == "Both morbidities") %>% 
    mutate(prop = n/sum(n)) %>% 
    pull(prop)
  
  age_over <- seq(0,100, 1)  
  outout <- data.frame(age = age_over,
                       life_exp_mean = map_dbl(age_over, fun_le_mean),
                       life_exp_inner_mean = map_dbl(age_over, fun_inner_le_mean),
                       life_exp_median = map_dbl(age_over, fun_le_median),
                       yll = map_dbl(age_over, fun_yll), 
                       ymm = map_dbl(age_over, fun_ymm),
                       age_get = prop_acquiring/sum(prop_acquiring),
                       age_get_wrong = prop_acquiring_wrong,
                       prop_acquiring = prop_acquiring,
                       type = paste0(metric1, "-", metric2))
  
  
  
  
  return(outout)
}



# Source a pre canned list of all condition combinations, where the prevalence 
# of each condition pair is captured to simulate most common pairs first
raw_LTC_counts <- DBI::dbGetQuery(con, "SELECT LTC_1, LTC_2, N FROM 
                  [].[].[tbl_Agg_Risk_Ratio_All]
                  WHERE ORG_NAME = 'England'")

dat_LTC_counts <- raw_LTC_counts %>%
  mutate(LTC_1 = factor(LTC_1)) %>%
  mutate(LTC_2 = factor(LTC_2, levels = levels(.$LTC_1))) %>%
  filter(as.numeric(LTC_1) < as.numeric(LTC_2)) %>%
  arrange(desc(N)) %>%
  mutate(LTC_1 = as.character(LTC_1)) %>%
  mutate(LTC_2 = as.character(LTC_2)) %>%
  mutate(completed = "FALSE")


# Do not simulate "severe" forms of LTCS (which are not included in analyses)  
dat_LTC_counts <- filter(dat_LTC_counts, LTC_1 == "Diabetes" | 
                           LTC_2 == "Diabetes") %>% 
  filter(!LTC_1 %in% c("Intermediate_Frailty_Risk_HFRS",
                       "Severe_COPD",
                       "High_Frailty_Risk_HFRS",
                       "Severe_Heart_Failure",
                       "Incurable_Cancer",
                       "End_Stage_Renal_Failure",
                       "Liver_Failure")) %>% 
  filter(!LTC_2 %in% c("Intermediate_Frailty_Risk_HFRS",
                       "Severe_COPD",
                       "High_Frailty_Risk_HFRS",
                       "Severe_Heart_Failure",
                       "Incurable_Cancer",
                       "End_Stage_Renal_Failure",
                       "Liver_Failure"))


# Iterate results generation over all condition pairs
for(i in 1:nrow(dat_LTC_counts)){
  print(i)
  gc()
  if(dat_LTC_counts[i,"completed"] == "FALSE"){
    tryCatch({
      out <- twovars(dat_LTC_counts[i,1], dat_LTC_counts[i,2])
      readr::write_csv(out, paste0("res_", dat_LTC_counts[i,1], "-",
                                   dat_LTC_counts[i,2], ".csv"))
      dat_LTC_counts[i,"completed"] <- "TRUE"
    },
    error = function(e){message('Error!!!')
      dat_LTC_counts[i,"completed"] <- "Error"})
    
  }
  saveRDS(dat_LTC_counts, 
          file = "dat_LTC_combos_all_v22_data_v21pop_MF_full.RDS")
}





