library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)
library(readr)


# Set results folder
result_set <-  "results_full"
#result_set <-  "results_male"
#result_set <-  "results_female"

base <- "./Male_Female_Markov_revised_transition_data/"

setwd(paste0(base, result_set))



inter_file_list <- list.files(pattern = "inter_mf_")

df <- data.frame(name=character(),
                 person_yll=numeric(), 
                 community_yll=numeric(), 
                 person_ymm=numeric(), 
                 community_ymm=numeric(), 
                 get=numeric(),
                 median_onset_age=numeric(),
                 quartile_onset_age=numeric(),
                 stringsAsFactors=FALSE) 

df_inter_ymm_yll <- data.frame(name=character(),
                               age=numeric(),
                               person_yll=numeric(), 
                               person_ymm=numeric(), 
                               stringsAsFactors=FALSE)

for(i in inter_file_list){
  
  file <- i
  
  out_name <-  gsub(x = file, "inter_mf_", "") %>% gsub(pattern = ".csv", 
                                                        replacement = "", .)
  
  temp_inter <- read.csv(file)
  temp_res <- readr::read_csv(paste0("res_", out_name, ".csv"))
  
  inter_c_prop_moving <- temp_inter %>% 
    filter(age_from + 1 == age_to) %>%
    filter(age_to <= 101) %>% 
    filter(from == "Neither morbidity" & to == "Both morbidities") %>% 
    pull(prop)
  
  
  
  inter_c_healthy <- temp_inter %>% filter(age_from == 0) %>% 
    filter(age_to <=100) %>% 
    filter(from == "Neither morbidity" & to == "Neither morbidity") %>% 
    pull(prop)
  
  inter_c_moving <- inter_c_prop_moving * c(1,inter_c_healthy)
  
  out_new <- select(temp_res, age, yll, ymm) %>% 
    mutate(age_get = inter_c_moving) 
  
  out_community_yll <- mutate(out_new, val = yll*age_get*1000) %>% 
    pull(val) %>% sum()
  out_community_ymm <- mutate(out_new, val = ymm*age_get*1000) %>% 
    pull(val) %>% sum()
  out_person_yll <- mutate(out_new, val = yll*age_get/sum(age_get)) %>%
    pull(val) %>% sum()
  out_person_ymm <- mutate(out_new, val = ymm*age_get/sum(age_get)) %>%
    pull(val) %>% sum()
  out_get <- sum(out_new$age_get)
  
  
  median_onset_age <- filter(temp_res, type == out_name) %>% 
    arrange(age) %>% 
    mutate(age_get = c(inter_c_moving/sum(inter_c_moving))) %>% 
    mutate(cum_get_prop = abs(cumsum(age_get)-0.5)) %>% 
    arrange(cum_get_prop) %>% 
    head(1) %>% 
    pull(age)

  
  df <- df %>% rbind(c(out_name, out_person_yll, out_community_yll, 
                       out_person_ymm, out_community_ymm,out_get, 
                       median_onset_age, quartile_onset_age))
  
  temp_inter_results <- out_new %>% 
    select(age, yll, ymm) %>% 
    mutate(out_name = out_name)
  
  df_inter_ymm_yll <- df_inter_ymm_yll %>% 
    rbind(temp_inter_results)
  
  
  if(nrow(df) %% 10 == 0){print(df)}
}

colnames(df) <- c("type", "mean_person_yll", "community_1000_yll",
                  "mean_person_ymm", "community_1000_ymm", 
                  "proportion_acquiring", "median_onset_age", 
                  "quartile_onset_age")
df <- df %>% mutate(mean_person_yll = as.numeric(mean_person_yll),
                    community_1000_yll = as.numeric(community_1000_yll),
                    mean_person_ymm = as.numeric(mean_person_ymm),
                    community_1000_ymm = as.numeric(community_1000_ymm),
                    proportion_acquiring = as.numeric(proportion_acquiring),
                    median_onset_age = as.numeric(median_onset_age), 
                    quartile_onset_age = as.numeric(quartile_onset_age))

df

df_clean <-  df %>% 
  mutate(type  = stringr::str_replace_all(type, "Diabetes", "")) %>% 
  mutate(type  = stringr::str_replace_all(type, "-", "")) %>% 
  mutate(type  = stringr::str_replace_all(type, "_", " ")) %>% 
  filter(!type %in% c("Multiple Sclerosis", "Sarcoidosis", "Autism", 
                      "Cystic Fibrosis", "Sickle Cell Disease")) %>% 
  mutate(type = if_else(type == "Parkinsons Disease", 
                        "Parkinson's Disease", 
                        type))

df_top10s <- rbind(
  df_clean %>% arrange(desc(mean_person_yll)) %>% 
    select(type, mean_person_yll) %>% 
    head(10) %>% 
    mutate(variable = "mean_person_yll") %>% 
    rename(value = mean_person_yll)
  ,
  df_clean %>% arrange(desc(community_1000_yll)) %>% 
    select(type, community_1000_yll) %>% 
    head(10) %>% 
    mutate(variable = "community_1000_yll") %>% 
    rename(value = community_1000_yll) 
  ,
  df_clean %>% arrange(desc(mean_person_ymm)) %>% 
    select(type, mean_person_ymm) %>% 
    head(10) %>% 
    mutate(variable = "mean_person_ymm") %>% 
    rename(value = mean_person_ymm)   
  ,
  df_clean %>% arrange(desc(community_1000_ymm)) %>% 
    select(type, community_1000_ymm) %>% 
    head(10) %>% 
    mutate(variable = "community_1000_ymm") %>% 
    rename(value = community_1000_ymm)        
  ,
  df_clean %>% arrange(desc(proportion_acquiring)) %>% 
    select(type, proportion_acquiring) %>% 
    head(10) %>% 
    mutate(variable = "proportion_acquiring") %>% 
    rename(value = proportion_acquiring) 
)

setwd(base)
write_csv(df_top10s, paste0("top10_",result_set, ".csv"))
