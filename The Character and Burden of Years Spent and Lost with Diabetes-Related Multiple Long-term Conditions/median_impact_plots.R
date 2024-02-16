library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)
library(readr)


result_list <-  c("results_female", "results_male", "results_full")

for(result_set in result_list){
  
  base <- "./Male_Female_Markov_revised_transition_data/"
  
  setwd(paste0(base, result_set))
  
  tbl <-
    list.files(pattern = "res_") %>% 
    map_df(~read_csv(.))
  
  out <- filter(tbl, grepl("Diabetes", type)) %>% 
    filter(!grepl("Frailty", type)) %>% 
    mutate(type  = stringr::str_replace_all(type, "Diabetes", "")) %>% 
    mutate(type  = stringr::str_replace_all(type, "-", "")) %>% 
    mutate(type  = stringr::str_replace_all(type, "_", " ")) %>% 
    filter(!type %in% c("Multiple Sclerosis", "Sarcoidosis", "Autism", 
                        "Cystic Fibrosis", "Sickle Cell Disease")) %>% 
    mutate(type = if_else(type == "Parkinsons Disease", 
                          "Parkinson's Disease", 
                          type))
  
  
  
  get_median_age_data <- function(group){
    filter(out, type == group) %>% 
      arrange(age) %>% 
      mutate(cum_get_prop = abs(cumsum(age_get)-0.5)) %>% 
      arrange(cum_get_prop) %>% 
      head(1)
  }
  
  out_median_data <- map_df(unique(out$type), get_median_age_data) 
  
  
  
  p_median <- out_median_data %>% select(age, life_exp_mean, ymm, yll, type) %>% 
    mutate(type = as.factor(type)) %>% 
    mutate(type = forcats::fct_reorder(type, age)) %>% 
    mutate(life_exp_w_mm = age + ymm) %>% 
    mutate(life_exp_wo_mm = age + ymm + yll) %>% 
    select(age, life_exp_w_mm, life_exp_wo_mm, type) %>% 
    pivot_longer(life_exp_w_mm:life_exp_wo_mm, names_to = "metric") %>% 
    ggplot(aes(y = type, xmin = age, xmax = value, colour = metric)) + 
    geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
    scale_x_continuous(name = "Onset Age of Condition Pair", limits = c(58,92), 
                       minor_breaks = seq(58, 92, 2)) + 
    scale_y_discrete(name = "Diabetes + ") +
    scale_colour_discrete(name = NULL, 
                          labels = c("Life years in those with conditions", 
                                     "Life years in general population")) +
    theme_minimal() +
    theme(legend.position = "top")
  
  
  
  setwd(base)
  ggsave(paste0("plot_median_", result_set, ".png"), plot = p_median, 
         width = 10*0.85, height = 7*0.85)
  saveRDS(p_median, paste0("plot_RDS_median_", result_set, ".RDS"))
}
