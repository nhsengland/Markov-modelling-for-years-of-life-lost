get_age_transition_matrix <- function(df, i){
  x <- df %>% 
    filter(age == i) %>% 
    pivot_wider(names_from = to, values_from = prop) %>% 
    select(-age, -from) %>% 
    as.matrix()
  
  # create 1 year transition probabilities from month transitions
  y <- x %^% 12 #infix operator from expm library
  
  rownames(y) <- colnames(x)
  return(y)
}



matrixproduct <- function(a, b){
  a %*% b #infix operator from expm library
}



calculate_transition_probabilities <- function(df){
  df %>% 
    group_by(age, from) %>% 
    mutate(row_total = sum(n)) %>% 
    mutate(prop = n/row_total) %>% 
    as.data.frame() %>% 
    select(-n, -row_total)
} 




revert_df_trans <- function(array){
  rname <- rownames(array)
  out <- array %>% 
    as.data.frame() %>% 
    mutate(from = rname) %>% 
    relocate(from)
  rownames(out) <- NULL
  return(out)
}


generate_all_table <- function(df){
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
    pivot_longer(cols = `Neither morbidity`:`Death`, names_to = "to", values_to = "prop")
}
