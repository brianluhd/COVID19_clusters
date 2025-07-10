# Define a function to call svymean and unweighted count
getSummary <- function(varformula, byformula, design){
  #' Calculates mean and se
  #' 
  #' @details Adapted from NHANES analysis [example code](https://wwwn.cdc.gov/nchs/data/Tutorials/Code/DB303_R.r).
  
  # Get mean, stderr, and unweighted sample size
  c <- svyby(varformula, byformula, design, unwtd.count) 
  p <- svyby(varformula, byformula, design, svymean) 
  outSum <- left_join(select(c,-se), p) 
  return(outSum)
}


# Define a function to call svymean and unweighted count
get_cont_summary <- function(varformula, byformula, design){
  #' Calculates mean and 95% CI
  #' 
  #' Reformats output of getSummary to mean and 95% CI for continuous features.
  
  # calculate mean and ci
  c <- svyby(varformula, byformula, design, unwtd.count ) 
  p <- svyby(varformula, byformula, design, svymean ) 
  ci <- confint(p)
  outSum <- left_join(select(c,-se), p)  %>% bind_cols(ci)
  
  var_name <- as.character(varformula)[2]
  
  df <- outSum %>%
    as_tibble() %>%
    select(-starts_with("se."), -counts) %>% # remove se columns
    mutate(across(where(is.numeric), ~ sprintf("%.1f", .x))) %>%
    mutate("{var_name}" := as.character(glue::glue("{.[[var_name]]} ({.[['2.5 %']]}, {.[['97.5 %']]})"))) %>%
    select(-se, -`2.5 %`, -`97.5 %`) %>%
    pivot_longer(cols = all_of(var_name), names_to = "Features") %>% # transpose
    pivot_wider(names_from = as.character(byformula)[2], 
                values_from = value)
  
  return(df)
}

get_cat_summary <- function(varformula, byformula, design){
  #' Calculates percentages
  #' 
  #' Reformats output of getSummary to percentages for categorical features.
  
  outSum <- getSummary(varformula, byformula, design)
  
  # reformate output to tibble
  var_name <- as.character(varformula)[2]
  
  df <- outSum %>%
    as_tibble() %>%
    select(-starts_with("se.", ignore.case = F), -counts) %>% # remove se columns
    mutate(across(starts_with(var_name), ~ (round(.x*100)))) %>% # convert proportion to percent and round to integer
    pivot_longer(cols = where(is.numeric), names_to = "Features") %>% # transpose
    pivot_wider(names_from = as.character(byformula)[2], 
                values_from = value) %>%
    mutate(Value = str_remove(Features, var_name),
           Features = str_extract(Features, var_name)) %>% # clean up names
    relocate(Value, .after = Features) %>%
    mutate(across(where(is.numeric), as.character))
  
  return(df)
}

get_survey_or <- function(response, dependent_variables, design){
  #' Calculates logistic regression
  #' 
  #' Calculates logistic regression based on given response, dependent variables, and survey design.
  
  or_formula <- reformulate(dependent_variables, response = response)
  
  or_regression <- svyglm(or_formula, design = design, family = quasibinomial)
  
  regression_ci <- exp(confint(or_regression))
  
  out <- as_tibble(summary(or_regression)$coefficients, rownames = "variables") %>%
    mutate(OR = exp(Estimate)) %>%
    bind_cols(regression_ci)
  
  return(out)
}