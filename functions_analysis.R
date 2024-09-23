continuous_feature <- function(l){
  #' Continuous features
  #' 
  #' Check if a feature is continuous or categorical.
  #' 
  #' @param l Vector of values for the feature to be determined
  
  if (class(l) == "numeric"){
    if ((max(l, na.rm = T) - min(l, na.rm = T) != 1) &
        (max(l, na.rm = T) - min(l, na.rm = T) != 0)){
      cont = T
    } else {
      cont = F
    }
  } else {
    cont = F
  }
  return(cont)
}

find_units <- function(df, features){
  #' Find units
  #' 
  #' Find the unit of each feature.
  #' 
  #' @param df Dataframe containing the data
  #' @param features list of features for which to find units
  #' 
  #' The unit of a feature must be labeled with the suffix "_units"
  #' after the name of the feature.
  
  # construct string for regex match 
  regex_match <- paste0("^", features, "_units$") %>% paste(collapse = "|")
  
  # regex match to find columns containing units for each feature
  unit_list <- df %>% 
    select(matches(regex_match)) %>% 
    unique() %>% as.list() # extract as list
  unit_list <- lapply(unit_list, na.omit) # remove where unit is NA from no available data
  unit_list <- lapply(unit_list, unique) # remove duplicates
  unit_list <- unlist(unit_list) # collapse into single vector
  
  # remove units suffix to get back name of each feature
  names(unit_list) <- str_replace_all(names(unit_list), "_units", "")
  
  # convert named list to tibble
  unit_tibble <- as_tibble(unit_list, rownames = "Features")
  unit_tibble <- unit_tibble %>% rename(Value = value)
  
  # add back in features with no found units
  all_features <- tibble(Features = features) %>%
    left_join(unit_tibble, by = join_by(Features))
  
  
  return(all_features)
}

cont_stats <- function(df, features, category, covariates = NA){
  #' Univariate statistics for continuous features
  #' 
  #' Perform Student's T test or type 2 ANOVA as appropriate for continuous features.
  #' P-values are rounded to the nearest 4 decimal points.
  #' 
  #' @param df Dataframe containing the data
  #' @param features list of features on which to perform the univariate analysis
  #' @param category name of the feature denoting the grouping of the features for analysis
  #' 
  #' P-values less than 0.0001 are denoted as <0.0001.
  
  p <- c()
  adjusted_p <- c()
  vif <- c()
  for (f in features){
    # unadjusted
    ## construct formula with dynamic variables
    fo <- reformulate(category, f)
    
    if (length(unique(df[[category]])) > 2){ # perform ANOVA if more than 2 groups
      anova_result <- Anova(lm(fo, data = df), type = 2)
      new_p <- anova_result$`Pr(>F)`[1]
    } else { # else perform student's t test
      ttest_result <- do.call("t.test", list(fo, substitute(df)))
      new_p <- ttest_result$`p.value`
    }
    
    ## structure p values
    if (new_p < 0.0001){
      new_p <- "<0.0001"
    } else {
      new_p <- sprintf("%.4f", new_p)
    }
    
    ## append to p values vector
    p <- append(p, new_p)
    
    # adjusted
    if (!all(is.na(covariates))){
      ## construct formula with dynamic variables, exclude covariate if same as exposure
      terms <- c(category, covariates[covariates != f]) # isolate terms for regression model to count number of terms later
      fo <- reformulate(terms, f)
      anova_result <- Anova(lm(fo, data = df), type = 2)
      new_adjusted_p <- anova_result$`Pr(>F)`[1]
      
      ## structure p values
      if (new_adjusted_p < 0.0001){
        new_adjusted_p <- "<0.0001"
      } else {
        new_adjusted_p <- sprintf("%.4f", new_adjusted_p)
      }
      
      ## append to p values vector
      adjusted_p <- append(adjusted_p, new_adjusted_p)
      
      # calculate max vif
      v <- vif(lm(fo, data = df)) # VIF output is different depending on degrees of freedom of covariates
      if (is.matrix(v)){ # if df >1 for any covariates, outputs matrix
        model_vif <- max(v[,"GVIF"])
      } else { # if df = 1 for all, outputs vector
        model_vif <- max(v)
      }
      vif <- append(vif, max(model_vif))
    }
  }
  
  # paste category name to keep track of p values
  p_value_name <- paste(stringr::str_to_title(category), "p-value")
  adjusted_p_value_name <- paste(stringr::str_to_title(category), "adjusted p-value")
  
  # construct output tibble
  stats_out <- tibble(Features = features, 
                      {{p_value_name}} := p,
                      {{adjusted_p_value_name}} := adjusted_p,
                      `Max VIF` = vif)
  
  return(stats_out)
}

reorder_columns <- function(df){
  #' Reorder columns
  #' 
  #' Reorder the columns as Feature, Summary, then Value.
  #' If Race is present, then White is listed first, followed by
  #' Black or African American, then Other.
  #' If Gender if present, Male is listed first.
  #' 
  df <- df %>% relocate(Features) %>%
    relocate(Summary, .after = Features) %>%
    relocate(Value, .after = Summary)
  
  if (all(c("White", "Other") %in% colnames(df))){
    df <- df %>% relocate(Other, .after = White )
  }
  if (all(c("White", "Black or African American") %in% colnames(df))){
    df <- df %>% relocate(`Black or African American`, .after = White )
  }
  if (all(c("Male", "Female") %in% colnames(df))){
    df <- df %>% relocate(Female, .after = Male )
  }
  if (all(c("SAID", "SIDD", "SIRD", "MOD", "MARD") %in% colnames(df))){
    df <- df %>% relocate(MARD, .after = SAID) %>%
      relocate(MOD, .after = SAID) %>%
      relocate(SIRD, .after = SAID) %>%
      relocate(SIDD, .after = SAID)
  }
  if ("Other" %in% colnames(df)){
    df <- df %>% relocate(any_of(c("SAID", "SIDD", "SIRD", "MOD", "MARD")), .after = Other)
  }
  
  return(df)
}

conf_mean <- function(x, level = 0.95, string = T, digits = 1){
  model <- lm(x~1)
  conf_int <- confint(model, level = level)
  if (string){
    conf_int <- sprintf(paste0("%.", digits, "f, ", "%.", digits, "f"), conf_int[1], conf_int[2])
  }
  return(conf_int)
}

calc_cont_summary <- function(df, features, categories, covariates = covariates, distribution_stat = "sd", digits = 1, ...){
  #' Descriptive statistics for continuous features
  #' 
  #' Calculate mean and standard deviation for continuous features.
  #' Values are rounded to the nearest 1 decimal point.
  #' 
  #' @param df Dataframe containing the data
  #' @param features list of features on which to calculate mean and standard deviation
  #' @param categories name of the feature denoting the grouping of the continuous features
  #' 
  #' Values are returned as mean (sd). Units for each feature were identified
  #' with `find_units`. Univariate statistics were calculated with `cont_stats`.
  
  output_summary <- tibble(Features = features,
                           Summary = case_when(distribution_stat == "sd" ~ "mean (SD)",
                                               distribution_stat == "ci" ~ "mean (95% CI)"))
  # find units
  unit_tibble <- find_units(df, features)
  
  # add units to tibble
  output_summary <- output_summary %>% left_join(unit_tibble, by = join_by(Features))
  for (c in categories){
    # calculate mean and standard deviations
    cat_mean <- df %>% group_by(.data[[c]]) %>% summarize_at(features, mean, na.rm = T)
    
    if (distribution_stat == "sd"){
      cat_sd <- df %>% group_by(.data[[c]]) %>% summarize_at(features, sd, na.rm = T)
    } else if (distribution_stat == "ci"){
      cat_conf <- df %>% group_by(.data[[c]]) %>% summarize_at(features, conf_mean, string = T, digits = digits)
    }
    
    
    # extract columns as list and format float as single digit strings as mean (sd)
    col_list <- list()
    col_list[[c]] <- cat_mean[[c]]
    for (i in 2:ncol(cat_mean)){
      col <- colnames(cat_mean)[i]
      
      if (distribution_stat == "sd"){
        col_list[[col]] <- sprintf(paste0("%.", digits, "f ", "(%.", digits, "f)"), cat_mean[[col]], cat_sd[[col]])
      } else if (distribution_stat == "ci"){
        col_list[[col]] <- sprintf(paste0("%.", digits, "f (%s)"), cat_mean[[col]], cat_conf[[col]])
      }
    }
    
    # convert list of mean (sd) to tibble
    cat_summary <- data.frame(col_list, row.names = cat_mean[[c]]) %>% # construct as dataframe to use rownames as names in each category
      t() %>% # flip rows and columns to get desired format, with features as rows
      as_tibble(rownames = "Features") # convert to tibble
    
    # perform ttest/anova
    p_tibble <- cont_stats(df, features, c, covariates)
    
    # append p values
    cat_summary <- cat_summary %>% left_join(p_tibble, by = join_by(Features))
    
    # join with outoput summary
    output_summary <- output_summary %>% left_join(cat_summary, by = join_by(Features))
    
    # calculate overall mean and distribution
    cont_mean <- unlist(lapply(df[features], mean, na.rm = T)) 
    if (distribution_stat == "sd"){
      cont_sd <- unlist(lapply(df[features], sd, na.rm = T)) 
      cont_overall <- sprintf(paste0("%.", digits, "f ", "(%.", digits, "f)"), cont_mean, cont_sd)
    } else if (distribution_stat == "ci"){
      cont_ci <- unlist(lapply(df[features], conf_mean, digits = digits))
      cont_overall <- sprintf(paste0("%.", digits, "f (%s)"), cont_mean, cont_ci)
    }
    
    # append to output
    output_summary <- output_summary %>%
      mutate(Overall = cont_overall) %>%
      relocate(Overall, .after = Value)
  }
  
  
  return(output_summary)
}

calc_cat_summary <- function(df, features, categories, covariates = NA){
  #' Descriptive and univariate statistics for categorical features
  #' 
  #' Calculate count and percentage for categorical features.
  #' Values are rounded to the nearest integer. Chi-squared test of independence 
  #' or Fisher's exact test were performed as appropriate for each feature.
  #' 
  #' @param df Dataframe containing the data
  #' @param features list of features on which to calculate mean and standard deviation
  #' @param categories list of features denoting the grouping of the categorical features
  #' 
  #' Values are returned as count (%). Units for each feature were identified
  #' with `find_units`.
  
  output_summary <- tibble(Features = features[1])
  
  for (f in features){
    f_summary <- tibble(Features = f,
                        Value = as.character(unique(df[[f]]))) %>%
      arrange(desc(Value)) %>% # show yes first
      drop_na() # remove na rows
    
    # calculate counts and percents
    counts <- table(df[f])
    counts <- counts[sort(names(counts), decreasing = T)] # show yes first)
    percents <- counts/sum(counts)*100
    
    f_summary <- f_summary %>% 
      mutate(Overall = sprintf("%.0f (%.0f)", counts, percents)) %>%
      relocate(Overall, .after = Value)
    for (c in categories){
      # create contingency table
      count_table <- df %>% 
        drop_na(all_of(f)) %>% # remove NA feature values
        count(across(all_of(f)), across(all_of(c))) %>% # count produces long form data
        pivot_wider(names_from = all_of(c), values_from = n) %>% # pivot to wide data for contingency table
        replace(is.na(.), 0)
      # calculate column-wise proportions
      prop <- prop.table(as.matrix(count_table[-1]), margin = 2) %>% as_tibble()
      prop <- prop * 100
      
      # append proportion to counts
      col_list <- list()
      for (i in 2:ncol(count_table)){
        col <- colnames(count_table)[i]
        col_list[[col]] <- sprintf("%.0f (%.0f)", count_table[[col]], prop[[col]])
      }
      
      cat_summary <- as_tibble(col_list) %>%
        mutate(Value = as.character(pull(count_table[1])))
      
      # unadjusted
      ## perform chi-squared test
      x <- chisq.test(count_table[-1])
      
      ## perform fisher exact if any expected counts are less than 5
      if (any(x$expected < 5)){
        x = fisher.test(count_table[-1])
        warning(paste0("Expected count less than 5 in ", f, " by ", c, ", Fisher's exact test was used instead."))
      }
      
      ## structure p values
      p <- x$p.value
      if (p < 0.0001){
        p <- "<0.0001"
      } else {
        p <- sprintf("%.4f", p)
      }
      
      ## append p-value
      p_value_name <- paste(stringr::str_to_title(c), "p-value")
      cat_summary <- cat_summary %>% mutate({{p_value_name}} := p)
      
      # adjusted
      if (!all(is.na(covariates))){
        ## determine if covariate is categorical or continuous
        mask <- c()
        for (co in covariates){
          mask <- append(mask, continuous_feature(df[[co]]))
        }
        cont_covariates <- covariates[mask]
        cat_covariates <- covariates[!mask]
        
        ## append factor name to categorical covariates
        cat_factor <- paste(cat_covariates, "factor", sep = "_")
        
        ## create factor for binary categories if doesn't exist
        for (i in 1:length(cat_covariates)){
          if (!(cat_factor[i] %in% colnames(df))){
            df[[cat_factor[i]]] <- factor(df[[cat_covariates[i]]])
            if (is.numeric(df[[cat_covariates[i]]])){
              df[[cat_factor[i]]] <- relevel(df[[cat_factor[i]]], ref = "0")
            }
          }
        }
  
        ## combine back cat and continuous covariates
        factor_covariates <- c(cat_factor, cont_covariates)
        
        ## create factor for binary response if doesn't exist
        response_factor <- paste(f, "factor", sep = "_")
        if (!(response_factor %in% colnames(df))){
          df[[response_factor]] <- factor(df[[f]])
          if (is.numeric(df[[f]])){
            df[[response_factor]] <- relevel(df[[response_factor]], ref = "0")
          }
        }
        
        ## create factor for binary category if doesn't exist
        cat_factor <- paste(c, "factor", sep = "_")
        if (!(cat_factor %in% colnames(df))){
          df[[cat_factor]] <- factor(df[[c]])
          if (is.numeric(df[[c]])){
            df[[cat_factor]] <- relevel(df[[cat_factor]], ref = "0")
          }
        }
        
        # declare p-value name
        adjusted_p_value_name <- paste(stringr::str_to_title(c), "adjusted p-value")
        
        # logistic regression if exposure is binary
        if (length(unique(df[[f]])) == 2){
          terms <- c(cat_factor, factor_covariates[factor_covariates != response_factor]) # isolate terms for regression model to count number of terms later
          fo <- reformulate(terms, response_factor)
          logit <- summary(glm(fo, data = df, family = "binomial"))
          logit_stats <- logit$coefficients %>% as_tibble(rownames = "Value") %>%
            rename({{adjusted_p_value_name}} := `Pr(>|z|)`) %>%
            select(Value, {{adjusted_p_value_name}}) %>%
            filter(Value == c | grepl(cat_factor, .$Value)) %>%
            mutate(Value = levels(df[[response_factor]])[2],
                   {{adjusted_p_value_name}} := ifelse(.[[adjusted_p_value_name]] < 0.0001, "<0.0001", sprintf("%.4f", .[[adjusted_p_value_name]]))) 
          
          # add max vif
          v = vif(glm(fo, data = df, family = "binomial"))
          if (is.vector(v)){ # if df = 1 for all, outputs vector
            logit_stats <- logit_stats %>% mutate(`Max VIF` = max(v))
            
          } else if (is.matrix(v)) { # if df >1 for any covariates, outputs matrix
            logit_stats <- logit_stats %>% mutate(`Max VIF` = max(v[,"GVIF"]))
          } else {
            logit_stats <- logit_stats %>% mutate(`Max VIF` = NA)
          }
          # append p-value
          cat_summary <- cat_summary %>% left_join(logit_stats, by = join_by(Value))

        } else if (length(unique(df[[f]])) > 2){ # if exposure is multinomial
          fo <- reformulate(c(cat_factor, factor_covariates[factor_covariates != response_factor]), response_factor)
          logit <- nnet::multinom(fo, data = df)
          logit_stats <- broom::tidy(logit) %>% 
            filter(grepl(c, .$term)) %>% select(`y.level`, `p.value`) %>% 
            rename(Value = `y.level`, {{adjusted_p_value_name}} := `p.value`) %>%
            mutate({{adjusted_p_value_name}} := ifelse(.[[adjusted_p_value_name]] < 0.0001, "<0.0001", sprintf("%.4f", .[[adjusted_p_value_name]])),
                   `Max VIF` = NA) # placeholder for VIF, since can't simply use vif function on multinomial regression
          
          # append p-value
          cat_summary <- cat_summary %>% left_join(logit_stats, by = join_by(Value))
        } else {
          adjusted_p <- NA
        }
      }
      
      # join to f_summary
      f_summary <- left_join(f_summary, cat_summary, by = join_by(Value))
    
    }
    
    # join each feature summary to output
    if (f == features[1]){
      output_summary <- left_join(output_summary, f_summary, by = join_by(Features))
    } else {
      output_summary <- rbind(output_summary, f_summary)
    }
      
  }
  
  # add Summary column to indicate counts, convert numerical binary to yes/no
  output_summary <- output_summary %>% mutate(Summary = "count (%)",
                                              Value = ifelse(.$Value == "1", "Yes", 
                                                             ifelse(.$Value == "0", "No", .$Value)))
  
  return(output_summary)
}


get_feature_summary <- function(df, features, categories, covariates = NA, ...){
  #' Aggregate descriptive and univariate statistics for features
  #' 
  #' Calculate mean and standard deviation for continuous features rounded to 1 decimal place.
  #' Calculate count and percentage for categorical features rounded to the nearest integer.
  #' 
  #' @param df Dataframe containing the data
  #' @param features list of features on which to calculate mean and standard deviation
  #' @param categories list of features denoting the grouping of the features to aggregate
  #' 
  #' Continuous and categorical features are determined using `continuous feature`.
  #' Statistics for continuous features were calculated using `calc_cont_summary`.
  #' Statistics for categorical features were calculated using `calc_cat_summary`.
  #' Aggregated columns were reordered using `reorder_columns`.
  #' 
  
  # identify and separate continuous vs categorical features
  mask <- c()
  for (f in features){
    mask <- append(mask, continuous_feature(df[[f]]))
  }
  cont_features <- features[mask]
  cat_features <- features[!mask]
  
  
  cat_summary <- NULL
  cont_summary <- NULL
  # calculate summaries
  if (length(cont_features) > 0){

    
    cont_summary <- calc_cont_summary(df, cont_features, categories, covariates, ...)
    
  }
  if (length(cat_features) > 0){
    cat_summary <- calc_cat_summary(df, cat_features, categories, covariates)
  }
  
  # combine continuous and categorical summaries
  if (!is.null(cat_summary) & !is.null(cont_summary)){
    combined <- full_join(cont_summary, cat_summary)
  } else if (is.null(cat_summary)){
    combined <- cont_summary
  } else if (is.null(cont_summary)){
    combined <- cat_summary
  }
  
  # round VIF to 2 digits if present
  if ("Max VIF" %in% colnames(combined)){
    combined <- combined %>% mutate(`Max VIF` = round(`Max VIF`, 2))
  }
  
  # re-order rows to original feature order as supplied
  combined$Features <- factor(combined$Features, levels = features)
  combined <- combined %>% arrange(Features)
  
  # re-order columns
  combined <- reorder_columns(combined)
  
  return(combined)
}

logit_table <- function(df, response, covariates){
  #' Logistic regression summary table
  #' 
  #' Produces summary table of logistic regression.
  #' 
  #' @param df dataframe containing the input data
  #' @param response response variable
  #' @param covariates covariates for regression
  #' 
  #' The response variable is printed as columns and covariates are printed as rows. 
  #' Categorical coviarates must be stored as factors and named as "[covariate]_factor" in `df`.
  
  # construct summary table
  logit_summary <- get_feature_summary(df, features = covariates, categories = response)
  
  # determine if covariate is categorical or continuous
  mask <- c()
  for (f in covariates){
    mask <- append(mask, continuous_feature(df[[f]]))
  }
  cont_covariates <- covariates[mask]
  cat_covariates <- covariates[!mask]
  
  # append factor name to categorical covariates
  cat_factor <- paste(cat_covariates, "factor", sep = "_")
  
  # create factor for binary categories if doesn't exist
  for (i in 1:length(cat_covariates)){
    if (!(cat_factor[i] %in% colnames(df))){
      df[[cat_factor[i]]] <- factor(df[[cat_covariates[i]]])
      if (is.numeric(df[[cat_covariates[i]]])){
        df[[cat_factor[i]]] <- relevel(df[[cat_factor[i]]], ref = "0")
      }
    }
  }
  
  # construct regression formula
  fo <- reformulate(c(cat_factor, cont_covariates), response)
  
  model <- do.call("glm", list(fo, substitute(df), family = "binomial"))
  
  # extract p-values
  p <- summary(model)$coefficients %>% 
    as_tibble(rownames = "comparison") %>% 
    select(comparison, "Pr(>|z|)") %>% 
    rename("Logit p-value" = "Pr(>|z|)") %>%
    mutate("Logit p-value" = ifelse(.$"Logit p-value" < 0.0001, "<0.0001", sprintf("%0.4f", .$"Logit p-value")))
  
  # extract coefficients
  co <- exp(coef(model))
  
  # extract confidence interval
  conf <- exp(confint.default(model))
  
  # construct odds ratio tibble
  odds <- tibble(comparison = names(co),
                 OR = sprintf("%0.2f (%0.2f, %0.2f)", co, conf[,1], conf[,2])) %>%
    left_join(p, by = join_by(comparison)) 
  
  # separate continuous and categorical covariates since need to join on different columns
  cat_odds <- odds %>% filter(grepl("_factor", .$comparison)) %>% 
    mutate(Features = str_replace_all(.$comparison, "(.*)_factor.*", "\\1"), # extract feature names
           comparison = str_replace_all(.$comparison, ".*_factor(.)", "\\1")) %>% # remove "factor" from comparison name
    mutate(comparison = replace(comparison, comparison == "1", "Yes")) %>% # replace 1 with "Yes"
    right_join(select(logit_summary, Features, Summary, Value), by = c("comparison" = "Value", "Features" = "Features")) %>%
    rename(Value = comparison) %>% # rename to match logit_summary colnames
    drop_na(OR)
  cont_odds <- odds %>% filter(!grepl("_factor", .$comparison)) %>%
    right_join(select(logit_summary, Features, Summary, Value), by = c("comparison" = "Features")) %>%
    rename(Features = comparison) %>% # rename to match logit_summary colnames
    drop_na(OR)
  
  odds <- rbind(cat_odds, cont_odds)
  
  # join odds to summary
  logit_summary <- left_join(logit_summary, odds, by = join_by(Features, Summary, Value))
  
  # rename p-value columns
  p_value_name <- paste(stringr::str_to_title(response), "p-value")
  logit_summary <- logit_summary %>% rename(`Unadjusted p-value` = all_of(p_value_name),
                                            `OR (95% CI)` = OR)
  # calculate VIF if more than 1 covariates
  if (length(covariates) > 1){
    vif_tibble <- vif(model) %>% as_tibble(rownames = "Features") %>%
      mutate(Features = str_replace_all(.$Features, "(.*)_factor", "\\1"))
    
    if ("value" %in% colnames(vif_tibble)){ # from if df = 1 for all covariates, vif returns list
      vif_tibble <- vif_tibble %>% rename(VIF = value)
    }
    
    # add VIF
    logit_summary <- left_join(logit_summary, vif_tibble, by = join_by(Features)) 
  }
  
  
  return(logit_summary)
}

anova_vif <- function(model, type = 2){
  #' Anova table with VIF
  #' 
  #' Appends VIF with Anova table
  #' 
  #' @param model lm model
  #' 
  #' VIF and Anova are calculated using the `car` package.
  
  anova_table <- as_tibble(Anova(model, type = type), rownames = "Covariates")
  
  vif_list <- as_tibble(vif(model), rownames = "Covariates") %>%
    rename(VIF = value)
  
  anova_table <- anova_table %>% left_join(vif_list, by = join_by(Covariates))
  
  return(anova_table)
}

feature_cor <- function(df, response, covariates, 
                        p_cutoff = 0.05, 
                        background_inc = c(), 
                        similar_measure = c(), 
                        demo = c(), 
                        collider = c(), 
                        duplicated = c(), 
                        strict = T){
  
  #' Correlation table
  #' 
  #' Calculates correlations for covariate selection in regression model
  #' 
  #' @param df dataframe
  #' @param response response variable
  #' @param covariates covariates to test
  #' @param p_cutff correlation p-value cutoff to determine significance
  #' @param background_inc list of covariates to include based on background knowledge
  #' @param similar_measure list of covariates that are similar to/same as the response variable and should be excluded
  #' @param demo list of demographic covariates to be included
  #' @param collider list of covariates that are colliders to be excluded
  #' @param duplicated list of covariates that are similar to/same as another included covariate and should be excluded
  #' @param strict boolean value indicating how to deal with significantly correlated covariates that are not include by
  #'  background knowlege. If True, only covariates specified by background knowlege will be included. If False, significantly 
  #'  correlated covariates that are not excluded by background knowlege will be included.
  
  p_value <- c()
  pearson <- c()
  covariates <- covariates[covariates != response]
  for (f in covariates){
    cor_res <- cor.test(df[[response]], df[[f]], use = "complete.obs")
    pearson <- append(pearson, cor_res$estimate)
    p_value <- append(p_value, cor_res$p.value)
  }
  
  out <- tibble(feature = covariates,
                pearson = pearson,
                "p-value" = p_value)
  
  out <- out %>% 
    # always include demographic covariates
    mutate(include = ifelse(.$feature %in% demo, 1, 0),
           reason = ifelse(.$feature %in% demo, "demographic", NA)) %>% 
    # include covariates based on background knowledge 
    mutate(include = ifelse(.$feature %in% background_inc, 1, .$include),
           reason = ifelse(.$feature %in% background_inc, "background knowledge", .$reason)) %>%     
    # exclude covariates with no correlation
    mutate(include = ifelse((.$`p-value` >= p_cutoff) & (.$include == 0), 0, .$include),
           reason = ifelse((.$`p-value` >= p_cutoff) & (.$include == 0), "no correlation", .$reason)) %>% 
    # exclude covariates that are the same/similar measure as the response
    mutate(include = ifelse(.$feature %in% similar_measure, 0, .$include),
           reason = ifelse(.$feature %in% similar_measure, "similar measure as response", .$reason)) %>% 
    # exclude covariates that are the same/similar measure as another covariate
    mutate(include = ifelse(.$feature %in% duplicated, 0, .$include),
           reason = ifelse(.$feature %in% duplicated, "similar measure as another covariate", .$reason)) %>% 
    # exclude colliders
    mutate(include = ifelse(.$feature %in% collider, 0, .$include),
           reason = ifelse(.$feature %in% collider, "collider", .$reason))
  
  if (strict){
    out <- out %>% mutate(include = ifelse(is.na(.$reason), 0, .$include),
                          reason =  ifelse(is.na(.$reason), "no background knowledge reason to include", .$reason))
  } else {
    out <- out %>% mutate(include = ifelse(is.na(.$reason), 1, .$include),
                          reason =  ifelse(is.na(.$reason), "no background knowledge reason to exclude", .$reason))
  }
  
  return(out)
}
