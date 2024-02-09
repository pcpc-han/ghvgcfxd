cie_bw <- function(outcome, exposure, contvar, catvar, data) {
    covariates <- c(contvar, catvar)  # All covariates
    variables_of_interest <- c(exposure, outcome, contvar, catvar)
    dt_subset <- data[complete.cases(data[variables_of_interest]), ]
    print(paste("Number of obs:", nrow(dt_subset)))
    
    # initial full model
    full_formula <- paste(outcome, "~", exposure, 
                          "+", paste(contvar, collapse = "+"), 
                          "+", paste(sapply(catvar, function(x) paste0("factor(", x, ")")), collapse = "+"))
    full_model <- lm(as.formula(full_formula), data = dt_subset)
    initial_adj_coef <- coef(full_model)[exposure]
    
    model_list <- list(full_model = full_model)
    
    while (length(covariates) > 0) { 
        
        CIE <- data.frame(Variable = character(), Change = numeric(), stringsAsFactors = FALSE)
        
        for (cov in covariates) {
            
            # Create a temporary formula without the current covariate
            temp_covariates <- setdiff(covariates, cov)
            temp_formula_parts <- c(exposure, temp_covariates[!(temp_covariates %in% catvar)], 
                                    sapply(temp_covariates[temp_covariates %in% catvar], function(x) paste0("factor(", x, ")")))
            temp_formula_str <- paste(outcome, "~", paste(temp_formula_parts, collapse = " + ")) %>% print
            
            # Fit the model
            temp_model <- lm(as.formula(temp_formula_str), data = dt_subset)
            tidy(temp_model, conf.int = TRUE)[2, ] %>% kable %>% print()
            temp_adj_coef <- coef(temp_model)[exposure]
            change_in_estimate <- 100 * (temp_adj_coef - initial_adj_coef) / initial_adj_coef
            
            # Store CIE
            CIE <- rbind(CIE, data.frame(Variable = cov, Change = change_in_estimate))
        }
        
        # Find covariate with least CIE to remove
        CIE %>% print()
        if (nrow(CIE) == 0) break 
        min_change <- CIE[which.min(abs(CIE$Change)), ]
        covariates <- setdiff(covariates, min_change$Variable)
        
        # Update model list
        updated_formula_parts <- c(exposure, covariates[!(covariates %in% catvar)], 
                                   sapply(covariates[covariates %in% catvar], function(x) paste0("factor(", x, ")")))
        updated_formula_str <- paste(outcome, "~", paste(updated_formula_parts, collapse = " + "))
        updated_model <- lm(as.formula(updated_formula_str), data = dt_subset)
        initial_adj_coef <- coef(updated_model)[exposure]
        model_list[[paste("Model without", min_change$Variable)]] <- updated_model
        
        cat("Removed:", min_change$Variable, "with a change of", min_change$Change, "% in the exposure coefficient.\n")
    }
}
