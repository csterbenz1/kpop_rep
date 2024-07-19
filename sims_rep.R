### Packages
library(MASS)
library(tidyverse)
library(survey)
#devtools::install_github("csterbenz1/KBAL", ref = "cat_kernel")
library(kbal)
library(parallel)
library(knitr)
library(glmnet)

##### The following scrip will run the simulations presented in section 5.2
#### NB: the following script is computationally intensive and timely to run
#### NB2: all simulation in Appendix D are produced in sims_res_rep.Rmd using the file saved from this script; users who do not with to rerun the timely simulations but reproduce the figures and tables from the simulation results may simply run the results .Rmd
#The file proceeds as follows
#0) Preliminary parameters to control arguments in the simulations
#1) Fit lasso selection model whose probabilities are later used to take Bernoulli draws to create a sample in each iteration
#2) Specify outcome model (identical in variables to selection model), add random noise, and adjust coefficients in an automated fashion to produce probability-like outcomes for the P(vote= Democrat)
#3) Create population targets and functions for computing standard errors according to various estimators
#4) run sims in series (note that running in parallel leads to memory issues and results in actually slower run time that running in series)


###### SET PARAMS  ###############
set.seed(9345876)
#change bath to data as needed
path_data = paste0(getwd(), "/generated_data/")
save_path = paste0(getwd(), "/weights/")
#loads original prepared data by authors
#if data prepared manually with data_prep_rep.R then should be named just cces_prepared.rds
cces_filename = "cces_prepared_orig.rds"

options(dplyr.print_max = 1e9)

####### Parameters for kpop results:
#ebal tolerance and max iterations for kpop
tolerance = 1e-4
maxit = 500
#adjust these both for runtime if needed
#increment in svd(K) dims in biasbound minimization search
increment = 5
#specify min/max svd(K) dims in biasbound min search
min_num_dims = NULL
max_num_dims = NULL
#to run just nonkpop methods
eval_kpop = TRUE
#to run kpop+mf
kpop_constraints = TRUE
SAVE = TRUE #save .Rdata results?
#can adjust accordingly to machine for adequate number of sims
#NB: this is the central determinant of runtime, can decrease to save time
nsims = 10

##### Central Params to adjust for Model Specification
#sd(y)*noise; 1-> r^2 = .5; sqrt(2) -> r^2 = .33; 1/2*sqrt(2) -> r^2 = .66;
noise = 1 
#adjusts sample size by dividing p(S) by scalar pS_denom (i.e. pS = plogis(XBeta)/pS_denom)
pS_denom = 60

###################### Formulas ################
formula_rake_demos_noeduc <- ~recode_age_bucket + recode_female + recode_race +
    recode_region + recode_pid_3way
#updated to include 6 way edu
formula_rake_demos_weduc <- ~recode_age_bucket + recode_female +
    recode_race + recode_region + recode_educ + recode_pid_3way

formula_rake_all_vars <- ~recode_age_bucket + recode_female +
    recode_race + recode_region + recode_pid_3way + recode_educ +
    recode_income_5way + recode_relig_6way + recode_born + recode_attndch_4way

create_targets <- function (target_design, target_formula) {
    target_mf <- model.frame(target_formula, model.frame(target_design))
    target_mm <- model.matrix(target_formula, target_mf)
    wts <- weights(target_design)
    
    return(colSums(target_mm * wts) / sum(wts))
}

### Post-stratification function
## For now assumes that strata variable is already created and in
## the data set and called "strata"‚
postStrat <- function(survey, pop_counts, pop_w_col, strata_pass, warn = T) {
    survey_counts <- survey %>%
        group_by(!!as.symbol(strata_pass)) %>%
        summarize(n = n()) %>%
        ungroup() %>%
        mutate(w_survey = n / sum(n))
    
    pop_counts <- pop_counts %>%
        rename(w_pop = matches(pop_w_col))
    
    if(warn == T & nrow(survey_counts) !=  nrow(pop_counts)) {
        missing_strat = pop_counts[! (( pop_counts[, strata_pass]%>% pull()) %in% (survey_counts[, strata_pass]%>% pull() )), strata_pass]
        warning(paste("Strata in Pop not found in Sample. Dropping", 
                      sum(pop_counts[(pop_counts[, strata_pass] %>% pull()) %in% 
                                         (missing_strat %>% pull()),"n" ]), 
                      "empty cells\n"), immediate.  =T )
    } 
    post_strat <- pop_counts %>%
        left_join(survey_counts, by = strata_pass) %>%
        filter(!is.na(w_survey)) %>%
        ## Normalizes back to 1 after dropping
        ## empty cells
        mutate(w_pop = w_pop * 1/sum(w_pop),
               w = w_pop / w_survey) %>%
        dplyr::select(!!as.symbol(strata_pass), w)
    
    survey <- survey %>%
        left_join(post_strat)
    
    return(survey)
}

#ensure that sample/selection model is not so severe that sample has no units in strata important to the selection model
check_sample <- function(sample, selection_model) {
    check = model.matrix(selection_model, data = sample)
    check = colSums(check)
    fail = check[which(check ==0)]
    fail_bin = length(check[which(check ==0)]) > 0 
    check_prop = check/nrow(sample)
    
    return(list(samp_prop = check_prop,
                fail  = fail, 
                fail_bin = fail_bin,
                counts = check))
}

#ensure that sample/selection model is not so severe that sample has no units in strata important to the selection and outcome model
check_sample_outcome <- function(sample, pop = cces, selection_model, interaction_cols, interaction_cols_2 = NULL) {
    
    vars = all.vars(selection_model)
    var = NULL
    counts = NULL
    prop = NULL
    u_outcome = NULL
    
    #uninteracted variables
    run_counts <- function(data) {
        for(i in 1:length(vars)) {
            t = data %>% group_by_at(vars[i]) %>%
                summarise(n = n(), 
                          avg_outcome = mean(outcome)) %>%
                mutate(prop = round(n/nrow(sample),4))
            var = c(var, as.character(t[,1] %>% pull()))
            counts = c(counts, t$n)
            prop = c(prop,  t$prop)
            u_outcome = c(u_outcome, t$avg_outcome)
        }
        #interactions
        t = suppressMessages(data %>% group_by_at(interaction_cols) %>% 
                                 summarise(n = n(),
                                           avg_outcome = mean(outcome)) %>%
                                 mutate(prop = round(n/nrow(sample), 4)))
        interaction = apply(t, 1,  function(r) paste(r[1],r[2], collapse = "_"))
        counts = data.frame(var  = c(var, interaction),
                            n = c(counts, t$n), 
                            prop = c(prop, t$prop),
                            avg_outcome = c(u_outcome, t$avg_outcome))
        
        if(!is.null(interaction_cols_2)) {
            t2 = suppressMessages(data %>% group_by_at(interaction_cols_2) %>% 
                                      summarise(n = n(),
                                                avg_outcome = mean(outcome)) %>%
                                      mutate(prop = round(n/nrow(sample), 4)))
            interaction = apply(t2, 1,  function(r) paste(r[1],r[2], collapse = "_"))
            append = cbind(data.frame(var = interaction), t2[, - c(1,2)])
            counts = rbind(counts, append)
        }
    }
    samp_counts = run_counts(sample)
    pop_counts = run_counts(pop)
    
    if(nrow(pop_counts) > nrow(samp_counts)) {
        samp_counts = rbind(samp_counts, data.frame(var = pop_counts$var[!(pop_counts$var %in% samp_counts$var)],
                                                    n = rep(0,sum(!(pop_counts$var %in% samp_counts$var))),
                                                    prop = rep(0,sum(!(pop_counts$var %in% samp_counts$var))),
                                                    avg_outcome = rep(-99,sum(!(pop_counts$var %in% samp_counts$var)))))
    }
    
    
    fail = sum(samp_counts$n == 0)
    bad = sum(samp_counts$prop <= 0.05)
    bad_strata = data.frame(strata = as.character(samp_counts$var[samp_counts$prop <= 0.05]), prop = samp_counts$prop[samp_counts$prop <= 0.05])
    v_bad = sum(samp_counts$prop <= 0.01)
    v_bad_strata = data.frame(strata = as.character(samp_counts$var[samp_counts$prop <= 0.01]), prop = samp_counts$prop[samp_counts$prop <= 0.01])
    counts$var[samp_counts$prop <= 0.01]
    
    samp_counts = samp_counts %>% mutate(leq_5pp = as.numeric(prop <= 0.05),
                                         leq_1pp = as.numeric(prop <= 0.01))
    
    
    return(list(samp_counts = samp_counts,
                fail = fail, 
                bad = bad, 
                v_bad = v_bad, 
                bad_strata = bad_strata,
                v_bad_strata = v_bad_strata))
}

#checks the outcome is not beyond the support (probability range 0-1)
check_outcome <- function(outcome) {
    beyond_support = (min(outcome) <0 | max(outcome) > 1)
    return(beyond_support)
}

#automated search for outcome coefficients that produce outcomes in probability range for all units
bound_outcome <- function(outcome, coefs, increment = 1, increment_intercept = .01, noise, cces_expanded, silent = T) {
    denom = 10
    fail = check_outcome(outcome)
    while(fail) {
        denom = denom + increment
        if(!silent) { cat(denom, ": ") }
        coefs_use = coefs/denom
        outcome = cces_expanded %*% coefs_use
        
        outcome = outcome + rnorm(nrow(cces_expanded), mean = 0, sd = sd(outcome)*noise)
        summary(outcome)
        if(max(outcome) <=1 & min(outcome) <0 ) {
            coefs[1] = coefs[1] + increment_intercept
            if(!silent) { cat("\nmoving intercept up", coefs[1],  "\n") }
            denom = denom - increment
        }
        if(max(outcome) >1 & min(outcome) >=0 ) {
            coefs[1] = coefs[1] - increment_intercept
            if(!silent) { cat("\nmoving intercept down", coefs[1],  "\n") }
            denom = denom - increment
        }
        fail = check_outcome(outcome)
        if(!silent) { cat(round(min(outcome),2), round(max(outcome),2),  "\n") }
    }
    if(!silent) { cat(paste("Min denom:", denom)) }
    return(list(outcome = outcome, coefs = coefs_use, denom = denom))
    
}


############# Load Data #####################
#these data have been cleaned already see app_modeled for how it was done
### Load Target Data
cces <- readRDS(paste0(path_data,cces_filename))
cces$recode_age_bucket = as.character(cces$recode_age_bucket)
cces$recode_age_3way= as.character(cces$recode_age_3way)

##################### LASSO: Selection #############################

selection_model = as.formula(~recode_pid_3way + recode_female + recode_age_bucket 
                             + recode_educ_3way 
                             + recode_race
                             + recode_born 
                             + recode_born:recode_age_bucket
                             + recode_pid_3way:recode_age_bucket
)
inter = c("recode_pid_3way", "recode_age_bucket")
inter_2 = c("recode_born", "recode_age_bucket")
cces_expanded = model.matrix(selection_model, data = cces)
coefs = matrix(NA,nrow = ncol(cces_expanded), ncol =1 )
rownames(coefs) = colnames(cces_expanded)
coefs
#(p(S)) for negative bias select non dem voters
coefs[,1] = c(-2, #intercept -5 w race
              2, #selection of indep pos
              2, #selection of R pos
              .5, #male
              .15, #36-50,
              .2, #51-64,
              .2, #65+,
              .7, #college
              -1 , #post-grad
              .5,#hispanic
              .3,#other
              .7,#white
              2, #bornagain
              1,#bornagain x 36-50
              1.5, #bornagain x 51-64
              2, #bornagain x 65+
              .3,#ind x 36-50
              .5, #rep x 36-50,
              1, #ind x 51-64,
              1, #rep x 51-64,
              -.2, #ind x 65+
              2 #rep x 65+
)

xbeta = cces_expanded %*% coefs
p_include = plogis(xbeta)
p_include = p_include/pS_denom
sum(p_include)


#################### DESIGN OUTCOME MODEL ##################
#p(D)
#start perfectly negatively correlated
coefs_outcome = -coefs
#scale out as a starting point before adjusting to be in probability range
coefs_outcome = coefs_outcome*1.5
cor(p_include, cces_expanded %*% coefs_outcome)
#adjust starting point of intercept so we're in ball part of prob range
#so the following adjustment of coefs doesn't take forever
coefs_outcome[1] = 25
cor(p_include, cces_expanded %*% coefs_outcome)

cat(paste("Adding sd(outcome)*",round(noise, 3), "\n")) 
#reset seed so noise added is directly controllable/replicable without running in order above
set.seed(1383904)
bound = bound_outcome(outcome = cces_expanded %*% coefs_outcome,
                      coefs = coefs_outcome,
                      cces_expanded = cces_expanded,
                      noise = noise, silent = F)
coefs_outcome = bound$coefs
xbeta_outcome = bound$outcome
if(check_outcome(xbeta_outcome)) {
    warning("Outcome beyond prob support for some units when noise is added",
            immediate. = T)
}

cat(paste("Mean outcome w/noise is", round(mean(xbeta_outcome)*100,3), "\n"))
cat(paste("Range of outcome w/noise is\n"))
cat(paste(summary(xbeta_outcome), "\n"))
s = summary(lm(update(selection_model, xbeta_outcome ~ .),data = cces))
R2_outcome = s$adj.r.squared
cat(paste("R^2 outcome is", round(s$adj.r.squared,3), "\n"))
cat(paste("Mean scaled outcome (target) is", round(mean(xbeta_outcome)*100,3)))
cat(paste("\nCorr of sampling prob and outcome ", round(cor(xbeta_outcome, p_include),3)))
cces$outcome = xbeta_outcome


######### Make STRATA variable in CCES ############
#post-stratification gets true selection model
formula_ps <- selection_model
cces <- bind_cols(cces, cces %>%
                      unite("strata", all.vars(formula_ps), remove = FALSE) %>%
                      dplyr::select(strata))


#################### Targets ###################
cces_svy <- suppressWarnings(svydesign(ids = ~1, data = cces))
margin_sim = svymean(~outcome, cces_svy)[1]* 100
margin_sim
targets_rake_demos_noeduc <- create_targets(cces_svy,
                                            formula_rake_demos_noeduc)
targets_rake_demos_weduc <- create_targets(cces_svy, formula_rake_demos_weduc)
targets_rake_all_vars <- create_targets(cces_svy,
                                        formula_rake_all_vars)
targets_demo_truth <- create_targets(cces_svy, selection_model)


## Make table of Population Counts for post-stratification for manual ps function
cces_counts <- cces %>%
    group_by(strata) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    mutate(w = n / sum(n, na.rm = TRUE))

########################### PREP SIMS ##################################

est_mean <- function(outcome, design) {
    svymean(as.formula(paste0("~", outcome)), design, na.rm = TRUE)[1]
}

#########################################
############## variance calc ###########
## Variance functions
var_fixed <- function(Y, weights, pop_size) {
    ## note: needs weights that sum to population total
    if(round(sum(weights)) != pop_size) { weights = weights*pop_size/sum(weights)}
    return(Hmisc::wtd.var(Y, weights))
}

## kott (14) (under poisson)
var_quasi <- function(weights, residuals, pop_size) {
    #moving from kott 14 w sum w =N to weights that sum to 1 + var of total to var of mean:
    #sum^n (w_i^2 - 1/N_pop w_i)e_i^2 
    return(sum((weights^2 - (weights / pop_size))*residuals^2))
}

## kott (15) linearization
var_linear <- function(weights, residuals, sample_size) {
    #moving from kott 15 w sum w =N to weights that sum to 1 + var of total to var of mean:
    # n/(n-1) sum^n (w_i*e_i)^2 - (1/n-1) [sum(w_i*e_i)^2] *using this nonapprox version (deriv in my notes)
    # approx = sum^n (w_i*e_i)^2 - (1/n) sum(w_i*e_i)^2
    n = sample_size
    return((n/(n-1))*sum((weights * residuals)^2) - 1/(n-1) * sum(weights * residuals)^2)
}

## chad
var_chad <- function(weights, residuals) {
    return(sum(weights^2 * residuals^2))
}

## calculate all variances
calc_SEs <- function(Y, residuals, pop_size, weights, sample_size) {
    if(round(sum(weights)) != 1 ) {
        weights = weights/sum(weights)
    }
    return(data.frame(SE_fixed = sqrt(var_fixed(Y, weights, pop_size) / length(Y)),
                      SE_quasi = sqrt(var_quasi(weights, residuals, pop_size)),
                      SE_linear = sqrt(var_linear(weights, residuals, sample_size)),
                      SE_chad = sqrt(var_chad(weights, residuals))))
}

########################### RUN SIMS ##################################

#double check params are as expected:sum(p_include)
SAVE
eval_kpop
kpop_constraints
#save some memory by deleting some larger objects for the parallel run
pryr::mem_used()
rm(bound, cces_expanded, s, xbeta, xbeta_outcome)
pryr::mem_used()
set.seed(1293908934)
system.time({
    for(nsim in 1:nsims)  {
        #
        if(nsim == 1) {
            sims = NULL
            SEs = NULL
            residuals_out = NULL
            dropped_cells_out = NULL
            sample_out = NULL 
            weights = NULL
            out = NULL
        }
        
        system.time({
            cat(paste("==============================================  SIM:",nsim, "/", nsims, "=", nsim/nsims,
                      "============================================== \n"))
            
            sample <- rbinom(nrow(cces), 1, p_include)
            survey_sim <- cces[sample == 1, ]
            
            survey_design <- suppressWarnings(svydesign(ids = ~1, data = survey_sim))
            
            ########################### check sample ##########################
            #checks to see if we got a sample that has no units in a strata in the selection model 
            #this would indicate our selection probabilities are too small
            check_s = check_sample(survey_sim, selection_model)
            bad_sample = check_s$fail_bin
            #tracks the counts of units in each strata in the selection model
            check_2 = check_sample_outcome(survey_sim, pop = cces, 
                                           selection_model, interaction_cols = inter, interaction_cols_2 = inter_2)
            check_nums = c(leq_5pp = check_2$bad,
                           leq_1pp = check_2$v_bad, 
                           fail = check_2$fail)
            s = survey_sim %>% group_by(recode_pid_3way,recode_female, recode_age_bucket,
                                        recode_educ_3way) %>% count() %>%
                mutate(n_s = round(n/nrow(survey_sim), 3))
            c = cces %>% group_by(recode_pid_3way, recode_female, recode_age_bucket,
                                  recode_educ_3way) %>% count() %>%
                mutate(n_c = round(n/nrow(cces), 3))
            count = nrow(c) - nrow(s)
            
            ############################################
            ## Unweighted estimate
            ############################################
            
            unweighted <- est_mean("outcome", survey_design)
            ############################################
            ## Sample size
            ############################################\
            n <- sum(sample) 
            
            ############################################
            ## Raking on demographics (no education)
            ############################################
            rake_demos_noeduc_svyd <- try(calibrate(design = survey_design,
                                                    formula = formula_rake_demos_noeduc,
                                                    population = targets_rake_demos_noeduc,
                                                    calfun = "raking"), silent = T)
            
            rake_demos_noeduc <- tryCatch(est_mean("outcome", rake_demos_noeduc_svyd), error = function(e) NA)
            
            #SEs
            if(class(rake_demos_noeduc_svyd)[1] == "try-error") {
                rake_demos_noeduc_se <- data.frame(SE_fixed = NA,
                                                   SE_quasi = NA, 
                                                   SE_linear = NA, 
                                                   SE_chad = NA) 
                names(rake_demos_noeduc_se) = paste0("rake_demos_noeduc_", names(rake_demos_noeduc_se))
                rake_demos_noeduc_se_SVY = data.frame(rake_demos_noeduc_se_SVY  = NA)
            } else {
                residuals = residuals(lm(update(formula_rake_demos_noeduc, outcome ~ .), 
                                         data = rake_demos_noeduc_svyd$variables))
                
                res_rake_demos_noeduc = data.frame(min = min(residuals), 
                                                   perc_25 = quantile(residuals, .25), 
                                                   mean = mean(residuals),
                                                   perc_75 = quantile(residuals, .75),
                                                   var = var(residuals))
                
                rake_demos_noeduc_se <- calc_SEs(Y = rake_demos_noeduc_svyd$variables$outcome, 
                                                 residuals = residuals, 
                                                 pop_size = nrow(cces), 
                                                 sample_size = sum(sample),
                                                 weights = weights(rake_demos_noeduc_svyd))
                names(rake_demos_noeduc_se) = paste0("rake_demos_noeduc_", names(rake_demos_noeduc_se))
                
                
                rake_demos_noeduc_se_SVY = data.frame(rake_demos_noeduc_se_SVY = data.frame(svymean(~outcome, 
                                                                                                    rake_demos_noeduc_svyd, 
                                                                                                    na.rm = TRUE))[1,2])
                
            }
            ############################################
            #### Raking on demographics (with education)
            ############################################
            
            rake_demos_weduc_svyd <- try(calibrate(design = survey_design,
                                                   formula = formula_rake_demos_weduc,
                                                   population = targets_rake_demos_weduc,
                                                   calfun = "raking"), silent = T)
            
            rake_demos_weduc <- tryCatch(est_mean("outcome", rake_demos_weduc_svyd), 
                                         error = function(e) NA)
            
            
            #SEs
            if(class(rake_demos_weduc_svyd)[1] == "try-error") {
                rake_demos_weduc_se <- data.frame(SE_fixed = NA,
                                                  SE_quasi = NA, 
                                                  SE_linear = NA, 
                                                  SE_chad = NA) 
                names(rake_demos_weduc_se) = paste0("rake_demos_weduc_", names(rake_demos_weduc_se))
                
            } else {
                residuals = residuals(lm(update(formula_rake_demos_weduc, outcome ~ .), 
                                         data = rake_demos_weduc_svyd$variables))
                res_rake_demos_wedu = data.frame(min = min(residuals), 
                                                 perc_25 = quantile(residuals, .25), 
                                                 mean = mean(residuals),
                                                 perc_75 = quantile(residuals, .75),
                                                 var = var(residuals))
                rake_demos_weduc_se <- calc_SEs(Y = rake_demos_weduc_svyd$variables$outcome, 
                                                residuals = residuals, 
                                                pop_size = nrow(cces),
                                                sample_size = sum(sample),
                                                weights = weights(rake_demos_weduc_svyd))
                names(rake_demos_weduc_se) = paste0("rake_demos_weduc_", names(rake_demos_weduc_se))
                
            }
            
            ############################################
            #### Raking on everything
            ############################################
            rake_all_svyd <- try(calibrate(design = survey_design,
                                           formula = formula_rake_all_vars,
                                           population = targets_rake_all_vars,
                                           calfun = "raking"))
            
            rake_all <- tryCatch(est_mean("outcome", rake_all_svyd), 
                                 error = function(e) NA)
            #SEs
            if(class(rake_all_svyd)[1] == "try-error") {
                rake_all_se <- data.frame(SE_fixed = NA,
                                          SE_quasi = NA, 
                                          SE_linear = NA, 
                                          SE_chad = NA) 
                names(rake_all_se) = paste0("rake_all_", names(rake_all_se))
                
            } else {
                residuals = residuals(lm(update(formula_rake_all_vars, outcome ~ .), 
                                         data = rake_all_svyd$variables))
                res_rake_all = data.frame(min = min(residuals), 
                                          perc_25 = quantile(residuals, .25), 
                                          mean = mean(residuals),
                                          perc_75 = quantile(residuals, .75),
                                          var = var(residuals))
                rake_all_se <- calc_SEs(Y = rake_all_svyd$variables$outcome, 
                                        residuals = residuals, 
                                        pop_size = nrow(cces), 
                                        sample_size = sum(sample),
                                        weights = weights(rake_all_svyd))
                names(rake_all_se) = paste0("rake_all_", names(rake_all_se))
                
            }
            
            ############################################
            ## Post-stratification: Truth
            ############################################
            
            #track empty cells:
            #this subsets cces strata to only those in survey_sim
            missing_strata <- unique(cces$strata)[!(unique(cces$strata) %in%
                                                        unique(survey_sim$strata))]
            cat(round(length(missing_strata)/ length(unique(cces$strata)),3),
                "% cces original strata missing from sample, ",
                " and", cces %>% filter(strata %in% missing_strata) %>% summarise(n()) %>% pull(), "/", nrow(cces), "units\n" )
            dropped_cells = cces %>% filter(strata %in% missing_strata) %>% group_by(strata) %>% count()
            #dropped_cells = data.frame(sum = sum(dropped_cells$n), strata = paste(dropped_cells$strata, collapse = " | "))
            dropped_cells_pass = sum(dropped_cells$n)
            
            post_stratification_svyd = svydesign(~1, data = postStrat(survey_sim, 
                                                                      cces_counts, "w", 
                                                                      strata_pass = "strata", 
                                                                      warn = F),
                                                 weights = ~w)
            
            post_stratification <- est_mean("outcome", post_stratification_svyd)
            
            #SEs can use svy pacakge
            post_stratification_se <- data.frame(post_strat_SE_svy = data.frame(svymean(~outcome, post_stratification_svyd,
                                                                                        na.rm = TRUE))[1,2])
            
            
            ############################################
            #### Raking on true model
            ############################################
            #error catching just a stop gap to prevent from breaking sims run if truth doesnt converge
            rake_truth_svyd <- try(calibrate(design = survey_design,
                                             formula = selection_model,
                                             population = targets_demo_truth,
                                             maxit = 100,
                                             calfun = "raking"), silent = T)
            
            if(class(rake_truth_svyd)[1] == "try-error") {
                rake_truth_svyd <- try(calibrate(design = survey_design,
                                                 formula = selection_model,
                                                 population = targets_demo_truth,
                                                 calfun = "raking",
                                                 maxit = 100,
                                                 epsilon = .009), silent = T)
            }
            
            if(class(rake_truth_svyd)[1] == "try-error") {
                
                rake_truth_se <- data.frame(SE_fixed = NA,
                                            SE_quasi = NA, 
                                            SE_linear = NA, 
                                            SE_chad = NA) 
                names(rake_truth_se) = paste0("rake_truth_", names(rake_truth_se))
                
            } else {
                lambdas <- 10^seq(3, -2, by = -.1)
                x <- model.matrix(update(selection_model, outcome ~ .),
                                  data = rake_truth_svyd$variables)[, -1]
                fit <- glmnet(x, 
                              rake_truth_svyd$variables$outcome, alpha = 0, lambda = lambdas)
                cv_fit <- cv.glmnet(x, rake_truth_svyd$variables$outcome, alpha = 0, lambda = lambdas)
                opt_lambda <- cv_fit$lambda.min
                fit <- cv_fit$glmnet.fit
                
                residuals = rake_truth_svyd$variables$outcome - predict(fit, s = opt_lambda, newx = x)
                res_rake_truth = data.frame(min = min(residuals), 
                                            perc_25 = quantile(residuals, .25), 
                                            mean = mean(residuals),
                                            perc_75 = quantile(residuals, .75),
                                            var = var(residuals))
                rake_truth_se <- tryCatch(calc_SEs(Y = rake_truth_svyd$variables$outcome,
                                                   residuals = residuals,
                                                   pop_size = nrow(cces),
                                                   sample_size = sum(sample),
                                                   weights = weights(rake_truth_svyd)), error = function(e) NA)
                names(rake_truth_se) = paste0("rake_truth_", names(rake_truth_se))
            }
            
            rake_truth <- tryCatch(est_mean("outcome", rake_truth_svyd), 
                                   error = function(e) NA)
            
            #HT and Hayek
            p_sample <- as.matrix((p_include[sample==1]))
            
            #HT
            ht_truth = sum((cces[sample ==1, "outcome"]/p_sample))/nrow(cces) 
            #Hayek
            hayek_truth = sum((cces[sample ==1, "outcome"]/p_sample))/sum(1/p_sample)
            
            ############################################
            ## Kpop: Categorical Data + b = argmax V(K)
            ############################################
            if(eval_kpop) {
                # Select the covariates for use in Kbal: updated cat data no cont age
                #one-hot coded for cat kernel
                kbal_data <- bind_rows(survey_sim %>% dplyr::select(recode_age_bucket,
                                                                    recode_female,
                                                                    recode_race,
                                                                    recode_region,
                                                                    recode_pid_3way,
                                                                    recode_educ,
                                                                    recode_income_5way,
                                                                    recode_relig_6way,
                                                                    recode_born,
                                                                    recode_attndch_4way),
                                       cces %>% dplyr::select(recode_age_bucket,
                                                              recode_female,
                                                              recode_race,
                                                              recode_region,
                                                              recode_pid_3way,
                                                              recode_educ,
                                                              recode_income_5way,
                                                              recode_relig_6way,
                                                              recode_born,
                                                              recode_attndch_4way))
                
                kbal_data_sampled <- c(rep(1, nrow(survey_sim)), rep(0, nrow(cces)))
                
                #### DEFAULT ######
                #rstudio_para_cat(c("nsim: ", nsim, " DEFAULT"))
                cat(paste("nsim:", nsim, "DEFAULT", "\n"))
                kbal_est <- kbal(allx=kbal_data,
                                 sampled = kbal_data_sampled,
                                 cat_data = TRUE,
                                 incrementby = increment,
                                 meanfirst = FALSE,
                                 ebal.tol = tolerance,
                                 ebal.maxit = maxit,
                                 minnumdims = min_num_dims,
                                 maxnumdims = max_num_dims,
                                 sampledinpop = FALSE,
                                 fullSVD = TRUE)
                
                kpop_svyd <- svydesign(~1, data = survey_sim,
                                       weights = kbal_est$w[kbal_data_sampled ==1])
                
                kpop <- est_mean("outcome", kpop_svyd)
                b_kpop = kbal_est$b
                #save memory by saving only the svd to re use
                svdK = kbal_est$svdK 
                numdims = kbal_est$numdims
                biasbound_r = kbal_est$biasbound_ratio
                biasbound = kbal_est$biasbound_opt
                
                ##### Kpop SEs
                kpop <- tryCatch(est_mean("outcome", kpop_svyd), error = function(e) NA)
                
                x <- as.matrix(data.frame(kbal_dims = kbal_est$svdK$v[, 1:kbal_est$numdims]))
                cv_fit <- cv.glmnet(x, kpop_svyd$variables$outcome, alpha = 0)
                lambda_pass = cv_fit$lambda.1se
                
                residuals = kpop_svyd$variables$outcome - predict(cv_fit$glmnet.fit,
                                                                  s = lambda_pass, newx = x)
                res_kpop = data.frame(min = min(residuals), 
                                      perc_25 = quantile(residuals, .25), 
                                      mean = mean(residuals),
                                      perc_75 = quantile(residuals, .75),
                                      var = var(residuals))
                kpop_se <- tryCatch(calc_SEs(Y = kpop_svyd$variables$outcome,
                                             residuals = residuals,
                                             pop_size = nrow(cces),
                                             sample_size = sum(sample),
                                             weights = weights(kpop_svyd)), error = function(e) NA)
                
                if(length(kpop_se) == 1) {
                    kpop_se <- data.frame(SE_fixed = NA, 
                                          SE_quasi = NA, 
                                          SE_linear = NA, 
                                          SE_chad = NA)
                }
                names(kpop_se) = tryCatch(paste0("kpop_", names(kpop_se)), error = function(e) NA)
                
                #CONVERGED
                dist_record = data.frame(t(kbal_est$dist_record))
                min_converged = dist_record[which.min(dist_record[dist_record$Ebal.Convergence ==1,"BiasBound"]), "Dims"]
                
                rm(kbal_est, residuals, x, cv_fit)
                
                #CONVERGED
                #### CONVG ####
                cat(paste("nsim:", nsim, " CONV", "\n"))
                if(is.null(min_converged) | length(min_converged) ==0) {
                    kpop_svyd_conv <- "dn converge"
                    kpop_conv <- "dn converge"
                    
                    numdims_conv = "dn converge"
                    biasbound_r_conv = "dn converge"
                    biasbound_conv = "dn converge"
                    kpop_conv_se = data.frame(SE_fixed = NA, 
                                              SE_quasi = NA, 
                                              SE_linear = NA, 
                                              SE_chad = NA)
                    res_kpop_conv = data.frame(min = NA, 
                                               perc_25 = NA, 
                                               mean =NA,
                                               perc_75 = NA,
                                               var = NA)
                    
                } else {
                    kbal_est_conv <- kbal(allx=kbal_data,
                                          K.svd = svdK,
                                          sampled = kbal_data_sampled,
                                          numdims = min_converged,
                                          ebal.tol = tolerance,
                                          ebal.maxit = maxit,
                                          minnumdims = min_num_dims,
                                          maxnumdims = max_num_dims,
                                          scale_data = FALSE,
                                          drop_MC = FALSE,
                                          incrementby = increment,
                                          meanfirst = FALSE,
                                          sampledinpop = FALSE,
                                          ebal.convergence = TRUE)
                    kpop_svyd_conv <- svydesign(~1, data = survey_sim,
                                                weights = kbal_est_conv$w[kbal_data_sampled ==1])
                    kpop_conv <- est_mean("outcome", kpop_svyd_conv)
                    
                    numdims_conv = kbal_est_conv$numdims
                    biasbound_r_conv = kbal_est_conv$biasbound_ratio
                    biasbound_conv = kbal_est_conv$biasbound_opt
                    
                    #SEs
                    x <- as.matrix(data.frame(kbal_dims = kbal_est_conv$svdK$v[, 1:kbal_est_conv$numdims]))
                    cv_fit <- cv.glmnet(x, kpop_svyd_conv$variables$outcome, alpha = 0)
                    fit <- cv_fit$glmnet.fit
                    lambda_pass = cv_fit$lambda.1se
                    residuals = kpop_svyd_conv$variables$outcome - predict(cv_fit$glmnet.fit, 
                                                                           s = lambda_pass, 
                                                                           newx = x)
                    res_kpop_conv = data.frame(min = min(residuals), 
                                               perc_25 = quantile(residuals, .25), 
                                               mean = mean(residuals),
                                               perc_75 = quantile(residuals, .75),
                                               var = var(residuals))
                    kpop_conv_se <- tryCatch(calc_SEs(Y = kpop_svyd_conv$variables$outcome,
                                                      residuals = residuals,
                                                      pop_size = nrow(cces),
                                                      sample_size = sum(sample),
                                                      weights = weights(kpop_svyd_conv)), error = function(e) NA)
                    if(length(kpop_conv_se) == 1) {
                        kpop_conv_se <- data.frame(SE_fixed = NA, 
                                                   SE_quasi = NA, 
                                                   SE_linear = NA, 
                                                   SE_chad = NA)
                    }
                    names(kpop_conv_se) = tryCatch(paste0("kpop_conv_", names(kpop_conv_se)), error = function(e) NA)
                    #KRLS SEs are exactly the same for coverged
                    rm(kbal_est_conv, residuals, x, cv_fit) 
                }
                
                
                ####### MF #######
                #rstudio_para_cat(c("nsim: ", nsim, " + aMF"))
                cat(paste("nsim:", nsim, "aMEANFIRST", "\n"))
                
                if(kpop_constraints) {
                    #########demos constraint method:
                    cat(paste("nsim:", nsim, "CONSTR (Demos)", "\n"))
                    kbal_demos_est <- kbal(K.svd = svdK,
                                           allx=kbal_data,
                                           cat_data = TRUE,
                                           sampled = kbal_data_sampled,
                                           ebal.tol = tolerance,
                                           ebal.maxit = maxit,
                                           minnumdims = min_num_dims,
                                           maxnumdims = max_num_dims,
                                           scale_data = FALSE,
                                           drop_MC = FALSE,
                                           incrementby = increment,
                                           meanfirst = TRUE,
                                           mf_columns = all.vars(formula_rake_demos_noeduc),
                                           sampledinpop = FALSE)
                    
                    kpop_demos_svyd <- svydesign(~1, data = survey_sim, 
                                                 weights = kbal_demos_est$w[kbal_data_sampled ==1])
                    
                    kpop_demos <- est_mean("outcome", kpop_demos_svyd)
                    
                    numdims_demos = kbal_demos_est$numdims
                    mf_appended_dims_demos = kbal_demos_est$meanfirst_dims
                    if(is.null(numdims_demos)) {
                        numdims_demos = c(NA) 
                        kpop_demos_se <- data.frame(SE_fixed = NA, 
                                                    SE_quasi = NA, 
                                                    SE_linear = NA, 
                                                    SE_chad = NA)
                    } else {
                        V <-  data.frame(kbal_dims = kbal_demos_est$svdK$v[, c(1:kbal_demos_est$numdims)])
                        X <- as.matrix(cbind(kbal_demos_est$appended_constraint_cols[kbal_data_sampled==1, ], V))
                        
                        cv_fit <- cv.glmnet(X, kpop_demos_svyd$variables$outcome, alpha = 0, 
                                            penalty.factor = c(rep(0, ncol(kbal_demos_est$appended_constraint_cols)), rep(1, kbal_demos_est$numdims)))
                        
                        lambda_pass = cv_fit$lambda.1se
                        residuals =  kpop_demos_svyd$variables$outcome - predict(cv_fit$glmnet.fit,
                                                                                 s = lambda_pass, 
                                                                                 newx = X)
                        res_kpop_demos = data.frame(min = min(residuals), 
                                                    perc_25 = quantile(residuals, .25), 
                                                    mean = mean(residuals),
                                                    perc_75 = quantile(residuals, .75),
                                                    var = var(residuals))
                        
                        kpop_demos_se <- tryCatch(calc_SEs(Y = kpop_demos_svyd$variables$outcome,
                                                           residuals = residuals,
                                                           pop_size = nrow(cces),
                                                           sample_size = sum(sample),
                                                           weights = weights(kpop_demos_svyd)), 
                                                  error = function(e) NA)
                        if(length(kpop_demos_se) == 1) {
                            kpop_demos_se <- data.frame(SE_fixed = NA, 
                                                        SE_quasi = NA, 
                                                        SE_linear = NA, 
                                                        SE_chad = NA)
                        }
                        names(kpop_demos_se) = tryCatch(paste0("kpop_demos_", names(kpop_demos_se)),
                                                        error = function(e) NA)
                    }
                    biasbound_r_demos = kbal_demos_est$biasbound_ratio
                    biasbound_demos = kbal_demos_est$biasbound_opt
                    
                    rm(kbal_demos_est, residuals, X, V, cv_fit)
                    
                    
                    #########demos + educ constraint method:
                    cat(paste("nsim:", nsim, "CONSTR (D+Edu)", "\n"))
                    kbal_demos_wedu_est <- kbal(K.svd = svdK,
                                                allx=kbal_data,
                                                cat_data = TRUE,
                                                sampled = kbal_data_sampled,
                                                ebal.tol = tolerance,
                                                ebal.maxit = maxit,
                                                minnumdims = min_num_dims,
                                                maxnumdims = max_num_dims,
                                                scale_data = FALSE,
                                                drop_MC = FALSE,
                                                incrementby = increment,
                                                meanfirst = TRUE,
                                                mf_columns = all.vars(formula_rake_demos_weduc),
                                                sampledinpop = FALSE)
                    kpop_demos_wedu_svyd <- svydesign(~1, data = survey_sim, 
                                                      weights = kbal_demos_wedu_est$w[kbal_data_sampled ==1])
                    
                    kpop_demos_wedu <- est_mean("outcome", kpop_demos_wedu_svyd)
                    
                    numdims_demos_wedu = kbal_demos_wedu_est$numdims
                    mf_appended_dims_demos_wedu = kbal_demos_wedu_est$meanfirst_dims
                    if(is.null(numdims_demos_wedu)) {
                        numdims_demos_wedu = c(NA)
                        kpop_demos_wedu_se <- data.frame(SE_fixed = NA, 
                                                         SE_quasi = NA, 
                                                         SE_linear = NA, 
                                                         SE_chad = NA)
                    } else {
                        V <-  data.frame(kbal_dims = kbal_demos_wedu_est$svdK$v[, c(1:kbal_demos_wedu_est$numdims)])
                        X <- as.matrix(cbind(kbal_demos_wedu_est$appended_constraint_cols[kbal_data_sampled==1, ], V))
                        
                        cv_fit <- cv.glmnet(X, kpop_demos_wedu_svyd$variables$outcome, alpha = 0,
                                            penalty.factor = c(rep(0, ncol(kbal_demos_wedu_est$appended_constraint_cols)),
                                                               rep(1, kbal_demos_wedu_est$numdims)))
                        
                        lambda_pass = cv_fit$lambda.1se
                        residuals =  kpop_demos_wedu_svyd$variables$outcome - predict(cv_fit$glmnet.fit,
                                                                                      s = lambda_pass, 
                                                                                      newx = X)
                        res_kpop_demos_wedu = data.frame(min = min(residuals), 
                                                         perc_25 = quantile(residuals, .25), 
                                                         mean = mean(residuals),
                                                         perc_75 = quantile(residuals, .75),
                                                         var = var(residuals))
                        kpop_demos_wedu_se <- tryCatch(calc_SEs(Y = kpop_demos_wedu_svyd$variables$outcome,
                                                                residuals = residuals,
                                                                pop_size = nrow(cces),
                                                                sample_size = sum(sample),
                                                                weights = weights(kpop_demos_wedu_svyd)), 
                                                       error = function(e) NA)
                        if(length(kpop_demos_wedu_se) == 1) {
                            kpop_demos_wedu_se <- data.frame(SE_fixed = NA, 
                                                             SE_quasi = NA, 
                                                             SE_linear = NA, 
                                                             SE_chad = NA)
                        }
                        names(kpop_demos_wedu_se) = tryCatch(paste0("kpop_demos_wedu_", names(kpop_demos_wedu_se)),
                                                             error = function(e) NA)
                        
                    }
                    biasbound_r_demos_wedu = kbal_demos_wedu_est$biasbound_ratio
                    biasbound_demos_wedu = kbal_demos_wedu_est$biasbound_opt
                    
                    rm(kbal_demos_wedu_est, residuals, X, V, cv_fit)
                    
                    
                    #########all constraint method:
                    cat(paste("nsim:", nsim, "CONSTR (All)", "\n"))
                    kbal_all_est <- kbal(K.svd = svdK,
                                         allx=kbal_data,
                                         cat_data = TRUE,
                                         sampled = kbal_data_sampled,
                                         ebal.tol = tolerance,
                                         ebal.maxit = maxit,
                                         minnumdims = min_num_dims,
                                         maxnumdims = max_num_dims,
                                         scale_data = FALSE,
                                         drop_MC = FALSE,
                                         incrementby = increment,
                                         meanfirst = TRUE,
                                         sampledinpop = FALSE)
                    kpop_all_svyd <- svydesign(~1, data = survey_sim, 
                                               weights = kbal_all_est$w[kbal_data_sampled ==1])
                    
                    kpop_all <- est_mean("outcome", kpop_all_svyd)
                    
                    numdims_all = kbal_all_est$numdims
                    mf_appended_dims_all = kbal_all_est$meanfirst_dims
                    if(is.null(numdims_all)) {
                        numdims_all = c(NA)
                        numdims_all_se <- data.frame(SE_fixed = NA, 
                                                     SE_quasi = NA, 
                                                     SE_linear = NA, 
                                                     SE_chad = NA)
                    } else {
                        V <-  data.frame(kbal_dims = kbal_all_est$svdK$v[, c(1:kbal_all_est$numdims)])
                        X <- as.matrix(cbind(kbal_all_est$appended_constraint_cols[kbal_data_sampled==1, ], V))
                        
                        cv_fit <- cv.glmnet(X, kpop_all_svyd$variables$outcome, alpha = 0, 
                                            penalty.factor = c(rep(0, ncol(kbal_all_est$appended_constraint_cols)),
                                                               rep(1, kbal_all_est$numdims)))
                        
                        lambda_pass = cv_fit$lambda.1se
                        residuals =  kpop_all_svyd$variables$outcome - predict(cv_fit$glmnet.fit,
                                                                               s = lambda_pass, 
                                                                               newx = X)
                        res_kpop_all = data.frame(min = min(residuals), 
                                                  perc_25 = quantile(residuals, .25), 
                                                  mean = mean(residuals),
                                                  perc_75 = quantile(residuals, .75),
                                                  var = var(residuals))
                        kpop_all_se <- tryCatch(calc_SEs(Y = kpop_all_svyd$variables$outcome,
                                                         residuals = residuals,
                                                         pop_size = nrow(cces),
                                                         sample_size = sum(sample),
                                                         weights = weights(kpop_all_svyd)), 
                                                error = function(e) NA)
                        if(length(kpop_demos_wedu_se) == 1) {
                            kpop_demos_wedu_se <- data.frame(SE_fixed = NA, 
                                                             SE_quasi = NA, 
                                                             SE_linear = NA, 
                                                             SE_chad = NA)
                        }
                        names(kpop_all_se) = tryCatch(paste0("kpop_all_", names(kpop_all_se)),
                                                      error = function(e) NA)
                        
                    }
                    
                    biasbound_r_all = kbal_all_est$biasbound_ratio
                    biasbound_all = kbal_all_est$biasbound_opt
                    
                    rm(kbal_all_est, residuals, X,V, cv_fit)
                } else {
                    kpop_demos <- NA
                    numdims_demos = c(NA) 
                    kpop_demos_se <- data.frame(SE_fixed = NA, 
                                                SE_quasi = NA, 
                                                SE_linear = NA, 
                                                SE_chad = NA)
                    res_kpop_demos = data.frame(min = NA, 
                                                perc_25 = NA, 
                                                mean = NA,
                                                perc_75 = NA,
                                                var = NA)
                    biasbound_r_demos = NA
                    biasbound_demos = NA
                    
                    kpop_demos_wedu <- NA
                    numdims_demos_wedu = c(NA) 
                    kpop_demos_wedu_se <- data.frame(SE_fixed = NA, 
                                                     SE_quasi = NA, 
                                                     SE_linear = NA, 
                                                     SE_chad = NA)
                    res_kpop_demos_wedu =  data.frame(min = NA, 
                                                      perc_25 = NA, 
                                                      mean = NA,
                                                      perc_75 = NA,
                                                      var = NA)
                    biasbound_r_demos_wedu = NA
                    biasbound_demos_wedu = NA
                    
                    kpop_all <- NA
                    numdims_all = c(NA) 
                    kpop_all_se <- data.frame(SE_fixed = NA, 
                                              SE_quasi = NA, 
                                              SE_linear = NA, 
                                              SE_chad = NA)
                    res_kpop_all =  data.frame(min = NA, 
                                               perc_25 = NA, 
                                               mean = NA,
                                               perc_75 = NA,
                                               var = NA)
                    biasbound_r_all = NA
                    biasbound_all = NA
                    
                    #for weights
                    kpop_demos_svyd = survey_design
                    kpop_demos_wedu_svyd = survey_design
                    kpop_all_svyd = survey_design
                }
                
                
                rm(svdK)
                
                ##### return
                kpop_res = list()
                kpop_res$sims = data.frame(b_out = b_kpop,
                                           kpop = kpop,
                                           kpop_conv = kpop_conv,
                                           kpop_demos = kpop_demos,
                                           kpop_demos_wedu = kpop_demos_wedu,
                                           kpop_all = kpop_all,
                                           bb = biasbound,
                                           bbr = biasbound_r,
                                           bb_conv = biasbound_conv,
                                           bbr_conv = biasbound_r_conv,
                                           bb_demos = biasbound_demos,
                                           bbr_demos = biasbound_r_demos,
                                           bb_demos_wedu = biasbound_demos_wedu,
                                           bbr_demos_wedu = biasbound_r_demos_wedu,
                                           bb_all = biasbound_all,
                                           bbr_all = biasbound_r_all,
                                           numdims,
                                           numdims_conv,
                                           numdims_demos,
                                           numdims_demos_wedu,
                                           numdims_all,
                                           mf_appended_dims_demos, 
                                           mf_appended_dims_demos_wedu, 
                                           mf_appended_dims_all)
                
                #Standard Errors:
                kpop_res$SEs = data.frame(rake_demos_noeduc_se,
                                          rake_demos_noeduc_se_SVY,
                                          rake_demos_weduc_se,
                                          rake_all_se,
                                          post_stratification_se,
                                          rake_truth_se,
                                          kpop_se,
                                          kpop_conv_se,
                                          kpop_demos_se,
                                          kpop_demos_wedu_se,
                                          kpop_all_se)
                
                #weights
                kpop_res$weights = list(b = b_kpop,
                                        kpop_w = weights(kpop_svyd),
                                        kpop_w_conv = weights(kpop_svyd_conv),
                                        kpop_demos_w = weights(kpop_demos_svyd),
                                        kpop_demos_wedu_w = weights(kpop_demos_wedu_svyd),
                                        kpop_all_w = weights(kpop_all_svyd))
                
                #residuals
                kpop_res$residuals = rbind(b = b_kpop,
                                           kpop = res_kpop,
                                           kpop_conv = res_kpop_conv,
                                           kpop_demos = res_kpop_demos,
                                           kpop_demos_wedu = res_kpop_demos_wedu,
                                           kpop_all = res_kpop_all,
                                           rake_truth = res_rake_truth,
                                           rake_demos = res_rake_demos_noeduc,
                                           rake_demos_wedu = res_rake_demos_wedu,
                                           rake_all = res_rake_all)
                
                rm(kpop_svyd, 
                   kpop_svyd_conv,
                   kpop_demos_svyd,
                   kpop_demos_wedu_svyd, kpop_all_svyd)
                
            }
            
            ############################################ OUTPUT
            #out = list()
            if(eval_kpop) {
                sims = rbind(sims,
                                 data.frame(nsim,
                                 n,
                                 unweighted,
                                 rake_demos_noeduc,
                                 rake_demos_weduc,
                                 rake_all,
                                 post_stratification,
                                 rake_truth,
                                 ht_truth, 
                                 hayek_truth,
                                 kpop_res$sims))
                
                SEs = rbind(SEs, data.frame(kpop_res$SEs))
                weights = rbind(weights, data.frame(kpop_res$weights))
                residuals_out = rbind(residuals_out, kpop_res$residuals)
                dropped_cells_out = c(dropped_cells_out, dropped_cells_pass)
                
                sample_out = rbind(sample_out, 
                                   c(drop_ps = count, 
                               bad_sample = bad_sample, 
                               check = check_nums))
               #rm(kpop_res)
                
            } else {
                
                sims = rbind(sims, data.frame(nsim,
                                      n,
                                      unweighted,
                                      rake_demos_noeduc,
                                      rake_demos_weduc,
                                      rake_all,
                                      post_stratification,
                                      rake_truth,
                                      ht_truth, 
                                      hayek_truth))
                
                SEs = rbind(SEs, data.frame(rake_demos_noeduc_se,
                                     rake_demos_weduc_se,
                                     rake_demos_noeduc_se_SVY,
                                     rake_all_se,
                                     rake_truth_se,
                                     post_stratification_se))
                
                
                residuals_out = rbind(residuals_out,
                                  rbind(b = NULL,
                                      rake_truth = res_rake_truth,
                                      rake_demos = res_rake_demos_noeduc,
                                      rake_demos_wedu = res_rake_demos_wedu ,
                                      rake_all = res_rake_all)
                )
                
                dropped_cells_out = c(dropped_cells_out, dropped_cells_pass)
                
                sample_out = rbind(sample_out, c(drop_ps = count, 
                               bad_sample = bad_sample, 
                               check = check_nums))
            } 
            
            out$sims = sims
            out$SEs = SEs
            out$weights = weights
            out$residuals = residuals_out
            out$dropped_cells = dropped_cells_out
            out$sample = sample_out 
        })
        
    }
    
})


################################## Clean Sims ##################
outcome = cces$outcome
sims = out

if(SAVE) {
    save(sims, outcome, tolerance, maxit, increment, min_num_dims, max_num_dims, noise, R2_outcome, 
         eval_kpop,
         coefs, coefs_outcome, selection_model, p_include, pS_denom,
         file = paste0(save_path, "./sims_kpop_",
                       Sys.Date(),
                       "_nsims", nsims,
                       ".RData"))
}
