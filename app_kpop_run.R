### Packages
library(tidyverse)
library(survey)
library(parallel)
#latest version:
#devtools::install_github("csterbenz1/KBAL", ref = "cat_kernel") 
library(kbal)
library(glmnet)


#the following file runs kpop for each specification in the Application Section 6 regarding 2016 Election 
#given the size of the data, this is time consuming (svd of a very large matrix is computationally intensive)
#note that it does not create or analyze results but simply runs kpop, but simply saves the results to load in the results file app_res.Rmd

#set appropriate path for prepared cces and pew data files, created in data_prep.R
path_data= paste0(getwd(),"/generated_data/" ) 
#path where weights will be saved
path_save = paste0(getwd(), "/weights/")
#to improve runtime and prevent issues, saveMEM will toggle removing the large kpop result objects from the environment
#after the main results have been extracted
saveMEM = TRUE

#Parameters to control details of kbal() run
tolerance = 1e-6
maxit = 500
increment = 1
min_num_dims = 1
max_num_dims = 500
manual_lambda = FALSE 
#T=lambda as that which minimizes cverror in residualization; F= 1 sd from min choice
lambda_min = FALSE 


############################ Raking Specifications ##############################
formula_rake_demos_noeduc <- ~recode_age_bucket + recode_female + recode_race +
    recode_region + recode_pid_3way

#updated to include 6 way edu
formula_rake_demos_weduc <- ~recode_age_bucket + recode_female +
    recode_race + recode_region + recode_educ + recode_pid_3way

formula_rake_all_vars <- ~recode_age_bucket + recode_female +
    recode_race + recode_region + recode_pid_3way + recode_educ +
    recode_income_5way + recode_relig_6way + recode_born + recode_attndch_4way

############################ Various Standard Error Functions ######################
var_fixed <- function(Y, weights, pop_size) {
    ## note: needs weights that sum to population total
    if(round(sum(weights)) != pop_size) { weights = weights*pop_size/sum(weights)}
    return(Hmisc::wtd.var(Y, weights))
}

## kott (14) (under poisson)
var_quasi <- function(weights, residuals, pop_size) {
    #moving from kott 14 w/sum(w) =N to weights that sum to 1 + var of total to var of mean:
    #sum^n (w_i^2 - 1/N_pop w_i)e_i^2 
    return(sum((weights^2 - (weights / pop_size))*residuals^2))
}

## kott (15) linearization
var_linear <- function(weights, residuals, sample_size) {
    #moving from kott 14 w sum w =N to weights that sum to 1 + var of total to var of mean:
    # n/(n-1) sum^n (w_i*e_i)^2 - (1/n-1) [sum^n] *using this for now
    # approx = sum^n (w_i*e_i)^2 - (1/n) [sum^n]
    n = sample_size
    return((n/(n-1))*sum((weights * residuals)^2) - 1/(n-1) * sum(weights * residuals)^2)
}

## Chad's version
var_chad <- function(weights, residuals) {
    return(sum(weights^2 * residuals^2))
}

## Calculate all variances
calc_SEs <- function(Y, residuals, pop_size, weights, sample_size) {
    if(round(sum(weights)) != 1 ) {
        weights = weights/sum(weights)
    }
    return(data.frame(SE_fixed = sqrt(var_fixed(Y, weights, pop_size) / length(Y)),
                      SE_quasi = sqrt(var_quasi(weights, residuals, pop_size)),
                      SE_linear = sqrt(var_linear(weights, residuals, sample_size)),
                      SE_chad = sqrt(var_chad(weights, residuals))))
}

############################## Load Data ##########################
## SURVEY DATA (PEW)
### Load: pew_prepared_orig.rds is authors original file
#if prepared manually with data_prep_rep.R then should be named pew_prepared.rds
pew <- readRDS(paste0(path_data, "pew_prepared_orig.rds"))
## AUXILIARY INFORMATION (CCES)
### Load cces_prepared_orig.rds is authors original file
#if prepared manually with data_prep_rep.R then should be named cces_prepared.rds
cces <- readRDS(paste0(path_data, "cces_prepared_orig.rds"))
#check
#sum to N
#no Nas in outcome
sum(cces$commonweight_vv_post)
sum(is.na(cces$recode_vote_2016))

############################## Define Outcome and Setup Data ##########################
#outcome is projected cces modeled vote difference created in data_prep.R
#it is the modeled difference in probability of voting dem versus rep (p(D) - p(R)) in a 3way lasso multinomial
pew = pew %>% mutate(outcome = diff_cces_on_pew)
cces = cces %>% mutate(outcome = diff_cces_on_cces)

#all available variables
kbal_data <- bind_rows(pew %>% dplyr::select(recode_age_bucket,
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
kbal_data_sampled <- c(rep(1, nrow(pew)), rep(0, nrow(cces)))

############################## Run: KPOP ####################################################
#################### 1) kpop default ###########################
kbal_est <- kbal(allx=kbal_data,
                 sampled = kbal_data_sampled,
                 cat_data = TRUE,
                 incrementby = increment,
                 meanfirst = FALSE,
                 ebal.tol = tolerance,
                 ebal.maxit = maxit,
                 population.w = cces$commonweight_vv_post,
                 minnumdims = min_num_dims,
                 maxnumdims = max_num_dims,
                 sampledinpop = FALSE,
                 fullSVD = TRUE)
#kbal() object is very large, save just estimates and diagnostic measures
kpop_svyd <- svydesign(~1, data = pew,
                       weights = kbal_est$w[kbal_data_sampled ==1])
kpop <- svymean(~outcome, kpop_svyd,na.rm = TRUE)*100
b_kpop = kbal_est$b
l1_orig = ifelse(!is.null(kbal_est$L1_orig),kbal_est$L1_orig, NA)
l1 = ifelse(!is.null(kbal_est$L1_opt),kbal_est$L1_opt, NA)
numdims = kbal_est$numdims
biasbound_r = kbal_est$biasbound_ratio
biasbound = kbal_est$biasbound_opt

#save memory by saving only the svd to re use and rm kbal object
svdK = kbal_est$svdK 

##### Kpop SEs
x <- as.matrix(data.frame(kbal_dims = kbal_est$svdK$v[, 1:kbal_est$numdims]))
cv_fit <- cv.glmnet(x, kpop_svyd$variables$outcome, alpha = 0)
#use 1SE within min CV error lamba choice
lambda_pass =  cv_fit$lambda.1se
residuals = kpop_svyd$variables$outcome - predict(cv_fit$glmnet.fit,
                                                  s = lambda_pass, newx = x)
res_kpop = data.frame(min = min(residuals), 
                      perc_25 = quantile(residuals, .25), 
                      mean = mean(residuals),
                      perc_75 = quantile(residuals, .75),
                      var = var(residuals))
kpop_se <- calc_SEs(Y = kpop_svyd$variables$outcome,
                    residuals = residuals,
                    pop_size = nrow(cces),
                    sample_size = nrow(pew),
                    weights = weights(kpop_svyd))
rownames(kpop_se) = "kpop"

################### 2) kpop + forcing ebal convergence ########################
dist_record = data.frame(t(kbal_est$dist_record))
min_converged = dist_record[which.min(dist_record[dist_record$Ebal.Convergence ==1,"BiasBound"]), "Dims"]

if(saveMEM) { rm(kbal_est) }

if(is.null(min_converged) | length(min_converged) ==0) {
    kpop_conv_svyd <- "dn converge"
    kpop_conv <- "dn converge"
    
    numdims_conv = "dn converge"
    biasbound_r_conv = "dn converge"
    biasbound_conv = "dn converge"
    kpop_conv_se = data.frame(SE_fixed = NA, 
                              SE_quasi = NA, 
                              SE_linear = NA, 
                              SE_chad = NA)
    
} else {
    kbal_est_conv <- kbal(allx=kbal_data,
                          sampled = kbal_data_sampled,
                          cat_data = TRUE,
                          K.svd = svdK,
                          numdims = min_converged,
                          ebal.tol = tolerance,
                          ebal.maxit = maxit,
                          minnumdims = min_num_dims,
                          maxnumdims = max_num_dims,
                          population.w = cces$commonweight_vv_post,
                          incrementby = increment,
                          meanfirst = FALSE,
                          sampledinpop = FALSE,
                          ebal.convergence = TRUE)
    #save estiamtes and diagnostics
    kpop_conv_svyd <- svydesign(~1, data = pew,
                                weights = kbal_est_conv$w[kbal_data_sampled ==1])
    kpop_conv <- svymean(~outcome, kpop_conv_svyd,na.rm = TRUE)*100
    numdims_conv = kbal_est_conv$numdims
    biasbound_r_conv = kbal_est_conv$biasbound_ratio
    biasbound_conv = kbal_est_conv$biasbound_opt
    l1_conv = ifelse(!is.null(kbal_est_conv$L1_opt), kbal_est_conv$L1_opt, NA)
    
    ### Kpop + Convergeced SEs
    x <- as.matrix(data.frame(kbal_dims = kbal_est_conv$svdK$v[, 1:kbal_est_conv$numdims]))
    cv_fit <- cv.glmnet(x, kpop_conv_svyd$variables$outcome, alpha = 0)
    lambda_pass =  cv_fit$lambda.1se
    fit <- cv_fit$glmnet.fit
    residuals = kpop_conv_svyd$variables$outcome - predict(cv_fit$glmnet.fit, 
                                                           s = lambda_pass, 
                                                           newx = x)
    res_kpop_conv = data.frame(min = min(residuals), 
                               perc_25 = quantile(residuals, .25), 
                               mean = mean(residuals),
                               perc_75 = quantile(residuals, .75),
                               var = var(residuals))
    kpop_conv_se <- calc_SEs(Y = kpop_conv_svyd$variables$outcome,
                             residuals = residuals,
                             pop_size = nrow(cces),
                             sample_size = nrow(pew),
                             weights = weights(kpop_conv_svyd))
    rownames(kpop_conv_se) = "kpop_conv"
    
    if(saveMEM) { rm(kbal_est_conv)  }
    
}

############### 3) kpop + MF (all) ############################
kbal_demos_est <- kbal(K.svd = svdK,
                       allx=kbal_data,
                       cat_data = TRUE,
                       sampled = kbal_data_sampled,
                       ebal.tol = tolerance,
                       ebal.maxit = maxit,
                       minnumdims = min_num_dims,
                       maxnumdims = max_num_dims,
                       population.w = cces$commonweight_vv_post,
                       incrementby = increment,
                       meanfirst = TRUE,
                       mf_columns = all.vars(formula_rake_demos_noeduc),
                       sampledinpop = FALSE)
#save estimates and diagnostics
kpop_demos_svyd <- svydesign(~1, data = pew, 
                             weights = kbal_demos_est$w[kbal_data_sampled ==1])

kpop_demos <- svymean(~outcome, kpop_demos_svyd, na.rm = TRUE)*100
numdims_demos = kbal_demos_est$numdims
l1_demos = ifelse(!is.null(kbal_demos_est$L1_opt), kbal_demos_est$L1_opt, NA)
biasbound_r_demos = kbal_demos_est$biasbound_ratio
biasbound_demos = kbal_demos_est$biasbound_opt

#check to ensure we found convergence with this MF constraint added
if(is.null(numdims_demos)) {
    numdims_demos = c(NA) 
    kpop_demos_se <- data.frame(SE_fixed = NA, 
                                SE_quasi = NA, 
                                SE_linear = NA, 
                                SE_chad = NA)
} else {
    #residuals for SEs in MF only regularize on dimensions of svd(K)
    V <-  data.frame(kbal_dims = kbal_demos_est$svdK$v[, c(1:kbal_demos_est$numdims)])
    X <- as.matrix(cbind(kbal_demos_est$appended_constraint_cols[kbal_data_sampled==1, ], V))
    
    cv_fit <- cv.glmnet(X, kpop_demos_svyd$variables$outcome, alpha = 0,
                        penalty.factor = c(rep(0, ncol(kbal_demos_est$appended_constraint_cols)), 
                                           rep(1, kbal_demos_est$numdims)))
    
    lambda_pass = cv_fit$lambda.1se
    residuals =  kpop_demos_svyd$variables$outcome - predict(cv_fit$glmnet.fit,
                                                             s = lambda_pass, 
                                                             newx = X)
    res_kpop_demos = data.frame(min = min(residuals), 
                                perc_25 = quantile(residuals, .25), 
                                mean = mean(residuals),
                                perc_75 = quantile(residuals, .75),
                                var = var(residuals))
    
    kpop_demos_se <- calc_SEs(Y = kpop_demos_svyd$variables$outcome,
                                       residuals = residuals,
                                       pop_size = nrow(cces),
                                       sample_size = nrow(pew),
                                       weights = weights(kpop_demos_svyd))
    rownames(kpop_demos_se) = "kpop_demos"
}

if(saveMEM) { rm(kbal_demos_est) }

################# 4) kpop + MF (demos + educ) ##########################
kbal_demos_wedu_est <- kbal(K.svd = svdK,
                            allx=kbal_data,
                            cat_data = TRUE,
                            sampled = kbal_data_sampled,
                            ebal.tol = tolerance,
                            ebal.maxit = maxit,
                            minnumdims = min_num_dims,
                            maxnumdims = max_num_dims,
                            population.w = cces$commonweight_vv_post,
                            incrementby = increment,
                            meanfirst = TRUE,
                            mf_columns = all.vars(formula_rake_demos_weduc),
                            sampledinpop = FALSE)
#save estimates and diagnostics
kpop_demos_wedu_svyd <- svydesign(~1, data = pew, 
                                  weights = kbal_demos_wedu_est$w[kbal_data_sampled ==1])
kpop_demos_wedu <- svymean(~outcome, kpop_demos_wedu_svyd, na.rm = TRUE)*100
numdims_demos_wedu = kbal_demos_wedu_est$numdims
l1_demos_weduc = ifelse(!is.null(kbal_demos_wedu_est$L1_opt), kbal_demos_wedu_est$L1_opt, NA)
biasbound_r_demos_wedu = kbal_demos_wedu_est$biasbound_ratio
biasbound_demos_wedu = kbal_demos_wedu_est$biasbound_opt

#check we found convergence with this MF constraint added
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
    kpop_demos_wedu_se <- calc_SEs(Y = kpop_demos_wedu_svyd$variables$outcome,
                                            residuals = residuals,
                                            pop_size = nrow(cces),
                                            sample_size = nrow(pew),
                                            weights = weights(kpop_demos_wedu_svyd))
    rownames(kpop_demos_wedu_se) = "kpop_demos_wedu"
}

if(saveMEM) { rm(kbal_demos_wedu_est)  }

################ 5) kpop + MF (demos + educ) ###########################
#########KPOP + all constraint method:
kbal_all_est <- kbal(K.svd = svdK,
                     allx=kbal_data,
                     cat_data = TRUE,
                     sampled = kbal_data_sampled,
                     ebal.tol = tolerance,
                     ebal.maxit = maxit,
                     minnumdims = min_num_dims,
                     maxnumdims = max_num_dims,
                     population.w = cces$commonweight_vv_post,
                     incrementby = increment,
                     meanfirst = TRUE,
                     mf_columns = all.vars(formula_rake_all_vars),
                     sampledinpop = FALSE)

kpop_all_svyd <- svydesign(~1, data = pew, 
                           weights = kbal_all_est$w[kbal_data_sampled ==1])
kpop_all <- svymean(~outcome, kpop_all_svyd, na.rm = TRUE)*100
numdims_all = kbal_all_est$numdims
l1_all = ifelse(!is.null(kbal_all_est$L1_opt),  kbal_all_est$L1_opt, NA)
biasbound_r_all = kbal_all_est$biasbound_ratio
biasbound_all = kbal_all_est$biasbound_opt

#check that we found convergecne for this specification of MF contstraints
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
    kpop_all_se <- alc_SEs(Y = kpop_all_svyd$variables$outcome,
                                     residuals = residuals,
                                     pop_size = nrow(cces),
                                     sample_size = nrow(pew),
                                     weights = weights(kpop_all_svyd))
    rownames(kpop_all_se) = "kpop_all"
}

if(saveMEM) {
    rm(kbal_all_est)
    rm(svdK)
}

########### Compile results of all Kpop runs and save#############
out = list()
out$est = data.frame(b_out,
                     tolerance = tolerance, 
                     maxit = maxit,
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
                     l1_orig ,
                     l1,
                     l1_conv,
                     l1_demos,
                     l1_demos_weduc,
                     l1_all)

#Standard Errors:
out$SEs = rbind(kpop_se,
                kpop_conv_se,
                kpop_demos_se,
                kpop_demos_wedu_se,
                kpop_all_se)

#weights
out$weights = list(b = b_out,
                   kpop_w = weights(kpop_svyd),
                   kpop_w_conv = weights(kpop_conv_svyd),
                   kpop_demos_w = weights(kpop_demos_svyd),
                   kpop_demos_wedu_w = weights(kpop_demos_wedu_svyd),
                   kpop_all_w = weights(kpop_all_svyd))

#residuals
out$residuals = rbind(b = b_out,
                      kpop_res = res_kpop,
                      kpop_w_conv = res_kpop_conv,
                      kpop_demos_w = res_kpop_demos,
                      kpop_demos_wedu_w = res_kpop_demos_wedu,
                      kpop_all_w = res_kpop_all)

###################################### 
#save output
save(out, tolerance, maxit, POPW,min_num_dims, max_num_dims, increment,TEST,
     file = paste0(path_save, "app_update_",
                   Sys.Date(), ".Rdata"))