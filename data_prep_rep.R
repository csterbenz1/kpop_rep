rm(list = ls())
library(MASS)
library(tidyverse)
library(survey)
library(srvyr)
library(glmnet)
library(caret)
library(foreign)
library(dataverse)

## Data Prep File proceeds in 4 main sections
## 1. Raw election results cleaning
## 2. Raw Pew pre-election wave cleaning
## 3. Raw CCES post-election wave cleaning
## 4. 3way regularized multinomial modeled outcome creation for use in 2016 Application (section 6)

#Output is all dumbped in:
sink("data_prep_out.Rout")
cat("Loaded packages are: \n")
(.packages())

#can adjust working_path as appropriate here
working_path = getwd()
cat(paste0("Cleaned and Prepared data will be saved in:\n\t", working_path, "/generated_data"))

### region definitions
northeast_states = c("Massachusetts", "Pennsylvania", "Connecticut",
                     "New Hampshire", "New York", "New Jersey", "Maine",
                     "Rhode Island", "Vermont")
midwest_states = c("Minnesota", "Iowa", "Nebraska", "Ohio", "Indiana",
                   "Wisconsin", "Kansas", "Missouri", "Illinois", "Michigan",
                   "South Dakota", "North Dakota") 
south_states = c("Delaware", "Texas", "North Carolina", "Georgia", "Alabama",
                 "West Virginia", "Maryland", "District of Columbia",
                 "Kentucky", "South Carolina", "Tennessee", "Louisiana",
                 "Oklahoma", "Florida", "Arkansas", "Virginia", "Mississippi") 
west_states = c("California", "Oregon", "Washington", "Hawaii", "New Mexico",
                "Alaska", "Colorado", "Utah", "Wyoming", "Arizona", "Nevada",
                "Montana", "Idaho")


################################################################################
#### 1. ELECTION RESULTS #######################################################
################################################################################

cat("Loading raw Election data\n")
if(!file.exists(paste0(working_path, "/data/1976-2020-president.dta"))) {
    if (!file.exists(paste0(working_path, "/data"))) {
        dir.create(file.path(paste0(working_path, "/data")))
    }
    cat("Downloading Election Data from dataverse and writing to local bin \n")
    writeBin(get_file_by_name(
        filename = "1976-2020-president.tab",
        dataset  = "doi:10.7910/DVN/42MVDX",
        server   = "dataverse.harvard.edu"
    ), paste0(working_path, "/data/1976-2020-president.dta"))
}
election_2016 <- tryCatch(read.dta(paste0(working_path, "/data/1976-2020-president.dta")), 
                          error = function(e) NA)
if(is.na(election_2016)) {
    warning("\nError loading from dataverse (may require dataverse login), try navigating to website using DOI and manually downloading to \"data\" folder as you will also later need to with the Pew data below and rerun", immediate. = TRUE)
    election <- read.csv( "./data/1976-2020-president.csv")
    
}
election_2016 <- election %>%
    filter(year == 2016) %>%
    select(state, party_simplified, candidatevotes, totalvotes) %>% 
    group_by(state) %>%
    summarize(demtotal = sum(candidatevotes[party_simplified == "DEMOCRAT"], na.rm = TRUE),
              reptotal = sum(candidatevotes[party_simplified == "REPUBLICAN"],
                             na.rm = TRUE),
              totalvotes = unique(totalvotes)) %>%
    ungroup()

cat(paste0("Saving cleaned 2016 Election results as: \"election.rds\" in \n\t", working_path, "/generated_data"))
saveRDS(election_2016, paste0(working_path, "/generated_data/election.rds"))



################################################################################
#### 2. Pew CLEANING ###########################################################
################################################################################

#Pew survey data is public, but does require creating an account; it can be accessed here
#https://www.pewresearch.org/dataset/october-2016-political-survey

### load original Pew .sav file downloaded from above link
#late october pew
cat("Loading raw Pew data\n")
pew <- read.spss("./data/Oct16 public.sav", to.data.frame=TRUE)

## Start by looking at missingness (fraction):
cat("Missingness in Pew data (fraction):\n")
lapply(pew[, c("age", "sex", "racethn", "state", "party", "educ2")], 
       function(x) (sum(is.na(x)) + sum(grepl("Don't know/Refused", x))) / length(x))

cat("Recoding/Cleaning basic demographics and outcome in Pew data\n")
## Outcome Variable:
#Q10 asked if you voted today who would you vote for among trump, clinton, stein, and johnson, 
#10a asked if they lean for these candidates if responded dnk
#Q11 forced a horse race between dem and rep (if you were forced to choose among)
#Q11 is only asked if they answer a third candidate for q10 or dnk
#Q11horse2 is q11 + q10 + q10a: trump: 816+42+70 = 928 and clinton: 944+28+80 = 1052 
pew <- pew %>%
    mutate(recode_vote_2016 =
               case_when(str_detect(Q11HORSE2, "Trump") ~ "Republican",
                         str_detect(Q11HORSE2, "Clinton") ~ "Democrat",
                         is.na(Q11HORSE2) ~ NA_character_,
                         TRUE ~ "Other"
               ))

## First recodes
pew <- pew %>% mutate(
  # age
  recode_age = ifelse(age == "Don't know/Refused (VOL.)", 
                      NA, 
                      #loads as a factor, undo that
                      as.numeric(as.character(age))),

  # gender
  recode_female = case_when(sex == "Female" ~ "Female",
                            TRUE ~ "Male"),
  # race/ethnicity
  # Note: combines missing with Other
  recode_race = case_when(racethn == "White, non-Hisp" ~ "White",
                          racethn == "Black, non-Hisp" ~ "Black",
                          racethn == "Hispanic" ~ "Hispanic",
                          TRUE ~ "Other"),
  # region (note there are no cases that do not fall into these groups such that we fill TRUE ~ South)
  recode_region = case_when(state %in% northeast_states ~ "Northeast",
                            state %in% west_states ~ "West",
                            state %in% midwest_states ~ "Midwest",
                            state %in% south_states ~ "South",
                            TRUE ~ "South"),
  
  # party -- combines refused + no preference with Independents
  recode_pid_3way = case_when( party == "Democrat" ~ "Dem",
                               party == "Republican" ~ "Rep",
                               TRUE ~ "Ind"),
  
  # state
  recode_inputstate = state,
  
  # education -- leaving missing for now (14 Don't now)
  recode_educ = factor(case_when( 
    educ2 == "Less than high school (Grades 1-8 or no formal schooling) " ~ "No HS",
    educ2 == "High school incomplete (Grades 9-11 or Grade 12 with NO diploma)" ~ "No HS",
    educ2 == "High school graduate (Grade 12 with diploma or GED certificate)" ~ "High school graduate",
    educ2 == "Some college, no degree (includes some community college)" ~ "Some college",
    educ2 == "Two year associate degree from a college or university" ~ "2-year",
    educ2 == "Four year college or university degree/Bachelor's degree (e.g., BS, BA, AB)" ~ "4-year",
    educ2 == "Some postgraduate or professional schooling, no postgraduate degree" ~ "Post-grad",
    educ2 == "Postgraduate or professional degree, including master's, doctorate, medical or law degree" ~ "Post-grad",
    TRUE ~ NA_character_), 
    levels = c("No HS", "High school graduate", "Some college", "2-year", "4-year", "Post-grad"))
)



#### adding indicator for any missingness
#note that region and inputstate and pid and female should not ever have any missing
pew <- pew %>% 
  mutate(missing = ifelse(is.na(recode_age) | is.na(recode_educ), 1, 0))

## modeling age
cat("Modeling age/income to delete missingness in Pew data\n")
## listwise delete missing on age/education
#45 observations dropped= 41 missing on age + 4 missing on educ
age_model <- glm(recode_age ~ recode_female + recode_race + recode_region + 
                   recode_pid_3way + recode_educ + child,
                 family = Gamma(),
                 data = pew)


#reports 41 observations dropped = those missing on age
age_model_no_educ <- glm(recode_age ~ recode_female + recode_race + 
                           recode_region + recode_pid_3way + child,
                 family = Gamma(),
                 data = pew)

## recode age and age buckets
pew <- pew %>%
  mutate(recode_age = case_when(
            ## if age is not missing, keep current age
            !is.na(recode_age) ~ recode_age,
            ## if age is missing but not education, use age_model
            is.na(recode_age) & !is.na(recode_educ) ~ round(predict(age_model, ., type = "response"), 0),
            ## if age and education are missing, use age_model_no_educ
            is.na(recode_educ) ~ round(predict(age_model_no_educ, ., type = "response"), 0)
  ),
  
  ## four way age bucket; should be no na's at this point
  recode_age_bucket = factor(case_when( recode_age <= 35 ~ "18 to 35",
                                        recode_age <= 50 ~ "36 to 50",
                                        recode_age < 65 ~ "51 to 64",
                                        !is.na(recode_age) ~ "65+"),
                             levels = c("18 to 35", "36 to 50", "51 to 64", "65+")),
  
  ## three way age bucket; should be no NAs at this point
  recode_age_3way = case_when( recode_age <= 50 ~ "a_18to50",
                               recode_age < 65 ~ "b_51to64",
                               !is.na(recode_age) ~ "c_65")
  )

## education model (ordered logit)
#drops 14 missing on educ
educ_model <- polr(recode_educ ~ recode_female + recode_race + recode_region + 
                     recode_pid_3way + recode_age + child,
                   data = pew)

pew <- pew %>%
  ## impute missing educ bucket
  mutate(recode_educ = case_when(
            ## if not missing, use current education
            !is.na(recode_educ) ~ recode_educ,
            ## if missing, use educ_model prediction
            is.na(recode_educ) ~ predict(educ_model, newdata =., type = "class")
  ),
  recode_educ_3way = factor(case_when(
      recode_educ %in% c("No HS", "High school graduate", "Some college") ~ "No College",
      recode_educ %in% c("2-year", "4-year") ~ "College",
      TRUE ~ "Post-grad"), 
  levels = c("No College", "College", "Post-grad"))
)  
  
####### Recode Variables: income, born again, church attendance,
### Pew
cat("Recoding/Cleaning additional variables in Pew data\n")
pew <- pew %>% mutate(
    recode_relig = factor(case_when(
        #say Unitarians are protestant (only 4 of these), Christians are protestant (269)
        relig == "Protestant (Baptist, Methodist, Non-denominational, Lutheran...)" ~ "Protestant",
        relig == "Christian (VOL.)" ~ "Protestant",
        relig == "Unitarian (Universalist) (VOL.)" ~ "Protestant",
        relig == "Jewish (Judaism)" ~ "Jewish",
        relig == "Roman Catholic (Catholic)" ~ "Catholic",
        relig == "Nothing in particular" ~ "Nothing in particular",
        relig == "Agnostic (not sure if there is a God)" ~ "Agnostic",
        relig == "Atheist (do not believe in God)" ~ "Atheist",
        relig == "Buddhist" ~ "Buddhist",
        relig == "Muslim (Islam)" ~ "Muslim",
        relig == "Mormon (Church of Jesus Christ of Latter-day Saints/LDS)" ~ "Mormon",
        relig == "Hindu" ~ "Hindu",
        relig == "Orthodox (Greek, Russian, or some other orthodox church)" ~ "Orthodox",
        #these people were asked a follow up if considered christian in $chr
        #out of 41 (31) of these: 11 (8) said yes Christ, 29 (22) no, 1 refused
        relig == "Something else (SPECIFY:______)" ~ "Something else",
        #those who respond don't know/refuse =NA (36) (NB: also asked if consider christian in chr)
        #out of these 50 (36): 30 (24) said yes christian, 11 (5) said no, 9 (7) refused
        relig == "Don't Know/Refused (VOL.)" & chr == "Yes" ~ "Protestant",
        #leaves only 20 (12) missing
        TRUE ~ NA_character_),
        levels = c("Protestant", "Jewish", "Catholic", "Nothing in particular", 
                   "Agnostic", "Atheist", "Buddhist", "Muslim", "Mormon", "Hindu",
                   "Orthodox", "Something else")),
    
    #supposedly only asked: if chr = yes or if relig = 1-4,13
    #(Protestant, Catholic, Mormon, Orthodox, Christian)
    #though looking at the NAs that's not totally true some in these groups still have NAs
    recode_born = case_when(born == "Yes, would" ~ "Yes",
                            born == "No, would not" ~ "No",
                            #46 refused but we think it's safe to say these people are not
                            born == "Don't know/Refused (VOL.)" ~ "No",
                            #NA's (589) these are mainly non-christian denominations
                            #so I'd say safe to say they are not born again
                            TRUE ~ "No"),
    
    #the names are actually ok here
    #leaving the 35 (28) that response don't know as they are 
    recode_attndch = attend, 
    
    recode_income = factor(case_when(
        income == "Less than $10,000" ~ "<10k",
        income == "10 to under $20,000" ~ "10-20k",
        income == "20 to under $30,000" ~ "20-30k",
        income == "30 to under $40,000" ~ "30-40k",
        income == "40 to under $50,000" ~ "40-50k",
        income == "50 to under $75,000" ~ "50-100k",
        income == "75 to under $100,000" ~ "50-100k",
        income == "100 to under $150,000 [OR]" ~ "100-150k",
        income == "$150,000 or more" ~ ">150k",
        #222 refused
        income == "[VOL. DO NOT READ] Don't know/Refused" ~ "Prefer not to say",
        #there should be no NAs, but as precaution lump them with prefer not to say
        TRUE ~ NA_character_),
        levels = c("<10k", "10-20k", "20-30k", "30-40k", "40-50k", "50-100k",
                   "100-150k",">150k", "Prefer not to say")),
    recode_income_5way = factor(case_when(
        income == "Less than $10,000" ~ "<20k",
        income == "10 to under $20,000" ~ "<20k",
        income == "20 to under $30,000" ~ "20-50k",
        income == "30 to under $40,000" ~ "20-50k",
        income == "40 to under $50,000" ~ "20-50k",
        income == "50 to under $75,000" ~ "50-100k",
        income == "75 to under $100,000" ~ "50-100k",
        income == "100 to under $150,000 [OR]" ~ "100-150k",
        income == "$150,000 or more" ~ ">150k",
        #222 refused
        income == "[VOL. DO NOT READ] Don't know/Refused" ~ "Prefer not to say",
        #there should be no NAs
        TRUE ~ NA_character_),
        levels = c("<20k", "20-50k", "50-100k", "100-150k",">150k", "Prefer not to say")),
    
    recode_income_4way = factor(case_when(
        income == "Less than $10,000" ~ "<50k",
        income == "10 to under $20,000" ~ "<50k",
        income == "20 to under $30,000" ~ "<50k",
        income == "30 to under $40,000" ~ "<50k",
        income == "40 to under $50,000" ~ "<50k",
        income == "50 to under $75,000" ~ "50-100k",
        income == "75 to under $100,000" ~ "50-100k",
        income == "100 to under $150,000 [OR]" ~ "100-150k",
        income == "$150,000 or more" ~ ">150k",
        #222 refused
        income == "[VOL. DO NOT READ] Don't know/Refused" ~ "Prefer not to say",
        #there should be no NAs
        TRUE ~ NA_character_),
        levels = c("<50k", "50-100k", "100-150k",">150k", "Prefer not to say")),
    
    recode_income_3way = factor(case_when(
        income == "Less than $10,000" ~ "<50k",
        income == "10 to under $20,000" ~ "<50k",
        income == "20 to under $30,000" ~ "<50k",
        income == "30 to under $40,000" ~ "<50k",
        income == "40 to under $50,000" ~ "<50k",
        income == "50 to under $75,000" ~ "50-150k",
        income == "75 to under $100,000" ~ "50-150k",
        income == "100 to under $150,000 [OR]" ~ "50-150k",
        income == "$150,000 or more" ~ ">150k",
        #222 refused
        income == "[VOL. DO NOT READ] Don't know/Refused" ~ "Prefer not to say",
        #there should be no NAs
        TRUE ~ NA_character_),
        levels = c("<50k", "50-150k",">150k","Prefer not to say"))
)


pew <- pew %>% mutate(
    #NB: this is actually 5way not 6, a naming mistake that is fixed in subsequent labels
    recode_relig_6way = factor(case_when(
        #clean up by saying those that are christian and not catholic is protestant
        recode_relig == "Protestant" ~ "Protestant",
        recode_relig == "Jewish" ~ "Jewish",
        recode_relig == "Catholic" ~ "Catholic",
        recode_relig == "Nothing in particular" ~ "Not religious",
        recode_relig == "Agnostic" ~ "Not religious",
        recode_relig == "Atheist" ~ "Not religious",
        recode_relig == "Buddhist" ~ "Other",
        recode_relig == "Muslim" ~ "Other",
        recode_relig == "Hindu" ~ "Other",
        recode_relig == "Mormon" ~ "Other",
        recode_relig == "Orthodox" ~ "Other",
        recode_relig == "Something else" & chr == "Yes"  ~ "Protestant",
        recode_relig == "Something else" & chr == "No"  ~ "Other",
        recode_relig == "Something else" & chr == "Don't know/Refused (VOL.)"  ~ "Not religious",
        #20 remaining NAs as Not religious
        TRUE ~ "Not religious"),
        levels = c("Protestant", "Jewish", "Catholic","Not religious", "Other")),

    recode_attndch_4way = factor(case_when(
        recode_attndch == "Once a week" ~ "Weekly",
        recode_attndch == "More than once a week" ~ "Weekly",
        recode_attndch == "Once or twice a month" ~ "Monthly",
        recode_attndch == "A few times a year" ~ "Yearly",
        recode_attndch == "Seldom" ~ "Never",
        recode_attndch == "Never" ~ "Never",
        recode_attndch == " Don't know/Refused (VOL.)" ~ "Never",
        #the 35 who refused say never go
        TRUE ~ "Never"),
        levels = c("Weekly","Monthly", "Yearly", "Never"))
)

####Interactions
cat("Creating a few interactions in PEW\n")
pew = pew %>% mutate(
    #12 levels
    recode_pid_race = as.factor(paste(recode_pid_3way, 
                                      recode_race, sep = ", ")),
    #9 levels
    recode_educ_pid = as.factor(paste(recode_educ_3way,
                                      recode_pid_3way, sep = ", ")),
    #36 levels
    recode_educ_pid_race = as.factor(paste(recode_educ_3way,
                                           recode_pid_3way, 
                                           recode_race, sep = ", ")),
    ## education 3 way bucket among white voters, no split among non-white
    recode_educ_wh_3way = factor(case_when(
        recode_race != "White" ~ "No Split",
        recode_educ %in% c("No HS", "High school graduate", "Some college") ~ "No College",
        recode_educ %in% c("2-year", "4-year") ~ "College",
        TRUE ~ "Post-grad"), 
        levels = c("No Split", "No College", "College", "Post-grad")),
    
    recode_race_reg_wh_educ = as.factor(paste(recode_region, 
                                           case_when(recode_race == "White" ~ paste(recode_race, recode_educ, sep = ", "),
                                                     TRUE ~ recode_race), sep = ", ")), 
    recode_race_educ_reg = as.factor(paste(recode_race, recode_region, 
                                           recode_educ_3way, sep = ", ")),
    #this has 4 levels
    recode_midwest_wh_edu = factor(case_when(
        (recode_race != "White" | recode_region != "Midwest") ~ "No Split", 
        TRUE ~ as.character(recode_educ_3way)), 
        levels = c("No Split", "No College", "College", "Post-grad")),
    #this has 13 levels
    recode_midwest_edu_race = as.factor(
        case_when(recode_region == "Midwest" ~ paste(recode_region,
                                                     recode_race,
                                                     recode_educ_3way, 
                                                     sep = ", "),
                  TRUE ~ "No Split")),
    #consider each individual age a "bucket"
    recode_age_factor = factor(case_when(recode_age <92 ~ as.character(recode_age), 
                                         TRUE ~ "92+")) 
)


## Remove those who say they definitely will not vote 
#463 NAs (all also have NA vote choice) and 46 don't plan to vote (some have vote pref but if they refuse then fine to drop)
#(5 of those that are NA here are missing also on age/educ so we get from 45 missing to 40 missing total)
cat("Dropping Pew respondents who said they definitely do not plan to vote\n")
pew <- pew %>%
  filter(plan1 %in% c("Plan to vote", "Already voted"))


################################################################################
#### 3. CCES CLEANING ##########################################################
################################################################################


## Get the 2016 CCES data from dataverse
## located at `https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi%3A10.7910/DVN/GDF6Z0`.
##
## Ansolabehere, Stephen; Schaffner, Brian F., 2017, "CCES Common Content, 2016", 
## https://doi.org/10.7910/DVN/GDF6Z0, Harvard Dataverse, V4, 
## UNF:6:WhtR8dNtMzReHC295hA4cg== [fileUNF]
## Download data `CCES16_Common_OUTPUT_Feb2018_VV.tab` as .dta

cat("Loading raw CCES data\n")
if(!file.exists(paste0(working_path, "/data/CCES16_Common_Content.dta"))) {
    if (!file.exists(paste0(working_path, "/data"))) {
        dir.create(file.path(paste0(working_path, "/data")))
    }
    cat("Downloading CCES Data from dataverse and writing to local bin (large, may take some time) \n")
    writeBin(get_file_by_name(
        filename = "CCES16_Common_OUTPUT_Feb2018_VV.tab",
        dataset  = "doi:10.7910/DVN/GDF6Z0",
        server   = "dataverse.harvard.edu"
    ), paste0(working_path, "/data/CCES16_Common_Content.dta"))
}
cces <- read.dta(paste0(working_path, "/data/CCES16_Common_Content.dta"))

# Common Content Weight -- commonweight_vv_post (post wave with vote validation)
### Drop invalid cases
cces <- cces %>%
    filter((CC16_401 == "I definitely voted in the General Election.") &
               !is.na(commonweight_vv_post)) 

#what's missing in CCES? just a few party id's
cat("Missingness in CCES data (raw count):\n")
lapply(cces[, c("birthyr", "inputstate_post", "gender", "race", "pid3", "educ")], 
       function(x) (sum(is.na(x)) + sum(grepl("Skipped", x)) + sum(grepl("Not Asked", x))))


cat("Recoding/Cleaning basic demographics and outcome in CCES data\n")
## Make outcome:
# Horserace (post) question 410a
cces <- cces %>%
    mutate( recode_vote_2016 = case_when(
        str_detect(CC16_410a, "Democrat") ~ "Democrat",
        str_detect(CC16_410a, "Republican") ~ "Republican",
        (as.numeric(CC16_410a) < 6 | as.numeric(CC16_410a) == 8) ~ "Other",
        #NAs = did not vote + not sure; (these are later dropped)
        #360: 229 = I'm not sure, 50 = NA,
        #81 = I didn't vote in this election (presumably only voted local since said voted in CC16_401)
        TRUE ~ NA_character_) )


## Demographic recodes
cces$inputstate_post = droplevels(cces$inputstate_post)

cces <- cces %>% mutate(
  recode_id = 1:nrow(cces),
  
  # age
  recode_age = 2016 - as.numeric(as.character(birthyr)),
  #4way
  recode_age_bucket = factor(case_when( recode_age <= 35 ~ "18 to 35",
                          recode_age <= 50 ~ "36 to 50",
                          recode_age < 65 ~ "51 to 64",
                          !is.na(recode_age) ~ "65+"), 
                          levels = c("18 to 35", "36 to 50", "51 to 64", "65+")),
  
  recode_age_3way = case_when( recode_age <= 50 ~ "a_18to50",
                          recode_age < 65 ~ "b_51to64",
                          !is.na(recode_age) ~ "c_65"),
  
  # gender
  recode_female = case_when(gender == "Female" ~ "Female",
                            TRUE ~ "Male"),
  
  # race: this means Mixed, Asian, Other, Native American, Middle Eastern all in Other
  recode_race = ifelse(race %in% c("White", "Black", "Hispanic"), as.character(race), "Other"),
  
  # region: there should be no states that do not fall in one of these categories (no NAs)
  recode_region = case_when(inputstate_post %in% northeast_states ~ "Northeast",
                            inputstate_post %in% west_states ~ "West",
                            inputstate_post %in% midwest_states ~ "Midwest",
                            inputstate_post %in% south_states ~ "South",
                            TRUE ~ "South"),
  
  # party: note 13 missing on pid are categorized as Indep
  recode_pid_3way = case_when( str_detect(pid3, "Democrat") ~ "Dem",
                        str_detect(pid3, "Republican") ~ "Rep",
                        TRUE ~ "Ind"),

  # educ
  recode_educ = factor(educ, levels = c("No HS", "High school graduate", 
                                        "Some college", "2-year", "4-year", 
                                        "Post-grad")),

  recode_educ_3way = factor(case_when(
      recode_educ %in% c("No HS", "High school graduate", "Some college") ~ "No College",
      recode_educ %in% c("2-year", "4-year") ~ "College",
      TRUE ~ "Post-grad"), 
      levels = c("No College", "College", "Post-grad")),
  
  # state
  recode_inputstate = inputstate_post
)

cces$recode_inputstate <- droplevels(cces$recode_inputstate)
cces$recode_educ <- droplevels(cces$recode_educ)

#### Recoding Variables:
cat("Recoding/Cleaning additional variables in CCES data\n")
cces <- cces %>% mutate(
    #if responded don't know/refused 
    #unitarians are protestant (only 4 of these), Christains are protestant (269)
    recode_relig = factor(case_when(
        religpew == "Protestant" ~ "Protestant",
        religpew == "Jewish" ~ "Jewish",
        religpew == "Roman Catholic" ~ "Catholic",
        religpew == "Nothing in particular" ~ "Nothing in particular",
        religpew == "Agnostic" ~ "Agnostic",
        religpew == "Atheist" ~ "Atheist",
        religpew == "Buddhist" ~ "Buddhist",
        religpew == "Muslim" ~ "Muslim",
        religpew == "Mormon" ~ "Mormon",
        religpew == "Hindu" ~ "Hindu",
        religpew == "Eastern or Greek Orthodox" ~ "Orthodox",
        religpew == "Something else" ~ "Something else",
        #NAs remain NAs 38
        TRUE ~ NA_character_),
        levels = c("Protestant", "Jewish", "Catholic", "Nothing in particular", 
                   "Agnostic", "Atheist", "Buddhist", "Muslim", "Mormon", "Hindu",
                   "Orthodox", "Something else")),
    
    #22 NAs
    recode_born = case_when(pew_bornagain == "No" ~ "No",
                            pew_bornagain == "Yes" ~ "Yes", 
                            #NAs (22) assume to be not
                            TRUE ~ "No"),
    
    #names ok; NAs (23) and don't know (319) 
    recode_attndch = pew_churatd,
    
    recode_income = factor(case_when(
        faminc == "Less than $10,000" ~ "<10k",
        faminc == "$10,000 - $19,999" ~ "10-20k",
        faminc == "$20,000 - $29,999" ~ "20-30k",
        faminc == "$30,000 - $39,999" ~ "30-40k",
        faminc == "$40,000 - $49,999" ~ "40-50k",
        faminc == "$50,000 - $59,999" ~ "50-100k",
        faminc == "$60,000 - $69,999" ~ "50-100k",
        faminc == "$70,000 - $79,999" ~ "50-100k",
        faminc == "$80,000 - $99,999" ~ "50-100k",
        faminc == "$100,000 - $119,999" ~ "100-150k",
        faminc == "$120,000 - $149,999" ~ "100-150k",
        faminc == "$150,000 - $199,999" ~ ">150k",
        faminc == "$200,000 - $249,999" ~ ">150k",
        faminc == "$250,000 - $349,999" ~ ">150k",
        faminc == "$350,000 - $499,999" ~ ">150k",
        faminc == "$500,000 or more" ~ ">150k",
        faminc == "$150,000 or more" ~ ">150k",
        #4788 do not say
        faminc == "Prefer not to say" ~ "Prefer not to say",
        #9 NAs say they are also prefer not to say
        TRUE ~ "Prefer not to say"),
        levels = c("<10k", "10-20k", "20-30k", "30-40k", "40-50k", "50-100k",
                   "100-150k",">150k","Prefer not to say")),
    
    recode_income_5way = factor(case_when(
        faminc == "Less than $10,000" ~ "<20k",
        faminc == "$10,000 - $19,999" ~ "<20k",
        faminc == "$20,000 - $29,999" ~ "20-50k",
        faminc == "$30,000 - $39,999" ~ "20-50k",
        faminc == "$40,000 - $49,999" ~ "20-50k",
        faminc == "$50,000 - $59,999" ~ "50-100k",
        faminc == "$60,000 - $69,999" ~ "50-100k",
        faminc == "$70,000 - $79,999" ~ "50-100k",
        faminc == "$80,000 - $99,999" ~ "50-100k",
        faminc == "$100,000 - $119,999" ~ "100-150k",
        faminc == "$120,000 - $149,999" ~ "100-150k",
        faminc == "$150,000 - $199,999" ~ ">150k",
        faminc == "$200,000 - $249,999" ~ ">150k",
        faminc == "$250,000 - $349,999" ~ ">150k",
        faminc == "$350,000 - $499,999" ~ ">150k",
        faminc == "$500,000 or more" ~ ">150k",
        faminc == "$150,000 or more" ~ ">150k",
        #4788 do not say
        faminc == "Prefer not to say" ~ "Prefer not to say",
        #9 NAs
        TRUE ~ "Prefer not to say"),
        levels = c("<20k", "20-50k", "50-100k", "100-150k",">150k","Prefer not to say")),
    recode_income_4way = factor(case_when(
        faminc == "Less than $10,000" ~ "<50k",
        faminc == "Less than $10,000" ~ "<50k",
        faminc == "$10,000 - $19,999" ~ "<50k",
        faminc == "$20,000 - $29,999" ~ "<50k",
        faminc == "$30,000 - $39,999" ~ "<50k",
        faminc == "$40,000 - $49,999" ~ "<50k",
        faminc == "$50,000 - $59,999" ~ "50-100k",
        faminc == "$60,000 - $69,999" ~ "50-100k",
        faminc == "$70,000 - $79,999" ~ "50-100k",
        faminc == "$80,000 - $99,999" ~ "50-100k",
        faminc == "$100,000 - $119,999" ~ "100-150k",
        faminc == "$120,000 - $149,999" ~ "100-150k",
        faminc == "$150,000 - $199,999" ~ ">150k",
        faminc == "$200,000 - $249,999" ~ ">150k",
        faminc == "$250,000 - $349,999" ~ ">150k",
        faminc == "$350,000 - $499,999" ~ ">150k",
        faminc == "$500,000 or more" ~ ">150k",
        faminc == "$150,000 or more" ~ ">150k",
        #4788 do not say
        faminc == "Prefer not to say" ~ "Prefer not to say",
        #9 NAs
        TRUE ~ "Prefer not to say"),
        levels = c("<50k", "50-100k", "100-150k",">150k", "Prefer not to say")),
    
    recode_income_3way = factor(case_when(
        faminc == "Less than $10,000" ~ "<50k",
        faminc == "Less than $10,000" ~ "<50k",
        faminc == "$10,000 - $19,999" ~ "<50k",
        faminc == "$20,000 - $29,999" ~ "<50k",
        faminc == "$30,000 - $39,999" ~ "<50k",
        faminc == "$40,000 - $49,999" ~ "<50k",
        faminc == "$50,000 - $59,999" ~ "50-150k",
        faminc == "$60,000 - $69,999" ~ "50-150k",
        faminc == "$70,000 - $79,999" ~ "50-150k",
        faminc == "$80,000 - $99,999" ~ "50-150k",
        faminc == "$100,000 - $119,999" ~ "50-150k",
        faminc == "$120,000 - $149,999" ~ "50-150k",
        faminc == "$150,000 - $199,999" ~ ">150k",
        faminc == "$200,000 - $249,999" ~ ">150k",
        faminc == "$250,000 - $349,999" ~ ">150k",
        faminc == "$350,000 - $499,999" ~ ">150k",
        faminc == "$500,000 or more" ~ ">150k",
        faminc == "$150,000 or more" ~ ">150k",
        #4788 do not say
        faminc == "Prefer not to say" ~ "Prefer not to say",
        #9 NAs
        TRUE ~ "Prefer not to say"),
        levels = c("<50k", "50-150k",">150k","Prefer not to say"))
)

cces <- cces %>% mutate(
    recode_relig_6way = factor(case_when(
        recode_relig == "Protestant" ~ "Protestant",
        recode_relig == "Jewish" ~ "Jewish",
        recode_relig == "Catholic" ~ "Catholic",
        recode_relig == "Nothing in particular" ~ "Not religious",
        recode_relig == "Agnostic" ~ "Not religious",
        recode_relig == "Atheist" ~ "Not religious",
        recode_relig == "Buddhist" ~ "Other",
        recode_relig == "Muslim" ~ "Other",
        recode_relig == "Hindu" ~ "Other",
        recode_relig == "Mormon" ~ "Other",
        recode_relig == "Orthodox" ~ "Other",
        #there's a lot of these (2k) 
        #quite of few of them are v religious going to church weekly (300)
        #when asked what type of protestant many identified so call them protestant
        recode_relig == "Something else" & !is.na(religpew_protestant) ~ "Protestant",
        recode_relig == "Something else" & !is.na(religpew_protestant) ~ "Other",
        #38 NAs to not religious
        TRUE ~ "Not religious"),
        levels = c("Protestant", "Jewish", "Catholic","Not religious", "Other")),
    
    recode_attndch_4way = factor(case_when(
        recode_attndch == "Once a week" ~ "Weekly",
        recode_attndch == "More than once a week" ~ "Weekly",
        recode_attndch == "Once or twice a month" ~ "Monthly",
        recode_attndch == "A few times a year" ~ "Yearly",
        recode_attndch == "Seldom" ~ "Never",
        recode_attndch == "Never" ~ "Never",
        recode_attndch == "Don't know" ~ "Never",
        #the 35 who refused say never go
        TRUE ~ "Never"),
        levels = c("Weekly","Monthly", "Yearly", "Never"))
)


####Interactions
cat("Creating a few interactions in CCES\n")
cces = cces %>% mutate(
    #12 levels
    recode_pid_race = as.factor(paste(recode_pid_3way, 
                                      recode_race, sep = ", ")),
    #9 levels
    recode_educ_pid = as.factor(paste(recode_educ_3way,
                                      recode_pid_3way, sep = ", ")),
    #36 levels
    recode_educ_pid_race = as.factor(paste(recode_educ_3way,
                                           recode_pid_3way, 
                                           recode_race, sep = ", ")),
    ## education 3 way bucket among white voters, no split among non-white
    recode_educ_wh_3way = factor(case_when(
        recode_race != "White" ~ "No Split",
        recode_educ %in% c("No HS", "High school graduate", "Some college") ~ "No College",
        recode_educ %in% c("2-year", "4-year") ~ "College",
        TRUE ~ "Post-grad"), 
        levels = c("No Split", "No College", "College", "Post-grad")),
    
    recode_race_reg_wh_educ = as.factor(paste(recode_region, 
                                           case_when(recode_race == "White" ~ paste(recode_race, recode_educ, sep = ", "),
                                                                        TRUE ~ recode_race), sep = ", ")), 
    recode_race_educ_reg = as.factor(paste(recode_race, recode_region, 
                                           recode_educ_3way, sep = ", ")),
    #this has 4 levels
    recode_midwest_wh_edu = factor(case_when(
        (recode_race != "White" | recode_region != "Midwest") ~ "No Split", 
        TRUE ~ as.character(recode_educ_3way)), 
        levels = c("No Split", "No College", "College", "Post-grad")),
    #this has 13 levels
    recode_midwest_edu_race = as.factor(
        case_when(recode_region == "Midwest" ~ paste(recode_region,
                                                     recode_race,
                                                     recode_educ_3way, 
                                                     sep = ", "),
                  TRUE ~ "No Split")),
   #consider each individual age a "bucket"
    recode_age_factor = factor(case_when(recode_age <92 ~ as.character(recode_age), 
                                         TRUE ~ "92+")) 
)

#drop 360 respondents who do not report their vote choice
cat("Dropping CCES respondents who did not report their vote choice\n")
cces <- cces %>%
    filter(!is.na(recode_vote_2016)) %>% 
    #readjusting weights so they sum to N
    mutate(commonweight_vv_post = commonweight_vv_post/ mean(commonweight_vv_post))

    



################################################################################
#### 4. MODELING 3WAY VOTE CHOCIE OUTCOME  ######################################
################################################################################

cat("Creating 3way modeled vote choice outcome in CCES for application\n")

## Model outcome: 3way multinomial regularized with a lasso penalty and with population weights from CCES
#Train Model: train percent 80%
percent <- sample.int(n = nrow(cces),  size = round((nrow(cces)/10)*8))
stack_data_train <-  data.frame(bind_rows(pew, cces[percent,]), 
                                S = c(rep(1, nrow(pew)), rep(0, nrow(cces[percent,]))))
#base model
C_model = as.formula(~ recode_female + 
                         recode_pid_3way + 
                         recode_educ +
                         recode_region +
                         recode_income_5way + 
                         recode_relig_6way + 
                         recode_born + 
                         recode_race +
                         recode_attndch_4way +
                         recode_age +
                         I(recode_age^2) + 
                         recode_female:recode_pid_3way + 
                         recode_age:recode_pid_3way)

mod <- model.matrix(C_model, data = stack_data_train)

################ CCES Training Model ################
#lambda from 10 fold default CV
cat("Fitting Lasso Regularized multinomial using 80% CCES data as training set\n")
lasso_lambda <- cv.glmnet(x= mod[stack_data_train$S == 0,-1], 
                          y = cces$recode_vote_2016[percent],
                          alpha = 1,
                          family = "multinomial",
                          weights = cces$commonweight_vv_post[percent],
                          intercept = TRUE)

lasso_cces <- glmnet(x= mod[stack_data_train$S == 0, -1], 
                     y = cces$recode_vote_2016[percent],
                     alpha = 1,
                     lambda = lasso_lambda$lambda.1se,
                     family = "multinomial",
                     weights = cces$commonweight_vv_post[percent],
                     intercept = TRUE)


###### Eval Performance on Test Set ##############
cat("Evaluting model fit on remaining 20% CCES test set\n")
stack_data_test <- data.frame(bind_rows(pew, cces[-percent,]), 
                              S = c(rep(1, nrow(pew)), rep(0, nrow(cces[-percent,]))))

mod_test <- model.matrix(C_model, data = stack_data_test)
lasso_test = predict(lasso_cces,
                     s = lasso_lambda$lambda.1se,
                     type = "response",
                     newx = mod_test[stack_data_test$S==0,-1])

#designate predicted vote choice as max prob outcome to check accuracy
mod_test_plur <- apply(lasso_test,1, which.max)
mod_test_out <- factor(case_when(mod_test_plur-1 == 2 ~ 
                                     "Republican",
                                 mod_test_plur-1 == 1 ~ "Other", 
                                 mod_test_plur-1 == 0 ~ "Democrat"),
                       levels = c("Democrat", "Other", "Republican"))
#overall correct
cat("Percent Correct Overall: CCES model projected on CCES (Testing Data)\n")
sum(mod_test_out == cces$recode_vote_2016[-percent])/length(cces$recode_vote_2016[-percent])
conf_cces_test <- confusionMatrix(as.factor(mod_test_out),
                                  as.factor(cces$recode_vote_2016[-percent]))
cat("Confusion Matrix: CCES model projected on CCES (Testing Data)\n [Table C.2 in Appendix C]")
print(t(t(conf_cces_test$table)/colSums(conf_cces_test$table)))

################### Model on Pew ############
stack_data_all <- data.frame(bind_rows(pew, cces), 
                             S = c(rep(1, nrow(pew)), rep(0, nrow(cces))))

mod_all <- model.matrix(C_model, data = stack_data_all)

########## Create and Save 3way Outcome from projection of CCES Model #########
## Project CCES model on Pew
cat("Projecting CCES Model onto Pew data\n")
lasso_cces_on_pew = predict(lasso_cces,
                            s = lasso_lambda$lambda.1se,
                            type = "response",
                            newx = mod_all[stack_data_all$S==1,-1])

## Project CCES model on full CCES data
cat("Projecting CCES Model onto full CCES data\n")
lasso_cces_on_cces = predict(lasso_cces,
                             s = lasso_lambda$lambda.1se,
                             type = "response",
                             newx = mod_all[stack_data_all$S==0,-1])

#save model's predicted probabilities of voting for each of the three vote choices
#Create outcome variable for application in Section 6 as:
#the difference of the probability of voting Dem - probability of voting Rep 
cat("Creating vote difference outcome p(D) - p(R) in Pew and CCES\n")
cces <- cces %>%
    mutate(mod_cces_on_cces_pD = lasso_cces_on_cces[,"Democrat",], 
           mod_cces_on_cces_pR = lasso_cces_on_cces[,"Republican",],
           mod_cces_on_cces_pO = lasso_cces_on_cces[,"Other",]) %>% 
    #main outcome variable used in the application created here as "diff_cces_on_cces"
    mutate(diff_cces_on_cces = mod_cces_on_cces_pD -
               mod_cces_on_cces_pR, 
           margin_cces_on_cces = (mod_cces_on_cces_pD -
                                      mod_cces_on_cces_pR)/( 
                                          mod_cces_on_cces_pD +
                                              mod_cces_on_cces_pR)
    )
pew <- pew %>%
    mutate(mod_cces_on_pew_pD = lasso_cces_on_pew[,"Democrat",], 
           mod_cces_on_pew_pR = lasso_cces_on_pew[,"Republican",], 
           mod_cces_on_pew_pO = lasso_cces_on_pew[,"Other",]) %>% 
    #main outcome variable used in the application created here as "diff_cces_on_pew"
    mutate(diff_cces_on_pew = mod_cces_on_pew_pD - mod_cces_on_pew_pR, 
           margin_cces_on_pew = (mod_cces_on_pew_pD -
                                     mod_cces_on_pew_pR)/( 
                                         mod_cces_on_pew_pD +
                                             mod_cces_on_pew_pR)
    )

cat(paste0("Saving recoded CCES and Pew data as \"cces_prepared.rds\" and \"pew_prepared.rds\" in: \n\t", working_path, "/generated_data"))
saveRDS(cces, file = paste0(working_path, "/generated_data/", "cces_prepared.rds"))
saveRDS(pew, file = paste0(working_path, "/generated_data/","pew_prepared.rds"))
sink()
file.show("data_prep_out.Rout")
