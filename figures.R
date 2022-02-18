

# title: "lifetime_smoking_epigenetic_aging"

# author: "Eric Klopack"



# NOTE: sample variance-covariance matrices (sample.cov.g), 
# sample means (sample.mean.g), and asymptotic variance-covariance matrix 
# (Gamma.g) were produced from publicly available data using the following 
# code:



# vbsdesign <- svydesign(id = ~secu,
#                        weights = ~weight_all,
#                        strata = ~stratum,
#                        data = mrs_weight,
#                        nest = T)
# 
# for(i in clocks_list_adj) {
#   
#   # Names of the observed variables
#   ov.names <- c("mort20", "cancr1314", "hibp1314", "lung1314", "heart1314", 
#                 i, "py_gt_18", "c_smoke", "par_smoke_1", 
#                 "par_smoke_2", "age", "race_black", "race_hispanic", 
#                 "race_other", "gender_female", "education_0_11", "education_12", 
#                 "education_13_15", "log_wealth", "bmi_over", "bmi_obese1", 
#                 "bmi_obese2", "binge_1_4", "binge_5_p")
#   
#   # The MP-inverse duplication matrix is handy for removing redundancy
#   Dplus <- lavaan::lav_matrix_duplication_ginv(length(ov.names))
#   # Create a formula that includes all observed variables for svymean
#   ov.formula <- as.formula(paste("~", paste(ov.names, collapse="+")))
#   
#   
#   sample.cov.g <- as.matrix(svyvar(ov.formula, design=vbsdesign, na.rm=TRUE))
#   Gamma.cov.g <- attr(sample.cov.g, "var")
#   Gamma.cov.g <- Dplus %*% Gamma.cov.g %*% t(Dplus)
#   
#   sample.mean.g <- svymean(ov.formula, design=vbsdesign, na.rm=TRUE)
#   Gamma.mean.g <- attr(sample.mean.g, "var")
#   
#   Gamma.g <- lavaan::lav_matrix_bdiag(Gamma.mean.g, Gamma.cov.g)
#   Gamma.g <- Gamma.g * 2978
#   
#   attr(sample.cov.g, "var") <- NULL
#   tmp  <- as.vector(sample.mean.g)
#   names(tmp) <- names(sample.mean.g)
#   sample.mean.g <- tmp
# }



library(lavaan)

clocks_list_adj <- c('horvath_adj', 'hannum_adj', 'levine_adj', 'grim_adj',
                     'mpoa38_adj')

for (i in clocks_list_adj) {
  
  model_sem <- paste0(
    '
mort20 ~
kmort*', i, ' +
ymort*py_gt_18 + 
smort*c_smoke +
c1mort*par_smoke_1 + c2mort*par_smoke_2 +
age + race_black + race_hispanic + race_other + gender_female + 
education_0_11 + education_12 + education_13_15 + log_wealth +
bmi_over + bmi_obese1 + bmi_obese2 +
binge_1_4 + binge_5_p

cancr1314 ~
kcancr*', i, ' +
ycancr*py_gt_18 + 
scancr*c_smoke +
c1cancr*par_smoke_1 + c2cancr*par_smoke_2 +
age + race_black + race_hispanic + race_other + gender_female + 
education_0_11 + education_12 + education_13_15 + log_wealth +
bmi_over + bmi_obese1 + bmi_obese2 +
binge_1_4 + binge_5_p

hibp1314 ~
khibp*', i, ' +
yhibp*py_gt_18 + 
shibp*c_smoke +
c1hibp*par_smoke_1 + c2hibp*par_smoke_2 +
age + race_black + race_hispanic + race_other + gender_female + 
education_0_11 + education_12 + education_13_15 + log_wealth +
bmi_over + bmi_obese1 + bmi_obese2 +
binge_1_4 + binge_5_p

lung1314 ~
klung*', i, ' +
ylung*py_gt_18 + 
slung*c_smoke +
c1lung*par_smoke_1 + c2lung*par_smoke_2 +
age + race_black + race_hispanic + race_other + gender_female + 
education_0_11 + education_12 + education_13_15 + log_wealth +
bmi_over + bmi_obese1 + bmi_obese2 +
binge_1_4 + binge_5_p

heart1314 ~
kheart*', i, ' +
yheart*py_gt_18 + 
sheart*c_smoke +
c1heart*par_smoke_1 + c2heart*par_smoke_2 +
age + race_black + race_hispanic + race_other + gender_female + 
education_0_11 + education_12 + education_13_15 + log_wealth +
bmi_over + bmi_obese1 + bmi_obese2 +
binge_1_4 + binge_5_p



', i, ' ~ 
p*py_gt_18 + 
g*c_smoke +
m1*par_smoke_1 + m2*par_smoke_2

py_gt_18 ~ o*c_smoke +
b1*par_smoke_1 + b2*par_smoke_2 +
age + race_black + race_hispanic + race_other + gender_female + 
education_0_11 + education_12 + education_13_15 + log_wealth +
bmi_over + bmi_obese1 + bmi_obese2 +
binge_1_4 + binge_5_p

c_smoke ~
r1*par_smoke_1 + r2*par_smoke_2 +
age + race_black + race_hispanic + race_other + gender_female + 
education_0_11 + education_12 + education_13_15 + log_wealth +
bmi_over + bmi_obese1 + bmi_obese2 +
binge_1_4 + binge_5_p

')
  
  sample.mean.g <- readRDS(paste0('data/mean_', i, '.rds'))
  sample.cov.g <- readRDS(paste0('data/cov_', i, '.rds'))
  Gamma.g <- readRDS(paste0('data/nacov_', i, '.rds'))
  
  model <- sem(model_sem, sample.mean = sample.mean.g, 
    sample.cov = sample.cov.g, NACOV = Gamma.g, sample.nobs = 2978, 
    estimator = 'MLM', ordered = c('c_smoke', 'cancr1314', 'hibp1314', 
    'lung1314', 'heart1314', 'mort20'))
  
  summary(model, standardized = T)
  
}
