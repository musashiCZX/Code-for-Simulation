
options(tinytex.verbose = TRUE)
library(metafor)#library(xlsx) # did not load on my computer
library(tidyverse)
library(metaviz)
library(lme4)
library(nlme)
library(reshape2)

ctrl <- lmeControl(opt="optim", msMaxIter=100)




library(tidyverse)
library(copula)
library(MASS)






set.seed(34532)


#----------------------------------------------------------------------------------------------------------------------------
#                       Data Generation Machine  
#----------------------------------------------------------------------------------------------------------------------------



Simulation <- function(nstudy = 6, ngroup = as.data.frame(matrix(data = replicate(nstudy, 10), 
                                                                 nrow = nstudy, ncol = 2)), 
                       beta0 = replicate(nstudy, 40), treatment_effect, baseline_effect, 
                       tau, sigma_ik = rep(1,nstudy),
                       sigma_YB = rep(4,nstudy), imbalance = replicate(nstudy, 0)){
  #nstudy denotes the number of studies
  #ngroup denotes the number of patients in each study, which should be a dataframe with 2 columns.
  #And the row number is equal to the nstudy.
  #beta0 denotes the mean baseline score for the patients in the control group.
  #treatment_effect denotes the treatment effect which we are interested in.
  #study effect is the random effect which allows the various baseline scores among groups.
  #baseline_effect is the effect of the baseline scores.
  #The beta0 can combine with the study effect.
  #Imbalance denotes the mean difference between the groups.
  
  ntot <- sum(ngroup)
  study <- vector()
  treatment <- vector()
  data.meta <- data.frame(study,treatment)
  #We allow setting the correlation manually. If nothing is input, we generate a random correlation dataset.
  treatment_effect_real <- treatment_effect + rnorm(nstudy, 0, tau)
  
  #We generate the data for each study.
  for (i in 1:nstudy){
    
    #The study Index.
    study <- rep(i,sum(ngroup[i,]))
    #The treatment Index
    treatment <- c(rep(1,ngroup[i,1]),rep(0,ngroup[i,2]))
    YB_T <- rnorm(ngroup[i,1], mean = beta0[i] + imbalance[i], sd = sigma_YB[i])
    YB_C <- rnorm(ngroup[i,2], mean = beta0[i], sd = sigma_YB[i])
    sd_baseline <- sigma_YB[i]
    #sd_follow_up <- treatment_effect_real[i]*sd_baseline/correlation[i] #We do not need this. Remove.
    j <- 1
    #j <- 1 is for the treatment group.
    residual_1 <- rnorm(ngroup[i,1], 0, sigma_ik[1,i])
    
    #We calculate the mean, sd. Then apply the quantile function on u1 and u2 to make them two correlated normal sets. 
    YF_1 <- beta0[i] + baseline_effect*(YB_T) + treatment_effect_real[i]*j + residual_1#rnorm(ngroup[i,1], 0, sigma_ik[1,i])# - mean(YB[1:(ngroup[i,1] + ngroup[i,2])])) 
    #beta0 should be changed to beta0i.
   # print(residual_1)
    #rtemp1 <- mvrnorm(n = ngroup[i,1], mu = c(mean_baseline_i, mean_follow_up_1),
    #                   Sigma = matrix(data = c(sd_baseline^2, 
    #                                            correlation[i]*sd_baseline*sd_follow_up , 
    #                                            correlation[i]*sd_baseline*sd_follow_up,
    #                                            sd_follow_up^2), 2, 2), empirical = TRUE)
    
    #YB_temp1 <- #qnorm(u1[,1], mean = mean_baseline_i, sd = sd_baseline)
    #YF_temp1 <- #qnorm(u1[,2], mean = mean_follow_up, sd = sd_follow_up)
    #temp1 <- cbind(YB_temp1, YF_temp1)
    residual_2 <- rnorm(ngroup[i,2], 0, sigma_ik[2,i])
    
    j <- 0
    #j <- 0 is for the control group.
    YF_0 <- beta0[i] + baseline_effect*(YB_C) + treatment_effect_real[i]*j + residual_2#rnorm(ngroup[i,2], 0, sigma_ik[2,i]) #- 
    # mean(YB[1:(ngroup[i,1] + ngroup[i,2])]))
  #  print(residual_2)
    Y <- cbind(c(YB_T, YB_C), c(YF_1, YF_0))
    
    
    
    #rtemp0 <- mvrnorm(n = ngroup[i,2], mu = c(mean_baseline_i, mean_follow_up_0),
    #                  Sigma = matrix(data = c(sd_baseline^2, 
    #                                          correlation[i]*sd_baseline*sd_follow_up , 
    #                                          correlation[i]*sd_baseline*sd_follow_up,
    #                                          sd_follow_up^2), 2, 2), empirical = TRUE)
    
    #YB_temp0 <- #qnorm(u2[,1], mean = mean_baseline_i, sd = sd_baseline)
    #YF_temp0 <- #qnorm(u2[,2], mean = mean_follow_up, sd = sd_follow_up)
    #temp0 <- cbind(YB_temp0, YF_temp0)
    
    
    
    data.meta <- rbind(data.meta, cbind(study, treatment, Y))
  }
  names(data.meta) <- c('study', 'treatment', 'YB', 'YF')
  data.meta$YChange <- data.meta$YF - data.meta$YB
  #data.meta$study <- as.factor(data.meta$study)
  #data.meta$treatment <- as.factor(data.meta$treatment)
  #treatment <- rep(c(rep(1,ngroup),rep(0,ngroup)), nstudy)
  
  #data.meta$YB <- rnorm(ntot, mean_baseline, sigma_YB)
  
  return(list(simulated_data = data.meta,
              ngroup = ngroup,
              nstudy = nstudy,
              treatment = treatment_effect,
              treatment_effect_real = treatment_effect_real,
              tau = tau,
              sigma_ik = sigma_ik))
  
}





#----------------------------------------------------------------------------------------------------------------------------
#                       Perform meta-analysis models based on the original data
#----------------------------------------------------------------------------------------------------------------------------



simulation_process <- function(data, data.meta1){

results_list <- as.data.frame(matrix(nrow = 5, ncol = 6))
row.names(results_list) <- c('change score model', 
              'final score model', 
              'Trowman Method', 
              'pseudo base model', 
              'pseudo full model')
names(results_list) <- c('estimate', 'se', 'ci.lb', 'ci.ub', 'tau', 'abs error with the true value')
#Transform the data to a wide style.
dat_agency <- data[data$group == 0,]

names(dat_agency) <- paste0(names(dat_agency),'_',0)
dat <- cbind(data[data$group == 1,],dat_agency)





#Calculate the change values.
dat$ChangeTreatment <- dat$MeanBaseline - dat$MeanPostBaseline

dat$ChangeControl <- dat$MeanBaseline_0 - dat$MeanPostBaseline_0

dat$ChangeSdTreatment <- sqrt((dat$sdBaseline)^2 + 
                                (dat$sdPostBaseline)^2 - 
                                2*dat$sdBaseline*dat$sdPostBaseline*dat$Correlation)

dat$ChangeSdControl <- sqrt((dat$sdPostBaseline_0)^2 + 
                              (dat$sdBaseline_0)^2 - 
                              2*dat$sdPostBaseline_0*dat$sdBaseline_0*dat$Correlation_0)

#print(dat)





#Fit the random effect meta-analysis model. The final model and the Change model.
meta_Follow <- escalc(measure = 'MD',
                      m1i = MeanPostBaseline, m2i = MeanPostBaseline_0, 
                      sd1i = sdPostBaseline, sd2i = sdPostBaseline_0, 
                      n1i = NCFB, n2i = NCFB_0, data = dat)
meta_Change <- escalc(measure = 'MD', 
                      m1i = ChangeControl, m2i = ChangeTreatment,
                      sd1i = ChangeSdControl, sd2i = ChangeSdTreatment,
                      n1i = NCFB_0, n2i = NCFB, data = dat)


resFollow <- rma(measure = 'MD',
                 m1i = MeanPostBaseline, m2i = MeanPostBaseline_0, 
                 sd1i = sdPostBaseline, sd2i = sdPostBaseline_0, 
                 n1i = NCFB, n2i = NCFB_0, data = dat,verbose=F, digits=5, control=list(maxiter=1000,stepadj=0.5))

resChange <- rma(measure = 'MD', 
                 m1i = MeanChange, m2i = MeanChange_0,
                 sd1i = sdChange, sd2i = sdChange_0,
                 n1i = NCFB, n2i = NCFB_0, data = dat,verbose=F, digits=5, control=list(maxiter=1000,stepadj=0.5))





#Make the forest plot to visualize the results.
forest(resChange,slab = dat$Study,
       xlab = paste('The effect size, weight, and confidence interval',
                    '\n','with respect to each study are listed on the right side.'),
       col = 2, addcred = T, 
       showweights = T, header = T, width = 1, cex = .7,
       main = 'The Forest plot of the Change scores model')





#List the estimate values of the change model.

results_list['change score model',] <- c(resChange$b, 
                                         resChange$se, 
                                         resChange$ci.lb, 
                                         resChange$ci.ub, 
                                         sqrt(resChange$tau2), 
                                         (resChange$b - data.meta1$treatment)^2)




#Identical with the change model. We apply the same execution to final model.
forest(resFollow, slab = dat$Study, 
       xlab = paste('Identical with the last plot. The right side is',
                    '\n', 'the effect size, weight, and confidence interval respectively.'),       
       col = 2, addcred = T, 
       showweights = T, header = T, width = 1, cex = .7,
       main = 'The forest plot of final scores model')

results_list['final score model',] = c(resFollow$b, 
                             resFollow$se, 
                             resFollow$ci.lb, 
                             resFollow$ci.ub, 
                             sqrt(resFollow$tau2), 
                             (resFollow$b - data.meta1$treatment)^2)

#The Recovered ANCOVA




dat$sd_pooling_hat_baseline <- sqrt(((dat$NCFB - 1)*(dat$sdBaseline)^2 + (dat$NCFB_0 - 1)*(dat$sdBaseline_0)^2)/(dat$NCFB + dat$NCFB_0 - 2))
dat$sd_pooling_hat_postbaseline <- sqrt(((dat$NCFB - 1)*(dat$sdPostBaseline)^2 + (dat$NCFB_0 - 1)*(dat$sdPostBaseline_0)^2)/(dat$NCFB + dat$NCFB_0 - 2))
dat$rpooling <- with(dat, ((NCFB_0 - 1)*Correlation_0*sdBaseline_0*sdPostBaseline_0 + (NCFB - 1)*Correlation*sdBaseline*sdPostBaseline)/((NCFB + NCFB_0 - 2)*sd_pooling_hat_baseline*sd_pooling_hat_postbaseline))
dat$beta_hat <- dat$rpooling*dat$sd_pooling_hat_postbaseline/dat$sd_pooling_hat_baseline

dat$change_ANCOVA <- dat$MeanPostBaseline - dat$MeanPostBaseline_0 - (dat$MeanBaseline - dat$MeanBaseline_0)*dat$beta_hat

dat$sdPostBaseline_ANCOVA <- sqrt(((dat$sd_pooling_hat_postbaseline)^2)*(1/dat$NCFB + 1/dat$NCFB_0)*(1 - dat$rpooling^2))

resANCOVA <- rma(yi = change_ANCOVA, sei = sdPostBaseline_ANCOVA, 
                 data = dat,verbose=F, digits=5, control=list(maxiter=1000,stepadj=0.5))
forest(resANCOVA, slab = dat$Study, 
       xlab = paste('Identical with the last plot. The right side is',
                    '\n', 'the effect size, weight, and confidence interval respectively.'),       
       col = 2, addcred = T, 
       showweights = T, header = T, width = 1, cex = 0.7,
       main = 'The forest plot of Recovered ANCOVA model')
results_list['Recovered ANCOVA',] = c(resANCOVA$b, 
                                     resANCOVA$se, 
                                     resANCOVA$ci.lb, 
                                     resANCOVA$ci.ub, 
                                     sqrt(resANCOVA$tau2), 
                                     (resANCOVA$b - data.meta1$treatment)^2)



#Fit by lm and weight by size of each group.
Trowman <- lm(MeanPostBaseline ~ group + MeanBaseline,
              data = data, weights = NCFB)
estimate = summary(Trowman)$coefficients['group','Estimate'] 
     se = summary(Trowman)$coefficients['group','Std. Error']
results_list['Trowman Method',] = c(estimate, 
     se,
     estimate - 1.96*se,
     estimate + 1.96*se,
     NA,
     (estimate - data.meta1$treatment)^2)
#We have to  calculate the CI manually.  
 
  



#The simulation process of the pseudo IPD. 
#We use the mean, standard deviation, correlation, and sample size
#To simulate the pseudo data for each group in each study.
Simulation_IPD <- function(meanB, sdB, meanF, sdF, Correlation, n,group,name){
  Yi1 <- rnorm(n)
  Yi2 <- rnorm(n)
  Yi1 <- (Yi1 - mean(Yi1))/sd(Yi1)
  Yi2 <- (Yi2 - mean(Yi2))/sd(Yi2)
  r_star <- cor(Yi1,Yi2)
  model <- lm(Yi2~Yi1)
  beta_hat <- model$coefficient
  residual_hat <- model$residual
  Yi3 <- (Yi1 * Correlation) + residual_hat*sqrt(1-Correlation^2)/sqrt(1-r_star^2)
  YBi <- Yi1*sdB + meanB
  YFi <- Yi3*sdF + meanF
  result <- data.frame(Baseline = YBi,
                       
                       Follow_up = YFi,
                       group = group,
                       Study = name,ID = 1:n)

  return(result)
}





#Apply the simulation function on the data. Name it 'pseudo_IPD'.
pseudo_IPD <- NULL
ntot <- nrow(data)
for (i in 1:ntot) 
{ pseudo_IPD <- rbind(pseudo_IPD,  
        Simulation_IPD( meanB=data[i,'MeanBaseline'],
             sdB=data[i,'sdBaseline'],
             meanF= data[i,'MeanPostBaseline'], 
             sdF = data[i,'sdPostBaseline'],
             Correlation = data[i,'Correlation'],
             n= data[i,'NCFB'],
             group = data[i, 'group'],
             name =data[i,'Study']))


}
# Center the Baseline and 'group'
pseudo_IPD$MeanBaseline <- ave(pseudo_IPD$Baseline,pseudo_IPD$Study)
pseudo_IPD$Baseline_centered <- pseudo_IPD$Baseline - pseudo_IPD$MeanBaseline
pseudo_IPD$group_centered <- as.numeric(pseudo_IPD$group) - 0.5


original_IPD <- data.meta1$simulated_data
original_IPD$MeanBaseline <- ave(original_IPD$YB,original_IPD$study)
original_IPD$Baseline_centered <- original_IPD$YB - original_IPD$MeanBaseline
original_IPD$group_centered <- as.numeric(original_IPD$treatment) - 0.5
names(original_IPD) <- c("Study", 
                         'group', 
                         'Baseline', 
                         'Follow_up',
                         'Change_Score', 
                         'MeanBaseline',
                         'Baseline_centered',
                         'group_centered')

# Utilize the 'aggregate' function to extract the information of the pseudo IPD.
# Check the mean and standard deviationof y1 and y2, and correlation y1, y2
check <- cbind(aggregate(Baseline~group+Study, data=pseudo_IPD, mean), 
              aggregate(Follow_up~group+Study, data=pseudo_IPD, mean)[3],
              aggregate(Baseline~group+Study, data=pseudo_IPD, sd)[3],
              aggregate(Follow_up~group+Study, data=pseudo_IPD, sd)[3],
              aggregate(ID~group+Study, data = pseudo_IPD,max)[3],
              as.vector(cbind(by(pseudo_IPD, pseudo_IPD[,c("group","Study")], 
                                 function(x) {cor(x$Baseline,x$Follow_up)}))))



names(check) <- c('group',
                  'Study', 
                  'MeanBaseline',
                  'MeanPostBaseline',
                  'sdBaseline',
                  'sdPostBaseline',
                  'NCFB',
                  'Correlation')








#Since the dataset 'check' is identical with the aggregate data.
#We use it in the standard meta-analysis.
check_agency <- check[check$group == 0,]

names(check_agency) <- paste0(names(check_agency),'_',0)
check_dat <- cbind(check[check$group == 1,],check_agency)





check_dat$ChangeTreatment <- 
  check_dat$MeanBaseline - check_dat$MeanPostBaseline

check_dat$ChangeControl <- 
  check_dat$MeanBaseline_0 - check_dat$MeanPostBaseline_0

check_dat$ChangeSdTreatment <- 
  sqrt((check_dat$sdBaseline)^2 + (check_dat$sdPostBaseline)^2 - 
         2*check_dat$sdBaseline*check_dat$sdPostBaseline*check_dat$Correlation)

check_dat$ChangeSdControl <- 
  sqrt((check_dat$sdPostBaseline_0)^2 + (check_dat$sdBaseline_0)^2 - 
         2*check_dat$sdPostBaseline_0*check_dat$sdBaseline_0*check_dat$Correlation_0)
#print(check_dat)




# Center the 'group'
#pseudo_IPD$group_centered <- as.numeric(pseudo_IPD$group) - 0.5




meta_pseudo_Follow <- escalc(measure = 'MD', 
                             m1i = MeanPostBaseline, 
                             m2i = MeanPostBaseline_0, 
                             sd1i = sdPostBaseline, 
                             sd2i = sdPostBaseline_0, 
                             n1i = NCFB, 
                             n2i = NCFB_0, 
                             data = check_dat)




meta_pseudo_change <- escalc(measure = 'MD', 
                             m1i = ChangeControl, 
                             m2i = ChangeTreatment, 
                             sd1i = ChangeSdControl, 
                             sd2i = ChangeSdTreatment, 
                             n1i = NCFB_0, 
                             n2i = NCFB, 
                             data = check_dat)


resFollow.pseudo <- rma(yi = yi, vi = vi, data = meta_pseudo_Follow,verbose=F, digits=5, control=list(maxiter=1000,stepadj=0.5))

resChange.pseudo <- rma(yi = yi, vi = vi, data = meta_pseudo_change,verbose=F, digits=5, control=list(maxiter=1000,stepadj=0.5),
                        knha = T)
forest(resChange.pseudo,slab = meta_pseudo_change$Study,
       main = 'Forest plot of the change score model for pseudo-IPD',
       xlab = 'The estimation is identical with the original data change score model',
       col = 2, addcred = T, 
       showweights = T, header = T, width = 1, cex = .8)




forest(resFollow.pseudo,slab = meta_pseudo_Follow$Study, 
       main = 'Forest plot of the final score model for pseudo-IPD',
       xlab = 'The estimation is identical with the original data final score model',
       col = 2, addcred = T, 
       showweights = T, header = T, width = 1, cex = .8)





#This is a function for extracting information from the lme.
#This function comes from Papadimitropoulou's code. 
groupeffect <- function(results)
{
  groupeff   <- summary(results)$tTable["group",]
# Add 95% CI extracted fom nlme using t-distribution
CI          <- intervals(results, which="fixed")
df          <- data.frame(CI$fixed)
names       <- c(names(groupeff),"95% CI lwb", "95% CI upb") 
groupeff    <- c(groupeff, df["group",]$lower, df["group",]$upper)
names(groupeff) <- names  
#print(groupeff)
# Estimated tau2  
tau2 <-  VarCorr((results)) [1,]
#print("Tau2 is")
#print(tau2)
return(as.numeric(c(groupeff, tau2)))
}



#The weighted Group Model. We regardless the information from the Baseline.
ctrl <- lmeControl(opt="optim", msMaxIter=100)
model_base <- lme(Follow_up ~ group, data = pseudo_IPD, 
                  random = ~-1 + group_centered|Study, 
                  control = ctrl,
                  weights = varIdent(form = ~1|Study))
model_base_results <- groupeffect(model_base)
results_list['pseudo base model',] <- c(model_base_results[c(1,2,6,7,9)], 
                             (model_base_results[1] - data.meta1$treatment)^2)


pseudo_IPD$arm <- paste0(pseudo_IPD$Study, '_', pseudo_IPD$group)
original_IPD$arm <- paste0(original_IPD$Study, '_', original_IPD$group)
#The Full model with the information of Baseline, Study, and group.
model_1 <- lme(Follow_up ~ Baseline_centered*as.factor(Study)
               + group, 
               data = pseudo_IPD,random = ~-1 + group_centered|Study, 
               control = ctrl)#,
               #weights = varIdent(form = ~Study|arm), method = 'REML')
model_full_results <- groupeffect(model_1)
results_list['pseudo full model',] <- c(model_full_results[c(1,2,6,7,9)], 
                                        (model_full_results[1] - data.meta1$treatment)^2)


model_original <- lme(Follow_up ~ Baseline_centered*as.factor(Study)
                      + group, 
                      data = original_IPD,random = ~-1 + group_centered|Study, 
                      control = ctrl)#, weights = varIdent(form = ~Study|arm), method = 'REML')
model_original_results <- groupeffect(model_original)
results_list['original full model',] <- c(model_original_results[c(1,2,6,7,9)], 
                                          (model_original_results[1] - data.meta1$treatment)^2)

#print(results_list['original full model',] - results_list['pseudo full model',])

Studies <- unique(pseudo_IPD$Study)
Two_stage_meta_analysis_data <- data.frame()
for (name_meta in Studies){
  data_study <- pseudo_IPD[pseudo_IPD$Study == name_meta,]
  lm_study <- lm(Follow_up ~ Baseline + group, data = data_study)
  record_study <- c(summary(lm_study)$coefficients['group', 'Estimate'],
                    summary(lm_study)$coefficients['group', 'Std. Error']^2,
                    nrow(data_study))
  Two_stage_meta_analysis_data <- rbind(Two_stage_meta_analysis_data, record_study)
}
Two_stage_meta_analysis_data <- cbind(Two_stage_meta_analysis_data, Studies)
names(Two_stage_meta_analysis_data) <- c('yi', 'vi', 'ni', 'studyi')

Two_stage_meta_analysis <- rma(yi = yi, vi = vi, slab = studyi,# random = ~ 1|studyi,
                                  data = Two_stage_meta_analysis_data,verbose=F, digits=5, 
                                  control=list(maxiter=1000,stepadj=0.5),
                                  knha = T)


results_list['pseudo two-stage model',] <-  c(Two_stage_meta_analysis$b, 
                                             Two_stage_meta_analysis$se, 
                                             Two_stage_meta_analysis$ci.lb, 
                                             Two_stage_meta_analysis$ci.ub, 
                                            sqrt(Two_stage_meta_analysis$tau2), 
                                            (Two_stage_meta_analysis$b - data.meta1$treatment)^2)



return(list(true_value = data.meta1$treatment,
            results = results_list))
}

#----------------------------------------------------------------------------------------------------------------------------
#                       Parameters tuning: We can vary the parameters globally  
#----------------------------------------------------------------------------------------------------------------------------

nstudy_exp <- 8
ngroup_exp <- 20
residual_exp <- replicate(nstudy_exp,c(16,4))
sd_baseline_exp <- 20
imbalance_exp <- 5



#----------------------------------------------------------------------------------------------------------------------------
#                       Prepare the blank data.frames to storage the simulation results
#----------------------------------------------------------------------------------------------------------------------------



## 
correlation_various_results <- vector(mode = 'list')
mean_baseline_results <- vector(mode = 'list')
baseline_effect_results <- vector(mode = 'list')
treatment_effect_results <- vector(mode = 'list')
baseline_effect_results <- vector(mode = 'list')



#Replicate the process and curve the results:

set.seed(3456)

#The candidate treatment effects.

Assess_Range <- c(50)
plot_number <- 1:length(Assess_Range)
Repeat_times <- 200

baseline_effect_results <- vector(mode = 'list')


##############################################The data.frame to store the error.

Change_Score_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Change_Score_modelaggregate_error) <- Assess_Range

Final_Score_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Final_Score_modelaggregate_error) <- Assess_Range

Recover_ANCOVA_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Recover_ANCOVA_modelaggregate_error) <- Assess_Range

Trowman_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Trowman_modelaggregate_error) <- Assess_Range

Pseudo_Base_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Pseudo_Base_modelaggregate_error) <- Assess_Range

Pseudo_Full_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Pseudo_Full_modelaggregate_error) <- Assess_Range

Pseudo_Two_Stage_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Pseudo_Two_Stage_modelaggregate_error) <- Assess_Range

Original_Full_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Original_Full_modelaggregate_error) <- Assess_Range

#########################################The data.frame to store the estimation


Change_Score_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Change_Score_modelaggregate_estimation) <- Assess_Range

Final_Score_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Final_Score_modelaggregate_estimation) <- Assess_Range

Recover_ANCOVA_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Recover_ANCOVA_modelaggregate_estimation) <- Assess_Range

Trowman_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Trowman_modelaggregate_estimation) <- Assess_Range

Pseudo_Base_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Pseudo_Base_modelaggregate_estimation) <- Assess_Range

Pseudo_Full_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Pseudo_Full_modelaggregate_estimation) <- Assess_Range

Pseudo_Two_Stage_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Pseudo_Two_Stage_modelaggregate_estimation) <- Assess_Range

Original_Full_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Original_Full_modelaggregate_estimation) <- Assess_Range

#########################################The tau


Change_Score_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Change_Score_modelaggregate_tau) <- Assess_Range

Final_Score_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Final_Score_modelaggregate_tau) <- Assess_Range

Recover_ANCOVA_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Recover_ANCOVA_modelaggregate_tau) <- Assess_Range

Trowman_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Trowman_modelaggregate_tau) <- Assess_Range

Pseudo_Base_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Pseudo_Base_modelaggregate_tau) <- Assess_Range

Pseudo_Full_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Pseudo_Full_modelaggregate_tau) <- Assess_Range

Pseudo_Two_Stage_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Pseudo_Two_Stage_modelaggregate_tau) <- Assess_Range

Original_Full_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Original_Full_modelaggregate_tau) <- Assess_Range

#######################################################The Sd

Change_Score_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Change_Score_modelaggregate_sd) <- Assess_Range

Final_Score_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Final_Score_modelaggregate_sd) <- Assess_Range

Recover_ANCOVA_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Recover_ANCOVA_modelaggregate_sd) <- Assess_Range

Trowman_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Trowman_modelaggregate_sd) <- Assess_Range

Pseudo_Base_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Pseudo_Base_modelaggregate_sd) <- Assess_Range

Pseudo_Full_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Pseudo_Full_modelaggregate_sd) <- Assess_Range

Pseudo_Two_Stage_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Pseudo_Two_Stage_modelaggregate_sd) <- Assess_Range

Original_Full_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range)))
names(Original_Full_modelaggregate_sd) <- Assess_Range


#----------------------------------------------------------------------------------------------------------------------------
#                       The main simulation experiment: We generate the original IPD and aggregate to the AD
#                       Then fit the meta-analysis models.
#----------------------------------------------------------------------------------------------------------------------------




for (item in 1:Repeat_times){
for (treatment_effect in Assess_Range){

data.meta1 <- Simulation(treatment_effect = -treatment_effect, 
                         baseline_effect = 0.5, tau = 13, sigma_ik = residual_exp,
                                               sigma_YB = replicate(nstudy, sd_baseline_exp), beta0 = replicate(nstudy, 40),
                                               nstudy = nstudy_exp, ngroup = as.data.frame(matrix(data = replicate(2*nstudy_exp, ngroup_exp), 
                                                                                              nrow = nstudy_exp, ncol = 2)))#, imbalance = replicate(nstudy, imbalance_exp))
#  Simulation(treatment_effect = treatment_effect, 
#                         baseline_effect = baseline_effect, tau = 4, sigma_ik = rep(1,8),
#                         sigma_YB = sigma_YB, mean_baseline = mean_baseline,correlation = correlation, 
#                         nstudy = nstudy, ngroup = ngroup)


aggregate_meta1 <- cbind(aggregate(YB~study+treatment, data=data.meta1$simulated_data, mean), 
                         aggregate(YB~study+treatment, data=data.meta1$simulated_data, sd)[3],
                         aggregate(YF~study+treatment, data=data.meta1$simulated_data, mean)[3],
                         
                         aggregate(YF~study+treatment, data=data.meta1$simulated_data, sd)[3], 
                         aggregate(YChange~study+treatment, data = data.meta1$simulated_data, mean)[3],
                         aggregate(YChange~study+treatment, data = data.meta1$simulated_data, sd)[3],
                         c(data.meta1$ngroup[,2], data.meta1$ngroup[,1]))
aggregate_meta1$group <- aggregate_meta1$treatment
ID <- aggregate_meta1$study
aggregate_meta1 <- cbind(ID, aggregate_meta1)
#aggregate_meta1$study <- as.factor(aggregate_meta1$study)

Correlation <- data.meta1$simulated_data %>%
  group_by(treatment, study) %>%
  summarise(Correlation = cor(YF, YB))

aggregate_meta1$Correlation <- Correlation$Correlation

aggregate_meta1$treatment <- NULL
names(aggregate_meta1) <- c("ID", 'Study', 
                            'MeanBaseline','sdBaseline',
                            'MeanPostBaseline','sdPostBaseline', 
                            'MeanChange', 'sdChange', 
                            'NCFB', 'group', 'Correlation')
data <- aggregate_meta1
results <- simulation_process(data, data.meta1)

treatment_effect_results[[paste0('treatment_effect = ',treatment_effect)]] <- results$results

#baseline_effect_results[[paste0('baseline_effect = ',baseline_effect)]] <- results$results

#mean_baseline_results[[paste0('mean_baseline = ',mean_baseline)]] <- results$results

#correlation_various_results[[paste0('correlation = ',correlation)]] <- results$results
}
  for (len in 1:length(treatment_effect_results)){
      Change_Score_modelaggregate_error[item,len] <- treatment_effect_results[[len]]['change score model','abs error with the true value']
      Final_Score_modelaggregate_error[item,len] <- treatment_effect_results[[len]]['final score model','abs error with the true value']
      Recover_ANCOVA_modelaggregate_error[item,len] <- treatment_effect_results[[len]]['Recovered ANCOVA','abs error with the true value']
      Trowman_modelaggregate_error[item,len] <- treatment_effect_results[[len]]['Trowman Method', 'abs error with the true value']
      Pseudo_Base_modelaggregate_error[item,len] <- treatment_effect_results[[len]]['pseudo base model', 'abs error with the true value']
      Pseudo_Full_modelaggregate_error[item,len] <- treatment_effect_results[[len]]['pseudo full model', 'abs error with the true value']
      Pseudo_Two_Stage_modelaggregate_error[item,len] <- treatment_effect_results[[len]]['pseudo two-stage model', 'abs error with the true value']
      Original_Full_modelaggregate_error[item,len] <- treatment_effect_results[[len]]['original full model', 'abs error with the true value']
      
      
      Change_Score_modelaggregate_estimation[item,len] <- treatment_effect_results[[len]]['change score model','estimate']
      Final_Score_modelaggregate_estimation[item,len] <- treatment_effect_results[[len]]['final score model','estimate']
      Recover_ANCOVA_modelaggregate_estimation[item,len] <- treatment_effect_results[[len]]['Recovered ANCOVA','estimate']
      Trowman_modelaggregate_estimation[item,len] <- treatment_effect_results[[len]]['Trowman Method', 'estimate']
      Pseudo_Base_modelaggregate_estimation[item,len] <- treatment_effect_results[[len]]['pseudo base model', 'estimate']
      Pseudo_Full_modelaggregate_estimation[item,len] <- treatment_effect_results[[len]]['pseudo full model', 'estimate']
      Pseudo_Two_Stage_modelaggregate_estimation[item,len] <- treatment_effect_results[[len]]['pseudo two-stage model', 'estimate']
      Original_Full_modelaggregate_estimation[item,len] <- treatment_effect_results[[len]]['original full model', 'estimate']
      
      Change_Score_modelaggregate_tau[item,len] <- treatment_effect_results[[len]]['change score model','tau']
      Final_Score_modelaggregate_tau[item,len] <- treatment_effect_results[[len]]['final score model','tau']
      Recover_ANCOVA_modelaggregate_tau[item,len] <- treatment_effect_results[[len]]['Recovered ANCOVA','tau']
      Trowman_modelaggregate_tau[item,len] <- treatment_effect_results[[len]]['Trowman Method', 'tau']
      Pseudo_Base_modelaggregate_tau[item,len] <- treatment_effect_results[[len]]['pseudo base model', 'tau']
      Pseudo_Full_modelaggregate_tau[item,len] <- treatment_effect_results[[len]]['pseudo full model', 'tau']
      Pseudo_Two_Stage_modelaggregate_tau[item,len] <- treatment_effect_results[[len]]['pseudo two-stage model', 'tau']
      Original_Full_modelaggregate_tau[item,len] <- treatment_effect_results[[len]]['original full model', 'tau']
      
      Change_Score_modelaggregate_sd[item,len] <- treatment_effect_results[[len]]['change score model','se']
      Final_Score_modelaggregate_sd[item,len] <- treatment_effect_results[[len]]['final score model','se']
      Recover_ANCOVA_modelaggregate_sd[item,len] <- treatment_effect_results[[len]]['Recovered ANCOVA','se']
      Trowman_modelaggregate_sd[item,len] <- treatment_effect_results[[len]]['Trowman Method', 'se']
      Pseudo_Base_modelaggregate_sd[item,len] <- treatment_effect_results[[len]]['pseudo base model', 'se']
      Pseudo_Full_modelaggregate_sd[item,len] <- treatment_effect_results[[len]]['pseudo full model', 'se']
      Pseudo_Two_Stage_modelaggregate_sd[item,len] <- treatment_effect_results[[len]]['pseudo two-stage model', 'se']
      Original_Full_modelaggregate_sd[item,len] <- treatment_effect_results[[len]]['original full model', 'se']
  }
  
}


#----------------------------------------------------------------------------------------------------------------------------
#                       Visualize the errors trend between models
#                       Visualize the estimates distribution with boxplot
#----------------------------------------------------------------------------------------------------------------------------


par(mfrow = c(1,1))
estimation_results <- vector(mode = 'list')
for (i in 1:length(Assess_Range)){
  estimation_results[[paste0('treatment_effect = ', Assess_Range[i])]] <- 
    cbind(Change_Score_modelaggregate_estimation[,i],
          Final_Score_modelaggregate_estimation[,i],
          Recover_ANCOVA_modelaggregate_estimation[,i],
          Trowman_modelaggregate_estimation[,i],
          Pseudo_Base_modelaggregate_estimation[,i],
          Pseudo_Full_modelaggregate_estimation[,i],
          Pseudo_Two_Stage_modelaggregate_estimation[,i],
          Original_Full_modelaggregate_estimation[,i])
        
  boxplot.matrix(as.matrix(estimation_results[[paste0('treatment_effect = ', Assess_Range[i])]]),
                 xaxt = 'n',
                 xlab = 'Models',
                 ylab = 'estimation distribution',
                 main = paste0('distribution of the estimation for treatment effect = ', 
                               Assess_Range[i]),
                 col = c('red','blue','yellow','green','brown','orange','purple','lightblue'))
  axis(1, at = 1:8, labels = c('change', 'final', 'RA', 'Trowman', 'PB', 'PF', 'PTS','OF'))
}
par(mfrow = c(1,1))
##Make a table with Mean, se, tau, ci, abs error with the true values. And the mean and sd of them.

for (plt in plot_number){
  plot(x = 1:Repeat_times, y = Change_Score_modelaggregate_error[,paste0(Assess_Range[plt])], type = 'b', pch = 20, col = 'red', xlab = 'index', ylab = 'errors')
  lines(x = 1:Repeat_times, y = Final_Score_modelaggregate_error[,paste0(Assess_Range[plt])], pch = 20, col = 'blue', type = 'b')
  lines(x = 1:Repeat_times, y = Recover_ANCOVA_modelaggregate_error[,paste0(Assess_Range[plt])], pch = 20, col = 'yellow', type = 'b')
  lines(x = 1:Repeat_times, y = Trowman_modelaggregate_error[,paste0(Assess_Range[plt])], pch = 20, col = 'green', type = 'b')
  lines(x = 1:Repeat_times, y = Pseudo_Base_modelaggregate_error[,paste0(Assess_Range[plt])], pch = 20, col = 'brown', type = 'b')
  lines(x = 1:Repeat_times, y = Pseudo_Full_modelaggregate_error[,paste0(Assess_Range[plt])], pch = 20, col = 'orange', type = 'b')
  lines(x = 1:Repeat_times, y = Pseudo_Two_Stage_modelaggregate_error[,paste0(Assess_Range[plt])], pch = 20, col = 'purple', type = 'b')
  lines(x = 1:Repeat_times, y = Original_Full_modelaggregate_error[,paste0(Assess_Range[plt])], pch = 20, col = 'lightblue', type = 'b')
  
  legend('topright', legend = c('change model', 
                                'final model',
                                'Recovered ANCOVA',
                                'Trowman Model',
                                'pseudo base model',
                                'pseudo full model',
                                'pseudo two-stage model',
                                'original full model'), 
         col = c('red','blue','yellow','green','brown','orange','purple','lightblue'), pch = 20)
}

results_treatments <- vector(mode = 'list')
blank <- as.data.frame(matrix(NA, nrow = 8, ncol = 4))
names(blank) <- c('estimation', 'sd', 'tau', 'error')
row.names(blank) <- c('change score model', 
                      'final score model',
                      'Trowman Method', 
                      'pseudo base model',
                      'pseudo full model', 
                      'Recovered ANCOVA',
                      'pseudo two-stage model',
                      'original full model')
for (treatment in Assess_Range){
  results_treatments[[paste0('treatment_effect = ',treatment)]] <- blank
}

for (treatment in Assess_Range){
  
  ##Add value to the final table:
  Change_Model_Error_Per_simulation <- apply(Change_Score_modelaggregate_error, 2, mean)
  Change_Model_Estimation_Per_simulation <- apply(Change_Score_modelaggregate_estimation, 2, mean)
  Change_Model_Tau_Per_simulation <- apply(Change_Score_modelaggregate_tau, 2, mean)
  Change_Model_Sd_Per_simulation <- apply(Change_Score_modelaggregate_sd, 2, mean)
  
  Final_Model_Error_Per_simulation <- apply(Final_Score_modelaggregate_error, 2, mean)
  Final_Model_Estimation_Per_simulation <- apply(Final_Score_modelaggregate_estimation, 2, mean)
  Final_Model_Tau_Per_simulation <- apply(Final_Score_modelaggregate_tau, 2, mean)
  Final_Model_Sd_Per_simulation <- apply(Final_Score_modelaggregate_sd, 2, mean)
  
  Recover_ANCOVA_Error_Per_simulation <- apply(Recover_ANCOVA_modelaggregate_error, 2, mean)
  Recover_ANCOVA_Estimation_Per_simulation <- apply(Recover_ANCOVA_modelaggregate_estimation, 2, mean)
  Recover_ANCOVA_Tau_Per_simulation <- apply(Recover_ANCOVA_modelaggregate_tau, 2, mean)
  Recover_ANCOVA_Sd_Per_simulation <- apply(Recover_ANCOVA_modelaggregate_sd, 2, mean)
  
  Trowman_Error_Per_simulation <- apply(Trowman_modelaggregate_error, 2, mean)
  Trowman_Estimation_Per_simulation <- apply(Trowman_modelaggregate_estimation, 2, mean)
  Trowman_Tau_Per_simulation <- apply(Trowman_modelaggregate_tau, 2, mean)
  Trowman_Sd_Per_simulation <- apply(Trowman_modelaggregate_sd, 2, mean)
  
  Pseudo_Base_Error_Per_simulation <- apply(Pseudo_Base_modelaggregate_error, 2, mean)
  Pseudo_Base_Estimation_Per_simulation <- apply(Pseudo_Base_modelaggregate_estimation, 2, mean)
  Pseudo_Base_Tau_Per_simulation <- apply(Pseudo_Base_modelaggregate_tau, 2, mean)
  Pseudo_Base_Sd_Per_simulation <- apply(Pseudo_Base_modelaggregate_sd, 2, mean)
  
  Pseudo_Full_Error_Per_simulation <- apply(Pseudo_Full_modelaggregate_error, 2, mean)
  Pseudo_Full_Estimation_Per_simulation <- apply(Pseudo_Full_modelaggregate_estimation, 2, mean)
  Pseudo_Full_Tau_Per_simulation <- apply(Pseudo_Full_modelaggregate_tau, 2, mean)
  Pseudo_Full_Sd_Per_simulation <- apply(Pseudo_Full_modelaggregate_sd, 2, mean)
  
  Pseudo_Two_Stage_Error_Per_simulation <- apply(Pseudo_Two_Stage_modelaggregate_error, 2, mean)
  Pseudo_Two_Stage_Estimation_Per_simulation <- apply(Pseudo_Two_Stage_modelaggregate_estimation, 2, mean)
  Pseudo_Two_Stage_Tau_Per_simulation <- apply(Pseudo_Two_Stage_modelaggregate_tau, 2, mean)
  Pseudo_Two_Stage_Sd_Per_simulation <- apply(Pseudo_Two_Stage_modelaggregate_sd, 2, mean)
  
  Original_Full_Error_Per_simulation <- apply(Original_Full_modelaggregate_error, 2, mean)
  Original_Full_Estimation_Per_simulation <- apply(Original_Full_modelaggregate_estimation, 2, mean)
  Original_Full_Tau_Per_simulation <- apply(Original_Full_modelaggregate_tau, 2, mean)
  Original_Full_Sd_Per_simulation <- apply(Original_Full_modelaggregate_sd, 2, mean)
  
  
  
  results_treatments[[paste0('treatment_effect = ',treatment)]]['change score model', 'estimation'] <- 
    Change_Model_Estimation_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['change score model', 'error'] <- 
    Change_Model_Error_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['change score model', 'tau'] <- 
    Change_Model_Tau_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['change score model', 'sd'] <- 
    Change_Model_Sd_Per_simulation[[paste0(treatment)]]
  
  results_treatments[[paste0('treatment_effect = ',treatment)]]['final score model', 'estimation'] <- 
    Final_Model_Estimation_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['final score model', 'error'] <- 
    Final_Model_Error_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['final score model', 'tau'] <- 
    Final_Model_Tau_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['final score model', 'sd'] <- 
    Final_Model_Sd_Per_simulation[[paste0(treatment)]]
  
  results_treatments[[paste0('treatment_effect = ',treatment)]]['Recovered ANCOVA', 'estimation'] <- 
    Recover_ANCOVA_Estimation_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['Recovered ANCOVA', 'error'] <- 
    Recover_ANCOVA_Error_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['Recovered ANCOVA', 'tau'] <- 
    Recover_ANCOVA_Tau_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['Recovered ANCOVA', 'sd'] <- 
    Recover_ANCOVA_Sd_Per_simulation[[paste0(treatment)]]
  
  results_treatments[[paste0('treatment_effect = ',treatment)]]['Trowman Method', 'estimation'] <- 
    Trowman_Estimation_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['Trowman Method', 'error'] <- 
    Trowman_Error_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['Trowman Method', 'tau'] <- 
    Trowman_Tau_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['Trowman Method', 'sd'] <- 
    Trowman_Sd_Per_simulation[[paste0(treatment)]]
  
  results_treatments[[paste0('treatment_effect = ',treatment)]]['pseudo base model', 'estimation'] <- 
    Pseudo_Base_Estimation_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['pseudo base model', 'error'] <- 
    Pseudo_Base_Error_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['pseudo base model', 'tau'] <- 
    Pseudo_Base_Tau_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['pseudo base model', 'sd'] <- 
    Pseudo_Base_Sd_Per_simulation[[paste0(treatment)]]
  
  results_treatments[[paste0('treatment_effect = ',treatment)]]['pseudo full model', 'estimation'] <- 
    Pseudo_Full_Estimation_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['pseudo full model', 'error'] <- 
    Pseudo_Full_Error_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['pseudo full model', 'tau'] <- 
    Pseudo_Full_Tau_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['pseudo full model', 'sd'] <- 
    Pseudo_Full_Sd_Per_simulation[[paste0(treatment)]]
  
  results_treatments[[paste0('treatment_effect = ',treatment)]]['pseudo two-stage model', 'estimation'] <- 
    Pseudo_Two_Stage_Estimation_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['pseudo two-stage model', 'error'] <- 
    Pseudo_Two_Stage_Error_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['pseudo two-stage model', 'tau'] <- 
    Pseudo_Two_Stage_Tau_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['pseudo two-stage model', 'sd'] <- 
    Pseudo_Two_Stage_Sd_Per_simulation[[paste0(treatment)]]
  
  results_treatments[[paste0('treatment_effect = ',treatment)]]['original full model', 'estimation'] <- 
    Original_Full_Estimation_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['original full model', 'error'] <- 
    Original_Full_Error_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['original full model', 'tau'] <- 
    Original_Full_Tau_Per_simulation[[paste0(treatment)]]
  results_treatments[[paste0('treatment_effect = ',treatment)]]['original full model', 'sd'] <- 
    Original_Full_Sd_Per_simulation[[paste0(treatment)]]
  
}

results_treatments_sd <- vector(mode = 'list')
blank_sd <- as.data.frame(matrix(NA, nrow = 8, ncol = 4))
names(blank_sd) <- c('estimation', 'sd', 'tau', 'error')
row.names(blank_sd) <- c('change score model',
                         'final score model',
                         'Trowman Method',
                         'pseudo base model',
                         'pseudo full model',
                         'Recovered ANCOVA',
                         'pseudo two-stage model',
                         'original full model')
for (treatment in Assess_Range){
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]] <- blank_sd
}

for (treatment in Assess_Range){
  
  ##Add value to the final table:
  Change_Model_Error_Per_simulation_sd <- apply(Change_Score_modelaggregate_error, 2, sd)
  Change_Model_Estimation_Per_simulation_sd <- apply(Change_Score_modelaggregate_estimation, 2, sd)
  Change_Model_Tau_Per_simulation_sd <- apply(Change_Score_modelaggregate_tau, 2, sd)
  Change_Model_Sd_Per_simulation_sd <- apply(Change_Score_modelaggregate_sd, 2, sd)
  
  Final_Model_Error_Per_simulation_sd <- apply(Final_Score_modelaggregate_error, 2, sd)
  Final_Model_Estimation_Per_simulation_sd <- apply(Final_Score_modelaggregate_estimation, 2, sd)
  Final_Model_Tau_Per_simulation_sd <- apply(Final_Score_modelaggregate_tau, 2, sd)
  Final_Model_Sd_Per_simulation_sd <- apply(Final_Score_modelaggregate_sd, 2, sd)
  
  Recover_ANCOVA_Error_Per_simulation_sd <- apply(Recover_ANCOVA_modelaggregate_error, 2, sd)
  Recover_ANCOVA_Estimation_Per_simulation_sd <- apply(Recover_ANCOVA_modelaggregate_estimation, 2, sd)
  Recover_ANCOVA_Tau_Per_simulation_sd <- apply(Recover_ANCOVA_modelaggregate_tau, 2, sd)
  Recover_ANCOVA_Sd_Per_simulation_sd <- apply(Recover_ANCOVA_modelaggregate_sd, 2, sd)
  
  Trowman_Error_Per_simulation_sd <- apply(Trowman_modelaggregate_error, 2, sd)
  Trowman_Estimation_Per_simulation_sd <- apply(Trowman_modelaggregate_estimation, 2, sd)
  Trowman_Tau_Per_simulation_sd <- apply(Trowman_modelaggregate_tau, 2, sd)
  Trowman_Sd_Per_simulation_sd <- apply(Trowman_modelaggregate_sd, 2, sd)
  
  Pseudo_Base_Error_Per_simulation_sd <- apply(Pseudo_Base_modelaggregate_error, 2, sd)
  Pseudo_Base_Estimation_Per_simulation_sd <- apply(Pseudo_Base_modelaggregate_estimation, 2, sd)
  Pseudo_Base_Tau_Per_simulation_sd <- apply(Pseudo_Base_modelaggregate_tau, 2, sd)
  Pseudo_Base_Sd_Per_simulation_sd <- apply(Pseudo_Base_modelaggregate_sd, 2, sd)
  
  Pseudo_Full_Error_Per_simulation_sd <- apply(Pseudo_Full_modelaggregate_error, 2, sd)
  Pseudo_Full_Estimation_Per_simulation_sd <- apply(Pseudo_Full_modelaggregate_estimation, 2, sd)
  Pseudo_Full_Tau_Per_simulation_sd <- apply(Pseudo_Full_modelaggregate_tau, 2, sd)
  Pseudo_Full_Sd_Per_simulation_sd <- apply(Pseudo_Full_modelaggregate_sd, 2, sd)
  
  Pseudo_Two_Stage_Error_Per_simulation_sd <- apply(Pseudo_Two_Stage_modelaggregate_error, 2, sd)
  Pseudo_Two_Stage_Estimation_Per_simulation_sd <- apply(Pseudo_Two_Stage_modelaggregate_estimation, 2, sd)
  Pseudo_Two_Stage_Tau_Per_simulation_sd <- apply(Pseudo_Two_Stage_modelaggregate_tau, 2, sd)
  Pseudo_Two_Stage_Sd_Per_simulation_sd <- apply(Pseudo_Two_Stage_modelaggregate_sd, 2, sd)
  
  Original_Full_Error_Per_simulation_sd <- apply(Original_Full_modelaggregate_error, 2, sd)
  Original_Full_Estimation_Per_simulation_sd <- apply(Original_Full_modelaggregate_estimation, 2, sd)
  Original_Full_Tau_Per_simulation_sd <- apply(Original_Full_modelaggregate_tau, 2, sd)
  Original_Full_Sd_Per_simulation_sd <- apply(Original_Full_modelaggregate_sd, 2, sd)
  
  
  
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['change score model', 'estimation'] <- 
    Change_Model_Estimation_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['change score model', 'error'] <- 
    Change_Model_Error_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['change score model', 'tau'] <- 
    Change_Model_Tau_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['change score model', 'sd'] <- 
    Change_Model_Sd_Per_simulation_sd[[paste0(treatment)]]
  
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['final score model', 'estimation'] <- 
    Final_Model_Estimation_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['final score model', 'error'] <- 
    Final_Model_Error_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['final score model', 'tau'] <- 
    Final_Model_Tau_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['final score model', 'sd'] <- 
    Final_Model_Sd_Per_simulation_sd[[paste0(treatment)]]
  
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['Recovered ANCOVA', 'estimation'] <- 
    Recover_ANCOVA_Estimation_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['Recovered ANCOVA', 'error'] <- 
    Recover_ANCOVA_Error_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['Recovered ANCOVA', 'tau'] <- 
    Recover_ANCOVA_Tau_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['Recovered ANCOVA', 'sd'] <- 
    Recover_ANCOVA_Sd_Per_simulation_sd[[paste0(treatment)]]
  
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['Trowman Method', 'estimation'] <- 
    Trowman_Estimation_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['Trowman Method', 'error'] <- 
    Trowman_Error_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['Trowman Method', 'tau'] <- 
    Trowman_Tau_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['Trowman Method', 'sd'] <- 
    Trowman_Sd_Per_simulation_sd[[paste0(treatment)]]
  
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['pseudo base model', 'estimation'] <- 
    Pseudo_Base_Estimation_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['pseudo base model', 'error'] <- 
    Pseudo_Base_Error_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['pseudo base model', 'tau'] <- 
    Pseudo_Base_Tau_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['pseudo base model', 'sd'] <- 
    Pseudo_Base_Sd_Per_simulation_sd[[paste0(treatment)]]
  
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['pseudo full model', 'estimation'] <- 
    Pseudo_Full_Estimation_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['pseudo full model', 'error'] <- 
    Pseudo_Full_Error_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['pseudo full model', 'tau'] <- 
    Pseudo_Full_Tau_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['pseudo full model', 'sd'] <- 
    Pseudo_Full_Sd_Per_simulation_sd[[paste0(treatment)]]
  
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['pseudo two-stage model', 'estimation'] <- 
    Pseudo_Two_Stage_Estimation_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['pseudo two-stage model', 'error'] <- 
    Pseudo_Two_Stage_Error_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['pseudo two-stage model', 'tau'] <- 
    Pseudo_Two_Stage_Tau_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['pseudo two-stage model', 'sd'] <- 
    Pseudo_Two_Stage_Sd_Per_simulation_sd[[paste0(treatment)]]
  
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['original full model', 'estimation'] <- 
    Original_Full_Estimation_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['original full model', 'error'] <- 
    Original_Full_Error_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['original full model', 'tau'] <- 
    Original_Full_Tau_Per_simulation_sd[[paste0(treatment)]]
  results_treatments_sd[[paste0('treatment_effect = ',treatment)]]['original full model', 'sd'] <- 
    Original_Full_Sd_Per_simulation_sd[[paste0(treatment)]]
  
}


#----------------------------------------------------------------------------------------------------------------------------
#                        The summary results 
#----------------------------------------------------------------------------------------------------------------------------


results_treatments_final <- vector(mode = 'list')
for (i in 1:length(results_treatments)){
  results_treatments_final[[paste0('treatment_effect = ', Assess_Range[i])]] <- 
    cbind(results_treatments[[i]], results_treatments_sd[[i]])
  names(results_treatments_final[[paste0('treatment_effect = ', Assess_Range[i])]]) <- 
    c('estimation', 'sd', 'tau', 'error', 'estimation_sd', 'sd_sd', 'tau_sd', 'error_sd')
}

















#----------------------------------------------------------------------------------------------------------------------------
#                       Prepare the blank data.frames to storage the simulation results for 
#                       the experiment with different baseline effect
#----------------------------------------------------------------------------------------------------------------------------






set.seed(123213123)

#The candidate treatment effects.

Assess_Range_Baseline <- c(0.2,0.5,0.8)
plot_number_Baseline <- 1:length(Assess_Range_Baseline)
Repeat_times_Baseline <- 200

baseline_effect_results <- vector(mode = 'list')


##############################################The data.frame to store the error.

Baseline_Change_Score_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Change_Score_modelaggregate_error) <- Assess_Range_Baseline

Baseline_Final_Score_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Final_Score_modelaggregate_error) <- Assess_Range_Baseline

Baseline_Recover_ANCOVA_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Recover_ANCOVA_modelaggregate_error) <- Assess_Range_Baseline

Baseline_Trowman_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Trowman_modelaggregate_error) <- Assess_Range_Baseline

Baseline_Pseudo_Base_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Pseudo_Base_modelaggregate_error) <- Assess_Range_Baseline

Baseline_Pseudo_Full_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Pseudo_Full_modelaggregate_error) <- Assess_Range_Baseline

Baseline_Pseudo_Two_Stage_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Pseudo_Two_Stage_modelaggregate_error) <- Assess_Range_Baseline


#########################################The data.frame to store the estimation


Baseline_Change_Score_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Change_Score_modelaggregate_estimation) <- Assess_Range_Baseline

Baseline_Final_Score_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Final_Score_modelaggregate_estimation) <- Assess_Range_Baseline

Baseline_Recover_ANCOVA_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Recover_ANCOVA_modelaggregate_estimation) <- Assess_Range_Baseline

Baseline_Trowman_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Trowman_modelaggregate_estimation) <- Assess_Range_Baseline

Baseline_Pseudo_Base_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Pseudo_Base_modelaggregate_estimation) <- Assess_Range_Baseline

Baseline_Pseudo_Full_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Pseudo_Full_modelaggregate_estimation) <- Assess_Range_Baseline

Baseline_Pseudo_Two_Stage_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Pseudo_Two_Stage_modelaggregate_estimation) <- Assess_Range_Baseline
#########################################The tau


Baseline_Change_Score_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Change_Score_modelaggregate_tau) <- Assess_Range_Baseline

Baseline_Final_Score_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Final_Score_modelaggregate_tau) <- Assess_Range_Baseline

Baseline_Recover_ANCOVA_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Recover_ANCOVA_modelaggregate_tau) <- Assess_Range_Baseline

Baseline_Trowman_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Trowman_modelaggregate_tau) <- Assess_Range_Baseline

Baseline_Pseudo_Base_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Pseudo_Base_modelaggregate_tau) <- Assess_Range_Baseline

Baseline_Pseudo_Full_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Pseudo_Full_modelaggregate_tau) <- Assess_Range_Baseline

Baseline_Pseudo_Two_Stage_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Pseudo_Two_Stage_modelaggregate_tau) <- Assess_Range_Baseline
#######################################################The Sd

Baseline_Change_Score_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Change_Score_modelaggregate_sd) <- Assess_Range_Baseline

Baseline_Final_Score_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Final_Score_modelaggregate_sd) <- Assess_Range_Baseline

Baseline_Recover_ANCOVA_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Recover_ANCOVA_modelaggregate_sd) <- Assess_Range_Baseline

Baseline_Trowman_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Trowman_modelaggregate_sd) <- Assess_Range_Baseline

Baseline_Pseudo_Base_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Pseudo_Base_modelaggregate_sd) <- Assess_Range_Baseline

Baseline_Pseudo_Full_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Pseudo_Full_modelaggregate_sd) <- Assess_Range_Baseline

Baseline_Pseudo_Two_Stage_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times_Baseline, ncol = length(Assess_Range_Baseline)))
names(Baseline_Pseudo_Two_Stage_modelaggregate_sd) <- Assess_Range_Baseline



#----------------------------------------------------------------------------------------------------------------------------
#                       The different baseline effect experiment. We can set the test baseline by 
#                       changing the 'Assess_range_Baseline'
#----------------------------------------------------------------------------------------------------------------------------



for (item in 1:Repeat_times_Baseline){
  for (baseline_effect in Assess_Range_Baseline){
    
    data.meta1 <- Simulation(treatment_effect = -50, 
                             baseline_effect = baseline_effect, tau = 13, sigma_ik = residual_exp,
                             sigma_YB = replicate(nstudy_exp,sd_baseline_exp),beta0  = replicate(nstudy_exp, 40),
                             nstudy = nstudy_exp, ngroup = as.data.frame(matrix(data = replicate(2*nstudy_exp, ngroup_exp), 
                                                                       nrow = nstudy_exp, ncol = 2)))#, imbalance = replicate(nstudy_exp, imbalance_exp))
    #  Simulation(treatment_effect = treatment_effect, 
    #                         baseline_effect = baseline_effect, tau = 4, sigma_ik = rep(1,8),
    #                         sigma_YB = sigma_YB, mean_baseline = mean_baseline,correlation = correlation, 
    #                         nstudy = nstudy, ngroup = ngroup)
    

    
    aggregate_meta1 <- cbind(aggregate(YB~study+treatment, data=data.meta1$simulated_data, mean), 
                             aggregate(YB~study+treatment, data=data.meta1$simulated_data, sd)[3],
                             aggregate(YF~study+treatment, data=data.meta1$simulated_data, mean)[3],
                             
                             aggregate(YF~study+treatment, data=data.meta1$simulated_data, sd)[3], 
                             aggregate(YChange~study+treatment, data = data.meta1$simulated_data, mean)[3],
                             aggregate(YChange~study+treatment, data = data.meta1$simulated_data, sd)[3],
                             c(data.meta1$ngroup[,2], data.meta1$ngroup[,1]))
    aggregate_meta1$group <- aggregate_meta1$treatment
    ID <- aggregate_meta1$study
    aggregate_meta1 <- cbind(ID, aggregate_meta1)
    #aggregate_meta1$study <- as.factor(aggregate_meta1$study)
    
    Correlation <- data.meta1$simulated_data %>%
      group_by(treatment, study) %>%
      summarise(Correlation = cor(YF, YB))
    
    aggregate_meta1$Correlation <- Correlation$Correlation
    
    aggregate_meta1$treatment <- NULL
    names(aggregate_meta1) <- c("ID", 'Study', 
                                'MeanBaseline','sdBaseline',
                                'MeanPostBaseline','sdPostBaseline', 
                                'MeanChange', 'sdChange', 
                                'NCFB', 'group', 'Correlation')
    data <- aggregate_meta1
    results <- simulation_process(data, data.meta1)
    
    baseline_effect_results[[paste0('baseline effect = ',baseline_effect)]] <- results$results
    
    #baseline_effect_results[[paste0('baseline_effect = ',baseline_effect)]] <- results$results
    
    #mean_baseline_results[[paste0('mean_baseline = ',mean_baseline)]] <- results$results
    
    #correlation_various_results[[paste0('correlation = ',correlation)]] <- results$results
  }
  for (len in 1:length(baseline_effect_results)){
    Baseline_Change_Score_modelaggregate_error[item,len] <- baseline_effect_results[[len]]['change score model','abs error with the true value']
    Baseline_Final_Score_modelaggregate_error[item,len] <- baseline_effect_results[[len]]['final score model','abs error with the true value']
    Baseline_Recover_ANCOVA_modelaggregate_error[item,len] <- baseline_effect_results[[len]]['Recovered ANCOVA','abs error with the true value']
    Baseline_Trowman_modelaggregate_error[item,len] <- baseline_effect_results[[len]]['Trowman Method', 'abs error with the true value']
    Baseline_Pseudo_Base_modelaggregate_error[item,len] <- baseline_effect_results[[len]]['pseudo base model', 'abs error with the true value']
    Baseline_Pseudo_Full_modelaggregate_error[item,len] <- baseline_effect_results[[len]]['pseudo full model', 'abs error with the true value']
    Baseline_Pseudo_Two_Stage_modelaggregate_error[item,len] <- baseline_effect_results[[len]]['pseudo two-stage model', 'abs error with the true value']
    
    Baseline_Change_Score_modelaggregate_estimation[item,len] <- baseline_effect_results[[len]]['change score model','estimate']
    Baseline_Final_Score_modelaggregate_estimation[item,len] <- baseline_effect_results[[len]]['final score model','estimate']
    Baseline_Recover_ANCOVA_modelaggregate_estimation[item,len] <- baseline_effect_results[[len]]['Recovered ANCOVA','estimate']
    Baseline_Trowman_modelaggregate_estimation[item,len] <- baseline_effect_results[[len]]['Trowman Method', 'estimate']
    Baseline_Pseudo_Base_modelaggregate_estimation[item,len] <- baseline_effect_results[[len]]['pseudo base model', 'estimate']
    Baseline_Pseudo_Full_modelaggregate_estimation[item,len] <- baseline_effect_results[[len]]['pseudo full model', 'estimate']
    Baseline_Pseudo_Two_Stage_modelaggregate_estimation[item,len] <- baseline_effect_results[[len]]['pseudo two-stage model', 'estimate']
    
    Baseline_Change_Score_modelaggregate_tau[item,len] <- baseline_effect_results[[len]]['change score model','tau']
    Baseline_Final_Score_modelaggregate_tau[item,len] <- baseline_effect_results[[len]]['final score model','tau']
    Baseline_Recover_ANCOVA_modelaggregate_tau[item,len] <- baseline_effect_results[[len]]['Recovered ANCOVA','tau']
    Baseline_Trowman_modelaggregate_tau[item,len] <- baseline_effect_results[[len]]['Trowman Method', 'tau']
    Baseline_Pseudo_Base_modelaggregate_tau[item,len] <- baseline_effect_results[[len]]['pseudo base model', 'tau']
    Baseline_Pseudo_Full_modelaggregate_tau[item,len] <- baseline_effect_results[[len]]['pseudo full model', 'tau']
    Baseline_Pseudo_Two_Stage_modelaggregate_tau[item,len] <- baseline_effect_results[[len]]['pseudo two-stage model', 'tau']
    
    Baseline_Change_Score_modelaggregate_sd[item,len] <- baseline_effect_results[[len]]['change score model','se']
    Baseline_Final_Score_modelaggregate_sd[item,len] <- baseline_effect_results[[len]]['final score model','se']
    Baseline_Recover_ANCOVA_modelaggregate_sd[item,len] <- baseline_effect_results[[len]]['Recovered ANCOVA','se']
    Baseline_Trowman_modelaggregate_sd[item,len] <- baseline_effect_results[[len]]['Trowman Method', 'se']
    Baseline_Pseudo_Base_modelaggregate_sd[item,len] <- baseline_effect_results[[len]]['pseudo base model', 'se']
    Baseline_Pseudo_Full_modelaggregate_sd[item,len] <- baseline_effect_results[[len]]['pseudo full model', 'se']
    Baseline_Pseudo_Two_Stage_modelaggregate_sd[item,len] <- baseline_effect_results[[len]]['pseudo two-stage model', 'se']
  }
  
}

##Make a table with Mean, se, tau, ci, abs error with the true values. And the mean and sd of them.

for (plt in plot_number_Baseline){
  plot(x = 1:Repeat_times_Baseline, y = Baseline_Change_Score_modelaggregate_error[,paste0(Assess_Range_Baseline[plt])], type = 'b', pch = 20, col = 'red',xlab = 'index', ylab = 'errors')
  lines(x = 1:Repeat_times_Baseline, y = Baseline_Final_Score_modelaggregate_error[,paste0(Assess_Range_Baseline[plt])], pch = 20, col = 'blue', type = 'b')
  lines(x = 1:Repeat_times_Baseline, y = Baseline_Recover_ANCOVA_modelaggregate_error[,paste0(Assess_Range_Baseline[plt])], pch = 20, col = 'yellow', type = 'b')
  lines(x = 1:Repeat_times_Baseline, y = Baseline_Trowman_modelaggregate_error[,paste0(Assess_Range_Baseline[plt])], pch = 20, col = 'green', type = 'b')
  lines(x = 1:Repeat_times_Baseline, y = Baseline_Pseudo_Base_modelaggregate_error[,paste0(Assess_Range_Baseline[plt])], pch = 20, col = 'brown', type = 'b')
  lines(x = 1:Repeat_times_Baseline, y = Baseline_Pseudo_Full_modelaggregate_error[,paste0(Assess_Range_Baseline[plt])], pch = 20, col = 'orange', type = 'b')
  lines(x = 1:Repeat_times_Baseline, y = Baseline_Pseudo_Two_Stage_modelaggregate_error[,paste0(Assess_Range_Baseline[plt])], pch = 20, col = 'purple', type = 'b')
  
  legend('topright', legend = c('change model', 
                                'final model', 
                                'Recovered ANCOVA',
                                'Trowman Model', 
                                'pseudo base model', 
                                'pseudo full model',
                                'pseudo two-stage model'), col = c('red', 'blue',
                                                                   'yellow', 'green', 
                                                                   'brown', 'orange',
                                                                   'purple'), pch = 20)
}

results_baselines <- vector(mode = 'list')
blank_baseline <- as.data.frame(matrix(NA, nrow = 7, ncol = 4))
names(blank_baseline) <- c('estimation', 'sd', 'tau', 'error')
row.names(blank_baseline) <- c('change score model', 
                               'final score model',
                               'Trowman Method', 
                               'pseudo base model',
                               'pseudo full model',
                               'Recovered ANCOVA',
                               'pseudo two-stage model')
for (baseline in Assess_Range_Baseline){
  results_baselines[[paste0('baseline_effect = ',baseline)]] <- blank_baseline
}

for (baseline in Assess_Range_Baseline){
  
  ##Add value to the final table:
  Baseline_Change_Model_Error_Per_simulation <- apply(Baseline_Change_Score_modelaggregate_error, 2, mean)
  Baseline_Change_Model_Estimation_Per_simulation <- apply(Baseline_Change_Score_modelaggregate_estimation, 2, mean)
  Baseline_Change_Model_Tau_Per_simulation <- apply(Baseline_Change_Score_modelaggregate_tau, 2, mean)
  Baseline_Change_Model_Sd_Per_simulation <- apply(Baseline_Change_Score_modelaggregate_sd, 2, mean)
  
  Baseline_Final_Model_Error_Per_simulation <- apply(Baseline_Final_Score_modelaggregate_error, 2, mean)
  Baseline_Final_Model_Estimation_Per_simulation <- apply(Baseline_Final_Score_modelaggregate_estimation, 2, mean)
  Baseline_Final_Model_Tau_Per_simulation <- apply(Baseline_Final_Score_modelaggregate_tau, 2, mean)
  Baseline_Final_Model_Sd_Per_simulation <- apply(Baseline_Final_Score_modelaggregate_sd, 2, mean)
  
  Baseline_Recover_ANCOVA_Error_Per_simulation <- apply(Baseline_Recover_ANCOVA_modelaggregate_error, 2, mean)
  Baseline_Recover_ANCOVA_Estimation_Per_simulation <- apply(Baseline_Recover_ANCOVA_modelaggregate_estimation, 2, mean)
  Baseline_Recover_ANCOVA_Tau_Per_simulation <- apply(Baseline_Recover_ANCOVA_modelaggregate_tau, 2, mean)
  Baseline_Recover_ANCOVA_Sd_Per_simulation <- apply(Baseline_Recover_ANCOVA_modelaggregate_sd, 2, mean)
  
  Baseline_Trowman_Error_Per_simulation <- apply(Baseline_Trowman_modelaggregate_error, 2, mean)
  Baseline_Trowman_Estimation_Per_simulation <- apply(Baseline_Trowman_modelaggregate_estimation, 2, mean)
  Baseline_Trowman_Tau_Per_simulation <- apply(Baseline_Trowman_modelaggregate_tau, 2, mean)
  Baseline_Trowman_Sd_Per_simulation <- apply(Baseline_Trowman_modelaggregate_sd, 2, mean)
  
  Baseline_Pseudo_Base_Error_Per_simulation <- apply(Baseline_Pseudo_Base_modelaggregate_error, 2, mean)
  Baseline_Pseudo_Base_Estimation_Per_simulation <- apply(Baseline_Pseudo_Base_modelaggregate_estimation, 2, mean)
  Baseline_Pseudo_Base_Tau_Per_simulation <- apply(Baseline_Pseudo_Base_modelaggregate_tau, 2, mean)
  Baseline_Pseudo_Base_Sd_Per_simulation <- apply(Baseline_Pseudo_Base_modelaggregate_sd, 2, mean)
  
  Baseline_Pseudo_Full_Error_Per_simulation <- apply(Baseline_Pseudo_Full_modelaggregate_error, 2, mean)
  Baseline_Pseudo_Full_Estimation_Per_simulation <- apply(Baseline_Pseudo_Full_modelaggregate_estimation, 2, mean)
  Baseline_Pseudo_Full_Tau_Per_simulation <- apply(Baseline_Pseudo_Full_modelaggregate_tau, 2, mean)
  Baseline_Pseudo_Full_Sd_Per_simulation <- apply(Baseline_Pseudo_Full_modelaggregate_sd, 2, mean)
 
  Baseline_Pseudo_Two_Stage_Error_Per_simulation <- apply(Baseline_Pseudo_Two_Stage_modelaggregate_error, 2, mean)
  Baseline_Pseudo_Two_Stage_Estimation_Per_simulation <- apply(Baseline_Pseudo_Two_Stage_modelaggregate_estimation, 2, mean)
  Baseline_Pseudo_Two_Stage_Tau_Per_simulation <- apply(Baseline_Pseudo_Two_Stage_modelaggregate_tau, 2, mean)
  Baseline_Pseudo_Two_Stage_Sd_Per_simulation <- apply(Baseline_Pseudo_Two_Stage_modelaggregate_sd, 2, mean) 
  
  
  
  
  
  results_baselines[[paste0('baseline_effect = ',baseline)]]['change score model', 'estimation'] <- 
    Baseline_Change_Model_Estimation_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['change score model', 'error'] <- 
    Baseline_Change_Model_Error_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['change score model', 'tau'] <- 
    Baseline_Change_Model_Tau_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['change score model', 'sd'] <- 
    Baseline_Change_Model_Sd_Per_simulation[[paste0(baseline)]]
  
  results_baselines[[paste0('baseline_effect = ',baseline)]]['final score model', 'estimation'] <- 
    Baseline_Final_Model_Estimation_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['final score model', 'error'] <- 
    Baseline_Final_Model_Error_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['final score model', 'tau'] <- 
    Baseline_Final_Model_Tau_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['final score model', 'sd'] <- 
    Baseline_Final_Model_Sd_Per_simulation[[paste0(baseline)]]
  
  results_baselines[[paste0('baseline_effect = ',baseline)]]['Recovered ANCOVA', 'estimation'] <- 
    Baseline_Recover_ANCOVA_Estimation_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['Recovered ANCOVA', 'error'] <- 
    Baseline_Recover_ANCOVA_Error_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['Recovered ANCOVA', 'tau'] <- 
    Baseline_Recover_ANCOVA_Tau_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['Recovered ANCOVA', 'sd'] <- 
    Baseline_Recover_ANCOVA_Sd_Per_simulation[[paste0(baseline)]]
  
  results_baselines[[paste0('baseline_effect = ',baseline)]]['Trowman Method', 'estimation'] <- 
    Baseline_Trowman_Estimation_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['Trowman Method', 'error'] <- 
    Baseline_Trowman_Error_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['Trowman Method', 'tau'] <- 
    Baseline_Trowman_Tau_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['Trowman Method', 'sd'] <- 
    Baseline_Trowman_Sd_Per_simulation[[paste0(baseline)]]
  
  results_baselines[[paste0('baseline_effect = ',baseline)]]['pseudo base model', 'estimation'] <- 
    Baseline_Pseudo_Base_Estimation_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['pseudo base model', 'error'] <- 
    Baseline_Pseudo_Base_Error_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['pseudo base model', 'tau'] <- 
    Baseline_Pseudo_Base_Tau_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['pseudo base model', 'sd'] <- 
    Baseline_Pseudo_Base_Sd_Per_simulation[[paste0(baseline)]]
  
  results_baselines[[paste0('baseline_effect = ',baseline)]]['pseudo full model', 'estimation'] <- 
    Baseline_Pseudo_Full_Estimation_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['pseudo full model', 'error'] <- 
    Baseline_Pseudo_Full_Error_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['pseudo full model', 'tau'] <- 
    Baseline_Pseudo_Full_Tau_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['pseudo full model', 'sd'] <- 
    Baseline_Pseudo_Full_Sd_Per_simulation[[paste0(baseline)]]
  
  results_baselines[[paste0('baseline_effect = ',baseline)]]['pseudo two-stage model', 'estimation'] <- 
    Baseline_Pseudo_Two_Stage_Estimation_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['pseudo two-stage model', 'error'] <- 
    Baseline_Pseudo_Two_Stage_Error_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['pseudo two-stage model', 'tau'] <- 
    Baseline_Pseudo_Two_Stage_Tau_Per_simulation[[paste0(baseline)]]
  results_baselines[[paste0('baseline_effect = ',baseline)]]['pseudo two-stage model', 'sd'] <- 
    Baseline_Pseudo_Two_Stage_Sd_Per_simulation[[paste0(baseline)]]
  
}

results_baselines_sd <- vector(mode = 'list')
blank_sd <- as.data.frame(matrix(NA, nrow = 7, ncol = 4))
names(blank_sd) <- c('estimation', 'sd', 'tau', 'error')
row.names(blank_sd) <- c('change score model',
                         'final score model',
                         'Trowman Method',
                         'pseudo base model',
                         'pseudo full model', 
                         'Recovered ANCOVA',
                         'pseudo two-stage model')
for (baseline in Assess_Range_Baseline){
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]] <- blank_sd
}

for (baseline in Assess_Range_Baseline){
  
  ##Add value to the final table:
  Baseline_Change_Model_Error_Per_simulation_sd <- apply(Baseline_Change_Score_modelaggregate_error, 2, sd)
  Baseline_Change_Model_Estimation_Per_simulation_sd <- apply(Baseline_Change_Score_modelaggregate_estimation, 2, sd)
  Baseline_Change_Model_Tau_Per_simulation_sd <- apply(Baseline_Change_Score_modelaggregate_tau, 2, sd)
  Baseline_Change_Model_Sd_Per_simulation_sd <- apply(Baseline_Change_Score_modelaggregate_sd, 2, sd)
  
  Baseline_Final_Model_Error_Per_simulation_sd <- apply(Baseline_Final_Score_modelaggregate_error, 2, sd)
  Baseline_Final_Model_Estimation_Per_simulation_sd <- apply(Baseline_Final_Score_modelaggregate_estimation, 2, sd)
  Baseline_Final_Model_Tau_Per_simulation_sd <- apply(Baseline_Final_Score_modelaggregate_tau, 2, sd)
  Baseline_Final_Model_Sd_Per_simulation_sd <- apply(Baseline_Final_Score_modelaggregate_sd, 2, sd)
  
  Baseline_Recover_ANCOVA_Error_Per_simulation_sd <- apply(Baseline_Recover_ANCOVA_modelaggregate_error, 2, sd)
  Baseline_Recover_ANCOVA_Estimation_Per_simulation_sd <- apply(Baseline_Recover_ANCOVA_modelaggregate_estimation, 2, sd)
  Baseline_Recover_ANCOVA_Tau_Per_simulation_sd <- apply(Baseline_Recover_ANCOVA_modelaggregate_tau, 2, sd)
  Baseline_Recover_ANCOVA_Sd_Per_simulation_sd <- apply(Baseline_Recover_ANCOVA_modelaggregate_sd, 2, sd)
  
  Baseline_Trowman_Error_Per_simulation_sd <- apply(Baseline_Trowman_modelaggregate_error, 2, sd)
  Baseline_Trowman_Estimation_Per_simulation_sd <- apply(Baseline_Trowman_modelaggregate_estimation, 2, sd)
  Baseline_Trowman_Tau_Per_simulation_sd <- apply(Baseline_Trowman_modelaggregate_tau, 2, sd)
  Baseline_Trowman_Sd_Per_simulation_sd <- apply(Baseline_Trowman_modelaggregate_sd, 2, sd)
  
  Baseline_Pseudo_Base_Error_Per_simulation_sd <- apply(Baseline_Pseudo_Base_modelaggregate_error, 2, sd)
  Baseline_Pseudo_Base_Estimation_Per_simulation_sd <- apply(Baseline_Pseudo_Base_modelaggregate_estimation, 2, sd)
  Baseline_Pseudo_Base_Tau_Per_simulation_sd <- apply(Baseline_Pseudo_Base_modelaggregate_tau, 2, sd)
  Baseline_Pseudo_Base_Sd_Per_simulation_sd <- apply(Baseline_Pseudo_Base_modelaggregate_sd, 2, sd)
  
  Baseline_Pseudo_Full_Error_Per_simulation_sd <- apply(Baseline_Pseudo_Full_modelaggregate_error, 2, sd)
  Baseline_Pseudo_Full_Estimation_Per_simulation_sd <- apply(Baseline_Pseudo_Full_modelaggregate_estimation, 2, sd)
  Baseline_Pseudo_Full_Tau_Per_simulation_sd <- apply(Baseline_Pseudo_Full_modelaggregate_tau, 2, sd)
  Baseline_Pseudo_Full_Sd_Per_simulation_sd <- apply(Baseline_Pseudo_Full_modelaggregate_sd, 2, sd)
  
  Baseline_Pseudo_Two_Stage_Error_Per_simulation_sd <- apply(Baseline_Pseudo_Two_Stage_modelaggregate_error, 2, sd)
  Baseline_Pseudo_Two_Stage_Estimation_Per_simulation_sd <- apply(Baseline_Pseudo_Two_Stage_modelaggregate_estimation, 2, sd)
  Baseline_Pseudo_Two_Stage_Tau_Per_simulation_sd <- apply(Baseline_Pseudo_Two_Stage_modelaggregate_tau, 2, sd)
  Baseline_Pseudo_Two_Stage_Sd_Per_simulation_sd <- apply(Baseline_Pseudo_Two_Stage_modelaggregate_sd, 2, sd)
  
  
  
  
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['change score model', 'estimation'] <- 
    Baseline_Change_Model_Estimation_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['change score model', 'error'] <- 
    Baseline_Change_Model_Error_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['change score model', 'tau'] <- 
    Baseline_Change_Model_Tau_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['change score model', 'sd'] <- 
    Baseline_Change_Model_Sd_Per_simulation_sd[[paste0(baseline)]]
  
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['final score model', 'estimation'] <- 
    Baseline_Final_Model_Estimation_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['final score model', 'error'] <- 
    Baseline_Final_Model_Error_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['final score model', 'tau'] <- 
    Baseline_Final_Model_Tau_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['final score model', 'sd'] <- 
    Baseline_Final_Model_Sd_Per_simulation_sd[[paste0(baseline)]]
  
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['Recovered ANCOVA', 'estimation'] <- 
    Baseline_Recover_ANCOVA_Estimation_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['Recovered ANCOVA', 'error'] <- 
    Baseline_Recover_ANCOVA_Error_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['Recovered ANCOVA', 'tau'] <- 
    Baseline_Recover_ANCOVA_Tau_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['Recovered ANCOVA', 'sd'] <- 
    Baseline_Recover_ANCOVA_Sd_Per_simulation_sd[[paste0(baseline)]]
  
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['Trowman Method', 'estimation'] <- 
    Baseline_Trowman_Estimation_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['Trowman Method', 'error'] <- 
    Baseline_Trowman_Error_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['Trowman Method', 'tau'] <- 
    Baseline_Trowman_Tau_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['Trowman Method', 'sd'] <- 
    Baseline_Trowman_Sd_Per_simulation_sd[[paste0(baseline)]]
  
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['pseudo base model', 'estimation'] <- 
    Baseline_Pseudo_Base_Estimation_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['pseudo base model', 'error'] <- 
    Baseline_Pseudo_Base_Error_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['pseudo base model', 'tau'] <- 
    Baseline_Pseudo_Base_Tau_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['pseudo base model', 'sd'] <- 
    Baseline_Pseudo_Base_Sd_Per_simulation_sd[[paste0(baseline)]]
  
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['pseudo full model', 'estimation'] <- 
    Baseline_Pseudo_Full_Estimation_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['pseudo full model', 'error'] <- 
    Baseline_Pseudo_Full_Error_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['pseudo full model', 'tau'] <- 
    Baseline_Pseudo_Full_Tau_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['pseudo full model', 'sd'] <- 
    Baseline_Pseudo_Full_Sd_Per_simulation_sd[[paste0(baseline)]]
  
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['pseudo two-stage model', 'estimation'] <- 
    Baseline_Pseudo_Two_Stage_Estimation_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['pseudo two-stage model', 'error'] <- 
    Baseline_Pseudo_Two_Stage_Error_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['pseudo two-stage model', 'tau'] <- 
    Baseline_Pseudo_Two_Stage_Tau_Per_simulation_sd[[paste0(baseline)]]
  results_baselines_sd[[paste0('baseline_effect = ',baseline)]]['pseudo two-stage model', 'sd'] <- 
    Baseline_Pseudo_Two_Stage_Sd_Per_simulation_sd[[paste0(baseline)]]
  
}
results_baselines_final <- vector(mode = 'list')
for (i in 1:length(results_baselines)){
  results_baselines_final[[paste0('baseline_effect = ', Assess_Range_Baseline[i])]] <- 
    cbind(results_baselines[[i]], results_baselines_sd[[i]])
  names(results_baselines_final[[paste0('baseline_effect = ', Assess_Range_Baseline[i])]]) <- 
    c('estimation', 'sd', 'tau', 'error', 'estimation_sd', 'sd_sd', 'tau_sd', 'error_sd')
}



#----------------------------------------------------------------------------------------------------------------------------
#                       Visualize the errors trend between models
#                       Visualize the estimates distribution with boxplot
#----------------------------------------------------------------------------------------------------------------------------





treatment_mse_table <- data.frame(matrix(nrow = 7, ncol = 1),
                                  row.names = c('change score model', 
                                                    'final score model', 
                                                    'Trowman Method', 
                                                    'pseudo base model', 
                                                    'pseudo full model',
                                                    'Recovered ANCOVA',
                                                    'pseudo two-stage model'))
names(treatment_mse_table) <- 'error'

Change_Score_model_error_sum <- 0
Final_Score_model_error_sum <- 0
Trowman_model_error_sum <- 0
Pseudo_Base_model_error_sum <- 0
Pseudo_Full_model_error_sum <- 0
Recovered_ANCOVA_model_error_sum <- 0
Pseudo_Two_Stage_model_error_sum <- 0

for (i in length(results_treatments_final)){
  Change_Score_model_error_sum <- Change_Score_model_error_sum + results_treatments_final[[i]]['change score model', 'error']
  Final_Score_model_error_sum <- Final_Score_model_error_sum + results_treatments_final[[i]]['final score model', 'error']
  Trowman_model_error_sum <- Trowman_model_error_sum + results_treatments_final[[i]]['Trowman Method', 'error']
  Pseudo_Base_model_error_sum <- Pseudo_Base_model_error_sum + results_treatments_final[[i]]['pseudo base model', 'error']
  Pseudo_Full_model_error_sum <- Pseudo_Full_model_error_sum + results_treatments_final[[i]]['pseudo full model', 'error']
  Recovered_ANCOVA_model_error_sum <- Recovered_ANCOVA_model_error_sum + results_treatments_final[[i]]['Recovered ANCOVA', 'error']
  Pseudo_Two_Stage_model_error_sum <- Pseudo_Two_Stage_model_error_sum + results_treatments_final[[i]]['pseudo two-stage model', 'error']
}
treatment_mse_table['change score model', 'error'] <- Change_Score_model_error_sum/7
treatment_mse_table['final score model', 'error'] <- Final_Score_model_error_sum/7
treatment_mse_table['Trowman Method', 'error'] <- Trowman_model_error_sum/7
treatment_mse_table['pseudo base model', 'error'] <- Pseudo_Base_model_error_sum/7
treatment_mse_table['pseudo full model', 'error'] <- Pseudo_Full_model_error_sum/7
treatment_mse_table['Recovered ANCOVA', 'error'] <- Recovered_ANCOVA_model_error_sum/7
treatment_mse_table['pseudo full model', 'error'] <- Pseudo_Two_Stage_model_error_sum /7





####################################################################################
#MSE of different methods for various baseline effect.




baseline_mse_table <- data.frame(matrix(nrow = 7, ncol = 1),
                                  row.names = c('change score model', 
                                                'final score model', 
                                                'Trowman Method', 
                                                'pseudo base model', 
                                                'pseudo full model',
                                                'Recovered ANCOVA',
                                                'pseudo two-stage model'))
names(baseline_mse_table) <- 'error'
Change_Score_model_error_sum <- 0
Final_Score_model_error_sum <- 0
Trowman_model_error_sum <- 0
Pseudo_Base_model_error_sum <- 0
Pseudo_Full_model_error_sum <- 0
Recovered_ANCOVA_model_error_sum <- 0
Pseudo_Two_Stage_model_error_sum <- 0

for (i in length(results_baselines_final)){
  Change_Score_model_error_sum <- Change_Score_model_error_sum + results_baselines_final[[i]]['change score model', 'error']
  Final_Score_model_error_sum <- Final_Score_model_error_sum + results_baselines_final[[i]]['final score model', 'error']
  Trowman_model_error_sum <- Trowman_model_error_sum + results_baselines_final[[i]]['Trowman Method', 'error']
  Pseudo_Base_model_error_sum <- Pseudo_Base_model_error_sum + results_baselines_final[[i]]['pseudo base model', 'error']
  Pseudo_Full_model_error_sum <- Pseudo_Full_model_error_sum + results_baselines_final[[i]]['pseudo full model', 'error']
  Recovered_ANCOVA_model_error_sum <- Recovered_ANCOVA_model_error_sum + results_baselines_final[[i]]['Recovered ANCOVA', 'error']
  Pseudo_Two_Stage_model_error_sum <- Pseudo_Two_Stage_model_error_sum + results_baselines_final[[i]]['pseudo two-stage model', 'error']
}
baseline_mse_table['change score model', 'error'] <- Change_Score_model_error_sum/7
baseline_mse_table['final score model', 'error'] <- Final_Score_model_error_sum/7
baseline_mse_table['Trowman Method', 'error'] <- Trowman_model_error_sum/7
baseline_mse_table['pseudo base model', 'error'] <- Pseudo_Base_model_error_sum/7
baseline_mse_table['pseudo full model', 'error'] <- Pseudo_Full_model_error_sum/7
baseline_mse_table['Recovered ANCOVA', 'error'] <- Recovered_ANCOVA_model_error_sum/7
baseline_mse_table['pseudo two-stage model', 'error'] <- Pseudo_Two_Stage_model_error_sum/7


#par(mfrow = c(2,2))
Baseline_estimation_results <- vector(mode = 'list')
for (i in 1:length(Assess_Range_Baseline)){
  Baseline_estimation_results[[paste0('baseline_effect = ', Assess_Range_Baseline[i])]] <- 
    cbind(Baseline_Change_Score_modelaggregate_estimation[,i],
          Baseline_Final_Score_modelaggregate_estimation[,i],
          Baseline_Recover_ANCOVA_modelaggregate_estimation[,i],
          Baseline_Trowman_modelaggregate_estimation[,i],
          Baseline_Pseudo_Base_modelaggregate_estimation[,i],
          Baseline_Pseudo_Full_modelaggregate_estimation[,i],
          Baseline_Pseudo_Two_Stage_modelaggregate_estimation[,i])
  boxplot.matrix(as.matrix(Baseline_estimation_results[[paste0('baseline_effect = ', Assess_Range_Baseline[i])]]),
                 xaxt = 'n',
                 xlab = 'Models',
                 ylab = 'estimation distribution',
                 main = paste0('distribution of the estimation for baseline effect = ', 
                               Assess_Range_Baseline[i]),
                 col = c('red','blue','yellow','green','brown','orange', 'purple'))
  axis(1, at = 1:7, labels = c('change', 'final', 'RA', 'Trowman', 'PB', 'PF', 'PTS'))
}
par(mfrow = c(1,1))














#----------------------------------------------------------------------------------------------------------------------------
#                       Prepare the blank data.frames to storage the simulation results for 
#                       the experiment with different random effect
#----------------------------------------------------------------------------------------------------------------------------



correlation_various_results <- vector(mode = 'list')
mean_baseline_results <- vector(mode = 'list')
baseline_effect_results <- vector(mode = 'list')
tau_effect_results <- vector(mode = 'list')




#Replicate the process and curve the results:

set.seed(3456)

#The candidate tau effects.

Assess_Range_Tau <- seq(6, 20, 7)
plot_number <- 1:length(Assess_Range_Tau)
Repeat_times <- 200

baseline_effect_results <- vector(mode = 'list')


##############################################The data.frame to store the error.

Tau_Change_Score_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Change_Score_modelaggregate_error) <- Assess_Range_Tau

Tau_Final_Score_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Final_Score_modelaggregate_error) <- Assess_Range_Tau

Tau_Recover_ANCOVA_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Recover_ANCOVA_modelaggregate_error) <- Assess_Range_Tau

Tau_Trowman_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Trowman_modelaggregate_error) <- Assess_Range_Tau

Tau_Pseudo_Base_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Pseudo_Base_modelaggregate_error) <- Assess_Range_Tau

Tau_Pseudo_Full_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Pseudo_Full_modelaggregate_error) <- Assess_Range_Tau

Tau_Pseudo_Two_Stage_modelaggregate_error <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Pseudo_Two_Stage_modelaggregate_error) <- Assess_Range_Tau
#########################################The data.frame to store the estimation


Tau_Change_Score_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Change_Score_modelaggregate_estimation) <- Assess_Range_Tau

Tau_Final_Score_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Final_Score_modelaggregate_estimation) <- Assess_Range_Tau

Tau_Recover_ANCOVA_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Recover_ANCOVA_modelaggregate_estimation) <- Assess_Range_Tau

Tau_Trowman_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Trowman_modelaggregate_estimation) <- Assess_Range_Tau

Tau_Pseudo_Base_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Pseudo_Base_modelaggregate_estimation) <- Assess_Range_Tau

Tau_Pseudo_Full_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Pseudo_Full_modelaggregate_estimation) <- Assess_Range_Tau

Tau_Pseudo_Two_Stage_modelaggregate_estimation <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Pseudo_Two_Stage_modelaggregate_estimation) <- Assess_Range_Tau

#########################################The tau


Tau_Change_Score_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Change_Score_modelaggregate_tau) <- Assess_Range_Tau

Tau_Final_Score_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Final_Score_modelaggregate_tau) <- Assess_Range_Tau

Tau_Recover_ANCOVA_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Recover_ANCOVA_modelaggregate_tau) <- Assess_Range_Tau

Tau_Trowman_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Trowman_modelaggregate_tau) <- Assess_Range_Tau

Tau_Pseudo_Base_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Pseudo_Base_modelaggregate_tau) <- Assess_Range_Tau

Tau_Pseudo_Full_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Pseudo_Full_modelaggregate_tau) <- Assess_Range_Tau

Tau_Pseudo_Two_Stage_modelaggregate_tau <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Pseudo_Two_Stage_modelaggregate_tau) <- Assess_Range_Tau

#######################################################The Sd

Tau_Change_Score_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Change_Score_modelaggregate_sd) <- Assess_Range_Tau

Tau_Final_Score_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Final_Score_modelaggregate_sd) <- Assess_Range_Tau

Tau_Recover_ANCOVA_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Recover_ANCOVA_modelaggregate_sd) <- Assess_Range_Tau

Tau_Trowman_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Trowman_modelaggregate_sd) <- Assess_Range_Tau

Tau_Pseudo_Base_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Pseudo_Base_modelaggregate_sd) <- Assess_Range_Tau

Tau_Pseudo_Full_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Pseudo_Full_modelaggregate_sd) <- Assess_Range_Tau

Tau_Pseudo_Two_Stage_modelaggregate_sd <- as.data.frame(matrix(NA, nrow = Repeat_times, ncol = length(Assess_Range_Tau)))
names(Tau_Pseudo_Two_Stage_modelaggregate_sd) <- Assess_Range_Tau



#----------------------------------------------------------------------------------------------------------------------------
#                       The different random effect experiment: 
#                       Change 'Assess_Range_Tau' to set the tested random effect 
#----------------------------------------------------------------------------------------------------------------------------



for (item in 1:Repeat_times){
  for (tau_effect in Assess_Range_Tau){
    
    data.meta1 <- Simulation(treatment_effect = -50, 
                             baseline_effect = 0.5, tau = tau_effect, sigma_ik = residual_exp,
                             sigma_YB = replicate(nstudy,sd_baseline_exp), beta0 = replicate(nstudy, 40),
                             nstudy = nstudy_exp, ngroup = as.data.frame(matrix(data = replicate(nstudy_exp, ngroup_exp), 
                                                                       nrow = nstudy, ncol = 2)), imbalance = replicate(nstudy_exp, imbalance_exp))
    #  Simulation(tau_effect = tau_effect, 
    #                         baseline_effect = baseline_effect, tau = 4, sigma_ik = rep(1,8),
    #                         sigma_YB = sigma_YB, mean_baseline = mean_baseline,correlation = correlation, 
    #                         nstudy = nstudy, ngroup = ngroup)
    
    
    aggregate_meta1 <- cbind(aggregate(YB~study+treatment, data=data.meta1$simulated_data, mean), 
                             aggregate(YB~study+treatment, data=data.meta1$simulated_data, sd)[3],
                             aggregate(YF~study+treatment, data=data.meta1$simulated_data, mean)[3],
                             
                             aggregate(YF~study+treatment, data=data.meta1$simulated_data, sd)[3], 
                             aggregate(YChange~study+treatment, data = data.meta1$simulated_data, mean)[3],
                             aggregate(YChange~study+treatment, data = data.meta1$simulated_data, sd)[3],
                             c(data.meta1$ngroup[,2], data.meta1$ngroup[,1]))
    aggregate_meta1$group <- aggregate_meta1$treatment
    ID <- aggregate_meta1$study
    aggregate_meta1 <- cbind(ID, aggregate_meta1)
   # aggregate_meta1$study <- as.factor(aggregate_meta1$study)
    
    Correlation <- data.meta1$simulated_data %>%
      group_by(treatment, study) %>%
      summarise(Correlation = cor(YF, YB))
    
    aggregate_meta1$Correlation <- Correlation$Correlation
    
    aggregate_meta1$treatment <- NULL
    names(aggregate_meta1) <- c("ID", 'Study', 
                                'MeanBaseline','sdBaseline',
                                'MeanPostBaseline','sdPostBaseline', 
                                'MeanChange', 'sdChange', 
                                'NCFB', 'group', 'Correlation')
    data <- aggregate_meta1
    results <- simulation_process(data, data.meta1)
    
    tau_effect_results[[paste0('tau_effect = ',tau_effect)]] <- results$results
    
    #baseline_effect_results[[paste0('baseline_effect = ',baseline_effect)]] <- results$results
    
    #mean_baseline_results[[paste0('mean_baseline = ',mean_baseline)]] <- results$results
    
    #correlation_various_results[[paste0('correlation = ',correlation)]] <- results$results
  }
  for (len in 1:length(tau_effect_results)){
    Tau_Change_Score_modelaggregate_error[item,len] <- tau_effect_results[[len]]['change score model','abs error with the true value']
    Tau_Final_Score_modelaggregate_error[item,len] <- tau_effect_results[[len]]['final score model','abs error with the true value']
    Tau_Recover_ANCOVA_modelaggregate_error[item,len] <- tau_effect_results[[len]]['Recovered ANCOVA','abs error with the true value']
    Tau_Trowman_modelaggregate_error[item,len] <- tau_effect_results[[len]]['Trowman Method', 'abs error with the true value']
    Tau_Pseudo_Base_modelaggregate_error[item,len] <- tau_effect_results[[len]]['pseudo base model', 'abs error with the true value']
    Tau_Pseudo_Full_modelaggregate_error[item,len] <- tau_effect_results[[len]]['pseudo full model', 'abs error with the true value']
    Tau_Pseudo_Two_Stage_modelaggregate_error[item,len] <- tau_effect_results[[len]]['pseudo two-stage model', 'abs error with the true value']
    
    
    Tau_Change_Score_modelaggregate_estimation[item,len] <- tau_effect_results[[len]]['change score model','estimate']
    Tau_Final_Score_modelaggregate_estimation[item,len] <- tau_effect_results[[len]]['final score model','estimate']
    Tau_Recover_ANCOVA_modelaggregate_estimation[item,len] <- tau_effect_results[[len]]['Recovered ANCOVA','estimate']
    Tau_Trowman_modelaggregate_estimation[item,len] <- tau_effect_results[[len]]['Trowman Method', 'estimate']
    Tau_Pseudo_Base_modelaggregate_estimation[item,len] <- tau_effect_results[[len]]['pseudo base model', 'estimate']
    Tau_Pseudo_Full_modelaggregate_estimation[item,len] <- tau_effect_results[[len]]['pseudo full model', 'estimate']
    Tau_Pseudo_Two_Stage_modelaggregate_estimation[item,len] <- tau_effect_results[[len]]['pseudo two-stage model', 'estimate']
    
    
    Tau_Change_Score_modelaggregate_tau[item,len] <- tau_effect_results[[len]]['change score model','tau']
    Tau_Final_Score_modelaggregate_tau[item,len] <- tau_effect_results[[len]]['final score model','tau']
    Tau_Recover_ANCOVA_modelaggregate_tau[item,len] <- tau_effect_results[[len]]['Recovered ANCOVA','tau']
    Tau_Trowman_modelaggregate_tau[item,len] <- tau_effect_results[[len]]['Trowman Method', 'tau']
    Tau_Pseudo_Base_modelaggregate_tau[item,len] <- tau_effect_results[[len]]['pseudo base model', 'tau']
    Tau_Pseudo_Full_modelaggregate_tau[item,len] <- tau_effect_results[[len]]['pseudo full model', 'tau']
    Tau_Pseudo_Two_Stage_modelaggregate_tau[item,len] <- tau_effect_results[[len]]['pseudo two-stage model', 'tau']
    
    
    Tau_Change_Score_modelaggregate_sd[item,len] <- tau_effect_results[[len]]['change score model','se']
    Tau_Final_Score_modelaggregate_sd[item,len] <- tau_effect_results[[len]]['final score model','se']
    Tau_Recover_ANCOVA_modelaggregate_sd[item,len] <- tau_effect_results[[len]]['Recovered ANCOVA','se']
    Tau_Trowman_modelaggregate_sd[item,len] <- tau_effect_results[[len]]['Trowman Method', 'se']
    Tau_Pseudo_Base_modelaggregate_sd[item,len] <- tau_effect_results[[len]]['pseudo base model', 'se']
    Tau_Pseudo_Full_modelaggregate_sd[item,len] <- tau_effect_results[[len]]['pseudo full model', 'se']
    Tau_Pseudo_Two_Stage_modelaggregate_sd[item,len] <- tau_effect_results[[len]]['pseudo two-stage model', 'se']
  }
  
}

##Make a table with Mean, se, tau, ci, abs error with the true values. And the mean and sd of them.

for (plt in plot_number){
  plot(x = 1:Repeat_times, y = Tau_Change_Score_modelaggregate_error[,paste0(Assess_Range_Tau[plt])], type = 'b', pch = 20, col = 'red', xlab = 'index', ylab = 'errors')
  lines(x = 1:Repeat_times, y = Tau_Final_Score_modelaggregate_error[,paste0(Assess_Range_Tau[plt])], pch = 20, col = 'blue', type = 'b')
  lines(x = 1:Repeat_times, y = Tau_Recover_ANCOVA_modelaggregate_error[,paste0(Assess_Range_Tau[plt])], pch = 20, col = 'yellow', type = 'b')
  lines(x = 1:Repeat_times, y = Tau_Trowman_modelaggregate_error[,paste0(Assess_Range_Tau[plt])], pch = 20, col = 'green', type = 'b')
  lines(x = 1:Repeat_times, y = Tau_Pseudo_Base_modelaggregate_error[,paste0(Assess_Range_Tau[plt])], pch = 20, col = 'brown', type = 'b')
  lines(x = 1:Repeat_times, y = Tau_Pseudo_Full_modelaggregate_error[,paste0(Assess_Range_Tau[plt])], pch = 20, col = 'orange', type = 'b')
  lines(x = 1:Repeat_times, y = Tau_Pseudo_Two_Stage_modelaggregate_error[,paste0(Assess_Range_Tau[plt])], pch = 20, col = 'purple', type = 'b')
  
  legend('topright', legend = c('change model',
                                'final model',
                                'Recovered ANCOVA',
                                'Trowman Model', 
                                'pseudo base model', 
                                'pseudo full model',
                                'pseudo two-stage model'),
         col = c('red', 'blue','yellow', 'green', 'brown','orange', 'purple'), pch = 20)
}

results_taus <- vector(mode = 'list')
blank <- as.data.frame(matrix(NA, nrow = 7, ncol = 4))
names(blank) <- c('estimation', 'sd', 'tau', 'error')
row.names(blank) <- c('change score model', 
                      'final score model',
                      'Trowman Method',
                      'pseudo base model',
                      'pseudo full model', 
                      'Recovered ANCOVA', 
                      'pseudo two-stage model')
for (tau in Assess_Range_Tau){
  results_taus[[paste0('tau_effect = ',tau)]] <- blank
}

for (tau in Assess_Range_Tau){
  
  ##Add value to the final table:
  Tau_Change_Model_Error_Per_simulation <- apply(Tau_Change_Score_modelaggregate_error, 2, mean)
  Tau_Change_Model_Estimation_Per_simulation <- apply(Tau_Change_Score_modelaggregate_estimation, 2, mean)
  Tau_Change_Model_Tau_Per_simulation <- apply(Tau_Change_Score_modelaggregate_tau, 2, mean)
  Tau_Change_Model_Sd_Per_simulation <- apply(Tau_Change_Score_modelaggregate_sd, 2, mean)
  
  Tau_Final_Model_Error_Per_simulation <- apply(Tau_Final_Score_modelaggregate_error, 2, mean)
  Tau_Final_Model_Estimation_Per_simulation <- apply(Tau_Final_Score_modelaggregate_estimation, 2, mean)
  Tau_Final_Model_Tau_Per_simulation <- apply(Tau_Final_Score_modelaggregate_tau, 2, mean)
  Tau_Final_Model_Sd_Per_simulation <- apply(Tau_Final_Score_modelaggregate_sd, 2, mean)
  
  Tau_Recover_ANCOVA_Error_Per_simulation <- apply(Tau_Recover_ANCOVA_modelaggregate_error, 2, mean)
  Tau_Recover_ANCOVA_Estimation_Per_simulation <- apply(Tau_Recover_ANCOVA_modelaggregate_estimation, 2, mean)
  Tau_Recover_ANCOVA_Tau_Per_simulation <- apply(Tau_Recover_ANCOVA_modelaggregate_tau, 2, mean)
  Tau_Recover_ANCOVA_Sd_Per_simulation <- apply(Tau_Recover_ANCOVA_modelaggregate_sd, 2, mean)
  
  Tau_Trowman_Error_Per_simulation <- apply(Tau_Trowman_modelaggregate_error, 2, mean)
  Tau_Trowman_Estimation_Per_simulation <- apply(Tau_Trowman_modelaggregate_estimation, 2, mean)
  Tau_Trowman_Tau_Per_simulation <- apply(Tau_Trowman_modelaggregate_tau, 2, mean)
  Tau_Trowman_Sd_Per_simulation <- apply(Tau_Trowman_modelaggregate_sd, 2, mean)
  
  Tau_Pseudo_Base_Error_Per_simulation <- apply(Tau_Pseudo_Base_modelaggregate_error, 2, mean)
  Tau_Pseudo_Base_Estimation_Per_simulation <- apply(Tau_Pseudo_Base_modelaggregate_estimation, 2, mean)
  Tau_Pseudo_Base_Tau_Per_simulation <- apply(Tau_Pseudo_Base_modelaggregate_tau, 2, mean)
  Tau_Pseudo_Base_Sd_Per_simulation <- apply(Tau_Pseudo_Base_modelaggregate_sd, 2, mean)
  
  Tau_Pseudo_Full_Error_Per_simulation <- apply(Tau_Pseudo_Full_modelaggregate_error, 2, mean)
  Tau_Pseudo_Full_Estimation_Per_simulation <- apply(Tau_Pseudo_Full_modelaggregate_estimation, 2, mean)
  Tau_Pseudo_Full_Tau_Per_simulation <- apply(Tau_Pseudo_Full_modelaggregate_tau, 2, mean)
  Tau_Pseudo_Full_Sd_Per_simulation <- apply(Tau_Pseudo_Full_modelaggregate_sd, 2, mean)
  
  Tau_Pseudo_Two_Stage_Error_Per_simulation <- apply(Tau_Pseudo_Two_Stage_modelaggregate_error, 2, mean)
  Tau_Pseudo_Two_Stage_Estimation_Per_simulation <- apply(Tau_Pseudo_Two_Stage_modelaggregate_estimation, 2, mean)
  Tau_Pseudo_Two_Stage_Tau_Per_simulation <- apply(Tau_Pseudo_Two_Stage_modelaggregate_tau, 2, mean)
  Tau_Pseudo_Two_Stage_Sd_Per_simulation <- apply(Tau_Pseudo_Two_Stage_modelaggregate_sd, 2, mean)
  
  
  
  
  results_taus[[paste0('tau_effect = ',tau)]]['change score model', 'estimation'] <- 
    Tau_Change_Model_Estimation_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['change score model', 'error'] <- 
    Tau_Change_Model_Error_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['change score model', 'tau'] <- 
    Tau_Change_Model_Tau_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['change score model', 'sd'] <- 
    Tau_Change_Model_Sd_Per_simulation[[paste0(tau)]]
  
  results_taus[[paste0('tau_effect = ',tau)]]['final score model', 'estimation'] <- 
    Tau_Final_Model_Estimation_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['final score model', 'error'] <- 
    Tau_Final_Model_Error_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['final score model', 'tau'] <- 
    Tau_Final_Model_Tau_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['final score model', 'sd'] <- 
    Tau_Final_Model_Sd_Per_simulation[[paste0(tau)]]
  
  results_taus[[paste0('tau_effect = ',tau)]]['Recovered ANCOVA', 'estimation'] <- 
    Tau_Recover_ANCOVA_Estimation_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['Recovered ANCOVA', 'error'] <- 
    Tau_Recover_ANCOVA_Error_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['Recovered ANCOVA', 'tau'] <- 
    Tau_Recover_ANCOVA_Tau_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['Recovered ANCOVA', 'sd'] <- 
    Tau_Recover_ANCOVA_Sd_Per_simulation[[paste0(tau)]]
  
  results_taus[[paste0('tau_effect = ',tau)]]['Trowman Method', 'estimation'] <- 
    Tau_Trowman_Estimation_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['Trowman Method', 'error'] <- 
    Tau_Trowman_Error_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['Trowman Method', 'tau'] <- 
    Tau_Trowman_Tau_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['Trowman Method', 'sd'] <- 
    Tau_Trowman_Sd_Per_simulation[[paste0(tau)]]
  
  results_taus[[paste0('tau_effect = ',tau)]]['pseudo base model', 'estimation'] <- 
    Tau_Pseudo_Base_Estimation_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['pseudo base model', 'error'] <- 
    Tau_Pseudo_Base_Error_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['pseudo base model', 'tau'] <- 
    Tau_Pseudo_Base_Tau_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['pseudo base model', 'sd'] <- 
    Tau_Pseudo_Base_Sd_Per_simulation[[paste0(tau)]]
  
  results_taus[[paste0('tau_effect = ',tau)]]['pseudo full model', 'estimation'] <- 
    Tau_Pseudo_Full_Estimation_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['pseudo full model', 'error'] <- 
    Tau_Pseudo_Full_Error_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['pseudo full model', 'tau'] <- 
    Tau_Pseudo_Full_Tau_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['pseudo full model', 'sd'] <- 
    Tau_Pseudo_Full_Sd_Per_simulation[[paste0(tau)]]
  
  results_taus[[paste0('tau_effect = ',tau)]]['pseudo two-stage model', 'estimation'] <- 
    Tau_Pseudo_Two_Stage_Estimation_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['pseudo two-stage model', 'error'] <- 
    Tau_Pseudo_Two_Stage_Error_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['pseudo two-stage model', 'tau'] <- 
    Tau_Pseudo_Two_Stage_Tau_Per_simulation[[paste0(tau)]]
  results_taus[[paste0('tau_effect = ',tau)]]['pseudo two-stage model', 'sd'] <- 
    Tau_Pseudo_Two_Stage_Sd_Per_simulation[[paste0(tau)]]
  
}

results_taus_sd <- vector(mode = 'list')
blank_sd <- as.data.frame(matrix(NA, nrow = 7, ncol = 4))
names(blank_sd) <- c('estimation', 'sd', 'tau', 'error')
row.names(blank_sd) <- c('change score model', 
                         'final score model',
                         'Trowman Method', 
                         'pseudo base model',
                         'pseudo full model', 
                         'Recovered ANCOVA',
                         'pseudo two-stage model')
for (tau in Assess_Range_Tau){
  results_taus_sd[[paste0('tau_effect = ',tau)]] <- blank_sd
}

for (tau in Assess_Range_Tau){
  
  ##Add value to the final table:
  Tau_Change_Model_Error_Per_simulation_sd <- apply(Tau_Change_Score_modelaggregate_error, 2, sd)
  Tau_Change_Model_Estimation_Per_simulation_sd <- apply(Tau_Change_Score_modelaggregate_estimation, 2, sd)
  Tau_Change_Model_Tau_Per_simulation_sd <- apply(Tau_Change_Score_modelaggregate_tau, 2, sd)
  Tau_Change_Model_Sd_Per_simulation_sd <- apply(Tau_Change_Score_modelaggregate_sd, 2, sd)
  
  Tau_Final_Model_Error_Per_simulation_sd <- apply(Tau_Final_Score_modelaggregate_error, 2, sd)
  Tau_Final_Model_Estimation_Per_simulation_sd <- apply(Tau_Final_Score_modelaggregate_estimation, 2, sd)
  Tau_Final_Model_Tau_Per_simulation_sd <- apply(Tau_Final_Score_modelaggregate_tau, 2, sd)
  Tau_Final_Model_Sd_Per_simulation_sd <- apply(Tau_Final_Score_modelaggregate_sd, 2, sd)
  
  Tau_Recover_ANCOVA_Error_Per_simulation_sd <- apply(Tau_Recover_ANCOVA_modelaggregate_error, 2, sd)
  Tau_Recover_ANCOVA_Estimation_Per_simulation_sd <- apply(Tau_Recover_ANCOVA_modelaggregate_estimation, 2, sd)
  Tau_Recover_ANCOVA_Tau_Per_simulation_sd <- apply(Tau_Recover_ANCOVA_modelaggregate_tau, 2, sd)
  Tau_Recover_ANCOVA_Sd_Per_simulation_sd <- apply(Tau_Recover_ANCOVA_modelaggregate_sd, 2, sd)
  
  Tau_Trowman_Error_Per_simulation_sd <- apply(Tau_Trowman_modelaggregate_error, 2, sd)
  Tau_Trowman_Estimation_Per_simulation_sd <- apply(Tau_Trowman_modelaggregate_estimation, 2, sd)
  Tau_Trowman_Tau_Per_simulation_sd <- apply(Tau_Trowman_modelaggregate_tau, 2, sd)
  Tau_Trowman_Sd_Per_simulation_sd <- apply(Tau_Trowman_modelaggregate_sd, 2, sd)
  
  Tau_Pseudo_Base_Error_Per_simulation_sd <- apply(Tau_Pseudo_Base_modelaggregate_error, 2, sd)
  Tau_Pseudo_Base_Estimation_Per_simulation_sd <- apply(Tau_Pseudo_Base_modelaggregate_estimation, 2, sd)
  Tau_Pseudo_Base_Tau_Per_simulation_sd <- apply(Tau_Pseudo_Base_modelaggregate_tau, 2, sd)
  Tau_Pseudo_Base_Sd_Per_simulation_sd <- apply(Tau_Pseudo_Base_modelaggregate_sd, 2, sd)
  
  Tau_Pseudo_Full_Error_Per_simulation_sd <- apply(Tau_Pseudo_Full_modelaggregate_error, 2, sd)
  Tau_Pseudo_Full_Estimation_Per_simulation_sd <- apply(Tau_Pseudo_Full_modelaggregate_estimation, 2, sd)
  Tau_Pseudo_Full_Tau_Per_simulation_sd <- apply(Tau_Pseudo_Full_modelaggregate_tau, 2, sd)
  Tau_Pseudo_Full_Sd_Per_simulation_sd <- apply(Tau_Pseudo_Full_modelaggregate_sd, 2, sd)
  
  Tau_Pseudo_Two_Stage_Error_Per_simulation_sd <- apply(Tau_Pseudo_Two_Stage_modelaggregate_error, 2, sd)
  Tau_Pseudo_Two_Stage_Estimation_Per_simulation_sd <- apply(Tau_Pseudo_Two_Stage_modelaggregate_estimation, 2, sd)
  Tau_Pseudo_Two_Stage_Tau_Per_simulation_sd <- apply(Tau_Pseudo_Two_Stage_modelaggregate_tau, 2, sd)
  Tau_Pseudo_Two_Stage_Sd_Per_simulation_sd <- apply(Tau_Pseudo_Two_Stage_modelaggregate_sd, 2, sd)
  
  
  
  
  results_taus_sd[[paste0('tau_effect = ',tau)]]['change score model', 'estimation'] <- 
    Tau_Change_Model_Estimation_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['change score model', 'error'] <- 
    Tau_Change_Model_Error_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['change score model', 'tau'] <- 
    Tau_Change_Model_Tau_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['change score model', 'sd'] <- 
    Tau_Change_Model_Sd_Per_simulation_sd[[paste0(tau)]]
  
  results_taus_sd[[paste0('tau_effect = ',tau)]]['final score model', 'estimation'] <- 
    Tau_Final_Model_Estimation_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['final score model', 'error'] <- 
    Tau_Final_Model_Error_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['final score model', 'tau'] <- 
    Tau_Final_Model_Tau_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['final score model', 'sd'] <- 
    Tau_Final_Model_Sd_Per_simulation_sd[[paste0(tau)]]
  
  results_taus_sd[[paste0('tau_effect = ',tau)]]['Recovered ANCOVA', 'estimation'] <- 
    Tau_Recover_ANCOVA_Estimation_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['Recovered ANCOVA', 'error'] <- 
    Tau_Recover_ANCOVA_Error_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['Recovered ANCOVA', 'tau'] <- 
    Tau_Recover_ANCOVA_Tau_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['Recovered ANCOVA', 'sd'] <- 
    Tau_Recover_ANCOVA_Sd_Per_simulation_sd[[paste0(tau)]]
  
  results_taus_sd[[paste0('tau_effect = ',tau)]]['Trowman Method', 'estimation'] <- 
    Tau_Trowman_Estimation_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['Trowman Method', 'error'] <- 
    Tau_Trowman_Error_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['Trowman Method', 'tau'] <- 
    Tau_Trowman_Tau_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['Trowman Method', 'sd'] <- 
    Tau_Trowman_Sd_Per_simulation_sd[[paste0(tau)]]
  
  results_taus_sd[[paste0('tau_effect = ',tau)]]['pseudo base model', 'estimation'] <- 
    Tau_Pseudo_Base_Estimation_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['pseudo base model', 'error'] <- 
    Tau_Pseudo_Base_Error_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['pseudo base model', 'tau'] <- 
    Tau_Pseudo_Base_Tau_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['pseudo base model', 'sd'] <- 
    Tau_Pseudo_Base_Sd_Per_simulation_sd[[paste0(tau)]]
  
  results_taus_sd[[paste0('tau_effect = ',tau)]]['pseudo full model', 'estimation'] <- 
    Tau_Pseudo_Full_Estimation_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['pseudo full model', 'error'] <- 
    Tau_Pseudo_Full_Error_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['pseudo full model', 'tau'] <- 
    Tau_Pseudo_Full_Tau_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['pseudo full model', 'sd'] <- 
    Tau_Pseudo_Full_Sd_Per_simulation_sd[[paste0(tau)]]
  
  results_taus_sd[[paste0('tau_effect = ',tau)]]['pseudo two-stage model', 'estimation'] <- 
    Tau_Pseudo_Two_Stage_Estimation_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['pseudo two-stage model', 'error'] <- 
    Tau_Pseudo_Two_Stage_Error_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['pseudo two-stage model', 'tau'] <- 
    Tau_Pseudo_Two_Stage_Tau_Per_simulation_sd[[paste0(tau)]]
  results_taus_sd[[paste0('tau_effect = ',tau)]]['pseudo two-stage model', 'sd'] <- 
    Tau_Pseudo_Two_Stage_Sd_Per_simulation_sd[[paste0(tau)]]
  
}
results_taus_final <- vector(mode = 'list')
for (i in 1:length(results_taus)){
  results_taus_final[[paste0('tau_effect = ', Assess_Range_Tau[i])]] <- 
    cbind(results_taus[[i]], results_taus_sd[[i]])
  names(results_taus_final[[paste0('tau_effect = ', Assess_Range_Tau[i])]]) <- 
    c('estimation', 'sd', 'tau', 'error', 'estimation_sd', 'sd_sd', 'tau_sd', 'error_sd')
}


#----------------------------------------------------------------------------------------------------------------------------
#                       Visualize the errors trend between models
#                       Visualize the estimates distribution with boxplot
#----------------------------------------------------------------------------------------------------------------------------



par(mfrow = c(3,2))
Tau_estimation_results <- vector(mode = 'list')
for (i in 1:length(Assess_Range_Tau)){
  Tau_estimation_results[[paste0('Tau_effect = ', Assess_Range_Tau[i])]] <- 
    cbind(Tau_Change_Score_modelaggregate_estimation[,i],
          Tau_Final_Score_modelaggregate_estimation[,i],
          Tau_Recover_ANCOVA_modelaggregate_estimation[,i],
          Tau_Trowman_modelaggregate_estimation[,i],
          Tau_Pseudo_Base_modelaggregate_estimation[,i],
          Tau_Pseudo_Full_modelaggregate_estimation[,i],
          Tau_Pseudo_Two_Stage_modelaggregate_estimation[,i])
  boxplot.matrix(as.matrix(Tau_estimation_results[[paste0('Tau_effect = ', Assess_Range_Tau[i])]]),
                 xaxt = 'n',
                 xlab = 'Models',
                 ylab = 'estimation distribution',
                 main = paste0('distribution of the estimation for baseline effect = ', 
                               Assess_Range_Tau[i]),
                 col = c('red','blue','yellow','green','brown','orange','purple'))
  axis(1, at = 1:7, labels = c('change model', 'final model', 'RA', 'Trowman', 'PB', 'PF','PTS'))
}
par(mfrow = c(1,1))


