## Code for paper 'A robust Bayesian bias-adjusted random effects model in evidence synthesis'.
## Plot the meta-analysis

## setwd("")

library('rjags')
library('runjags')
library('metafor')
library('gtools')
library('Matrix')

## Data 
## n_studies <- 4  ## number of studies
obs <- matrix(c(10,80,5,17,16,41,16,44), byrow = TRUE, nrow = 4, ncol = 2)
sample_size <- matrix(c(201,298,40,40,122,122,172,170), byrow = TRUE, nrow = 4, ncol = 2) 
parameters <- c("beta", "sigma2_theta", "mu", "OR", "delta", "p")


colnames(obs) = c('Control', 'Treatment')
colnames(sample_size) = c('Control', 'Treatment')
rownames(obs) = c('REFLEX', 'WA16291', 'DANCER', 'SERENE')
rownames(sample_size) = c('REFLEX', 'WA16291', 'DANCER', 'SERENE')

## run JAGS and summarize posterior samples

## Unadjusted bias model

unadjusted_bias_model <- run.jags('bias_adjusted_model.txt', 
                                  data = list(n_studies = 4,
                                              N = sample_size,
                                              obs = obs,
                                              q = c(1,1,1,1),
                                              mu_mu = 0,
                                              sigma2_mu = 10,
                                              mu_beta = 0,
                                              sigma2_beta = 10,
                                              alpha = 0.01,
                                              lambda = 0.01),
                                  monitor = parameters, 
                                  n.chains = 2, burnin = 5000, 
                                  sample = 15000, method = 'rjags')

unadjusted_bias_model

#######################################################################################
## Forestplot
## Log Odds ratio
study_name = c('REFLEX', 'WA16291', 'DANCER', 'SERENE', 'Overall effect')


estimated_mean_unadjusted_model = c(unadjusted_bias_model$summary$statistics['delta[1]','Mean'],
                                    unadjusted_bias_model$summary$statistics['delta[2]','Mean'],
                                    unadjusted_bias_model$summary$statistics['delta[3]','Mean'],
                                    unadjusted_bias_model$summary$statistics['delta[4]','Mean'], 
                                    unadjusted_bias_model$summary$statistics['mu','Mean'])

estimated_lower_ci_unadjusted_model = c(unadjusted_bias_model$summary$quantiles['delta[1]', '2.5%'],
                                        unadjusted_bias_model$summary$quantiles['delta[2]', '2.5%'],
                                        unadjusted_bias_model$summary$quantiles['delta[3]', '2.5%'],
                                        unadjusted_bias_model$summary$quantiles['delta[4]', '2.5%'],
                                        unadjusted_bias_model$summary$quantiles['mu', '2.5%'])

estimated_upper_ci_unadjusted_model = c(unadjusted_bias_model$summary$quantiles['delta[1]', '97.5%'],
                                        unadjusted_bias_model$summary$quantiles['delta[2]', '97.5%'],
                                        unadjusted_bias_model$summary$quantiles['delta[3]', '97.5%'],
                                        unadjusted_bias_model$summary$quantiles['delta[4]', '97.5%'],
                                        unadjusted_bias_model$summary$quantiles['mu', '97.5%'])


## Forestplot of unadjusted model
forest(x = estimated_mean_unadjusted_model, 
       ci.lb = estimated_lower_ci_unadjusted_model, 
       ci.ub = estimated_upper_ci_unadjusted_model,
       xlab = 'Log-Odds Ratio',
       slab = study_name, 
       pch = c(rep(3,4),18),
       psize = 1.5,
       col = c(rep("black",4),"red"),
       header = FALSE,
       alim=c(-0.1,2.5),
       ylim = c(0,13),
       rows = c(10,8,6,4,2))
text(-0.97,12, "Study",     pos=4, font=2, cex=.8)
text(3.4,12, "Estimate [95% PI]",     pos=4, font=2, cex=.8)

#########################################################
## Domain 1 and 2
## Forestplot of unadjusted model
forest(x = estimated_mean_unadjusted_model[1:4], 
       ci.lb = estimated_lower_ci_unadjusted_model[1:4], 
       ci.ub = estimated_upper_ci_unadjusted_model[1:4],
       xlab = 'Log-Odds Ratio',
       slab = study_name[1:4], 
       pch = rep(3,4), 
       psize = 1.5,
       refline = 0,
       col = rep("black",4),
       header = FALSE,
       ylim = c(0,14), 
       alim=c(-0.1,2.5),
       rows = c(11,9,7,5),
       main = "Domain 1 - 2")
text(-0.97,13, "Study", pos=4, font=2, cex=.8)
text(3.31,14, "Estimate [95% PI]", pos=4, font=2, cex=.8)
text(2.63,13, "LB-UB Estimate [LB p2.5; UB p97.5]", col = 'blue', pos=4, font=2, cex=.8)

polygon(x = c(estimated_lower_ci_unadjusted_model[5], estimated_mean_unadjusted_model[5], 
              estimated_upper_ci_unadjusted_model[5], estimated_mean_unadjusted_model[5],
              estimated_lower_ci_unadjusted_model[5]), 
        c(2.5,3,2.5,2,2.5), col = 'black', border = 'black')
text(-0.98, 2.5, labels = paste0('OVERALL EFFECT'), col = "black", pos = 4, cex = 1)

text(3.27, 2.5, labels = paste0(round(estimated_mean_unadjusted_model[5], 2), ' ',
                           '[',round(estimated_lower_ci_unadjusted_model[5], 2),', ',
                           round(estimated_upper_ci_unadjusted_model[5],2),']'), 
     col = "black", pos=4, cex=1)


## Result of RBA domain 1 and 2
#load('results_domain_1.RData')

RBA_lower_mean_adjusted_model_D1 = c(results_domain_1$output_minimum_expectation$summary[1]$statistics['delta[1]', 'Mean'],
                                  results_domain_1$output_minimum_expectation$summary[1]$statistics['delta[2]', 'Mean'],
                                  results_domain_1$output_minimum_expectation$summary[1]$statistics['delta[3]', 'Mean'],
                                  results_domain_1$output_minimum_expectation$summary[1]$statistics['delta[4]', 'Mean'],
                                  results_domain_1$output_minimum_expectation$summary[1]$statistics['mu', 'Mean'])

RBA_upper_mean_adjusted_model_D1 = c(results_domain_1$output_maximum_expectation$summary[1]$statistics['delta[1]', 'Mean'],
                                     results_domain_1$output_maximum_expectation$summary[1]$statistics['delta[2]', 'Mean'],
                                     results_domain_1$output_maximum_expectation$summary[1]$statistics['delta[3]', 'Mean'],
                                     results_domain_1$output_maximum_expectation$summary[1]$statistics['delta[4]', 'Mean'],
                                     results_domain_1$output_maximum_expectation$summary[1]$statistics['mu', 'Mean'])

RBA_lower_ci_adjusted_model_D1 = c(results_domain_1$output_minimum_expectation$summary$quantiles['delta[1]', '2.5%'],
                                   results_domain_1$output_minimum_expectation$summary$quantiles['delta[2]', '2.5%'],
                                   results_domain_1$output_minimum_expectation$summary$quantiles['delta[3]', '2.5%'],
                                   results_domain_1$output_minimum_expectation$summary$quantiles['delta[4]', '2.5%'],
                                   results_domain_1$output_minimum_expectation$summary$quantiles['mu', '2.5%'])

RBA_upper_ci_adjusted_model_D1 = c(results_domain_1$output_maximum_expectation$summary$quantiles['delta[1]', '97.5%'],
                                   results_domain_1$output_maximum_expectation$summary$quantiles['delta[2]', '97.5%'],
                                   results_domain_1$output_maximum_expectation$summary$quantiles['delta[3]', '97.5%'],
                                   results_domain_1$output_maximum_expectation$summary$quantiles['delta[4]', '97.5%'],
                                   results_domain_1$output_maximum_expectation$summary$quantiles['mu', '97.5%'])

## Adding the results of RBA to the foresplot
# Study 1 REFLEX
polygon(x = c(RBA_lower_mean_adjusted_model_D1[1], RBA_upper_mean_adjusted_model_D1[1],RBA_upper_mean_adjusted_model_D1[1], 
              RBA_lower_mean_adjusted_model_D1[1]), c(9.85,9.85,10.15,10.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_D1[1], RBA_upper_ci_adjusted_model_D1[1]), c(10,10), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_D1[1], RBA_lower_ci_adjusted_model_D1[1]), c(9.85,10.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_D1[1], RBA_upper_ci_adjusted_model_D1[1]), c(9.85,10.15), border = 'blue')
text(2.98, 10, labels = paste0('(', min(round(RBA_lower_mean_adjusted_model_D1[1], 2),
                                    round(RBA_upper_mean_adjusted_model_D1[1],2)),'-',
                                max(round(RBA_lower_mean_adjusted_model_D1[1], 2),
                                    round(RBA_upper_mean_adjusted_model_D1[1],2)),')', ' ',
                            '[',round(RBA_lower_ci_adjusted_model_D1[1], 2),', ',
                                round(RBA_upper_ci_adjusted_model_D1[1],2),']'), col = "blue", pos=4, cex=1)

# Study 2 WA16291
polygon(x = c(RBA_lower_mean_adjusted_model_D1[2], RBA_upper_mean_adjusted_model_D1[2],
              RBA_upper_mean_adjusted_model_D1[2], 
              RBA_lower_mean_adjusted_model_D1[2]), c(7.85,7.85,8.15,8.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_D1[2], RBA_upper_ci_adjusted_model_D1[2]), c(8,8), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_D1[2], RBA_lower_ci_adjusted_model_D1[2]), c(7.85,8.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_D1[2], RBA_upper_ci_adjusted_model_D1[2]), c(7.85,8.15), border = 'blue')
text(2.98, 8, labels = paste0('(', min(round(RBA_lower_mean_adjusted_model_D1[2], 2),
                                    round(RBA_upper_mean_adjusted_model_D1[2],2)),'-',
                           max(round(RBA_lower_mean_adjusted_model_D1[2], 2),
                               round(RBA_upper_mean_adjusted_model_D1[2],2)), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D1[2], 2),', ',
                           round(RBA_upper_ci_adjusted_model_D1[2],2),']'), col = "blue", pos=4, cex=1)


# Study 3 DANCER 
polygon(x = c(RBA_lower_mean_adjusted_model_D1[3], RBA_upper_mean_adjusted_model_D1[3],
              RBA_upper_mean_adjusted_model_D1[3], 
              RBA_lower_mean_adjusted_model_D1[3]), c(5.85,5.85,6.15,6.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_D1[3], RBA_upper_ci_adjusted_model_D1[3]), c(6,6), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_D1[3], RBA_lower_ci_adjusted_model_D1[3]), c(5.85,6.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_D1[3], RBA_upper_ci_adjusted_model_D1[3]), c(5.85,6.15), border = 'blue')
text(2.98, 6, labels = paste0('(',min(round(RBA_lower_mean_adjusted_model_D1[3], 2),
                                   round(RBA_upper_mean_adjusted_model_D1[3],2)),'0','-',
                           max(round(RBA_lower_mean_adjusted_model_D1[3], 2),
                               round(RBA_upper_mean_adjusted_model_D1[3],2)), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D1[3], 2),', ',
                           round(RBA_upper_ci_adjusted_model_D1[3],2),']'), col = "blue", pos=4, cex=1)


#Study 4 SERENE
polygon(x = c(RBA_lower_mean_adjusted_model_D1[4], RBA_upper_mean_adjusted_model_D1[4],
              RBA_upper_mean_adjusted_model_D1[4], 
              RBA_lower_mean_adjusted_model_D1[4]), c(3.85,3.85,4.15,4.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_D1[4], RBA_upper_ci_adjusted_model_D1[4]), c(4,4), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_D1[4], RBA_lower_ci_adjusted_model_D1[4]), c(3.85,4.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_D1[4], RBA_upper_ci_adjusted_model_D1[4]), c(3.85,4.15), border = 'blue')
text(2.98, 4, labels = paste0('(',min(round(RBA_lower_mean_adjusted_model_D1[4], 2),
                                   round(RBA_upper_mean_adjusted_model_D1[4],2)),'0','-',
                               max(round(RBA_lower_mean_adjusted_model_D1[4], 2),
                               round(RBA_upper_mean_adjusted_model_D1[4],2)), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D1[4], 2),', ',
                           round(RBA_upper_ci_adjusted_model_D1[4],2),']'), col = "blue", pos=4, cex=1)


polygon(x = c(RBA_lower_ci_adjusted_model_D1[5], RBA_lower_mean_adjusted_model_D1[5], 
              RBA_lower_mean_adjusted_model_D1[5], RBA_lower_ci_adjusted_model_D1[5]), 
        c(1,1.5,0.5,1), col = 'blue', border = 'blue')

polygon(x = c(RBA_lower_mean_adjusted_model_D1[5], RBA_upper_mean_adjusted_model_D1[5], 
              RBA_upper_mean_adjusted_model_D1[5], RBA_lower_mean_adjusted_model_D1[5]), 
        c(1.5,1.5,0.5,0.5), col = 'aquamarine4', border = 'aquamarine4')

polygon(x = c(RBA_upper_mean_adjusted_model_D1[5], RBA_upper_ci_adjusted_model_D1[5], 
              RBA_upper_mean_adjusted_model_D1[5], RBA_upper_mean_adjusted_model_D1[5]), 
        c(1.5,1,0.5,1.5), col = 'blue', border = 'blue')

text(2.98, 1, labels = paste0('(',round(RBA_lower_mean_adjusted_model_D1[5], 2),'-',
                             round(RBA_upper_mean_adjusted_model_D1[5],2), ')', ' ',
                             '[',round(RBA_lower_ci_adjusted_model_D1[5], 2),', ',
                             round(RBA_upper_ci_adjusted_model_D1[5],2),']'), col = "blue", pos=4, cex=1)



############################################################################################################
############################################################################################################
#########################################################
## Domain 3
## Forestplot of unadjusted model
forest(x = estimated_mean_unadjusted_model[1:4], 
       ci.lb = estimated_lower_ci_unadjusted_model[1:4], 
       ci.ub = estimated_upper_ci_unadjusted_model[1:4],
       xlab = 'Log-Odds Ratio',
       slab = study_name[1:4], 
       pch = rep(3,4), 
       psize = 1.5,
       refline = 0,
       col = rep("black",4),
       header = FALSE,
       ylim = c(0,14), 
       alim=c(-0.1,2.5),
       rows = c(11,9,7,5),
       main = "Domain 3")
text(-0.97,13, "Study", pos=4, font=2, cex=.8)
text(3.31,14, "Estimate [95% PI]", pos=4, font=2, cex=.8)
text(2.63,13, "LB-UB Estimate [LB p2.5; UB p97.5]", col = 'blue', pos=4, font=2, cex=.8)


polygon(x = c(estimated_lower_ci_unadjusted_model[5], estimated_mean_unadjusted_model[5], 
              estimated_upper_ci_unadjusted_model[5], estimated_mean_unadjusted_model[5],
              estimated_lower_ci_unadjusted_model[5]), 
        c(2.5,3,2.5,2,2.5), col = 'black', border = 'black')
text(-0.98, 2.5, labels = paste0('OVERALL EFFECT'), col = "black", pos = 4, cex = 1)

text(3.27, 2.5, labels = paste0(round(estimated_mean_unadjusted_model[5], 2), ' ',
                                '[',round(estimated_lower_ci_unadjusted_model[5], 2),', ',
                                round(estimated_upper_ci_unadjusted_model[5],2),']'), 
     col = "black", pos=4, cex=1)

## Result of RBA domain 3
#load('results_domain_3.RData')

RBA_lower_mean_adjusted_model_D3 = c(results_domain_3$output_minimum_expectation$summary[1]$statistics['delta[1]', 'Mean'],
                                     results_domain_3$output_minimum_expectation$summary[1]$statistics['delta[2]', 'Mean'],
                                     results_domain_3$output_minimum_expectation$summary[1]$statistics['delta[3]', 'Mean'],
                                     results_domain_3$output_minimum_expectation$summary[1]$statistics['delta[4]', 'Mean'],
                                     results_domain_3$output_minimum_expectation$summary[1]$statistics['mu', 'Mean'])

RBA_upper_mean_adjusted_model_D3 = c(results_domain_3$output_maximum_expectation$summary[1]$statistics['delta[1]', 'Mean'],
                                     results_domain_3$output_maximum_expectation$summary[1]$statistics['delta[2]', 'Mean'],
                                     results_domain_3$output_maximum_expectation$summary[1]$statistics['delta[3]', 'Mean'],
                                     results_domain_3$output_maximum_expectation$summary[1]$statistics['delta[4]', 'Mean'],
                                     results_domain_3$output_maximum_expectation$summary[1]$statistics['mu', 'Mean'])

RBA_lower_ci_adjusted_model_D3 = c(results_domain_3$output_minimum_expectation$summary$quantiles['delta[1]', '2.5%'],
                                   results_domain_3$output_minimum_expectation$summary$quantiles['delta[2]', '2.5%'],
                                   results_domain_3$output_minimum_expectation$summary$quantiles['delta[3]', '2.5%'],
                                   results_domain_3$output_minimum_expectation$summary$quantiles['delta[4]', '2.5%'],
                                   results_domain_3$output_minimum_expectation$summary$quantiles['mu', '2.5%'])

RBA_upper_ci_adjusted_model_D3 = c(results_domain_3$output_maximum_expectation$summary$quantiles['delta[1]', '97.5%'],
                                   results_domain_3$output_maximum_expectation$summary$quantiles['delta[2]', '97.5%'],
                                   results_domain_3$output_maximum_expectation$summary$quantiles['delta[3]', '97.5%'],
                                   results_domain_3$output_maximum_expectation$summary$quantiles['delta[4]', '97.5%'],
                                   results_domain_3$output_maximum_expectation$summary$quantiles['mu', '97.5%'])

## Adding the results of RBA to the foresplot
# Study 1 REFLEX
polygon(x = c(RBA_lower_mean_adjusted_model_D3[1], RBA_upper_mean_adjusted_model_D3[1],
              RBA_upper_mean_adjusted_model_D3[1], 
              RBA_lower_mean_adjusted_model_D3[1]), c(9.85,9.85,10.15,10.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_D3[1], RBA_upper_ci_adjusted_model_D3[1]), c(10,10), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_D3[1], RBA_lower_ci_adjusted_model_D3[1]), c(9.85,10.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_D3[1], RBA_upper_ci_adjusted_model_D3[1]), c(9.85,10.15), border = 'blue')
text(2.98, 10, labels = paste0('(', min(round(RBA_lower_mean_adjusted_model_D3[1], 2),
                                     round(RBA_upper_mean_adjusted_model_D3[1],2)),'-',
                            max(round(RBA_lower_mean_adjusted_model_D3[1], 2),
                                round(RBA_upper_mean_adjusted_model_D3[1],2)),')', ' ',
                            '[',round(RBA_lower_ci_adjusted_model_D3[1], 2),', ',
                            round(RBA_upper_ci_adjusted_model_D3[1],2),']'), col = "blue", pos=4, cex=1)

# Study 2 WA16291
polygon(x = c(RBA_lower_mean_adjusted_model_D3[2], RBA_upper_mean_adjusted_model_D3[2],
              RBA_upper_mean_adjusted_model_D3[2], 
              RBA_lower_mean_adjusted_model_D3[2]), c(7.85,7.85,8.15,8.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_D3[2], RBA_upper_ci_adjusted_model_D3[2]), c(8,8), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_D3[2], RBA_lower_ci_adjusted_model_D3[2]), c(7.85,8.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_D3[2], RBA_upper_ci_adjusted_model_D3[2]), c(7.85,8.15), border = 'blue')
text(2.98, 8, labels = paste0('(', min(round(RBA_lower_mean_adjusted_model_D3[2], 2),
                                    round(RBA_upper_mean_adjusted_model_D3[2],2)),'-',
                           max(round(RBA_lower_mean_adjusted_model_D3[2], 2),
                               round(RBA_upper_mean_adjusted_model_D3[2],2)), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D3[2], 2),', ',
                           round(RBA_upper_ci_adjusted_model_D3[2],2),']'), col = "blue", pos=4, cex=1)


# Study 3 DANCER 
polygon(x = c(RBA_lower_mean_adjusted_model_D3[3], RBA_upper_mean_adjusted_model_D3[3],
              RBA_upper_mean_adjusted_model_D3[3], 
              RBA_lower_mean_adjusted_model_D3[3]), c(5.85,5.85,6.15,6.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_D3[3], RBA_upper_ci_adjusted_model_D3[3]), c(6,6), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_D3[3], RBA_lower_ci_adjusted_model_D3[3]), c(5.85,6.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_D3[3], RBA_upper_ci_adjusted_model_D3[3]), c(5.85,6.15), border = 'blue')
text(2.98, 6, labels = paste0('(',min(round(RBA_lower_mean_adjusted_model_D3[3], 2),
                                   round(RBA_upper_mean_adjusted_model_D3[3],2)),'-',
                           max(round(RBA_lower_mean_adjusted_model_D3[3], 2),
                               round(RBA_upper_mean_adjusted_model_D3[3],2)), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D3[3], 2),', ',
                           round(RBA_upper_ci_adjusted_model_D3[3],2),'0', ']'), col = "blue", pos=4, cex=1)


#Study 4 SERENE
polygon(x = c(RBA_lower_mean_adjusted_model_D3[4], RBA_upper_mean_adjusted_model_D3[4],
              RBA_upper_mean_adjusted_model_D3[4], 
              RBA_lower_mean_adjusted_model_D3[4]), c(3.85,3.85,4.15,4.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_D3[4], RBA_upper_ci_adjusted_model_D3[4]), c(4,4), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_D3[4], RBA_lower_ci_adjusted_model_D3[4]), c(3.85,4.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_D3[4], RBA_upper_ci_adjusted_model_D3[4]), c(3.85,4.15), border = 'blue')
text(2.98, 4, labels = paste0('(',min(round(RBA_lower_mean_adjusted_model_D3[4], 2),
                                   round(RBA_upper_mean_adjusted_model_D3[4],2)),'-',
                           max(round(RBA_lower_mean_adjusted_model_D3[4], 2),
                               round(RBA_upper_mean_adjusted_model_D3[4],2)), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D3[4], 2),', ',
                           round(RBA_upper_ci_adjusted_model_D3[4],2),']'), col = "blue", pos=4, cex=1)


polygon(x = c(RBA_lower_ci_adjusted_model_D3[5], RBA_lower_mean_adjusted_model_D3[5], 
              RBA_lower_mean_adjusted_model_D3[5], RBA_lower_ci_adjusted_model_D3[5]), 
        c(1,1.5,0.5,1), col = 'blue', border = 'blue')

polygon(x = c(RBA_lower_mean_adjusted_model_D3[5], RBA_upper_mean_adjusted_model_D3[5], 
              RBA_upper_mean_adjusted_model_D3[5], RBA_lower_mean_adjusted_model_D3[5]), 
        c(1.5,1.5,0.5,0.5), col = 'aquamarine4', border = 'aquamarine4')

polygon(x = c(RBA_upper_mean_adjusted_model_D3[5], RBA_upper_ci_adjusted_model_D3[5], 
              RBA_upper_mean_adjusted_model_D3[5], RBA_upper_mean_adjusted_model_D3[5]), 
        c(1.5,1,0.5,1.5), col = 'blue', border = 'blue')

text(2.98, 1, labels = paste0('(',round(RBA_lower_mean_adjusted_model_D3[5], 2),'-',
                           round(RBA_upper_mean_adjusted_model_D3[5],2), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D3[5], 2),', ',
                           round(RBA_upper_ci_adjusted_model_D3[5],2),']'), col = "blue", pos=4, cex=1)


############################################################################################################
############################################################################################################
#########################################################
## Domain 4
## Forestplot of unadjusted model
forest(x = estimated_mean_unadjusted_model[1:4], 
       ci.lb = estimated_lower_ci_unadjusted_model[1:4], 
       ci.ub = estimated_upper_ci_unadjusted_model[1:4],
       xlab = 'Log-Odds Ratio',
       slab = study_name[1:4], 
       pch = rep(3,4), 
       psize = 1.5,
       refline = 0,
       col = rep("black",4),
       header = FALSE,
       ylim = c(0,14), 
       alim=c(-0.1,2.5),
       rows = c(11,9,7,5),
       main = "Domain 4")
text(-0.97,13, "Study", pos=4, font=2, cex=.8)
text(3.31,14, "Estimate [95% PI]", pos=4, font=2, cex=.8)
text(2.63,13, "LB-UB Estimate [LB p2.5; UB p97.5]", col = 'blue', pos=4, font=2, cex=.8)

polygon(x = c(estimated_lower_ci_unadjusted_model[5], estimated_mean_unadjusted_model[5], 
              estimated_upper_ci_unadjusted_model[5], estimated_mean_unadjusted_model[5],
              estimated_lower_ci_unadjusted_model[5]), 
        c(2.5,3,2.5,2,2.5), col = 'black', border = 'black')
text(-0.98, 2.5, labels = paste0('OVERALL EFFECT'), col = "black", pos = 4, cex = 1)

text(3.27, 2.5, labels = paste0(round(estimated_mean_unadjusted_model[5], 2), ' ',
                                '[',round(estimated_lower_ci_unadjusted_model[5], 2),', ',
                                round(estimated_upper_ci_unadjusted_model[5],2),']'), 
     col = "black", pos=4, cex=1)


## Result of RBA domain 4
#load('results_domain_4.RData')

RBA_lower_mean_adjusted_model_D4 = c(results_domain_4$output_minimum_expectation$summary[1]$statistics['delta[1]', 'Mean'],
                                     results_domain_4$output_minimum_expectation$summary[1]$statistics['delta[2]', 'Mean'],
                                     results_domain_4$output_minimum_expectation$summary[1]$statistics['delta[3]', 'Mean'],
                                     results_domain_4$output_minimum_expectation$summary[1]$statistics['delta[4]', 'Mean'],
                                     results_domain_4$output_minimum_expectation$summary[1]$statistics['mu', 'Mean'])

RBA_upper_mean_adjusted_model_D4 = c(results_domain_4$output_maximum_expectation$summary[1]$statistics['delta[1]', 'Mean'],
                                     results_domain_4$output_maximum_expectation$summary[1]$statistics['delta[2]', 'Mean'],
                                     results_domain_4$output_maximum_expectation$summary[1]$statistics['delta[3]', 'Mean'],
                                     results_domain_4$output_maximum_expectation$summary[1]$statistics['delta[4]', 'Mean'],
                                     results_domain_4$output_maximum_expectation$summary[1]$statistics['mu', 'Mean'])

RBA_lower_ci_adjusted_model_D4 = c(results_domain_4$output_minimum_expectation$summary$quantiles['delta[1]', '2.5%'],
                                   results_domain_4$output_minimum_expectation$summary$quantiles['delta[2]', '2.5%'],
                                   results_domain_4$output_minimum_expectation$summary$quantiles['delta[3]', '2.5%'],
                                   results_domain_4$output_minimum_expectation$summary$quantiles['delta[4]', '2.5%'],
                                   results_domain_4$output_minimum_expectation$summary$quantiles['mu', '2.5%'])

RBA_upper_ci_adjusted_model_D4 = c(results_domain_4$output_maximum_expectation$summary$quantiles['delta[1]', '97.5%'],
                                   results_domain_4$output_maximum_expectation$summary$quantiles['delta[2]', '97.5%'],
                                   results_domain_4$output_maximum_expectation$summary$quantiles['delta[3]', '97.5%'],
                                   results_domain_4$output_maximum_expectation$summary$quantiles['delta[4]', '97.5%'],
                                   results_domain_4$output_maximum_expectation$summary$quantiles['mu', '97.5%'])

## Adding the results of RBA to the foresplot
# Study 1 REFLEX
polygon(x = c(RBA_lower_mean_adjusted_model_D4[1], RBA_upper_mean_adjusted_model_D4[1],
              RBA_upper_mean_adjusted_model_D4[1], 
              RBA_lower_mean_adjusted_model_D4[1]), c(9.85,9.85,10.15,10.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_D4[1], RBA_upper_ci_adjusted_model_D4[1]), c(10,10), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_D4[1], RBA_lower_ci_adjusted_model_D4[1]), c(9.85,10.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_D4[1], RBA_upper_ci_adjusted_model_D4[1]), c(9.85,10.15), border = 'blue')
text(2.98, 10, labels = paste0('(', min(round(RBA_lower_mean_adjusted_model_D4[1], 2),
                                     round(RBA_upper_mean_adjusted_model_D4[1],2)), '0', '-',
                            max(round(RBA_lower_mean_adjusted_model_D4[1], 2),
                                round(RBA_upper_mean_adjusted_model_D4[1],2)),')', ' ',
                            '[',round(RBA_lower_ci_adjusted_model_D4[1], 2),', ',
                            round(RBA_upper_ci_adjusted_model_D4[1],2),']'), col = "blue", pos=4, cex=1)

# Study 2 WA16291
polygon(x = c(RBA_lower_mean_adjusted_model_D4[2], RBA_upper_mean_adjusted_model_D4[2],
              RBA_upper_mean_adjusted_model_D4[2], 
              RBA_lower_mean_adjusted_model_D4[2]), c(7.85,7.85,8.15,8.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_D4[2], RBA_upper_ci_adjusted_model_D4[2]), c(8,8), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_D4[2], RBA_lower_ci_adjusted_model_D4[2]), c(7.85,8.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_D4[2], RBA_upper_ci_adjusted_model_D4[2]), c(7.85,8.15), border = 'blue')
text(2.98, 8, labels = paste0('(', min(round(RBA_lower_mean_adjusted_model_D4[2], 2),
                                    round(RBA_upper_mean_adjusted_model_D4[2],2)),'-',
                           max(round(RBA_lower_mean_adjusted_model_D4[2], 2),
                               round(RBA_upper_mean_adjusted_model_D4[2],2)), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D4[2], 2),', ',
                           round(RBA_upper_ci_adjusted_model_D4[2],2),']'), col = "blue", pos=4, cex=1)


# Study 3 DANCER 
polygon(x = c(RBA_lower_mean_adjusted_model_D4[3], RBA_upper_mean_adjusted_model_D4[3],
              RBA_upper_mean_adjusted_model_D4[3], 
              RBA_lower_mean_adjusted_model_D4[3]), c(5.85,5.85,6.15,6.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_D4[3], RBA_upper_ci_adjusted_model_D4[3]), c(6,6), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_D4[3], RBA_lower_ci_adjusted_model_D4[3]), c(5.85,6.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_D4[3], RBA_upper_ci_adjusted_model_D4[3]), c(5.85,6.15), border = 'blue')
text(2.98, 6, labels = paste0('(',min(round(RBA_lower_mean_adjusted_model_D4[3], 2),
                                   round(RBA_upper_mean_adjusted_model_D4[3],2)),'-',
                           max(round(RBA_lower_mean_adjusted_model_D4[3], 2),
                               round(RBA_upper_mean_adjusted_model_D4[3],2)), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D4[3], 2),', ',
                           round(RBA_upper_ci_adjusted_model_D4[3],2),']'), col = "blue", pos=4, cex=1)


#Study 4 SERENE
polygon(x = c(RBA_lower_mean_adjusted_model_D4[4], RBA_upper_mean_adjusted_model_D4[4],
              RBA_upper_mean_adjusted_model_D4[4], 
              RBA_lower_mean_adjusted_model_D4[4]), c(3.85,3.85,4.15,4.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_D4[4], RBA_upper_ci_adjusted_model_D4[4]), c(4,4), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_D4[4], RBA_lower_ci_adjusted_model_D4[4]), c(3.85,4.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_D4[4], RBA_upper_ci_adjusted_model_D4[4]), c(3.85,4.15), border = 'blue')
text(2.98, 4, labels = paste0('(',min(round(RBA_lower_mean_adjusted_model_D4[4], 2),
                                   round(RBA_upper_mean_adjusted_model_D4[4],2)),'-',
                           max(round(RBA_lower_mean_adjusted_model_D4[4], 2),
                               round(RBA_upper_mean_adjusted_model_D4[4],2)), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D4[4], 2), '0' ,', ',
                           round(RBA_upper_ci_adjusted_model_D4[4],2),']'), col = "blue", pos=4, cex=1)


polygon(x = c(RBA_lower_ci_adjusted_model_D4[5], RBA_lower_mean_adjusted_model_D4[5], 
              RBA_lower_mean_adjusted_model_D4[5], RBA_lower_ci_adjusted_model_D4[5]), 
        c(1,1.5,0.5,1), col = 'blue', border = 'blue')

polygon(x = c(RBA_lower_mean_adjusted_model_D4[5], RBA_upper_mean_adjusted_model_D4[5], 
              RBA_upper_mean_adjusted_model_D4[5], RBA_lower_mean_adjusted_model_D4[5]), 
        c(1.5,1.5,0.5,0.5), col = 'aquamarine4', border = 'aquamarine4')

polygon(x = c(RBA_upper_mean_adjusted_model_D4[5], RBA_upper_ci_adjusted_model_D4[5], 
              RBA_upper_mean_adjusted_model_D4[5], RBA_upper_mean_adjusted_model_D4[5]), 
        c(1.5,1,0.5,1.5), col = 'blue', border = 'blue')

text(2.98, 1, labels = paste0('(',round(RBA_lower_mean_adjusted_model_D4[5], 2),'-',
                           round(RBA_upper_mean_adjusted_model_D4[5],2), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D4[5], 2),', ',
                           round(RBA_upper_ci_adjusted_model_D4[5],2),']'), col = "blue", pos=4, cex=1)


############################################################################################################
############################################################################################################
#########################################################
## Domain 5 and 6

## Forestplot of unadjusted model
forest(x = estimated_mean_unadjusted_model[1:4], 
       ci.lb = estimated_lower_ci_unadjusted_model[1:4], 
       ci.ub = estimated_upper_ci_unadjusted_model[1:4],
       xlab = 'Log-Odds Ratio',
       slab = study_name[1:4], 
       pch = rep(3,4), 
       psize = 1.5,
       refline = 0,
       col = rep("black",4),
       header = FALSE, 
       ylim = c(0,14), 
       alim=c(-0.1,2.5),
       rows = c(11,9,7,5),
       main = "Domain 5 - 6")
text(-0.97,13, "Study", pos=4, font=2, cex=.8)
text(3.31,14, "Estimate [95% PI]", pos=4, font=2, cex=.8)
text(2.63,13, "LB-UB Estimate [LB p2.5; UB p97.5]", col = 'blue', pos=4, font=2, cex=.8)


polygon(x = c(estimated_lower_ci_unadjusted_model[5], estimated_mean_unadjusted_model[5], 
              estimated_upper_ci_unadjusted_model[5], estimated_mean_unadjusted_model[5],
              estimated_lower_ci_unadjusted_model[5]), 
        c(2.5,3,2.5,2,2.5), col = 'black', border = 'black')
text(-0.98, 2.5, labels = paste0('OVERALL EFFECT'), col = "black", pos = 4, cex = 1)

text(3.27, 2.5, labels = paste0(round(estimated_mean_unadjusted_model[5], 2), ' ',
                                '[',round(estimated_lower_ci_unadjusted_model[5], 2),', ',
                                round(estimated_upper_ci_unadjusted_model[5],2),']'), 
     col = "black", pos=4, cex=1)

## Result of RBA domain 5 and 6
#load('results_domain_5.RData')

RBA_lower_mean_adjusted_model_D5 = c(results_domain_5$output_minimum_expectation$summary[1]$statistics['delta[1]', 'Mean'],
                                     results_domain_5$output_minimum_expectation$summary[1]$statistics['delta[2]', 'Mean'],
                                     results_domain_5$output_minimum_expectation$summary[1]$statistics['delta[3]', 'Mean'],
                                     results_domain_5$output_minimum_expectation$summary[1]$statistics['delta[4]', 'Mean'],
                                     results_domain_5$output_minimum_expectation$summary[1]$statistics['mu', 'Mean'])

RBA_upper_mean_adjusted_model_D5 = c(results_domain_5$output_maximum_expectation$summary[1]$statistics['delta[1]', 'Mean'],
                                     results_domain_5$output_maximum_expectation$summary[1]$statistics['delta[2]', 'Mean'],
                                     results_domain_5$output_maximum_expectation$summary[1]$statistics['delta[3]', 'Mean'],
                                     results_domain_5$output_maximum_expectation$summary[1]$statistics['delta[4]', 'Mean'],
                                     results_domain_5$output_maximum_expectation$summary[1]$statistics['mu', 'Mean'])

RBA_lower_ci_adjusted_model_D5 = c(results_domain_5$output_minimum_expectation$summary$quantiles['delta[1]', '2.5%'],
                                   results_domain_5$output_minimum_expectation$summary$quantiles['delta[2]', '2.5%'],
                                   results_domain_5$output_minimum_expectation$summary$quantiles['delta[3]', '2.5%'],
                                   results_domain_5$output_minimum_expectation$summary$quantiles['delta[4]', '2.5%'],
                                   results_domain_5$output_minimum_expectation$summary$quantiles['mu', '2.5%'])

RBA_upper_ci_adjusted_model_D5 = c(results_domain_5$output_maximum_expectation$summary$quantiles['delta[1]', '97.5%'],
                                   results_domain_5$output_maximum_expectation$summary$quantiles['delta[2]', '97.5%'],
                                   results_domain_5$output_maximum_expectation$summary$quantiles['delta[3]', '97.5%'],
                                   results_domain_5$output_maximum_expectation$summary$quantiles['delta[4]', '97.5%'],
                                   results_domain_5$output_maximum_expectation$summary$quantiles['mu', '97.5%'])

## Adding the results of RBA to the foresplot
# Study 1 REFLEX
polygon(x = c(RBA_lower_mean_adjusted_model_D5[1], RBA_upper_mean_adjusted_model_D5[1],
              RBA_upper_mean_adjusted_model_D5[1], 
              RBA_lower_mean_adjusted_model_D5[1]), c(9.85,9.85,10.15,10.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_D5[1], RBA_upper_ci_adjusted_model_D5[1]), c(10,10), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_D5[1], RBA_lower_ci_adjusted_model_D5[1]), c(9.85,10.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_D5[1], RBA_upper_ci_adjusted_model_D5[1]), c(9.85,10.15), border = 'blue')
text(2.98, 10, labels = paste0('(', min(round(RBA_lower_mean_adjusted_model_D5[1], 2),
                                     round(RBA_upper_mean_adjusted_model_D5[1],2)),'-',
                            max(round(RBA_lower_mean_adjusted_model_D5[1], 2),
                                round(RBA_upper_mean_adjusted_model_D5[1],2)),')', ' ',
                            '[',round(RBA_lower_ci_adjusted_model_D5[1], 2),', ',
                            round(RBA_upper_ci_adjusted_model_D5[1],2),']'), col = "blue", pos=4, cex=1)

# Study 2 WA16291
polygon(x = c(RBA_lower_mean_adjusted_model_D5[2], RBA_upper_mean_adjusted_model_D5[2],
              RBA_upper_mean_adjusted_model_D5[2], 
              RBA_lower_mean_adjusted_model_D5[2]), c(7.85,7.85,8.15,8.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_D5[2], RBA_upper_ci_adjusted_model_D5[2]), c(8,8), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_D5[2], RBA_lower_ci_adjusted_model_D5[2]), c(7.85,8.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_D5[2], RBA_upper_ci_adjusted_model_D5[2]), c(7.85,8.15), border = 'blue')
text(2.98, 8, labels = paste0('(', min(round(RBA_lower_mean_adjusted_model_D5[2], 2),
                                    round(RBA_upper_mean_adjusted_model_D5[2],2)),'-',
                           max(round(RBA_lower_mean_adjusted_model_D5[2], 2),
                               round(RBA_upper_mean_adjusted_model_D5[2],2)), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D5[2], 2),', ',
                           round(RBA_upper_ci_adjusted_model_D5[2],2),']'), col = "blue", pos=4, cex=1)


# Study 3 DANCER 
polygon(x = c(RBA_lower_mean_adjusted_model_D5[3], RBA_upper_mean_adjusted_model_D4[3],
              RBA_upper_mean_adjusted_model_D5[3], 
              RBA_lower_mean_adjusted_model_D5[3]), c(5.85,5.85,6.15,6.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_D5[3], RBA_upper_ci_adjusted_model_D5[3]), c(6,6), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_D5[3], RBA_lower_ci_adjusted_model_D5[3]), c(5.85,6.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_D5[3], RBA_upper_ci_adjusted_model_D5[3]), c(5.85,6.15), border = 'blue')
text(2.98, 6, labels = paste0('(',min(round(RBA_lower_mean_adjusted_model_D5[3], 2),
                                   round(RBA_upper_mean_adjusted_model_D5[3],2)),'-',
                           max(round(RBA_lower_mean_adjusted_model_D5[3], 2),
                               round(RBA_upper_mean_adjusted_model_D5[3],2)), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D5[3], 2),'0' , ', ',
                           round(RBA_upper_ci_adjusted_model_D5[3],2),']'), col = "blue", pos=4, cex=1)


#Study 4 SERENE
polygon(x = c(RBA_lower_mean_adjusted_model_D5[4], RBA_upper_mean_adjusted_model_D5[4],
              RBA_upper_mean_adjusted_model_D5[4], 
              RBA_lower_mean_adjusted_model_D5[4]), c(3.85,3.85,4.15,4.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_D5[4], RBA_upper_ci_adjusted_model_D5[4]), c(4,4), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_D5[4], RBA_lower_ci_adjusted_model_D5[4]), c(3.85,4.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_D5[4], RBA_upper_ci_adjusted_model_D5[4]), c(3.85,4.15), border = 'blue')
text(2.98, 4, labels = paste0('(',min(round(RBA_lower_mean_adjusted_model_D5[4], 2),
                                   round(RBA_upper_mean_adjusted_model_D5[4],2)),'-',
                           max(round(RBA_lower_mean_adjusted_model_D5[4], 2),
                               round(RBA_upper_mean_adjusted_model_D5[4],2)), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D5[4], 2),', ',
                           round(RBA_upper_ci_adjusted_model_D5[4],2),']'), col = "blue", pos=4, cex=1)


polygon(x = c(RBA_lower_ci_adjusted_model_D5[5], RBA_lower_mean_adjusted_model_D5[5], 
              RBA_lower_mean_adjusted_model_D5[5], RBA_lower_ci_adjusted_model_D5[5]), 
        c(1,1.5,0.5,1), col = 'blue', border = 'blue')

polygon(x = c(RBA_lower_mean_adjusted_model_D5[5], RBA_upper_mean_adjusted_model_D5[5], 
              RBA_upper_mean_adjusted_model_D5[5], RBA_lower_mean_adjusted_model_D5[5]), 
        c(1.5,1.5,0.5,0.5), col = 'aquamarine4', border = 'aquamarine4')

polygon(x = c(RBA_upper_mean_adjusted_model_D5[5], RBA_upper_ci_adjusted_model_D5[5], 
              RBA_upper_mean_adjusted_model_D5[5], RBA_upper_mean_adjusted_model_D5[5]), 
        c(1.5,1,0.5,1.5), col = 'blue', border = 'blue')

text(2.98, 1, labels = paste0('(',round(RBA_lower_mean_adjusted_model_D5[5], 2),'-',
                           round(RBA_upper_mean_adjusted_model_D5[5],2), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D5[5], 2),', ',
                           round(RBA_upper_ci_adjusted_model_D5[5],2),']'), col = "blue", pos=4, cex=1)


############################################################################################################
############################################################################################################
#########################################################
## All domains 
## Forestplot of unadjusted model
forest(x = estimated_mean_unadjusted_model[1:4], 
       ci.lb = estimated_lower_ci_unadjusted_model[1:4], 
       ci.ub = estimated_upper_ci_unadjusted_model[1:4],
       xlab = 'Log-Odds Ratio',
       slab = study_name[1:4], 
       pch = rep(3,4), 
       psize = 1.5,
       refline = 0,
       col = rep("black",4),
       header = FALSE, 
       ylim = c(0,14), 
       alim=c(-0.1,2.5),
       rows = c(11,9,7,5),
       main = "All domains")
text(-0.97,13, "Study", pos=4, font=2, cex=.8)
text(3.31,14, "Estimate [95% PI]", pos=4, font=2, cex=.8)
text(2.63,13, "LB-UB Estimate [LB p2.5; UB p97.5]", col = 'blue', pos=4, font=2, cex=.8)


polygon(x = c(estimated_lower_ci_unadjusted_model[5], estimated_mean_unadjusted_model[5], 
              estimated_upper_ci_unadjusted_model[5], estimated_mean_unadjusted_model[5],
              estimated_lower_ci_unadjusted_model[5]), 
        c(2.5,3,2.5,2,2.5), col = 'black', border = 'black')
text(-0.98, 2.5, labels = paste0('OVERALL EFFECT'), col = "black", pos = 4, cex = 1)

text(3.27, 2.5, labels = paste0(round(estimated_mean_unadjusted_model[5], 2), ' ',
                                '[',round(estimated_lower_ci_unadjusted_model[5], 2),', ',
                                round(estimated_upper_ci_unadjusted_model[5],2),']'), 
     col = "black", pos=4, cex=1)


## Result of RBA all domains
#load('results_all_domains.RData')

RBA_lower_mean_adjusted_model_all_domains = c(results_all_domains$output_minimum_expectation$summary[1]$statistics['delta[1]', 'Mean'],
                                              results_all_domains$output_minimum_expectation$summary[1]$statistics['delta[2]', 'Mean'],
                                              results_all_domains$output_minimum_expectation$summary[1]$statistics['delta[3]', 'Mean'],
                                              results_all_domains$output_minimum_expectation$summary[1]$statistics['delta[4]', 'Mean'],
                                              results_all_domains$output_minimum_expectation$summary[1]$statistics['mu', 'Mean'])

RBA_upper_mean_adjusted_model_all_domains = c(results_all_domains$output_maximum_expectation$summary[1]$statistics['delta[1]', 'Mean'],
                                              results_all_domains$output_maximum_expectation$summary[1]$statistics['delta[2]', 'Mean'],
                                              results_all_domains$output_maximum_expectation$summary[1]$statistics['delta[3]', 'Mean'],
                                              results_all_domains$output_maximum_expectation$summary[1]$statistics['delta[4]', 'Mean'],
                                              results_all_domains$output_maximum_expectation$summary[1]$statistics['mu', 'Mean'])

RBA_lower_ci_adjusted_model_all_domains = c(results_all_domains$output_minimum_expectation$summary$quantiles['delta[1]', '2.5%'],
                                            results_all_domains$output_minimum_expectation$summary$quantiles['delta[2]', '2.5%'],
                                            results_all_domains$output_minimum_expectation$summary$quantiles['delta[3]', '2.5%'],
                                            results_all_domains$output_minimum_expectation$summary$quantiles['delta[4]', '2.5%'],
                                            results_all_domains$output_minimum_expectation$summary$quantiles['mu', '2.5%'])

RBA_upper_ci_adjusted_model_all_domains = c(results_all_domains$output_maximum_expectation$summary$quantiles['delta[1]', '97.5%'],
                                            results_all_domains$output_maximum_expectation$summary$quantiles['delta[2]', '97.5%'],
                                            results_all_domains$output_maximum_expectation$summary$quantiles['delta[3]', '97.5%'],
                                            results_all_domains$output_maximum_expectation$summary$quantiles['delta[4]', '97.5%'],
                                            results_all_domains$output_maximum_expectation$summary$quantiles['mu', '97.5%'])

## Adding the results of RBA to the foresplot
# Study 1 REFLEX
polygon(x = c(RBA_lower_mean_adjusted_model_all_domains[1], RBA_upper_mean_adjusted_model_all_domains[1],
              RBA_upper_mean_adjusted_model_all_domains[1], 
              RBA_lower_mean_adjusted_model_all_domains[1]), c(9.85,9.85,10.15,10.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_all_domains[1], RBA_upper_ci_adjusted_model_all_domains[1]), c(10,10), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_all_domains[1], RBA_lower_ci_adjusted_model_all_domains[1]), c(9.85,10.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_all_domains[1], RBA_upper_ci_adjusted_model_all_domains[1]), c(9.85,10.15), border = 'blue')
text(2.98, 10, labels = paste0('(', min(round(RBA_lower_mean_adjusted_model_all_domains[1], 2),
                                     round(RBA_upper_mean_adjusted_model_all_domains[1],2)),'-',
                            max(round(RBA_lower_mean_adjusted_model_all_domains[1], 2),
                                round(RBA_upper_mean_adjusted_model_all_domains[1],2)),'0',')', ' ',
                            '[',round(RBA_lower_ci_adjusted_model_all_domains[1], 2),', ',
                            round(RBA_upper_ci_adjusted_model_all_domains[1],2),']'), col = "blue", pos=4, cex=1)

# Study 2 WA16291
polygon(x = c(RBA_lower_mean_adjusted_model_all_domains[2], RBA_upper_mean_adjusted_model_all_domains[2],
              RBA_upper_mean_adjusted_model_all_domains[2], 
              RBA_lower_mean_adjusted_model_all_domains[2]), c(7.85,7.85,8.15,8.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_all_domains[2], RBA_upper_ci_adjusted_model_all_domains[2]), c(8,8), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_all_domains[2], RBA_lower_ci_adjusted_model_all_domains[2]), c(7.85,8.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_all_domains[2], RBA_upper_ci_adjusted_model_all_domains[2]), c(7.85,8.15), border = 'blue')
text(2.98, 8, labels = paste0('(', min(round(RBA_lower_mean_adjusted_model_all_domains[2], 2),
                                    round(RBA_upper_mean_adjusted_model_all_domains[2],2)),'-',
                           max(round(RBA_lower_mean_adjusted_model_all_domains[2], 2),
                               round(RBA_upper_mean_adjusted_model_all_domains[2],2)), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_all_domains[2], 2),'0',', ',
                           round(RBA_upper_ci_adjusted_model_all_domains[2],2),']'), col = "blue", pos=4, cex=1)


# Study 3 DANCER 
polygon(x = c(RBA_lower_mean_adjusted_model_all_domains[3], RBA_upper_mean_adjusted_model_all_domains[3],
              RBA_upper_mean_adjusted_model_all_domains[3], 
              RBA_lower_mean_adjusted_model_all_domains[3]), c(5.85,5.85,6.15,6.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_all_domains[3], RBA_upper_ci_adjusted_model_all_domains[3]), c(6,6), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_all_domains[3], RBA_lower_ci_adjusted_model_all_domains[3]), c(5.85,6.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_all_domains[3], RBA_upper_ci_adjusted_model_all_domains[3]), c(5.85,6.15), border = 'blue')
text(2.98, 6, labels = paste0('(',min(round(RBA_lower_mean_adjusted_model_all_domains[3], 2),
                                   round(RBA_upper_mean_adjusted_model_all_domains[3],2)),'-',
                           max(round(RBA_lower_mean_adjusted_model_all_domains[3], 2),
                               round(RBA_upper_mean_adjusted_model_all_domains[3],2)), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_all_domains[3], 2),', ',
                           round(RBA_upper_ci_adjusted_model_all_domains[3],2),']'), col = "blue", pos=4, cex=1)


#Study 4 SERENE
polygon(x = c(RBA_lower_mean_adjusted_model_all_domains[4], RBA_upper_mean_adjusted_model_all_domains[4],
              RBA_upper_mean_adjusted_model_all_domains[4], 
              RBA_lower_mean_adjusted_model_all_domains[4]), c(3.85,3.85,4.15,4.15), col = 'aquamarine4', border = 'aquamarine4')
polygon(x = c(RBA_lower_ci_adjusted_model_all_domains[4], RBA_upper_ci_adjusted_model_all_domains[4]), c(4,4), border = 'blue')
polygon(x = c(RBA_lower_ci_adjusted_model_all_domains[4], RBA_lower_ci_adjusted_model_all_domains[4]), c(3.85,4.15), border = 'blue')
polygon(x = c(RBA_upper_ci_adjusted_model_all_domains[4], RBA_upper_ci_adjusted_model_all_domains[4]), c(3.85,4.15), border = 'blue')
text(2.98, 4, labels = paste0('(',min(round(RBA_lower_mean_adjusted_model_all_domains[4], 2),
                                   round(RBA_upper_mean_adjusted_model_all_domains[4],2)),'-',
                           max(round(RBA_lower_mean_adjusted_model_all_domains[4], 2),
                               round(RBA_upper_mean_adjusted_model_all_domains[4],2)), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_all_domains[4], 2),', ',
                           round(RBA_upper_ci_adjusted_model_all_domains[4],2),']'), col = "blue", pos=4, cex=1)


polygon(x = c(RBA_lower_ci_adjusted_model_all_domains[5], RBA_lower_mean_adjusted_model_all_domains[5], 
              RBA_lower_mean_adjusted_model_all_domains[5], RBA_lower_ci_adjusted_model_all_domains[5]), 
        c(1,1.5,0.5,1), col = 'blue', border = 'blue')

polygon(x = c(RBA_lower_mean_adjusted_model_all_domains[5], RBA_upper_mean_adjusted_model_all_domains[5], 
              RBA_upper_mean_adjusted_model_all_domains[5], RBA_lower_mean_adjusted_model_all_domains[5]), 
        c(1.5,1.5,0.5,0.5), col = 'aquamarine4', border = 'aquamarine4')

polygon(x = c(RBA_upper_mean_adjusted_model_all_domains[5], RBA_upper_ci_adjusted_model_all_domains[5], 
              RBA_upper_mean_adjusted_model_all_domains[5], RBA_upper_mean_adjusted_model_all_domains[5]), 
        c(1.5,1,0.5,1.5), col = 'blue', border = 'blue')

text(2.98, 1, labels = paste0('(',round(RBA_lower_mean_adjusted_model_all_domains[5], 2),'-',
                           round(RBA_upper_mean_adjusted_model_all_domains[5],2), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_all_domains[5], 2),'0',', ',
                           round(RBA_upper_ci_adjusted_model_all_domains[5],2),']'), col = "blue", pos=4, cex=1)


############################################################################################################
############################################################################################################
#########################################################
## Overall effect per domains 

forest(x = estimated_mean_unadjusted_model[1:4], 
       ci.lb = estimated_lower_ci_unadjusted_model[1:4], 
       ci.ub = estimated_upper_ci_unadjusted_model[1:4],
       xlab = 'Log-Odds Ratio',
       slab = study_name[1:4], 
       pch = rep(3,4), 
       psize = 1.5,
       refline = 0,
       col = rep("white",4),
       header = FALSE, 
       ylim = c(0,14), 
       alim=c(-0.1,2.5),
       rows = c(11,9, 7,5),
       main = "Overall effect")
text(-0.97,13, "Bias Domain", pos=4, font=2, cex=.8)
text(3.31,14, "Estimate [95% PI]", pos=4, font=2, cex=.8)
text(2.63,13, "LB-UB Estimate [LB p2.5; UB p97.5]", col = 'blue', pos=4, font=2, cex=.8)


## Adding the results of RBA to the foresplot
# Domain 1-2

polygon(x = c(RBA_lower_ci_adjusted_model_D1[5], RBA_lower_mean_adjusted_model_D1[5], 
              RBA_lower_mean_adjusted_model_D1[5], RBA_lower_ci_adjusted_model_D1[5]), 
        c(11,11.5,10.5,11), col = 'blue', border = 'blue')

polygon(x = c(RBA_lower_mean_adjusted_model_D1[5], RBA_upper_mean_adjusted_model_D1[5], 
              RBA_upper_mean_adjusted_model_D1[5], RBA_lower_mean_adjusted_model_D1[5]), 
        c(11.5,11.5,10.5,10.5), col = 'aquamarine4', border = 'aquamarine4')

polygon(x = c(RBA_upper_mean_adjusted_model_D1[5], RBA_upper_ci_adjusted_model_D1[5], 
              RBA_upper_mean_adjusted_model_D1[5], RBA_upper_mean_adjusted_model_D1[5]), 
        c(11.5,11,10.5,11.5), col = 'blue', border = 'blue')

text(2.98, 11, labels = paste0('(',round(RBA_lower_mean_adjusted_model_D1[5], 2),'-',
                           round(RBA_upper_mean_adjusted_model_D1[5],2), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D1[5], 2),', ',
                           round(RBA_upper_ci_adjusted_model_D1[5],2),']'), col = "blue", pos=4, cex=1)

text(-0.95, 11, labels = paste0('1-2'), col = "black", pos = 4, cex = 1)

##############
# Domain 3

polygon(x = c(RBA_lower_ci_adjusted_model_D3[5], RBA_lower_mean_adjusted_model_D3[5], 
              RBA_lower_mean_adjusted_model_D3[5], RBA_lower_ci_adjusted_model_D3[5]), 
        c(9,9.5,8.5,9), col = 'blue', border = 'blue')

polygon(x = c(RBA_lower_mean_adjusted_model_D3[5], RBA_upper_mean_adjusted_model_D3[5], 
              RBA_upper_mean_adjusted_model_D3[5], RBA_lower_mean_adjusted_model_D3[5]), 
        c(9.5,9.5,8.5,8.5), col = 'aquamarine4', border = 'aquamarine4')

polygon(x = c(RBA_upper_mean_adjusted_model_D3[5], RBA_upper_ci_adjusted_model_D3[5], 
              RBA_upper_mean_adjusted_model_D3[5], RBA_upper_mean_adjusted_model_D3[5]), 
        c(9.5,9,8.5,9.5), col = 'blue', border = 'blue')

text(2.98, 9, labels = paste0('(',round(RBA_lower_mean_adjusted_model_D3[5], 2),'-',
                           round(RBA_upper_mean_adjusted_model_D3[5],2), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D3[5], 2),', ',
                           round(RBA_upper_ci_adjusted_model_D3[5],2),']'), col = "blue", pos=4, cex=1)

text(-0.95, 9, labels = paste0('3'), col = "black", pos = 4, cex = 1)


# Domain 4

polygon(x = c(RBA_lower_ci_adjusted_model_D4[5], RBA_lower_mean_adjusted_model_D4[5], 
              RBA_lower_mean_adjusted_model_D4[5], RBA_lower_ci_adjusted_model_D4[5]), 
        c(7,7.5,6.5,7), col = 'blue', border = 'blue')

polygon(x = c(RBA_lower_mean_adjusted_model_D4[5], RBA_upper_mean_adjusted_model_D4[5], 
              RBA_upper_mean_adjusted_model_D4[5], RBA_lower_mean_adjusted_model_D4[5]), 
        c(7.5,7.5,6.5,6.5), col = 'aquamarine4', border = 'aquamarine4')

polygon(x = c(RBA_upper_mean_adjusted_model_D4[5], RBA_upper_ci_adjusted_model_D4[5], 
              RBA_upper_mean_adjusted_model_D4[5], RBA_upper_mean_adjusted_model_D4[5]), 
        c(7.5,7,6.5,7.5), col = 'blue', border = 'blue')

text(2.98, 7, labels = paste0('(',round(RBA_lower_mean_adjusted_model_D4[5], 2),'-',
                           round(RBA_upper_mean_adjusted_model_D4[5],2), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D4[5], 2),', ',
                           round(RBA_upper_ci_adjusted_model_D4[5],2),']'), col = "blue", pos=4, cex=1)


text(-0.95, 7, labels = paste0('4'), col = "black", pos = 4, cex = 1)

# Domain 5-6

polygon(x = c(RBA_lower_ci_adjusted_model_D5[5], RBA_lower_mean_adjusted_model_D5[5], 
              RBA_lower_mean_adjusted_model_D5[5], RBA_lower_ci_adjusted_model_D5[5]), 
        c(5,5.5,4.5,5), col = 'blue', border = 'blue')

polygon(x = c(RBA_lower_mean_adjusted_model_D5[5], RBA_upper_mean_adjusted_model_D5[5], 
              RBA_upper_mean_adjusted_model_D5[5], RBA_lower_mean_adjusted_model_D5[5]), 
        c(5.5,5.5,4.5,4.5), col = 'aquamarine4', border = 'aquamarine4')

polygon(x = c(RBA_upper_mean_adjusted_model_D5[5], RBA_upper_ci_adjusted_model_D5[5], 
              RBA_upper_mean_adjusted_model_D5[5], RBA_upper_mean_adjusted_model_D5[5]), 
        c(5.5,5,4.5,5.5), col = 'blue', border = 'blue')

text(2.98, 5, labels = paste0('(',round(RBA_lower_mean_adjusted_model_D5[5], 2),'-',
                           round(RBA_upper_mean_adjusted_model_D5[5],2), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_D5[5], 2),', ',
                           round(RBA_upper_ci_adjusted_model_D5[5],2),']'), col = "blue", pos=4, cex=1)


text(-0.95, 5, labels = paste0('5-6'), col = "black", pos = 4, cex = 1)


# Domain all
polygon(x = c(RBA_lower_ci_adjusted_model_all_domains[5], RBA_lower_mean_adjusted_model_all_domains[5], 
              RBA_lower_mean_adjusted_model_all_domains[5], RBA_lower_ci_adjusted_model_all_domains[5]), 
        c(3,3.5,2.5,3), col = 'blue', border = 'blue')

polygon(x = c(RBA_lower_mean_adjusted_model_all_domains[5], RBA_upper_mean_adjusted_model_all_domains[5], 
              RBA_upper_mean_adjusted_model_all_domains[5], RBA_lower_mean_adjusted_model_all_domains[5]), 
        c(3.5,3.5,2.5,2.5), col = 'aquamarine4', border = 'aquamarine4')

polygon(x = c(RBA_upper_mean_adjusted_model_all_domains[5], RBA_upper_ci_adjusted_model_all_domains[5], 
              RBA_upper_mean_adjusted_model_all_domains[5], RBA_upper_mean_adjusted_model_all_domains[5]), 
        c(3.5,3,2.5,3.5), col = 'blue', border = 'blue')

text(2.98, 3, labels = paste0('(',round(RBA_lower_mean_adjusted_model_all_domains[5], 2),'-',
                           round(RBA_upper_mean_adjusted_model_all_domains[5],2), ')', ' ',
                           '[',round(RBA_lower_ci_adjusted_model_all_domains[5], 2),'0',', ',
                           round(RBA_upper_ci_adjusted_model_all_domains[5],2),']'), col = "blue", pos=4, cex=1)

text(-0.95, 3, labels = paste0('All'), col = "black", pos = 4, cex = 1)

#####################################
## unadjusted model
polygon(x = c(estimated_lower_ci_unadjusted_model[5], estimated_mean_unadjusted_model[5], 
              estimated_upper_ci_unadjusted_model[5], estimated_mean_unadjusted_model[5],
              estimated_lower_ci_unadjusted_model[5]), 
        c(1,1.5,1,0.5,1), col = 'black', border = 'black')
text(-0.95, 1, labels = paste0('None'), col = "black", pos = 4, cex = 1)

text(3.26, 1, labels = paste0(round(estimated_mean_unadjusted_model[5], 2), ' ',
                                '[',round(estimated_lower_ci_unadjusted_model[5], 2),', ',
                                round(estimated_upper_ci_unadjusted_model[5],2),']'), 
     col = "black", pos=4, cex=1)
