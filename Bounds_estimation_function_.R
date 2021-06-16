library('rjags')
library('runjags')
library('forestplot')
library('gtools')
library('Matrix')

## Data 
## n_studies <- 4  ## number of studies
obs <- matrix(c(10,80,5,17,16,41,16,44), byrow = TRUE, nrow = 4, ncol = 2)
sample_size <- matrix(c(201,298,40,40,122,122,172,170), byrow = TRUE, nrow = 4, ncol = 2) 
#parameters <- c("beta", "theta", "sigma2_theta", "mu", "OR", "delta", "p")
parameters <- c("beta", "sigma2_theta", "mu", "OR", "delta", "p")

colnames(obs) = c('Control', 'Treatment')
colnames(sample_size) = c('Control', 'Treatment')
rownames(obs) = c('REFLEX', 'WA16291', 'DANCER', 'SERENE')
rownames(sample_size) = c('REFLEX', 'WA16291', 'DANCER', 'SERENE')


## Functions 
## q_values must be a dataframe with q (study quality) combinations and 4 columns q1, q2, q3 q4

## this function estimates bounds on expectations
estimate_bounds_finite_set <- function(n_studies, sample_size, obs, mu_phi, q_values, 
                                       mu_mu = 0, sigma2_mu = 10, mu_beta = 0, sigma2_beta = 10, 
                                       alpha = 0.01, lambda = 0.01, 
                                       n_chains = 2, n_iter = 20000, threshold = 1, percentile = 0.05){
  
  num_q_values <- dim(q_values)[1]
  
  output <- run.jags('bias_adjusted_model_150121.txt', 
                     data <- list(n_studies = n_studies, N = sample_size, obs = obs,
                                  q = c(q_values[1,1], q_values[1,2], q_values[1,3], q_values[1,4]),
                                  mu_mu = mu_mu, sigma2_mu = sigma2_mu, mu_beta = mu_beta, sigma2_beta = sigma2_beta,
                                  alpha = alpha, lambda = lambda),
                     monitor =  c("beta", "sigma2_theta", "mu", "OR", "delta", "p"), 
                     n.chains = n_chains, burnin = n_iter/2, 
                     sample = n_iter, method = 'rjags')
  
  ## Expectation
  min_mu <- output$summary[1]$statistics['mu', 'Mean'] 
  max_mu <- min_mu
  
  argmin_q_expectation <- c(q_values[1,1], q_values[1,2], q_values[1,3], q_values[1,4])
  argmax_q_expectation <- argmin_q_expectation
  
  out_min_mu <- output
  out_max_mu <- out_min_mu
  
  ## Exceedance probability
  min_exceeding_prob <- mean(as.mcmc(output, var = 'mu') >= threshold) 
  max_exceeding_prob <- min_exceeding_prob
  
  argmin_q_exceeding_prob <- c(q_values[1,1], q_values[1,2], q_values[1,3], q_values[1,4])
  argmax_q_exceeding_prob <- argmin_q_exceeding_prob
  
  out_min_exceeding_prob <- output
  out_max_exceeding_prob <- out_min_exceeding_prob
  
  # Percentile of expectation
  min_percentile <- quantile(as.mcmc(output, var = 'mu'), probs = percentile) 
  max_percentile <- min_percentile
  
  argmin_q_percentile <- c(q_values[1,1], q_values[1,2], q_values[1,3], q_values[1,4])
  argmax_q_percentile <- argmin_q_percentile
  
  out_min_percentile <- output
  out_max_percentile <- out_min_percentile
  
  for(i in 2:num_q_values){
    print(i)
    out <- run.jags('bias_adjusted_model_150121.txt', 
                    data <- list(n_studies = n_studies, N = sample_size, obs = obs,
                                 q = c(q_values[i,1], q_values[i,2], q_values[i,3], q_values[i,4]),
                                 mu_mu = mu_mu, sigma2_mu = sigma2_mu, mu_beta = mu_beta, sigma2_beta = sigma2_beta,
                                 alpha = alpha, lambda = lambda),
                    monitor =  c("beta", "sigma2_theta", "mu", "OR", "delta", "p"), 
                    n.chains = n_chains, burnin = n_iter/2, 
                    sample = n_iter, method = 'rjags')
    
    ## Expectation
    if(out$summary[1]$statistics['mu', 'Mean'] < min_mu){
      min_mu <- out$summary[1]$statistics['mu', 'Mean']
      argmin_q_expectation <- c(q_values[i,1], q_values[i,2], q_values[i,3], q_values[i,4])
      out_min_mu <- out
    }
    if(out$summary[1]$statistics['mu', 'Mean'] > max_mu){
      max_mu <- out$summary[1]$statistics['mu', 'Mean']
      argmax_q_expectation <- c(q_values[i,1], q_values[i,2], q_values[i,3], q_values[i,4])
      out_max_mu <- out
    }
    
    ## Exceedance probability
    temp_prob = mean(as.mcmc(out, var = 'mu') >= threshold)
    
    if(temp_prob < min_exceeding_prob){
      min_exceeding_prob <- temp_prob
      argmin_q_exceeding_prob <- c(q_values[i,1], q_values[i,2], q_values[i,3], q_values[i,4])
      out_min_exceeding_prob <- out
    }
    if(temp_prob > max_exceeding_prob){
      max_exceeding_prob <- temp_prob
      argmax_q_exceeding_prob <- c(q_values[i,1], q_values[i,2], q_values[i,3], q_values[i,4])
      out_max_exceeding_prob <- out
    }
    
    ## Percentile of expectation
    temp_percentile_expectation = quantile(as.mcmc(out, var = 'mu'), probs = percentile)
    
    if(temp_percentile_expectation < min_percentile){
      min_percentile <- temp_percentile_expectation
      argmin_q_percentile <- c(q_values[i,1], q_values[i,2], q_values[i,3], q_values[i,4])
      out_min_percentile <- out
    }
    if(temp_percentile_expectation > max_percentile){
      max_percentile <- temp_percentile_expectation
      argmax_q_percentile <- c(q_values[i,1], q_values[i,2], q_values[i,3], q_values[i,4])
      out_max_percentile <- out
    }
    
  }
  
  return(list(bounds = list(expectation = c(minimum = min_mu, maximum = max_mu),
                            exceeding_prob = c(min_exceeding_prob, max_exceeding_prob),
                            percentile = c(min_percentile, max_percentile)),
              q = list (expectation = c(argmin_q_expectation = argmin_q_expectation, 
                                        argmax_q_expectation = argmax_q_expectation),
                        exceeding_prob = c(argmin_q_exceeding_prob = argmin_q_exceeding_prob, 
                                           argmax_q_exceeding_prob = argmax_q_exceeding_prob),
                        percentile = c(argmin_q_percentile = argmin_q_percentile, 
                                            argmax_q_percentile = argmax_q_percentile)),
              output_minimum_expectation = out_min_mu,
              output_maximum_expectation = out_max_mu, 
              output_minimum_exceeding_prob = out_min_exceeding_prob,
              output_maximum_exceeding_prob = out_max_exceeding_prob,
              out_minimum_percentile = out_min_percentile,
              out_maximum_percentile = out_max_percentile))
}


#########################################################################################################
#########################################################################################################

## Domain 5 and 6
## low_val_domain_5 = 0.6
## high_val_domain_5 = 1
## q_values_domain_5 = matrix(rep(seq(low_val_domain_5, high_val_domain_5, 0.1), 4), ncol = 4)
## colnames(q_values_domain_5) = c('Q1','Q2', 'Q3','Q4')
## q_values_domain_5 = as.data.frame(q_values_domain_5)

results_domain_5 <- estimate_bounds_finite_set(n_studies = 4, sample_size = sample_size, obs = obs, 
                                              q_values = q_values_domain_5, 
                                              mu_mu = 0, sigma2_mu = 10, mu_beta = 0, sigma2_beta = 10, 
                                              alpha = 0.01, lambda = 0.01, 
                                              n_chains = 2, n_iter = 20000, threshold = 1, percentile = 0.05)

save(results_domain_5, file = 'results_domain_5.RData')

##############################################################
##############################################################

## Domain 1 and 2
## low_val_domain_1 = 0.1
## high_val_domain_1 = 1
## space = length(seq(0.1,1,0.1))
## q_values_domain_1 = matrix(rep(c(0,0,0,0), 10000), ncol = 4)


vals = seq(0.1,1,0.1)
len_vals = length(vals)
n = 1
for(i in 1:len_vals){
  for(j in 1:len_vals){
    for(k in 1:len_vals){
      for(l in 1:len_vals){
        q_values_domain_1[n,] = c(vals[i], vals[j], vals[k], vals[l])
        n = n + 1
      }  
    }
  }
}

## colnames(q_values_domain_1) = c('Q1','Q2', 'Q3','Q4')
## q_values_domain_1 = as.data.frame(q_values_domain_1)

results_domain_1 <- estimate_bounds_finite_set(n_studies = 4, sample_size = sample_size, obs = obs, 
                                              q_values = q_values_domain_1, 
                                              mu_mu = 0, sigma2_mu = 10, mu_beta = 0, sigma2_beta = 10, 
                                              alpha = 0.01, lambda = 0.01, 
                                              n_chains = 2, n_iter = 20000, threshold = 1, percentile = 0.05)

save(results_domain_1, file = 'results_domain_1.RData')

#############################################################################
#############################################################################

## Domain 4
## the inequality constraint region Ax < b with x = (x1,x2) is transformed into 
## Ax = b with x = (x1,x2,x3,x4,x5,x6)  
A <- matrix(c(1,-1,1,0,0,0,
              0, 1,0,1,0,0,
              0,-1,0,0,1,0,
             -1, 0,0,0,0,1), nrow = 4, ncol = 6, byrow = TRUE)
b <- c(0, 1, -0.6, -0.1)

# combinations of system equations with 3 variables
extended_number_variables <- 6
number_of_equations <- 4

comb <- combinations(extended_number_variables, number_of_equations, repeats.allowed = FALSE)
# number of combinations
n_comb <- dim(comb)[1]

## Find the points that define the convex set. 
solutions <- matrix(0, nrow = n_comb, ncol = extended_number_variables)

for(i in 1:n_comb){
  rank_A <- rankMatrix(A[, c(comb[i,])])[1]
  rank_Ab <- rankMatrix(cbind(A[, c(comb[i,])],b))[1]
  nrows <- dim(A)[1]
  if(rank_A == rank_Ab){
    if(rank_A == nrows){
      solutions[i, c(comb[i,])] <- solve(A[, c(comb[i,])], b)
    }
    if(rank_A < nrows){
      solutions[i, c(comb[i,])] <- 'many'
    }
  }
  else{
    solutions[i,] <- NaN
  }
}

solutions


## Now, we use the variables we are interested and remove the NaN solutions
nrows <- dim(solutions[,1:2])[1]
n_NaN<- rep(0,nrows)

for(i in 1:nrows){
  if(is.nan(as.numeric(solutions[i,1]))){
    n_NaN[i] <- i 
  }
  else{
    n_NaN[i] <- 0
  }
}

## These are the vertices of the polyhedron. 
## We are interested in the values inside the polyhedron
solutions_unique <- unique(solutions[-c(n_NaN),1:2])
solutions_unique

## selecting the vertices that satisfies the conditions
solutions_unique_conditions <- c()

#### Revisar

for(i in 1:dim(solutions_unique)[1]){
  q_1 = solutions_unique[i,][1]
  q_2 = solutions_unique[i,][2]
  if(q_1 >= 0.1 & q_1 <= q_2 & q_2 >= 0.6  & q_2 <= 1){
    solutions_unique_conditions <- rbind(solutions_unique_conditions, solutions_unique[i,])
  }
}
solutions_unique_conditions

## number of vertices
num_alphas <- dim(solutions_unique_conditions)[1]

all_points <- c()
## points in the convex set
## \sum_{i = 1}^{nrows} alpha_i * pi  with \sum_{i = 1}^{nrows} alpha_i = 1
##

alpha_1 <- seq(0,1, by = 0.1)
L1 <- length(alpha_1)

for(i in 1:L1){
  alpha_2 <- seq(0,1 - alpha_1[i], by = 0.1)
  L2 <- length(alpha_2)
  for(j in 1:L2){
    alpha_3 <- seq(0,1 - alpha_1[i] - alpha_2[j], by = 0.1)
    L3 <- length(alpha_3)
    for(k in 1:L3){
      alpha_4 <- seq(0,1 - alpha_1[i] - alpha_2[j] - alpha_3[k], by = 0.1)
      L4 <- length(alpha_4)
      for(m in 1:L4){
        alpha_5 <- seq(0,1 - alpha_1[i] - alpha_2[j] - alpha_3[k] - alpha_4[m], by = 0.1)
        L5 <- length(alpha_5)
        for(l in 1:L5){
          alpha_6 <- seq(0,1 - alpha_1[i] - alpha_2[j] - alpha_3[k] - alpha_4[m] - alpha_5[l], by = 0.1)
          L6 <- length(alpha_6)
          for(n in 1:L6){
            alpha_7 <- seq(0,1 - alpha_1[i] - alpha_2[j] - alpha_3[k] - alpha_4[m] - alpha_5[l] - alpha_6[n], by = 0.1)
            L7 <- length(alpha_7)  
            for(v in 1:L7){
              alpha_8 <- 1 - alpha_1[i] - alpha_2[j] - alpha_3[k] - alpha_4[m] - alpha_5[l] - alpha_6[n] - alpha_7[v]
              
              p <- alpha_1[i]* as.numeric(solutions_unique[1,]) + alpha_2[j]* as.numeric(solutions_unique[2,]) +
                alpha_3[k]* as.numeric(solutions_unique[3,]) + alpha_4[m]* as.numeric(solutions_unique[4,]) +
                alpha_5[l]* as.numeric(solutions_unique[5,]) + alpha_6[n]* as.numeric(solutions_unique[6,]) +
                alpha_7[v]* as.numeric(solutions_unique[7,]) + alpha_8 * as.numeric(solutions_unique[8,])
              
              if(p[1] >= 0.6 & p[2] > 0 & p[3] > 0){
                all_points <- rbind(all_points, p)
              }
            }
          }
        }
      }
    }
  }
}

all_points <- as.data.frame(all_points, col.rows = NULL)
colnames(all_points) <- c('Q1','Q3','Q4')

## Remove duplicate points
my_all_points <- unique(all_points) 
