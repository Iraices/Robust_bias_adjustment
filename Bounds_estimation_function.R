library('rjags')
library('runjags')
library('forestplot')
library('gtools')
library('Matrix')

## Data 
obs <- matrix(c(10,80,5,17,16,41,16,44), byrow = TRUE, nrow = 4, ncol = 2)
sample_size <- matrix(c(201,298,40,40,122,122,172,170), byrow = TRUE, nrow = 4, ncol = 2) 
parameters <- c("beta", "sigma2_theta", "mu", "OR", "delta", "p")

colnames(obs) = c('Control', 'Treatment')
colnames(sample_size) = c('Control', 'Treatment')
rownames(obs) = c('REFLEX', 'WA16291', 'DANCER', 'SERENE')
rownames(sample_size) = c('REFLEX', 'WA16291', 'DANCER', 'SERENE')


## Functions 
## q_values must be a dataframe with q (study quality) combinations and 4 columns q1, q2, q3 q4

## this function estimates bounds on expectation and percentile of the overall effect
estimate_bounds_finite_set <- function(n_studies, sample_size, obs, mu_phi, q_values, 
                                       mu_mu = 0, sigma2_mu = 10, mu_beta = 0, sigma2_beta = 10, 
                                       alpha = 0.01, lambda = 0.01, 
                                       n_chains = 2, n_iter = 20000, threshold = 1, percentile = 0.05){
  
  num_q_values <- dim(q_values)[1]
  
  output <- run.jags('bias_adjusted_model.txt', 
                     data <- list(n_studies = n_studies, N = sample_size, obs = obs,
                                  q = c(q_values[1,1], q_values[1,2], q_values[1,3], q_values[1,4]),
                                  mu_mu = mu_mu, sigma2_mu = sigma2_mu, mu_beta = mu_beta, sigma2_beta = sigma2_beta,
                                  alpha = alpha, lambda = lambda),
                     monitor =  c("beta", "sigma2_theta", "mu", "OR", "delta", "p"), 
                     n.chains = n_chains, burnin = n_iter/2, 
                     sample = n_iter, method = 'rjags')
  
  ## Expectation of the overall effect (mu)
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
  
  # Percentile of overall effect (mu)
  min_percentile <- quantile(as.mcmc(output, var = 'mu'), probs = percentile) 
  max_percentile <- min_percentile
  
  argmin_q_percentile <- c(q_values[1,1], q_values[1,2], q_values[1,3], q_values[1,4])
  argmax_q_percentile <- argmin_q_percentile
  
  out_min_percentile <- output
  out_max_percentile <- out_min_percentile
  
  # 2.5th Percentile of overall effect (mu)
  min_2_5th_percentile <- quantile(as.mcmc(output, var = 'mu'), probs = 0.025) 
  max_2_5th_percentile <- min_2_5th_percentile
  
  argmin_q_2_5th_percentile <- c(q_values[1,1], q_values[1,2], q_values[1,3], q_values[1,4])
  argmax_q_2_5th_percentile <- argmin_q_2_5th_percentile
  
  out_min_2_5th_percentile <- output
  out_max_2_5th_percentile <- out_min_2_5th_percentile
  
  # 97.5th Percentile of overall effect (mu)
  min_97_5th_percentile <- quantile(as.mcmc(output, var = 'mu'), probs = 0.975) 
  max_97_5th_percentile <- min_97_5th_percentile
  
  argmin_q_97_5th_percentile <- c(q_values[1,1], q_values[1,2], q_values[1,3], q_values[1,4])
  argmax_q_97_5th_percentile <- argmin_q_97_5th_percentile
  
  out_min_97_5th_percentile <- output
  out_max_97_5th_percentile <- out_min_97_5th_percentile
  
  for(i in 2:num_q_values){
    print(i)
    out <- run.jags('bias_adjusted_model.txt', 
                    data <- list(n_studies = n_studies, N = sample_size, obs = obs,
                                 q = c(q_values[i,1], q_values[i,2], q_values[i,3], q_values[i,4]),
                                 mu_mu = mu_mu, sigma2_mu = sigma2_mu, mu_beta = mu_beta, sigma2_beta = sigma2_beta,
                                 alpha = alpha, lambda = lambda),
                    monitor =  c("beta", "sigma2_theta", "mu", "OR", "delta", "p"), 
                    n.chains = n_chains, burnin = n_iter/2, 
                    sample = n_iter, method = 'rjags')
    
    ## Expectation overall effect (mu)
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
    
    ## Percentile of overall effect (mu)
    temp_percentile = quantile(as.mcmc(out, var = 'mu'), probs = percentile)
    
    if(temp_percentile < min_percentile){
      min_percentile <- temp_percentile
      argmin_q_percentile <- c(q_values[i,1], q_values[i,2], q_values[i,3], q_values[i,4])
      out_min_percentile <- out
    }
    if(temp_percentile > max_percentile){
      max_percentile <- temp_percentile
      argmax_q_percentile <- c(q_values[i,1], q_values[i,2], q_values[i,3], q_values[i,4])
      out_max_percentile <- out
    }
    ## 2.5 th Percentile of overall effect (mu)
    temp_2_5th_percentile = quantile(as.mcmc(out, var = 'mu'), probs = 0.025)
    
    if(temp_2_5th_percentile < min_2_5th_percentile){
      min_2_5th_percentile <- temp_2_5th_percentile
      argmin_q_2_5th_percentile <- c(q_values[i,1], q_values[i,2], q_values[i,3], q_values[i,4])
      out_min_2_5th_percentile <- out
    }
    if(temp_2_5th_percentile > max_2_5th_percentile){
      max_2_5th_percentile <- temp_2_5th_percentile
      argmax_q_2_5th_percentile <- c(q_values[i,1], q_values[i,2], q_values[i,3], q_values[i,4])
      out_max_2_5th_percentile <- out
    }
    ## 97.5 th Percentile of overall effect (mu)
    temp_97_5th_percentile = quantile(as.mcmc(out, var = 'mu'), probs = 0.975)
    
    if(temp_97_5th_percentile < min_97_5th_percentile){
      min_97_5th_percentile <- temp_97_5th_percentile
      argmin_q_97_5th_percentile <- c(q_values[i,1], q_values[i,2], q_values[i,3], q_values[i,4])
      out_min_97_5th_percentile <- out
    }
    if(temp_97_5th_percentile > max_97_5th_percentile){
      max_97_5th_percentile <- temp_97_5th_percentile
      argmax_q_97_5th_percentile <- c(q_values[i,1], q_values[i,2], q_values[i,3], q_values[i,4])
      out_max_97_5th_percentile <- out
    }
    
  }
  
  return(list(bounds = list(expectation = c(minimum = min_mu, maximum = max_mu),
                            exceeding_prob = c(min_exceeding_prob, max_exceeding_prob),
                            percentile = c(min_percentile, max_percentile),
                            percentile_2_5th = c(min_2_5th_percentile, max_2_5th_percentile),
                            percentile_97_5th = c(min_97_5th_percentile, max_97_5th_percentile)),
              q = list (expectation = c(argmin_q_expectation = argmin_q_expectation, 
                                        argmax_q_expectation = argmax_q_expectation),
                        exceeding_prob = c(argmin_q_exceeding_prob = argmin_q_exceeding_prob, 
                                           argmax_q_exceeding_prob = argmax_q_exceeding_prob),
                        percentile = c(argmin_q_percentile = argmin_q_percentile, 
                                            argmax_q_percentile = argmax_q_percentile),
                        percentile_2_5th = c(argmin_q_2_5th_percentile = argmin_q_2_5th_percentile,
                                             argmax_q_2_5th_percentile = argmax_q_2_5th_percentile),
                        percentile_97_5th = c(argmin_q_97_5th_percentile = argmin_q_97_5th_percentile,
                                              argmax_q_97_5th_percentile = argmax_q_97_5th_percentile)),
              output_minimum_expectation = out_min_mu,
              output_maximum_expectation = out_max_mu, 
              output_minimum_exceeding_prob = out_min_exceeding_prob,
              output_maximum_exceeding_prob = out_max_exceeding_prob,
              output_minimum_percentile = out_min_percentile,
              output_maximum_percentile = out_max_percentile,
              output_minimum_2_5th_percentile = out_min_2_5th_percentile,
              output_maximum_2_5th_percentile = out_max_2_5th_percentile,
              output_minimum_97_5th_percentile = out_min_97_5th_percentile,
              output_maximum_97_5th_percentile = out_max_97_5th_percentile))
}


#########################################################################################################
#########################################################################################################

## Domain 5 and 6
low_val_domain_5 = 0.5
high_val_domain_5 = 0.95
q_values_domain_5 = matrix(rep(seq(low_val_domain_5, high_val_domain_5, 0.05), 4), ncol = 4)
colnames(q_values_domain_5) = c('Q1','Q2', 'Q3','Q4')
q_values_domain_5 = as.data.frame(q_values_domain_5)

results_domain_5 <- estimate_bounds_finite_set(n_studies = 4, sample_size = sample_size, obs = obs, 
                                              q_values = q_values_domain_5, 
                                              mu_mu = 0, sigma2_mu = 10, mu_beta = 0, sigma2_beta = 10, 
                                              alpha = 0.01, lambda = 0.01, 
                                              n_chains = 2, n_iter = 20000, threshold = 1, percentile = 0.05)

# save(results_domain_5, file = 'results_domain_5.RData')

##############################################################
##############################################################

## Domain 1 and 2
low_val_domain_1 = 0.1
high_val_domain_1 = 0.95
q_values_domain_1 = matrix(rep(c(0,0,0,0), 10000), ncol = 4)


vals = seq(0.1,0.95, length.out = 10)
len_vals = length(vals) #10
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

colnames(q_values_domain_1) = c('Q1','Q2', 'Q3','Q4')
q_values_domain_1 = as.data.frame(q_values_domain_1)

results_domain_1 <- estimate_bounds_finite_set(n_studies = 4, sample_size = sample_size, obs = obs, 
                                              q_values = q_values_domain_1, 
                                              mu_mu = 0, sigma2_mu = 10, mu_beta = 0, sigma2_beta = 10, 
                                              alpha = 0.01, lambda = 0.01, 
                                              n_chains = 2, n_iter = 20000, threshold = 1, percentile = 0.05)

# save(results_domain_1, file = 'results_domain_1.RData')

#############################################################################
#############################################################################

## Domain 4
## The region
## q_1 <= q_2 
## q_1 >= 0.1
## 0.5 <= q_2 <= 0.95 
# The vertices of the region are 
p_1 = c(0.1, 0.95)
p_2 = c(0.1, 0.5)
p_3 = c(0.5, 0.5)
p_4 = c(0.95, 0.95)
## We are interested in the points inside the region

points_inside_domain4 <- c()
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
      alpha_4 <- 1 - alpha_1[i] - alpha_2[j] - alpha_3[k]
      
      p <- alpha_1[i] * p_1 + alpha_2[j] * p_2 + alpha_3[k] * p_3 + alpha_4 * p_4
      if(p[1] >= 0.1 & p[1] <= p[2] & p[2] >= 0.5 & p[2] <= 0.95){
      
        points_inside_domain4 <- rbind(points_inside_domain4, p)
      }
    }
  }
}  
points_inside_domain4


points_inside_domain4 <- as.data.frame(points_inside_domain4, col.rows = NULL)
## Add q_3 and q_4
points_inside_domain4$V3 = points_inside_domain4$V2
points_inside_domain4$V4 = points_inside_domain4$V2
q_values_domain_4 <- unique(points_inside_domain4) 
q_values_domain_4
colnames(q_values_domain_4) <- c('Q1','Q2','Q3','Q4')

## 286 points

results_domain_4 <- estimate_bounds_finite_set(n_studies = 4, sample_size = sample_size, obs = obs, 
                                               q_values = q_values_domain_4, 
                                               mu_mu = 0, sigma2_mu = 10, mu_beta = 0, sigma2_beta = 10, 
                                               alpha = 0.01, lambda = 0.01, 
                                               n_chains = 2, n_iter = 20000, threshold = 1, percentile = 0.05)

# save(results_domain_4, file = 'results_domain_4.RData')

#############################################################################
#############################################################################
## Domain 3
## The region
## q_3 <= q_1 
## q_4 <= q_1 
## q_3 >= 0.1
## q_4 >= 0.1
## 0.5 <= q_1 <= 0.95 

## Ax = b
A <- matrix(c(-1,1,0,1,0,0,0,0,0,
              -1,0,1,0,1,0,0,0,0,
              1,0,0,0,0,1,0,0,0,
              -1,0,0,0,0,0,1,0,0,
              0,-1,0, 0,0,0,0,1,0,
              0,0,-1,0,0,0,0,0,1), nrow = 6, ncol = 9, byrow = TRUE)

b <- c(0, 0, 0.95,-0.5, -0.1, -0.1)

# combinations of system equations with 6 variables
extended_number_variables <- 9
number_of_equations <- 6

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

nrows <- dim(solutions)[1]
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
solutions_unique <- unique(solutions[-c(n_NaN),])

solutions_unique_conditions <- c()

for(i in 1:dim(solutions_unique)[1]){
  q_1 = solutions_unique[i,][1]
  q_3 = solutions_unique[i,][2]
  q_4 = solutions_unique[i,][3]
  if(q_1 >= 0.5 & q_1 <= 0.95 & q_3 >= 0.1 & q_3 <= q_1 & q_4 >= 0.1 & q_4 <= q_1){
    solutions_unique_conditions <- rbind(solutions_unique_conditions, solutions_unique[i,])
  }
}
solutions_unique_conditions[,c(1,2,3)]


# The vertices of the polyhedron are 
## (q_1, q_3, q_4)
v_1 = c(0.5, 0.1, 0.1)
v_2 = c(0.95, 0.1, 0.1)
v_3 = c(0.5, 0.1, 0.5)
v_4 = c(0.95, 0.1, 0.95)
v_5 = c(0.5, 0.5, 0.1)
v_6 = c(0.95, 0.95, 0.1)
v_7 = c(0.5, 0.5, 0.5)
v_8 = c(0.95, 0.95, 0.95)
## We are interested in the points inside the region

points_inside_domain3 <- c()
## points in the convex set
## \sum_{i = 1}^{nrows} alpha_i * pi  with \sum_{i = 1}^{nrows} alpha_i = 1
##

alpha_1 <- seq(0,1, by = 0.2)
L1 <- length(alpha_1)

for(i in 1:L1){
  alpha_2 <- seq(0,1 - alpha_1[i], by = 0.2)
  L2 <- length(alpha_2)
  for(j in 1:L2){
    alpha_3 <- seq(0,1 - alpha_1[i] - alpha_2[j], by = 0.2)
    L3 <- length(alpha_3)
    for(k in 1:L3){
      alpha_4 <- seq(0,1 - alpha_1[i] - alpha_2[j] - alpha_3[k], by = 0.2)
      L4 <- length(alpha_4)
      for(m in 1:L4){
        alpha_5 <- seq(0,1 - alpha_1[i] - alpha_2[j] - alpha_3[k] - alpha_4[m], by = 0.2)
        L5 <- length(alpha_5)
        for(l in 1:L5){
          alpha_6 <- seq(0,1 - alpha_1[i] - alpha_2[j] - alpha_3[k] - alpha_4[m] - alpha_5[l], by = 0.2)
          L6 <- length(alpha_6)
          for(n in 1:L6){
            alpha_7 <- seq(0,1 - alpha_1[i] - alpha_2[j] - alpha_3[k] - alpha_4[m] - alpha_5[l] - alpha_6[n], by = 0.2)
            L7 <- length(alpha_7)  
            for(v in 1:L7){
              alpha_8 <- 1 - alpha_1[i] - alpha_2[j] - alpha_3[k] - alpha_4[m] - alpha_5[l] - alpha_6[n] - alpha_7[v]
              
              p <- alpha_1[i]* v_1 + alpha_2[j]* v_2 + alpha_3[k]* v_3 + alpha_4[m]* v_4 +
                    alpha_5[l]* v_5 + alpha_6[n]* v_6 + alpha_7[v]* v_7 + alpha_8 * v_8
              
              if(p[1] >= 0.5 & p[1] <= 0.95 & p[2] <= p[1] & p[3] <= p[1] & p[2] >= 0.1 & p[3] >= 0.1){
                points_inside_domain3 <- rbind(points_inside_domain3, p)
              }
            }
          }
        }
      }
    }
  }
}

points_inside_domain3 <- as.data.frame(points_inside_domain3, col.rows = NULL)
colnames(points_inside_domain3) <- c('Q1','Q3','Q4')
points_inside_domain3$Q2 = points_inside_domain3$Q1 
points_inside_domain3 = points_inside_domain3[,c('Q1','Q2','Q3','Q4')]

## Remove duplicate points
q_values_domain_3 <- unique(points_inside_domain3) 

## 736 points

results_domain_3 <- estimate_bounds_finite_set(n_studies = 4, sample_size = sample_size, obs = obs, 
                                               q_values = q_values_domain_3, 
                                               mu_mu = 0, sigma2_mu = 10, mu_beta = 0, sigma2_beta = 10, 
                                               alpha = 0.01, lambda = 0.01, 
                                               n_chains = 2, n_iter = 20000, threshold = 1, percentile = 0.05)

# save(results_domain_3, file = 'results_domain_3.RData')


#############################################################################
#############################################################################
## All Domains
## The region
## q_1 <= q_2 
## q_3 <= q_2 
## q_3 = q_4
## q_1 >= 0.1
## q_3 >= 0.1
## 0.1 <= q_2 <= 0.95


## Ax = b
A <- matrix(c( 1, -1,  0, 1, 0, 0, 0, 0, 0,
              -1,  0,  0, 0, 1, 0, 0, 0, 0,
               0, -1,  1, 0, 0, 1, 0, 0, 0,
               0,  0, -1, 0, 0, 0, 1, 0, 0,
               0,  1,  0, 0, 0, 0, 0, 1, 0,
               0, -1,  0, 0, 0, 0, 0, 0, 1), nrow = 6, ncol = 9, byrow = TRUE)

b <- c(0, -0.1, 0, -0.1, 0.95, -0.1)

# combinations of system equations with 6 variables
extended_number_variables <- 9
number_of_equations <- 6

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

nrows <- dim(solutions)[1]
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
solutions_unique <- unique(solutions[-c(n_NaN),])

solutions_unique_conditions <- c()

for(i in 1:dim(solutions_unique)[1]){
  q_1 = solutions_unique[i,][1]
  q_2 = solutions_unique[i,][2]
  q_3 = solutions_unique[i,][3]
  if(q_1 >= 0.1 & q_1 <= q_2 & q_3 >= 0.1 & q_3 <= q_2 & q_2 >= 0.1 & q_2 <= 0.95){
    solutions_unique_conditions <- rbind(solutions_unique_conditions, solutions_unique[i,])
  }
}
solutions_unique_conditions[,c(1,2,3)]


# The vertices of the polyhedron are 
## (q_1, q_2, q_3)
z_1 = c(0.1, 0.1, 0.1)
z_2 = c(0.1, 0.95, 0.1)
z_3 = c(0.1, 0.95, 0.95)
z_4 = c(0.95, 0.95, 0.1)
z_5 = c(0.95, 0.95, 0.95)

## We are interested in the points inside the region

points_inside_all_domains <- c()
## points in the convex set
## \sum_{i = 1}^{nrows} alpha_i * zi  with \sum_{i = 1}^{nrows} alpha_i = 1
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
        alpha_5 <- 1 - alpha_1[i] - alpha_2[j] - alpha_3[k] - alpha_4[m] 
        
        p <- alpha_1[i]* z_1 + alpha_2[j]* z_2 + alpha_3[k]* z_3 + alpha_4[m]* z_4 + alpha_5 * z_5
              
        if(p[1] >= 0.1 & p[1] <= p[2] & p[3] >= 0.1 & p[3] <= p[2] & p[2] >= 0.1 & p[2] <= 0.95){
          points_inside_all_domains <- rbind(points_inside_all_domains, p)

        }
      }
    }
  }
}

points_inside_all_domains <- as.data.frame(points_inside_all_domains, col.rows = NULL)
colnames(points_inside_all_domains) <- c('Q1','Q2','Q3')
points_inside_all_domains$Q4 = points_inside_all_domains$Q3 
points_inside_all_domains = points_inside_all_domains[,c('Q1','Q2','Q3','Q4')]

## Remove duplicate points
q_values_all_domains <- unique(points_inside_all_domains) 

## 839 points

results_all_domains <- estimate_bounds_finite_set(n_studies = 4, sample_size = sample_size, obs = obs, 
                                               q_values = q_values_all_domains, 
                                               mu_mu = 0, sigma2_mu = 10, mu_beta = 0, sigma2_beta = 10, 
                                               alpha = 0.01, lambda = 0.01, 
                                               n_chains = 2, n_iter = 20000, threshold = 1, percentile = 0.05)

# save(results_all_domains, file = 'results_all_domains.RData')


