# Loading packages
library(readr)  
library(cmdstanr)
library(posterior)
library(bayesplot)
library(tidyverse)
library(MCMCpack)
library(glue)
library(loo)
library(fabricatr)
library(matrixStats)

#Loading the data
data = read_tsv("../../data/PCAWG_SV_with_only_clust_and_11_annotations_noNA.tsv")

### retrieving arguments from command line ###

# Getting arguments from command line
args <- commandArgs(trailingOnly = TRUE)

# Getting selected variables (i.e features)
var_args <- strsplit(args[1], ",")[[1]]
var_select <- vector(length = length(var_args))

for (i in 1:length(var_select)){
  var_select[i] = as.numeric(var_args[i])
}

# Getting number of signatures
K = as.numeric(args[2])

# Getting number of iterations
iterations = as.numeric(args[3])

# Getting refresh rate
refresh_rate = as.numeric(args[4])

# Getting folder with wanted results
folder <- args[5]

# Getting number of bins
n_bins <- as.numeric(args[6])

# Getting number of features from selected variables
num_f = length(var_select) - 1

### Preprocessing the data ###

#Removing patients with less than 20 observations 
data_count <- data %>% 
  group_by(sample) %>% 
  mutate(count = n())

data = data[-which(data_count$count < 20),]

#Transforming sample ID into unique integers
data$sample = as.numeric(as.factor(data.matrix(data[,1])))

#Transforming sv_class to discrete variable with following mapping (DEL: 1, DUP: 2, INV: 3)
data$sv_class = as.numeric(as.factor(data$sv_class))

#Transforming length into 5 bins based on fixed breaks
data$length = cut(data$length, breaks = c(0, 5000, 50000 , 500000, 5000000, max(data$length)), labels = c(1,2,3,4,5))

#Transforming clustering to (non-cluster: 1, cluster: 2)
data$clustering = data$clustering + 1

#Transforming the following features into a n_bins factor variable based on quantiles.
data$telomeredist = split_quantile(data$telomeredist, type = n_bins)
data$GC_content = split_quantile(data$GC_content, type = n_bins)
data$replication_time = split_quantile(data$replication_time, type = n_bins)
data$NCO = split_quantile(data$NCO, type = n_bins)
data$DSB = split_quantile(data$DSB, type = n_bins)
data$deltaDSB = split_quantile(data$deltaDSB, type = n_bins)

#Transofrming the last 5 features into 3 bins based on following mapping (0: 1, 0.5:2, 1:3)
data$fragile_sites = as.numeric(as.factor(data$fragile_sites))
data$LINE = as.numeric(as.factor(data$LINE))
data$SINE = as.numeric(as.factor(data$SINE))
data$Simple_repeat = as.numeric(as.factor(data$Simple_repeat))
data$LTR = as.numeric(as.factor(data$LTR))

#Selecting only the variables specified in the arguments
data = data[, var_select] # Here you can filter rows

# Setting seed for reproducibility
set.seed(123) 

# Function to split data into train and test per patient
split_data <- function(data, train_frac = 0.8) {
  train_data <- data %>%
    group_by(sample) %>%
    group_map(~ {
      n <- nrow(.x)
      train_indices <- sample(n, size = floor(train_frac * n))
      list(train = .x[train_indices, ], test = .x[-train_indices, ])
    }, .keep = TRUE)
  
  # Combine results
  train_set <- bind_rows(lapply(train_data, function(x) x$train))
  test_set  <- bind_rows(lapply(train_data, function(x) x$test))

  if (nrow(train_set) + nrow(test_set) != nrow(data)){
    print("Data splitting malfunctioned")
    break
  }
  
  list(train = train_set, test = test_set)
}

# Split the data
splits <- split_data(data)

train_data <- splits$train
test_data <- splits$test

### Defining variables from the Stan script ###

# Initializing the Stan model based on number of features
path_to_stan_model <- glue("{num_f}F/MFLDA_{num_f}F.stan")
file <- file.path(path_to_stan_model)
mod <- cmdstan_model(file)

# Size parameters
M = nrow(unique(train_data[,1]))      # Number of patients
N = nrow(train_data)                  # Total number of mutations
print(glue("Number of obs in N_train: {N}"))

# Finding number of unique types in each feature
T_c <- vector(length = num_f)
for(i in 2:length(var_select)){
  T_c[i-1] = nrow(unique(data[,i]))
}

# Defining the observed values from the data
obs <- vector(mode = 'list', length = length(var_select))
for(i in 1:length(var_select)){
  obs[[i]] <- as.integer(unlist(train_data[,i]))
}

# Defining priors for the model 
priors <- vector(mode = 'list', length = (num_f + 1))
priors[[1]] <- c(rep(1,K)) # Alpha: prior for theta
for(i in 1:num_f){
  priors[[i+1]] <- c(rep(1,T_c[[i]][1])) # Beta_i: prior for phi_i distribution over types
}

#Creating function to initialize values
init_fun <- function(T_c, K, M, num_f){ 
  init_list = list(theta = rdirichlet(M, rep(1, K)))
  for(i in 1:num_f){
    init_name = paste0("phi_c", glue({i}))
    init_list[[init_name]] = rdirichlet(K, rep(1, T_c[i]))
  }
  return(list(init_list)) 
}

### Creating the MCMC runs ###

# Adding the constant variables to data list
data_list <- list(K = K, C = num_f, M = M, N = N, T_c = T_c, 
                  pat = obs[[1]], alpha = priors[[1]])

# Adding the feature depended variables to data list
for(i in 1:num_f){
  feature_name = paste0("obs_c", glue({i}))
  data_list[[feature_name]] = obs[[i+1]]
  
  prior_name = paste0("beta_c", glue({i}))
  data_list[[prior_name]] = priors[[i+1]]
}

# Sampling from the model
print(glue("Running model with nsignature: {K}"))
time_before <- Sys.time()
fit <- mod$sample(
  data = data_list,
  init = init_fun(T_c, K, M, num_f),
  seed = 476,
  chains = 1,
  iter_sampling = iterations,
  parallel_chains = 1,
  refresh = refresh_rate
)
time_after <- Sys.time()
used_time <- time_after - time_before
print(glue("Finished model with nsignature: {K}"))

# Defining the summary from the sample
summ <- fit$summary()

# Creating loo for the sample
loglik = fit$draws("log_lik")
r_eff = relative_eff(exp(loglik))
loo_result = loo(loglik, r_eff = r_eff)
waic_result = waic(loglik)

### Calculating the log likelihood on the test_data ###

# Finding posterior mean for the theta variables from the above sample
theta_start = 2
theta_end = M*K + 1
theta = matrix(summ$mean[theta_start:theta_end], nrow = M, ncol = K, byrow = FALSE)
col_names_signature = vector(length = K)
for(i in 1:K){
  col_names_signature[[i]] = glue("{i}")
}
theta_df = setNames(data.frame(theta), col_names_signature)

# Finding the posterior mean for the phi values
parameters = list()
previous_end = theta_end
for(i in 1:num_f){
  feature_start = previous_end + 1
  feature_end = (feature_start - 1) + T_c[i]*K
  parameter_matrix = matrix(summ$mean[feature_start: feature_end], nrow = K, ncol = T_c[i], byrow = FALSE)
  feature_name = glue("phi_c{var_select[i+1]}_df") 
  previous_end = feature_end
  if(var_select[i+1] == 2){
    parameters[[feature_name]] = setNames(data.frame(parameter_matrix), c("DEL", "DUP", "INV"))
  } else {
    parameter_names = vector(length = T_c[i])
    for(j in 1:T_c[i]){
      parameter_names[j] = glue("Level {j}")
    }
    parameters[[feature_name]] = setNames(data.frame(parameter_matrix), parameter_names)
  }
}

# Defining necessary variables
pat_test = as.integer(unlist(test_data[,1]))
N_test = nrow(test_data)
print(glue("Number of obs in N_test: {N_test}"))
log_lik = 0

# Calculating log likelihood for the test data
for(n in 1:N_test){
  gamma = vector(length=K)
  for(k in 1:K){
    gamma[k] = gamma[k] + log(theta_df[pat_test[n], k])
    for(c in 1:num_f){
      feature_name = glue("phi_c{var_select[c+1]}_df")
      phi_c = parameters[[feature_name]]
      obs_val = as.integer(unlist(test_data[n, c+1]))
      gamma[k] = gamma[k] + log(phi_c[k, obs_val])
    }
  }
  log_lik = log_lik + logSumExp(gamma)
}

print(glue("Log-Likelihood: {log_lik}"))

### Saving the results and arguments ###

# Define arguments
args_list = list(var_select = var_select, 
                 K = K, M = M, used_time = used_time, 
                 T_c = T_c, 
                 log_lik = log_lik, 
                 loo_result = loo_result, r_eff = r_eff,
                 waic_result = waic_result)

# Making function to generate paths
generate_path <- function(type, folder, var_select, K, iterations){
  path <- folder
  path <- paste0(path, "/")
  path <- paste0(path, type)
  path <- paste0(path, "_var-")

  for(i in var_select){
    path <- paste0(path, glue({i}))
  }

  path <- paste0(path, glue("_nsignature-{K}"))

  path <- paste0(path, glue("_niter-{iterations}"))

  return(path)
}

# Creating paths 
path_to_summary <- generate_path("summ", folder, var_select, K, iterations)
path_to_arguments <- generate_path("args", folder, var_select, K, iterations)
path_to_model <- generate_path("model", folder, var_select, K, iterations)

path_to_model = paste0(path_to_model, ".RDS")

# Save summary, arguments and model object 
write.csv(summ, paste0(path_to_summary, ".csv"))
saveRDS(args_list, paste0(path_to_arguments, ".Rdata"))
fit$save_object(file = path_to_model)


