# MFLDA
Detecting Patterns of Somatic Structural Variation Using Multi-Feature Latent Dirichlet Allocation (MFLDA) 

The `MFLDA_with_CV.R` script uses the `cmdstanr` package to sample from the posteior as defined in the stan script `MFLDA_7F.stan`. The posterior in `MFLDA_7F.stan` is defined as the posterior of the MFLDA model with 7 features (i.e. mutation type, length, cluster, telomere distance, replication time, fragile
site and LINE) as presented in the article. The number of features can be changed by chaning the stan script. `MFLDA_with_CV.R` returns a summary containing the summary statistics from the MCMC sampling (e.g. posterior mean, rhat, ess, etc.) for each parameter. It further saves a list of the arguments used for the specific sample and the model object used for sampling. 

`MFLDA_with_CV.R` takes multiple arguments when called.
1. First the selected features needs to be defined. Here, these are given as their corresponding column number in the dataset seperated by comma. OBS. you always have to select column 1, this is expected to be the patient ID and is used in the model fit.
2. Then the number of signatures needs to be specified.
3. Then the number of iterations
4. Then the refresh rate for the model fit (how often the model show progression while fitting)
5. Then the subfolder in which the result is going to be stored. 
6. Lastly you have to specify how many bins you want each continous feature to be discretized into.

The below example is used to call MFLDA with 7 features, 9 signatures, 700 iterations, 20 refresh rate, "results" as results subfolder, and 5 bins

```Rscript MCLDA_with_CV.R 1,2,3,4,5,7,11,12 9 700 20 results 5```

`MFLDA_with_CV.R` further requires the path to the data. This is defined in the script in line 14. It expects the data to be a table in a tsv file with each row representing a mutation, the first column being unique patient IDs and the remaning columns to be the features of the mutations. If the features in the dataset differ from the ones in the dataset from the article the discretization step in the data preprocessing needs to be changed accordingly. 


SPØRGSMÅL??
1. Skal det script der vise hvordan jeg har fundet significant sammenhæng mellem funde signaturer og driver mutationer/cancer typer også med?
2. Vil det ikke kræve at vi også har data med, da det er meget specifikt for den data vi har tilgængelig? Ellers skal der i hvertfald nok skrives lidt om i koden.
3. Skal vi have data med? Ellers skal jeg nok skrive lidt mere specifikt hvilke features vi havde i vores data. 
