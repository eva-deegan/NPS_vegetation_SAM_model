# Load packages
package.list <- c("tidyverse", 
                  "jagsUI", #jags wrapper
                  "coda",
                  "MCMCvis") #mcmc output
for(i in package.list){library(i, character.only = T)}

# Load data ---------------------------------------------------------------
# Capture command line arguments for running on high performance computing cluster
args <- commandArgs(trailingOnly = TRUE) 
park_index <- as.numeric(args[1])
covariate <- args[2]
photosynthesis_type <- args[3]


park_list <- c("1ARCH","2BLCA_CURE","3CANY","4CARE",
               "5DINO","6CHIR","7CORO","8GICL","9SAGE","10SAGW","11TONT","12WHSA")


#run one park at a time
selected_park <- park_list[park_index]
cat("Running JAGS model for", selected_park, "\n")
  
#input_file <- file.path(insert prepped data file path) 
park_data <- readRDS(input_file)



#alpha set to 1 for dirch prior
alpha=c(1,1,1,1)
#add to data list
park_data2 = c(park_data, list(alpha))
names(park_data2)[15]=c("alpha")
data_list <- park_data2
#verifying matrix format
data_list[["y"]] <- as.matrix(data_list[["y"]])

# Parameters to monitor ------------------------------------------------------

params <- c(
  'deviance',
  'b0.species',
  'mu.species',
  'p0.species',
  'p.overall',
  'b0.plot',
  'b',
  #re-weighted seasonal weights
  'wWT.bar',
  'wST.bar',
  'wWP.bar',
  'wSP.bar',
  #OG weights
  'wST',
  'wWT',
  'wSP',
  'wWP',
  #first season re-normalized weights
  #based on the number of months surveyed
  'wS1T.hat',
  'wW1T.hat',
  'wS1P.hat',
  'wW1P.hat',
  #antecedent variables
  'mean.AntWintP',
  'mean.AntSummP',
  'sd.antSummP',
  'sd.antWintP',
  'mean.AntWintT',
  'mean.AntSummT',
  'sd.antSummT',
  'sd.antWintT',
  #variation
  'sig.b', 
  'sig.plot',
  'sig.species',
  'mu.b', 
  #goodness of fit
  'Dsum')
# inits --------------------------------------------------------------
#input_file = file.path(insert initials RDS file path)
inits <- readRDS(input_file)


# JAGS model --------------------------------------------------------------
#model = file.path(insert file path to model code)
set.seed(12345)
mod <- jagsUI::jags(data = data_list,
                    inits = inits,
                    model.file = model,
                    parameters.to.save = params,
                    parallel = TRUE,
                    n.chains = 6,
                    n.burnin = 10000,
                    n.iter = 500000,
                    n.thin = 10,
                    DIC = TRUE)


#save as an R data object
current_time <- Sys.time()
timestamp_str <- format(current_time, "%Y%m%d_%H%M%S")

print(timestamp_str)

#saveRDS(mod, file = file.path(insert file path where model should be saved))


# Check convergence -------------------------------------------------------
print("Checking convergence...")
#MCMCtrace(object = mod$samples, 
#  pdf = TRUE, 
 # wd = insert directory to save mcmc plots
#  filename = insert name of pdf file for mcmc plots
#  type = "both",
#  open_pdf  = FALSE)
print("Convergence check complete.")
# Get RHat per parameter ------------------------------------------------
print("Getting Rhat values...")
mcmc <- mod$samples
gel <- gelman.diag(mcmc, multivariate = FALSE)
gell_df <- as.data.frame(gel$psrf)
gell_df$Parameter <- rownames(gell_df)
colnames(gell_df) <- c("Point_est", "Rhat", "Parameter")

print("% of parameters converged (under 1.2)")
pct <- round(length(which(gell_df$Rhat <= 1.2)) / nrow(gell_df) * 100, 2)
print(paste0(pct, "%"))



