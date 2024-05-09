############################################################
#
# IPCW BN on Flchain
# 
#
############################################################

########################
#
# Libraries
#
########################

library(survival)
library(dplyr)
library(bnlearn)
library(readr)
library(caret)
library(ggsurvfit)

# Set a seed for reproducibility
set.seed(123)

########################
#
# Loading Data
#
########################
# Load the data
data(flchain, package="survival")
df <- flchain
df$flc.grp <- as.factor(df$flc.grp)
df$death <- as.numeric(df$death)
df = df[-27,]
df = df[,-c(3,7,8,11)]
dataset <- df
rows <- numeric(nrow(dataset))
########################
#
# Getting IPWC Weights
#
########################

# Initialise a vector to hold the survival probabilities
observedTimes <- rows
observedDeaths <- rows
head(dataset)

tau=4600 # Censored past this

# Getting the observedTime and observedDeath
for(i in 1:nrow(dataset)) {
  obs = dataset[i,]
  observedTime = obs$futime
  observedDeath = obs$death
  if (observedTime>tau)
    # If they died or are alive past tau, they are censored at tau.
    # They are observed alive at tao and their time is tau.
  {
    observedTime=tau
    observedDeath=0
  }
  if (observedTime<tau && obs$death==0)
    # If they are alive before tau but left the study they are censored at their futime.
    # They are observed alive at their futime and their observed time is their futime.
  {
    observedTime=observedTime
    observedDeath=0
  }
  if (observedTime<tau && obs$death==1)
    # If they are dead before tau, they are not censored.
    # They are observed dead at their futime and their observed time is their futime.
  {
    observedTime=observedTime
    observedDeath=1
  }
  observedTimes[i] <- observedTime
  observedDeaths[i] <- observedDeath
}

# Add the observedTimes and observedDeaths as a new column
dataset$observedTimes <- observedTimes
dataset$observedDeaths <- observedDeaths

# Find survival probabilities
survivalObj <- Surv(time = dataset$futime, event = dataset$death)
KMSurvModel <- survfit(survivalObj ~ dataset$flc.grp, data = dataset, type="kaplan-meier")

# Get the number of groups
num_groups <- length(levels(dataset$flc.grp))
# Initialize lists to store the survival times and probabilities for each group
survivalTimesList <- vector("list", num_groups)
survivalProbsList <- vector("list", num_groups)

currentStart = 1
# Loop over each group
for (i in 1:num_groups) {
  # Extract survival times and probabilities
  survivalTimesList[[i]] <- KMSurvModel$time[currentStart:(currentStart - 1 + KMSurvModel[["strata"]][paste0("dataset$flc.grp=", i)])]
  survivalProbsList[[i]] <- KMSurvModel$surv[currentStart:(currentStart - 1 + KMSurvModel[["strata"]][paste0("dataset$flc.grp=", i)])]
  currentStart = currentStart + KMSurvModel[["strata"]][paste0("dataset$flc.grp=", i)]
  
  # Add a time point of 0 with a survival probability of 1 at the start
  survivalTimesList[[i]] <- c(0, survivalTimesList[[i]])
  survivalProbsList[[i]] <- c(1, survivalProbsList[[i]])
}

# Initialize a vector to store the survival probabilities at the observed times
observedProbs <- numeric(nrow(dataset))

for (i in 1:nrow(dataset)) {
  # Get the group and observed time for this row
  group <- dataset$flc.grp[i]
  observedTime <- dataset$observedTimes[i]
  
  # Get the survival times and probabilities for this group
  survivalTimes <- survivalTimesList[[group]]
  survivalProbs <- survivalProbsList[[group]]
  
  # Find the largest survival time that is less than or equal to the observed time
  maxSurvivalTimeIndex <- max(which(survivalTimes <= observedTime))
  
  # Get the corresponding survival probability
  observedProbs[i] <- survivalProbs[maxSurvivalTimeIndex]
  if(dataset$observedDeath[i]==1)
    #If they are not censored their probability of them being uncensored (i.e., the event occurring) is 1
  {
    observedProbs[i] <- 1
  }
}

# Add the survival probabilities as a new column in the dataset
dataset$observedProbs <- observedProbs

#Getting the inverse probabilities for the weights
dataset <- dataset %>%
  mutate(weights = ifelse(observedTimes != 0,(1/observedProbs), 0))

########################
#
# BNs 5 fold cross validation
#
########################

#BN IPCW

datasetW = dataset[,c(1:5,9)]
weights <- (dataset$weights)

datasetW[] <- lapply(datasetW, function(x) as.numeric(x))
datasetW$sex <- as.factor(datasetW$sex)

# Resample the data according to the weights
resampled_samples <- datasetW[sample(nrow(datasetW), size = sum(weights)*100, replace = TRUE, prob = weights), ]

#Defining the DAG
dag <- model2network("[kappa][lambda][sex][flc.grp|sex:age:lambda:kappa][age][observedDeaths|sex:flc.grp:age]")

# Now fit the network with the resampled data
bnIPCW <- bn.fit(dag, data = resampled_samples)

# Predict the test data and calculate the mean squared prediction error
pred <- predict(bnIPCW, node = "observedDeaths", data = resampled_samples)
mspe <- mean((pred - resampled_samples$observedDeaths)^2)
mspe

# Create 5 folds
folds <- createFolds(resampled_samples$observedDeaths, k = 5)
total_cvmspe <- 0

for (i in 1:length(folds)) {
  # Split the data into training and test sets
  train_data <- resampled_samples[-folds[[i]],]
  test_data  <- resampled_samples[folds[[i]],]
  
  # Fit the model on the training data
  bn <- bn.fit(dag, data = train_data)
  
  # Predict the test data and calculate the mean squared prediction error
  pred <- predict(bn, node = "observedDeaths", data = test_data)
  cvmspe <- mean((pred - test_data$observedDeaths)^2)
  total_cvmspe <- total_cvmspe + cvmspe
}
x<-total_cvmspe/5
x

#IF no weights

bnNoWeights <- bn.fit(dag, data = datasetW)

# Predict the test data and calculate the mean squared prediction error
pred <- predict(bnNoWeights, node = "observedDeaths", data = datasetW)
mspe <- mean((pred - datasetW$observedDeaths)^2)
mspe

# Create 5 folds
folds <- createFolds(datasetW$observedDeaths, k = 5)
total_cvmspe <- 0
for (i in 1:length(folds)) {
  
  # Split the data into training and test sets
  train_data <- datasetW[-folds[[i]],]
  test_data  <- datasetW[folds[[i]],]
  
  # Fit the model on the training data
  bn <- bn.fit(dag, data = train_data)
  
  # Predict the test data and calculate the mean squared prediction error
  pred <- predict(bn, node = "observedDeaths", data = test_data)
  cvmspe <- mean((pred - test_data$observedDeaths)^2)
  total_cvmspe <- total_cvmspe + cvmspe
}
x<-total_cvmspe/5
x

#IF remove censoring

datasetNoCensored= dataset[dataset$observedDeaths == 1,]
datasetNoCensored <- datasetNoCensored[,c(1:5,9)]
datasetNoCensored[] <- lapply(datasetNoCensored, function(x) as.numeric(x))
datasetNoCensored$sex <- as.factor(datasetNoCensored$sex)
bnNoCensored <- bn.fit(dag, data = datasetNoCensored)

# Predict the test data and calculate the mean squared prediction error
pred <- predict(bnNoCensored, node = "observedDeaths", data = datasetW)
mspe <- mean((pred - datasetW$observedDeaths)^2)
mspe
########################
#
# Plotting
#
########################

pdf("DAG.pdf", width = 10, height = 10)
dag <- model2network("[Kappa][Lambda][Sex][FLC Group|Sex:Age:Lambda:Kappa][Age][Observed Deaths|Sex:FLC Group:Age]")
plot(dag)
dev.off()

#Plotting curves
pdf("kaplan-meier.pdf", width = 10, height = 10)
survfit(survivalObj ~ dataset$flc.grp, data = dataset, type="kaplan-meier") %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + add_confidence_interval()
dev.off()