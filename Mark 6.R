library(survival)
library(dplyr)
library(bnlearn)
library(readr)
library(caret)
library(ggsurvfit)

# Set a seed for reproducibility
set.seed(123)

#from paper
plotD3bn <- function(bn) {
  varNames <- nodes(bn)
  # Nodes should be zero indexed!
  links <- data.frame(arcs(bn)) %>%
    mutate(from = match(from, varNames)-1, to = match(to, varNames)-1, value = 1)
  
  nodes <- data.frame(name = varNames) %>%
    mutate(group = 1, size = 30)
  
  networkD3::forceNetwork(
    Links = links,  
    Nodes = nodes,
    Source = "from",
    Target = "to",
    Value = "value",
    NodeID = "name",
    Group = "group",
    fontSize = 20,
    zoom = TRUE,
    arrows = TRUE,
    bounded = TRUE,
    opacityNoHover = 1
  )
}

df <- read_csv("D:\\Trintiy\\Senior Sophister\\Semester 2\\Data Analysis\\Assignment\\flchain\\flchain.csv")
head(df)

#plot(df$death, df$futime)
#cor.test(df$death, df$futime)

# Convert character variables to factors
df$sex <- ifelse(df$sex == "M", 1, ifelse(df$sex == "F", 2, df$sex))
df$sex <- as.numeric(as.character(df$sex))
#df$chapter <- ifelse(df$death == 0, 'Not Dead', df$chapter)
#df$chapter <- as.factor(df$chapter)
df$flc.grp <- as.factor(df$flc.grp)
df$death <- as.numeric(df$death)
new = df[-27,]
df = new[,-c(1,8,9)]

dataset <- df
rows <- numeric(nrow(dataset))

# Initialise a vector to hold the survival probabilities
observedTimes <- rows
observedDeaths <- rows
ageBucketed <- rows
head(dataset)

tau=4500 # Censoredpast this

# Loop over each row in the dataset
for(i in 1:nrow(dataset)) {
  obs = dataset[i,]
  observedTime = obs$futime
  observedDeath = obs$death
  if (observedTime>tau)
  {
    observedTime=tau
    observedDeath=0
  }
  if (observedTime<tau && obs$death==0)
  {
    observedTime=0
    observedDeath=0
  }
  ageBucketed[i] <- cut(obs$age, breaks = seq(50, 100, by = 5), include.lowest = TRUE, right = FALSE)
  observedTimes[i] <- observedTime
  observedDeaths[i] <- observedDeath
}

# Add the survival probabilities as a new column
dataset$observedTimes <- observedTimes
dataset$observedDeaths <- observedDeaths
dataset$ageBucketed <- ageBucketed

survivalObj <- Surv(time = dataset$futime, event = dataset$death)
KMSurvModel <- survfit(survivalObj ~ dataset$flc.grp, data = dataset, type="kaplan-meier")
plot(KMSurvModel)
survfit(survivalObj ~ dataset$flc.grp, data = dataset, type="kaplan-meier") %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + add_confidence_interval()
  #+ add_risktable()

# Get the number of groups
num_groups <- length(levels(dataset$flc.grp))

# Initialize lists to store the survival times and probabilities for each group
survivalTimesList <- vector("list", num_groups)
survivalProbsList <- vector("list", num_groups)

currentStart = 1
KMSurvModel[["strata"]][["dataset$flc.grp=1"]]
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
}

# Add the survival probabilities as a new column in the dataset
dataset$observedProbs <- observedProbs

dataset <- dataset %>%
  mutate(weights = ifelse(observedTimes != 0,(1/observedProbs), 0))

##############
#BN IPWC
##############

datasetW = dataset[,c(2,4,5,6,12,8)]
weights <- dataset$weights

#datasetW$sex <- as.numeric(datasetW$sex)
#datasetW$ageBucketed <- as.numeric(datasetW$ageBucketed)

# Convert the rest to numeric
#datasetW$kappa <- as.numeric(datasetW$kappa)
#datasetW$lambda <- as.numeric(datasetW$lambda)
#datasetW$flc.grp <- as.numeric(datasetW$flc.grp)
#datasetW$death <- as.numeric(datasetW$death)

datasetW[] <- lapply(datasetW, function(x) as.numeric(x))

# Resample the data according to the weights
resampled_samples <- datasetW[sample(nrow(datasetW), size = sum(weights), replace = TRUE, prob = weights), ]

dag <- model2network("[kappa][lambda][sex][flc.grp|sex:ageBucketed:lambda:kappa][ageBucketed][death|sex:flc.grp:ageBucketed]")

# Now fit the network with the resampled data
bnIPWC <- bn.fit(dag, data = resampled_samples, method = "hard-em-g")

plotD3bn(bnIPWC)

# Create 5 folds
folds <- createFolds(resampled_samples$death, k = 5)
total_cvmspe <- 0
#dag <- model2network("[sex][flc.grp|sex:ageBucketed][ageBucketed|sex][death|sex:flc.grp:ageBucketed]")

for (i in 1:length(folds)) {
  # Split the data into training and test sets
  train_data <- resampled_samples[-folds[[i]],]
  test_data  <- resampled_samples[folds[[i]],]
  
  # Fit the model on the training data
  bnIPWC <- bn.fit(dag, data = train_data)
  
  # Predict the test data and calculate the mean squared prediction error
  pred <- predict(bnIPWC, node = "death", data = test_data, method = "exact")
  cvmspe <- mean((pred - test_data$death)^2)
  total_cvmspe <- total_cvmspe + cvmspe
}
x<-total_cvmspe/5
x
################
#IF no weights
################

oldData <- dataset[,c(2,4,5,6,12,8)]
oldData[] <- lapply(oldData, function(x) as.numeric(as.character(x)))
bnNoWeights <- bn.fit(dag, data = oldData)

#plot(bnNoWeights)
plotD3bn(bnNoWeights)

# Create 5 folds
folds <- createFolds(oldData$death, k = 5)
total_cvmspe <- 0
#dag <- model2network("[sex][flc.grp][ageBucketed|sex:flc.grp][death|sex:flc.grp:ageBucketed]")

for (i in 1:length(folds)) {
  
  # Split the data into training and test sets
  train_data <- oldData[-folds[[i]],]
  test_data  <- oldData[folds[[i]],]
  
  # Fit the model on the training data
  bn <- bn.fit(dag, data = train_data)

  # Predict the test data and calculate the mean squared prediction error
  pred <- predict(bn, node = "death", data = test_data, method = "exact")
  cvmspe <- mean((pred - test_data$death)^2)
  total_cvmspe <- total_cvmspe + cvmspe
}
x<-total_cvmspe/5
x

######################
#IF remove censoring
######################
datasetNoCensored= dataset[dataset['weights'] != 0,]
datasetNoCensored <- datasetNoCensored[,c(2,4,5,6,12,8)]
datasetNoCensored[] <- lapply(datasetNoCensored, function(x) as.numeric(as.character(x)))

bnNoCensored <- bn.fit(dag, data = datasetNoCensored)

#plot(bnNoCensored)
plotD3bn(bnNoCensored)

# Create 5 folds
folds <- createFolds(datasetNoCensored$death, k = 5)
total_cvmspe <- 0
#dag <- model2network("[sex][flc.grp][ageBucketed|sex:flc.grp][death|sex:flc.grp:ageBucketed]")

for (i in 1:length(folds)) {
  
  # Split the data into training and test sets
  train_data <- datasetNoCensored[-folds[[i]],]
  test_data  <- datasetNoCensored[folds[[i]],]
  
  # Fit the model on the training data
  bn <- bn.fit(dag, data = train_data)
  
  # Predict the test data and calculate the mean squared prediction error
  pred <- predict(bn, node = "death", data = test_data, method = "exact")
  cvmspe <- mean((pred - test_data$death)^2)
  total_cvmspe <- total_cvmspe + cvmspe
}
x<-total_cvmspe/5
x


###################
#Comparing
###################
bn_structure_IPCW = bn.net(bnIPWC)
bn_structure_NoCensored = bn.net(bnNoCensored)
bn_structure_NoWeights = bn.net(bnNoWeights)

compare(bn_structure, bn_structure_NoCensored, arcs = TRUE)
hamming(bn_structure_IPCW,bn_structure_NoCensored)
