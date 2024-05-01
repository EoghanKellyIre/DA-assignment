library(survival)
library(dplyr)
library(bnlearn)

dataset = myeloid

# The myeloid data is from R
# id
# subject identifier, 1-646
# 
# trt
# treatment arm A or B
# 
# futime
# time to death or last follow-up
# 
# death
# 1 if futime is a death, 0 for censoring
# 
# txtime
# time to hematropetic stem cell transplant
# 
# crtime
# time to complete response
# 
# rltime
# time to relapse of disease

###########################

#Going Through the first datapoint

###########################

firstObs <- dataset[1,]

#event Time is time to hematopoietic stem cell transplant
#censor Time is time to death or last follow-up

#Set a maximum observed time of tau
tau = 200

# If txtime is NA, it’s replaced with a very large number.The event did not occur during the study period
if(is.na(firstObs[,6]))
{
  firstObs[,6] = 100000000000
}

# If treatment occured it's value would be less then the time of start of study.
# Therefore the observed time would be after treatment
observedTime = min(firstObs[,6], firstObs[,4])

# threshold of interest set. study up to 200 days
observedTime = min(observedTime, tau)

# If the minimum of tau and the event time is greater than the censor time, it
# means that both the event time and tau are beyond the censor time. In other words,
# the event did not occur, and tau did not reach within the observation period.
# Therefore, we don’t have complete information for this individual up to tau, and
# their data is considered censored before tau. Setting it to 0 ensures that these individuals do
# not contribute to the estimation of the survival function beyond their censor time
if((min(tau,firstObs[,6]) >  firstObs[,4]))
{
  observedTime = 0
}

#creating a survival object with the time to event and if censored (death=0)
survivalObj <- Surv(time = firstObs[,4], event = firstObs[,5])

#fitting a Kaplan-Meier survival model to the data
KMSurvModel <- survfit(survivalObj ~ firstObs[,2], data = firstObs, type="kaplan-meier")
plot(KMSurvModel)

#Finds the survival probability at the observed time
givenSurvival <- summary(KMSurvModel, times = observedTime)
survivalProbability <- givenSurvival[[6]]

###########################

#Applying to data set

###########################

# Initialise a vector to hold the survival probabilities
observedTimes <- numeric(nrow(dataset))
head(dataset)

tau=500 #was 200 last paper idk if better or not tbh

# Loop over each row in the myeloid dataset
for(i in 1:nrow(dataset)) {
  obs <- dataset[i,]
  eventTime = dataset[i, 6]
  if(is.na(eventTime)) {
    eventTime = 100000000000
  }
  
  observedTime = min(eventTime, dataset[i, 4])
  observedTime = min(observedTime, tau)
  
  if((min(tau, eventTime) > dataset[i, 4])) {
    observedTime = 0
  }
  
  observedTimes[i] <- observedTime
}

# Add the survival probabilities as a new column in the myeloid dataset
dataset$observedTimes <- observedTimes

survivalObj <- Surv(time = dataset[,4], event = dataset[,5])
KMSurvModel <- survfit(survivalObj ~ dataset[,2], data = dataset, type="kaplan-meier")
plot(KMSurvModel)


# Extract survival times and probabilities
survivalTimesA <- KMSurvModel$time[1:KMSurvModel[["strata"]][["dataset[, 2]=A"]]]
survivalProbsA <- KMSurvModel$surv[1:KMSurvModel[["strata"]][["dataset[, 2]=A"]]]
# Add a time point of 0 with a survival probability of 1 at the start
survivalTimesA <- c(0, survivalTimesA)
survivalProbsA <- c(1, survivalProbsA)

# Extract survival times and probabilities
survivalTimesB <- KMSurvModel$time[-(1:KMSurvModel[["strata"]][["dataset[, 2]=A"]])]
survivalProbsB <- KMSurvModel$surv[-(1:KMSurvModel[["strata"]][["dataset[, 2]=A"]])]
# Add a time point of 0 with a survival probability of 1 at the start
survivalTimesB <- c(0, survivalTimesB)
survivalProbsB <- c(1, survivalProbsB)


# Initialize a vector to store the survival probabilities at the observed times
observedProbs <- numeric(nrow(dataset))

# Loop over each observed time
for(i in 1:nrow(dataset)) {
  observedTime <- dataset$observedTimes[i]
  
  # Extract the survival probability at this time
  if(dataset$trt[i]=="A")
  {
    # Find the largest survival time that is less than or equal to the observed time
    intervalIndex <- max(which(survivalTimesA <= observedTime))
    observedProbs[i] <- survivalProbsA[intervalIndex]
  }
  else
  {
    intervalIndex <- max(which(survivalTimesB <= observedTime))
    observedProbs[i] <- survivalProbsB[intervalIndex]
  }
}

# Add the survival probabilities as a new column in the dataset
dataset$observedProbs <- observedProbs

dataset <- dataset %>%
  mutate(weights = ifelse(observedTimes != 0,(1/observedProbs), 0))

###########################

#Moving forward

###########################
datasetBefore <- dataset

dataset$sex <- ifelse(dataset$sex == "m", 1, 2)
dataset$trt <- ifelse(dataset$trt == "A", 1, 2)
dataset <- dataset[,-c(1, 10, 9)]

weights <- dataset$weights +1
dataset$trt <- weights*dataset$trt
dataset$sex <- weights*dataset$sex
dataset <- dataset[,-c(3, 5, 6, 7, 8) ]

dataset <- dataset[c('sex', 'death','trt')]

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

bn <- hc(dataset)
plot(bn)
fitted <- bn.fit(bn, dataset)
plotD3bn(bn)

