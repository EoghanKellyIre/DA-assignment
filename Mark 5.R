library(survival)
library(dplyr)
library(bnlearn)
library(readr)

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
library(ggsurvfit)
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

############################
#Old

###############
datasetW = dataset[,c(2,6,12,8, 14)]
write.csv(datasetW, file = "dataFromR.csv")
datasetW = dataset[,c(2,6,12,8)]
weights <- dataset$weights +1
datasetW$ageBucketed <- weights*datasetW$ageBucketed
datasetW$sex <- weights*datasetW$sex
  
bn <- hc(datasetW)

plot(bn)

# Convert all integer columns to numeric
datasetW[] <- lapply(datasetW, function(x) as.numeric(as.character(x)))
fitted <- bn.fit(bn, datasetW, weights=weights)

plotD3bn(bn)


# Resample the data according to the weights
resampled_samples <- samples[sample(nrow(samples), size = sum(weights), replace = TRUE, prob = weights), ]

# Now fit the network with the resampled data
bn <- bn.fit(bn, data = resampled_samples)



##############################
# New Code
#############################



pvec_zi_given_alive = numeric(ncol(datasetW)-1)

for (i in 1:(3))
{
  for (level in 1:(2)) #Sex
  {
      sum = 0 
      for (j in 1:nrow(datasetW))
      {
        if( datasetW[j, 4] == 0 )
        {
          if(datasetW[j,i]  == level){ sum = sum + 1}
        }
      } 
  }
}

# change the numbers to suit whatever col the nodes correspond to 
nodes = c(1,2,3) # sex, grp, age
numCategNode = c(2,2,10)


### corresp to Eq 14
P_ziGivenDoA <- function(i, statusDoA)
{
  e = statusDoA # dead or alive
  sum1 = 0
  
  pVec = zeros(numCategNode[i])
  
  # get numerator
  for (level in 1:numCategNode[i])
  {
    sum1 = 0
    for (j in 1:n)
    {
      if((datasetW[j,4] == e)  && (datasetW[j,i] == e)){ sum1 = sum1 + datasetW[j,5]}
    }
    pvec[level] <- sum1 
  }
  
  # get denominator
  sum2 = 0
  for (j in 1:n)
  {
    if (datasetW[j,i] == e) { sum2 = sum2 + datasetW[j,5]} 
  }   
  
  pvec <- pvec / sum2
  return(pvec)
}


library(reticulate)
use_condaenv("DAProject", required = TRUE)

pomegranate <- import("pomegranate")
BayesianNetwork <- pomegranate$BayesianNetwork$from_samples

model <- pomegranate$BayesianNetwork$from_samples(datasetW, weights=weights, algorithm = 'exact')

# Resample the data according to the weights
resampled_samples <- datasetW[sample(nrow(datasetW), size = sum(weights), replace = TRUE, prob = weights)]

bn <- hc(datasetW)

# Now fit the network with the resampled data
bn <- bn.fit(bn, data = resampled_samples)

plotD3bn(bn)



###############################
#Other New

#########################


datasetW = dataset[,c(2,6,12,8)]
weights <- dataset$weights

# Resample the data according to the weights
resampled_samples <- datasetW[sample(nrow(datasetW), size = sum(weights), replace = TRUE, prob = weights), ]

datasetW[] <- lapply(datasetW, function(x) as.numeric(as.character(x)))

dag <- model2network("[sex][flc.grp][ageBucketed|sex:flc.grp][death|sex:flc.grp:ageBucketed]")

# Now fit the network with the resampled data
bn <- bn.fit(dag, data = resampled_samples)

plotD3bn(bn)

