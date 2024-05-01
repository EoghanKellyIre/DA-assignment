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
df <- df[,-1]

# Convert character variables to factors
df$sex <- ifelse(df$sex == "M", 1, ifelse(df$sex == "F", 2, df$sex))
df$sex <- as.numeric(as.character(df$sex))
#df$chapter <- ifelse(df$death == 0, 'Not Dead', df$chapter)
#df$chapter <- as.factor(df$chapter)
df$flc.grp <- as.factor(df$flc.grp)
df$death <- as.numeric(df$death)

dataset <- df

# Initialise a vector to hold the survival probabilities
observedTimes <- numeric(nrow(dataset))
head(dataset)

#tau=500 #was 200 last paper idk if better or not tbh

# Loop over each row in the dataset
for(i in 1:nrow(dataset)) {
  observedTime = dataset[i, 9]
  observedTimes[i] <- observedTime
}

# Add the survival probabilities as a new column in the myeloid dataset
dataset$observedTimes <- observedTimes

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

datasetW = dataset[,-c(3:5, 7:9, 11:14)]
weights <- dataset$weights +1
datasetW$age <- weights*datasetW$age
datasetW$sex <- weights*datasetW$sex
  
bn <- hc(datasetW)

plot(bn)

# Convert all integer columns to numeric
#datasetW[] <- lapply(datasetW, function(x) as.numeric(as.character(x)))
fitted <- bn.fit(bn, datasetW)

plotD3bn(bn)
