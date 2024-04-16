library(bnlearn)
library(tidyverse)
library(lubridate)
library(networkD3)
library(readr)
library(plyr)


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
df$sex <- ifelse(df$sex == "M", 1, ifelse(df$sex == "F", 0, df$sex))
df$sex <- as.numeric(as.character(df$sex))
df$chapter <- ifelse(df$death == 0, 'Not Dead', df$chapter)
df$chapter <- as.factor(df$chapter)

# Remove rows with NA values
df <- na.omit(df)

df <- data.frame(lapply(df, function(x) {
  if (is.double(x)) as.numeric(x) else x
}))

bn <- hc(df)
bn1 <- mmhc(df)
bn2 <- tabu(df)
bn3 <- pc.stable(df)


plot(bn)
#plot(bn1)
plot(bn2)
plot(bn3)



fitted <- bn.fit(bn, df)
#fitted1 <- bn.fit(bn1, df)
fitted2 <- bn.fit(bn2, df)
fitted3 <- bn.fit(bn3, df)
fitted4 <- bn.fit(bn4, df)

plotD3bn(bn)

cv.bnlearn(df, bn, loss="mse", k=1)



actual = df[9,]
predicted_age <- predict(fitted, node = "futime", data = actual)
predicted_age1 <- predict(fitted2, node = "futime", data = actual)
predicted_age2 <- predict(fitted3, node = "futime", data = actual)
predicted_age
predicted_age1
predicted_age2
actual[9]

##############
#hc
##############

# Split the data into training and test sets
set.seed(123)  # for reproducibility
train_indices1 <- sample(1:nrow(df), nrow(df) * 0.8)
train_set1 <- df[train_indices1, ]
test_set1 <- df[-train_indices1, ]

# Fit the model on the training set
bn1_train1 <- hc(train_set1)
fitted1_train <- bn.fit(bn1_train1, train_set1)

# Predict on the test set
predicted_age1 <- predict(fitted1_train, node = "age", data = test_set1)

# Compute prediction error
error1 <- mean((predicted_age1 - test_set1$age)^2)

##############
#pc.stable
##############

# Split the data into training and test sets
train_indices2 <- sample(1:nrow(df), nrow(df) * 0.8)
train_set2 <- df[train_indices2, ]
test_set2 <- df[-train_indices2, ]

# Fit the model on the training set
bn1_train2 <- pc.stable(train_set2)
fitted1_train2 <- bn.fit(bn1_train2, train_set2)

# Predict on the test set
predicted_age2 <- predict(fitted1_train2, node = "age", data = test_set)

# Compute prediction error
error2 <- mean((predicted_age2 - test_set2$age)^2)

##############
#tabu
##############

train_indices3 <- sample(1:nrow(df), nrow(df) * 0.8)
train_set3 <- df[train_indices3, ]
test_set3 <- df[-train_indices3, ]

# Fit the model on the training set
bn1_train3 <- tabu(train_set3)
fitted1_train3 <- bn.fit(bn1_train3, train_set3)

# Predict on the test set
predicted_age3 <- predict(fitted1_train3, node = "age", data = test_set)

# Compute prediction error
error3 <- mean((predicted_age3 - test_set3$age)^2)

##############
#tabu
##############

train_indices3 <- sample(1:nrow(df), nrow(df) * 0.8)
train_set3 <- df[train_indices3, ]
test_set3 <- df[-train_indices3, ]

# Fit the model on the training set
bn1_train3 <- tabu(train_set3)
fitted1_train3 <- bn.fit(bn1_train3, train_set3)

# Predict on the test set
predicted_age3 <- predict(fitted1_train3, node = "age", data = test_set)

# Compute prediction error
error3 <- mean((predicted_age3 - test_set3$age)^2)

##############
#mmhc
##############

train_indices3 <- sample(1:nrow(df), nrow(df) * 0.8)
train_set3 <- df[train_indices3, ]
test_set3 <- df[-train_indices3, ]

# Fit the model on the training set
bn1_train3 <- mmhc(train_set3)
fitted1_train3 <- bn.fit(bn1_train3, train_set3)

# Predict on the test set
predicted_age3 <- predict(fitted1_train3, node = "age", data = test_set)

# Compute prediction error
error3 <- mean((predicted_age3 - test_set3$age)^2)

##############
#rsmax2
##############

train_indices3 <- sample(1:nrow(df), nrow(df) * 0.8)
train_set3 <- df[train_indices3, ]
test_set3 <- df[-train_indices3, ]

# Fit the model on the training set
bn1_train3 <- rsmax2(train_set3)
fitted1_train3 <- bn.fit(bn1_train3, train_set3)

# Predict on the test set
predicted_age3 <- predict(fitted1_train3, node = "age", data = test_set)

# Compute prediction error
error3 <- mean((predicted_age3 - test_set3$age)^2)


##############
#h2pc
##############

train_indices3 <- sample(1:nrow(df), nrow(df) * 0.8)
train_set3 <- df[train_indices3, ]
test_set3 <- df[-train_indices3, ]

# Fit the model on the training set
bn1_train3 <- h2pc(train_set3)
fitted1_train3 <- bn.fit(bn1_train3, train_set3)

# Predict on the test set
predicted_age3 <- predict(fitted1_train3, node = "age", data = test_set)

# Compute prediction error
error3 <- mean((predicted_age3 - test_set3$age)^2)


##############
#mmpc
##############

train_indices3 <- sample(1:nrow(df), nrow(df) * 0.8)
train_set3 <- df[train_indices3, ]
test_set3 <- df[-train_indices3, ]

# Fit the model on the training set
bn1_train3 <- mmhc(train_set3)
fitted1_train3 <- bn.fit(bn1_train3, train_set3)

# Predict on the test set
predicted_age3 <- predict(fitted1_train3, node = "futime", data = test_set)

# Compute prediction error
error3 <- mean((predicted_age3 - test_set3$futime )^2)

