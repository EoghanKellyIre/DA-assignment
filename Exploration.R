# Load the necessary library
library(survival)
library(corrplot)

# Load the data
data(flchain, package="survival")

# View the first few rows of the data
head(flchain)

# Summary of the data
summary(flchain)

######################
#Distributions
######################
# Plotting age distribution
hist(flchain$age, main = "Age Distribution", xlab = "Age", col = "lightblue")
# Gender distribution
barplot(table(flchain$sex), main = "Gender Distribution", xlab = "Gender", ylab = "Count", col = c("pink", "lightblue"))
# Boxplot of creatinine levels
boxplot(flchain$creatinine, main = "Boxplot of Creatinine Levels", ylab = "Creatinine Levels")
# Survival status distribution
barplot(table(flchain$death), main = "Survival Status Distribution", xlab = "Survival Status", ylab = "Count", names.arg = c("Alive", "Dead"), col = c("green", "red"))

######################
#Distributions of death/leave
######################
# Filter the data for those who died
death_data <- flchain[flchain$death == 1,]
hist(death_data$futime, main = "Time Until Death", xlab = "Days from Enrollment Until Death", col = "lightblue")
# Filter the data for those who not dead
alive_data <- flchain[flchain$death == 0,]
#alive_data <- alive_data[alive_data$futime > 4600,]
hist(alive_data$futime, main = "Time Until leave study", xlab = "Days from Enrollment Until leave", col = "lightblue")

######################
#Corr plot
######################
# Select some variables for the correlation plot
df<-flchain[,-c(11)]
df$sex <- ifelse(df$sex == "M", 1, ifelse(df$sex == "F", 2, df$sex))
selected_vars <- df

# Compute the correlation matrix
cor_matrix <- cor(selected_vars, use = "complete.obs")

# Create the correlation plot
corrplot(cor_matrix, method = "circle")
