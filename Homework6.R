# Load the survival library
library(survival)

# Reading and preprocessing the data
data = read.table("/home/raf/Survival\ analysis/Oscar.txt", header = TRUE, sep = "\t")
data['Age'] <- data['Final']-data['Birth']

# The survival object
surv_obj <- Surv(data$Age, data$Alive)
# Plotting the overall survival curve
km_fit <- survfit(surv_obj ~ 1)
mean_surv_time <- integrate(function(x) predict(km_fit, x, type = "surv"), 0, max(data$Final))$value

plot(km_fit, xlab = "Time (years)", ylab = "Survival Probability", main = "Kaplan-Meier Survival Curve")


# The log-rank test
group <- factor(data$Wins)
survdiff(surv_obj ~ group)