# For ordering data
library(dplyr)
# Loading the survival package
library(survival)

# Q1

# Interval censored data
d <- data.frame(
  start = c(0, 1, 1, 2, 4),
   end = c(3, 2, 3, 5, 6)
)
surv_obj <- with(d, Surv(start, end, type = "interval2"))
surv_fit <- survfit(surv_obj ~ 1)
print(surv_fit)
plot(surv_fit$time, surv_fit$surv)


# Q3

data(bladder1)
