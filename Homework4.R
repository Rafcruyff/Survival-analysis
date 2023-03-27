# For ordering data
library(dplyr)
# Loading the leukemia data 
library(KMsurv)
# For the KM curve
library(survival)
data(drug6mp)
d <- drug6mp
?drug6mp

# The Kaplan Meier estimator without using the survival package.

d <- arrange(d, t1)
for (i in 1:21) {
  d$at_risk[i] <- sum(d$t1 >= d$t1[i])
}
d <- mutate(d, KM = cumprod(1-1/at_risk))

# The Kaplan Meier estimator using the survival package.

KM_surv_object = Surv(d$t1,d$relapse)
KM_surv = survfit(KM_surv_object ~ 1)
plot(KM_surv)
