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

t <- 1:max(d$t1)
dfKM <- data.frame(t=t) %>%
  mutate(at_risk = 0,
         events = 0,
         censor = 0)

for (i in t){
  dfKM$at_risk[i] <- sum(d$t1 >= i)
  dfKM$events[i] <- nrow(d[d$t1 == i & d$relapse == 1,])
  dfKM$censor[i] <- nrow(d[d$t1 == i & d$relapse == 0,])
}

z <- qnorm(.975)

dfKM <- mutate(dfKM,
               KM = cumprod(1-events/at_risk),
               KM_var = sqrt(KM^2*cumsum(events/(at_risk*(at_risk-events)))),
               KM_min = KM^(exp(KM_var/(KM*log(KM)))),
               KM_plus = KM^(exp(-KM_var/(KM*log(KM)))))



plot(dfKM$KM, type = 's')
lines(t, dfKM$KM_min, type = "s", lty = 2)
lines(t, dfKM$KM_plus, type = "s", lty = 2)

# The Kaplan Meier estimator using the survival package.

KM_surv_object = Surv(d$t1,d$relapse)
KM_surv = survfit(KM_surv_object ~ 1)
plot(KM_surv)
