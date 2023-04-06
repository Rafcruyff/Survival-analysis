# For ordering data
library(dplyr)
# Loading the leukemia data 
library(KMsurv)
# For the KM curve
library(survival)
data(drug6mp)
d <- drug6mp
?drug6mp


# Q1

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
               KM_min = KM^(exp(z*KM_var/(KM*log(KM)))),
               KM_plus = KM^(exp(-z*KM_var/(KM*log(KM)))))



plot(dfKM$KM, type = 's')
lines(t, dfKM$KM_min, type = "s", lty = 2)
lines(t, dfKM$KM_plus, type = "s", lty = 2)

# The Kaplan Meier estimator using the survival package.

KM_surv_object = Surv(d$t1,d$relapse)
KM_surv = survfit(KM_surv_object ~ 1, conf.type = 'log-log')
plot(KM_surv)

# Q2

# For ordering data
library(dplyr)
# Loading the leukemia data 
library(KMsurv)
# For the KM curve
library(survival)
data("kidney")
kidney1 <- kidney[kidney$type==1,]
kidney2 <- kidney[kidney$type==2,]

time_points <- sort(unique(kidney$time))
time_points1 <- sort(unique(kidney1$time))
time_points2 <- sort(unique(kidney2$time))

d <- data.frame(t = time_points)

for (i in 1:length(time_points)){
  time_point = time_points[i]
  d$Y[i] <- sum(kidney$time >= time_point)
  d$d[i] <- sum(kidney[kidney$time == time_point, "delta"])
  d$Y1[i] <- sum(kidney1$time >= time_point)
  d$d1[i] <- sum(kidney1[kidney1$time == time_point, "delta"])
  d$Y2[i] <- sum(kidney2$time >= time_point)
  d$d2[i] <- sum(kidney2[kidney2$time == time_point, "delta"])
}

d <- mutate(d,
            L = (Y1*Y2)/(Y1+Y2),
            int1 = L/Y1*d1,
            int2 = L/Y2*d2,
            int3 = L^2/(Y1*Y2)*(d1+d2))

Z1 = sum(d$int1, na.rm = TRUE)-sum(d$int2, na.rm = TRUE)
V11 = sum(d$int3, na.rm = TRUE)
U1 = Z1/sqrt(V11)
p_value <- 2 * pnorm(-abs(U1))

# Now with the package

surv_obj <- Surv(time = kidney$time, event = kidney$delta)
group <- factor(kidney$type)
survdiff(surv_obj ~ group)


# Q3

# For ordering data
library(dplyr)
# Loading the leukemia data 
library(KMsurv)
# For the KM curve
library(survival)
data("bmt")
kidney1 <- kidney[kidney$type==1,]
kidney2 <- kidney[kidney$type==2,]

time_points <- sort(unique(kidney$time))
time_points1 <- sort(unique(kidney1$time))
time_points2 <- sort(unique(kidney2$time))

d <- data.frame(t = time_points)

for (i in 1:length(time_points)){
  time_point = time_points[i]
  d$Y[i] <- sum(kidney$time >= time_point)
  d$d[i] <- sum(kidney[kidney$time == time_point, "delta"])
  d$Y1[i] <- sum(kidney1$time >= time_point)
  d$d1[i] <- sum(kidney1[kidney1$time == time_point, "delta"])
  d$Y2[i] <- sum(kidney2$time >= time_point)
  d$d2[i] <- sum(kidney2[kidney2$time == time_point, "delta"])
}

d <- mutate(d,
            L = (Y1*Y2)/(Y1+Y2),
            int1 = L/Y1*d1,
            int2 = L/Y2*d2,
            int3 = L^2/(Y1*Y2)*(d1+d2))

Z1 = sum(d$int1, na.rm = TRUE)-sum(d$int2, na.rm = TRUE)
V11 = sum(d$int3, na.rm = TRUE)
U1 = Z1/sqrt(V11)
p_value <- 2 * pnorm(-abs(U1))

# Now with the package

surv_obj <- Surv(time = kidney$time, event = kidney$delta)
group <- factor(kidney$type)
survdiff(surv_obj ~ group)
