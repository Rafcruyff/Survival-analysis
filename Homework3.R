# For ordering data
library(dplyr)
# Loading the leukemia data 
library(KMsurv)
data(drug6mp)
d <- drug6mp
?drug6mp
# Check with the mice package.
library(mice)
# Estimating the hazard rate of the placebo and 6mp patients using 
# the Nelson Aalen estimator. 
d <- arrange(d, t1) %>% 
  mutate(n_at_risk1 = 1,
         n_at_risk2 = 1)%>%
  group_by(t1) %>%
  mutate(num_relapse1 = sum(relapse))%>%
  group_by(t2) %>%
  mutate(num_relapse2 = sum(relapse))

for (i in 1:21) {
  d$n_at_risk1[i] <- sum(d$t1 >= d$t1[i])
  d$n_at_risk2[i] <- sum(d$t2 >= d$t2[i])
}

d <- group_by(d,t1) %>%
  mutate( 
    NA1 = num_relapse1/n_at_risk1) %>%
  group_by(t2) %>%
  mutate( 
    NA2 = num_relapse2/n_at_risk2)

group_b

par(mfrow = c(1, 2))
plot(d$t1, d$NA1, type = 's')
plot(d$t1, d$NA1mice, type = 's')

d <- arrange(d, t2 ,desc(relapse))
for (i in 1:21) {
  d$n_at_risk2[i] <- sum(d$t2 >= d$t2[i])
}

d <- mutate(d, 
            NA2 = cumsum(relapse/n_at_risk2),
            NA2mice = nelsonaalen(d, t2, relapse))


par(mfrow = c(1, 2))
plot(d$t2, d$NA2, type = 's')
plot(d$t2, d$NA2mice, type = 's')
