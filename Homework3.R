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
         n_at_risk2 = 1) %>%
  group_by(t1) %>%
  mutate(num_relapse1 = sum(relapse))%>%
  group_by(t2) %>%
  mutate(num_relapse2 = sum(relapse)) %>%
  ungroup()

NA1mice = nelsonaalen(d, t1, relapse)


for (i in 1:21) {
  d$n_at_risk1[i] <- sum(d$t1 >= d$t1[i])
  d$n_at_risk2[i] <- sum(d$t2 >= d$t2[i])
}

d <- group_by(d,t1) %>%
  mutate( 
    NA1 = num_relapse1/n_at_risk1) %>%
  group_by(t2) %>%
  mutate( 
    NA2 = num_relapse2/n_at_risk2) %>%
  ungroup()


NA1 <- unique(d[,c('t1','NA1')]) %>%
  arrange('t1') %>%
  mutate(NA1 = cumsum(NA1))

d <- arrange(d, t2)
NA2mice = nelsonaalen(d, t2, relapse)

NA2 <- unique(d[,c('t2','NA2')]) %>%
  arrange(t2) %>%
  mutate(NA2 = cumsum(NA2))


par(mfrow = c(1, 2))
plot(NA1$t1, NA1$NA1, type = 's')
plot(d$t1, NA1mice, type = 's')

par(mfrow = c(1, 2))
plot(NA2$t2, NA2$NA2, type = 's')
plot(d$t2, NA2mice, type = 's')

t <-  0:max(c(NA1$t1,NA2$t2))
NA1_step <- rep(0,length(t))
NA2_step <-  rep(0,length(t))

sf <- data.frame(t = t, NA1_step = NA1_step, NA2_step = NA2_step)

for (t in NA1$t1){
  sf$NA1_step[which(sf$t >= t)] <- NA1$NA1[which(NA1$t1 == t)]
}

for (t in NA2$t2){
  sf$NA2_step[which(sf$t >= t)] <- NA2$NA2[which(NA2$t2 == t)]
}

plot(sf$t, sf$NA1_step, type = 's', col='red')
lines(sf$t, sf$NA2_step, type = 's', col = 'blue')

# For the first 10 months the cumulative hazards look the same, but after that 
# the cum hazard of the placebo group seems to diverge from that of the 6MP group.
# This means that the survival rate of the placebo group should decrease t.o.v. the 6MP group.

# We can check this with the hazard ratio.
sf <- mutate(sf, 
             haz_ratio = NA1_step/NA2_step)

plot(sf$t, sf$haz_ratio, type = 's', col='green')

# Inconclusive

# The Kaplan Meier estimators:

d <- arrange(d, t1) %>%
  mutate(KM1 = cumprod(1-1/n_at_risk1)) %>%
  arrange(d, t2) %>%
  mutate(KM2 = cumprod(1-1/n_at_risk2))

# and their variance estimators

d <- arrange(d, t1) %>%
  mutate(KM_var1 = KM1^2*cumprod(1/n_at_risk1^2)) %>%
  arrange(d, t2) %>%
  mutate(KM_var2 = KM2^2*cumprod(1/n_at_risk2^2))

# Let's plot them

par(mfrow = c(1, 2))
plot(d$t1, d$KM1)#, type = 's')
plot(d$t2, d$KM2)#, type = 's')
