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

# Plot them together
plot(NA1$t1, NA1$NA1, type = 's', col='red')
lines(NA2$t2, NA2$NA2, type = 's', col = 'blue')

# For the first 10 months the cumulative hazards look the same, but after that 
# the cum hazard of the placebo group seems to diverge from that of the 6MP group.
# This means that the survival rate of the placebo group should decrease t.o.v. the 6MP group.

# We can check this with the hazard ratio.
t <-  0:max(c(NA1$t1,NA2$t2))
NA1_step <- rep(0,length(t))
NA2_step <-  rep(0,length(t))

sf <- data.frame(t = t, NA1_step = NA1_step, NA2_step = NA2_step)

for (t in t){
  
  }
