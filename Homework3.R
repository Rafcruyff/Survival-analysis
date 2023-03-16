# Loading discfrail and survival
library(discfrail)
library(survival)
function(nelsonaalen_npdf)

# Loading the leukemia data 
library(KMsurv)
data("drug6mp")
# Data info
?drug6mp
# Compute the Nelson-Aalen estimator for both groups
# First we define a timeline
t = 1:350/10

# Placebo:
t_p = drug6mp['t1']
# Compute number at risk
y_p = rep(21,350)
for (i in t_p){
  y_p[which(y_p>i)] = y_p[which(y_p>i)]-1
}
  
NA_p = nelsonaalen_npdf(t_p)
