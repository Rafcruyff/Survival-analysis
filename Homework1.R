set.seed(1234)
# Sample size
n = 100
rate = 0.2
# Sample from Exp(.2)
x = rexp(n, rate)
# Sample from Unif(3,6)
cc = runif(n, min=3, max=6)
# Apply censoring to find the times
tt = pmin(x,cc)
# Find which times are not censored
d = (x <= cc)
# Find how many events there are
num_events = sum(d)
# Compute the maximum likelihood estimator
# of the uncensored times using log likelihood

# For this we take the product of f_X(t_i|lambda)f_Y(d_i) for all i as likelihood
# then 100log(lambda/3)-sum_(i=1)^100 t_i  will be log-likelihood.
# and MLE will be
MLE_u = n/sum(tt)

# Compute the maximum likelihood estimator
# of the uncensored times using log likelihood
MLE_c = num_events/sum(tt)
