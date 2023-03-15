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
# then 100log(lambda/3)-sum_(i=1)^100 x_i  will be log-likelihood.
# and MLE will be
MLE_u = n/sum(x)

# Compute the maximum likelihood estimator
# of the uncensored times using log likelihood
MLE_c = num_events/sum(tt)

# Make a plot of the log-likelihood. first make a fine grid on (0.1,0.3)
grid = 0.1 + (1:200)/1000
# compute the log-likelihoods
l = log(grid)*57-grid*sum(tt) 
# Plot it
plot(grid, l)
# Find the maximum log likelihood
l_max = max(l)
lambda_est = grid[which.max(l)]
# Find lambda_up and lambda_down.
l_deviate = l_max-1.92 
lambda_down = grid[min(which(l > l_deviate))]
lambda_up = grid[max(which(l > l_deviate))]
# Plot in the graph
abline(v = lambda_down, col = "darkgreen")
abline(v = lambda_up, col = "darkgreen")
abline(h = l_deviate, col = "red")

# Now we do the same for the uncensored data.

# compute the log-likelihoods for the uncensored data
l_u = n*log(grid/3)-grid*sum(x) 
# Plot it
plot(grid, l_u)
# Find the maximum log likelihood
l_u_max = max(l_u)
lambda_est = grid[which.max(l_u)]
# Find lambda_up and lambda_down.
l_u_deviate = l_u_max-1.92 
lambda_u_down = grid[min(which(l_u > l_u_deviate))]
lambda_u_up = grid[max(which(l_u > l_u_deviate))]
# Plot in the graph
abline(v = lambda_u_down, col = "darkgreen")
abline(v = lambda_u_up, col = "darkgreen")
abline(h = l_u_deviate, col = "red")

