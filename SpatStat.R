library(spatstat)
library(stpp)
library(dplyr)

# Simulate a realisation of the inhomogeneous Poisson process at [0, 10]2 with the following
# properties:
# the points form two groups
# the centres of the groups have random locations
# for each group, all distances between the points and the centre are less than 2.
# Repeat the simulations twice and produce two plots of simulated points.

# read this https://hpaulkeeler.com/simulating-an-inhomogeneous-poisson-point-process/

# Setup plotting window to show two plots side by side
par(mfrow= c(1,2))

# Generate random coordinates - 4 of them between 0-10 which will be used as centre coordinates
a <- runif(4, min = 0, max = 10)
a

# Define first intensity function for a poisson process
lambda1 <- function(x,y) {
  50 * as.numeric((abs(x - a[1]) < 2) & (abs(y - a[2]) < 2))
}
# Returns 50 if a point is within 2 units of a1,a2 otherwise 0
# So we are creating a square reigon centred at a1,a2 where points are likely to appear a intensity = 50 and 0 elsewhere

# Simulate and plot the first point patter
plot(rpoispp(lambda1, win = owin(c(0,10), c(0,10))))
# owin defines a 10x10 square window

# Now we can define the second intensity function by updating lambda1 for the other plot
lambda1 <- function(x, y) {
  50 * as.numeric((abs(x - a[3]) < 2) & (abs(y - a[4]) < 2))
}
# Defines a new reigon centred at a3,a4


# Simulate and plot the second point pattern
plot(rpoispp(lambda1, win = owin(c(0,10), c(0,10))), add = TRUE)
# Simualtes another point pattern in the second cluster, add = true overlays this new pattern on top of the other so you can see
# the new one in the same plot

# Now we can repeat the process for the second plot
# Generate the random numbers
a <- runif(4, min = 0, max = 10)
a
lambda1 <- function(x, y) {
  50 * as.numeric((abs(x - a[1]) < 2) & (abs(y - a[2]) < 2))
}
plot(rpoispp(lambda1, win=owin(c(0,10),c(0,10))))
lambda1 <- function(x, y) {
  50 * as.numeric((abs(x - a[3]) < 2) & (abs(y - a[4]) < 2))
}
plot(rpoispp(lambda1, win=owin(c(0,10),c(0,10))), add=TRUE)




  













