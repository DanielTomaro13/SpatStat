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





# Use the data mucosa from SPATSTAT. The data give the locations of the centres of two types
# of cells in a cross-section of the gastric mucosa of a rat. This is a marked point pattern. The
# marks have the levels ‘ECL’ (enterochromaffin) and ‘other’.
# a. Plot the data.
# b. c. d. e. Investigate the intensity of the dataset. Plot the estimated intensity.
# Plot the data points and a contour plot for the estimated intensity in the same figure.
# Simulate a Poisson process with the estimated intensity.
# Separate the data into the sub-patterns of points of types ‘ECL’ and ‘other’ and plot
# their intensities.
# f. Plot the cross-type pair correlation function for ‘ECL’ and ‘other’ marks. Interpret the
# plot.
# g. Remove the marks from the point pattern and perform the spatial Kolmogorov–
# Smirnov test for the uniform distribution of the x coordinate. Then test the uniform
# distribution of the y coordinate. Explain your results.

data(mucosa)
summary(mucosa)
str(mucosa)
plot(mucosa)

# Get the average intensity, the intensity is just the average number of points per unit of area we store it in lamb
lamb <- summary(mucosa)$intensity
lamb

# Divide the space into grid boxes, we will do it into 6 col and 3 rows (6x3) to count how many fall in each and save as Q
quadratcount(mucosa, nx = 6, ny = 3)
Q <- quadratcount(mucosa, nx = 6, ny = 3)

# Let us plot this now
plot(mucosa, cex = 0.5, pch = '+') # We plot mucosa as normal
plot(Q, add = TRUE, cex = 2) # We add the grids over the top of the original plot

# We can now estimate point density using Kernel Smoothing
den <- density(mucosa)
den
plot(den)
# density() gives a smooth estimate of point intensity over the area.
# Brighter/warmer areas = higher density of points.
den1 <- density(mucosa, sigma = 1) # custom smoothing sigma = bandwith controls how smooth the estimate is
plot(den1)
plot(mucosa, add = TRUE, cex = 1) # Add original points on top

contour(den1) # Draw contour lines on the smoothed map, which helps visualize reigons with equal intensity
plot(mucosa, add = TRUE, cex = 0.5) # Add points on top again

# Simulate poisson process with same intensity and window 
pei <- rpoispp(lamb, win = owin(c(0, 0.999653), c(0.0125351, 0.808533)))
pei
plot(pei)
# This generates a random point pattern with constant intensity = lamb in the same area.
# Helps compare to find out if, is the original data random or clustered/regular

# Plot density for each cell type ("ECL" vs "other")
plot(density(split(mucosa)))
# split(mucosa) separates the data by mark (cell type).
# This plots separate intensity maps for each type.

# Analyze interaction between cell types
p <- pcfcross(mucosa, "ECL", "other")
plot(p)
# pcfcross() computes the pair correlation function between ECL and other cells.
# This tells you if the two types are attracted to or repelled from each other at different distances.

# Remove marks to treat all cells as the same type
mucosa_unmarked <- unmark(mucosa)
# Now you’re treating it as an unmarked point pattern, ignoring cell type.

# Test spatial distribution using CDF
cdf.test(mucosa_unmarked, function(x, y) {x})
testx <- cdf.test(mucosa_unmarked, function(x, y) {x})
plot(testx)
# This is a spatial Kolmogorov-Smirnov test.
# It compares the distribution of x coordinates against what’s expected under Complete Spatial Randomness (CSR).
# plot(testx) shows the test result visually.
cdf.test(mucosa_unmarked, function(x, y) {y})
testy <- cdf.test(mucosa_unmarked, function(x, y) {y})
plot(testy)
# Same test, but for y coordinate.
# Helps detect if points are evenly spread across x or y dimensions, or show clustering or bias.






