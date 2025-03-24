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


# Now we will work with nbfires
data(nbfires)
str(nbfires)

fire_types <- nbfires$marks$fire.type
summary(fire_types)
plot(fire_types)
plot(nbfires)
# This extracts the fire type mark for each point and plots the counts and shows all fire locations on the map

# Seperating points by type
forest_points <- nbfires[nbfires$marks$fire.type == 'forest']
grass_points <- nbfires[nbfires$marks$fire.type == 'grass']
dump_points <- nbfires[nbfires$marks$fire.type == 'dump']
other_points <- nbfires[nbfires$marks$fire.type == 'other']

# Plotting the spatial density for each fire type
plot(density(forest_points))
plot(density(grass_points))
plot(density(dump_points))
plot(density(other_points))
# This uses kernel smoothing to estimate the spatial intensity or density of each type
# Higher intensity means the type is more common

# Empty space distances (distance from random locations to the nearest fire)
forest_points <- nbfires[nbfires$marks$fire.type == "forest"]
emp <- distmap(forest_points)
plot(emp, main = "Empty space distances")
plot(forest_points, add = TRUE)

# distmap() computes a map of empty space distances which is the distance from every location to the nearest forest fire
# plot(emp) shows this as a heatmap with darker colours closer to the fire
# The forest points plot adds fire locations on top of the map

# Analyze spatial distribution of forest fires
Fest(forest_points)
plot(Fest(forest_points))
# F function test examines the empty space distribution which is the distrances from random locations to the nearest point
# This tests for clustering vs randomness
Gest(forest_points)
plot(Gest(forest_points))
# G function test looks at the nearest neighbor distances between the actual points and tells you if points tend to be closer (cluster) or evenly spaced (regular)

# Removing makrs for spatial modeling
forest_points_no_marks <- unmark(forest_points)
plot(forest_points_no_marks)

# Removes fire type labels so we treat the pattern as unmarked and just spatial locations
# This is needed to fit some spatial models

# Fitting a spatial model
fit3rd <- ppm(forest_points_no_marks, ~polynom(x, y, 3))
plot(fit3rd, how = "image", se = FALSE)

# This fits a poisson point process model using a 3rd degree polynomial in x and y. This captures spatial trends in the intensity surface or how fire density
# changes across the reigon. The plot() shows the predicted intensity as a heatmanp and turns off standard error shading

# Let us now work with Spatio Temporal Point Processes using the STPP package

# Defining a spatio temporal intensity function
lbda1 <- function(x, y, t) {
  1000 * (1 - x)^2 * y^2 * t^3
}
# This defines a lambda function which is the intensity or how likely a point is to occur at each x,y,t location
# (1 - x)^2	Higher intensity near x = 0 (left side of x-axis)
# y^2	Higher intensity near y = 1 (top of y-axis)
# t^3	Intensity increases over time (later times = more points)
# 10000	Scales everything up to get more points
# more points appear over time, and mostly in the upper-left region of the space.

# Simulate Spatio Temporal Poisson Point Process
ipp1 <- rpp(lambda = lbda1)
# This simulates a spatio temporal poisson point process using the lambda I defined earlier
# The result ipp1 includes locations and times of all simulated points
# specifically ipp1$xyt contains the (x,y,t) coordinates
str(ipp1)
# This is the internal structure with $xyt the matrix of x,y,t values for each point
# $npoints is the number of points
# $lambda is the intensity function used

# Convert to a 3D point object
ipp2 <- as.3dpoints(ipp1$xyt)
# Converts the xyt matrix into a 3D point object that you can plot in 3D which is useful for visualizing points in space and time together

# Plotting the 3D spatio temporal point pattern
plot(ipp2, pch = 19, mark = TRUE)
# Shows the 3D scatter plot of points and sets settings for the points to make solid and add labels 

# Animate the points over time
animation(ipp1$xyt, runtime = 10)
# This animates the point pattern over time (t).
# runtime = 10 sets animation length to 10 seconds.
# shows how points appear progressively as time increases.

# Plotting the point pattern in 2D
plot(ipp1$xyt)
# This plots the 3D points as a 2D projection, usually plotting x and y.
# Time (t) isn't shown here, but it's still part of the data.

# Plot the distibution of event times
plot(density(ipp1$xyt[,3]))
# ipp1$xyt[,3] extracts the t (time) values from all the points.
# density() computes a kernel density estimate for the time values.
# The plot shows when events were most common over time.







