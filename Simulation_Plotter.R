# Generate simulated patient data and plot the associated
# sliding window point cloud and persistence curves. 

### HEADING
### This lines will only be executed once. They load the R-library TDA
# it will install it if needed.
if (!require(TDA)) install.packages("TDA", dependencies = TRUE)
library(TDA)
library(ggpubr)
library(geometry)
require(ggplot2)
source("TDAsignal.R")
source("persistentcurvemaker.R")
source("landscapemaker.R")
source("simulate_profiles.R")
source("rankfuncmaker.R")
source("curves_from_persist.R")

windowsize=2 

# Change the simulation parameters
lengthofab = 5
meanofab = 1
stddofab = 0.2
dimofstudy = 0
numprobes = 20
aber_start = 1
  
# Generate the test and control data
testprof = generate_test_profile(meanofab, stddofab, numprobes, 1, lengthofab)
axislim = max(testprof)
  
numpoints = length(testprof)
indicesaux=c((1:numpoints),(1:(windowsize-1)))
  
# Create sliding window point cloud
signalauxtest = c(testprof, testprof[1:(windowsize - 1)])
  
test_cloud = matrix(0, nrow = numpoints, ncol = windowsize)
for (i in 1:numpoints) {
  test_cloud[i,] = signalauxtest[i:(i + windowsize - 1)]
}
  
# Computes the bounds of the axes for the point clouds
maxpcx = max(test_cloud[, 1])
minpcx = min(test_cloud[, 1])
maxpcy = max(test_cloud[, 2])
minpcy = min(test_cloud[, 2])
  
# Computes the maximum distance between points in the test and control clouds
# takes the largest of the two divided by 3 for 3 as maximum, after that
# all simplices are included in the VR complex. 
testconvexhullpointindices = unique(c(convhulln(test_cloud)))
testconvexhull = test_cloud[testconvexhullpointindices, ]
maxtestdist = max(c(as.matrix(dist(testconvexhull,upper = TRUE))))

# Construct the betti, lifespan and landscape curves
TDAsignaloutputtest <- apply(matrix(testprof), 2, TDAsignal, dimstudy = dimofstudy)
diag_test = TDAsignaloutputtest[[1]][[3]]
maxfiltrationvalue = max(diag_test[,3])
par_discretization=seq(from=0,to=maxfiltrationvalue,by=0.01)
  
landscapes = landscape_from_persist(diag_test, par_discretization, dimofstudy, 4)
landscapes_2 = landscapes[2,]
landscapes_3 = landscapes[3,]
landscapes_4 = landscapes[4,]

list_profilestest <-lapply(TDAsignaloutputtest, persistentcurvemaker,
                             dimstudy = dimofstudy, par_discretization=par_discretization)
  
# initialize persistence curve matrix
bettiprofilematrixtest <- matrix(data=0, 
                                   nrow=1,
                                   ncol=length(par_discretization))
  
lifespanmatrixtest <- bettiprofilematrixtest

bettiprofilematrixtest[1,] = list_profilestest[[1]]$betti1

lifespanmatrixtest[1,] = list_profilestest[[1]]$lifespan
  
theme = theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=16, face = "bold")) + theme(axis.text = element_text(size=14, face = "bold")) + theme(axis.title = element_text(size=16, face = "bold")) + theme(legend.position=c(0.8, 0.9)) + theme(legend.text=element_text(size=12))
xlabel = xlab("Filtration Parameter")

# Plot the signals of the control and test patient data
signalx = seq(from=1,to=numprobes, by=1)
signaly2 = testprof
signaldf = data.frame(signalx, signaly2)
signalplot = ggplot(signaldf, aes(signalx)) +  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) + geom_point(aes(y=signaly2), colour="blue") + ggtitle("Simulated Copy Number Data") + xlab("Genes") + ylab("Copy Number") + theme
  
# Plot the test patient's sliding window point cloud
test_cloud_x = test_cloud[, 1]
test_cloud_y = test_cloud[, 2]
test_cloud_df = data.frame(test_cloud_x, test_cloud_y)
tcplot = ggplot(test_cloud_df, aes(test_cloud_x)) +  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) + geom_point(aes(y=test_cloud_y), colour = "blue") + ggtitle("Sliding Window Point Cloud of Simulated Data") + xlim(-axislim, axislim) + ylim(-axislim, axislim) + xlab("Copy Number") + ylab("Copy Number") + theme 

# Plot the persistence curves
y2 = bettiprofilematrixtest[1, ]
df = data.frame(par_discretization, y2)
betti_curveplot = ggplot(df, aes(par_discretization)) + geom_line(aes(y=y2), colour="blue", size = 1.5) + ggtitle("Betti Curve for Simulated Data") + xlabel + ylab("Connected Components") + theme
  
y1 = landscapes_2
y2 = landscapes_3
y3 = landscapes_4
df = data.frame(par_discretization, y1, y2, y3)
landscape_curveplot = ggplot(df, aes(par_discretization)) + geom_line(aes(y=y1, colour="Landscape 2"), size = 1.5) + geom_line(aes(y=y2, colour="Landscape 3"), size = 1.5) + geom_line(aes(y=y3, colour="Landscape 4"), size = 1.5) + scale_colour_manual("", 
                                                                                                                                                                                                                                       breaks = c("Landscape 2", "Landscape 3", "Landscape 4"),
                                                                                                                                                                                                                                       values = c("red", "blue", "gray"))  + ggtitle("Landscapes for Simulated Data") + xlabel + ylab("") + theme
  
y2 = lifespanmatrixtest[1, ]
df = data.frame(par_discretization, y2)
lifespan_curveplot = ggplot(df, aes(par_discretization)) + geom_line(aes(y=y2), colour="blue", size = 1.5) + ggtitle("Lifespan Curve for Simulated Data") + xlabel + ylab("Lifespan") + theme


# Output the initial data plot, sw point cloud, betti, lifespan and landscape
# plot
print(signalplot)
print(tcplot)
print(betti_curveplot)
print(lifespan_curveplot)
print(landscape_curveplot)


