library(CubSFS)
library(parallel)
library(tidyverse)

#3. Set the number of cores (use detectCores to use all)
numcores <- detectCores()

#4.Defining all the species specific stuff we'll need for the following code
# mu requires mutations in sites/generation
# genyr is generation time in years
# sfsfilename is the input sfs name
# outputfilename is what you want to call the output (pop size through time)
mu <- 2.9*10^(-9)
genyr <- 0.3
sfsfilename <- "West_MAFpop0.obs"
outputfilename <- "West"

write.table(rbind(cbind("Mutation",mu),
            cbind("Generations in years",genyr),
            cbind("Input file",sfsfilename),
            cbind("Output prefix",outputfilename)),
            paste(outputfilename,"_parameters.txt",sep=""),
            quote = FALSE,row.names = FALSE,col.names = FALSE)

#5. Overwrote the CValpha function to allow for 
# parallelization (as per the tutorial) (in file)
CValpha <- function (SFS.groups, nb.groups, knots, theta.prm, alpha, n.samples, 
                     Wmat, Vmat, is.folded, nb.time = 500, conv.lim = 10^{
                       -3
                     }, loc.solv = "COBYLA", control = list(maxeval = 20000, ftol_abs = 10^(-5), 
                                                            stopval = 0)) 
{
  fit <- mclapply(split(1:nb.groups, as.factor(1:nb.groups)), 
                helpCValpha, SFS.groups = SFS.groups, n.samples = n.samples, 
                knots = knots, theta.prm = theta.prm, alpha = alpha, 
                Wmat = Wmat, Vmat = Vmat, is.folded = is.folded, nb.time = nb.time, 
                conv.lim = conv.lim, loc.solv = loc.solv, control = control, mc.cores=numcores)
  MSE <- mean(unlist(fit))
  cat("alpha =", alpha, "and Mean GOF = ", MSE, "\n")
  return(MSE)
}

# also added an mclapply in Line 16 of estiamteCubSFS.R
#source("CubSFS/CubSFS/R/CValpha.R")
#source("CubSFS/CubSFS/R/estimateCubSFS.R")
#source("CubSFS/CubSFS/R/estimateAlpha.R")
## had an error that stopped AiCc of folded spectra
#source("CubSFS/CubSFS/R/AicCubSFS.R")

#6. Loading in the data
# object should be a named numeric vector, and include monomorphic sites 
# (for now). Note: has to be a rounded SFS: will not work with fractions
# e.g. after downprojecting
obsSFS <- read_table2(sfsfilename,skip = 1)
obsSFS <- as.matrix(obsSFS)[1,]
obsSFS <- round(obsSFS,digits=0)

#7. Describing the number of samples and number of sites
# Removing the mononorphic sites from the spectra
# Because our SFS is folded, removing the "top half" of the SFS
n.samples <- 28
total.nb.sites <- sum(obsSFS)
obsSFS <- obsSFS[-1]
obsSFS <- obsSFS[1:(length(obsSFS)/2)]

#8. Looking at the fit of the oberved to the expected SFS
# n.knots = number of time points demographic history carved up into
# t_m the furtherest back in time estimate in coalescent units
# Expected height of a genealogical tree is equal to 2 coalescent units
# (1 coalescent unit = 2Ne generations, where Ne = number of diploid indivs)
# alpha = smoothing parameter will be estimated for us if we leave this 
# as null

# Now setting seed and running it.
set.seed(98765)
# This step takes a while, and does not print off updates as it goes
starttime <- Sys.time()
resAlpha <- estimateCubSFS(obsSFS,n.samples,n.knots=29,t_m=0.25,is.folded=TRUE)
endttime <- Sys.time()
print("Running time for generating resAlpha")
endttime-starttime
# Resultant plots
resAlpha$CoalRatePlot
resAlpha$SFSPlot
resAlpha$qqPlot
# This gives you the estimate of alpha used (resAlpha$alpha$minimum)
resAlpha$alpha

# Saving the output just in case we have to regenerate pop size through time
# after comparison to other species
saveRDS(resAlpha, paste(outputfilename,"original.resAlpha.rds",sep=""))

old.total.nb.sites <- total.nb.sites

# test plotting
TimeInYears <- seq(0,1000000,by=1000)
resNt <- PlotCubSFS(n.samples,mu,old.total.nb.sites,genyr,resAlpha, abstime=TimeInYears,is.folded = TRUE)

TimeInYears <- seq(0,50000,by=1000)
resNt <- PlotCubSFS(n.samples,mu,old.total.nb.sites,genyr,resAlpha, abstime=TimeInYears, is.folded = TRUE)

TimeInYears <- seq(0,1000,by=100)
resNt <- PlotCubSFS(n.samples,mu,old.total.nb.sites,genyr,resAlpha, abstime=TimeInYears, is.folded=TRUE)

#9. Monomorphic site counts appear to be too high for these
# datasets (probably because the loci filtered out by MAF
# got dumped into the monomorphic category). Going to
# estimate what it actually should be based on expSFS

x <- (seq(1,length(rev(resAlpha$expSFS))))-1
model_data <- data.frame(x=x,y=resAlpha$expSF)

theta.0 <- min(model_data$y)*0.5
model.0 <- lm(log(y-theta.0) ~ x, data=model_data)
alpha.0 <- exp(coef(model.0)[1])
beta.0 <- coef(model.0)[2]
start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)  
model <- nls(y ~ (alpha * exp(beta * x) + theta), data=model_data, start=start)  
plot(model_data$x, model_data$y)
lines(model_data$x, predict(model,list(x=model_data$x)),col='skyblue',lwd=3)
total.nb.sites <- round(predict(model,list(x=-1))+sum(obsSFS))  

# test plotting
TimeInYears <- seq(0,1000000,by=1000)
resNt <- PlotCubSFS(n.samples,mu,total.nb.sites,genyr,resAlpha, abstime=TimeInYears,is.folded = TRUE)

TimeInYears <- seq(0,50000,by=1000)
resNt <- PlotCubSFS(n.samples,mu,total.nb.sites,genyr,resAlpha, abstime=TimeInYears, is.folded = TRUE)

TimeInYears <- seq(0,1000,by=100)
resNt <- PlotCubSFS(n.samples,mu,total.nb.sites,genyr,resAlpha, abstime=TimeInYears, is.folded=TRUE)

#10. Because of MAF filtering applied across all of these datasets
# counts for low frequency variants are much too low. Therefore,
# instead of doing the analyses on the obsSFS, we use the expected SFS
# extracted from the previous run for all sites < 0.01 in frequency
cutoff <- ceiling(0.01*n.samples)
obsSFS2 <- obsSFS
if (cutoff >= 1) {
  obsSFS2[1:cutoff] <- round(resAlpha$expSFS,digits=0)[1:cutoff]
  names(obsSFS2) <- names(resAlpha$SFS)
  set.seed(98765)
  starttime <- Sys.time()
  resAlpha2 <- estimateCubSFS(obsSFS2,n.samples,n.knots=29,t_m=0.25,is.folded=TRUE)
  endttime <- Sys.time()
  print("Running time for generating resAlpha2")
  endttime-starttime
} else {
  resAlpha2 <- resAlpha
}
# Resultant plots
resAlpha2$CoalRatePlot
resAlpha2$SFSPlot

resAlpha2$qqPlot

# This gives you the estimate of alpha (resAlpha$alpha$minimum)
resAlpha2$alpha

# Saving the output just in case we have to regenerate pop size through time
# after comparison to other species
saveRDS(resAlpha2, paste(outputfilename,"estimated.resAlpha2.rds",sep=""))

# test plotting
TimeInYears <- seq(0,1000000,by=1000)
resNt <- PlotCubSFS(n.samples,mu,total.nb.sites,genyr,resAlpha2, abstime=TimeInYears,is.folded = TRUE)

TimeInYears <- seq(0,50000,by=1000)
resNt <- PlotCubSFS(n.samples,mu,total.nb.sites,genyr,resAlpha2, abstime=TimeInYears, is.folded = TRUE)

TimeInYears <- seq(0,1000,by=100)
resNt <- PlotCubSFS(n.samples,mu,total.nb.sites,genyr,resAlpha2, abstime=TimeInYears, is.folded=TRUE)

#11. Bootstrapping to get uncertainty estimates going back in time
# Following code generates the bootstrap estimates for 'tweaked MAF' SFS
set.seed(98765)
# Number of bootstrap replicates. Due to time constraints, suggest 10, not 50 (they run sequentially)
nb.boot <- 10
# Resampling sites for the bootstrap replicates
boot <- rmultinom(nb.boot,total.nb.sites,c(total.nb.sites-sum(obsSFS2),obsSFS2)/total.nb.sites)
# Stripping the monomorphic sites off the top of the bootstrap replicates
boot <- boot[-1,]

# This will generate demographic reconstructions for each bootstrap
# If it gets stuck and doesn't finish in a "decent" amount of time
# compared to the estimate, probably need to kill and resimulate
# that replicate (probably b/c of an issue similar to: 
# https://github.com/facebook/prophet/issues/842)
# You'll know if it has stalled, because CPU use will go almost to 0 for rsession
print("Estimated running time for generating all boots")
(endttime-starttime)*nb.boot
starttime <- Sys.time()
boot.res <- list()
boot.res[[1]] <- estimateCubSFS(boot[,1],n.samples=n.samples,n.knots=29,t_m=0.25,is.folded=TRUE)
endttime <- Sys.time()
print("This rep took ")
endttime-starttime

# We are doing a loop rather than apply because some boots can cause the optimisation 
# to hang for whatever reason. This way if you exit the loop, you still get some results
# saved rather than losing ALL of the boots if one of the latter replicates doesn't work
# out. See above for figuring out if you've got stuck.
for (i in 2:nb.boot) {
  print(paste("Replicate ",i," is starting at:",sep=""))
  print(Sys.time())
  print("If it takes much longer than:")
  print(endttime-starttime)
  print("You may need to kill and resimulate that replicate")
  boot.res[[i]] <- tempboot <- estimateCubSFS(boot[,i],n.samples=n.samples,n.knots=29,t_m=0.25,is.folded=TRUE)
}

# If you do need to manually redo some replicates
# probably easiest just to resimulate the boots
# and then manually run them e.g.
#boot.res[[9]] <- estimateCubSFS(boot[,1],n.samples=n.samples,n.knots=29,t_m=0.25,is.folded=TRUE)
#boot.res[[10]] <- estimateCubSFS(boot[,2],n.samples=n.samples,n.knots=29,t_m=0.25,is.folded=TRUE)
#boot.res[[8]] <- estimateCubSFS(boot[,3],n.samples=n.samples,n.knots=29,t_m=0.25,is.folded=TRUE)

#10. Saving ResAlpha and boot.res just in case we want to redo the following
# steps if the timescale/y-axis for any of the other species needs tweaking
saveRDS(boot.res, paste(outputfilename,".boot.res.rds",sep=""))

#11. Can iteratively plot, modifying this to capture whole pop increase
# suggest 250,000 and 50,000 as options - can tweak this if any of the
# species show any far different trends over those periods
TimeInYears <- seq(0,1000000,by=1000)
resNt <- PlotCubSFS(n.samples,mu,total.nb.sites,genyr,resAlpha2,boot.res=boot.res, abstime=TimeInYears,is.folded = TRUE)


#12. Writing out the 50,000 time scale data so we can pull it back 
# in and plot all the species at the same time after they are done
write.table(resNt$Ntdf,paste(outputfilename,"popsize.txt",sep=""),row.names = FALSE, quote=FALSE)


TimeInYears <- seq(0,200000,by=1000)
resNt <- PlotCubSFS(n.samples,mu,total.nb.sites,genyr,resAlpha2,boot.res=boot.res, abstime=TimeInYears, is.folded = TRUE)

write.table(resNt$Ntdf,paste(outputfilename,"popsize_200.txt",sep=""),row.names = FALSE, quote=FALSE)


# Final look see at 1,000
TimeInYears <- seq(0,1000,by=100)
resNt <- PlotCubSFS(n.samples,mu,total.nb.sites,genyr,resAlpha2,boot.res=boot.res, abstime=TimeInYears, is.folded=TRUE)

#13. Saving the AicC for the actual and estimated SFS
write.table(cbind(rbind("Original SFS:",
  AicCubSFS(n.samples, SFS.res=resAlpha, is.folded=TRUE)),
  rbind("Partly estimated SFS:",
  AicCubSFS(n.samples, SFS.res=resAlpha2, is.folded=TRUE))),
  paste(outputfilename,"aicc.txt",sep=""),row.names = FALSE, quote=FALSE)

sessionInfo()



data <- read_table2("Westpopsize.txt")

ggplot(data) +
  geom_ribbon(aes(x=time,ymin=LCL,ymax=UCL,fill='#b5b0b0'),alpha=0.2,color="black",size=0.1) +
  theme_bw(base_size = 4) + scale_y_log10(expand=c(0,0)) + ylab("Effective population size") + xlab ("Time (years)") +
  theme(axis.title=element_text(face="bold"), legend.title = element_text(face="bold")) + theme(aspect.ratio = 1) + 
  theme(legend.key.size = unit(0.5,"line")) +
geom_line(aes(x=time,y=Nt),size=0.5) +
  geom_line(aes(x=time,y=median),size=0.5,linetype="dashed") + xlab ("Time (years)") +
  theme(axis.title=element_text(face="bold"))


ggsave("West.pdf",width=8.7,height=17.8,units="cm")



data <- read_table2("Westpopsize_200.txt")

ggplot(data) +
  geom_ribbon(aes(x=time,ymin=LCL,ymax=UCL,fill='#b5b0b0'),alpha=0.2,color="black",size=0.1) +
  theme_bw(base_size = 4) + scale_y_log10(expand=c(0,0)) + ylab("Effective population size") + xlab ("Time (years)") +
  theme(axis.title=element_text(face="bold"), legend.title = element_text(face="bold")) + theme(aspect.ratio = 1) + 
  theme(legend.key.size = unit(0.5,"line")) +
geom_line(aes(x=time,y=Nt),size=0.5) +
  geom_line(aes(x=time,y=median),size=0.5,linetype="dashed") + xlab ("Time (years)") +
  theme(axis.title=element_text(face="bold"))


ggsave("West_200.pdf",width=8.7,height=17.8,units="cm")


