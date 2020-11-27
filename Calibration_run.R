
#load the needed pakages
library(rstan)
library(rstanarm)
library(raster)
library(MASS)
library(viridis)
library(Metrics)
library(RColorBrewer)
library(plotrix)
library(dplyr)
library(plyr)
library("bayesplot")
library(sROC)
library(waffle)
library(ggsci)
library(wesanderson)
library(truncnorm)
library(Hmisc)
library(ggplot2)
library(ggridges)
library(RColorBrewer)
library(modeest)
library(hrbrthemes)


#q0 prior list
#from Joffre G.I. Agren, D. Gillon, and E. Bosatta, R, Richard Joffre, Gi Ågren, Dominique Gillon, and Ernesto Bosatta. 2001. “Organic Matter Quality in Ecological Studies: Theory Meets Experiment.” Oikos 93: 451–58. https://doi.org/10.1034/j.1600-0706.2001.930310.x.
joffreQ<- c(0.931,
            1.000, 0.965,
            1.079,
            0.941, 1.009,
            1.062, 1.018, 1.067,
            1.045, 1.115,
            0.976, 1.035,
            0.965, 0.966, 0.968, 0.969,
            0.998, 1.000, 1.035, 1.045, 1.085, 1.148, 1.150, 1.178, 1.178, 1.185, 1.191, 1.195, 1.216, 1.222, 1.229, 1.268, 1.356, 1.368, 1.388, 1.465,
            0.998, 1.025, 1.048,
            1.091, 1.119, 1.147,
            1.108, 1.168,
            0.97)

## a few control options for running the models
rstan_options(mc.cores = parallel::detectCores()) # in case you want to run multicore
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


##rescaling functions to be used later on
range01 <- function(x){
  (x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))
  }
rangeMinMax <- function(x, a, b){
  (b-a)*(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))+a
}

##column and row mean functions to be used later ot
 colMax <- function (colData) {
   apply(colData, MARGIN=c(2), max)
 }
 colMin <- function (colData) {
   apply(colData, MARGIN=c(2), min)
 }

 ## Function to add an alpha value to a colour
 add.alpha <- function(col, alpha=1){
   if(missing(col))
     stop("Please provide a vector of colours.")
   apply(sapply(col, col2rgb)/255, 2,
         function(x)
           rgb(x[1], x[2], x[3], alpha=alpha))
 }


## Read the data and declare some objects used later on
 calibration_data<-read.csv("./Data/Dataset_with_coordinates.csv")
 Tarasov_data<-read.csv("./Data/Tarasova.csv")

tree_palette=c("limegreen","deepskyblue", "firebrick1", "darkorange")
parms_palette=brewer.pal(8, "Set1")
decay_pch<-c(21,22,23,24,25,19,1)
decay_pch<-seq(1:8)

legend_decay_class<-c("snag, < 1/3 of tree height broken",
                      "snag, > 1/3 of tree height broken, but standing part >1.3.m",
                      "log, fallen with roots",
                      "log, broken below 1.3 m",
                      "log, man made",
                      "log, butt or other stem part left during logging",
                      "Tarasov et al.")


parameters_list<-c(expression(beta),
                   expression(eta[11]),
                   expression(q[0]),
                   expression(e[0]),
                  expression(f[c]),
                  expression(Delta),
                  expression(paste(epsilon[T[max]])),
                  expression(paste(epsilon[u[0]])))

#plot all the data
hist(calibration_data$cmass.omass)
plot(2003-calibration_data$Year.died, calibration_data$cmass.omass,
     pch=calibration_data$X1.pine.2.spruce.3.birch..Tree.species, col=tree_palette[calibration_data$X1.pine.2.spruce.3.birch..Tree.species])

points(Tarasov_data$t, Tarasov_data$y, pch=16)




## Working on the dataset preparation

#isolate the spruce and pines
no.birch<-calibration_data$X1.pine.2.spruce.3.birch..Tree.species!=3
calibration_data_filtered<-calibration_data[no.birch,]
calibration_data_filtered<-calibration_data
#omit the nas in the target variables
na_list<-(is.na(calibration_data_filtered$Year.died) | is.na(calibration_data_filtered$cmass.omass))
calibration_data_filtered<-calibration_data_filtered[!na_list,]

#plot only spruce and pine
hist(calibration_data_filtered$cmass.omass)
plot(2003-calibration_data_filtered$Year.died, calibration_data_filtered$cmass.omass)

#calculate the maximum deviation of the mass loss by age group
age_vector<-calibration_data_filtered$Age.when.tree.died..Tree.age
sd_vector<-c()
for(i in 1:length(unique(age_vector))){
  sd_vector[i]<- sd(calibration_data_filtered$cmass.omass[calibration_data_filtered$Age.when.tree.died..Tree.age==unique(age_vector)[i]], na.rm = T)
}


#wood density modifier
modifier<-rangeMinMax(calibration_data_filtered$Original.wood.density, 0.8,1.2)

#assemble the calibration data for STAN
Lat<-calibration_data_filtered$N_vector
u0_calc=c((0.0855+0.0157*(50.6-0.768*Lat)), Tarasov_data$u) # calculating u0 outside STAN
calib_data<-list(N=dim(calibration_data_filtered)[1]+dim(Tarasov_data)[1],
                 time_gone=c(2003-calibration_data_filtered$Year.died, Tarasov_data$t),
                 mass_left=c(calibration_data_filtered$cmass.omass, Tarasov_data$y),
                 u0=u0_calc,
                 max_sd=max(sd_vector, na.rm=T),
                 tmax=c(calibration_data_filtered$Stem.diameter, Tarasov_data$tmax),
                 scaled_d=rangeMinMax(c(calibration_data_filtered$Original.wood.density, Tarasov_data$dens), 0.8,1.2))

#assemble the new classifiers
tree_class<-c(calibration_data_filtered$X1.pine.2.spruce.3.birch..Tree.species, rep(4, dim(Tarasov_data)[1]))
decay_class<-c(calibration_data_filtered$X1.snag..less.than.1.3.of.tree.height.broken.2.snag..more.than.1.3.of.tree.height.broken..but.standing.part.higher.than.1.3.m...3.log..fallen.with.roots.4.log..broken.below.1.3.m.5.log..man.made.6.log..butt.or.other.stem.part.left.during.logging...Appearance, rep(7, dim(Tarasov_data)[1]))
decay_class_remapped<- mapvalues(as.factor(decay_class), from = c("1","2","3","4","5","6","7"), to = c("1", "2","3","3","3","3","4"))
legend_decay_class_remapped<-c("snag, < 1/3 broken",
                                "snag",
                                "log",
                                "Tarasov et al.")



#doublecheck the data
which(u0_calc<0)
plot(calib_data$time_gone, calib_data$mass_left)
is.na(calib_data)






##################################################################################
#######################  GENERAL CALIBRATION  ####################################
##################################################################################

#setting the number of iterations
iterations=10000

#running the Stan model and creating the calibration objects
fit <- stan(file = 'dec_model_Q_general.stan', data = calib_data,  chains = 4, iter = iterations, thin=10,  cores=8, control = list(adapt_delta = 0.9))
print(fit)
str(fit)
posteriors<-as.data.frame(fit)
names<-colnames(posteriors)


## plot parameters posteriors and priors

#define the priors
prior_realizations=iterations
beta        = rnorm(prior_realizations, 7,0.7*0.15);
eta_11      = rnorm(prior_realizations, 0.36,0.36*0.15);
q0          = rnorm(prior_realizations,  mean(joffreQ), sd(joffreQ));
e0          = rnorm(prior_realizations, 0.25,0.25*0.15);
fc          = rnorm(prior_realizations, 0.5,0.5*0.15);
delay       = runif(prior_realizations, 0,10);
tmax_error  = rnorm(prior_realizations, 1,1*0.15);
u0_error    = rnorm(prior_realizations, 1,1*0.15);
priors<-list(beta,
             eta_11,
             q0,
             e0,
             fc,
             delay,
             tmax_error,
             u0_error)

posteriors_table<-mat.or.vec(dim(posteriors)[2], 3)

png("chain_test.png")
par(mfrow=c(2,2))
plot(posteriors[,1])
plot(posteriors[,2])
plot(posteriors[,3])
plot(posteriors[,4])
dev.off()

#png("Parameters.png", width=3000, height=2500, res=300)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Parameters.png", width=3000, height=2500, res=300)
par(mfrow=c(3,3))
for(i in 1:(length(names)-2)){
  dist<-posteriors[,i]
  posteriors_table[i, 1]<-mean(posteriors[,i])
  posteriors_table[i, 2:3]<-quantile(posteriors[,i], c(0.025, 0.975))
  plot(density(as.numeric(dist)), main=parameters_list[i], col=parms_palette[i], xlim=range(c(density(as.numeric(dist))$x), density(priors[[i]])$x))
  polygon(density(priors[[i]]),lty=2, col=add.alpha("grey",0.8))
  polygon(density(as.numeric(dist)), col=add.alpha(parms_palette[i],0.8))
}
dev.off()

colnames(posteriors_table)=c("mean", "min", "max")
rownames(posteriors_table)=colnames(posteriors)

write.csv(posteriors_table, "General_parameters.csv")


# thinning the posteriors to plot the simulations
thinning=250
resampling_vector<-seq(from=1, to=iterations*0.1, by=thinning)
resampled_posteriors<-posteriors[resampling_vector,]


# plot the simulations
time_sim=120
time_simulation=seq(from=0, to=time_sim)
time_simulation1=time_simulation

#create the table to hold the posterior probability desnity for each of the data points
RMSE_table<-mat.or.vec(length(calib_data[[2]]), dim(resampled_posteriors)[1])
SS_table<-mat.or.vec(length(calib_data[[2]]), dim(resampled_posteriors)[1])

mass_left_simulated_list<-list()
for(j in 1:length(calib_data[[2]])){

    mass_left_simulated<-mat.or.vec(dim(resampled_posteriors)[1], time_sim+1)

    for( i in 1: dim(resampled_posteriors)[1]){

      zeta = (1-resampled_posteriors$e0[i])/(resampled_posteriors$beta[i]*resampled_posteriors$eta_11[i]*resampled_posteriors$e0[i]);
      alpha = resampled_posteriors$fc[i]*resampled_posteriors$beta[i]*resampled_posteriors$eta_11[i]*u0_calc[j]*resampled_posteriors$u0_error[i]*resampled_posteriors$q0[i]^resampled_posteriors$beta[i];

      #imålement the IF cycle for tmax with a vector
      selection_vector<-(time_simulation<(calib_data$tmax[i]*resampled_posteriors$tmax_error[i]))

      #for all the times when time<tmax
      mass_left_simulated[i,selection_vector] <- ((2/(calib_data$tmax[i]*resampled_posteriors$tmax_error[i]))*(1/(alpha*(1-zeta)))*((1+alpha*time_simulation[selection_vector][selection_vector])^(1-zeta)-
                                                               (1-(time_simulation[selection_vector]/(calib_data$tmax[i]*resampled_posteriors$tmax_error[i]))))+
                              ((2/(calib_data$tmax[i]*resampled_posteriors$tmax_error[i])^2)*(1/(alpha^2*(1-zeta)*(2-zeta)))*(1-(1+alpha*time_simulation[selection_vector])^(2-zeta)))+
                              (1-(time_simulation[selection_vector]/(calib_data$tmax[i]*resampled_posteriors$tmax_error[i])))^2)
      #for all the times when time>=tmax
      mass_left_simulated[i,!selection_vector] <- (2/(calib_data$tmax[i]*resampled_posteriors$tmax_error[i]))*(1/(alpha*(1-zeta)))*(1+alpha*time_simulation[!selection_vector])^(1-zeta)+
                              ((2/((calib_data$tmax[i]*resampled_posteriors$tmax_error[i])^2))*(1/(alpha^2*(1-zeta)*(2-zeta)))*(((1+alpha*(time_simulation[!selection_vector]-(calib_data$tmax[i]*resampled_posteriors$tmax_error[i])))^(2-zeta))-((1+alpha*time_simulation[!selection_vector])^(2-zeta))))

      #calculate the RMSE for that point
      RMSE_table[j,i]<-rmse(mass_left_simulated[i, calib_data$time_gone[j]],calib_data$mass_left[j])
      tmax=calib_data$tmax[i]*resampled_posteriors$tmax_error[i]
      SS_table[j,i]<-((1/alpha)*(1/(zeta-1)))+(tmax/3)


      }

    mass_left_simulated_list[[j]]<-mass_left_simulated
}





#build the big table
mass_left_simulated_table<-mass_left_simulated_list[[1]]

for(j in 1:length(calib_data[[2]])){
  mass_left_simulated_table<-rbind(mass_left_simulated_table,mass_left_simulated_list[[j]])
}
dim(mass_left_simulated_table)

plot(mass_left_simulated_table[,30])

#plotting the lines
#png("Lines_plot.png", width=3500, height=2500, res=300)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Lines_plot.png", width=3500, height=2500, res=300)
plot(2003-calibration_data_filtered$Year.died, calibration_data_filtered$cmass.omass, ylim=c(0,1.03), col="red", xlab="time", ylab="mass remaining", xlim=c(0,120), xaxs="i",yaxs="i", cex=0)
for(j in 1:dim(calibration_data_filtered)[1]){
  for( i in 1: dim(resampled_posteriors)[1]){
    lines( mass_left_simulated_list[[j]][i,])
  }

}
points(calib_data$time_gone, calib_data$mass_left, ylim=c(0,1.1),
       col=tree_palette[tree_class],#[calibration_data_filtered$X1.pine.2.spruce.3.birch..Tree.species],
       bg=add.alpha(tree_palette[calibration_data_filtered$X1.pine.2.spruce.3.birch..Tree.species],0.3),
       pch=decay_pch[decay_class], #calibration_data_filtered$Decay.class+20,
       xlim=c(0,65))
legend("topright", c(legend_decay_class), pch=seq(1:6), bty="n")
legend("bottomleft", c("Pine","Spruce", "Birch", "Spruce (Tarasov et al.)"), pch=21, col=tree_palette, pt.bg=add.alpha(tree_palette,0.3), bty="n")
dev.off()

#plotting the lines
#png("Range_plot.png", width=3500, height=2500, res=300)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Range_plot.png", width=3500, height=2500, res=300)
plot(2003-calibration_data_filtered$Year.died, calibration_data_filtered$cmass.omass, ylim=c(-0.2,1.03), col="red", xlab="Years", ylab="mass remaining", xlim=c(0,120),  xaxs="i",yaxs="i", cex=0)
polygon( c(time_simulation, rev(time_simulation)), c(colMax(mass_left_simulated_table),rev(colMin(mass_left_simulated_table))), col=add.alpha("darkorange",0.4), border = add.alpha("darkorange",0.9))
points(calib_data$time_gone, calib_data$mass_left,
       col=tree_palette[tree_class],#[calibration_data_filtered$X1.pine.2.spruce.3.birch..Tree.species],
       bg=add.alpha(tree_palette[tree_class],0.3),
       pch=decay_pch[decay_class])
legend("topright", c(legend_decay_class), pch=seq(1:8), bty="n")
legend("bottomleft", c("Pine","Spruce", "Birch", "Spruce (Tarasov et al.)"), pch=21, col=tree_palette, pt.bg=add.alpha(tree_palette,0.3), bty="n")
abline(h=0, lty=2)
dev.off()

plot(calib_data$time_gone, calib_data$mass_left)

x_dens<-c()
y_dens<-c()
for(i in 1:dim(mass_left_simulated_table)[2]){
  y_dens<-c(y_dens,mass_left_simulated_table[,i])
  x_dens<-c(x_dens, rep(i, length(mass_left_simulated_table[,i])))
  }

length(x_dens)
length(y_dens)

#plot(x_dens, y_dens)

raster_plot<-kde2d(x_dens, y_dens, n=180, lims = c(c(0,time_sim), c(-0.1,1.3)))


#png("Density_plot.png", width=3500, height=2500, res=300)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Density_plot.png", width=3500, height=2500, res=300)
image(raster_plot, ylim=c(0,1.03), xlim=c(0,120), c(rev(magma(130))), xlab="Years", ylab="Fraction of initial mass")
points(calib_data$time_gone, calib_data$mass_left, ylim=c(0,1.1),
       col=tree_palette[tree_class],#[calibration_data_filtered$X1.pine.2.spruce.3.birch..Tree.species],
       bg=add.alpha(tree_palette[calibration_data_filtered$X1.pine.2.spruce.3.birch..Tree.species],0.4),
       pch=decay_pch[decay_class], #calibration_data_filtered$Decay.class+20,
       xlim=c(0,65))
legend("topright", c(legend_decay_class), pch=seq(1:6), bty="n", text.col="black", col="black", cex=1.3)
legend("bottomright", c("Pine","Spruce", "Birch", "Spruce (Tarasov et al.)"),text.col="black", cex=1.3, pch=21, col=tree_palette, pt.bg=add.alpha(tree_palette,0.2), bty="n")
#abline(h=0, lty=2, col="white")
dev.off()


#png("Density_plot_presentation.png", width=3500, height=2000, res=300)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Density_plot_presentation.png", width=3500, height=2000, res=300)
image(raster_plot, ylim=c(0,1.03), xlim=c(0,120), c(rev(viridis(130))))
points(calib_data$time_gone, calib_data$mass_left, ylim=c(0,1.1),
       col=tree_palette[tree_class],#[calibration_data_filtered$X1.pine.2.spruce.3.birch..Tree.species],
       bg=add.alpha(tree_palette[calibration_data_filtered$X1.pine.2.spruce.3.birch..Tree.species],0.4),
       pch=decay_pch[decay_class], #calibration_data_filtered$Decay.class+20,
       xlim=c(0,65))
legend("topright", c(legend_decay_class), pch=seq(1:6), bty="n", text.col="black", col="black", cex=1.3)
legend("bottomright", c("Pine","Spruce", "Birch", "Spruce (Tarasov et al.)"),text.col="black", cex=1.3, pch=21, col=tree_palette, pt.bg=add.alpha(tree_palette,0.2), bty="n")
#abline(h=0, lty=2, col="white")

dev.off()



##################################################################################
########################  HSY sensitivity  #######################################
##################################################################################

  thinning=1
  resampling_vector.long<-seq(from=1, to=iterations*0.1, by=thinning)
  resampled_posteriors.long<-posteriors[resampling_vector.long,]


  # rebuild a RMSE table
  time_sim=100
  time_simulation=seq(from=0, to=time_sim)
  if(max(calib_data$time_gone)>time_sim){print("Attention, the simulation time is shorter than the longest data point age")}

  #create the table to hold the posterior probability desnity for each of the data points
  RMSE_table.long<-mat.or.vec(length(calib_data[[2]]), dim(resampled_posteriors.long)[1])
  SS_table.long<-mat.or.vec(length(calib_data[[2]]), dim(resampled_posteriors.long)[1])

  mass_left_simulated_list.long<-list()
  for(j in 1:length(calib_data[[2]])){

    mass_left_simulated.long<-mat.or.vec(dim(resampled_posteriors.long)[1], time_sim+1)

    for( i in 1: dim(resampled_posteriors.long)[1]){

      zeta = (1-resampled_posteriors.long$e0[i])/(resampled_posteriors.long$beta[i]*resampled_posteriors.long$eta_11[i]*resampled_posteriors.long$e0[i]);
      alpha = resampled_posteriors.long$fc[i]*resampled_posteriors.long$beta[i]*resampled_posteriors.long$eta_11[i]*u0_calc[j]*resampled_posteriors.long$u0_error[i]*resampled_posteriors.long$q0[i]^resampled_posteriors.long$beta[i]; # alpha is recalculated every time because of varying u0 in different sites

      #imålement the IF cycle for tmax with a vector
      selection_vector<-(time_simulation<(calib_data$tmax[j]*resampled_posteriors.long$tmax_error[i]))


      #for all the times when time<tmax
      mass_left_simulated.long[i,selection_vector] <- ((2/(calib_data$tmax[i]*resampled_posteriors.long$tmax_error[i]))*(1/(alpha*(1-zeta)))*((1+alpha*time_simulation[selection_vector][selection_vector])^(1-zeta)-
                                                                                                 (1-(time_simulation[selection_vector]/(calib_data$tmax[i]*resampled_posteriors.long$tmax_error[i]))))+
                                                    ((2/(calib_data$tmax[i]*resampled_posteriors.long$tmax_error[i])^2)*(1/(alpha^2*(1-zeta)*(2-zeta)))*(1-(1+alpha*time_simulation[selection_vector])^(2-zeta)))+
                                                    (1-(time_simulation[selection_vector]/(calib_data$tmax[i]*resampled_posteriors.long$tmax_error[i])))^2)
      #for all the times when time>=tmax
      mass_left_simulated.long[i,!selection_vector] <- (2/(calib_data$tmax[i]*resampled_posteriors.long$tmax_error[i]))*(1/(alpha*(1-zeta)))*(1+alpha*time_simulation[!selection_vector])^(1-zeta)+
        ((2/((calib_data$tmax[i]*resampled_posteriors.long$tmax_error[i])^2))*(1/(alpha^2*(1-zeta)*(2-zeta)))*(((1+alpha*(time_simulation[!selection_vector]-(calib_data$tmax[i]*resampled_posteriors.long$tmax_error[i])))^(2-zeta))-((1+alpha*time_simulation[!selection_vector])^(2-zeta))))

      #calculate the RMSE for that point
      RMSE_table.long[j,i]<-rmse(mass_left_simulated.long[i, calib_data$time_gone[j]],calib_data$mass_left[j])
      SS_table.long[j,i]<-((1/alpha)*(1/(zeta-1)))+((calib_data$tmax[i]*resampled_posteriors.long$tmax_error[i])/3)


    }

    mass_left_simulated_list.long[[j]]<-mass_left_simulated.long
  }


  dim(RMSE_table.long)

  mlv(SS_table.long, method = "Parzen")
  min(SS_table.long)
  max(SS_table.long)


  #decay.class<-calibration_data_filtered$X1.snag..less.than.1.3.of.tree.height.broken.2.snag..more.than.1.3.of.tree.height.broken..but.standing.part.higher.than.1.3.m...3.log..fallen.with.roots.4.log..broken.below.1.3.m.5.log..man.made.6.log..butt.or.other.stem.part.left.during.logging...Appearance
  #tree.class<-calibration_data_filtered$X1.pine.2.spruce.3.birch..Tree.species
  decay.class<-decay_class
  tree.class<-tree_class
  decay.tree.class<-interaction(decay.class, tree.class)
  RMSE_table_classed<-as.data.frame(cbind(RMSE_table.long, decay.tree.class))
  RMSE_table_classed$decay.tree.class<-decay.tree.class

  length(decay.class)
  length(tree.class)

  class.levels<-levels(as.factor(RMSE_table_classed$decay.tree.class))

  decay.class.text<-legend_decay_class[decay.class]
  unique(decay.class.text)
  tree.class.text<-c("Pine","Spruce","Birch", "Spruce (Tarasov)")[tree.class]
  unique(tree.class.text)
  decay.tree.class.text<-interaction(decay.class.text, tree.class.text)
  unique(decay.tree.class.text)
  class.levels.text<-(unique(as.factor(decay.tree.class.text)))


  which((decay.class.text==levels(as.factor(decay.class.text))[5]) & (tree.class.text==levels(as.factor(tree.class.text))[1]))
  which(decay.tree.class.text==class.levels.text[5])
  which(is.na(decay.tree.class.text==class.levels.text[j]))

  #create the table to hold the results
  HSY_table<-mat.or.vec((dim(resampled_posteriors)[2]-1), length(class.levels.text))
  rownames(HSY_table)<-colnames(resampled_posteriors)[1:(length(colnames(resampled_posteriors))-1)]
  colnames(HSY_table)<-class.levels.text

  for(j in 1:length(class.levels.text)){
      RMSE.table.subclass<-RMSE_table_classed[decay.tree.class.text==class.levels.text[j],]
      dim(RMSE.table.subclass)

      if(dim(RMSE.table.subclass)[1]>1){
      RMSE.mean.subclass<-colMeans(RMSE.table.subclass[,1:((iterations*0.1)/thinning)], na.rm=T)
      quantiles.subclass<-quantile(RMSE.mean.subclass, c(0.1,0.9), na.rm=T)
      which.subclass<-which(RMSE.mean.subclass<quantiles.subclass[1])
      bin1<-resampled_posteriors.long[RMSE.mean.subclass<quantiles.subclass[1],]
      bin2<-resampled_posteriors.long[!RMSE.mean.subclass<quantiles.subclass[1],]

      
      
      for(i in 1:(dim(resampled_posteriors)[2]-1)){
        resampled<-gdata::resample(bin2[,i],100)
        HSY_table[i,j]<-ks.test(bin1[,i],resampled)$statistic
      }
      }else{
        HSY_table[,j]<-rep(NA, (dim(resampled_posteriors)[2]-1))
        }

  }






  plot(ecdf(bin1$beta), col="blue", main=NA, pch=NA, verticals=T)
  plot(ecdf(bin2$beta), col="red", add=T)


  #png(filename = "HSY_CDF.png", height=1500, width=3000, res=300)
  png(filename = "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/HSY_CDF.png", height=1500, width=3000, res=300)
  par(mfrow=c(2,4))
  for(i in 1:8){
  exclude_na1<-is.na(bin1[,i])
  exclude_na2<-is.na(bin2[,i])
  cdf1<-kCDF(bin1[,i][!exclude_na1])
  cdf2<-kCDF(bin2[,i][!exclude_na2])
  plot(cdf1$x, cdf1$Fhat, col="blue", main=parameters_list[i], type="l", ylab="CDF", xlab="Parameter value")
  lines(cdf2$x, cdf2$Fhat, col="red")
  if (i==1){legend("bottomright", c("fit <90% C.I.", "fit <90% C.I."), bty="n", col=c("blue", "red"), lty=1, cex=0.8)}
  }
  dev.off()





  colnames(HSY_table)<-class.levels.text

  HSY_table.scaled<-HSY_table
  for(i in 1:dim(HSY_table)[2]){
    if(!all(is.na(HSY_table[,i]))){
    HSY_table.scaled[,i]<-range01(as.numeric(HSY_table[,i]))
    }else{
      HSY_table.scaled[,i]<-rep(NA, dim(HSY_table.scaled)[1])
    }
  }

  HSY_table.ranked<-HSY_table
  for(i in 1:dim(HSY_table)[2]){
    if(!all(is.na(HSY_table[,i]))){
      HSY_table.ranked[,i]<-rank(as.numeric(HSY_table[,i]))
    }else{
      HSY_table.ranked[,i]<-rep(NA, dim(HSY_table.ranked)[1])
    }
  }



  #png(filename = "HSY_table_scaled.png", height=2800, width=3500, res=350)
  png(filename = "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/HSY_table_scaled.png", height=2800, width=3500, res=350)
  par(mar=c(30,5,2,2))
  color2D.matplot(HSY_table.scaled[1:length(parameters_list),], axes=F, xlab="", ylab="", show.values = 4, vcex = 1, extremes = c("white", "red"))
  axis(1, at = seq_len(ncol(HSY_table.scaled[1:length(parameters_list),])) - 0.5, las=2,
       labels = colnames(HSY_table.scaled[1:length(parameters_list),]), tick = FALSE, cex.axis = 1)
  axis(2, at = seq_len(nrow(HSY_table.scaled[1:length(parameters_list),])) -0.5,
       labels = (parameters_list), tick = FALSE, las = 1, cex.axis = 1)
  dev.off()





  HSY_table.col=HSY_table
  HSY_table.col[is.na(HSY_table.col)]<-0

  map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
  }

  cellcol<-map2color(HSY_table.col, rev(heat.colors(length(HSY_table.col))))
  #png(filename = "HSY_table.png", height=2800, width=3500, res=350)
  png(filename = "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/HSY_table.png", height=2800, width=3500, res=350)
  par(mar=c(25.5,5,2,2))
  color2D.matplot(HSY_table[1:length(parameters_list),], axes=F, xlab="", ylab="", show.values = 3, vcex = 1, extremes = c("white", "red"))#, cellcolors = cellcol)
  axis(1, at = seq_len(ncol(HSY_table[1:length(parameters_list),])) - 0.5, las=2,
       labels = colnames(HSY_table[1:length(parameters_list),]), tick = FALSE, cex.axis = 1)
  axis(2, at = seq_len(nrow(HSY_table[1:length(parameters_list),])) -0.5,
       labels = rev(parameters_list), tick = FALSE, las = 1, cex.axis = 1)
  dev.off()

  HSY_table.col=HSY_table.ranked
  HSY_table.col[is.na(HSY_table.ranked)]<-0
  cols<-unique(as.vector(HSY_table.col))

  #cellcol<-map2color(HSY_table.col, rev(heat.colors(length(cols))))
  cellcol<-map2color(HSY_table.col, c("white",rev(brewer.pal(length(cols), "RdYlBu"))))
  cellcol<-map2color(HSY_table.col,c("white",(viridis(length(cols)))))
  #png(filename = "HSY_table_ranked.png", height=2800, width=2500, res=350)
  png(filename = "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/HSY_table_ranked.png", height=2800, width=2500, res=350)
  par(mar=c(25.5,5,2,2))
  color2D.matplot(HSY_table.ranked[1:length(parameters_list),], axes=F, xlab="", ylab="", show.values = 1, vcex = 1, extremes = c("white", "red"))#, cellcolors = cellcol)
  axis(1, at = seq_len(ncol(HSY_table.ranked[1:length(parameters_list),])) - 0.5, las=2,
       labels = colnames(HSY_table.ranked[1:length(parameters_list),]), tick = FALSE, cex.axis = 1)
  axis(2, at = seq_len(nrow(HSY_table.ranked[1:length(parameters_list),])) -0.5,
       labels = rev(parameters_list), tick = FALSE, las = 1, cex.axis = 1)
  dev.off()



  #count the individuals
  dat <- data.frame(decay.tree.class.text)
  count<-dat %>%
    group_by(decay.tree.class.text) %>%
    summarise(no_rows = length(decay.tree.class.text))

  palette<-c(pal_jama("default")(7),pal_npg("nrc")(9))


  waffle_vec<-count$no_rows
  names(waffle_vec)<-count$decay.tree.class.text
  #png(filename = "Class.count.png", height=1500, width=3600, res=350)
  png(filename = "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Class.count.png", height=1500, width=3600, res=350)
  waffle(waffle_vec, rows=30, col=palette)
  dev.off()








##################################################################################
########################  CALIBRATION BY CLASSES  ################################
##################################################################################


#plotting the lines
#png("Range_plot_simple.png", width=4500, height=3500, res=300)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Range_plot_simple.png", width=4500, height=3500, res=300)

par(mfrow=c(2,2))
for(i in 1:4){
plot(2003-calibration_data_filtered$Year.died, calibration_data_filtered$cmass.omass, ylim=c(-0.2,1.03), col="red", xlab="time", ylab="mass remaining", xlim=c(0,120),  xaxs="i",yaxs="i", cex=0, main=legend_decay_class_remapped[i])
polygon( c(time_simulation1, rev(time_simulation1)), c(colMax(mass_left_simulated_table),rev(colMin(mass_left_simulated_table))), col=add.alpha("darkorange",0.4), border = add.alpha("darkorange",0.9))
selection<-decay_class_remapped==i
points(calib_data$time_gone[selection], calib_data$mass_left[selection],
       col=tree_palette[tree_class][selection],#[calibration_data_filtered$X1.pine.2.spruce.3.birch..Tree.species],
       bg=add.alpha(tree_palette[tree_class][selection],0.3),
        pch=decay_pch[decay_class_remapped][selection])
if(i!=4){legend("topright", c("Pine","Spruce", "Birch"), pch=21, col=tree_palette, pt.bg=add.alpha(tree_palette,0.3), bty="n")}
abline(h=0, lty=2)
legend("bottomleft", paste("(",LETTERS[i],")"), bty="n")
}
dev.off()



decay.class<-decay_class_remapped
tree.class<-tree_class
decay.tree.class<-interaction(decay.class, tree.class)
RMSE_table_classed<-as.data.frame(cbind(RMSE_table.long, decay.tree.class))
RMSE_table_classed$decay.tree.class<-decay.tree.class

length(decay.class)
length(tree.class)

class.levels<-levels(as.factor(RMSE_table_classed$decay.tree.class))

decay.class.text<-legend_decay_class_remapped[decay.class]
unique(decay.class.text)
tree.class.text<-c("Pine","Spruce","Birch", "Spruce (Tarasov)")[tree.class]
unique(tree.class.text)
decay.tree.class.text<-interaction(decay.class.text, tree.class.text)
unique(decay.tree.class.text)
class.levels.text<-(unique(as.factor(decay.tree.class.text)))

decay.tree.class.text[1]
decay.tree.class[1]


#count the individuals
count<-count(decay.tree.class.text)
palette<-c(pal_jama("default")(1),pal_npg("nrc")(9))

length(decay.tree.class.text)
sum(count$freq)

waffle_vec<-count$freq
names(waffle_vec)<-paste(count$x," ",waffle_vec)
#png(filename = "Class_count_reclassed.png", height=1500, width=3600, res=350)
png(filename = "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Class_count_reclassed.png", height=1500, width=3600, res=350)
waffle(waffle_vec, rows=32, col=palette)
dev.off()

sum(waffle_vec)

#create the table to hold the results
HSY_table_remapped<-mat.or.vec((dim(resampled_posteriors)[2]-2), length(class.levels.text))
rownames(HSY_table_remapped)<-colnames(resampled_posteriors)[1:(length(colnames(resampled_posteriors))-2)]
colnames(HSY_table_remapped)<-class.levels.text

for(j in 1:length(class.levels.text)){
  RMSE.table.subclass<-RMSE_table_classed[decay.tree.class.text==class.levels.text[j],]
  last<-dim(RMSE.table.subclass)[2]

  if(dim(RMSE.table.subclass)[1]>1){
    RMSE.mean.subclass<-colMeans(RMSE.table.subclass[,-last], na.rm=T)
    quantiles.subclass<-quantile(RMSE.mean.subclass, c(0.05,0.95), na.rm=T)
    which.subclass<-which(RMSE.mean.subclass<quantiles.subclass[1])
    bin1<-resampled_posteriors.long[RMSE.mean.subclass<quantiles.subclass[1],]
    bin2<-resampled_posteriors.long[!RMSE.mean.subclass<quantiles.subclass[1],]

    for(i in 1:(dim(resampled_posteriors)[2]-2)){
      HSY_table_remapped[i,j]<-ks.test(bin1[,i],bin2[,i])$statistic
    }
  }else{
    HSY_table_remapped[,j]<-rep(NA, (dim(resampled_posteriors)[2]-2))
  }

}


colMeans(RMSE.table.subclass, na.rm=T)


colnames(HSY_table_remapped)<-class.levels.text


rowMeans(HSY_table_remapped)
round(sort(rowMeans(HSY_table_remapped)),2)

HSY_table_remapped.scaled<-HSY_table_remapped
for(i in 1:dim(HSY_table_remapped)[2]){
  if(!all(is.na(HSY_table_remapped[,i]))){
    HSY_table_remapped.scaled[,i]<-range01(as.numeric(HSY_table_remapped[,i]))
  }else{
    HSY_table_remapped.scaled[,i]<-rep(NA, dim(HSY_table_remapped.scaled)[1])
  }
}

HSY_table_remapped.ranked<-HSY_table_remapped
for(i in 1:dim(HSY_table_remapped)[2]){
  if(!all(is.na(HSY_table_remapped[,i]))){
    HSY_table_remapped.ranked[,i]<-rank(as.numeric(HSY_table_remapped[,i]))
  }else{
    HSY_table_remapped.ranked[,i]<-rep(NA, dim(HSY_table_remapped.ranked)[1])
  }
}

rowMeans(HSY_table_remapped.ranked)

order<-order(rowMeans(HSY_table_remapped_ave))
rowMeans(HSY_table_remapped_ave)[order]

#png(filename = "HSY_table_remapped_scaled.png", height=2800, width=3500, res=350)
png(filename = "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/HSY_table_remapped_scaled.png", height=2800, width=3500, res=350)
par(mar=c(15,5,2,2))
color2D.matplot(HSY_table_remapped.scaled, axes=F, xlab="", ylab="", show.values = 4, vcex = 1, extremes = c("white", "red"))
axis(1, at = seq_len(ncol(HSY_table_remapped.scaled)) - 0.5, las=2,
     labels = colnames(HSY_table_remapped.scaled), tick = FALSE, cex.axis = 1)
axis(2, at = seq_len(nrow(HSY_table_remapped.scaled)) -0.5,
     labels = rev(parameters_list), tick = FALSE, las = 1, cex.axis = 1)
dev.off()





HSY_table_remapped.col=HSY_table_remapped
HSY_table_remapped.col[is.na(HSY_table_remapped.col)]<-0

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

colnames(HSY_table_remapped)
logs<-c(1,5,9)
snags<-c(2,6,8)
semisnag<-c(3,4,7)

HSY_table_remapped_ave<-mat.or.vec(dim(HSY_table_remapped)[1], 3)
HSY_table_remapped_ave[,1]<-rowMeans(HSY_table_remapped[,logs])
HSY_table_remapped_ave[,2]<-rowMeans(HSY_table_remapped[,snags])
HSY_table_remapped_ave[,3]<-rowMeans(HSY_table_remapped[,semisnag])

colnames(HSY_table_remapped_ave)<-c("logs", "snags", "snags <1/3 broken")
rownames(HSY_table_remapped_ave)<-rownames(HSY_table_remapped)

png(filename = "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/HSY_table_remapped_average.png", height=2000, width=3500, res=350)
barplot(t(HSY_table_remapped_ave), beside=T, col=c("chocolate1", "cadetblue1","chartreuse1"), las=2, ylab="Sensitivity", ylim=c(0,0.8), xpd=FALSE, names.arg=parameters_list)
legend("topright", c("logs", "snags >1/3 broken", "snags <1/3 broken"), bty="n", col=c("chocolate1", "cadetblue1","chartreuse1"), pch=16)
box()
dev.off()



colnames(HSY_table_remapped)
HSY_reshuffle_vec<-c(1,2,3,5,6,4,9,8,7,10)
HSY_table_remapped_reshuffled<-HSY_table_remapped[,HSY_reshuffle_vec]

cellcol<-map2color(HSY_table_remapped.col, rev(heat.colors(length(HSY_table_remapped.col))))
#png(filename = "HSY_table_remapped.png", height=2800, width=3500, res=350)
png(filename = "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/HSY_table_remapped.png", height=2800, width=3500, res=350)
par(mar=c(15,5,2,2))
color2D.matplot(HSY_table_remapped_reshuffled, axes=F, xlab="", ylab="", show.values = 3, vcex = 1, extremes = c("white", "red"))#, cellcolors = cellcol)
axis(1, at = seq_len(ncol(HSY_table_remapped_reshuffled)) - 0.5, las=2,
     labels = colnames(HSY_table_remapped_reshuffled), tick = FALSE, cex.axis = 1)
axis(2, at = seq_len(nrow(HSY_table_remapped_reshuffled)) -0.5,
     labels = rev(parameters_list), tick = FALSE, las = 1, cex.axis = 1)
dev.off()



#png(filename = "HSY_table_remapped_presentation.png", height=2000, width=3500, res=350)
png(filename = "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/HSY_table_remapped_presentation.png", height=2000, width=3500, res=350)
par(mar=c(14,5,2,2))
color2D.matplot(HSY_table_remapped, axes=F, xlab="", ylab="", show.values = 3, vcex = 1, extremes = c("white", "red"))#, cellcolors = cellcol)
axis(1, at = seq_len(ncol(HSY_table_remapped)) - 0.5, las=2,
     labels = colnames(HSY_table_remapped), tick = FALSE, cex.axis = 1)
axis(2, at = seq_len(nrow(HSY_table_remapped)) -0.5,
     labels = rev(parameters_list), tick = FALSE, las = 1, cex.axis = 1)
dev.off()


HSY_table_remapped.col=HSY_table_remapped.ranked
HSY_table_remapped.col[is.na(HSY_table_remapped.ranked)]<-0
cols<-unique(as.vector(HSY_table_remapped.col))

#cellcol<-map2color(HSY_table_remapped.col, rev(heat.colors(length(cols))))
cellcol<-map2color(HSY_table_remapped.col, c("white",rev(brewer.pal(length(cols), "RdYlBu"))))
cellcol<-map2color(HSY_table_remapped.col,c("white",(viridis(length(cols)))))
#png(filename = "HSY_table_remapped_ranked.png", height=2800, width=2500, res=350)
png(filename = "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/HSY_table_remapped_ranked.png", height=2800, width=2500, res=350)
par(mar=c(15,5,2,2))
color2D.matplot(HSY_table_remapped.ranked, axes=F, xlab="", ylab="", show.values = 1, vcex = 1, extremes = c("white", "red"))#, cellcolors = cellcol)
axis(1, at = seq_len(ncol(HSY_table_remapped.ranked)) - 0.5, las=2,
     labels = colnames(HSY_table_remapped.ranked), tick = FALSE, cex.axis = 1)
axis(2, at = seq_len(nrow(HSY_table_remapped.ranked)) -0.5,
     labels = rev(parameters_list), tick = FALSE, las = 1, cex.axis = 1)
dev.off()





rank(as.numeric(levels(decay.tree.class)))


decay.tree.class.remapped<-seq(1:16)[decay.tree.class]
decay.tree.class.numeric<-unique(as.numeric(decay.tree.class))
unique(decay.tree.class)
unique.decay.tree.class.text=unique(decay.tree.class.text)
levels(unique.decay.tree.class.text)[16]<-"Tarasov et al."



levels.decay.tree.class.text=levels(decay.tree.class.text)

decay.tree.class.palette<-c(pal_igv("default")(10))
decay.tree.class.palette.alpha<-add.alpha(decay.tree.class.palette,0.6)
#color palette for the new classifier
classes_palette<-c(brewer.pal(5, "Greens")[2:5],
                   brewer.pal(5, "Blues")[2:5],
                   brewer.pal(5, "Reds")[2:5],
                   rep("darkorange",4))

######local fit
calib_data_local<-list(N=dim(calibration_data_filtered)[1]+dim(Tarasov_data)[1],
                 time_gone=c(2003-calibration_data_filtered$Year.died, Tarasov_data$t),
                 mass_left=c(calibration_data_filtered$cmass.omass, Tarasov_data$y),
                 u0=u0_calc,
                 max_sd=max(sd_vector, na.rm=T),
                 tmax=c(calibration_data_filtered$Stem.diameter, Tarasov_data$tmax),
                 scaled_d=rangeMinMax(c(calibration_data_filtered$Original.wood.density, Tarasov_data$dens), 0.8,1.2),
                 decay_tree_class=as.numeric(decay.tree.class),
                 tree_class=as.numeric(tree.class),
                 delay_class=as.numeric(decay.class==2 | decay.class==1))

#running the Stan model and creating the calibration objects
start.time <- Sys.time()
fit_local <- stan(file = 'dec_model_Q_specific.stan', data = calib_data_local,  chains = 4, iter = iterations, cores=4,  control = list(adapt_delta = 0.9), thin=2)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

posteriors_local<-as.data.frame(fit_local )
names_local <-colnames(posteriors_local )


parameters_list_local<-c(expression(beta),
                         expression(eta[11]),
                         expression(paste(q[0]," Birch")),
                         expression(paste(q[0]," Spruce")),
                         expression(paste(q[0]," Pine")),
                         expression(paste(q[0]," Tarasov et al.")),
                         expression(e[0]),
                         expression(f[c]),
                         "Delay snags",
                         "Delay logs",
                         expression(paste(T[max]," error log.Pine ")),
                         expression(paste(T[max]," error snag.Pine")),
                         expression(paste(T[max]," error  snag, < 1/3 broken.Pine")),
                         expression(paste(T[max]," error snag, < 1/3 broken.Spruce")),
                         expression(paste(T[max]," error log.Spruce")),
                         expression(paste(T[max]," error snag.Spruce")),
                         expression(paste(T[max]," error snag, < 1/3 broken.Birch")),
                         expression(paste(T[max]," error snag.Birch")),
                         expression(paste(T[max]," error log.Birch")),
                         expression(paste(T[max]," error Tarasov et al.")),
                         expression(paste(u[0]," error log.Pine ")),
                         expression(paste(u[0]," error snag.Pine")),
                         expression(paste(u[0]," error  snag, < 1/3 broken.Pine")),
                         expression(paste(u[0]," error snag, < 1/3 broken.Spruce")),
                         expression(paste(u[0]," error log.Spruce")),
                         expression(paste(u[0]," error snag.Spruce")),
                         expression(paste(u[0]," error snag, < 1/3 broken.Birch")),
                         expression(paste(u[0]," error snag.Birch")),
                         expression(paste(u[0]," error log.Birch")),
                         expression(paste(u[0]," error Tarasov et al.")))




#diagnostics
print(fit_local)

available_mcmc(pattern = "_nuts_")

lp_cp <- log_posterior(fit_local)
head(lp_cp)

np_cp <- nuts_params(fit_local)
head(np_cp)

chain_palette<-wes_palette("Darjeeling2", 4)

#trace of mixing
posterior_cp <- as.array(fit_local)
#png("traceplot.png", height = 2600, width = 2000, res=300)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/traceplot.png", height = 2600, width = 2000, res=300)
par(mfrow=c(3,2))
plot(posterior_cp[,1,"beta"], type="l", col=chain_palette[1], ylab=parameters_list_local[1], ylim=range(posterior_cp[,1,"beta"])*c(0.98,1.02))
for(i in 2:4){lines(posterior_cp[,i,"beta"], col=chain_palette[i])}
legend("bottomleft", "(A)", bty="n")
legend("topright",c("Chain 1", "Chain 2", "Chain 3", "Chain 4"), col=chain_palette, lty=1, bty="n")
plot(posterior_cp[,1,"eta_11"], type="l", col=chain_palette[1], ylab=parameters_list_local[2], ylim=range(posterior_cp[,1,"eta_11"])*c(0.98,1.02))
for(i in 2:4){lines(posterior_cp[,i,"eta_11"], col=chain_palette[i])}
legend("bottomleft", "(B)", bty="n")
plot(posterior_cp[,1,"q0_p"], type="l", col=chain_palette[1], ylab=parameters_list_local[5], ylim=range(posterior_cp[,1,"q0_p"])*c(0.98,1.02))
legend("bottomleft", "(C)", bty="n")
for(i in 2:4){lines(posterior_cp[,i,"q0_p"], col=chain_palette[i])}
plot(posterior_cp[,1,"e0"], type="l", col=chain_palette[1], ylab=parameters_list_local[7], ylim=range(posterior_cp[,1,"e0"])*c(0.98,1.02))
legend("bottomleft", "(D)", bty="n")
for(i in 2:4){lines(posterior_cp[,i,"e0"], col=chain_palette[i])}
plot(posterior_cp[,1,"tmax_error1"], type="l", col=chain_palette[1], ylab=parameters_list_local[11], ylim=range(posterior_cp[,1,"tmax_error1"])*c(0.98,1.02))
legend("bottomleft", "(E)", bty="n")
for(i in 2:4){lines(posterior_cp[,i,"tmax_error1"], col=chain_palette[i])}
plot(posterior_cp[,1,"u0_error1"], type="l", col=chain_palette[1], ylab=parameters_list_local[21], ylim=range(posterior_cp[,1,"u0_error1"])*c(0.98,1.02))
legend("bottomleft", "(F)", bty="n")
for(i in 2:4){lines(posterior_cp[,i,"u0_error1"], col=chain_palette[i])}
dev.off()



#png("traceplot_presentation.png", height = 1400, width = 2800, res=300)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/traceplot_presentation.png", height = 1400, width = 2800, res=300)
par(mfrow=c(1,2))
plot(posterior_cp[,1,"beta"], type="l", col=chain_palette[1], ylab=parameters_list_local[1], ylim=range(posterior_cp[,1,"beta"])*c(0.98,1.02))
for(i in 2:4){lines(posterior_cp[,i,"beta"], col=chain_palette[i])}
legend("bottomleft", "(A)", bty="n")
legend("topright",c("Chain 1", "Chain 2", "Chain 3", "Chain 4"), col=chain_palette, lty=1, bty="n")
plot(posterior_cp[,1,"eta_11"], type="l", col=chain_palette[1], ylab=parameters_list_local[2], ylim=range(posterior_cp[,1,"eta_11"])*c(0.98,1.02))
for(i in 2:4){lines(posterior_cp[,i,"eta_11"], col=chain_palette[i])}
legend("bottomleft", "(B)", bty="n")
dev.off()

#energy
#png("MCMCenergy.png", height = 1800, width = 2000, res=300)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/MCMCenergy.png", height = 1800, width = 2000, res=300)
color_scheme_set("red")
mcmc_nuts_energy(np_cp)
dev.off()

#rhats
rhats <- rhat(fit_local)
print(rhats)

str(rhats)
color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats)

#png("Rhat.png", height = 2000, width = 2000, res=320)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Rhat.png", height = 2000, width = 2000, res=320)
par(mar=c(15,5,2,2))
barplot(rhats[c(-31,-32)], ylim=c(0.9995,1.1),  xpd=FALSE, las=2, names.arg=parameters_list_local, col="cadetblue")
mtext(expression(hat(R)), 2, padj=-3.7)
abline(h=1.05, lty=2, col="grey")
abline(h=1.1, lty=3, col="grey")
box()
dev.off()


ratios_cp <- neff_ratio(fit_local)
print(ratios_cp)

#autocorrelations

acbeta<-list()
for(i in 1:4){acbeta[[i]]<-acf(posterior_cp[,i,"beta"])}
aceta_11<-list()
for(i in 1:4){aceta_11[[i]]<-acf(posterior_cp[,i,"eta_11"])}
acq0_p<-list()
for(i in 1:4){acq0_p[[i]]<-acf(posterior_cp[,i,"q0_p"])}
ace_0<-list()
for(i in 1:4){ace_0[[i]]<-acf(posterior_cp[,i,"e0"])}
actmax_error1<-list()
for(i in 1:4){actmax_error1[[i]]<-acf(posterior_cp[,i,"tmax_error1"])}
acu0_error1<-list()
for(i in 1:4){acu0_error1[[i]]<-acf(posterior_cp[,i,"u0_error1"])}

#png("ACF.png", height = 2600, width = 2000, res=300)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/ACF.png", height = 2600, width = 2000, res=300)
par(mfrow=c(3,2))

lwd=2
plot(acbeta[[1]]$acf, type="l", ylim=c(-0.3,1), col=chain_palette[1], main=parameters_list_local[1], ylab="Autocorrelation", xlab="Lag", lwd=lwd)
lines(acbeta[[2]]$acf, lty=1, col=chain_palette[2], lwd=lwd)
lines(acbeta[[3]]$acf, lty=1, col=chain_palette[3], lwd=lwd)
lines(acbeta[[4]]$acf, lty=1, col=chain_palette[4], lwd=lwd)
legend("topright",c("Chain 1", "Chain 2", "Chain 3", "Chain 4"), col=chain_palette, lty=1, bty="n")
abline(h=0, col="grey")
legend("bottomright", "(A)", bty="n")

plot(aceta_11[[1]]$acf, type="l", ylim=c(-0.3,1), col=chain_palette[1], main=parameters_list_local[2], ylab="Autocorrelation", xlab="Lag", lwd=lwd)
lines(aceta_11[[2]]$acf, lty=1, col=chain_palette[2], lwd=lwd)
lines(aceta_11[[3]]$acf, lty=1, col=chain_palette[3], lwd=lwd)
lines(aceta_11[[4]]$acf, lty=1, col=chain_palette[4], lwd=lwd)
abline(h=0, col="grey")
legend("bottomright", "(B)", bty="n")

plot(acq0_p[[1]]$acf, type="l", ylim=c(-0.3,1), col=chain_palette[1], main=parameters_list_local[5], ylab="Autocorrelation", xlab="Lag", lwd=lwd)
lines(acq0_p[[2]]$acf, lty=1, col=chain_palette[2], lwd=lwd)
lines(acq0_p[[3]]$acf, lty=1, col=chain_palette[3], lwd=lwd)
lines(acq0_p[[4]]$acf, lty=1, col=chain_palette[4], lwd=lwd)
abline(h=0, col="grey")
legend("bottomright", "(C)", bty="n")

plot(ace_0[[1]]$acf, type="l", ylim=c(-0.3,1), col=chain_palette[1], main=parameters_list_local[7], ylab="Autocorrelation", xlab="Lag", lwd=lwd)
lines(ace_0[[2]]$acf, lty=1, col=chain_palette[2], lwd=lwd)
lines(ace_0[[3]]$acf, lty=1, col=chain_palette[3], lwd=lwd)
lines(ace_0[[4]]$acf, lty=1, col=chain_palette[4], lwd=lwd)
abline(h=0, col="grey")
legend("bottomright", "(D)", bty="n")

plot(actmax_error1[[1]]$acf, type="l", ylim=c(-0.3,1), col=chain_palette[1], main=parameters_list_local[11], ylab="Autocorrelation", xlab="Lag", lwd=lwd)
lines(actmax_error1[[2]]$acf, lty=1, col=chain_palette[2], lwd=lwd)
lines(actmax_error1[[3]]$acf, lty=1, col=chain_palette[3], lwd=lwd)
lines(actmax_error1[[4]]$acf, lty=1, col=chain_palette[4], lwd=lwd)
abline(h=0, col="grey")
legend("bottomright", "(E)", bty="n")

plot(acu0_error1[[1]]$acf, type="l", ylim=c(-0.3,1), col=chain_palette[1], main=parameters_list_local[21], ylab="Autocorrelation", xlab="Lag", lwd=lwd)
lines(acu0_error1[[2]]$acf, lty=1, col=chain_palette[2], lwd=lwd)
lines(acu0_error1[[3]]$acf, lty=1, col=chain_palette[3], lwd=lwd)
lines(acu0_error1[[4]]$acf, lty=1, col=chain_palette[4], lwd=lwd)
abline(h=0, col="grey")
legend("bottomright", "(F)", bty="n")

dev.off()




###### select posteriors where zeta>1
posteriors_local

zeta_vec<- (1-posteriors_local$e0)/(posteriors_local$beta*posteriors_local$eta_11*posteriors_local$e0)
zeta_select<-which(!zeta_vec<1)
length(zeta_select)

plot(density(zeta_vec))
plot(density(zeta_vec[zeta_select]))


posteriors_local_select<-posteriors_local[zeta_select,]



plot(density(posteriors_local_select$beta))
plot(density(posteriors_local$beta))
plot(density(posteriors_local_select$eta_11))
plot(density(posteriors_local$eta_11))
plot(density(posteriors_local_select$q0_s))





## plot parameters posteriors and priors

#define the priors
prior_realizations=30000
beta        = rnorm(prior_realizations, 7,0.7*0.15);
eta_11      = rnorm(prior_realizations, 0.36,0.36*0.15);
q0          = rnorm(prior_realizations, 1.101848,0.1275926);
e0          = rtruncnorm(prior_realizations, mean=0.3, sd=0.25, a=0, b=0.6);
fc          = rnorm(prior_realizations, 0.5,0.5*0.15);
delay       = runif(prior_realizations, 0,5);
tmax_error  = rnorm(prior_realizations, 1,1*0.15);
u0_error    = rnorm(prior_realizations, 1,1*0.15);

priors<-list(beta,
             eta_11,
             q0,
             e0,
             fc,
             delay,
             tmax_error,
             u0_error)


posteriors_local_table<-mat.or.vec(dim(posteriors_local)[2], 3)
rownames(posteriors_local_table)<-colnames(posteriors_local)

#png("ACF.png", height = 2600, width = 2000, res=300)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Parameters_local_1.png", height = 2600, width = 2000, res=300)
par(mfrow=c(2,2))
for(i in c(1,2,7,8)){
  dist<-posteriors_local[[i]]

  if(i==1 | i==2){
    prior<-priors[[i]]
    } else if (i==7) {
      prior<-priors[[4]]
    } else if (i==8) {
      prior<-priors[[5]]
    }

  posteriors_local_table[i, 1]<-mean(posteriors_local[[i]])
  posteriors_local_table[i, 2:3]<-quantile(posteriors_local[[i]], c(0.025, 0.975))
  plot(density(as.numeric(dist)), main=parameters_list_local[i], col="darkgreen", xlim=range(c(density(as.numeric(dist))$x), density(prior)$x), ylim=range(c(density(as.numeric(dist))$y), density(priors[[i]])$y))
  polygon(density(prior),lty=2, col=add.alpha("grey",0.8))
  polygon(density(as.numeric(dist)), col=add.alpha("darkgreen",0.3))
  legend("topleft", paste("(",c("A","B","","","","","C","D")[i],")"), bty="n")

}
dev.off()



#png("Parameters_local_2.png", width=4000, height=4000, res=380)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Parameters_local_2.png", width=4000, height=4000, res=380)
par(mfrow=c(2,2))
colnames(posteriors_local)
#plot q0
dist<-priors[[3]]
range_vector<-unlist(c(posteriors_local[,3:6]))
plot(density(as.numeric(dist)), main=parameters_list[3], xlim=range(c(range_vector, density(priors[[3]])$x)),
     ylim=c(0,15), col=NA, xlab="")
polygon(density(priors[[3]]),lty=2, col=add.alpha("lightgrey",0.8))

for(i in 1:4){
  dist<-posteriors_local[,2+i]
  posteriors_local_table[2+i, 1]<-mean(posteriors_local[[2+i]])
  posteriors_local_table[2+i, 2:3]<-quantile(posteriors_local[[2+i]], c(0.025, 0.975))
  dens=density(as.numeric(dist))
  which_max=which.max(dens$y)
  polygon(dens, col=add.alpha(tree_palette[i],0.2))
  text(dens$x[which_max], dens$y[which_max]*1.03, parameters_list_local[3:6][i], col=tree_palette[i],cex=0.9)
  }
legend("topleft", "(A)", bty="n")
legend("topright", c("prior distribution"), pch=16, col=c(add.alpha("lightgrey",0.8)), bty="n", cex=1.1)

#plot delays
delay_palette<-c("salmon", "saddlebrown")

dist<-priors[[6]]
range_vector<-unlist(c(posteriors_local[,9:10]))
plot(density(as.numeric(dist)), main=expression(Delta), xlim=range(c(density(as.numeric(range_vector))$x), density(priors[[6]])$x),
     ylim=c(0,7.3), col=NA, xlab="")
polygon(density(priors[[6]]),lty=2, col=add.alpha("lightgrey",0.8))

for(i in 1:2){
  dist<-posteriors_local[,8+i]
  posteriors_local_table[8+i, 1]<-mean(posteriors_local[[8+i]])
  posteriors_local_table[8+i, 2:3]<-quantile(posteriors_local[[8+i]], c(0.025, 0.975))
  dens=density(as.numeric(dist))
  which_max=which.max(dens$y)
  polygon(dens, col=add.alpha(delay_palette[i],0.2))
  text(dens$x[which_max], dens$y[which_max]*1.03, c("snag","log")[i], col=delay_palette[i],cex=0.9)
}
legend("topleft", "(B)", bty="n")


#plot tmax
dist<-priors[[7]]
range_vector<-unlist(c(posteriors_local[,11:20]))
plot(density(as.numeric(dist)), main=parameters_list[7], xlim=range(c(density(as.numeric(range_vector))$x), density(priors[[7]])$x),
     ylim=c(0, max(c(density(as.numeric(range_vector))$y), density(priors[[7]])$y))*1.8, col=NA, xlab="")
polygon(density(priors[[7]]),lty=2, col=add.alpha("lightgrey",0.8))

for(i in 1:length(decay.tree.class.palette)){
  dist<-posteriors_local[,10+i]
  posteriors_local_table[10+i, 1]<-mean(posteriors_local[[10+i]])
  posteriors_local_table[10+i, 2:3]<-quantile(posteriors_local[[10+i]], c(0.025, 0.975))
  dens=density(as.numeric(dist))
  which_max=which.max(dens$y)
  polygon(dens, col=add.alpha(decay.tree.class.palette[i],0.2))
  text(dens$x[which_max], dens$y[which_max]*1.03, unique.decay.tree.class.text[i], col=decay.tree.class.palette[i],cex=0.9)

}

legend("topleft", "(C)", bty="n")


#plot u0
dist<-priors[[8]]
range_vector<-unlist(c(posteriors_local[,21:30]))
plot(density(as.numeric(dist)), main=parameters_list[8], xlim=range(c(density(as.numeric(range_vector))$x), density(priors[[8]])$x),
     ylim=range(c(density(as.numeric(range_vector))$y), density(priors[[8]])$y)*1.4, col=NA, xlab="")
polygon(density(priors[[8]]),lty=2, col=add.alpha("lightgrey",0.8))

for(i in 1:length(decay.tree.class.palette)){
  dist<-posteriors_local[,20+i]
  posteriors_local_table[20+i, 1]<-mean(posteriors_local[[20+i]])
  posteriors_local_table[20+i, 2:3]<-quantile(posteriors_local[[20+i]], c(0.025, 0.975))
  dens=density(as.numeric(dist))
  which_max=which.max(dens$y)
  polygon(dens, col=add.alpha(decay.tree.class.palette[i],0.2))
  text(dens$x[which_max], dens$y[which_max]*1.03, unique.decay.tree.class.text[i], col=decay.tree.class.palette[i],cex=0.9)
}

legend("topleft", "(D)", bty="n")


dev.off()


#png("Parameters_local_1_presentation.png", width=4500, height=2500, res=350)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Parameters_local_1_presentation.png", width=4000, height=2500, res=350)
par(mfrow=c(2,2))
for(i in c(1,2,7,8)){
  dist<-posteriors_local[[i]]

  if(i==1 | i==2){
    prior<-priors[[i]]
  } else if (i==7) {
    prior<-priors[[4]]
  } else if (i==8) {
    prior<-priors[[5]]
  }

  posteriors_local_table[i, 1]<-mean(posteriors_local[[i]])
  posteriors_local_table[i, 2:3]<-quantile(posteriors_local[[i]], c(0.025, 0.975))
  plot(density(as.numeric(dist)), main=parameters_list_local[i], col="darkgreen", xlim=range(c(density(as.numeric(dist))$x), density(prior)$x), ylim=range(c(density(as.numeric(dist))$y), density(priors[[i]])$y), xlab="")
  polygon(density(prior),lty=2, col=add.alpha("lightgrey",0.8))
  polygon(density(as.numeric(dist)), col=add.alpha("darkgreen",0.3))
  legend("topleft", paste("(",c("A","B","","","","","C","D")[i],")"), bty="n")
  
  if(i==1){
    legend("topright", c("prior distribution", "posterior distribution"), pch=16, col=c(add.alpha("lightgrey",0.8), add.alpha("darkgreen",0.3)), bty="n", cex=1.1)
  }
  
  
}
dev.off()


#png("Parameters_local_presentation.png", width=4500, height=2500, res=350)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Parameters_local_presentation.png", width=3800, height=2200, res=350)
par(mfrow=c(2,2))
colnames(posteriors_local)
#plot q0
dist<-priors[[3]]
range_vector<-unlist(c(posteriors_local[,3:6]))
plot(density(as.numeric(dist)), main=parameters_list[3], xlim=range(c(range_vector, density(priors[[3]])$x)),
     ylim=c(0,17), col=NA)
polygon(density(priors[[3]]),lty=2, col=add.alpha("lightgrey",0.8))

for(i in 1:4){
  dist<-posteriors_local[,2+i]
  posteriors_local_table[2+i, 1]<-mean(posteriors_local[[2+i]])
  posteriors_local_table[2+i, 2:3]<-quantile(posteriors_local[[2+i]], c(0.025, 0.975))
  dens=density(as.numeric(dist))
  which_max=which.max(dens$y)
  polygon(dens, col=add.alpha(tree_palette[i],0.2))
  text(dens$x[which_max], dens$y[which_max]*1.03, parameters_list_local[3:6][i], col=tree_palette[i],cex=0.9)
}

#plot delays
delay_palette<-c("salmon", "saddlebrown")

dist<-priors[[6]]
range_vector<-unlist(c(posteriors_local[,9:10]))
plot(density(as.numeric(dist)), main=parameters_list[6], xlim=range(c(density(as.numeric(range_vector))$x), density(priors[[6]])$x),
     ylim=c(0,7.3), col=NA)
polygon(density(priors[[6]]),lty=2, col=add.alpha("lightgrey",0.8))

for(i in 1:2){
  dist<-posteriors_local[,8+i]
  posteriors_local_table[8+i, 1]<-mean(posteriors_local[[8+i]])
  posteriors_local_table[8+i, 2:3]<-quantile(posteriors_local[[8+i]], c(0.025, 0.975))
  dens=density(as.numeric(dist))
  which_max=which.max(dens$y)
  polygon(dens, col=add.alpha(delay_palette[i],0.2))
  text(dens$x[which_max], dens$y[which_max]*1.03, c("snag","log")[i], col=delay_palette[i],cex=0.9)
}


#plot tmax
dist<-priors[[7]]
range_vector<-unlist(c(posteriors_local[,11:20]))
plot(density(as.numeric(dist)), main=parameters_list[7], xlim=range(c(density(as.numeric(range_vector))$x), density(priors[[7]])$x),
     ylim=c(0, max(c(density(as.numeric(range_vector))$y), density(priors[[7]])$y))*1.8, col=NA)
polygon(density(priors[[7]]),lty=2, col=add.alpha("lightgrey",0.8))

for(i in 1:length(decay.tree.class.palette)){
  dist<-posteriors_local[,10+i]
  posteriors_local_table[10+i, 1]<-mean(posteriors_local[[10+i]])
  posteriors_local_table[10+i, 2:3]<-quantile(posteriors_local[[10+i]], c(0.025, 0.975))
  dens=density(as.numeric(dist))
  which_max=which.max(dens$y)
  polygon(dens, col=add.alpha(decay.tree.class.palette[i],0.2))
  text(dens$x[which_max], dens$y[which_max]*1.03, unique.decay.tree.class.text[i], col=decay.tree.class.palette[i],cex=0.9)

}

#plot u0
dist<-priors[[8]]
range_vector<-unlist(c(posteriors_local[,21:30]))
plot(density(as.numeric(dist)), main=parameters_list[8], xlim=range(c(density(as.numeric(range_vector))$x), density(priors[[8]])$x),
     ylim=range(c(density(as.numeric(range_vector))$y), density(priors[[8]])$y)*1.6, col=NA)
polygon(density(priors[[8]]),lty=2, col=add.alpha("lightgrey",0.8))

for(i in 1:length(decay.tree.class.palette)){
  dist<-posteriors_local[,20+i]
  posteriors_local_table[20+i, 1]<-mean(posteriors_local[[20+i]])
  posteriors_local_table[20+i, 2:3]<-quantile(posteriors_local[[20+i]], c(0.025, 0.975))
  dens=density(as.numeric(dist))
  which_max=which.max(dens$y)
  polygon(dens, col=add.alpha(decay.tree.class.palette[i],0.2))
  text(dens$x[which_max], dens$y[which_max]*1.03, unique.decay.tree.class.text[i], col=decay.tree.class.palette[i],cex=0.9)
}


dev.off()








colnames(posteriors_local_table)<-c("Mean", "Min", "Max")


write_tab<-cbind(rownames(round(posteriors_local_table[1:30,],2)),round(posteriors_local_table[1:30,],2))

colnames(write_tab)<-c("Parameter","Mean", "Min", "Max")
write_tab[,1]<-c("$\\beta$",
                  "$\\eta_{11}$",
                  "$q_{0,b}$",
                  "$q_{0,s}$",
                  "$q_{0,p}$",
                  "$q_{0,t}$",
                  "$e_{0}$",
                  "$f_c$",
                  "$\\Delta_S$",
                  "$\\Delta_L$",
                  "$\\epsilon_{T_{max,1}}$",
                  "$\\epsilon_{T_{max,2}}$",
                  "$\\epsilon_{T_{max,3}}$",
                  "$\\epsilon_{T_{max,4}}$",
                  "$\\epsilon_{T_{max,5}}$",
                  "$\\epsilon_{T_{max,6}}$",
                  "$\\epsilon_{T_{max,7}}$",
                  "$\\epsilon_{T_{max,8}}$",
                  "$\\epsilon_{T_{max,9}}$",
                  "$\\epsilon_{T_{max,10}}$",
                  "$\\epsilon_{u_{0,1}}$",
                  "$\\epsilon_{u_{0,2}}$",
                  "$\\epsilon_{u_{0,3}}$",
                  "$\\epsilon_{u_{0,4}}$",
                  "$\\epsilon_{u_{0,5}}$",
                  "$\\epsilon_{u_{0,6}}$",
                  "$\\epsilon_{u_{0,7}}$",
                  "$\\epsilon_{u_{0,8}}$",
                  "$\\epsilon_{u_{0,9}}$",
                  "$\\epsilon_{u_{0,10}}$")
write_tab<-as.data.frame(write_tab)
write.csv(write_tab, "Local_parameters.csv",  row.names=FALSE)
#write.csv(write_tab, "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Local_parameters.csv",  row.names=FALSE)

write_tab_global<-write_tab[c(1,2,7,8),]
write_tab_local<-write_tab[c(3,4,5,6,9:30),]

write.csv(write_tab_global, "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Global_parameters.csv",  row.names=FALSE, quote = FALSE)
write.csv(write_tab_local, "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Local_parameters.csv",  row.names=FALSE, quote = FALSE)




#png("Parameters_local_ridgeplot_q0.png", width=2000, height=2000, res=350)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Parameters_local_ridgeplot_q0.png", width=2000, height=2000, res=350)
#q0
ridgeplotting_number<-c(priors[[3]])
ridgeplotting_names<-c(rep("Prior", length(priors[[8]])))
# ridgeplotting_names<-c()
# ridgeplotting_number<-c()
for(i in 1:4){
  ridgeplotting_number<-c(ridgeplotting_number,posteriors_local[,2+i])
  ridgeplotting_names<-c(ridgeplotting_names, as.character(rep(c("Birch","Spruce","Pine",  "Tarasov et al.")[i], length(posteriors_local[,2+i]))))
}
ridgeplotting_table<-as.data.frame(cbind(as.factor(ridgeplotting_names), (ridgeplotting_number)))
colnames(ridgeplotting_table)<-c("classes","value")
ridgeplotting_table$classes<-as.factor(ridgeplotting_names)
levels(ridgeplotting_table$classes)
ridgeplotting_table$classes<-factor(ridgeplotting_table$classes,levels(ridgeplotting_table$classes)[c(3,1,2,4,5)])

# tmax
ggplot(ridgeplotting_table, aes(x = value, y =classes, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 0.9) +
  #scale_fill_viridis(alpha=1, begin=0, end=1, option = "C") +
  labs(title =parameters_list[3]) +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )
dev.off()

#png("Parameters_local_ridgeplot_delay.png", width=2000, height=1500, res=350)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Parameters_local_ridgeplot_delay.png", width=2000, height=1500, res=350)
#delay
ridgeplotting_number<-c(priors[[6]])
ridgeplotting_names<-c(rep("Prior", length(priors[[8]])))
# ridgeplotting_names<-c()
# ridgeplotting_number<-c()
for(i in 1:2){
  ridgeplotting_number<-c(ridgeplotting_number,posteriors_local[,8+i])
  ridgeplotting_names<-c(ridgeplotting_names, as.character(rep(c("snag","log")[i], length(posteriors_local[,8+i]))))
}
ridgeplotting_table<-as.data.frame(cbind(as.factor(ridgeplotting_names), (ridgeplotting_number)))
colnames(ridgeplotting_table)<-c("classes","value")
ridgeplotting_table$classes<-as.factor(ridgeplotting_names)
levels(ridgeplotting_table$classes)
ridgeplotting_table$classes<-factor(ridgeplotting_table$classes,levels(ridgeplotting_table$classes)[c(2,1,3)])

# tmax
ggplot(ridgeplotting_table, aes(x = value, y =classes, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 0.9) +
  #scale_fill_viridis(alpha=1, begin=0, end=1, option = "C") +
  labs(title =parameters_list[6]) +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )
dev.off()


#png("Parameters_local_ridgeplot_tmax.png", width=2000, height=3000, res=350)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Parameters_local_ridgeplot_tmax.png", width=2000, height=3000, res=350)
ridgeplotting_number<-c(priors[[7]])
ridgeplotting_names<-c(rep("Prior", length(priors[[7]])))
# ridgeplotting_names<-c()
# ridgeplotting_number<-c()
for(i in 1:length(decay.tree.class.palette)){
  ridgeplotting_number<-c(ridgeplotting_number,posteriors_local[,10+i])
  ridgeplotting_names<-c(ridgeplotting_names, as.character(rep(unique.decay.tree.class.text[i], length(posteriors_local[,10+i]))))
 }
ridgeplotting_table<-as.data.frame(cbind(as.factor(ridgeplotting_names), (ridgeplotting_number)))
colnames(ridgeplotting_table)<-c("classes","value")
ridgeplotting_table$classes<-as.factor(ridgeplotting_names)
levels(ridgeplotting_table$classes)
ridgeplotting_table$classes<-factor(ridgeplotting_table$classes,levels(ridgeplotting_table$classes)[c(4,1,2,3,5,6,7,8,9,10,11)])

# tmax
  ggplot(ridgeplotting_table, aes(x = value, y =classes, fill = ..x..)) +
    geom_density_ridges_gradient(scale = 0.9) +
    #scale_fill_viridis(alpha=1, begin=0, end=1, option = "C") +
    labs(title =parameters_list[7]) +
    theme_ipsum() +
    theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)
    )
  dev.off()


 # png("Parameters_local_ridgeplot_u0.png", width=2000, height=3000, res=350)
  png("Parameters_local_ridgeplot_u0.png", width=2000, height=3000, res=350)
  ridgeplotting_number<-c(priors[[8]])
  ridgeplotting_names<-c(rep("Prior", length(priors[[8]])))
  #ridgeplotting_number<-c()
  #ridgeplotting_names<-c()
  for(i in 1:length(decay.tree.class.palette)){
    ridgeplotting_number<-c(ridgeplotting_number,posteriors_local[,20+i])
    ridgeplotting_names<-c(ridgeplotting_names, as.character(rep(unique.decay.tree.class.text[i], length(posteriors_local[,20+i]))))
  }
  length(as.character(rep(unique.decay.tree.class.text[i], length(posteriors_local[,20+i]))))

  ridgeplotting_table<-as.data.frame(cbind(as.factor(ridgeplotting_names), (ridgeplotting_number)))
  colnames(ridgeplotting_table)<-c("classes","value")
  ridgeplotting_table$classes<-as.factor(ridgeplotting_names)
  levels(ridgeplotting_table$classes)
  ridgeplotting_table$classes<-factor(ridgeplotting_table$classes,levels(ridgeplotting_table$classes)[c(4,1,2,3,5,6,7,8,9,10,11)])

  # u0
  ggplot(ridgeplotting_table, aes(x = value, y =classes, fill = ..x..)) +
    geom_density_ridges_gradient(scale = 0.9) +
    #scale_fill_viridis(alpha=1, begin=0, end=1, option = "C") +
    labs(title =parameters_list[8]) +
    theme_ipsum() +
    theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 8)
    )
  dev.off()



  
  
  
    
##################################################################################
##################  CALIBRATION BY CLASSES - SIMULATION ##########################
##################################################################################

# thinning the posteriors to plot the simulations
thinning=50
resampling_vector<-seq(from=1, to=dim(posteriors_local)[1], by=thinning)
resampled_posteriors_local<-posteriors_local[resampling_vector,]
dim(resampled_posteriors_local)

# plot the simulations
time_sim=500
time_simulation=seq(from=0, to=time_sim)
time_simulation2=time_simulation


#create the table to hold the posterior probability desnity for each of the data points
RMSE_table_local<-mat.or.vec(length(calib_data_local$tree_class), dim(resampled_posteriors_local)[1])
Inflection_table_local<-mat.or.vec(length(calib_data_local$tree_class), dim(resampled_posteriors_local)[1])

#create the vector to hold the steady states
SS_table_local<-mat.or.vec(length(calib_data_local$tree_class), dim(resampled_posteriors_local)[1])

mass_left_simulated_list_local<-list()
for(j in 1:length(calib_data_local$tree_class)){

  mass_left_simulated_local<-mat.or.vec(dim(resampled_posteriors_local)[1], time_sim+1)

  for( i in 1: dim(resampled_posteriors_local)[1]){

    zeta = (1-resampled_posteriors_local$e0[i])/(resampled_posteriors_local$beta[i]*resampled_posteriors_local$eta_11[i]*resampled_posteriors_local$e0[i]);


    # if statements for q0
    if(calib_data_local$tree_class[j]==3){
      q0  = resampled_posteriors_local$q0_b[i]; # if wood is birch
    }else if(calib_data_local$tree_class[j]==2){
      q0  = resampled_posteriors_local$q0_s[i]; # if wood is spruce
    }else if(calib_data_local$tree_class[j]==1){
      q0  = resampled_posteriors_local$q0_p[i]; # if wood is pine
    }else if(calib_data_local$tree_class[j]==4){
      q0  = resampled_posteriors_local$q0_t[i];} # if the data is from Tarasov et al.



    # if statements for u0 error
    if(calib_data_local$decay_tree_class[j]==3){
      alpha = resampled_posteriors_local$fc[i]*resampled_posteriors_local$beta[i]*resampled_posteriors_local$eta_11[i]*u0_calc[i]*resampled_posteriors_local$u0_error1[i]*q0^resampled_posteriors_local$beta[i]; # alpha is recalculated every time because of varying u0 in different sites
      tmax_error  = resampled_posteriors_local$tmax_error1[i]*calib_data_local$scaled_d[j];
      }else if(calib_data_local$decay_tree_class[j]==2){
      alpha = resampled_posteriors_local$fc[i]*resampled_posteriors_local$beta[i]*resampled_posteriors_local$eta_11[i]*u0_calc[i]*resampled_posteriors_local$u0_error2[i]*q0^resampled_posteriors_local$beta[i]; # alpha is recalculated every time because of varying u0 in different sites
      tmax_error  = resampled_posteriors_local$tmax_error2[i]*calib_data_local$scaled_d[j];
      }else if(calib_data_local$decay_tree_class[j]==1){
      alpha = resampled_posteriors_local$fc[i]*resampled_posteriors_local$beta[i]*resampled_posteriors_local$eta_11[i]*u0_calc[i]*resampled_posteriors_local$u0_error3[i]*q0^resampled_posteriors_local$beta[i]; # alpha is recalculated every time because of varying u0 in different sites
      tmax_error  = resampled_posteriors_local$tmax_error3[i]*calib_data_local$scaled_d[j];
      }else if(calib_data_local$decay_tree_class[j]==5){
      alpha = resampled_posteriors_local$fc[i]*resampled_posteriors_local$beta[i]*resampled_posteriors_local$eta_11[i]*u0_calc[i]*resampled_posteriors_local$u0_error4[i]*q0^resampled_posteriors_local$beta[i]; # alpha is recalculated every time because of varying u0 in different sites
      tmax_error  = resampled_posteriors_local$tmax_error4[i]*calib_data_local$scaled_d[j];
      }else if(calib_data_local$decay_tree_class[j]==7){
      alpha = resampled_posteriors_local$fc[i]*resampled_posteriors_local$beta[i]*resampled_posteriors_local$eta_11[i]*u0_calc[i]*resampled_posteriors_local$u0_error5[i]*q0^resampled_posteriors_local$beta[i]; # alpha is recalculated every time because of varying u0 in different sites
      tmax_error  = resampled_posteriors_local$tmax_error5[i]*calib_data_local$scaled_d[j];
      }else if(calib_data_local$decay_tree_class[i]==6){
      alpha = resampled_posteriors_local$fc[i]*resampled_posteriors_local$beta[i]*resampled_posteriors_local$eta_11[i]*u0_calc[i]*resampled_posteriors_local$u0_error6[i]*q0^resampled_posteriors_local$beta[i]; # alpha is recalculated every time because of varying u0 in different sites
      tmax_error  = resampled_posteriors_local$tmax_error6[i]*calib_data_local$scaled_d[j];
      }else if(calib_data_local$decay_tree_class[j]==9){
      alpha = resampled_posteriors_local$fc[i]*resampled_posteriors_local$beta[i]*resampled_posteriors_local$eta_11[i]*u0_calc[i]*resampled_posteriors_local$u0_error7[i]*q0^resampled_posteriors_local$beta[i]; # alpha is recalculated every time because of varying u0 in different sites
      tmax_error  = resampled_posteriors_local$tmax_error7[i]*calib_data_local$scaled_d[j];
      }else if(calib_data_local$decay_tree_class[j]==10){
      alpha = resampled_posteriors_local$fc[i]*resampled_posteriors_local$beta[i]*resampled_posteriors_local$eta_11[i]*u0_calc[i]*resampled_posteriors_local$u0_error8[i]*q0^resampled_posteriors_local$beta[i]; # alpha is recalculated every time because of varying u0 in different sites
      tmax_error  = resampled_posteriors_local$tmax_error8[i]*calib_data_local$scaled_d[j];
      }else if(calib_data_local$decay_tree_class[j]==11){
      alpha = resampled_posteriors_local$fc[i]*resampled_posteriors_local$beta[i]*resampled_posteriors_local$eta_11[i]*u0_calc[i]*resampled_posteriors_local$u0_error9[i]*q0^resampled_posteriors_local$beta[i]; # alpha is recalculated every time because of varying u0 in different sites
      tmax_error  = resampled_posteriors_local$tmax_error9[i]*calib_data_local$scaled_d[j];
      }else if(calib_data_local$decay_tree_class[j]==16){
      alpha = resampled_posteriors_local$fc[i]*resampled_posteriors_local$beta[i]*resampled_posteriors_local$eta_11[i]*u0_calc[i]*resampled_posteriors_local$u0_error10[i]*q0^resampled_posteriors_local$beta[i]; # alpha is recalculated every time because of varying u0 in different sites
      tmax_error  = resampled_posteriors_local$tmax_error10[i]*calib_data_local$scaled_d[j];}

    # if statements fordelay
    if(calib_data_local$delay_class[j]==1){
      time_simulation_eff=time_simulation-resampled_posteriors_local$delayS[i];
    }else if(calib_data_local$delay_class[j]==0){
      time_simulation_eff=time_simulation-resampled_posteriors_local$delayL[i];
    }

    time_simulation_eff[time_simulation_eff<0]=0

    #imålement the IF cycle for tmax with a vector
    selection_vector<-(time_simulation_eff<(calib_data$tmax[i]*tmax_error))

    #for all the times when time<tmax
    mass_left_simulated_local[i,selection_vector] <- ((2/(calib_data$tmax[i]*tmax_error))*(1/(alpha*(1-zeta)))*((1+alpha*time_simulation_eff[selection_vector])^(1-zeta)-
                                                                                               (1-(time_simulation_eff[selection_vector]/(calib_data$tmax[i]*tmax_error))))+
                                                  ((2/(calib_data$tmax[i]*tmax_error)^2)*(1/(alpha^2*(1-zeta)*(2-zeta)))*(1-(1+alpha*time_simulation_eff[selection_vector])^(2-zeta)))+
                                                  (1-(time_simulation_eff[selection_vector]/(calib_data$tmax[i]*tmax_error)))^2)
    #for all the times when time>=tmax
    mass_left_simulated_local[i,!selection_vector] <- (2/(calib_data$tmax[i]*tmax_error))*(1/(alpha*(1-zeta)))*(1+alpha*time_simulation_eff[!selection_vector])^(1-zeta)+
      ((2/((calib_data$tmax[i]*tmax_error)^2))*(1/(alpha^2*(1-zeta)*(2-zeta)))*(((1+alpha*(time_simulation_eff[!selection_vector]-(calib_data$tmax[i]*tmax_error)))^(2-zeta))-((1+alpha*time_simulation_eff[!selection_vector])^(2-zeta))))

    #calculate the RMSE for that point
    RMSE_table_local[j,i]<-rmse(mass_left_simulated_local[i, calib_data_local$time_gone[j]],calib_data_local$mass_left[j])

    #find if the simulation has one inflection point
    Inflection_table_local[j,i]<-which.min(mass_left_simulated_local[i,])==length(mass_left_simulated_local[i,])

    #calculate the steady state for each run
    tmax=calib_data$tmax[i]*tmax_error
    SS_table_local[j,i]<-((1/alpha)*(1/(zeta-1)))+(tmax/3)


  }

  mass_left_simulated_list_local[[j]]<-mass_left_simulated_local
}


plot(density(SS_table_local[j,]))



###test if there are any simulation with inflection point
test<-which(Inflection_table_local==0)
if(length(test)<0){print("Attention!!!! There are some simulations which present an inflection point before the end")}


#build the LIST of big tables
# attention: one element of the list is the big table for that particular tree class!!!!
mass_left_simulated_local_listoflists<-list()
for(i in 1:length(decay.tree.class.numeric)){
selection_vector_which<-which(calib_data_local$decay_tree_class==decay.tree.class.numeric[i])
mass_left_simulated_table_local1<-mass_left_simulated_list_local[[selection_vector_which[1]]]
  for(j in 1:length(calib_data_local$decay_tree_class)){
    mass_left_simulated_table_local1<-rbind(mass_left_simulated_table_local1,mass_left_simulated_list_local[[selection_vector_which[j]]])
    }
  mass_left_simulated_local_listoflists[[i]]<-mass_left_simulated_table_local1
}

str(mass_left_simulated_local_listoflists)
dim(mass_left_simulated_table_local1)



#plotting the range
#png("Range_plot_byclass.png", width=2600, height=4000, res=360)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Range_plot_byclass.png", width=2600, height=4000, res=360)
par(mfrow=c(5,2))
decay_pch.reclassed<-seq(1:16)
for(i in 1:10){
  selection_vector<-calib_data_local$decay_tree_class==decay.tree.class.numeric[i]
  mass_left_simulated_table_local_selected<-mass_left_simulated_local_listoflists[[i]]
  plot(2003-calibration_data_filtered$Year.died, calibration_data_filtered$cmass.omass, ylim=c(0,1.03), col="red", xlab="Years", ylab="Fraction of initial mass",
       xlim=c(0,120),  xaxs="i",yaxs="i", cex=0, main=unique.decay.tree.class.text[i])
  polygon( c(time_simulation, rev(time_simulation)), c(colMax(mass_left_simulated_table_local_selected),rev(colMin(mass_left_simulated_table_local_selected))),
           col=add.alpha(unique(classes_palette[calib_data_local$decay_tree_class[selection_vector]]),0.4),
           border = add.alpha(unique(classes_palette[calib_data_local$decay_tree_class[selection_vector]]),0.9))
  points(calib_data_local$time_gone[selection_vector], calib_data_local$mass_left[selection_vector],
         col=classes_palette[calib_data_local$decay_tree_class[selection_vector]],
         bg=add.alpha(classes_palette[calib_data_local$decay_tree_class[selection_vector]],0.3),
         pch=calib_data_local$tree_class[selection_vector])
  #legend("topright", c(legend_decay_class), pch=seq(1:8), bty="n")
  #legend("bottomleft", c("Pine","Spruce", "Birch", "Spruce (Tarasov et al.)"), pch=21, col=tree_palette, pt.bg=add.alpha(tree_palette,0.3), bty="n")
  abline(h=0, lty=2)
  legend("topright", paste("(",LETTERS[i],")"), bty="n")

  }
dev.off()






#plotting the range
#png("Range_plot_byclass_presentation.png", width=4600, height=2000, res=360)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Range_plot_byclass_presentation.png", width=4600, height=2000, res=360)
par(mfrow=c(2,5))
decay_pch.reclassed<-seq(1:16)
for(i in 1:10){
  selection_vector<-calib_data_local$decay_tree_class==decay.tree.class.numeric[i]
  mass_left_simulated_table_local_selected<-mass_left_simulated_local_listoflists[[i]]
  plot(2003-calibration_data_filtered$Year.died, calibration_data_filtered$cmass.omass, ylim=c(0,1.03), col="red", xlab="time", ylab="mass remaining",
       xlim=c(0,120),  xaxs="i",yaxs="i", cex=0, main=unique.decay.tree.class.text[i])
  polygon( c(time_simulation, rev(time_simulation)), c(colMax(mass_left_simulated_table_local_selected),rev(colMin(mass_left_simulated_table_local_selected))),
           col=add.alpha(unique(classes_palette[calib_data_local$decay_tree_class[selection_vector]]),0.4),
           border = add.alpha(unique(classes_palette[calib_data_local$decay_tree_class[selection_vector]]),0.9))
  points(calib_data_local$time_gone[selection_vector], calib_data_local$mass_left[selection_vector],
         col=classes_palette[calib_data_local$decay_tree_class[selection_vector]],
         bg=add.alpha(classes_palette[calib_data_local$decay_tree_class[selection_vector]],0.3),
         pch=calib_data_local$tree_class[selection_vector])
  #legend("topright", c(legend_decay_class), pch=seq(1:8), bty="n")
  #legend("bottomleft", c("Pine","Spruce", "Birch", "Spruce (Tarasov et al.)"), pch=21, col=tree_palette, pt.bg=add.alpha(tree_palette,0.3), bty="n")
  abline(h=0, lty=2)
  legend("topright", paste("(",LETTERS[i],")"), bty="n")

}
dev.off()



SS_table_plot<-mat.or.vec(10,1000)
SS_table_values<-mat.or.vec(10,3)
for(i in 1:length(decay.tree.class.numeric)){
  selection_vector<-calib_data_local$decay_tree_class==decay.tree.class.numeric[i]
  quant<-quantile(SS_table_local[selection_vector,], c(0.05, 0.95))
  mode<-mlv(SS_table_local[selection_vector,])
  SS_table_plot[i,]<-sample(SS_table_local[selection_vector,], 1000)
  SS_table_values[i,1]<-mode
  SS_table_values[i,2:3]<-quant
}
#png("SS_plot_byclass.png", width=2600, height=3000, res=360)
png("./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/SS_plot_byclass.png", width=2600, height=3000, res=360)
par(mar=c(15,5,2,2))
boxplot(t(SS_table_plot), names=unique.decay.tree.class.text, las=2, ylim=c(0, 1850),
        ylab= expression(paste("C (T", ha^-1, ")")), col=classes_palette[as.numeric(unique.decay.tree.class.text)])

dev.off()

SS_table_values<-round(SS_table_values,0)
SS_table_values<-cbind(as.character(unique.decay.tree.class.text),SS_table_values)
colnames(SS_table_values)<-c("Decayclass","Mode","Min", "Max")
SS_table_values[3,1]<-"snag <1/3 broken.Pine"
SS_table_values[4,1]<-"snag <1/3 broken.Spruce"
SS_table_values[7,1]<-"snag <1/3 broken.Birch"
SS_table_values[,1]<-latexTranslate(SS_table_values[,1])




SS_table_values<-rbind(SS_table_values,c("General",round(mlv(SS_table),0),round(min(SS_table),0),round(max(SS_table),0)))

(SS_table_values)

    write.csv(SS_table_values, "SteadyStates.csv", row.names = F)
    write.csv(SS_table_values, "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/SteadyStates.csv", row.names = F, quote = FALSE)


mean(RMSE_table_local, na.rm=T)
mean(RMSE_table, na.rm=T)



names_vec<-levels(unique.decay.tree.class.text)[calib_data_local$decay_tree_class]

RMSE_local <-aggregate(rowMeans(RMSE_table_local), by=list(names_vec),
                       FUN=mean, na.rm=TRUE)
RMSE_global <-aggregate(rowMeans(RMSE_table), by=list(names_vec),
                       FUN=mean, na.rm=TRUE)

RMSE_aggregated<-rbind(RMSE_global[,2], RMSE_local[,2])

#png(filename = "RMSE.png", height =2500, width=3000, res=300)
png(filename = "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/RMSE.png", height =2500, width=3000, res=300)
par(mar=c(12,5,2,2))
barplot(RMSE_aggregated,beside=TRUE, names.arg=RMSE_global[,1], las=2, col=c("darkorange", "darkgreen"), angle=c(45,125), density=35, ylab="RMSE", ylim=c(0,0.45))
legend("topleft", c("General calibration", "Local calibration"),
       fill=c("darkorange", "darkgreen"), angle=c(45,125), density=35, bty="n")
box()

dev.off()

#png(filename = "RMSE_presentation.png", height =2500, width=4200, res=300)
png(filename = "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/RMSE_presentation.png", height =2500, width=4200, res=300)
par(mar=c(12,5,2,2))
barplot(RMSE_aggregated,beside=TRUE, names.arg=RMSE_global[,1], las=2, col=c("darkorange", "darkgreen"), angle=c(45,125), density=35, ylab="RMSE", ylim=c(0,0.45))
legend("topleft", c("General calibration", "Local calibration"),
       fill=c("darkorange", "darkgreen"), angle=c(45,125), density=35, bty="n")
box()

dev.off()

tab_predictions<-mat.or.vec(11,3)

tab_predictions[1,1]<-mean(mass_left_simulated[,50])
tab_predictions[2,1]<-mean(mass_left_simulated_local_listoflists[[1]][,50])
tab_predictions[3,1]<-mean(mass_left_simulated_local_listoflists[[2]][,50])
tab_predictions[4,1]<-mean(mass_left_simulated_local_listoflists[[3]][,50])
tab_predictions[5,1]<-mean(mass_left_simulated_local_listoflists[[4]][,50])
tab_predictions[6,1]<-mean(mass_left_simulated_local_listoflists[[5]][,50])
tab_predictions[7,1]<-mean(mass_left_simulated_local_listoflists[[6]][,50])
tab_predictions[8,1]<-mean(mass_left_simulated_local_listoflists[[7]][,50])
tab_predictions[9,1]<-mean(mass_left_simulated_local_listoflists[[8]][,50])
tab_predictions[10,1]<-mean(mass_left_simulated_local_listoflists[[9]][1,50])
tab_predictions[11,1]<-mean(mass_left_simulated_local_listoflists[[10]][1,50])

tab_predictions[1,2]<-mean(mass_left_simulated[,100])
tab_predictions[2,2]<-mean(mass_left_simulated_local_listoflists[[1]][,75])
tab_predictions[3,2]<-mean(mass_left_simulated_local_listoflists[[2]][,75])
tab_predictions[4,2]<-mean(mass_left_simulated_local_listoflists[[3]][,75])
tab_predictions[5,2]<-mean(mass_left_simulated_local_listoflists[[4]][,75])
tab_predictions[6,2]<-mean(mass_left_simulated_local_listoflists[[5]][,75])
tab_predictions[7,2]<-mean(mass_left_simulated_local_listoflists[[6]][,75])
tab_predictions[8,2]<-mean(mass_left_simulated_local_listoflists[[7]][,75])
tab_predictions[9,2]<-mean(mass_left_simulated_local_listoflists[[8]][,75])
tab_predictions[10,2]<-mean(mass_left_simulated_local_listoflists[[9]][1,75])
tab_predictions[11,2]<-mean(mass_left_simulated_local_listoflists[[10]][1,75])

tab_predictions[1,3]<-mean(mass_left_simulated[,120])
tab_predictions[2,3]<-mean(mass_left_simulated_local_listoflists[[1]][,120])
tab_predictions[3,3]<-mean(mass_left_simulated_local_listoflists[[2]][,120])
tab_predictions[4,3]<-mean(mass_left_simulated_local_listoflists[[3]][,120])
tab_predictions[5,3]<-mean(mass_left_simulated_local_listoflists[[4]][,120])
tab_predictions[6,3]<-mean(mass_left_simulated_local_listoflists[[5]][,120])
tab_predictions[7,3]<-mean(mass_left_simulated_local_listoflists[[6]][,120])
tab_predictions[8,3]<-mean(mass_left_simulated_local_listoflists[[7]][,120])
tab_predictions[9,3]<-mean(mass_left_simulated_local_listoflists[[8]][,120])
tab_predictions[10,3]<-mean(mass_left_simulated_local_listoflists[[9]][1,120])
tab_predictions[11,3]<-mean(mass_left_simulated_local_listoflists[[10]][1,120])


tab_predictions_diff<-mat.or.vec(10,3)
tab_predictions_diff[,1]<-tab_predictions[2:11,1]-tab_predictions[1,1]
tab_predictions_diff[,2]<-tab_predictions[2:11,2]-tab_predictions[1,2]
tab_predictions_diff[,3]<-tab_predictions[2:11,3]-tab_predictions[1,3]




mass_left_simulated_ave<-colMeans(mass_left_simulated)
year_list<-c(50,85,120)

tab_predictions_diff_MCMC<-mat.or.vec(10,3)
for(i in 1:3){
tab_predictions_diff_MCMC[1,i]<-mean((mass_left_simulated_local_listoflists[[1]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
tab_predictions_diff_MCMC[2,i]<-mean((mass_left_simulated_local_listoflists[[2]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
tab_predictions_diff_MCMC[3,i]<-mean((mass_left_simulated_local_listoflists[[3]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
tab_predictions_diff_MCMC[4,i]<-mean((mass_left_simulated_local_listoflists[[4]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
tab_predictions_diff_MCMC[5,i]<-mean((mass_left_simulated_local_listoflists[[5]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
tab_predictions_diff_MCMC[6,i]<-mean((mass_left_simulated_local_listoflists[[6]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
tab_predictions_diff_MCMC[7,i]<-mean((mass_left_simulated_local_listoflists[[7]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
tab_predictions_diff_MCMC[8,i]<-mean((mass_left_simulated_local_listoflists[[8]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
tab_predictions_diff_MCMC[9,i]<-mean((mass_left_simulated_local_listoflists[[9]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
tab_predictions_diff_MCMC[10,i]<-mean((mass_left_simulated_local_listoflists[[10]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
}

tab_predictions_diff_MCMC_min<-mat.or.vec(10,3)
for(i in 1:3){
  tab_predictions_diff_MCMC_min[1,i]<-min((mass_left_simulated_local_listoflists[[1]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_min[2,i]<-min((mass_left_simulated_local_listoflists[[2]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_min[3,i]<-min((mass_left_simulated_local_listoflists[[3]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_min[4,i]<-min((mass_left_simulated_local_listoflists[[4]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_min[5,i]<-min((mass_left_simulated_local_listoflists[[5]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_min[6,i]<-min((mass_left_simulated_local_listoflists[[6]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_min[7,i]<-min((mass_left_simulated_local_listoflists[[7]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_min[8,i]<-min((mass_left_simulated_local_listoflists[[8]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_min[9,i]<-min((mass_left_simulated_local_listoflists[[9]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_min[10,i]<-min((mass_left_simulated_local_listoflists[[10]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
}

tab_predictions_diff_MCMC_max<-mat.or.vec(10,3)
for(i in 1:3){
  tab_predictions_diff_MCMC_max[1,i]<-max((mass_left_simulated_local_listoflists[[1]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_max[2,i]<-max((mass_left_simulated_local_listoflists[[2]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_max[3,i]<-max((mass_left_simulated_local_listoflists[[3]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_max[4,i]<-max((mass_left_simulated_local_listoflists[[4]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_max[5,i]<-max((mass_left_simulated_local_listoflists[[5]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_max[6,i]<-max((mass_left_simulated_local_listoflists[[6]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_max[7,i]<-max((mass_left_simulated_local_listoflists[[7]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_max[8,i]<-max((mass_left_simulated_local_listoflists[[8]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_max[9,i]<-max((mass_left_simulated_local_listoflists[[9]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_max[10,i]<-max((mass_left_simulated_local_listoflists[[10]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
}

tab_predictions_diff_MCMC_sd<-mat.or.vec(10,3)
for(i in 1:3){
  tab_predictions_diff_MCMC_sd[1,i]<-sd((mass_left_simulated_local_listoflists[[1]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_sd[2,i]<-sd((mass_left_simulated_local_listoflists[[2]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_sd[3,i]<-sd((mass_left_simulated_local_listoflists[[3]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_sd[4,i]<-sd((mass_left_simulated_local_listoflists[[4]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_sd[5,i]<-sd((mass_left_simulated_local_listoflists[[5]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_sd[6,i]<-sd((mass_left_simulated_local_listoflists[[6]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_sd[7,i]<-sd((mass_left_simulated_local_listoflists[[7]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_sd[8,i]<-sd((mass_left_simulated_local_listoflists[[8]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_sd[9,i]<-sd((mass_left_simulated_local_listoflists[[9]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
  tab_predictions_diff_MCMC_sd[10,i]<-sd((mass_left_simulated_local_listoflists[[10]][,year_list[i]])-mass_left_simulated_ave[year_list[i]])
}




#png(filename = "Predictions.png", height =2600, width=2500, res=300)
png(filename = "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Predictions.png", height =2000, width=2500, res=300)
par(mar=c(5,15,2,2))
barplot<-barplot(t(tab_predictions_diff_MCMC)*100,beside=TRUE, names.arg=unique.decay.tree.class.text, las=2, horiz=TRUE,
        col=c("darkviolet", "firebrick", "chartreuse"), angle=c(40,60, 120), density=35, xlab="Difference in predictions, % C remaining", xlim=c(-40,40))
legend("topright", c("50 years", "85 years", "120 years"),fill=c("darkviolet", "firebrick", "chartreuse"), angle=c(40,60,120), density=35, bty="n")
box()
abline(v=0, lwd=1.5)
arrows(t(tab_predictions_diff_MCMC-tab_predictions_diff_MCMC_sd)*100, barplot, t(tab_predictions_diff_MCMC+tab_predictions_diff_MCMC_sd)*100, barplot, angle=90, code=3, length=0.04, col="gray40")
dev.off()


#png(filename = "Predictions_presentation.png", height =1600, width=2800, res=300)
png(filename = "./../../../../Dropbox/Apps/Overleaf/Modeling local persistence of coarse dead wood residuals in managed forests/Predictions_presentation.png", height =1600, width=2800, res=300)
par(mar=c(5,15,2,2))
barplot<-barplot(t(tab_predictions_diff_MCMC)*100,beside=TRUE, names.arg=c(as.character(RMSE_global[,1])), las=2, horiz=TRUE,
                 col=c("darkviolet", "firebrick", "chartreuse"), angle=c(40,60, 120), density=35, xlab="Difference in predictions, % C remaining", xlim=c(-40,40))
legend("topright", c("50 years", "85 years", "120 years"),fill=c("darkviolet", "firebrick", "chartreuse"), angle=c(40,60,120), density=35, bty="n")
box()
abline(v=0, lwd=1.5)
arrows(t(tab_predictions_diff_MCMC-tab_predictions_diff_MCMC_sd)*100, barplot, t(tab_predictions_diff_MCMC+tab_predictions_diff_MCMC_sd)*100, barplot, angle=90, code=3, length=0.04, col="gray40")
dev.off()




mean(abs(tab_predictions[1,1]-tab_predictions[,1]))
mean(abs(tab_predictions[2,1]-tab_predictions[,2]))
mean(abs(tab_predictions[3,1]-tab_predictions[,3]))



save.image("~/Documents/Q/Energimyndhigeten/Model_calibration/workspace.RData")
