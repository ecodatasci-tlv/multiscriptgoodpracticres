# ex1_m2.r
# 
#   Copyright 2014 Matthew V Talluto, Isabelle Boulangeat, Dominique Gravel
# 
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 3 of the License, or (at
#   your option) any later version.
#   
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   General Public License for more details.
# 
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
# 
# 
# run the integrated model


## set the RNG seed for reproducibility
## disable this for a true random run
# library(sp)
# library(raster)
# library(ggplot2)
# library(caTools)
# library(BIOMOD)
# library(ecospat)
# library(dismo)
# library(spaa)
# library(SDMTools)
# source("niche.overlap.functions.r")
# source("ex1_globals.r")
coef_data <- read.csv(file = paste(save_path,"physio_data.csv", sep=""))
# use model2 to compute priors for the metamodel
coef_data$tau <- 1/(coef_data$sd^2)

# factor by which to reduce the precision of the prior; is allows the user to decrease the confidence
# in the prior estimates
precisionReduction <- 10

## compute the integrated meta-model
#halophila_coordinates<-read.csv("coordinates.csv") #newdata is combined with the ones from eilat workshop
dups <- duplicated(halophila_coordinates[,]) #finding duplicates
halophila_coordinates <- halophila_coordinates[!dups,] #removing duplicates

#just for mediterranean:
halophila_med <- subset(halophila_coordinates, halophila_coordinates$region=="Mediterranean")#using only MED data
#turning the coordinates into binary form- 1 for presence
sp.occr_loc_med=data.frame(lon=halophila_med$lon,lat=halophila_med$lat)
sp_name_med=data.frame(Occurrence=rep(1,nrow(sp.occr_loc_med)))
species_t_med=SpatialPointsDataFrame(coords=sp.occr_loc_med,data=sp_name_med)

preds_r=stackOpen("env_data.stk")#opening recent envioronmental conditions
preds_d=stackOpen("env_depth_res.stk")
newext <- c(-10, 39, 30, 48)#coordinates to crop by
preds_r_med=crop(preds_r, newext)#conditions by croped region
preds_d_med=crop(preds_d, newext)
depth_meter = reclassify(preds_d_med,d_mask_mtx)
preds_d_m=raster::mask(preds_r_med,depth_meter)
preds_r_med=preds_d_m
preds_r_med <- aggregate(preds_r_med, fact=2, fun=mean)

preds_rm <- preds_r_med
rm_co=preds_rm
newext <- c(35, 37, 30, 33)
rm_co_df=crop(rm_co, newext)
data_rm<-as(rm_co_df, 'SpatialPolygons')
data_rest <- preds_r_med
rest_med_p <- raster::mask(data_rest,data_rm,inverse=TRUE)
preds_r_med<-rest_med_p
plot(preds_r_med)
#extract the coordinates from the binary form and removing all NA's
temp<-extract(preds_r_med,species_t_med)
temp3 <- as.data.frame(species_t_med)
temp3 <- data.frame(lon=temp3$lon, lat=temp3$lat, temp)
temp3<-na.omit(temp3)
temp4 <- data.frame(temp3$lon, temp3$lat)
sp_name_med=data.frame(Occurrence=rep(1,nrow(temp4)))
species_t_med=SpatialPointsDataFrame(coords=temp4,data=sp_name_med)

env_data<-as.data.frame(rasterToPoints(preds_r_med))
#pseudo_points <- data.frame(x=env_data[,1], y=env_data[,2])

presence_data <- data.frame(x=temp3$lon, y=temp3$lat)
#pseudo_points <- as.matrix(rbind(pseudo_points,presence_data))
#pa <- c(rep(0, times=nrow(env_data)),rep(1,times=nrow(presence_data)))

#pseudo <- pseudo.abs(pseudo_points,pa,strategy = "circles", distance = 2, nb.points=nrow(presence_data)*5, add.pres = FALSE)
#sp.absent <- as.data.frame(pseudo_points[pseudo,])

coef_b <- matrix(nrow = 0, ncol=5)
for (i in 1:100){
  
  bg.profile <- OCSVMprofiling(xy = presence_data, varstack = preds_r_med)
  bagr.radius <- backgroundRadius(xy = presence_data, background = bg.profile$absence,
                                  start = 0.249, by = 0.249, unit = "decimal degrees") 
  pseudo_TS <- pseudoAbsences(xy = presence_data,
                              background = bagr.radius,
                              exclusion.buffer = 0.249, prevalence = -1.5) # (-1.5) means 5 times of presence
  
  pseudo_all <- pseudo_TS$species1$PA01$km2070 #this is the maximum distance between presence and absence
  sp.absent <-subset(pseudo_all[,1:2], pseudo_all[,3] == 0)
  
  pred_pres <- extract(preds_r_med,species_t_med)
  pred_abs <- extract(preds_r_med, sp.absent)
  pred_pres <-cbind(pred_pres, presence=rep(1,nrow(pred_pres)))
  pred_abs <-cbind(pred_abs, presence=rep(0,nrow(pred_abs)))
  pred <- as.data.frame(rbind(pred_pres,pred_abs))

# run a simple bayesian logistic regression on the data
source("combined_model_selection.R")
# create the model
model1 <- jags.model(jags_model_name,inits = jags.inits ,data = m1Data, n.chains=settings$chains, n.adapt=settings$tuning )
update(model1, settings$burnin)

m1Results <- coda.samples(model1, m1Pars, settings$samples*settings$thin, thin=settings$thin)


#analyzing the data:
coef_b <- rbind(coef_b, as.matrix(m1Results))
}
#creating histograms for the coefficients:

for (i in 1:ncol(coef_b)){
  col_name <- colnames(coef_b)[i]
  png(filename = paste(save_path,col_name,"_hist_combined.png",sep = ""))
  hist(coef_b[,i], xlab = col_name, col = "blue", breaks = 50, main = paste("Histogram for coefficient",col_name))
  dev.off()
}


#creating a probability graph:

minsst_rng <- c(0:20)
maxsst_rng <- c(18:35)

data_sstmin <- data.frame(minsst_rng,rep(median(temp3$Temperature_Max)))
data_sstmax <- data.frame(rep(median(temp3$Temperature_Min)),maxsst_rng)
prob_sstmin_com <- na.omit(predictions(data_sstmin,coef_b))
prob_sstmax_com <- na.omit(predictions(data_sstmax,coef_b))

data_temp_rng <- expand.grid(minsst_rng,maxsst_rng)
prob_temp_rng <- na.omit(predictions(data_temp_rng,coef_b))
#for two predictors together:
prob_mean <- colMeans(prob_temp_rng)
SD_prob <- sapply(prob_temp_rng,function(x)sd(x))
m_SD <- cbind(data_temp_rng,prob_mean,SD_prob)
m_SD <- as.data.frame(m_SD)

prob_plot <- ggplot(m_SD, aes(x=m_SD$Var1, y=m_SD$Var2, z = m_SD$prob_mean)) +
  geom_raster(aes(fill = m_SD$prob_mean), interpolate = TRUE)+
  scale_fill_gradient(low = "#91C8F8", high = "#0E5797", name="Mean")+
  labs(x="Min SST", y="Max SST")
prob_plot
ggsave(filename = paste(save_path,"prob_mean_combined.png",sep = ""))

prob_plot <- ggplot(m_SD, aes(x=m_SD$Var1, y=m_SD$Var2, z = m_SD$SD_prob)) +
  geom_raster(aes(fill = m_SD$SD_prob), interpolate = TRUE)+
  scale_fill_gradient(low = "#91C8F8", high = "#0E5797",name="SD")+
  labs(x="Min SST", y="Max SST")
prob_plot
ggsave(filename = paste(save_path,"prob_sd_combined.png",sep = ""))


#for minsst:
prob_sstmin_mean_com <- colMeans(prob_sstmin_com)
SD_prob_sstmin_com <- sapply(prob_sstmin_com,function(x)sd(x))
percentiles <- sapply(prob_sstmin_com,function(x)quantile(x, probs = seq(0,1,0.25)))
percentiles <- percentiles[-c(1,3,5),]
percentiles<- t(percentiles)
m_SD_min_com <- cbind(minsst_rng,prob_sstmin_mean_com,SD_prob_sstmin_com,percentiles)
m_SD_min_com <- as.data.frame(m_SD_min_com)
names(m_SD_min_com)[4] <- paste("twentyfive_per")
names(m_SD_min_com)[5] <- paste("seventyfive_per")


ploting(m_SD_min_com,"The probability to find the species in a certain minimum sea surface temperature",
        expression(paste("Temperature (", degree~c,")")))
ggsave(filename = paste(save_path,"minsst_combined.png",sep = ""))

minsst_plot <- ggplot() +
  # blue plot for sst_min from combined model
  geom_point(data=m_SD_min_com, aes(x=m_SD_min_com[,1], y=m_SD_min_com[,2])) + 
  geom_line(data=m_SD_min_com, aes(x=m_SD_min_com[,1], y=m_SD_min_com[,2]), 
              colour="darkblue", size=1) +
  geom_ribbon(aes(x=m_SD_min_com[,1], ymin = m_SD_min_com[,4], ymax = m_SD_min_com[,5], fill="blue"),alpha = .25) +
  # red plot for sst_min from sdm regular
  geom_point(data=m_SD_min, aes(x=m_SD_min[,1], y=m_SD_min[,2])) + 
  geom_line(data=m_SD_min, aes(x=m_SD_min[,1], y=m_SD_min[,2]),
              colour="red", size=1)+
  geom_ribbon(aes(x=m_SD_min[,1],ymin = m_SD_min[,4], ymax = m_SD_min[,5], fill="red"),alpha = .25) +
  scale_fill_identity(name = 'Model', guide = 'legend',labels = c('Combined SDM','Regular SDM')) +
  labs(x=expression(paste("Temperature (", degree~c,")")), y="Probability")+
  ggtitle("Comparison between regular SDM and combined SDM")
minsst_plot
ggsave(filename = paste(save_path,"minsst_plot_combined.png",sep = ""))

#for maxsst:
prob_sstmax_mean_com <- colMeans(prob_sstmax_com)
SD_prob_sstmax_com <- sapply(prob_sstmax_com,function(x)sd(x))
percentiles <- sapply(prob_sstmax_com,function(x)quantile(x, probs = seq(0,1,0.25)))
percentiles <- percentiles[-c(1,3,5),]
percentiles<- t(percentiles)
m_SD_max_com <- cbind(maxsst_rng,prob_sstmax_mean_com,SD_prob_sstmax_com,percentiles)
m_SD_max_com <- as.data.frame(m_SD_max_com)
names(m_SD_max_com)[4] <- paste("twentyfive_per")
names(m_SD_max_com)[5] <- paste("seventyfive_per")


ploting(m_SD_max_com,"The probability to find the species in a certain maximum sea surface temperature",
        expression(paste("Temperature (", degree~c,")")))
ggsave(filename = paste(save_path,"maxsst_combined.png",sep = ""))

maxsst_plot <- ggplot() +
  # blue plot for sst_max from combined model
  geom_point(data=m_SD_max_com, aes(x=m_SD_max_com[,1], y=m_SD_max_com[,2])) + 
  geom_line(data=m_SD_max_com, aes(x=m_SD_max_com[,1], y=m_SD_max_com[,2]), 
              colour="blue", size=1) +
  geom_ribbon(aes(x=m_SD_max_com[,1],ymin = m_SD_max_com[,4], ymax = m_SD_max_com[,5], fill="blue"),alpha = .25) +
  # red plot for sst_max from sdm regular
  geom_point(data=m_SD_max, aes(x=m_SD_max[,1], y=m_SD_max[,2])) + 
  geom_line(data=m_SD_max, aes(x=m_SD_max[,1], y=m_SD_max[,2]), 
              colour="red", size=1)+
  geom_ribbon(aes(x=m_SD_max[,1],ymin = m_SD_max[,4], ymax = m_SD_max[,5], fill="red"),alpha = .25) +
  scale_fill_identity(name = 'Model', guide = 'legend',labels = c('Combined SDM','Regular SDM')) +
  labs(x=expression(paste("Temperature (", degree~c,")")), y="Probability")+
  ggtitle("Comparison between regular SDM and combined SDM")
maxsst_plot
ggsave(filename = paste(save_path,"maxsst_plot_combined.png",sep = ""))

#making probability and uncertainty maps:
predictors_coo_com <- rasterToPoints(preds_r_med)
predictors_coo_com <- as.data.frame(na.omit(predictors_coo_com))
names(predictors_coo_com)[1] <- paste("lon")
names(predictors_coo_com)[2] <- paste("lat")

predictors_data_com <- predictors_coo_com[,3:4]

prob_sp_com <- na.omit(predictions(predictors_data_com,coef_b))

prob_sp_mean_com <- colMeans(prob_sp_com)
max_prob_sp_mean_com <- max(prob_sp_mean_com)
data_m_p_c <- rasterFromXYZ(data.frame(x=predictors_coo_com$lon,y=predictors_coo_com$lat, z=prob_sp_mean_com))
png(filename = paste(save_path,"mean_map_present_combined.png",sep = ""))
plot(data_m_p_c, zlim=c(0,max_prob_sp_mean_com))
points(temp3)
dev.off()

SD_prob_sp_com <- apply(prob_sp_com, 2, sd)
data_sd_p_c <- rasterFromXYZ(data.frame(x=predictors_coo_com$lon,y=predictors_coo_com$lat, z=SD_prob_sp_com))
png(filename = paste(save_path,"SD_map_present_combined.png",sep = ""))
plot(data_sd_p_c, zlim=c(0,1))
dev.off()

diff_p_cr <- (prob_sp_mean_com-prob_sp_mean)
diff_p_cr_mean_map <- rasterFromXYZ(data.frame(x=predictors_coo_com$lon,y=predictors_coo_com$lat, z=diff_p_cr))
png(filename = paste(save_path,"mean_difference_map_present.png",sep = ""))
plot(diff_p_cr_mean_map)
dev.off()

####future prediction maps:
preds_fu=stackOpen("min_max_data_future.stk")#opening recent envioronmental conditions
preds_d=stackOpen("env_depth_res.stk")
new_preds_fu_med=resample(preds_fu, depth_meter, method="bilinear")
preds_fu = new_preds_fu_med
newext <- c(-10, 39, 30, 48)#coordinates to crop by
preds_fu_med=crop(preds_fu, newext)#conditions by croped region
preds_d_med=crop(preds_d, newext)
depth_meter = reclassify(preds_d_med,d_mask_mtx)
preds_d_m=raster::mask(new_preds_fu_med,depth_meter)
preds_fu_med=preds_d_m
preds_fu_med <- aggregate(preds_fu_med, fact=2, fun=mean)
plot(preds_fu_med)

predictors_coo_fu_com <- rasterToPoints(preds_fu_med)
predictors_coo_fu_com <- as.data.frame(na.omit(predictors_coo_fu_com))
names(predictors_coo_fu_com)[1] <- paste("lon")
names(predictors_coo_fu_com)[2] <- paste("lat")

predictors_data_fu_com <- predictors_coo_fu_com[,3:4]

prob_sp_fu_com <- na.omit(predictions(predictors_data_fu_com,coef_b))

prob_sp_mean_fu_com <- colMeans(prob_sp_fu_com)
max_prob_sp_mean_fu_com <- max(prob_sp_mean_fu_com)
data_m_f_c <- rasterFromXYZ(data.frame(x=predictors_coo_fu_com$lon,y=predictors_coo_fu_com$lat, z=prob_sp_mean_fu_com))
png(filename = paste(save_path,"mean_map_future_combined.png",sep = ""))
plot(data_m_f_c, zlim=c(0,max_prob_sp_mean_fu_com))
dev.off()

SD_prob_sp_fu_com <- apply(prob_sp_fu_com, 2, sd)
data_sd_f_c <- rasterFromXYZ(data.frame(x=predictors_coo_fu_com$lon,y=predictors_coo_fu_com$lat, z=SD_prob_sp_fu_com))
png(filename = paste(save_path,"SD_map_future_combined.png",sep = ""))
plot(data_sd_f_c, zlim=c(0,1))
dev.off()

diff_f_cr <- (prob_sp_mean_fu_com-prob_sp_mean_fu)
diff_f_cr_mean_map <- rasterFromXYZ(data.frame(x=predictors_coo_fu_com$lon,y=predictors_coo_fu_com$lat, z=diff_f_cr))
png(filename = paste(save_path,"mean_difference_map_future.png",sep = ""))
plot(diff_f_cr_mean_map)
dev.off()

## cross-validation of the model:
cross_val_c <- pred
AUC_mine_c <- rep(NA,100)
for (i in 1:100){
  train_rows <- sample(nrow(cross_val_c))
  train_m <- pred[train_rows[1:(nrow(cross_val_c)*0.7)],]
  test_m <- pred[train_rows[(nrow(cross_val_c)*0.7):nrow(cross_val_c)],]
  
    
    CData <- m1Data
    CData$sstmin<-train_m$Temperature_Min
    CData$sstmax<-train_m$Temperature_Max
    CData$presence <- train_m$presence
    
  
  # create the model
  model2 <- jags.model(jags_model_name, inits = jags.inits,data = CData, n.chains=settings$chains, n.adapt=settings$tuning )
  update(model2, settings$burnin)
  
  m2Results <- coda.samples(model2, m1Pars, settings$cv_samples*settings$thin, thin=settings$thin)
  
  coef_CV <- as.matrix(m2Results)
  CV_pred <- na.omit(predictions(test_m,coef_CV))
  CV_pred_mean <- colMeans(CV_pred)
  test_m$prediction <- CV_pred_mean
  
  test_m_pres <- subset(test_m, test_m$presence ==1)
  test_m_abs <- subset(test_m, test_m$presence ==0)
  th <- seq(1,0,by=-0.02)
  true_positive <- rep(NA,length(th))
  false_positive <- rep(NA,length(th))
  for (j in 1:length(th)){
    true_positive[j] <-sum(test_m_pres$prediction>th[j])/nrow(test_m_pres)
    false_positive[j] <-sum(test_m_abs$prediction>th[j])/nrow(test_m_abs)
  }
  
  validation <- rbind(true_positive,false_positive)
  validation <- as.data.frame(t(validation))
  cross_validation <- ggplot(validation, aes(x=validation[,2],y=validation[,1]))+
    xlim(0,1)+ylim(0,1)+
    geom_line(size=1, color = "blue")
  cross_validation <- cross_validation + labs(x="false_positive", y="true_positive")
  cross_validation
  ggsave(filename = paste(save_path,"ROC_combined.png",sep = ""))
  
  AUC_mine_c[i] <- trapz(validation$false_positive, validation$true_positive)
}
AUC_mean_c = mean(na.omit(AUC_mine_c))

#Cross-validation using BOYCE index:
temp3$presence <- rep(1, nrow(temp3))
absence_data <- cbind(sp.absent,pred_abs)
names(absence_data) <-c("lon","lat","Temperature_Min","Temperature_Max","presence")
cross_val_c_B<- rbind(temp3,absence_data)
B_Index_mine_c <- rep(NA,100)

for (i in 1:100){
  train_rows <- sample(nrow(cross_val_c_B))
  train_m <- cross_val_c_B[train_rows[1:(nrow(cross_val_c_B)*0.7)],]
  test_m <- cross_val_c_B[train_rows[(nrow(cross_val_c_B)*0.7):nrow(cross_val_c_B)],]
  
   
  CData <- m1Data
  CData$sstmin<-train_m$Temperature_Min
  CData$sstmax<-train_m$Temperature_Max
  CData$presence <- train_m$presence
  
  # create the model
  model2 <- jags.model(jags_model_name, inits = jags.inits,data = CData, n.chains=settings$chains, n.adapt=settings$tuning )
  update(model2, settings$burnin)
  
  m2Results <- coda.samples(model2, m1Pars, settings$cv_samples*settings$thin, thin=settings$thin)
  
  coef_CV <- as.matrix(m2Results)
  test_m = subset(test_m, test_m$presence==1)
  test_m_matrix <- as.matrix(test_m)
  # suitability_data <- raster(subset(test_m_matrix,select=c(lon,lat,prediction)))
  presence_data <- subset(test_m,select=c(lon,lat))
  presence_data_train <- subset(train_m,select=c(lon,lat))
  
  predictors_coo_BO <- rasterToPoints(preds_r_med)
  predictors_coo_BO <- as.data.frame(na.omit(predictors_coo_BO))
  names(predictors_coo_BO)[1] <- paste("lon")
  names(predictors_coo_BO)[2] <- paste("lat")
  
  predictors_data_BO <- predictors_coo_BO[,3:4]
  
  prob_sp_BO <- na.omit(predictions(predictors_data_BO,coef_CV))
  
  prob_sp_BO_mean <- colMeans(prob_sp_BO)
  data_m_p_BO <- rasterFromXYZ(data.frame(x=predictors_coo_BO$lon,y=predictors_coo_BO$lat, z=prob_sp_BO_mean))
  
  boyce_index <- ecospat.boyce(data_m_p_BO,presence_data, window.w = 0.2)
  B_Index_mine_c[i] <- boyce_index$Spearman.cor
  
}
BO_mean_c = mean(na.omit(B_Index_mine_c))

# cross-validation by kappa:
# kappa_val_c <- pred
# kappa_mine_c <- rep(NA,100)
# for (i in 1:100){
#   train_rows <- sample(nrow(kappa_val_c))
#   train_m <- pred[train_rows[1:(nrow(kappa_val_c)*0.7)],]
#   test_m <- pred[train_rows[(nrow(kappa_val_c)*0.7):nrow(kappa_val_c)],]
#   
#   
#   CData <- m1Data
#   CData$sstmin<-train_m$Temperature_Min
#   CData$sstmax<-train_m$Temperature_Max
#   CData$presence <- train_m$presence
#   
#   # create the model
#   model2 <- jags.model(jags_model_name, inits = jags.inits,data = CData, n.chains=settings$chains, n.adapt=settings$tuning )
#   update(model2, settings$burnin)
#   
#   m2Results <- coda.samples(model2, m1Pars, settings$cv_samples*settings$thin, thin=settings$thin)
#   
#   coef_CV_kap <- as.matrix(m2Results)
#   CV_pred_kap <- na.omit(predictions(test_m,coef_CV_kap))
#   CV_pred_mean_kap <- colMeans(CV_pred_kap)
#   test_m$prediction <- CV_pred_mean_kap
#   
#   mat <- confusion.matrix(test_m$presence,test_m$prediction, threshold = 0.2)
#   kappa_mine_c[i] <- Kappa(mat)
# }
# kap_mean_c = mean(na.omit(kappa_mine_c))

## cross-validation by region:
pres_points <- data.frame(temp3$lon, temp3$lat)
for (i in 1: nrow(pres_points)){
  if (pres_points$temp3.lon[i]>(-10) & pres_points$temp3.lon[i]<23 &
      pres_points$temp3.lat[i]>30 & pres_points$temp3.lat[i] < 37){
    pres_points$index[i] <- 1 
  } else if (pres_points$temp3.lon[i]>(-10) & pres_points$temp3.lon[i]<23 &
             pres_points$temp3.lat[i]>37 & pres_points$temp3.lat[i] < 48){
    pres_points$index[i] <- 2
  }else if (pres_points$temp3.lon[i]>(23) & pres_points$temp3.lon[i]<39 &
            pres_points$temp3.lat[i]>30 & pres_points$temp3.lat[i] < 37){
    pres_points$index[i] <- 3 # includes the levant area
  } else {pres_points$index[i] <- 4}
}

abs_points <- data.frame(sp.absent$x, sp.absent$y)# x=lon, y=lat
for (i in 1: nrow(abs_points)){
  if (abs_points$sp.absent.x[i]>(-10) & abs_points$sp.absent.x[i]<23 &
      abs_points$sp.absent.y[i]>30 & abs_points$sp.absent.y[i] < 37){
    abs_points$index[i] <- 1 
  } else if (abs_points$sp.absent.x[i]>(-10) & abs_points$sp.absent.x[i]<23 &
             abs_points$sp.absent.y[i]>37 & abs_points$sp.absent.y[i] < 48){
    abs_points$index[i] <- 2
  }else if (abs_points$sp.absent.x[i]>(23) & abs_points$sp.absent.x[i]<39 &
            abs_points$sp.absent.y[i]>30 & abs_points$sp.absent.y[i] < 37){
    abs_points$index[i] <- 3 # includes the levant area
  } else {abs_points$index[i] <- 4}
}
#test:
region_three_p <- subset(pres_points, pres_points$index==3)
region_three_p <- region_three_p[,1:2]
sp_name_med_cv=data.frame(Occurrence=rep(1,nrow(region_three_p)))
species_t_med_cv_p=SpatialPointsDataFrame(coords=region_three_p,data=sp_name_med_cv)

region_three_a <- subset(abs_points, abs_points$index==3)
region_three_a <- region_three_a[,1:2]
sp_name_med_cv=data.frame(Occurrence=rep(1,nrow(region_three_a)))
species_t_med_cv_a=SpatialPointsDataFrame(coords=region_three_a,data=sp_name_med_cv)

cv_pres <- extract(preds_r_med,species_t_med_cv_p)
cv_abs <- extract(preds_r_med, species_t_med_cv_a)
cv_pres <-cbind(cv_pres, presence=rep(1,nrow(cv_pres)))
cv_abs <-cbind(cv_abs, presence=rep(0,nrow(cv_abs)))
cv_test <- as.data.frame(rbind(cv_pres,cv_abs))

bo_test_p <- as.data.frame(cbind(cv_pres,region_three_p))
bo_test_a <- as.data.frame(cbind(cv_abs,region_three_a))
names(bo_test_a)<-c("Temperature_Min","Temperature_Max","presence","temp3.lon","temp3.lat")
bo_test <- rbind(bo_test_p,bo_test_a)

#train:
region_three_p <- subset(pres_points, pres_points$index!=3)
region_three_p <- region_three_p[,1:2]
sp_name_med_cv=data.frame(Occurrence=rep(1,nrow(region_three_p)))
species_t_med_cv_p=SpatialPointsDataFrame(coords=region_three_p,data=sp_name_med_cv)

region_three_a <- subset(abs_points, abs_points$index!=3)
region_three_a <- region_three_a[,1:2]
sp_name_med_cv=data.frame(Occurrence=rep(1,nrow(region_three_a)))
species_t_med_cv_a=SpatialPointsDataFrame(coords=region_three_a,data=sp_name_med_cv)

cv_pres <- extract(preds_r_med,species_t_med_cv_p)
cv_abs <- extract(preds_r_med, species_t_med_cv_a)
cv_pres <-cbind(cv_pres, presence=rep(1,nrow(cv_pres)))
cv_abs <-cbind(cv_abs, presence=rep(0,nrow(cv_abs)))
cv_train <- as.data.frame(rbind(cv_pres,cv_abs))

bo_train_p <- as.data.frame(cbind(cv_pres,region_three_p))
bo_train_a <- as.data.frame(cbind(cv_abs,region_three_a))
names(bo_train_a)<-c("Temperature_Min","Temperature_Max","presence","temp3.lon","temp3.lat")
bo_train <- rbind(bo_train_p,bo_train_a)

train <- cv_train
test <- cv_test
CData <- m1Data
CData$sstmin<-train$Temperature_Min
CData$sstmax<-train$Temperature_Max
CData$presence <- train$presence


# create the model
model2 <- jags.model(jags_model_name,inits = jags.inits, data = CData, n.chains=settings$chains, n.adapt=settings$tuning )
update(model2, settings$burnin)

m2Results <- coda.samples(model2, m1Pars, settings$cv_samples*settings$thin, thin=settings$thin)

coef_CV <- as.matrix(m2Results)
CV_pred <- na.omit(predictions(test,coef_CV))
CV_pred_mean <- colMeans(CV_pred)
test$prediction <- CV_pred_mean

test_pres <- subset(test, test$presence ==1)
test_abs <- subset(test, test$presence ==0)
th <- seq(1,0,by=-0.02)
true_positive <- rep(NA,length(th))
false_positive <- rep(NA,length(th))
for (j in 1:length(th)){
  true_positive[j] <-sum(test_pres$prediction>th[j])/nrow(test_pres)
  false_positive[j] <-sum(test_abs$prediction>th[j])/nrow(test_abs)
}

validation <- rbind(true_positive,false_positive)
validation <- as.data.frame(t(validation))
cross_validation <- ggplot(validation, aes(x=validation[,2],y=validation[,1]))+
  xlim(0,1)+ylim(0,1)+
  geom_line(size=1, color = "blue")
cross_validation <- cross_validation + labs(x="false_positive", y="true_positive")
cross_validation
ggsave(filename = paste(save_path,"ROC_combined.png",sep = ""))

AUC_hot_region_c <- trapz(validation$false_positive, validation$true_positive)

#by region using BOYCE index:
cross_val_B<- bo_train

train_m <- cross_val_B
test_m <- bo_test

CData <- m1Data
CData$sstmin<-train_m$Temperature_Min
CData$sstmax<-train_m$Temperature_Max
CData$presence <- train_m$presence

# create the model
model2 <- jags.model(jags_model_name,inits = jags.inits, data = CData, n.chains=settings$chains, n.adapt=settings$tuning )
update(model2, settings$burnin)

m2Results <- coda.samples(model2, m1Pars, settings$cv_samples*settings$thin, thin=settings$thin)

coef_CV <- as.matrix(m2Results)
test_m = subset(test_m, test_m$presence==1)
test_m_matrix <- as.matrix(test_m)
# suitability_data <- raster(subset(test_m_matrix,select=c(lon,lat,prediction)))
presence_data <- subset(test_m,select=c(temp3.lon,temp3.lat))
presence_data_train <- subset(train_m,select=c(temp3.lon,temp3.lat))

predictors_coo_BO <- rasterToPoints(preds_r_med)
predictors_coo_BO <- as.data.frame(na.omit(predictors_coo_BO))
names(predictors_coo_BO)[1] <- paste("lon")
names(predictors_coo_BO)[2] <- paste("lat")

predictors_data_BO <- predictors_coo_BO[,3:4]

prob_sp_BO <- predictions(predictors_data_BO,coef_CV)

prob_sp_BO_mean <- colMeans(prob_sp_BO)
data_m_p_BO <- rasterFromXYZ(data.frame(x=predictors_coo_BO$lon,y=predictors_coo_BO$lat, z=prob_sp_BO_mean))

boyce_index <- ecospat.boyce(data_m_p_BO,presence_data, window.w = 0.2)
BO_mean_region_c <- boyce_index$Spearman.cor

# by region using kappa:

# train_m <- cv_train
# test_m <- cv_test
# 
# 
# CData <- m1Data
# CData$sstmin<-train_m$Temperature_Min
# CData$sstmax<-train_m$Temperature_Max
# CData$presence <- train_m$presence
# 
# 
# # create the model
# model2 <- jags.model(jags_model_name,inits = jags.inits, data = CData, n.chains=settings$chains, n.adapt=settings$tuning )
# update(model2, settings$burnin)
# 
# m2Results <- coda.samples(model2, m1Pars, settings$cv_samples*settings$thin, thin=settings$thin)
# 
# coef_CV_kap <- as.matrix(m2Results)
# CV_pred_kap <- predictions(test_m,coef_CV_kap)
# CV_pred_mean_kap <- colMeans(CV_pred_kap)
# test_m$prediction <- CV_pred_mean_kap
# 
# mat <- confusion.matrix(test_m$presence,test_m$prediction, threshold = 0.2)
# kap_mean_region_c <- Kappa(mat)

# train <- rbind(temp_hot,temp_medium)
#   CData <- m1Data
#   CData$sstmin<-train$Temperature_Min
#   CData$sstmax<-train$Temperature_Max
#   CData$presence <- train$presence
#   
#   
#   # create the model
#   model2 <- jags.model(jags_model_name, inits = jags.inits,data = CData, n.chains=settings$chains, n.adapt=settings$tuning )
#   update(model2, settings$burnin)
#   
#   m2Results <- coda.samples(model2, m1Pars, settings$cv_samples*settings$thin, thin=settings$thin)
#   
#   coef_CV <- as.matrix(m2Results)
#   CV_pred <- na.omit(predictions(temp_hot,coef_CV))
#   CV_pred_mean <- colMeans(CV_pred)
#   temp_hot$prediction <- CV_pred_mean
#   
#   temp_cold_pres <- subset(temp_cold, temp_cold$presence ==1)
#   temp_cold_abs <- subset(temp_cold, temp_cold$presence ==0)
#   th <- seq(1,0,by=-0.02)
#   true_positive <- rep(NA,length(th))
#   false_positive <- rep(NA,length(th))
#   for (j in 1:length(th)){
#     true_positive[j] <-sum(temp_cold_pres$prediction>th[j])/nrow(temp_cold_pres)
#     false_positive[j] <-sum(temp_cold_abs$prediction>th[j])/nrow(temp_cold_abs)
#   }
#   
#   validation <- rbind(true_positive,false_positive)
#   validation <- as.data.frame(t(validation))
#   cross_validation <- ggplot(validation, aes(x=validation[,2],y=validation[,1]))+
#     xlim(0,1)+ylim(0,1)+
#     geom_line(size=1, color = "blue")
#   cross_validation <- cross_validation + labs(x="false_positive", y="true_positive")
#   cross_validation
#   ggsave(filename = paste(save_path,"ROC_combined.png",sep = ""))
#   
#   AUC_to_cold_c <- trapz(validation$false_positive, validation$true_positive)
  
# overlap calculations:
p_vs_f_prob_data <- brick(data_m_p_c,data_m_f_c)
p_vs_f_prob_data <- na.omit(as.data.frame(rasterToPoints(p_vs_f_prob_data)))
names(p_vs_f_prob_data)=c("lon","lat","present","future")

levant_coo = subset(p_vs_f_prob_data, p_vs_f_prob_data$lon>28 & p_vs_f_prob_data$lon<37 &
                      p_vs_f_prob_data$lat>31 & p_vs_f_prob_data$lat < 36)

west_coo = subset(p_vs_f_prob_data, p_vs_f_prob_data$lon>(-5) & p_vs_f_prob_data$lon<17 &
                    p_vs_f_prob_data$lat>35 & p_vs_f_prob_data$lat < 43)
for (i in 1: nrow(p_vs_f_prob_data)){
  if (p_vs_f_prob_data$lon[i]>28 & p_vs_f_prob_data$lon[i]<37 &
      p_vs_f_prob_data$lat[i]>31 & p_vs_f_prob_data$lat[i] < 36){
    p_vs_f_prob_data$color[i] <- "red" 
  } else if (p_vs_f_prob_data$lon[i]>(-5) & p_vs_f_prob_data$lon[i]<17 &
             p_vs_f_prob_data$lat[i]>35 & p_vs_f_prob_data$lat[i] < 43){
    p_vs_f_prob_data$color[i] <- "light green"
  } else {p_vs_f_prob_data$color[i] <- "blue"}
}

cor.test(p_vs_f_prob_data$present,p_vs_f_prob_data$future)
reg1 <- lm(p_vs_f_prob_data$future~p_vs_f_prob_data$present)
p_f_overlap_com <- ggplot(p_vs_f_prob_data, aes(x=p_vs_f_prob_data[,3],y=p_vs_f_prob_data[,4], color = p_vs_f_prob_data$color))+
  xlim(0,1)+ylim(0,1)+
  geom_point(size=1.5)+
  scale_color_manual(values=unique(p_vs_f_prob_data$color),name = 'Overlap',labels = c('All of the Mediterranean','Western Mediterranean','Levant region'))+
  stat_smooth(method = "lm", col = "black")
p_f_overlap_com <- p_f_overlap_com + labs(x="Present mean probabilities", y="Future mean probabilities")
p_f_overlap_com
ggsave(filename = paste(save_path,"present_future_overlap_combined.png",sep = ""))
cor_coeff_f_vs_p_com<-coef(reg1)
cor_coeff_f_vs_p_com


###### niche similarity:
# for the all MED:
niche_overlap_I_stat_c <- nicheOverlap(data_m_p_c,data_m_f_c)
niche_overlap_D_stat_c <- nicheOverlap(data_m_p_c,data_m_f_c, stat = 'D')
# for the Levant:
present_com_data <- data_m_p_c
newext <- c(28,37, 31, 36)#coordinates to crop by
present_com_data=crop(present_com_data, newext)#conditions by croped region

future_com_data <- data_m_f_c
newext <- c(28,37, 31, 36)#coordinates to crop by
future_com_data=crop(future_com_data, newext)#conditions by croped region

niche_overlap_I_stat_c_L <- nicheOverlap(present_com_data,future_com_data)
niche_overlap_D_stat_c_L <- nicheOverlap(present_com_data,future_com_data, stat = 'D')

# for the rest of the MED:
present_data_poly_com<-as(present_com_data, 'SpatialPolygons')
present_data_rest_com <- data_m_p_c
rest_med_p_com <- raster::mask(present_data_rest_com,present_data_poly_com,inverse=TRUE)

future_data_poly_com<-as(future_com_data, 'SpatialPolygons')
future_data_rest_com <- data_m_f_c
rest_med_f_com <- raster::mask(future_data_rest_com,future_data_poly_com,inverse=TRUE)

niche_overlap_I_stat_c_R <- nicheOverlap(rest_med_p_com,rest_med_f_com)
niche_overlap_D_stat_c_R <- nicheOverlap(rest_med_p_com,rest_med_f_com, stat = 'D')

# for west med:
present_com_data_w <- data_m_p_c
newext <- c(-5,17, 35, 43)#coordinates to crop by
present_com_data_w=crop(present_com_data_w, newext)#conditions by croped region

future_com_data_w <- data_m_f_c
newext <- c(-5,17, 35, 43)#coordinates to crop by
future_com_data_w=crop(future_com_data_w, newext)#conditions by croped region

niche_overlap_I_stat_c_w <- nicheOverlap(present_com_data_w,future_com_data_w)
niche_overlap_D_stat_c_w <- nicheOverlap(present_com_data_w,future_com_data_w, stat = 'D')

# overlap between sdm regular and sdm combined:
#####present:
p_vs_p_prob_data <- brick(data_m_p_c,data_m_p)
p_vs_p_prob_data <- na.omit(as.data.frame(rasterToPoints(p_vs_p_prob_data)))
names(p_vs_p_prob_data)=c("lon","lat","present_combined","present_regular")

levant_coo_p = subset(p_vs_p_prob_data, p_vs_p_prob_data$lon>28 & p_vs_p_prob_data$lon<37 &
                      p_vs_p_prob_data$lat>31 & p_vs_p_prob_data$lat < 36)

west_coo_p = subset(p_vs_p_prob_data, p_vs_p_prob_data$lon>(-5) & p_vs_p_prob_data$lon<17 &
                      p_vs_p_prob_data$lat>35 & p_vs_p_prob_data$lat < 43)
for (i in 1: nrow(p_vs_p_prob_data)){
  if (p_vs_p_prob_data$lon[i]>28 & p_vs_p_prob_data$lon[i]<37 &
      p_vs_p_prob_data$lat[i]>31 & p_vs_p_prob_data$lat[i] < 36){
    p_vs_p_prob_data$color[i] <- "red" 
  } else if (p_vs_p_prob_data$lon[i]>(-5) & p_vs_p_prob_data$lon[i]<17 &
             p_vs_p_prob_data$lat[i]>35 & p_vs_p_prob_data$lat[i] < 43){
    p_vs_p_prob_data$color[i] <- "light green"
  } else {p_vs_p_prob_data$color[i] <- "blue"}
}

cor.test(p_vs_p_prob_data$present_combined,p_vs_p_prob_data$present_regular)
reg1 <- lm(p_vs_p_prob_data$present_regular~p_vs_p_prob_data$present_combined)
p_p_overlap_com <- ggplot(p_vs_p_prob_data, aes(x=p_vs_p_prob_data[,3],y=p_vs_p_prob_data[,4], color = p_vs_p_prob_data$color))+
  xlim(0,1)+ylim(0,1)+
  geom_point(size=1.5)+
  scale_color_manual(values=unique(p_vs_p_prob_data$color),name = 'Overlap',labels = c('All of the Mediterranean','Western Mediterranean','Levant region'))+
  stat_smooth(method = "lm", col = "black")
p_p_overlap_com <- p_p_overlap_com + labs(x="Combined mean probabilities", y="Regular mean probabilities")
p_p_overlap_com
ggsave(filename = paste(save_path,"present_combined_regular_overlap.png",sep = ""))
cor_coeff_p_vs_p_com<-coef(reg1)
cor_coeff_p_vs_p_com

###### niche similarity:
niche_overlap_I_stat_p_com_vs_reg <- nicheOverlap(data_m_p_c,data_m_p)
niche_overlap_D_stat_p_com_vs_reg <- nicheOverlap(data_m_p_c,data_m_p, stat = 'D')
# for the Levant:
p_com_data <- data_m_p_c
newext <- c(28,37, 31, 36)#coordinates to crop by
p_com_data=crop(p_com_data, newext)#conditions by croped region

present_data <- data_m_p
newext <- c(28,37, 31, 36)#coordinates to crop by
present_data=crop(present_data, newext)#conditions by croped region

niche_overlap_I_stat_p_L_com_vs_reg <- nicheOverlap(p_com_data,present_data)
niche_overlap_D_stat_p_L_com_vs_reg <- nicheOverlap(p_com_data,present_data, stat = 'D')

# for the rest of the MED:
present_data_poly_com<-as(p_com_data, 'SpatialPolygons')
present_data_rest_com <- data_m_p_c
rest_med_present_com <- raster::mask(present_data_rest_com,present_data_poly_com,inverse=TRUE)

present_data_poly<-as(present_data, 'SpatialPolygons')
present_data_rest <- data_m_p
rest_med_present <- raster::mask(present_data_rest,present_data_poly,inverse=TRUE)

niche_overlap_I_stat_R_com_vs_reg <- nicheOverlap(rest_med_present_com,rest_med_present)
niche_overlap_D_stat_R_com_vs_reg <- nicheOverlap(rest_med_present_com,rest_med_present, stat = 'D')

# for west med:
p_com_data_w <- data_m_p_c
newext <- c(-5,17, 35, 43)#coordinates to crop by
p_com_data_w=crop(p_com_data_w, newext)#conditions by croped region

present_data_w <- data_m_p
newext <- c(-5,17, 35, 43)#coordinates to crop by
present_data_w=crop(present_data_w, newext)#conditions by croped region

niche_overlap_I_stat_p_w_com_vs_reg <- nicheOverlap(p_com_data_w,present_data_w)
niche_overlap_D_stat_p_w_com_vs_reg <- nicheOverlap(p_com_data_w,present_data_w, stat = 'D')

####future:
f_vs_f_prob_data <- brick(data_m_f_c,data_m_f)
f_vs_f_prob_data <- na.omit(as.data.frame(rasterToPoints(f_vs_f_prob_data)))
names(f_vs_f_prob_data)=c("lon","lat","future_combined","future_regular")

levant_coo_f = subset(f_vs_f_prob_data, f_vs_f_prob_data$lon>28 & f_vs_f_prob_data$lon<37 &
                      f_vs_f_prob_data$lat>31 & f_vs_f_prob_data$lat < 36)

west_coo_f = subset(f_vs_f_prob_data, f_vs_f_prob_data$lon>(-5) & f_vs_f_prob_data$lon<17 &
                      f_vs_f_prob_data$lat>35 & f_vs_f_prob_data$lat < 43)
for (i in 1: nrow(f_vs_f_prob_data)){
  if (f_vs_f_prob_data$lon[i]>28 & f_vs_f_prob_data$lon[i]<37 &
      f_vs_f_prob_data$lat[i]>31 & f_vs_f_prob_data$lat[i] < 36){
    f_vs_f_prob_data$color[i] <- "red" 
  } else if (f_vs_f_prob_data$lon[i]>(-5) & f_vs_f_prob_data$lon[i]<17 &
             f_vs_f_prob_data$lat[i]>35 & f_vs_f_prob_data$lat[i] < 43){
    f_vs_f_prob_data$color[i] <- "light green"
  } else {f_vs_f_prob_data$color[i] <- "blue"}
}

cor.test(f_vs_f_prob_data$future_combined,f_vs_f_prob_data$future_regular)
reg1 <- lm(f_vs_f_prob_data$future_regular~f_vs_f_prob_data$future_combined)
f_f_overlap_com <- ggplot(f_vs_f_prob_data, aes(x=f_vs_f_prob_data[,3],y=f_vs_f_prob_data[,4], color = f_vs_f_prob_data$color))+
  xlim(0,1)+ylim(0,1)+
  geom_point(size=1.5)+
  scale_color_manual(values=unique(f_vs_f_prob_data$color),name = 'Overlap',labels = c('All of the Mediterranean','Western Mediterranean','Levant region'))+
  stat_smooth(method = "lm", col = "black")
f_f_overlap_com <- f_f_overlap_com + labs(x="Combined mean probabilities", y="Regular mean probabilities")
f_f_overlap_com
ggsave(filename = paste(save_path,"future_combined_regular_overlap.png",sep = ""))
cor_coeff_f_vs_f_com<-coef(reg1)
cor_coeff_f_vs_f_com

###### niche similarity:
niche_overlap_I_stat_f_com_vs_reg <- nicheOverlap(data_m_f_c,data_m_f)
niche_overlap_D_stat_f_com_vs_reg <- nicheOverlap(data_m_f_c,data_m_f, stat = 'D')
# for the Levant:
future_com_data <- data_m_f_c
newext <- c(28,37, 31, 36)#coordinates to crop by
future_com_data=crop(future_com_data, newext)#conditions by croped region

future_data <- data_m_f
newext <- c(28,37, 31, 36)#coordinates to crop by
future_data=crop(future_data, newext)#conditions by croped region

niche_overlap_I_stat_c_L_com_vs_reg <- nicheOverlap(future_data,future_com_data)
niche_overlap_D_stat_c_L_com_vs_reg <- nicheOverlap(future_data,future_com_data, stat = 'D')

# for the rest of the MED:
future_data_poly_com<-as(future_com_data, 'SpatialPolygons')
future_data_rest_com <- data_m_f_c
rest_med_future_com <- raster::mask(future_data_rest_com,future_data_poly_com,inverse=TRUE)

future_data_poly<-as(future_data, 'SpatialPolygons')
future_data_rest <- data_m_f
rest_med_future <- raster::mask(future_data_rest,future_data_poly,inverse=TRUE)

niche_overlap_I_stat_c_R_com_vs_reg <- nicheOverlap(rest_med_future_com,rest_med_future)
niche_overlap_D_stat_c_R_com_vs_reg <- nicheOverlap(rest_med_future_com,rest_med_future, stat = 'D')

# for west med:
future_com_data_w <- data_m_f_c
newext <- c(-5,17, 35, 43)#coordinates to crop by
future_com_data_w=crop(future_com_data_w, newext)#conditions by croped region

future_data_w <- data_m_f
newext <- c(-5,17, 35, 43)#coordinates to crop by
future_data_w=crop(future_data_w, newext)#conditions by croped region

niche_overlap_I_stat_c_w_com_vs_reg <- nicheOverlap(future_data_w,future_com_data_w)
niche_overlap_D_stat_c_w_com_vs_reg <- nicheOverlap(future_data_w,future_com_data_w, stat = 'D')

save.image(file=paste(save_path,"workspace_combined_bayes.RData", sep = ""))
