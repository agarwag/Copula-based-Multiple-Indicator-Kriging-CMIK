#############################################################
######################### LOOADING LIBRARIES ###############
##############################################################
library(dplyr)
###########################################################################################
######################### data set for one station just to see the variabels ###############
###########################################################################################
# STATIONS DATA
Stations<-read.csv("/home/agarwag/Desktop/KAUST/p4/Europe climate data/stations_rr.txt")
names(Stations)
dim(Stations)
str(Stations)

# SOURCE DATA 
Source <- read.csv("/home/agarwag/Desktop/KAUST/p4/Europe climate data/ECA_nonblend_info_rr.txt", comment.char="#")
names(Source)
dim(Source)
str(Source)
Source$SOUID = as.integer(as.character(Source$SOUID))

Stations <- within(Stations, {
  dms <- do.call(rbind, strsplit(as.character(LAT), ":"))
  dms_sign <- substr(as.character(LAT),1,1)
  dec.lat <-  as.numeric(paste0(dms_sign,(abs(as.numeric(dms[,1])) + (as.numeric(dms[,2]) + as.numeric(dms[,3])/60)/60)))
  rm(dms); rm(dms_sign)
})
#converting degrees to decimal longitude
Stations <- within(Stations, {
  dms <- do.call(rbind, strsplit(as.character(LON), ":"))
  dms_sign <- substr(as.character(LON),1,1)
  dec.lon <-  as.numeric(paste0(dms_sign,(abs(as.numeric(dms[,1])) + (as.numeric(dms[,2]) + as.numeric(dms[,3])/60)/60)))
  rm(dms); rm(dms_sign)
})


latitude.station<-Stations$dec.lat
longitude.station<-Stations$dec.lon
elevation.station<-Stations$HGHT

#plot(longitude.station, latitude.station)

######################################################################################################
######################### Converting latitude and longitude in decimal for Source data ###############
######################################################################################################

#Problem: Dimesions of the latitude, longitude, elevation does not match if we look for the unique values
Source =  Source[complete.cases(Source),] # remove NAs

#converting degrees to decimal Latitude
Source <- within(Source, {
  dms <- do.call(rbind, strsplit(as.character(LAT), ":"))
  dms_sign <- substr(as.character(LAT),1,1)
  dec.lat <-  as.numeric(paste0(dms_sign,(abs(as.numeric(dms[,1])) + (as.numeric(dms[,2]) + as.numeric(dms[,3])/60)/60)))
  rm(dms); rm(dms_sign)
})
#converting degrees to decimal longitude
Source <- within(Source, {
  dms <- do.call(rbind, strsplit(as.character(LON), ":"))
  dms_sign <- substr(as.character(LON),1,1)
  dec.lon <-  as.numeric(paste0(dms_sign,(abs(as.numeric(dms[,1])) + (as.numeric(dms[,2]) + as.numeric(dms[,3])/60)/60)))
  rm(dms); rm(dms_sign)
})



#Source<-Source[!duplicated(Source[,c("dec.lat", "dec.lon")]),]

latitude.source<-unique(Source$dec.lat)
length(latitude.source)
longitude.source<-unique(Source$dec.lon)
length(longitude.source)
elevation.source<-unique(Source$HGHT)
length(elevation.source)

#plot(Source$dec.lon, Source$dec.lat)
######################################################################################################
######################### Extracting the data for a country ###############
######################################################################################################
# France: FR, Switzerland: CH, ITALY: IT, GERMANY: DE, UNIKTED KINGDOM: GB, spain: ES
#NOTE: We need to look for the source data because for the precipitation data there is seperate file for each staions and their names are as the source id
Source_country <-Source %>% filter(Source$CN=="ES")

#plot(Source_country$LON,Source_country$LAT)
## 207 obs for Spain
#creating the index id to look for the particular country data 
id.source=unique(Source_country$SOUID)
name.source.data<-noquote(paste("RR_SOUID",id.source,".txt",sep = ""))

#Sample name RR_SOUID100003
######################################################################################################
############ Selecting the Stations data for germany (seperate file for each station) ##################
######################################################################################################
#NOTE:
#(1): 1 April to 31 Marcg 2019 for all stations
#(2): We Take stations where we have more than 50% data are available

data<-NULL
for (i in 1:length(id.source)) {
  file.name<-name.source.data[i]
  aux0=paste("/home/agarwag/Desktop/KAUST/p4/Europe climate data/ECA_nonblend_rr/",file.name,sep="")
  lines<-readLines(aux0)
  mytable.i <- read.table(text = lines[-c(1:18)], sep = ',', header = T)
  date <- gsub("(.*G.{3}).*","\\1",mytable.i$DATE)
  dat1 <- data.frame(Year=as.numeric(substr(date,1,4)), Month=as.numeric(substr(date,5,6)), Day=as.numeric(substr(date,7,8)),stringsAsFactors=FALSE); 
  mytable.i<-cbind(mytable.i,dat1)
  # extracting the data from 2018 to 2019 and only autumn
  year<-2019
  data.year<-NULL
  for (y in 1: length(year)) {
    table.i<-mytable.i %>% filter(mytable.i$Year==year[y])
    data.y<- table.i %>% filter(table.i$Month==9 |table.i$Month==10| table.i$Month==11|table.i$Month==12)
    data.year<-rbind(data.year,data.y)
  }
  mytable.i<- data.year
  data<-rbind(data,mytable.i)
}

####### replacing the missing  values, -9999 by NA
# 
data$RR<-replace(data$RR,which(data$RR==-9999),NA)
head(data)
tail(data)
dim(data)


######################################################################################################
############ Merge tha data with Station date to get the latitude, logitude and elevation ############
######################################################################################################

#Note: 
# (1): We need to merge with Stations data for the latitude and longitude but the infromation in data frame Stations are less 
# (2): We can not merge with Source dataframe with the source id as the in that case we have two diffrent sttaiuons have same longitude and latitude 

final.data<-merge(data, Source_country)
head(final.data)
names(final.data)


data_daily = subset(final.data,DATE=="20190901")
#hist(data_daily$RR)

data_monthly <- aggregate(RR ~ SOUID+Month+dec.lat+dec.lon, data=final.data, mean)


data_nov = subset(data_monthly, Month==11)
data_nov = data_nov[data_nov$dec.lat>35,]
# main = "Spain precipitation (November 2019)",
hist(data_nov$RR, xlab="Precipitation (mm)",main="", cex.main=1.5, cex.lab=1.4)

ggplot() + geom_histogram(aes(x=data_nov$RR), breaks=hist(data_nov$RR, plot=F)$breaks, col="black")+
  xlab("Precipitation (mm)")+ ylab("Frequency")+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.5)))
summary(data_nov$RR)
length(data_nov$RR)
#data_nov_plot = data_nov[data_nov$dec.lat>35,]

hist(log(data_nov$RR), xlab="Precipitation (mm)",main="", cex.main=1.5, cex.lab=1.4)

library(ggplot2) 
library(maps)
library(maptools)

spain = map_data("world", region = "spain")
#plot(spain)

ggplot(spain, aes(long, lat ))+geom_path(aes(group=group))+ coord_map()  + ggtitle(" Spain Precipitation (Nov 2019)")+
  geom_point(data = data_nov, aes(dec.lon, dec.lat, colour = RR)) +scale_color_viridis_c("Precipitation (mm)") 


load("/home/agarwag/Desktop/KAUST/p4/map_spain.RData")
library(ggmap)
library(viridis)
ggmap(map_spain)+
  #scale_x_continuous(limits = c(4,16))+
  #scale_y_continuous(limits = c(47,55))+
  geom_point(data = data_nov, aes(dec.lon, dec.lat, colour = RR), size=2) +
  scale_color_viridis("Precipitation (mm)", option = "D") +
  labs(x="Longitude (degrees)", y="Latitude (degrees)") +
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1.1)),
        legend.title = element_text(size=rel(1.3)))


library(rgdal)
####### Transforming Longitude Latitude to Universal Transverse Mercator (UTM) projections ####
LongLatToUTM<-function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## assigning a CRS (Coordinate Reference System)
  
  res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep=''))) ## Transforming the CRS
  return(data.frame(x= res$X, y= res$Y))
  #return(as.data.frame(res))
}

#spain zone 30
loc_m <- LongLatToUTM(data_nov$dec.lon,data_nov$dec.lat,zone = 30)

head(loc_m)

data_nov$x = loc_m[,1]
data_nov$y = loc_m[,2]
library(geoR);library(fields);library(maps)

## split the data into 80% training and 20% testing dataset
#set.seed(124)
set.seed(123)
n_iters = 500
train_data_list = list(); test_data_list=list()

for (index in 1:n_iters) {
  dt = sort(sample(nrow(data_nov), nrow(data_nov)*.8, replace = F))
  
  train_data_list[[index]] = data_nov[dt,]
  test_data_list[[index]] = data_nov[-dt,]
}

### OR

## split the data intok groups for k-fold crossvalidation
#set.seed(123)
#n_iters = k_cv = 10
#ss = sample(1:k_cv, nrow(data_nov), replace=T)
#for (index in 1:k_cv) {
#  train_data_list[[index]] = data_nov[ss!=index,]
#  test_data_list[[index]] = data_nov[ss==index,]
#}


parallel_gp = list()
for (index in 1:n_iters) {
  
  train_data = train_data_list[[index]]
  test_data = test_data_list[[index]]
  
  train_data_nov = train_data[,"RR"]; train_loc = cbind(x = train_data$x, y= train_data$y)
  
  #### predict cdfs on the test data
  test_data_nov =  test_data[,"RR"]; test_loc = cbind(x= test_data$x, y= test_data$y)
  
  data_train_gp = data.frame(train_loc, precip = train_data_nov)
  
  data_test_gp = data.frame(test_loc, precip = test_data_nov)
  
  coordinates(data_train_gp) = ~x+y
  coordinates(data_test_gp) = ~x+y
  
  #plot(variogram(RR~1, data_train_gp))
  #Estimate parameters by maximum likelihood:
  ml <- likfit(data=data_train_gp$precip,coords=train_loc,
               fix.nugget=F,cov.model="exponential", 
               ini = c(1500, 266639),nugget=18)
  
  
  pred<-krige.conv(data= data_train_gp$precip, coords= train_loc, locations=test_loc,
                   krige=krige.control(cov.model="exponential", type.krige = "SK", beta = mean(train_data_nov),
                                       cov.pars=c(ml$sigmasq ,ml$phi),
                                       nugget= ml$nugget))
  
  parallel_gp[[index]] = list(tv = test_data_nov, predicted = pred$predict, var = pred$krige.var)
  
}


#### rmse, mad, pit and crps: gaussian prediction
library(scoringRules)

pit_RR_g =  vector(mode = "list", length = n_iters)
rmse_g= matrix(0, nrow = n_iters, ncol = length(test_data_list[[1]]$RR))
mad_g = matrix(0, nrow = n_iters, ncol = length(test_data_list[[1]]$RR))

crps_iter = vector();
crpsF_y_g = matrix(0, nrow = n_iters, ncol = length(test_data_list[[1]]$RR))
length_ci = vector();length_ci_g= vector(); alpha=0.1
for (i in 1:n_iters) {
  ### generate data from predictive distribution at all test locations
  tv = parallel_gp[[i]]$tv
  predicted = parallel_gp[[i]]$predicted
  var = parallel_gp[[i]]$var
  
  for (j  in 1:length(tv)) {
    pit_RR_g[[i]][j] = pnorm(tv[j], mean= predicted[j], sd = sqrt(var[j]))
    crps_iter[j] =   crps(y = tv[j], family = "normal",  mean= predicted[j], sd = sqrt(var[j]))
    
    lower_ci = qnorm(alpha/2, mean= predicted[j], sd = sqrt(var[j]))
    upper_ci = qnorm(1-alpha/2, mean= predicted[j], sd = sqrt(var[j]))
    length_ci[j] = upper_ci - lower_ci
  }
  crpsF_y_g[i,] =  crps_iter
  rmse_g[i,] = (tv-predicted)^2 ## RMSE
  mad_g[i,] = abs(tv-predicted) ## MAD
  length_ci_g[i] = mean(length_ci)
}
hist(unlist(pit_RR_g), freq=F, xlab = "Probability Integral Transform", main=paste("Gaussian "))
abline(h=1)    
mean(apply(rmse_g, 2, function(x) sqrt(mean(x))))
mean(apply(mad_g, 2, median))

mean(apply(crpsF_y_g,2,median))
mean(apply(crpsF_y_g,2,mean))
mean(length_ci_g)


#plot( parallel_gp[[5]]$tv, parallel_gp[[5]]$predicted)




### log transformation for Gaussian process

parallel_gp_log = list()
for (index in 1:n_iters) {
  
  train_data = train_data_list[[index]]
  test_data = test_data_list[[index]]
  
  train_data_nov_log = log(train_data[,"RR"]); train_loc = cbind(x = train_data$x, y= train_data$y)
  
  #### predict cdfs on the test data
  test_data_nov =  test_data[,"RR"]; test_loc = cbind(x= test_data$x, y= test_data$y)
  
  data_train_gp = data.frame(train_loc, precip = train_data_nov_log)
  data_test_gp = data.frame(test_loc, precip = test_data_nov)
  
  coordinates(data_train_gp) = ~x+y
  coordinates(data_test_gp) = ~x+y
  
  #plot(variogram(RR~1, data_train_gp))
  #Estimate parameters by maximum likelihood:
  ml <- likfit(data=data_train_gp$precip,coords=train_loc,
               fix.nugget=F,cov.model="exponential", 
               ini = c(1.5, 266639),nugget=0.08)
  
  
  pred<-krige.conv(data= data_train_gp$precip, coords= train_loc, locations=test_loc,
   krige=krige.control(cov.model="exponential", type.krige = "SK", beta = mean(train_data_nov_log),
                                       cov.pars=c(ml$sigmasq ,ml$phi),
                                       nugget= ml$nugget))
  
  ## back transforming the prdeiction and sd
  parallel_gp_log[[index]]=list(tv = test_data_nov, predicted = exp(pred$predict + pred$krige.var/2),
                                var = exp(2*pred$predict +pred$krige.var)*(exp(pred$krige.var)-1) )
  
}


#### rmse, mad, pit and crps: gaussian prediction
library(scoringRules)

pit_RR_g_log =  vector(mode = "list", length = n_iters)
rmse_g_log= matrix(0, nrow = n_iters, ncol = length(test_data_list[[1]]$RR))
mad_g_log = matrix(0, nrow = n_iters, ncol = length(test_data_list[[1]]$RR))

crps_iter = vector();
crpsF_y_g_log = matrix(0, nrow = n_iters, ncol = length(test_data_list[[1]]$RR))
length_ci = vector();length_ci_g= vector(); alpha=0.1
for (i in 1:n_iters) {
  ### generate data from predictive distribution at all test locations
  tv = parallel_gp_log[[i]]$tv
  predicted = parallel_gp_log[[i]]$predicted
  var = parallel_gp_log[[i]]$var
  
  for (j  in 1:length(tv)) {
    pit_RR_g_log[[i]][j] = pnorm(tv[j], mean= predicted[j], sd = sqrt(var[j]))
    crps_iter[j] =   crps(y = tv[j], family = "normal",  mean= predicted[j], sd = sqrt(var[j]))
    
    lower_ci = qnorm(alpha/2, mean= predicted[j], sd = sqrt(var[j]))
    upper_ci = qnorm(1-alpha/2, mean= predicted[j], sd = sqrt(var[j]))
    length_ci[j] = upper_ci - lower_ci
  }
  crpsF_y_g_log[i,] =  crps_iter
  rmse_g_log[i,] = (tv-predicted)^2 ## RMSE
  mad_g_log[i,] = abs(tv-predicted) ## MAD
  length_ci_g[i] = mean(length_ci)
}
hist(unlist(pit_RR_g_log), freq=F, xlab = "Probability Integral Transform",
     main=paste("Gaussian with transformation "), breaks = 10, ylim=c(0,3))
abline(h=1)    
mean(apply(rmse_g_log, 2, function(x) sqrt(mean(x))))
mean(apply(mad_g_log, 2, median))

mean(apply(crpsF_y_g_log,2,median))
mean(apply(crpsF_y_g_log,2,mean))
mean(length_ci_g)



### Indicator kriging


##################### Applying copula indicator kriging ######################3

library(copula)
library(fields)
library(sp)
library(mvtnorm)


#############################################################
#   Use Gaussian copula to estimate the covariance matrix   #
#                  Nonparametric marginal                   #
#############################################################

###Function: Pseudo neg logLikelihood of Gaussian copula
logL.Gcopula.Exp <- function(par,Y, coord){
  #For now, assume exponential correlation structure
  #so par is a scalar. 
  #Need modify the function to allow users to specify correlation structure.
  
  #Y: observed data, n-dimensional vector
  #coord: the coordinates (n*2 matrix)
  #par: theta in the exponential correlation function
  
  n = length(Y)
  
  # distance matrix
  dm <- as.matrix(dist(coord))
  
  # correlation matrix (exponential)
  R.m <- exp(-dm/par)
  params = R.m[lower.tri(R.m)]
  
  #myCop <- normalCopula(param=params, dim = n, dispstr = "un")  ## dim of coupla is equal to the number of spatial locations
  u = ecdf(Y)(Y)  ### extract the empirical distibution functions value
  
  #truncate u to avoid u equals exactly zero or 1
  u = pmax(0.001,pmin(u, 0.999))
  
  #use true marginal
  #u = pnorm(Y)
  
  logl <- dmvnorm(qnorm(u), mean=rep(0,n), sigma=R.m,log = TRUE) - sum(dnorm(qnorm(u), log=TRUE)) # log-density for X2(s) (denominator)
  
  #logl = log(dCopula(u, myCop))
  return(-logl)
}





### thresholds
#qseq = c(seq(0.05,0.5,length.out = 5),0.55, 0.6,0.65, seq(0.7, 0.975, length.out = 7)) ## good results

qseq = c(0.05, seq(0.1,0.7, by=0.1), seq(0.75,0.975, length.out = 9))
## choosing threshold
th = quantile(data_nov$RR, qseq )

#summary(data_nov$RR)

# Plot empirical cdf and thresholds
plot(ecdf(data_nov$RR), main="Global ECDF of RR conc.")
abline(v=th, lty=2)


start.time = Sys.time()
library(doParallel)
registerDoParallel(cores=50)
parallel = list(); tv_RR = list()
for (index in 1:n_iters) {
  
  train_data = train_data_list[[index]]
  test_data = test_data_list[[index]]
  
  train_data_nov = train_data[,"RR"]
  train_loc = cbind(x = train_data$x, y= train_data$y)
  
  th = quantile(train_data_nov, qseq )
  #plot(ecdf(train_data$RR), main="Global ECDF of RR conc.")
  #abline(v=th, lty=2)
  
  Y = train_data_nov
  coords = train_loc ### converting distances from m to km
  
  
  theta.hat <- optim(par = 1900,logL.Gcopula.Exp,Y=Y,coord=coords, method="L-BFGS-B",
                     lower=0.01,upper=600000)$par
  
  
  
  #qseq = seq(0.1,0.9, length.out = 9)
  #th = quantile(Y, qseq) 
  
  dist_mat =  as.matrix(dist(coords))
  h = dist_mat[lower.tri(dist_mat)]
  
  totC<-array(NA,dim = c(dim(dist_mat)[1],dim(dist_mat)[2],length(th),length(th)))
  
  myCop = vector(mode = "list", length = length(h))
  for (k in 1:length(h)) {
    myCop[[k]] =  normalCopula(param=exp(-h[k]/theta.hat), dim = 2, dispstr = "un")
  }
  
  b=rep(0, length(h))
  
  ## i from 
  e = rep(1:length(th), times=1:length(th))
  f = vector(); for(j in 1:length(th)) {f = c(f, 1:j)}
  ij = cbind(e,f)
  
  C = matrix(NA, nrow = nrow(dist_mat), ncol=ncol(dist_mat))
  Cij= foreach(z=1:dim(ij)[1]) %dopar% {    
    
    i = ij[z,1]
    j = ij[z,2] 
    #rhat: estimated correlation based on copula
    uhat_i = mean(Y< th[i]) # E[I(Y(s)<th_i)] = P[Y(s)<th_i] = F_Y(th_i)
    uhat_j = mean(Y< th[j])
    
    
    for(k in 1:length(h))
    {
      #estimated based on copula
      #myCop <- normalCopula(param=exp(-h[k]/theta.hat), dim = 2, dispstr = "un")
      b[k] = pCopula(c(uhat_i,uhat_j), myCop[[k]]) #P(U1<u, U2<u)
    }
    Cov = b - uhat_i*uhat_j
    myCop_diag <- normalCopula(param=exp(-0/theta.hat), dim = 2, dispstr = "un")
    Cov_diag = pCopula(c(uhat_i,uhat_j), myCop_diag) - uhat_i*uhat_j
    
    
    ## Construct the covariace matrix
    C[lower.tri(C)] <- Cov #put elements in lower triangular
    diag(C) = Cov_diag
    C = t(C) ## make it upper traingular
    C[lower.tri(C)] <- Cov
    
    
    totC[,,i,j] = C
  }
  
  
  #totC[,,2,1]=0;  totC[,,3,1]=0; totC[,,3,2]=0
  
  for (z in 1:dim(ij)[1]) {
    i = ij[z,1]; j = ij[z,2]
    totC[,,i,j] = Cij[[z]]
  }
  
  myC <-matrix(NA,nrow=dim(dist_mat)[1]*length(th),ncol=dim(dist_mat)[2]*length(th))
  
  for(i in 1:length(th))
  {
    for(j in 1:i)
    {
      myC[(dim(dist_mat)[1]*(i-1)+1):(dim(dist_mat)[1]*i),(dim(dist_mat)[2]*(j-1)+1):(dim(dist_mat)[2]*j)] = totC[,,i,j]
    }
  }
  
  for(i in 1:(length(th)-1))
  {
    for(j in (i+1):length(th))
    {
      myC[(dim(dist_mat)[1]*(i-1)+1):(dim(dist_mat)[1]*i),(dim(dist_mat)[2]*(j-1)+1):(dim(dist_mat)[2]*j)] = t(totC[,,j,i])
    }
  }
  
  
  
  #### predict cdfs on the test data
  test_data_nov =  test_data[,"RR"]
  test_loc = cbind(x= test_data$x, y= test_data$y)
  
  
  F_ind_list = list()
  inv_C = solve(myC)
  
  #dim(grd)[1]
  F_ind_list = foreach(l=1:dim(test_loc)[1]) %dopar% {
    
    newloc = data.frame(x= test_loc[l,1] ,y= test_loc[l,2]) 
    crosslocs<-rdist(newloc,train_loc)
    
    totC12<-array(NA,dim = c(dim(crosslocs)[1],dim(crosslocs)[2],length(th),length(th)))
    
    myCop = vector(mode="list", length = length(crosslocs))
    b = rep(0, length(crosslocs))
    
    for(k in 1:length(crosslocs))
    {
      myCop[[k]] <- normalCopula(param=exp(-crosslocs[k]/theta.hat), dim = 2, dispstr = "un")
    }
    
    for (i in 1:length(th)) {
      for (j in 1:length(th)) {
        
        uhat_i = mean(Y< th[i]) # E[I(Y(s)<th_i)] = P[Y(s)<th_i] = F_Y(th_i)
        uhat_j = mean(Y< th[j])
        
        for(k in 1:length(crosslocs))
        {
          #estimated based on copula
          #myCop <- normalCopula(param=exp(-crosslocs[k]/theta.hat), dim = 2, dispstr = "un")
          b[k] = pCopula(c(uhat_i,uhat_j), myCop[[k]]) #P(U1<u, U2<u)
        }
        
        Cov = b - uhat_i*uhat_j
        
        totC12[,,i,j] = Cov
        
        
      }
    }
    
    
    myC12 <-matrix(NA,nrow=dim(crosslocs)[1]*length(th),ncol=dim(crosslocs)[2]*length(th))
    
    
    for(i in 1:length(th))
    {
      for(j in 1:length(th))
      {
        myC12[(dim(crosslocs)[1]*(i-1)+1):(dim(crosslocs)[1]*i),(dim(crosslocs)[2]*(j-1)+1):(dim(crosslocs)[2]*j)] = totC12[,,i,j]
      }
    }
    
    
    
    ################################### Indicator variable data ################################
    data_ind = list(); mu_ind= vector()
    data_all = NULL
    for (i in 1:length(th)) {
      data_ind[[i]] = as.numeric(Y< th[i])
      data_all = c(data_all,data_ind[[i]])
 
    }
    mu_ind = qseq
    mu_all = rep(mu_ind, each= length(Y))

    ind_krig = mu_ind+ myC12%*%inv_C%*%(data_all - mu_all)
    
    F_ind = ind_krig
    
    
    ## Fix values at extremes, all values should lie between 0 and 1
    for (i in 1:length(F_ind)) {
      
      if(F_ind[i]<0){F_ind[1:i]=0}
      if(F_ind[i]>=1) {F_ind[i:length(F_ind)] = 1}
      
    }
    
    ############ Extrapolate to 0 and 1  ############
    xseq = th
    if(F_ind[length(xseq)]==1)
    {F_ind=F_ind} else { F_ind = c(F_ind,1)
    xseq = c(xseq, (max(data_nov$RR)+xseq[length(xseq)])/2) }
    
    if(F_ind[1]==0)
    {F_ind=F_ind} else {F_ind= c(0,F_ind)
    xseq = c( (min(data_nov$RR)+xseq[1])/2 ,xseq)}
    
    
    ## Order relation correction (Carr 1994)
    ##upward downward
    
    for (i in 1:(length(F_ind)-1)) {
      if(F_ind[i+1]< F_ind[i]){
        m = which(F_ind[(i+1):length(F_ind)]> F_ind[i])[1]
        
        for (j in (i+1):(i+m-1)) {
          F_ind[j] = F_ind[i] + ((F_ind[i+m]-F_ind[i])*(j-i))/m
        }
      }
    } 
    
    #F_ind_list[[l]] = cbind(xseq, F_ind)
    return(cbind(xseq, F_ind))
    
  }
  
  
  
  parallel[[index]]= F_ind_list
  tv_RR[[index]] = test_data_nov
  
}

end.time = Sys.time()
end.time  - start.time



####### Smooth distribution using the global distribution and rescaling
#lines(smoothF_RR[[1]][[18]])
smoothglobal = function(Fx, globalF){
  
  ecdf_RR = ecdf(globalF)
  ecdf_x = environment(ecdf_RR)$x; ecdf_y = environment(ecdf_RR)$y 
  
  Fsmooth= vector(); xsmooth= vector()
  for (t in 1:(nrow(Fx)[1]-1)  ) {
    
    xseq_1 = Fx[,1][t]; xseq_2 = Fx[,1][t+1]
    
    ind_x = (xseq_1 < ecdf_x) & (ecdf_x  <= xseq_2)
    
    a = Fx[,2][t]; b = Fx[,2][t+1]
    
    
    ff = c(ecdf_y[ind_x])
    if( (max(ff)- min(ff)) ==0 ) {fs=a} else fs = ((ff- min(ff))/ (max(ff)- min(ff)))*(b-a)+a
    x_ff = c(ecdf_x[ind_x])
    
    
    Fsmooth = c(Fsmooth, fs)
    xsmooth = c(xsmooth, x_ff)
  }
  Fsmooth =  c(Fsmooth, b)
  xsmooth = c(xsmooth, xseq_2)
  return(cbind(xsmooth, Fsmooth))
}

smoothtry = smoothglobal(Fx = parallel[[3]][[7]], globalF = data_nov$RR )
plot(parallel[[3]][[7]], ylab="F")
lines(smoothtry,  col="grey" )


L=200; 
crpsF_y = matrix(0, nrow = n_iters, ncol = length(test_data_list[[1]]$RR));  
smoothF_RR = list()
pit_RR =  vector(mode = "list", length = n_iters) ; rmse_RR=vector(); rmse_RR_mean = vector();

rmse_RR = matrix(0, nrow = n_iters, ncol = length(test_data_list[[1]]$RR))
mad_RR = matrix(0, nrow = n_iters, ncol = length(test_data_list[[1]]$RR))

x_seq = seq(min(data_nov$RR), max(data_nov$RR), length.out = 10000) #RR
x_mid = x_seq[-length(x_seq)] + diff(x_seq)/2
delta_x = diff(x_seq)[1]  

for (index in 1:n_iters) {
  smoothF_RR[[index]] =   lapply (parallel[[index]], function(a) smoothglobal(Fx=a, globalF = data_nov$RR) )
  median_pred_RR = as.numeric(lapply(smoothF_RR[[index]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
 
  quantile_pred_RR_l = matrix(0, nrow=L, ncol = length(median_pred_RR))
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_RR_l[l,] = as.numeric(lapply(smoothF_RR[[index]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  

  mean_pred_RR = apply(quantile_pred_RR_l, 2, mean)
  tv = tv_RR[[index]]
  rmse_RR[index,] = (tv-median_pred_RR)^2 ## RMSE
  rmse_RR_mean[index] = sqrt(mean((tv-mean_pred_RR)^2)) ## RMSE
  mad_RR[index,] = abs(tv-median_pred_RR) ## MAD
  
  crpsF_y_iter= vector()
  for (j  in 1:length(tv)) {
    
    pit_RR[[index]][j] =   mean(quantile_pred_RR_l[,j]<=tv[j])
    f_x = function(x,y) { {mean(quantile_pred_RR_l[,j]<=x) - ifelse( y <= x,1,0)}^2 }
    fx_mid = sapply(x_mid, function(x) f_x(x, tv[j]))
    crpsF_y_iter[j] = sum(delta_x*fx_mid, na.rm = T)
  }
  crpsF_y[index,] =  crpsF_y_iter
}

hist(unlist(pit_RR), freq=F, xlab = "Proability Integral Transform", main=paste("Gaussian Copula "))
abline(h=1)
mean(apply(crpsF_y,2,median))
mean(rmse_RR); mean(mad_RR)

mean(apply(rmse_RR, 2, function(x) sqrt(mean(x))))
mean(apply(mad_RR, 2, median))

mean(rmse_RR_mean)
mean(length_ci_gcop)

sd(crpsF_y)
sd(crpsF_y_g)


length(qseq)

###################################################################################################
###############################################################################################

library(gstat)
parallel_var= list(); tv_RR_var = list()
for (index in 1:n_iters) {
  
  
  train_data = train_data_list[[index]]
  test_data = test_data_list[[index]]
  
  train_data_nov = train_data[,"RR"]
  train_loc = cbind(x = train_data$x, y= train_data$y)
  
  th = quantile(train_data_nov, qseq )
  
  data_train_var = data.frame(train_loc,RR = train_data_nov)
  
  #### predict cdfs on the test data
  test_data_RR = test_data[,"RR"]
  test_loc = cbind(x = test_data$x, y= test_data$y)
  
  data_test_var = data.frame(test_loc, RR = test_data_RR)
  
  coordinates(data_train_var) = ~x+y
  coordinates(data_test_var) = ~x+y
  
  ### indicator cokriging
  
  quartz = th
  # These are the actual values of the quantiles
  
  RR.i <- gstat(id = "RR1", formula = I(RR < quartz[1]) ~ 1, data = data_train_var, 
                  nmax = 7, beta = qseq[1], set = list(order = 4, zero = 1e-05))
  RR.i <- gstat(RR.i, "RR2", formula = I(RR < quartz[2]) ~ 1, data = data_train_var, 
                  nmax = 7, beta = qseq[2])
  RR.i <- gstat(RR.i, "RR3", formula = I(RR < quartz[3]) ~ 1, data = data_train_var, 
                  nmax = 7, beta = qseq[3])
  RR.i <- gstat(RR.i, "RR4", formula = I(RR < quartz[4]) ~ 1, data = data_train_var, 
                  nmax = 7, beta = qseq[4])
  RR.i <- gstat(RR.i, "RR5", formula = I(RR < quartz[5]) ~ 1, data = data_train_var, 
                  nmax = 7, beta = qseq[5])
  RR.i <- gstat(RR.i, "RR6", formula = I(RR < quartz[6]) ~ 1, data = data_train_var, 
                  nmax = 7, beta = qseq[6])
  RR.i <- gstat(RR.i, "RR7", formula = I(RR < quartz[7]) ~ 1, data = data_train_var, 
                  nmax = 7, beta = qseq[7])
  RR.i <- gstat(RR.i, "RR8", formula = I(RR < quartz[8]) ~ 1, data = data_train_var, 
                  nmax = 7, beta = qseq[8])
  RR.i <- gstat(RR.i, "RR9", formula = I(RR < quartz[9]) ~ 1, data = data_train_var, 
                  nmax = 7, beta = qseq[9])
  RR.i <- gstat(RR.i, "RR10", formula = I(RR < quartz[10]) ~ 1, data = data_train_var, 
                  nmax = 7, beta = qseq[10])
  RR.i <- gstat(RR.i, "RR11", formula = I(RR < quartz[11]) ~ 1, data = data_train_var, 
                  nmax = 7, beta = qseq[11])
  RR.i <- gstat(RR.i, "RR12", formula = I(RR < quartz[12]) ~ 1, data = data_train_var, 
                  nmax = 7, beta = qseq[12])
  RR.i <- gstat(RR.i, "RR13", formula = I(RR < quartz[13]) ~ 1, data = data_train_var, 
                  nmax = 7, beta = qseq[13])
  RR.i <- gstat(RR.i, "RR14", formula = I(RR < quartz[14]) ~ 1, data = data_train_var, 
                  nmax = 7, beta = qseq[14])
  RR.i <- gstat(RR.i, "RR15", formula = I(RR < quartz[15]) ~ 1, data = data_train_var, 
                  nmax = 7, beta = qseq[15])
  RR.i <- gstat(RR.i, "RR16", formula = I(RR < quartz[16]) ~ 1, data = data_train_var, 
                nmax = 7, beta = qseq[16])
  RR.i <- gstat(RR.i, "RR17", formula = I(RR < quartz[17]) ~ 1, data = data_train_var, 
                nmax = 7, beta = qseq[17])
  
  RR.i <- gstat(RR.i, model = vgm(0.1, "Exp", 300000, 0.1), fill.all = T)
  
  
  RR.vg <- variogram(RR.i)
  #plot(RR.vg)
  RR.quartfit = fit.lmc(RR.vg, RR.i)
  
  Fvar_list= list()
  # now do the co-kriging
  for (l in 1:dim(data_test_var@coords)[1]) {
    
    cokrig.RRquart <- predict(RR.quartfit, data_test_var[l,data_test_var$coordinates])
    
    Fvar = c(cokrig.RRquart@data$RR1.pred,cokrig.RRquart@data$RR2.pred,cokrig.RRquart@data$RR3.pred,cokrig.RRquart@data$RR4.pred,cokrig.RRquart@data$RR5.pred,cokrig.RRquart@data$RR6.pred,cokrig.RRquart@data$RR7.pred,cokrig.RRquart@data$RR8.pred,
             cokrig.RRquart@data$RR9.pred,cokrig.RRquart@data$RR10.pred,cokrig.RRquart@data$RR11.pred,cokrig.RRquart@data$RR12.pred,cokrig.RRquart@data$RR13.pred,cokrig.RRquart@data$RR14.pred,cokrig.RRquart@data$RR15.pred,
             cokrig.RRquart@data$RR16.pred, cokrig.RRquart@data$RR17.pred)
    
    
    ############ Extrapolate to 0 and 1  ############
    xseq= quartz
    if(Fvar[length(quartz)]==1)
    {Fvar=Fvar} else {Fvar = c(Fvar,1)
    xseq = c(xseq,(max(data_train_var$RR)+xseq[length(xseq)])/2)}
    
    if(Fvar[1]==0)
    {Fvar=Fvar} else {Fvar= c(0,Fvar)
    xseq = c( (min(data_train_var$RR)+xseq[1])/2 ,xseq)}
    
    Fvar_list[[l]] = cbind(xseq, Fvar)
    
  }
  
  #parallel_var[[index]] = list(test_data_RR, median_pred_var)
  parallel_var[[index]] = Fvar_list
  tv_RR_var[[index]] = test_data_RR
  
}



####### Smooth distribution using the global distribution and rescaling

crpsF_y_var = matrix(0, nrow = n_iters, ncol = length(test_data_list[[1]]$RR)); crpsF_y_iter= vector(); smoothF_RR = list()
pit_RR_var =  vector(mode = "list", length = n_iters) ; rmse_RR_var=vector();
rmse_RR_mean_var = vector()
L=200; 
length_ci_var= vector()
mad_RR_var = matrix(0, nrow = n_iters, ncol = length(test_data_list[[1]]$RR))

x_seq = seq(min(data_nov$RR), max(data_nov$RR), length.out = 10000) #RR
x_mid = x_seq[-length(x_seq)] + diff(x_seq)/2
delta_x = diff(x_seq)[1]  

for (index in 1:n_iters) {
  smoothF_RR[[index]] =   lapply (parallel_var[[index]], function(a) smoothglobal(Fx=a, globalF = data_nov$RR) )
  median_pred_RR = as.numeric(lapply(smoothF_RR[[index]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
  quantile_pred_RR_l = matrix(0, nrow=L, ncol = length(median_pred_RR))
  
  lower_ci = as.numeric(lapply(smoothF_RR[[index]], function(a) {min(a[,1][which(a[,2]>=alpha/2)])}))
  upper_ci = as.numeric(lapply(smoothF_RR[[index]], function(a) {min(a[,1][which(a[,2]>=1-alpha/2)])}))
  length_ci_var[index] = mean(upper_ci-lower_ci)
  
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_RR_l[l,] = as.numeric(lapply(smoothF_RR[[index]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  
  mean_pred_RR = apply(quantile_pred_RR_l, 2, mean)
  tv = tv_RR_var[[index]]
  rmse_RR_var[index] = sqrt(mean((tv-median_pred_RR)^2)) ## RMSE
  rmse_RR_mean_var[index] = sqrt(mean((tv-mean_pred_RR)^2)) ## RMSE
  mad_RR_var[index,] = abs(tv-median_pred_RR) ## MAD
  
  for (j  in 1:length(tv)) {
    
    pit_RR_var[[index]][j] =   mean(quantile_pred_RR_l[,j]<=tv[j])
    f_x = function(x,y) { {mean(quantile_pred_RR_l[,j]<=x) - ifelse( y <= x,1,0)}^2 }
    fx_mid = sapply(x_mid, function(x) f_x(x, tv[j]))
    crpsF_y_iter[j] = sum(delta_x*fx_mid, na.rm = T)
  }
  crpsF_y_var[index,] =  crpsF_y_iter
}


hist(unlist(pit_RR_var), freq=F, main=paste("Variogram"), xlab="Probability Integral Transform")
abline(h=1)
mean(crpsF_y_var, na.rm = T)
mean(apply(crpsF_y_var,2,median))
mean(rmse_RR_var)
mean(apply(mad_RR_var,2, median))
mean(rmse_RR_mean_var)
mean(length_ci_var)



#############################################################
#   Use t-copula to estimate the covariance matrix   #
#                  Nonparametric marginal                   #
#############################################################

###Function: Pseudo neg logLikelihood of t-copula
logL.tcopula.Exp <- function(par,Y, coord){
  

  n = length(Y)
  
  # distance matrix
  dm <- as.matrix(dist(coord))
  
  # correlation matrix (exponential)
  R.m <- exp(-dm/par[1])
  params = R.m[lower.tri(R.m)]
  
  #myCop <- normalCopula(param=params, dim = n, dispstr = "un")  ## dim of coupla is equal to the number of spatial locations
  u = ecdf(Y)(Y)  ### extract the empirical distibution functions value
  
  #truncate u to avoid u equals exactly zero or 1
  u = pmax(0.001,pmin(u, 0.999))
  
  #use true marginal
  #u = pnorm(Y)
  #nu= par[2] #
  nu= 4
  logl <- dmvt(qt(u, df = nu), df=nu, sigma=R.m,log = TRUE) - sum(dt(qt(u, df=nu), df=nu, log=TRUE)) # log-density for X2(s) (denominator)
  
  #logl = log(dCopula(u, myCop))
  return(-logl)
}




## if we estimate df then it reaches close to gaussina copula results as estimate df is very high
start.time = Sys.time()
library(doParallel)
registerDoParallel(cores=50)
parallel_t = list(); tv_RR = list()
for (index in 1:n_iters) {
  
  train_data = train_data_list[[index]]
  test_data = test_data_list[[index]]
  
  train_data_nov = train_data[,"RR"]
  train_loc = cbind(x = train_data$x, y= train_data$y)
  
  th = quantile(train_data_nov, qseq )
  
  Y = train_data_nov
  coords = train_loc ### converting distances from m to km
  
  
  theta.hat_vec <- optim(par = c(152635),logL.tcopula.Exp,Y=Y,coord=coords, method="L-BFGS-B",
                         lower=c(0.01),upper=c(10^6))$par
  
  theta.hat = theta.hat_vec[1]
  nu = 4#round(theta.hat_vec[2],0) #4
  
  #qseq = seq(0.1,0.9, length.out = 9)
  #th = quantile(Y, qseq) 
  
  dist_mat =  as.matrix(dist(coords))
  h = dist_mat[lower.tri(dist_mat)]
  
  totC<-array(NA,dim = c(dim(dist_mat)[1],dim(dist_mat)[2],length(th),length(th)))
  
  myCop = vector(mode = "list", length = length(h))
  for (k in 1:length(h)) {
    myCop[[k]] <- tCopula(param=exp(-h[k]/theta.hat), dim = 2, dispstr = "un", df=nu)
  }
  
  b=rep(0, length(h))
  
  ## i from 
  e = rep(1:length(th), times=1:length(th))
  f = vector(); for(j in 1:length(th)) {f = c(f, 1:j)}
  ij = cbind(e,f)
  
  C = matrix(NA, nrow = nrow(dist_mat), ncol=ncol(dist_mat))
  Cij= foreach(z=1:dim(ij)[1]) %dopar% {    
    
    i = ij[z,1]
    j = ij[z,2] 
    #rhat: estimated correlation based on copula
    uhat_i = mean(Y< th[i]) # E[I(Y(s)<th_i)] = P[Y(s)<th_i] = F_Y(th_i)
    uhat_j = mean(Y< th[j])
    
    
    for(k in 1:length(h))
    {
      #estimated based on copula
      #myCop <- tCopula(param=exp(-h[k]/theta.hat), dim = 2, dispstr = "un", df=nu)
      b[k] = pCopula(c(uhat_i,uhat_j), myCop[[k]]) #P(U1<u, U2<u)
    }
    Cov = b - uhat_i*uhat_j
    myCop_diag <- tCopula(param=exp(-0/theta.hat), dim = 2, dispstr = "un", df=nu)
    Cov_diag = pCopula(c(uhat_i,uhat_j), myCop_diag) - uhat_i*uhat_j
    
    
    ## Construct the covariace matrix
    C[lower.tri(C)] <- Cov #put elements in lower triangular
    diag(C) = Cov_diag
    C = t(C) ## make it upper traingular
    C[lower.tri(C)] <- Cov
    
    
    totC[,,i,j] = C
  }
  
  
  #totC[,,2,1]=0;  totC[,,3,1]=0; totC[,,3,2]=0
  
  for (z in 1:dim(ij)[1]) {
    i = ij[z,1]; j = ij[z,2]
    totC[,,i,j] = Cij[[z]]
  }
  
  myC <-matrix(NA,nrow=dim(dist_mat)[1]*length(th),ncol=dim(dist_mat)[2]*length(th))
  
  for(i in 1:length(th))
  {
    for(j in 1:i)
    {
      myC[(dim(dist_mat)[1]*(i-1)+1):(dim(dist_mat)[1]*i),(dim(dist_mat)[2]*(j-1)+1):(dim(dist_mat)[2]*j)] = totC[,,i,j]
    }
  }
  
  for(i in 1:(length(th)-1))
  {
    for(j in (i+1):length(th))
    {
      myC[(dim(dist_mat)[1]*(i-1)+1):(dim(dist_mat)[1]*i),(dim(dist_mat)[2]*(j-1)+1):(dim(dist_mat)[2]*j)] = t(totC[,,j,i])
    }
  }
  
  
  
  #### predict cdfs on the test data
  test_data_nov =  test_data[,"RR"]
  test_loc = cbind(x= test_data$x, y= test_data$y)
  
  F_ind_list = list()
  inv_C = solve(myC)
  
  #dim(grd)[1]
  F_ind_list = foreach(l=1:dim(test_loc)[1]) %dopar% {
    
    newloc = data.frame(x= test_loc[l,1] ,y= test_loc[l,2]) 
    crosslocs<-rdist(newloc,train_loc)
    
    totC12<-array(NA,dim = c(dim(crosslocs)[1],dim(crosslocs)[2],length(th),length(th)))
    
    myCop = vector(mode="list", length = length(crosslocs))
    b = rep(0, length(crosslocs))
    
    for(k in 1:length(crosslocs))
    {
      myCop[[k]] <- tCopula(param=exp(-crosslocs[k]/theta.hat), dim = 2, dispstr = "un", df=nu)
    }
    
    
    for (i in 1:length(th)) {
      for (j in 1:length(th)) {
        
        uhat_i = mean(Y< th[i]) # E[I(Y(s)<th_i)] = P[Y(s)<th_i] = F_Y(th_i)
        uhat_j = mean(Y< th[j])
        
        b=vector()
        for(k in 1:length(crosslocs))
        {
          #estimated based on copula
          #myCop <- tCopula(param=exp(-crosslocs[k]/theta.hat), dim = 2, dispstr = "un", df=nu)
          b[k] = pCopula(c(uhat_i,uhat_j), myCop[[k]]) #P(U1<u, U2<u)
        }
        
        Cov = b - uhat_i*uhat_j
        
        totC12[,,i,j] = Cov
        
        
      }
    }
    
    
    myC12 <-matrix(NA,nrow=dim(crosslocs)[1]*length(th),ncol=dim(crosslocs)[2]*length(th))
    
    
    for(i in 1:length(th))
    {
      for(j in 1:length(th))
      {
        myC12[(dim(crosslocs)[1]*(i-1)+1):(dim(crosslocs)[1]*i),(dim(crosslocs)[2]*(j-1)+1):(dim(crosslocs)[2]*j)] = totC12[,,i,j]
      }
    }
    
    
    
    ################################### Indicator variable data ################################
    data_ind = list(); mu_ind= vector()
    data_all = NULL
    for (i in 1:length(th)) {
      data_ind[[i]] = as.numeric(Y< th[i])
      data_all = c(data_all,data_ind[[i]])
      
      # mu_ind[i] = mean(data_ind[[i]])
    }
    mu_ind = qseq
    mu_all = rep(mu_ind, each= length(Y))
    
    # inv_C = solve(myC)  outside the loop
    #ind_krig = myC12%*%inv_C%*%data_all
    ind_krig = mu_ind+ myC12%*%inv_C%*%(data_all - mu_all)
    
    F_ind = ind_krig
    
    
    ## Fix values at extremes, all values should lie between 0 and 1
    for (i in 1:length(F_ind)) {
      
      if(F_ind[i]<0){F_ind[1:i]=0}
      if(F_ind[i]>=1) {F_ind[i:length(F_ind)] = 1}
      
    }
    
    ############ Extrapolate to 0 and 1  ############
    xseq = th
    if(F_ind[length(xseq)]==1)
    {F_ind=F_ind} else { F_ind = c(F_ind,1)
    xseq = c(xseq, (max(data_nov$RR)+xseq[length(xseq)])/2) }
    
    if(F_ind[1]==0)
    {F_ind=F_ind} else {F_ind= c(0,F_ind)
    xseq = c( (min(data_nov$RR)+xseq[1])/2 ,xseq)}
    
    
    ## Order relation correction (Carr 1994)
    ##upward downward
    
    for (i in 1:(length(F_ind)-1)) {
      if(F_ind[i+1]< F_ind[i]){
        m = which(F_ind[(i+1):length(F_ind)]> F_ind[i])[1]
        
        for (j in (i+1):(i+m-1)) {
          F_ind[j] = F_ind[i] + ((F_ind[i+m]-F_ind[i])*(j-i))/m
        }
      }
    } 
    
    #F_ind_list[[l]] = cbind(xseq, F_ind)
    return(cbind(xseq, F_ind))
    
  }
  
  
  
  parallel_t[[index]]= F_ind_list
  tv_RR[[index]] = test_data_nov
  
}

end.time = Sys.time()
end.time  - start.time





####### Smooth distribution using the global distribution and rescaling


L=200; 
crpsF_y_t =matrix(0, nrow = n_iters, ncol = length(test_data_list[[1]]$RR));  smoothF_RR = list()
pit_RR_t =  vector(mode = "list", length = n_iters) ; 
rmse_RR_t=vector(); rmse_RR_mean_t = vector()
length_ci_tcop= vector()
mad_RR_t = matrix(0, nrow = n_iters, ncol = length(test_data_list[[1]]$RR))

x_seq = seq(min(data_nov$RR), max(data_nov$RR), length.out = 10000) #RR
x_mid = x_seq[-length(x_seq)] + diff(x_seq)/2
delta_x = diff(x_seq)[1]  


for (index in 1:n_iters) {
  smoothF_RR[[index]] =   lapply (parallel_t[[index]], function(a) smoothglobal(Fx=a, globalF = data_nov$RR) )
  median_pred_RR = as.numeric(lapply(smoothF_RR[[index]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
  quantile_pred_RR_l = matrix(0, nrow=L, ncol = length(median_pred_RR))
  
  lower_ci = as.numeric(lapply(smoothF_RR[[index]], function(a) {min(a[,1][which(a[,2]>=alpha/2)])}))
  upper_ci = as.numeric(lapply(smoothF_RR[[index]], function(a) {min(a[,1][which(a[,2]>=1-alpha/2)])}))
  length_ci_tcop[index] = mean(upper_ci-lower_ci)
  
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_RR_l[l,] = as.numeric(lapply(smoothF_RR[[index]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  
  mean_pred_RR = apply(quantile_pred_RR_l, 2, mean)
  tv = tv_RR[[index]]
  rmse_RR_t[index] = sqrt(mean((tv-median_pred_RR)^2)) ## RMSE
  rmse_RR_mean_t[index] = sqrt(mean((tv-mean_pred_RR)^2)) ## RMSE
  mad_RR_t[index,] = abs(tv-median_pred_RR) ## MAD
  
  crpsF_y_iter= vector()
  for (j  in 1:length(tv)) {
    
    pit_RR_t[[index]][j] =   mean(quantile_pred_RR_l[,j]<=tv[j])
    f_x = function(x,y) { {mean(quantile_pred_RR_l[,j]<=x) - ifelse( y <= x,1,0)}^2 }
    fx_mid = sapply(x_mid, function(x) f_x(x, tv[j]))
    crpsF_y_iter[j] = sum(delta_x*fx_mid, na.rm = T)
  }
  crpsF_y_t[index,] =  crpsF_y_iter
}

hist(unlist(pit_RR_t), freq=F, main=paste("t-Copula"), xlab="Probability Integral Transform")

abline(h=1)
mean(crpsF_y_t, na.rm = T)
mean(apply(crpsF_y_t,2,median))
mean(rmse_RR_t); 
mean(apply(mad_RR_t, 2, median))
mean(rmse_RR_mean_t)
mean(length_ci_tcop)

sd(crpsF_y_t)
sd(crpsF_y_g)




rbind(mean(apply(mad_g, 2, median)),mean(apply(mad_g_log, 2, median)),mean(apply(mad_RR_var, 2, median)),mean(apply(mad_RR, 2, median)),mean(apply(mad_RR_t, 2, median)))
rbind(mean(apply(crpsF_y_g, 2, mean)),mean(apply(crpsF_y_g_log, 2, mean)), 
      mean(apply(crpsF_y_var, 2, mean)),mean(apply(crpsF_y, 2, mean)),mean(apply(crpsF_y_t, 2, mean)))

rbind(mean(apply(crpsF_y_g, 2, median)), mean(apply(crpsF_y_g_log, 2, median)),
      mean(apply(crpsF_y_var, 2, median)),mean(apply(crpsF_y, 2, median)),mean(apply(crpsF_y_t, 2, median)))


##standard error
1.2533*rbind(mean(apply(mad_g, 2, sd)/sqrt(n_iters)),mean(apply(mad_g_log, 2, median)/sqrt(n_iters)),
      mean(apply(mad_RR_var, 2, median)/sqrt(n_iters)),mean(apply(mad_RR, 2, median)/sqrt(n_iters)),mean(apply(mad_RR_t, 2, median)/sqrt(n_iters)))
rbind(mean(apply(crpsF_y_g, 2, sd)/sqrt(n_iters)), mean(apply(crpsF_y_g_log, 2, median)/sqrt(n_iters)),
      mean(apply(crpsF_y_var, 2, sd)/sqrt(n_iters)),mean(apply(crpsF_y, 2, median)/sqrt(n_iters)),mean(apply(crpsF_y_t, 2, median)/sqrt(n_iters)))
1.2533*rbind(mean(apply(crpsF_y_g, 2, sd)/sqrt(n_iters)), mean(apply(crpsF_y_g_log, 2, median)/sqrt(n_iters)),
      mean(apply(crpsF_y_var, 2, sd)/sqrt(n_iters)),mean(apply(crpsF_y, 2, median)/sqrt(n_iters)),mean(apply(crpsF_y_t, 2, median)/sqrt(n_iters)))



####################### Figures for manuscript ############################

#  PIT

hist(unlist(pit_RR_g), freq=F, ylim=c(0,3), main=paste("Gaussian Kriging"), 
     xlab="Probability Integral Transform", breaks=10, cex.main=1.5, cex.lab=1.4)
abline(h=1)
hist(unlist(pit_RR_g_log), freq=F, ylim=c(0,3), main=paste("Gaussian Kriging (log transfromation)"), 
     xlab="Probability Integral Transform", breaks=10, cex.main=1.5, cex.lab=1.4)
abline(h=1)
hist(unlist(pit_RR_var), freq=F, ylim=c(0,3), main=paste("Varoigram"), 
     xlab="Probability Integral Transform", breaks=10, cex.main=1.5, cex.lab=1.4)
abline(h=1)
hist(unlist(pit_RR), freq=F, ylim=c(0,3), main=paste("Gaussian Copula"), 
     xlab="Probability Integral Transform", breaks=10, cex.main=1.5, cex.lab=1.4)
abline(h=1)
hist(unlist(pit_RR_t), freq=F, ylim=c(0,3), main=paste("t-Copula"), 
     xlab="Probability Integral Transform", breaks=10, cex.main=1.5, cex.lab=1.4)
abline(h=1)
hist(unlist(pit_RR_g_log), freq=F, ylim=c(0,3), main=paste("Gaussian with transfromation"), 
     xlab="Probability Integral Transform", breaks=10, cex.main=1.5, cex.lab=1.4)
abline(h=1)


## gplot
ggplot() + geom_histogram(aes(x=unlist(pit_RR_g), y=..density..), 
          breaks=hist(unlist(pit_RR_g), plot=F, breaks = 10)$breaks,  col="black")+ ylim(0,3)+
  xlab("Proability Integral Transform")+ ylab("Density")+
  ggtitle("GK")+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.8)),
        plot.title = element_text(size=rel(2), hjust = 0.5))+
  geom_hline(yintercept = 1)

ggplot() + geom_histogram(aes(x=unlist(pit_RR_g_log), y=..density..), 
            breaks=hist(unlist(pit_RR_g), plot=F, breaks = 10)$breaks,  col="black")+ ylim(0,3)+
  xlab("Proability Integral Transform")+ ylab("Density")+
  ggtitle("logGK")+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.8)),
        plot.title = element_text(size=rel(2), hjust = 0.5))+
  geom_hline(yintercept = 1)

ggplot() + geom_histogram(aes(x=unlist(pit_RR_var), y=..density..), 
                          breaks=hist(unlist(pit_RR_var), plot=F, breaks = 10)$breaks,  col="black")+ ylim(0,3)+
  xlab("Proability Integral Transform")+ ylab("Density")+
  ggtitle("VMIK")+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.5)),
        plot.title = element_text(size=rel(1.8), hjust = 0.5))+
  geom_hline(yintercept = 1)

ggplot() + geom_histogram(aes(x=unlist(pit_RR), y=..density..), 
                          breaks=hist(unlist(pit_RR), plot=F, breaks = 10)$breaks,  col="black")+ ylim(0,3)+
  xlab("Proability Integral Transform")+ ylab("Density")+
  ggtitle(expression(paste("CMIK"[G])))+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.8)),
        plot.title = element_text(size=rel(2), hjust = 0.5))+
  geom_hline(yintercept = 1)

ggplot() + geom_histogram(aes(x=unlist(pit_RR_t), y=..density..), 
                          breaks=hist(unlist(pit_RR_t), plot=F, breaks = 10)$breaks,  col="black")+ ylim(0,3)+
  xlab("Proability Integral Transform")+ ylab("Density")+
  ggtitle(expression(paste("CMIK"[t])))+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.8)),
        plot.title = element_text(size=rel(2), hjust = 0.5))+
  geom_hline(yintercept = 1)




















#### Applying the best method- tcopula #####

#### predict cdfs in a grid over study region
####Define the grid for interpolation
##Define the grid extent:

spain = map_data("world", region = "spain")

x.range <- range(spain$long)
y.range <- range(spain$lat)
grd <- as.matrix(expand.grid(x=seq(from=x.range[1], to=x.range[2], length.out = 100),
                             y=seq(from=y.range[1], to=y.range[2], length.out = 100) ))

spain_latlon_all=NULL
for (i in unique(spain$group)) {
  spain_group = subset(spain, group==i)
  spain_latlon = cbind(spain_group$long, spain_group$lat)
  spain_latlon_all = rbind(spain_latlon_all, spain_latlon, c(NA,NA))
}

library(mgcv)
pred_grd = grd[in.out(spain_latlon_all,grd),] 
dim(pred_grd)
plot(pred_grd)

pred_grd_m <- LongLatToUTM(pred_grd[,1],pred_grd[,2],zone = 30)




## Fit copula model to the data and predict in a grid to create map

Y = data_nov$RR
coords = loc_m ### converting distances from m to km


theta.hat_vec <- optim(par = c(1900),logL.tcopula.Exp,Y=Y,coord=coords, method="L-BFGS-B",
                       lower=0.01,upper=10^6)$par

theta.hat = theta.hat_vec[1]
nu = 4#round(theta.hat_vec[2],0)

#qseq = seq(0.1,0.9, length.out = 9)
#th = quantile(Y, qseq) 

dist_mat =  as.matrix(dist(coords))
h = dist_mat[lower.tri(dist_mat)]

totC<-array(NA,dim = c(dim(dist_mat)[1],dim(dist_mat)[2],length(th),length(th)))

myCop = vector(mode = "list", length = length(h))
for (k in 1:length(h)) {
  myCop[[k]] <- tCopula(param=exp(-h[k]/theta.hat), dim = 2, dispstr = "un", df=nu)
}

b=rep(0, length(h))

## i from 
e = rep(1:length(th), times=1:length(th))
f = vector(); for(j in 1:length(th)) {f = c(f, 1:j)}
ij = cbind(e,f)

C = matrix(NA, nrow = nrow(dist_mat), ncol=ncol(dist_mat))
Cij= foreach(z=1:dim(ij)[1]) %dopar% {    
  
  i = ij[z,1]
  j = ij[z,2] 
  #rhat: estimated correlation based on copula
  uhat_i = mean(Y< th[i]) # E[I(Y(s)<th_i)] = P[Y(s)<th_i] = F_Y(th_i)
  uhat_j = mean(Y< th[j])
  
  
  for(k in 1:length(h))
  {
    #estimated based on copula
    #myCop <- tCopula(param=exp(-h[k]/theta.hat), dim = 2, dispstr = "un", df=nu)
    b[k] = pCopula(c(uhat_i,uhat_j), myCop[[k]]) #P(U1<u, U2<u)
  }
  Cov = b - uhat_i*uhat_j
  myCop_diag <- tCopula(param=exp(-0/theta.hat), dim = 2, dispstr = "un", df=nu)
  Cov_diag = pCopula(c(uhat_i,uhat_j), myCop_diag) - uhat_i*uhat_j
  
  
  ## Construct the covariace matrix
  C[lower.tri(C)] <- Cov #put elements in lower triangular
  diag(C) = Cov_diag
  C = t(C) ## make it upper traingular
  C[lower.tri(C)] <- Cov
  
  
  totC[,,i,j] = C
}


#totC[,,2,1]=0;  totC[,,3,1]=0; totC[,,3,2]=0

for (z in 1:dim(ij)[1]) {
  i = ij[z,1]; j = ij[z,2]
  totC[,,i,j] = Cij[[z]]
}

myC <-matrix(NA,nrow=dim(dist_mat)[1]*length(th),ncol=dim(dist_mat)[2]*length(th))

for(i in 1:length(th))
{
  for(j in 1:i)
  {
    myC[(dim(dist_mat)[1]*(i-1)+1):(dim(dist_mat)[1]*i),(dim(dist_mat)[2]*(j-1)+1):(dim(dist_mat)[2]*j)] = totC[,,i,j]
  }
}

for(i in 1:(length(th)-1))
{
  for(j in (i+1):length(th))
  {
    myC[(dim(dist_mat)[1]*(i-1)+1):(dim(dist_mat)[1]*i),(dim(dist_mat)[2]*(j-1)+1):(dim(dist_mat)[2]*j)] = t(totC[,,j,i])
  }
}





F_ind_list = list()
inv_C = solve(myC)

#dim(grd)[1]
F_ind_list = foreach(l=1:dim(pred_grd_m)[1]) %dopar% {
  
  newloc = data.frame(x= pred_grd_m[l,1] ,y= pred_grd_m[l,2]) 
  crosslocs<-rdist(newloc,loc_m)
  
  totC12<-array(NA,dim = c(dim(crosslocs)[1],dim(crosslocs)[2],length(th),length(th)))
  
  myCop = vector(mode="list", length = length(crosslocs))
  b = rep(0, length(crosslocs))
  
  for(k in 1:length(crosslocs))
  {
    myCop[[k]] <- tCopula(param=exp(-crosslocs[k]/theta.hat), dim = 2, dispstr = "un", df=nu)
  }
  
  
  for (i in 1:length(th)) {
    for (j in 1:length(th)) {
      
      uhat_i = mean(Y< th[i]) # E[I(Y(s)<th_i)] = P[Y(s)<th_i] = F_Y(th_i)
      uhat_j = mean(Y< th[j])
      
      b=vector()
      for(k in 1:length(crosslocs))
      {
        #estimated based on copula
        #myCop <- tCopula(param=exp(-crosslocs[k]/theta.hat), dim = 2, dispstr = "un", df=nu)
        b[k] = pCopula(c(uhat_i,uhat_j), myCop[[k]]) #P(U1<u, U2<u)
      }
      
      Cov = b - uhat_i*uhat_j
      
      totC12[,,i,j] = Cov
      
      
    }
  }
  
  
  myC12 <-matrix(NA,nrow=dim(crosslocs)[1]*length(th),ncol=dim(crosslocs)[2]*length(th))
  
  
  for(i in 1:length(th))
  {
    for(j in 1:length(th))
    {
      myC12[(dim(crosslocs)[1]*(i-1)+1):(dim(crosslocs)[1]*i),(dim(crosslocs)[2]*(j-1)+1):(dim(crosslocs)[2]*j)] = totC12[,,i,j]
    }
  }
  
  
  
  ################################### Indicator variable data ################################
  data_ind = list(); mu_ind= vector()
  data_all = NULL
  for (i in 1:length(th)) {
    data_ind[[i]] = as.numeric(Y< th[i])
    data_all = c(data_all,data_ind[[i]])
    
    # mu_ind[i] = mean(data_ind[[i]])
  }
  mu_ind = qseq
  mu_all = rep(mu_ind, each= length(Y))
  
  # inv_C = solve(myC)  outside the loop
  #ind_krig = myC12%*%inv_C%*%data_all
  ind_krig = mu_ind+ myC12%*%inv_C%*%(data_all - mu_all)
  
  F_ind = ind_krig
  
  
  ## Fix values at extremes, all values should lie between 0 and 1
  for (i in 1:length(F_ind)) {
    
    if(F_ind[i]<0){F_ind[i]=0}
    if(F_ind[i]>=1) {F_ind[i:length(F_ind)] = 1}
    
  }
  
  ############ Extrapolate to 0 and 1  ############
  xseq = th
  if(F_ind[length(xseq)]==1)
  {F_ind=F_ind} else { F_ind = c(F_ind,1)
  xseq = c(xseq, (max(data_nov$RR)+xseq[length(xseq)])/2) }
  
  if(F_ind[1]==0)
  {F_ind=F_ind} else {F_ind= c(0,F_ind)
  xseq = c( (min(data_nov$RR)+xseq[1])/2 ,xseq)}
  
  
  ## Order relation correction (Carr 1994)
  ##upward downward
  
  for (i in 1:(length(F_ind)-1)) {
    if(F_ind[i+1]< F_ind[i]){
      m = which(F_ind[(i+1):length(F_ind)]> F_ind[i])[1]
      
      for (j in (i+1):(i+m-1)) {
        F_ind[j] = F_ind[i] + ((F_ind[i+m]-F_ind[i])*(j-i))/m
      }
    }
  } 
  
  #F_ind_list[[l]] = cbind(xseq, F_ind)
  return(cbind(xseq, F_ind))
  
}





####### Smooth distribution using the global distribution and rescaling
#lines(smoothF_RR[[1]][[18]])

L=200; 

  smoothF_RR =   lapply (F_ind_list, function(a) smoothglobal(Fx=a, globalF = data_nov$RR) )
  median_pred_RR = as.numeric(lapply(smoothF_RR, function(a) {min(a[,1][which(a[,2]>=0.5)])}))
  quantile_pred_RR_l = matrix(0, nrow=L, ncol = length(median_pred_RR))
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_RR_l[l,] = as.numeric(lapply(smoothF_RR, function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  
  mean_pred_RR = apply(quantile_pred_RR_l, 2, mean)
  
  c=50
  F_c =   apply(quantile_pred_RR_l, 2, function(x) mean(x<=c))

  lowlim = uplim = length_gamma = matrix(0, nrow = 100, ncol = length(median_pred_RR)) 
  alpha = 0.1
  gamma = seq(0,alpha, length.out=100)
  for (i in 1:length(gamma)) {
    lowlim[i,] =  as.numeric(lapply(smoothF_RR, function(a) {min(a[,1][which(a[,2]>= gamma[i])])}))
    uplim[i,] =  as.numeric(lapply(smoothF_RR, function(a) {min(a[,1][which(a[,2]>=1-gamma[i])])}))
    length_gamma[i,] = uplim[i,]-lowlim[i,]
      }
  
  lim_index = apply(length_gamma,2, function(x) which.min(x) )
  
  lowlim_short = uplim_short = vector(length=length(median_pred_RR))
  for (j in 1:length(median_pred_RR)) {
    lowlim_short[j] =  lowlim[lim_index[j],j]
    uplim_short[j] = uplim[lim_index[j],j]
  }
  
  length_c = (uplim_short-lowlim_short)
  mean(length_ci_precip)
  summary(length_c)
  
  
  ## gaussian prediction for precipitation data
  #Estimate parameters by maximum likelihood:
  ml <- likfit(data=data_nov$RR,coords=coords,
               fix.nugget=F,cov.model="exponential", 
               ini = c(1500, 266639),nugget=18)
  
  
  pred<-krige.conv(data= data_nov$RR, coords= coords, locations=pred_grd_m,
                   krige=krige.control(cov.model="exponential", type.krige = "SK", beta = mean(data_nov$RR),
                                       cov.pars=c(ml$sigmasq ,ml$phi),
                                       nugget= ml$nugget))
  
  predicted_g = pred$predict; var_g = pred$krige.var
  ll_g = predicted_g+ qnorm(0.05,0,1)*sqrt(var_g); ul_g =predicted_g+ qnorm(0.95,0,1)*sqrt(var_g)
  length_g = ul_g-ll_g
  
library(viridis)
  spain = map_data("world", region = "spain")
  #plot(spain)
  pred_results = data.frame(loc=pred_grd, median= median_pred_RR, pexceed=1-F_c, 
                            length_ci_precip=length_c, length_g=length_g)
  
  ggplot(spain, aes(long, lat ))+geom_path(aes(group=group))+ coord_map()  + ggtitle(" Predicted Spain Precipitation (November 2019)")+
    geom_point(data = pred_results, aes(loc.x, loc.y, colour = median), size=2) +scale_color_viridis_c("Precipitation (mm)") +geom_path(aes(group=group))
  
  ggplot(spain, aes(long, lat ))+geom_path(aes(group=group))+ coord_map()  + ggtitle(" Probability of exceeding 100 mm")+
    geom_point(data = pred_results, aes(loc.x, loc.y, colour = pexceed), size=2) +scale_color_viridis_c("Probability") + geom_path(aes(group=group))

  p = ggmap(map_spain)  
  
  ## Prediction map
  p + ggtitle(" Predicted Spain Precipitation (November 2019)")+
    geom_point(data = pred_results, aes(loc.x, loc.y, colour = median),size=2, alpha=0.6) +
    scale_color_viridis("Precipitation (mm)")+ labs(x= "Longitude (degrees)", y="Latitude (degrees)")+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1.1)),
        legend.title = element_text(size=rel(1.3)),
        plot.title = element_text(size=rel(1.5)))
  
  ## Exceedance map
  p + ggtitle(" Probability of exceeding daily average 50 mm")+
   # geom_point(data = pred_results, aes(x=loc.x,y= loc.y, col = pexceed))+
   stat_contour_filled(data = pred_results, aes(x=loc.x,y= loc.y, z=pexceed), alpha=0.85)+
   labs(x= "Longitude (degrees)", y="Latitude (degrees)")+scale_fill_viridis_d(name="Probability")+ 
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1.1)),
        legend.title = element_text(size=rel(1.3)),
        plot.title = element_text(size=rel(1.5)))
        
  
  ## Uncertainity map: length of 90% ci
  p + ggtitle(" Length of 90% CI")+
    geom_point(data = pred_results, aes(loc.x, loc.y, colour = length_ci_precip-length_g),size=2, alpha=0.5) +
    scale_color_viridis("Length (mm)")+ labs(x= "Longitude (degrees)", y="Latitude (degrees)")
 
  p + ggtitle(" Length of 90% CI")+
    geom_point(data = pred_results, aes(loc.x, loc.y, colour = length_g),size=2, alpha=0.5) +
    scale_color_viridis("Length (mm)")+ labs(x= "Longitude (degrees)", y="Latitude (degrees)")
  
   
  

  
  

