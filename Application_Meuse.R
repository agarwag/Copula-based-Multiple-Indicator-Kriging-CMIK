library(sp)
library(gstat)
library(ggplot2)
data(meuse)
class(meuse)
head(meuse)
str(meuse)

meuse_sp = meuse
coordinates(meuse_sp) = ~x+y; 
proj4string(meuse_sp) <- CRS("+init=epsg:28992"); 
longlat <- spTransform(SpatialPoints(coordinates(meuse_sp),proj4string=meuse_sp@proj4string), 
                       CRS("+proj=longlat +datum=WGS84"))

meuse$Longitude = coordinates(longlat)[,1]; meuse$Latitude = coordinates(longlat)[,2]



## Spatial map of observed data
load("/home/agarwag/Desktop/KAUST/p4/map_meuse.RData")
library(ggmap)
library(viridis)
ggmap(map_meuse)+
  scale_x_continuous(limits = c(5.71,5.78))+
  scale_y_continuous(limits = c(50.95,50.995))+
  geom_point(data=meuse , aes(Longitude,Latitude, col = zinc), size=2) +
  scale_color_viridis("Zinc (ppm)")+
labs(x="Longitude (degrees)", y="Latitude (degrees)") +
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1.1)),
        legend.title = element_text(size=rel(1.3)))


#### Histogram of zinc concentrations
hist(meuse$zinc, xlab = "Zinc concentrations (ppm)", main="", cex.main=1.5, cex.lab=1.4)

ggplot() + geom_histogram(aes(x=meuse$zinc), breaks=hist(meuse$zinc, plot=F)$breaks, col="black")+
              xlab("Zinc concentrations (ppm)")+ ylab("Frequency")+
              theme(axis.text = element_text(size=rel(1.2)),
                    axis.title = element_text(size=rel(1.5)))
                    


hist(log(meuse$zinc))
####### Transforming Longitude Latitude to Universal Transverse Mercator (UTM) projections ####
LongLatToUTM<-function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## assigning a CRS (Coordinate Reference System)
  
  res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep=''))) ## Transforming the CRS
  return(data.frame(x= res$X, y= res$Y))
  #return(as.data.frame(res))
}

meuse_m <- LongLatToUTM(meuse$Longitude,meuse$Latitude,zone = 31)

## changing rdh coodinates to utm for meuse data analysis
meuse$x=meuse_m$x; meuse$y=meuse_m$y

library(geoR);library(fields);library(maps)

## split the data into 80% training and 20% testing dataset
#set.seed(124)
set.seed(126)
#set.seed(123)
n_iters=500
train_data_list = list(); test_data_list=list()

for (index in 1:n_iters) {
  
dt = sort(sample(nrow(meuse), nrow(meuse)*.8, replace = F))
  
train_data_list[[index]] = meuse[dt,]
test_data_list[[index]] = meuse[-dt,]
}

### OR

## split the data intok groups for k-fold crossvalidation
#set.seed(123)
#n_iters = k_cv = 10
#ss = sample(1:k_cv, nrow(meuse), replace=T)
#for (index in 1:k_cv) {
#  train_data_list[[index]] = meuse[ss!=index,]
#  test_data_list[[index]] = meuse[ss==index,]
#}


parallel_gp = list()
for (index in 1:n_iters) {
  
  train_data = train_data_list[[index]]
  test_data = test_data_list[[index]]
  
  train_meuse = train_data[,"zinc"]; train_loc = cbind(x = train_data$x, y= train_data$y)
  
  #### predict cdfs on the test data
  test_meuse =  test_data[,"zinc"]; test_loc = cbind(x= test_data$x, y= test_data$y)
  
  data_train_gp = data.frame(train_loc, zinc = train_meuse)
  
  data_test_gp = data.frame(test_loc, zinc = test_meuse)
  
  coordinates(data_train_gp) = ~x+y
  coordinates(data_test_gp) = ~x+y
  
  plot(variogram(zinc~1, data_train_gp))
  #Estimate parameters by maximum likelihood:
  ml <- likfit(data=data_train_gp$zinc,coords=train_loc,
               fix.nugget=F,cov.model="exponential", 
               ini = c(328896, 1193),nugget=13286)
  
  
  pred<-krige.conv(data= data_train_gp$zinc, coords= train_loc, locations=test_loc,
                   krige=krige.control(cov.model="exponential", type.krige = "SK", beta = mean(train_meuse),
                                       cov.pars=c(ml$sigmasq ,ml$phi),
                                       nugget= ml$nugget))
  
  parallel_gp[[index]] = list(tv = test_meuse, predicted = pred$predict, var = pred$krige.var)
  
}


#### rmse, mad, pit and crps: gaussian prediction
library(scoringRules)
crpsF_y_iter= vector()

pit_zinc_g =  vector(mode = "list", length = n_iters)
mad_g = rmse_g= crpsF_y_g= matrix(0, nrow = n_iters, ncol = length(test_data_list[[1]]$zinc))

crps_iter = vector(); crps_g = vector()
length_ci = vector();length_ci_g= vector(); alpha=0.1
for (i in 1:n_iters) {
  ### generate data from predictive distribution at all test locations
  tv = parallel_gp[[i]]$tv
  predicted = parallel_gp[[i]]$predicted
  var = parallel_gp[[i]]$var
  
  for (j  in 1:length(tv)) {
    pit_zinc_g[[i]][j] = pnorm(tv[j], mean= predicted[j], sd = sqrt(var[j]))
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
hist(unlist(pit_zinc_g), freq=F, xlab = "PIT", main=paste("Histogram of PIT "))
abline(h=1)    
mean(apply(rmse_g, 2, function(x) sqrt(mean(x))))
mean(apply(mad_g, 2, median))
mean(apply(crpsF_y_g,2,mean))
mean(apply(crpsF_y_g,2,median))

#plot( parallel_gp[[2]]$tv, parallel_gp[[2]]$predicted)



### log transformation for Gaussian process


parallel_gp_log = list()
for (index in 1:n_iters) {
  
  train_data = train_data_list[[index]]
  test_data = test_data_list[[index]]
  
  train_meuse = log(train_data[,"zinc"]); train_loc = cbind(x = train_data$x, y= train_data$y)
  
  #### predict cdfs on the test data
  test_meuse =  test_data[,"zinc"]; test_loc = cbind(x= test_data$x, y= test_data$y)
  
  data_train_log_gp = data.frame(train_loc, zinc = train_meuse)
  data_test_log_gp = data.frame(test_loc, zinc = test_meuse)
  
  coordinates(data_train_log_gp) = ~x+y
  coordinates(data_test_log_gp) = ~x+y
  
  #plot(variogram(zinc~1, data_train_log_gp))
  #Estimate parameters by maximum likelihood:
  ml <- likfit(data=data_train_log_gp$zinc,coords=train_loc,
               fix.nugget=F,cov.model="exponential", 
               ini = c(0.9, 500),nugget=0.04)
  
  
  pred<-krige.conv(data= data_train_log_gp$zinc, coords= train_loc, locations=test_loc,
        krige=krige.control(cov.model="exponential", type.krige = "SK", beta = mean(train_meuse),
                                       cov.pars=c(ml$sigmasq ,ml$phi),
                                       nugget= ml$nugget))
  
  parallel_gp_log[[index]] = list(tv = test_meuse, predicted = exp(pred$predict+0.5*pred$krige.var),
                                  var = exp(2*pred$predict +pred$krige.var)*(exp(pred$krige.var)-1))
 
}


#### rmse, mad, pit and crps: gaussian prediction
library(scoringRules)
crpsF_y_iter= vector()

pit_zinc_g_log =  vector(mode = "list", length = n_iters)
mad_g_log = rmse_g_log= crpsF_y_g_log= matrix(0, nrow = n_iters, ncol = length(test_data_list[[1]]$zinc))

length_ci = vector();length_ci_g= vector(); alpha=0.1
for (i in 1:n_iters) {
  ### generate data from predictive distribution at all test locations
  tv = parallel_gp_log[[i]]$tv
  predicted = parallel_gp_log[[i]]$predicted
  var = parallel_gp_log[[i]]$var
  
  for (j  in 1:length(tv)) {
    pit_zinc_g_log[[i]][j] = pnorm(tv[j], mean= predicted[j], sd = sqrt(var[j]))
    crps_iter[j] =   crps(y = tv[j], family = "normal",  mean= predicted[j], sd = sqrt(var[j]))
    #lower_ci = qnorm(alpha/2, mean= predicted[j], sd = sqrt(var[j]))
    #upper_ci = qnorm(1-alpha/2, mean= predicted[j], sd = sqrt(var[j]))
   # length_ci[j] = upper_ci - lower_ci
  }
  
  crpsF_y_g_log[i,] =  crps_iter
  rmse_g_log[i,] = (tv-predicted)^2 ## RMSE
  mad_g_log[i,] = abs(tv-predicted) ## MAD
  #length_ci_g[i] = mean(length_ci)
}
hist(unlist(pit_zinc_g_log), freq=F, xlab = "PIT", main=paste("Histogram of PIT "))
abline(h=1)    
mean(apply(rmse_g_log, 2, function(x) sqrt(mean(x))))
mean(apply(mad_g_log, 2, median))
mean(apply(crpsF_y_g_log,2,mean))

mean(apply(crpsF_y_g_log,2,median))

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





#meuse$TPrp = c(scale(meuse$TPrp))
qseq = c(seq(0.025,0.5,length.out = 5),0.55, 0.6,0.65, seq(0.7, 0.99, length.out = 7)) ## saved environment 

#qseq = c(seq(0.025,0.5,length.out = 5),0.55, 0.6,0.65, seq(0.7, 0.975, length.out = 9)) ## regularly used

#qseq = c(0.05,0.1625, 0.275, 0.3875,0.46, 0.5,0.55, 0.6,0.625, 0.65, seq(0.7, 0.975, length.out = 10))

#qseq = c(seq(0.1,0.8, by=0.1), seq(0.85,0.99,length.out = 4)) ## good pit
#seq = c(seq(0.025,0.45,length.out = 5),0.5, 0.55, 0.6,0.63, 0.675, seq(0.7, 0.99, length.out = 10))

th = quantile(meuse$zinc, qseq )

#summary(meuse$zinc)
## choosing threshold
#th = c(seq(0.2,3, length.out = 15))
plot(ecdf(meuse$zinc), main="Global ECDF of zinc conc.")
abline(v=th, lty=2)


start.time = Sys.time()
library(doParallel)
registerDoParallel(cores=50)
parallel = list(); tv_zinc = list()
for (index in 1:n_iters) {
  
  train_data = train_data_list[[index]]
  test_data = test_data_list[[index]]
  
  train_meuse = train_data[,"zinc"]
  train_loc = cbind(x = train_data$x, y= train_data$y)
  

  th = quantile(train_meuse, qseq )
  #plot(ecdf(train_meuse), main="Global ECDF of zinc conc.")
  #abline(v=th, lty=2)
  
  Y = train_meuse
  coords = train_loc 
  
 
  theta.hat <- optim(par = 400,logL.Gcopula.Exp,Y=Y,coord=coords, method="L-BFGS-B",
                     lower=0.01,upper=2000)$par
  

  
  dist_mat =  as.matrix(dist(coords))
  h = dist_mat[lower.tri(dist_mat)]
  
  totC<-array(NA,dim = c(dim(dist_mat)[1],dim(dist_mat)[2],length(th),length(th)))
  
  ## i from 
  e = rep(1:length(th), times=1:length(th))
  f = vector(); for(j in 1:length(th)) {f = c(f, 1:j)}
  ij = cbind(e,f)
  
  myCop = vector(mode = "list", length = length(h))
  for (k in 1:length(h)) {
    myCop[[k]] =  normalCopula(param=exp(-h[k]/theta.hat), dim = 2, dispstr = "un")
  }
  
  b=vector()
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
  test_meuse =  test_data[,"zinc"]
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
         # myCop <- normalCopula(param=exp(-crosslocs[k]/theta.hat), dim = 2, dispstr = "un")
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
    xseq = c(xseq, (max(meuse$zinc)+xseq[length(xseq)])/2) }
    
    if(F_ind[1]==0)
    {F_ind=F_ind} else {F_ind= c(0,F_ind)
    xseq = c( (min(meuse$zinc)+xseq[1])/2 ,xseq)}
    
    
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
  tv_zinc[[index]] = test_meuse
  
}

end.time = Sys.time()
end.time  - start.time




####### Smooth distribution using the global distribution and rescaling
#lines(smoothF_zinc[[1]][[18]])
smoothglobal = function(Fx, globalF){
  
  ecdf_zinc = ecdf(globalF)
  ecdf_x = environment(ecdf_zinc)$x; ecdf_y = environment(ecdf_zinc)$y 
  
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

smoothtry = smoothglobal(Fx = parallel[[1]][[2]], globalF = meuse$zinc )
plot(parallel[[1]][[2]], ylab="F")
lines(smoothtry,  col="grey" )

L=200; 
smoothF_zinc = list()
pit_zinc =  vector(mode = "list", length = n_iters) ; 
rmse_zinc=mad_zinc =crpsF_y= matrix(0, nrow = n_iters, ncol = length(test_data_list[[1]]$zinc))
length_ci_gcop = vector(); alpha=0.1


x_seq = seq(min(meuse$zinc), max(meuse$zinc), length.out = 10000) #zinc
x_mid = x_seq[-length(x_seq)] + diff(x_seq)/2
delta_x = diff(x_seq)[1]  

#parallel_sim[[1]][[i]]$Fsim[[1]]
for (index in 1:n_iters) {
  smoothF_zinc[[index]] =   lapply (parallel[[index]], function(a) smoothglobal(Fx=a, globalF = meuse$zinc) )
  median_pred_zinc = as.numeric(lapply(smoothF_zinc[[index]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
  quantile_pred_zinc_l = matrix(0, nrow=L, ncol = length(median_pred_zinc))
 
  lower_ci = as.numeric(lapply(smoothF_zinc[[index]], function(a) {min(a[,1][which(a[,2]>=alpha/2)])}))
  upper_ci = as.numeric(lapply(smoothF_zinc[[index]], function(a) {min(a[,1][which(a[,2]>=1-alpha/2)])}))
  length_ci_gcop[index] = mean(upper_ci-lower_ci)
  
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_zinc_l[l,] = as.numeric(lapply(smoothF_zinc[[index]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  
  #mean_pred_zinc = apply(quantile_pred_zinc_l, 2, mean)
  tv = tv_zinc[[index]]
  rmse_zinc[index,] = (tv-median_pred_zinc)^2 ## RMSE
 # rmse_RR_mean[index] = sqrt(mean((tv-mean_pred_RR)^2)) ## RMSE
  mad_zinc[index,] = abs(tv-median_pred_zinc)
  
  crpsF_y_iter= vector()
  for (j  in 1:length(tv)) {
    
    pit_zinc[[index]][j] =   mean(quantile_pred_zinc_l[,j]<=tv[j])
    f_x = function(x,y) { {mean(quantile_pred_zinc_l[,j]<=x) - ifelse( y <= x,1,0)}^2 }
    fx_mid = sapply(x_mid, function(x) f_x(x, tv[j]))
    crpsF_y_iter[j] = sum(delta_x*fx_mid, na.rm = T)
  }
  crpsF_y[index,] =  crpsF_y_iter
}

hist(unlist(pit_zinc), freq=F, main=paste("Histogram of PIT"), breaks = 10)
abline(h=1)
mean(apply(crpsF_y,2,mean))
mean(apply(rmse_zinc, 2, function(x) sqrt(mean(x))))
mean(apply(mad_zinc, 2, median))

sd(crpsF_y)
sd(crpsF_y_g)



#################################################################################################
###############################################################################################
###################################### Variogram Approach  #############################################
#qseq = c(seq(0.05,0.5,length.out = 5),0.55, 0.6,0.65, seq(0.7, 0.975, length.out = 7))
#qseq = c(seq(0.05,0.95,length.out = 15))

#qseq = seq(0.1,0.9, length.out = 10)
th = quantile(meuse$zinc, qseq )

summary(meuse$zinc)
## choosing threshold
#th = c(seq(150,1000, length.out = 11), seq(1100, 1600, length.out = 4))
plot(ecdf(meuse$zinc), main="Global ECDF of zinc concentrations")
abline(v=th, lty=2)
library(gstat)
set.seed(123)
parallel_var= list(); tv_zinc_var = list()
for (index in 1:n_iters) {
  
  
  train_data = train_data_list[[index]]
  test_data = test_data_list[[index]]
  
  train_meuse = train_data[,"zinc"]
  train_loc = cbind(x = train_data$x, y= train_data$y)
  
  th = quantile(train_meuse, qseq )
  
  data_train_var = data.frame(train_loc,zinc = train_meuse)
  
  #### predict cdfs on the test data
  test_data_zinc = test_data[,"zinc"]
  test_loc = cbind(x = test_data$x, y= test_data$y)
  
  data_test_var = data.frame(test_loc, zinc = test_data_zinc)
  
  coordinates(data_train_var) = ~x+y
  coordinates(data_test_var) = ~x+y
  
  ### indicator cokriging

  quartz = th
  # These are the actual values of the quantiles
  
  zinc.i <- gstat(id = "zinc1", formula = I(zinc < quartz[1]) ~ 1, data = data_train_var, 
                 nmax = 7, beta = qseq[1], set = list(order = 4, zero = 1e-05))
  zinc.i <- gstat(zinc.i, "zinc2", formula = I(zinc < quartz[2]) ~ 1, data = data_train_var, 
                 nmax = 7, beta = qseq[2])
  zinc.i <- gstat(zinc.i, "zinc3", formula = I(zinc < quartz[3]) ~ 1, data = data_train_var, 
                 nmax = 7, beta = qseq[3])
  zinc.i <- gstat(zinc.i, "zinc4", formula = I(zinc < quartz[4]) ~ 1, data = data_train_var, 
                 nmax = 7, beta = qseq[4])
  zinc.i <- gstat(zinc.i, "zinc5", formula = I(zinc < quartz[5]) ~ 1, data = data_train_var, 
                 nmax = 7, beta = qseq[5])
  zinc.i <- gstat(zinc.i, "zinc6", formula = I(zinc < quartz[6]) ~ 1, data = data_train_var, 
                 nmax = 7, beta = qseq[6])
  zinc.i <- gstat(zinc.i, "zinc7", formula = I(zinc < quartz[7]) ~ 1, data = data_train_var, 
                 nmax = 7, beta = qseq[7])
  zinc.i <- gstat(zinc.i, "zinc8", formula = I(zinc < quartz[8]) ~ 1, data = data_train_var, 
                 nmax = 7, beta = qseq[8])
  zinc.i <- gstat(zinc.i, "zinc9", formula = I(zinc < quartz[9]) ~ 1, data = data_train_var, 
                 nmax = 7, beta = qseq[9])
  zinc.i <- gstat(zinc.i, "zinc10", formula = I(zinc < quartz[10]) ~ 1, data = data_train_var, 
                 nmax = 7, beta = qseq[10])
  zinc.i <- gstat(zinc.i, "zinc11", formula = I(zinc < quartz[11]) ~ 1, data = data_train_var, 
                 nmax = 7, beta = qseq[11])
  zinc.i <- gstat(zinc.i, "zinc12", formula = I(zinc < quartz[12]) ~ 1, data = data_train_var, 
                 nmax = 7, beta = qseq[12])
  zinc.i <- gstat(zinc.i, "zinc13", formula = I(zinc < quartz[13]) ~ 1, data = data_train_var, 
                 nmax = 7, beta = qseq[13])
  zinc.i <- gstat(zinc.i, "zinc14", formula = I(zinc < quartz[14]) ~ 1, data = data_train_var, 
                 nmax = 7, beta = qseq[14])
  zinc.i <- gstat(zinc.i, "zinc15", formula = I(zinc < quartz[15]) ~ 1, data = data_train_var, 
                 nmax = 7, beta = qseq[15])
#  zinc.i <- gstat(zinc.i, "zinc16", formula = I(zinc < quartz[16]) ~ 1, data = data_train_var, 
#                  nmax = 7, beta = qseq[16])
#  zinc.i <- gstat(zinc.i, "zinc17", formula = I(zinc < quartz[17]) ~ 1, data = data_train_var, 
#                  nmax = 7, beta = qseq[17])
#  zinc.i <- gstat(zinc.i, "zinc18", formula = I(zinc < quartz[18]) ~ 1, data = data_train_var, 
#                  nmax = 7, beta = qseq[18])
#  zinc.i <- gstat(zinc.i, "zinc19", formula = I(zinc < quartz[19]) ~ 1, data = data_train_var, 
#                  nmax = 7, beta = qseq[19])
#  zinc.i <- gstat(zinc.i, "zinc20", formula = I(zinc < quartz[20]) ~ 1, data = data_train_var, 
#                  nmax = 7, beta = qseq[20])
  
  
  zinc.i <- gstat(zinc.i, model = vgm(0.1, "Exp", 500, 0.1), fill.all = T)
  
  
  zinc.vg <- variogram(zinc.i)
  #plot(zinc.vg)
  zinc.quartfit = fit.lmc(zinc.vg, zinc.i)
  
  Fvar_list= list()
  # now do the co-kriging
  for (l in 1:dim(data_test_var@coords)[1]) {
    
    cokrig.zincquart <- predict(zinc.quartfit, data_test_var[l,data_test_var$coordinates])
    
    Fvar = c(cokrig.zincquart@data$zinc1.pred,cokrig.zincquart@data$zinc2.pred,cokrig.zincquart@data$zinc3.pred,cokrig.zincquart@data$zinc4.pred,cokrig.zincquart@data$zinc5.pred,cokrig.zincquart@data$zinc6.pred,cokrig.zincquart@data$zinc7.pred,cokrig.zincquart@data$zinc8.pred,
             cokrig.zincquart@data$zinc9.pred,cokrig.zincquart@data$zinc10.pred,cokrig.zincquart@data$zinc11.pred,cokrig.zincquart@data$zinc12.pred,cokrig.zincquart@data$zinc13.pred,cokrig.zincquart@data$zinc14.pred,cokrig.zincquart@data$zinc15.pred) #,
            # cokrig.zincquart@data$zinc16.pred,cokrig.zincquart@data$zinc17.pred,cokrig.zincquart@data$zinc18.pred,cokrig.zincquart@data$zinc19.pred,cokrig.zincquart@data$zinc20.pred)
    
    
    ############ Extrapolate to 0 and 1  ############
    xseq= quartz
    if(Fvar[length(quartz)]==1)
    {Fvar=Fvar} else {Fvar = c(Fvar,1)
    xseq = c(xseq,(max(data_train_var$zinc)+xseq[length(xseq)])/2)}
    
    if(Fvar[1]==0)
    {Fvar=Fvar} else {Fvar= c(0,Fvar)
    xseq = c( (min(data_train_var$zinc)+xseq[1])/2 ,xseq)}
    
    Fvar_list[[l]] = cbind(xseq, Fvar)
    
  }

  #parallel_var[[index]] = list(test_data_zinc, median_pred_var)
  parallel_var[[index]] = Fvar_list
  tv_zinc_var[[index]] = test_data_zinc
  
}

####### Smooth distribution using the global distribution and rescaling
crpsF_y_var =vector();  crpsF_y_iter= vector(); smoothF_zinc = list()
pit_zinc_var =  vector(mode = "list", length = n_iters) ; 
rmse_zinc_var= mad_zinc_var= crpsF_y_var=matrix(0, nrow = n_iters, ncol = length(test_data_list[[1]]$zinc))
L=200; 
length_ci_var= vector()


x_seq = seq(min(meuse$zinc), max(meuse$zinc), length.out = 10000) #zinc
x_mid = x_seq[-length(x_seq)] + diff(x_seq)/2
delta_x = diff(x_seq)[1]  

#parallel_var_sim[[1]][[i]]$Fsim[[1]]
for (index in 1:n_iters) {
  smoothF_zinc[[index]] =   lapply (parallel_var[[index]], function(a) smoothglobal(Fx=a, globalF = meuse$zinc) )
  median_pred_zinc = as.numeric(lapply(smoothF_zinc[[index]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
  
  lower_ci = as.numeric(lapply(smoothF_zinc[[index]], function(a) {min(a[,1][which(a[,2]>=alpha/2)])}))
  upper_ci = as.numeric(lapply(smoothF_zinc[[index]], function(a) {min(a[,1][which(a[,2]>=1-alpha/2)])}))
  length_ci_var[index] = mean(upper_ci-lower_ci)
  
  quantile_pred_zinc_l = matrix(0, nrow=L, ncol = length(median_pred_zinc))
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_zinc_l[l,] = as.numeric(lapply(smoothF_zinc[[index]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  
  mean_pred_zinc = apply(quantile_pred_zinc_l, 2, mean)
  tv = tv_zinc_var[[index]]
  rmse_zinc_var[index,] = (tv-median_pred_zinc)^2 ## RMSE
  mad_zinc_var[index,] = abs(tv-median_pred_zinc) ## MAD
  
  for (j  in 1:length(tv)) {
    
    pit_zinc_var[[index]][j] =   mean(quantile_pred_zinc_l[,j]<=tv[j])
    f_x = function(x,y) { {mean(quantile_pred_zinc_l[,j]<=x) - ifelse( y <= x,1,0)}^2 }
    fx_mid = sapply(x_mid, function(x) f_x(x, tv[j]))
    crpsF_y_iter[j] = sum(delta_x*fx_mid, na.rm = T)
  }
  crpsF_y_var[index,] =  crpsF_y_iter
}

hist(unlist(pit_zinc_var), freq=F, main=paste("Histogram of PIT"))
abline(h=1)
mean(apply(crpsF_y_var,2,mean))
mean(apply(rmse_zinc_var,2, mean))
mean(apply(mad_RR_var,2, median))




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
  R.m <- exp(-dm/par)
  params = R.m[lower.tri(R.m)]
  
  #myCop <- normalCopula(param=params, dim = n, dispstr = "un")  ## dim of coupla is equal to the number of spatial locations
  u = ecdf(Y)(Y)  ### extract the empirical distibution functions value
  
  #truncate u to avoid u equals exactly zero or 1
  u = pmax(0.001,pmin(u, 0.999))
  
  #use true marginal
  #u = pnorm(Y)
  nu= 4
  logl <- dmvt(qt(u, df = nu), df=nu, sigma=R.m,log = TRUE) - sum(dt(qt(u, df=nu), df=nu, log=TRUE)) # log-density for X2(s) (denominator)
  
  #logl = log(dCopula(u, myCop))
  return(-logl)
}




#meuse$TPrp = c(scale(meuse$TPrp))
#qseq = c(seq(0.05,0.5,length.out = 5),0.55, 0.6,0.65, seq(0.7, 0.975, length.out = 7))

#qseq = c(0.05,0.1625, 0.275, 0.3875,0.46, 0.5,0.55, 0.6,0.625, 0.65, seq(0.7, 0.975, length.out = 10))

#qseq = seq(0.1,0.9, length.out = 15)
th = quantile(meuse$zinc, qseq )


## choosing threshold
plot(ecdf(meuse$zinc), main="Global ECDF of zinc conc.")
abline(v=th, lty=2)


start.time = Sys.time()
library(doParallel)
registerDoParallel(cores=50)
parallel_t = list(); tv_zinc = list()
for (index in 1:n_iters) {
  
  train_data = train_data_list[[index]]
  test_data = test_data_list[[index]]
  
  train_meuse = train_data[,"zinc"]
  train_loc = cbind(x = train_data$x, y= train_data$y)
  
  th = quantile(train_meuse, qseq )
  
  Y = train_meuse
  coords = train_loc ### converting distances from m to km
  
  
  theta.hat_vec <- optim(par = c(400),logL.tcopula.Exp,Y=Y,coord=coords, method="L-BFGS-B",
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
     # myCop <- tCopula(param=exp(-h[k]/theta.hat), dim = 2, dispstr = "un", df=nu)
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
  test_meuse =  test_data[,"zinc"]
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
      
      if(F_ind[i]<0){F_ind[i]=0}
      if(F_ind[i]>=1) {F_ind[i:length(F_ind)] = 1}
      
    }
    
    ############ Extrapolate to 0 and 1  ############
    xseq = th
    if(F_ind[length(xseq)]==1)
    {F_ind=F_ind} else { F_ind = c(F_ind,1)
    xseq = c(xseq, (max(meuse$zinc)+xseq[length(xseq)])/2) }
    
    if(F_ind[1]==0)
    {F_ind=F_ind} else {F_ind= c(0,F_ind)
    xseq = c( (min(meuse$zinc)+xseq[1])/2 ,xseq)}
    
    
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
  tv_zinc[[index]] = test_meuse
  
}

end.time = Sys.time()
end.time  - start.time



####### Smooth distribution using the global distribution and rescaling

L=200; 
crpsF_y_t =vector();  smoothF_zinc = list()
pit_zinc_t =  vector(mode = "list", length = n_iters) ; 
rmse_zinc_t = mad_zinc_t =crpsF_y_t= matrix(0, nrow = n_iters, ncol = length(test_data_list[[1]]$zinc))
length_ci_tcop = vector()

x_seq = seq(min(meuse$zinc), max(meuse$zinc), length.out = 10000) #zinc
x_mid = x_seq[-length(x_seq)] + diff(x_seq)/2
delta_x = diff(x_seq)[1]  

for (index in 1:n_iters) {
  smoothF_zinc[[index]] =   lapply (parallel_t[[index]], function(a) smoothglobal(Fx=a, globalF = meuse$zinc) )
  median_pred_zinc = as.numeric(lapply(smoothF_zinc[[index]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
  
  lower_ci = as.numeric(lapply(smoothF_zinc[[index]], function(a) {min(a[,1][which(a[,2]>=alpha/2)])}))
  upper_ci = as.numeric(lapply(smoothF_zinc[[index]], function(a) {min(a[,1][which(a[,2]>=1-alpha/2)])}))
  length_ci_tcop[index] = mean(upper_ci-lower_ci)
  
  quantile_pred_zinc_l = matrix(0, nrow=L, ncol = length(median_pred_zinc))
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_zinc_l[l,] = as.numeric(lapply(smoothF_zinc[[index]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  
  mean_pred_zinc = apply(quantile_pred_zinc_l, 2, mean)
  tv = tv_zinc[[index]]
  rmse_zinc_t[index,] = (tv-median_pred_zinc)^2 ## RMSE
  mad_zinc_t[index,] = abs(tv-median_pred_zinc) ## MAD
  
  crpsF_y_iter= vector()
  for (j  in 1:length(tv)) {
    
    pit_zinc_t[[index]][j] =   mean(quantile_pred_zinc_l[,j]<=tv[j])
    f_x = function(x,y) { {mean(quantile_pred_zinc_l[,j]<=x) - ifelse( y <= x,1,0)}^2 }
    fx_mid = sapply(x_mid, function(x) f_x(x, tv[j]))
    crpsF_y_iter[j] = sum(delta_x*fx_mid, na.rm = T)
  }
  crpsF_y_t[index,] =  crpsF_y_iter
}

hist(unlist(pit_zinc_t), freq=F, main=paste("Histogram of PIT"), breaks = 10)
abline(h=1)

mean(apply(crpsF_y_t,2,mean))
mean(apply(rmse_zinc_t,2, mean)); 
mean(apply(mad_zinc_t, 2, median))

sd(crpsF_y_t)
sd(crpsF_y_g)







##### average length of ci and coverage probability #################
alpha = 0.1
### gaussian kriging
length_ci_g_zinc = vector(); cp_g_zinc = vector()

for (i in 1:n_iters) {
  tv = parallel_gp[[i]]$tv
  predicted = parallel_gp[[i]]$predicted
  var = parallel_gp[[i]]$var
  
  ll = qnorm(alpha/2, mean = predicted, sd = sqrt(var))
  ul = qnorm(1-alpha/2, mean = predicted, sd =sqrt(var))
  length_ci_g_zinc[i] = mean(ul-ll)
  cp_g_zinc[i] = sum(tv > ll & tv < ul)/length(tv)
}
mean(cp_g_zinc)
mean(length_ci_g_zinc)


length_ci_logg_zinc = vector(); cp_logg_zinc = vector()
for (i in 1:n_iters) {
  tv = parallel_gp_log[[i]]$tv
  predicted = parallel_gp_log[[i]]$predicted
  var = parallel_gp_log[[i]]$var
  
  ll = qnorm(alpha/2, mean = predicted, sd = sqrt(var))
  ul = qnorm(1-alpha/2, mean = predicted, sd =sqrt(var))
  length_ci_logg_zinc[i] = mean(ul-ll)
  cp_logg_zinc[i] = sum(tv > ll & tv < ul)/length(tv)
}
mean(cp_logg_zinc)
mean(length_ci_logg_zinc)

## variogram
############## average length ci and coverage proabability
length_ci_var_zinc = vector(); cp_var_zinc= vector()
for (i in 1:n_iters) {
  smoothF_zinc[[i]] =   lapply (parallel_var[[index]], function(a) smoothglobal(Fx=a, globalF = meuse$zinc) )
  
  ll = as.numeric(lapply(smoothF_zinc[[i]], function(a) {min(a[,1][which(a[,2]>= alpha/2)])}))
  ul = as.numeric(lapply(smoothF_zinc[[i]], function(a) {min(a[,1][which(a[,2]>= 1- alpha/2)])}))
  
  tv = tv_zinc_var[[i]]
  length_ci_var_zinc[i] = mean(ul-ll)
  cp_var_zinc[i] = sum(tv > ll & tv < ul)/length(tv)
  
}
mean(length_ci_var_zinc)
mean(cp_var_zinc)

## Gaussian copula
length_ci_zinc = vector(); cp_zinc= vector()
for (i in 1:n_iters) {
  smoothF_zinc[[index]] =   lapply (parallel[[index]], function(a) smoothglobal(Fx=a, globalF = meuse$zinc) )
  
  ll = as.numeric(lapply(smoothF_zinc[[i]], function(a) {min(a[,1][which(a[,2]>= alpha/2)])}))
  ul = as.numeric(lapply(smoothF_zinc[[i]], function(a) {min(a[,1][which(a[,2]>= 1- alpha/2)])}))
  
  tv = tv_zinc[[i]]
  length_ci_zinc[i] = mean(ul-ll)
  cp_zinc[i] = sum(tv > ll & tv < ul)/length(tv)
  
}
mean(length_ci_zinc)
mean(cp_zinc)
### t copula
length_ci_t_zinc = vector(); cp_t_zinc= vector()
for (i in 1:n_iters) {
  smoothF_zinc[[i]] =   lapply (parallel_t[[i]], function(a) smoothglobal(Fx=a, globalF = meuse$zinc) )
  
  ll = as.numeric(lapply(smoothF_zinc[[i]], function(a) {min(a[,1][which(a[,2]>= alpha/2)])}))
  ul = as.numeric(lapply(smoothF_zinc[[i]], function(a) {min(a[,1][which(a[,2]>= 1- alpha/2)])}))
  
  tv = tv_zinc[[i]]
  length_ci_t_zinc[i] = mean(ul-ll)
  cp_t_zinc[i] = sum(tv > ll & tv < ul)/length(tv)
  
}
mean(length_ci_t_zinc)
mean(cp_t_zinc)


## figures and tables for manuscript

hist(unlist(pit_zinc_g), freq=F, ylim=c(0,2.2), main=paste("Gaussian Kriging"), 
     xlab="Probability Integral Transform", breaks=10, cex.main=1.5, cex.lab=1.4)
abline(h=1)
hist(unlist(pit_zinc_g_log), freq=F, ylim=c(0,2.2), main=paste("Gaussian Kriging (log transformation)"), 
     xlab="Probability Integral Transform", breaks=10, cex.main=1.5, cex.lab=1.4)
abline(h=1)
hist(unlist(pit_zinc_var), freq=F, ylim=c(0,2.2), main=paste("Variogram"), 
     xlab="Probability Integral Transform", breaks=10, cex.main=1.5, cex.lab=1.4)
abline(h=1)
hist(unlist(pit_zinc), freq=F,ylim=c(0,2.2), main=paste("Gaussian copula"), 
     xlab="Probability Integral Transform", breaks=10, cex.main=1.5, cex.lab=1.4)
abline(h=1)
hist(unlist(pit_zinc_t), freq=F, ylim=c(0,2.2), main=paste("t-copula"), 
     xlab="Probability Integral Transform", breaks=10, cex.main=1.5, cex.lab=1.4)
abline(h=1)





mean(apply(mad_g, 2, median)); mean(apply(mad_g_log, 2, median)); mean(apply(mad_zinc_var, 2, median)) ;mean(apply(mad_zinc, 2, median)) ; mean(apply(mad_zinc_t, 2, median)) 

mean(apply(crpsF_y_g,2,mean));mean(apply(crpsF_y_g_log,2,mean)); mean(apply(crpsF_y_var,2,mean)); mean(apply(crpsF_y,2,mean)); mean(apply(crpsF_y_t,2,mean))

mean(apply(crpsF_y_g,2,median)); mean(apply(crpsF_y_g_log,2,median)); mean(apply(crpsF_y_var,2,median)); mean(apply(crpsF_y,2,median)); mean(apply(crpsF_y_t,2,median))



##standard error
1.2533*rbind(mean(apply(mad_g, 2, sd)/sqrt(n_iters)),mean(apply(mad_g_log, 2, median)/sqrt(n_iters)),
             mean(apply(mad_zinc_var, 2, median)/sqrt(n_iters)),mean(apply(mad_zinc, 2, median)/sqrt(n_iters)),mean(apply(mad_zinc_t, 2, median)/sqrt(n_iters)))
rbind(mean(apply(crpsF_y_g, 2, sd)/sqrt(n_iters)), mean(apply(crpsF_y_g_log, 2, median)/sqrt(n_iters)),
      mean(apply(crpsF_y_var, 2, sd)/sqrt(n_iters)),mean(apply(crpsF_y, 2, median)/sqrt(n_iters)),mean(apply(crpsF_y_t, 2, median)/sqrt(n_iters)))
1.2533*rbind(mean(apply(crpsF_y_g, 2, sd)/sqrt(n_iters)), mean(apply(crpsF_y_g_log, 2, median)/sqrt(n_iters)),
             mean(apply(crpsF_y_var, 2, sd)/sqrt(n_iters)),mean(apply(crpsF_y, 2, median)/sqrt(n_iters)),mean(apply(crpsF_y_t, 2, median)/sqrt(n_iters)))




## average length and ci  

## average length ci and coverage ###
#length
mean(length_ci_g_zinc);mean(length_ci_logg_zinc);mean(length_ci_var_zinc);mean(length_ci_zinc);mean(length_ci_t_zinc)

# coverage
mean(cp_g_zinc);mean(cp_logg_zinc);mean(cp_var_zinc);mean(cp_zinc);mean(cp_t_zinc)













################# Application using t-Copula kriging method ############
data("meuse.area")
meuse_area_sp = as.data.frame(meuse.area)
coordinates(meuse_area_sp) = ~x+y;
proj4string(meuse_area_sp) <- CRS("+init=epsg:28992")
longlat.grid <- spTransform(SpatialPoints(coordinates(meuse_area_sp),proj4string=meuse_area_sp@proj4string), 
                            CRS("+proj=longlat +datum=WGS84"))

meuse.area_lonlat = cbind(Longitude = coordinates(longlat.grid)[,1],
                          Latitude = coordinates(longlat.grid)[,2])

x.range <- range(meuse.area_lonlat[,1])
y.range <- range(meuse.area_lonlat[,2])
grd <- as.matrix(expand.grid(x=seq(from=x.range[1], to=x.range[2], length.out = 100),
                             y=seq(from=y.range[1], to=y.range[2], length.out = 100) ))

plot(meuse.area_lonlat)
library(mgcv)
pred_grd_meuse = grd[in.out(meuse.area_lonlat,grd),] 
pred_grd_meuse_m <- LongLatToUTM(pred_grd_meuse[,1],pred_grd_meuse[,2],zone = 31)


library(copula)
library(fields)
library(sp)
library(mvtnorm)


###Function: Pseudo neg logLikelihood of t-copula
logL.tcopula.Exp <- function(par,Y, coord){
  
  
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
  nu= 4
  logl <- dmvt(qt(u, df = nu), df=nu, sigma=R.m,log = TRUE) - sum(dt(qt(u, df=nu), df=nu, log=TRUE)) # log-density for X2(s) (denominator)
  
  #logl = log(dCopula(u, myCop))
  return(-logl)
}




th = quantile(meuse$zinc, qseq )


## choosing threshold
plot(ecdf(meuse$zinc), main="Global ECDF of zinc conc.")
abline(v=th, lty=2)


start.time = Sys.time()
library(doParallel)
registerDoParallel(cores=50)
parallel = list(); tv_zinc = list()

  
  
  Y = meuse$zinc
  coords = cbind(meuse$x, meuse$y) ### converting distances from m to km
  
  
  theta.hat_vec <- optim(par = c(400),logL.tcopula.Exp,Y=Y,coord=coords, method="L-BFGS-B",
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
      # myCop <- tCopula(param=exp(-h[k]/theta.hat), dim = 2, dispstr = "un", df=nu)
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
  
  
  
  
  #### predict cdfs on the meuse grid 
  meuse_grid = data.frame(x=pred_grd_meuse_m$x, y=pred_grd_meuse_m$y)
  
  F_ind_list = list()
  inv_C = solve(myC)
  
  #dim(grd)[1]
  F_ind_list = foreach(l=1:dim(meuse_grid)[1]) %dopar% {
    
    newloc = data.frame(x= meuse_grid[l,1] ,y= meuse_grid[l,2]) 
    crosslocs<-rdist(newloc,coords)
    
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
    data_ind = list()
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
    xseq = c(xseq, (max(meuse$zinc)+xseq[length(xseq)])/2) }
    
    if(F_ind[1]==0)
    {F_ind=F_ind} else {F_ind= c(0,F_ind)
    xseq = c( (min(meuse$zinc)+xseq[1])/2 ,xseq)}
    
    
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
  

  
end.time = Sys.time()
end.time  - start.time



####### Smooth distribution using the global distribution and rescaling

  smoothglobal = function(Fx, globalF){
  
  ecdf_zinc = ecdf(globalF)
  ecdf_x = environment(ecdf_zinc)$x; ecdf_y = environment(ecdf_zinc)$y 
  
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

smoothtry = smoothglobal(Fx = F_ind_list[[20]], globalF = meuse$zinc )
plot(F_ind_list[[20]], ylab="F")
lines(smoothtry,  col="grey" )

L=200; 



  smoothF_zinc =   lapply (F_ind_list, function(a) smoothglobal(Fx=a, globalF = meuse$zinc) )
  median_pred_zinc = as.numeric(lapply(F_ind_list, function(a) {min(a[,1][which(a[,2]>=0.5)])}))
  quantile_pred_zinc_l = matrix(0, nrow=L, ncol = length(median_pred_zinc))
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_zinc_l[l,] = as.numeric(lapply(smoothF_zinc, function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  
  mean_pred_zinc = apply(quantile_pred_zinc_l, 2, mean)
  
  
  c=500
  F_c =   apply(quantile_pred_zinc_l, 2, function(x) mean(x<=c))
 
#The evident structure here is that zinc concentration is larger close to the river Meuse banks

## c=500 ppm
  
  pred_grd_meuse = as.data.frame(pred_grd_meuse); colnames(pred_grd_meuse) = c("Longitude", "Latitude")
  pred_grd_meuse$median_zinc = median_pred_zinc
  pred_grd_meuse$p_exceed = 1-F_c
  
  ### Predicted values
  ggplot(pred_grd_meuse) + geom_point(aes(Longitude,Latitude, col = median_zinc))+coord_equal()+
    scale_color_viridis_c("Zinc (ppm)")
  
  # Probability of exceedance
  ggplot(pred_grd_meuse) + geom_point(aes(Longitude,Latitude, col = p_exceed))+coord_equal() +
    scale_color_viridis_c("Probability")+ ggtitle("Probability of exceedance 1000 ppm") 
    
  
  
  p_meuse = ggmap(map_meuse)+
    scale_x_continuous(limits = c(5.71,5.78))+
    scale_y_continuous(limits = c(50.95,50.995))
  
  ## Prediction map
  p_meuse + ggtitle(" Predicted zinc concentrations (ppm)")+
    geom_point(data = pred_grd_meuse, aes(Longitude, Latitude, colour = median_zinc),size=1, alpha=0.8) +
    scale_color_viridis("Zinc (ppm)")+ labs(x= "Longitude (degrees)", y="Latitude (degrees)")+
    theme(axis.text = element_text(size=rel(1.2)),
          axis.title = element_text(size=rel(1.5)),
          legend.text = element_text(size=rel(1.1)),
          legend.title = element_text(size=rel(1.3)),
          plot.title = element_text(size=rel(1.5)))
  
  ## Exceedance map
  p_meuse + ggtitle(" Probability of exceeding 500 ppm zinc concentrations")+
    # geom_point(data = pred_results, aes(x=loc.x,y= loc.y, col = pexceed))+
    stat_contour_filled(data = pred_grd_meuse, aes(x=Longitude,y=Latitude, z=p_exceed), alpha=0.85)+
    labs(x= "Longitude (degrees)", y="Latitude (degrees)")+scale_fill_viridis_d(name="Probability")+ 
    theme(axis.text = element_text(size=rel(1.2)),
          axis.title = element_text(size=rel(1.5)),
          legend.text = element_text(size=rel(1.1)),
          legend.title = element_text(size=rel(1.3)),
          plot.title = element_text(size=rel(1.5)))
  
  head(pred_grd_meuse)
 
  #The meuse data set provided by package sp is a data set comprising of four
  #heavy metals measured in the top soil in a flood plain along the river Meuse,
  #along with a handful of covariates. The process governing heavy metal distribution seems that polluted sediment is carried by the river, and mostly deposited
  #close to the river bank, and areas with low elevation.
  
  
  
  plot(meuse.grid$x, meuse.grid$y)
plot(pred_grd_meuse_m)
plot(meuse.area)  
pred_grd_meuse_m[1,]
