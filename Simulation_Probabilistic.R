
############ Simulation Study: tukey g and h random fields ################
###########################################################################
library(copula)
library(fields)
library(sp)
library(mvtnorm)

### Dtgh1: function to simulate Tukey g and h random fields

Dtgh1 = function(nx=10, theta=0.5, g, h)
{
  #generate data from gaussian field with exponential correlation not too close to each other
  #nx: no of grids, sample size=nx^2
  #theta: strength of autocorrelation, larger theta means higher autocorrelation
  n = nx^2
  grid = expand.grid(1:nx, 1:nx)
  xy = matrix(runif(n*2, -0.4, 0.4), n, 2)
  loc = (grid-0.5+xy)/nx
  x = loc[,1]
  y = loc[,2]
  
  
  # distance matrix
  dm <- as.matrix(dist(loc))
  # create a correlation structure (exponential)
  omega1 <- exp(-dm/theta)
  # create an autocorrelated random field
  z <- t(chol(omega1))%*%rnorm(n)
  
  ## Apply Tukey g and h transformation
  if(g==0){tgh = z*exp(h*z^2/2)}
  else {tgh = (exp(g*z)-1)*exp(h*z^2/2)/g}
  #quilt.plot(loc, tgh, nx=50, ny=50, zlim = c(-4,4))
  out = data.frame(x, y, tgh)
  return(list(simulation=out, Sigma22=omega1))
}




### Simulate data from Tukey g and h
set.seed(16739)
theta = 0.25
sim_size = 100

g_values = c(0.5, 0.75, 0)
h_values = c(0.3, 0, 0.5)
train_data_list =  vector(mode = "list", length = 3)
test_data_list =  vector(mode = "list", length = 3)
data_list = vector(mode = "list", length = 3)

n=100
dt = sort(sample(n, n*.8, replace = F))

for (c in 1:length(g_values)) {
  
  for(i in 1:sim_size){
    
    g= g_values[c]
    h= h_values[c]
    
    dat = dat2 = Dtgh1(nx=10, theta=theta, g, h)
    dat = dat2[[1]]
    
    
    
    data_list[[c]][[i]] = dat
    train_data_list[[c]][[i]] = dat[dt,]
    test_data_list[[c]][[i]] = dat[-dt,]
  }
}  


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
  
  myCop <- normalCopula(param=params, dim = n, dispstr = "un")
  u = ecdf(Y)(Y)  ### extract the empirical distibution functions value
  
  #truncate u to avoid u equals exactly zero or 1
  u = pmax(0.001,pmin(u, 0.999))
  
  #use true marginal
  #u = pnorm(Y)
  
  logl = log(dCopula(u, myCop))
  return(-logl)
}




library(doParallel)
registerDoParallel(cores=50)

start.time= Sys.time()
parallel_sim = list(); dat_gh= list()
for (c in 1:length(g_values)) {
  
  if(c==1){ qseq = c(seq(0.025,0.7, length.out = 10) , seq(0.8,0.975, length.out = 5))## for skew and ht
  } else if(c==2) { qseq = c(seq(0.05,0.7, length.out = 10) , seq(0.8,0.975, length.out = 5)) ## for positively skewed)
  }  else { qseq = c(0.025,seq(0.05,0.95, length.out = 13),0.975)} ## for heavy tail
  
  parallel_sim[[c]] = foreach(s=1:sim_size)  %dopar% {
    
  
    dat = data_list[[c]][[s]]
    train_data = train_data_list[[c]][[s]]
    test_data = test_data_list[[c]][[s]]
 
    Y = train_data$tgh
    coord = train_data[,1:2]
    
    coordinates(train_data) = ~x+y
    
    theta.hat <- optim(par=0.25,logL.Gcopula.Exp,Y=Y,coord=coord, method="L-BFGS-B",
                       lower=0.01,upper=1)$par

  ################################# Mutiple thresholds ###########################################
    th = quantile(Y, qseq) 
    
    #plot(ecdf(dat$tgh), main="Global ECDF of RR conc.")
    #abline(v=th, lty=2)
    
    dist_mat =  as.matrix(dist(coord))
    h = dist_mat[lower.tri(dist_mat)]
    
    totC<-array(NA,dim = c(dim(dist_mat)[1],dim(dist_mat)[2],length(th),length(th)))
    
    myCop = vector(mode = "list", length = length(h))
    for (k in 1:length(h)) {
      myCop[[k]] =  normalCopula(param=exp(-h[k]/theta.hat), dim = 2, dispstr = "un")
    }
    
    for (i in 1:length(th)) {
      for (j in 1:i) {
        
        #rhat: estimated correlation based on copula
        uhat_i = mean(Y< th[i]) # E[I(Y(s)<th_i)] = P[Y(s)<th_i] = F_Y(th_i)
        uhat_j = mean(Y< th[j])
        
        b=vector()
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
        C = matrix(NA, nrow = nrow(dist_mat), ncol=ncol(dist_mat))
        C[lower.tri(C)] <- Cov #put elements in lower triangular
        diag(C) = Cov_diag
        C = t(C) ## make it upper traingular
        C[lower.tri(C)] <- Cov
        
        
        totC[,,i,j] = C
      }
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
    
    inv_C = solve(myC)
    
    Fsim_gh_list = list()
    library(fields)
    for (l  in 1:dim(test_data)[1]) {
      
      newloc= test_data[l,1:2]
      crosslocs<-rdist(newloc,coord)
      
      totC12<-array(NA,dim = c(dim(crosslocs)[1],dim(crosslocs)[2],length(th),length(th)))
      
      for(k in 1:length(crosslocs))
      {
        myCop[[k]] <- normalCopula(param=exp(-crosslocs[k]/theta.hat), dim = 2, dispstr = "un")
      }
      for (i in 1:length(th)) {
        for (j in 1:length(th)) {
          
          uhat_i = mean(Y< th[i]) # E[I(Y(s)<th_i)] = P[Y(s)<th_i] = F_Y(th_i)
          uhat_j = mean(Y< th[j])
          
          b=vector()
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
      data_ind = list()
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
      xseq = c(xseq, (max(dat$tgh)+xseq[length(xseq)])/2) }
      
      if(F_ind[1]==0)
      {F_ind=F_ind} else {F_ind= c(0,F_ind)
      xseq = c( (min(dat$tgh)+xseq[1])/2 ,xseq)}
      
      ## Order relation correction (Carr 1994)
      ##upward downward
      
      for (i in 1:(length(F_ind)-1)) {
        if(F_ind[i+1]< F_ind[i]){
          m = which(F_ind[(i+1):length(F_ind)] >= F_ind[i])[1]
          
          for (j in (i+1):(i+m-1)) {
            F_ind[j] = F_ind[i] + ((F_ind[i+m]-F_ind[i])*(j-i))/m
          }
        }
      }
      
      Fsim_gh_list[[l]] = cbind(xseq, F_ind)
    }
    return(list(Fsim =Fsim_gh_list, tv= test_data[,3], data_tgh= dat$tgh ))
  }
}
end.time= Sys.time()
end.time-start.time



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



L=200; 
crpsF_y_a = matrix(0, nrow = sim_size, ncol = length(test_data_list[[1]][[1]]$tgh))
pit_sim_a =  vector(mode = "list", length = sim_size) ; rmse_sim_a=vector(); 
rmse_sim_mean_a = vector(); mad_sim_a=vector();  smoothF_gh = list()

### Marginal dsitribution: Skewed and ht distribution
c=1
x_seq = seq(min(unlist(data_list[[1]])), max(unlist(data_list[[1]])), length.out = 10000) #skew and ht
x_mid = x_seq[-length(x_seq)] + diff(x_seq)/2
delta_x = diff(x_seq)[1]

for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim[[c]][[i]]$Fsim, function(a) smoothglobal(Fx=a, globalF = parallel_sim[[c]][[i]]$data_tgh) )
  median_pred_sim_a = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
 
  quantile_pred_sim_l = matrix(0, nrow=L, ncol = length(median_pred_sim_a))
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_sim_l[l,] = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  
  tv_gh = parallel_sim[[c]][[i]]$tv
  
  crpsF_y_iter= vector()
  for (j  in 1:length(tv_gh)) {
    
    pit_sim_a[[i]][j] =   mean(quantile_pred_sim_l[,j]<=tv_gh[j])
    f_x = function(x,y) { {mean(quantile_pred_sim_l[,j]<=x) - ifelse( y <= x,1,0)}^2 }
    fx_mid = sapply(x_mid, function(x) f_x(x, tv_gh[j]))
    crpsF_y_iter[j] = sum(delta_x*fx_mid, na.rm = T)
  }
  crpsF_y_a[i,] =  crpsF_y_iter
}

hist(unlist(pit_sim_a), freq=F, xlab = "Proability Integral Transform", main=paste("Gaussian Copula (g=0.4, h=0.2) "), breaks = 10)
abline(h=1)
mean(apply(crpsF_y_a, 2, median))
mean(rmse_sim_a); mean(mad_sim_a)


### Marginal dsitribution: Skewed distribution (g=0.5, h=0)
L=200; 
crpsF_y_b = matrix(0, nrow = sim_size, ncol = length(test_data_list[[1]][[1]]$tgh))
smoothF_gh = list()
pit_sim_b =  vector(mode = "list", length = sim_size) ; rmse_sim_b=vector(); 
rmse_sim_mean_b = vector(); mad_sim_b=vector()

c=2

x_seq = seq(min(unlist(data_list[[2]])), max(unlist(data_list[[2]])), length.out = 10000) #skewed
x_mid = x_seq[-length(x_seq)] + diff(x_seq)/2
delta_x = diff(x_seq)[1]


for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim[[c]][[i]]$Fsim, function(a) smoothglobal(Fx=a, globalF = parallel_sim[[c]][[i]]$data_tgh) )
  median_pred_sim_b = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
  
  quantile_pred_sim_l = matrix(0, nrow=L, ncol = length(median_pred_sim_b))
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_sim_l[l,] = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  

  tv_gh = parallel_sim[[c]][[i]]$tv
  
  crpsF_y_iter= vector()
  for (j  in 1:length(tv_gh)) {
    
    pit_sim_b[[i]][j] =   mean(quantile_pred_sim_l[,j]<=tv_gh[j])
    f_x = function(x,y) { {mean(quantile_pred_sim_l[,j]<=x) - ifelse( y <= x,1,0)}^2 }
    fx_mid = sapply(x_mid, function(x) f_x(x, tv_gh[j]))
    crpsF_y_iter[j] = sum(delta_x*fx_mid, na.rm = T)
  }
  crpsF_y_b[i,] =  crpsF_y_iter
}

hist(unlist(pit_sim_b), freq=F, xlab = "Proability Integral Transform", main=paste("Gaussian Copula (g=0.5, h=0) "), breaks = 10)
abline(h=1)
mean(apply(crpsF_y_b, 2, median))
mean(rmse_sim_b); mean(mad_sim_b)


### Marginal distribution: Heavy tailed distribution (g=0, h=0.5)
L=200; 
crpsF_y_c = matrix(0, nrow = sim_size, ncol = length(test_data_list[[1]][[1]]$tgh))
smoothF_gh = list()
pit_sim_c =  vector(mode = "list", length = sim_size) ; rmse_sim_c=vector(); 
rmse_sim_mean_c = vector(); mad_sim_c=vector()

c=3

x_seq = seq(min(unlist(data_list[[3]])), max(unlist(data_list[[3]])), length.out = 10000) #ht
x_mid = x_seq[-length(x_seq)] + diff(x_seq)/2
delta_x = diff(x_seq)[1]

for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim[[c]][[i]]$Fsim, function(a) smoothglobal(Fx=a, globalF = parallel_sim[[c]][[i]]$data_tgh) )
  median_pred_sim_c = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
  
  quantile_pred_sim_l = matrix(0, nrow=L, ncol = length(median_pred_sim_c))
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_sim_l[l,] = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  
  tv_gh = parallel_sim[[c]][[i]]$tv
  
  crpsF_y_iter= vector()
  for (j  in 1:length(tv_gh)) {
    
    pit_sim_c[[i]][j] =   mean(quantile_pred_sim_l[,j]<=tv_gh[j])
    f_x = function(x,y) { {mean(quantile_pred_sim_l[,j]<=x) - ifelse( y <= x,1,0)}^2 }
    fx_mid = sapply(x_mid, function(x) f_x(x, tv_gh[j]))
    crpsF_y_iter[j] = sum(delta_x*fx_mid, na.rm = T)
  }
  crpsF_y_c[i,] =  crpsF_y_iter
}

hist(unlist(pit_sim_c), freq=F, xlab = "Proability Integral Transform", main=paste("Gaussian Copula (g=0,h=0.5) "), breaks = 10)
abline(h=1)
mean(apply(crpsF_y_c, 2, median))
mean(rmse_sim_c); mean(mad_sim_c)








##################### Gaussian prediction for Tukey g and h process ####################################
## assuming that the distribution is gaussian
library(fields)
library(geoR)
library(scoringRules)


gaussian_pred = list(); gaussian_var = list(); tv = list(); parallel_g= list()
for (c in 1:length(g_values)) {
  
  for(i in 1:sim_size){
    
    train_data = train_data_list[[c]][[i]]
    test_data = test_data_list[[c]][[i]]
    
    train_loc = train_data[,1:2]
    test_loc = test_data[,1:2]
    
    if(c==1){
    ml <- likfit(data=train_data$tgh,coords=train_loc,
                   fix.nugget=F,cov.model="exponential",
                   ini.cov.pars  = c(0.84, 0.1), nugget = 0.01)
    }
    else if( c==2){
    ml <- likfit(data=train_data$tgh,coords=train_loc,
                 fix.nugget=F,cov.model="exponential",
                 ini.cov.pars  = c(0.98, 0.23), nugget = 0.01)
    }
    else if(c==3){
    ml <- likfit(data=train_data$tgh,coords=train_loc,
                 fix.nugget=F,cov.model="exponential",
                 ini.cov.pars  = c(4, 0.1), nugget = 0.01)
  }

    coordinates(train_data) = ~x+y
    coordinates(test_data) = ~x+y
    
    pred<-krige.conv(data= train_data$tgh, coords= train_loc, locations=test_loc,
                     krige=krige.control(cov.model="exponential", type.krige = "SK", beta = mean(train_data$tgh),
                                         cov.pars=c(ml$sigmasq ,ml$phi),
                                         nugget= ml$nugget))

    #alpha = 0.05
    #l = qnorm(alpha/2,conditional_mean,sqrt(conditional_variance))
    #u = qnorm(1-alpha/2,conditional_mean,sqrt(conditional_variance))
    gaussian_pred[[i]] = list(tv = test_data$tgh, predicted = pred$predict, var = pred$krige.var)
    
}
  parallel_g[[c]] = gaussian_pred
}


#### pit and crps


### Gaussian marginal distribution (g=0.4,h=0.2)
c=1
pit_g_a =  vector(mode = "list", length = sim_size)
rmse_g_a= vector(); mad_g_a = vector()
crps_iter = vector(); crps_g_a = matrix(0, nrow = sim_size, ncol = length(test_data_list[[1]][[1]]$tgh))
for (i in 1:sim_size) {
  tv = parallel_g[[c]][[i]]$tv
  predicted = parallel_g[[c]][[i]]$predicted
  var = parallel_g[[c]][[i]]$var
  
  for (j  in 1:length(tv)) {
    pit_g_a[[i]][j] = pnorm(tv[j], mean= predicted[j], sd = sqrt(var[j]))
    crps_iter[j] =   crps(y = tv[j], family = "normal",  mean= predicted[j], sd = sqrt(var[j]))
  }
  crps_g_a[i,] =  crps_iter
  rmse_g_a[i] = sqrt(mean((tv-predicted)^2)) ## RMSE
  mad_g_a[i] = mean(abs(tv-predicted)) ## MAD
}
hist(unlist(pit_g_a), freq=F, main="Gaussian prediction (g=0.4, h=0.2)", breaks=10)
abline(h=1)    
mean(rmse_g_a); mean(mad_g_a)
mean(apply(crps_g_a,2, mean))



### Skewed marginal distribution (g=0.5,h=0)
c=2
pit_g_b =  vector(mode = "list", length = sim_size)
rmse_g_b= vector(); mad_g_b = vector()
crps_iter = vector(); crps_g_b = matrix(0, nrow = sim_size, ncol = length(test_data_list[[1]][[1]]$tgh))
for (i in 1:sim_size) {
  tv = parallel_g[[c]][[i]]$tv
  predicted = parallel_g[[c]][[i]]$predicted
  var = parallel_g[[c]][[i]]$var
  
  for (j  in 1:length(tv)) {
    pit_g_b[[i]][j] = pnorm(tv[j], mean= predicted[j], sd = sqrt(var[j]))
    crps_iter[j] =   crps(y = tv[j], family = "normal",  mean= predicted[j], sd = sqrt(var[j]))
  }
  crps_g_b[i,] =  crps_iter
  rmse_g_b[i] = sqrt(mean((tv-predicted)^2)) ## RMSE
  mad_g_b[i] = mean(abs(tv-predicted)) ## MAD
}
hist(unlist(pit_g_b), freq=F, main="Gaussian prediction (g=0.5, h=0)", breaks=10)
abline(h=1)    
mean(rmse_g_b); mean(mad_g_b)
mean(apply(crps_g_b, 2, median))





### Heavy tailed marginal distribution (g=0,h=0.5)
c=3
pit_g_c =  vector(mode = "list", length = sim_size)
rmse_g_c= vector(); mad_g_c = vector()
crps_iter = vector(); crps_g_c = matrix(0, nrow = sim_size, ncol = length(test_data_list[[1]][[1]]$tgh))
for (i in 1:sim_size) {
  tv = parallel_g[[c]][[i]]$tv
  predicted = parallel_g[[c]][[i]]$predicted
  var = parallel_g[[c]][[i]]$var
  
  for (j  in 1:length(tv)) {
    pit_g_c[[i]][j] = pnorm(tv[j], mean= predicted[j], sd = sqrt(var[j]))
    crps_iter[j] =   crps(y = tv[j], family = "normal",  mean= predicted[j], sd = sqrt(var[j]))
  }
  crps_g_c[i,] =  crps_iter
  rmse_g_c[i] = sqrt(mean((tv-predicted)^2)) ## RMSE
  mad_g_c[i] = mean(abs(tv-predicted)) ## MAD
}
hist(unlist(pit_g_c), freq=F, main="Gaussian prediction (g=0, h=0.5)", breaks = 10)
abline(h=1)    
mean(rmse_g_c); mean(mad_g_c)
mean(apply(crps_g_c,2, median))


######################################## t copula ###########################################

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
  nu= par[2] #
  #nu= 4
  logl <- dmvt(qt(u, df = nu), df=nu, sigma=R.m,log = TRUE) - sum(dt(qt(u, df=nu), df=nu, log=TRUE)) # log-density for X2(s) (denominator)
  
  #logl = log(dCopula(u, myCop))
  return(-logl)
}



library(doParallel)
registerDoParallel(cores=50)

start.time= Sys.time()
parallel_sim_t = list(); dat_gh= list()
for (c in 1:length(g_values)) {
  
  if(c==1){ qseq = c(seq(0.025,0.7, length.out = 10) , seq(0.8,0.975, length.out = 5))## for skew and ht
  } else if(c==2) { qseq = c(seq(0.05,0.7, length.out = 10) , seq(0.8,0.975, length.out = 5)) ## for positively skewed)
  }  else { qseq = c(0.025,seq(0.05,0.95, length.out = 13),0.975)} ## for heavy tail
  
  parallel_sim_t[[c]] = foreach(s=1:sim_size)  %dopar% {
    
    
    dat = data_list[[c]][[s]]
    train_data = train_data_list[[c]][[s]]
    test_data = test_data_list[[c]][[s]]
    
    Y = train_data$tgh
    coord = train_data[,1:2]
    
    coordinates(train_data) = ~x+y
    
    theta.hat_vec <- optim(par = c(0.25,10),logL.tcopula.Exp,Y=Y,coord=coord, method="L-BFGS-B",
                           lower=0.01,upper=30)$par
    
    theta.hat = theta.hat_vec[1]
    nu = round(theta.hat_vec[2],0)  #4
    
    #theta.hat; nu
    ################################# Mutiple thresholds ###########################################
    
    th = quantile(Y, qseq) 
    
    #plot(ecdf(dat$tgh), main="Global ECDF of RR conc.")
    #abline(v=th, lty=2)
    
    dist_mat =  as.matrix(dist(coord))
    h = dist_mat[lower.tri(dist_mat)]
    
    totC<-array(NA,dim = c(dim(dist_mat)[1],dim(dist_mat)[2],length(th),length(th)))
    
    myCop = vector(mode = "list", length = length(h))
    for (k in 1:length(h)) {
      myCop[[k]] <- tCopula(param=exp(-h[k]/theta.hat), dim = 2, dispstr = "un", df=nu)
    }
    
    for (i in 1:length(th)) {
      for (j in 1:i) {
        
        #rhat: estimated correlation based on copula
        uhat_i = mean(Y< th[i]) # E[I(Y(s)<th_i)] = P[Y(s)<th_i] = F_Y(th_i)
        uhat_j = mean(Y< th[j])
        
        b=vector()
        for(k in 1:length(h))
        {
          #estimated based on copula
          #myCop <- normalCopula(param=exp(-h[k]/theta.hat), dim = 2, dispstr = "un")
          b[k] = pCopula(c(uhat_i,uhat_j), myCop[[k]]) #P(U1<u, U2<u)
        }
        Cov = b - uhat_i*uhat_j
        myCop_diag <- tCopula(param=exp(-0/theta.hat), dim = 2, dispstr = "un", df=nu)
        Cov_diag = pCopula(c(uhat_i,uhat_j), myCop_diag) - uhat_i*uhat_j
        
        ## Construct the covariace matrix
        C = matrix(NA, nrow = nrow(dist_mat), ncol=ncol(dist_mat))
        C[lower.tri(C)] <- Cov #put elements in lower triangular
        diag(C) = Cov_diag
        C = t(C) ## make it upper traingular
        C[lower.tri(C)] <- Cov
        
        totC[,,i,j] = C
      }
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
    
    inv_C = solve(myC)
    
    Fsim_gh_list = list()
    library(fields)
    for (l  in 1:dim(test_data)[1]) {
      
      newloc= test_data[l,1:2]
      crosslocs<-rdist(newloc,coord)
      
      totC12<-array(NA,dim = c(dim(crosslocs)[1],dim(crosslocs)[2],length(th),length(th)))
      
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
      data_ind = list()
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
      xseq = c(xseq, (max(dat$tgh)+xseq[length(xseq)])/2) }
      
      if(F_ind[1]==0)
      {F_ind=F_ind} else {F_ind= c(0,F_ind)
      xseq = c( (min(dat$tgh)+xseq[1])/2 ,xseq)}
      
      ## Order relation correction (Carr 1994)
      ##upward downward
      
      for (i in 1:(length(F_ind)-1)) {
        if(F_ind[i+1]< F_ind[i]){
          m = which(F_ind[(i+1):length(F_ind)] >= F_ind[i])[1]
          
          for (j in (i+1):(i+m-1)) {
            F_ind[j] = F_ind[i] + ((F_ind[i+m]-F_ind[i])*(j-i))/m
          }
        }
      }
      
      Fsim_gh_list[[l]] = cbind(xseq, F_ind)
      
    }
    return(list(Fsim =Fsim_gh_list, tv= test_data[,3], data_tgh= dat$tgh ))
  }
}
end.time= Sys.time()
end.time-start.time



#### marginal distribution : Skew and ht

L=200; 
crpsF_y_t_a = matrix(0, nrow = sim_size, ncol = length(test_data_list[[1]][[1]]$tgh));  smoothF_gh = list()
pit_sim_t_a =  vector(mode = "list", length = sim_size) ; rmse_sim_t_a=vector(); mad_sim_t_a=vector()

c=1
x_seq = seq(min(unlist(data_list[[1]])), max(unlist(data_list[[1]])), length.out = 10000) #skew and ht
x_mid = x_seq[-length(x_seq)] + diff(x_seq)/2
delta_x = diff(x_seq)[1]  

for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim_t[[c]][[i]]$Fsim, function(a) smoothglobal(Fx=a, globalF = parallel_sim_t[[c]][[i]]$data_tgh) )
  median_pred_sim_t_a = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
  
  quantile_pred_sim_l = matrix(0, nrow=L, ncol = length(median_pred_sim_t_a))
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_sim_l[l,] = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  
  tv_gh = parallel_sim_t[[c]][[i]]$tv
  
  crpsF_y_iter= vector()
  for (j  in 1:length(tv_gh)) {
    
    pit_sim_t_a[[i]][j] =   mean(quantile_pred_sim_l[,j]<=tv_gh[j])
    f_x = function(x,y) { {mean(quantile_pred_sim_l[,j]<=x) - ifelse( y <= x,1,0)}^2 }
    fx_mid = sapply(x_mid, function(x) f_x(x, tv_gh[j]))
    crpsF_y_iter[j] = sum(delta_x*fx_mid, na.rm = T)
  }
  crpsF_y_t_a[i,] =  crpsF_y_iter
}

hist(unlist(pit_sim_t_a), freq=F, xlab = "Proability Integral Transform", main=paste("t-Copula (g=0, h=0)"), breaks=10)
abline(h=1)
mean(apply(crpsF_y_t_a, 2,median))
#mean(rmse_sim_t_a); mean(mad_sim_t_a)


#### marginal distribution : Skewed (g=0.5, h=0)

L=200; 
crpsF_y_t_b = matrix(0, nrow = sim_size, ncol = length(test_data_list[[1]][[1]]$tgh));  smoothF_gh = list()
pit_sim_t_b =  vector(mode = "list", length = sim_size) ; rmse_sim_t_b=vector(); mad_sim_t_b=vector()

c=2
x_seq = seq(min(unlist(data_list[[2]])), max(unlist(data_list[[2]])), length.out = 10000) #skew 
x_mid = x_seq[-length(x_seq)] + diff(x_seq)/2
delta_x = diff(x_seq)[1]  

for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim_t[[c]][[i]]$Fsim, function(a) smoothglobal(Fx=a, globalF = parallel_sim_t[[c]][[i]]$data_tgh) )
  median_pred_sim_t_b = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
  
  quantile_pred_sim_l = matrix(0, nrow=L, ncol = length(median_pred_sim_t_b))
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_sim_l[l,] = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  
  tv_gh = parallel_sim_t[[c]][[i]]$tv
  
  crpsF_y_iter= vector()
  for (j  in 1:length(tv_gh)) {
    
    pit_sim_t_b[[i]][j] =   mean(quantile_pred_sim_l[,j]<=tv_gh[j])
    f_x = function(x,y) { {mean(quantile_pred_sim_l[,j]<=x) - ifelse( y <= x,1,0)}^2 }
    fx_mid = sapply(x_mid, function(x) f_x(x, tv_gh[j]))
    crpsF_y_iter[j] = sum(delta_x*fx_mid, na.rm = T)
  }
  crpsF_y_t_b[i,] =  crpsF_y_iter
}

hist(unlist(pit_sim_t_b), freq=F, xlab = "Proability Integral Transform",
     main=paste("t-Copula (g=0.5, h=0)"), breaks=10)
abline(h=1)
mean(apply(crpsF_y_t_b, 2, median))
#mean(rmse_sim_t_b); mean(mad_sim_t_b)




#### marginal distribution : heavy tailed (g=0, h=0.5)

L=200; 
crpsF_y_t_c = matrix(0, nrow = sim_size, ncol = length(test_data_list[[1]][[1]]$tgh));  smoothF_gh = list()
pit_sim_t_c =  vector(mode = "list", length = sim_size) ; rmse_sim_t_c=vector(); mad_sim_t_c=vector()

c=3
x_seq = seq(min(unlist(data_list[[3]])), max(unlist(data_list[[3]])), length.out = 10000) # ht
x_mid = x_seq[-length(x_seq)] + diff(x_seq)/2
delta_x = diff(x_seq)[1]  

for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim_t[[c]][[i]]$Fsim, function(a) smoothglobal(Fx=a, globalF = parallel_sim_t[[c]][[i]]$data_tgh) )
  median_pred_sim_t_c = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
  
  quantile_pred_sim_l = matrix(0, nrow=L, ncol = length(median_pred_sim_t_c))
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_sim_l[l,] = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  
  tv_gh = parallel_sim_t[[c]][[i]]$tv
  
  crpsF_y_iter= vector()
  for (j  in 1:length(tv_gh)) {
    
    pit_sim_t_c[[i]][j] =   mean(quantile_pred_sim_l[,j]<=tv_gh[j])
    f_x = function(x,y) { {mean(quantile_pred_sim_l[,j]<=x) - ifelse( y <= x,1,0)}^2 }
    fx_mid = sapply(x_mid, function(x) f_x(x, tv_gh[j]))
    crpsF_y_iter[j] = sum(delta_x*fx_mid, na.rm = T)
  }
  crpsF_y_t_c[i,] =  crpsF_y_iter
}

hist(unlist(pit_sim_t_c), freq=F, xlab = "Proability Integral Transform", main=paste("t-Copula (g=0, h=0.5)"), breaks=10)
abline(h=1)
mean(apply(crpsF_y_t_c, 2, median))
#mean(rmse_sim_t_c); mean(mad_sim_t_c)


###################################### Variogram ###################################
Fsimv_gh_list = list(); tv_gh_list = list(); data_tgh_list= list()
parallel_sim_var= list()

library(gstat)
for (k in 1:length(g_values)) {
  
  if(c==1){ qseq = c(seq(0.025,0.7, length.out = 10) , seq(0.8,0.975, length.out = 5))## for skew and ht
  } else if(c==2) { qseq = c(seq(0.05,0.7, length.out = 10) , seq(0.8,0.975, length.out = 5)) ## for positively skewed)
  }  else { qseq = c(0.025,seq(0.05,0.95, length.out = 13),0.975)} ## for heavy tail
  
  for(i in 1:sim_size){
 
    dat = data_list[[k]][[i]]
    train_data = train_data_list[[k]][[i]]
    test_data = test_data_list[[k]][[i]]
    
    coordinates(train_data) = ~x+y
    coordinates(test_data) = ~x+y
    
    ### indicator cokriging
    quartz <- quantile(train_data$tgh, qseq)
    # These are the actual values of the quantiles
    
    tgh.i <- gstat(id = "tgh1", formula = I(tgh < quartz[1]) ~ 1, data = train_data, 
                   nmax = 7, beta = qseq[1], set = list(order = 4, zero = 1e-05))
    tgh.i <- gstat(tgh.i, "tgh2", formula = I(tgh < quartz[2]) ~ 1, data = train_data, 
                   nmax = 7, beta = qseq[2])
    tgh.i <- gstat(tgh.i, "tgh3", formula = I(tgh < quartz[3]) ~ 1, data = train_data, 
                   nmax = 7, beta = qseq[3])
    tgh.i <- gstat(tgh.i, "tgh4", formula = I(tgh < quartz[4]) ~ 1, data = train_data, 
                   nmax = 7, beta = qseq[4])
    tgh.i <- gstat(tgh.i, "tgh5", formula = I(tgh < quartz[5]) ~ 1, data = train_data, 
                   nmax = 7, beta = qseq[5])
    tgh.i <- gstat(tgh.i, "tgh6", formula = I(tgh < quartz[6]) ~ 1, data = train_data, 
                   nmax = 7, beta = qseq[6])
    tgh.i <- gstat(tgh.i, "tgh7", formula = I(tgh < quartz[7]) ~ 1, data = train_data, 
                   nmax = 7, beta = qseq[7])
    tgh.i <- gstat(tgh.i, "tgh8", formula = I(tgh < quartz[8]) ~ 1, data = train_data, 
                   nmax = 7, beta = qseq[8])
    tgh.i <- gstat(tgh.i, "tgh9", formula = I(tgh < quartz[9]) ~ 1, data = train_data, 
                   nmax = 7, beta = qseq[9])
    tgh.i <- gstat(tgh.i, "tgh10", formula = I(tgh < quartz[10]) ~ 1, data = train_data, 
                   nmax = 7, beta = qseq[10])
    tgh.i <- gstat(tgh.i, "tgh11", formula = I(tgh < quartz[11]) ~ 1, data = train_data, 
                   nmax = 7, beta = qseq[11])
    tgh.i <- gstat(tgh.i, "tgh12", formula = I(tgh < quartz[12]) ~ 1, data = train_data, 
                   nmax = 7, beta = qseq[12])
    tgh.i <- gstat(tgh.i, "tgh13", formula = I(tgh < quartz[13]) ~ 1, data = train_data, 
                   nmax = 7, beta = qseq[13])
    tgh.i <- gstat(tgh.i, "tgh14", formula = I(tgh < quartz[14]) ~ 1, data = train_data, 
                   nmax = 7, beta = qseq[14])
    tgh.i <- gstat(tgh.i, "tgh15", formula = I(tgh < quartz[15]) ~ 1, data = train_data, 
                   nmax = 7, beta = qseq[15])
    
    tgh.i <- gstat(tgh.i, model = vgm(psill=0.04, "Exp", range =  0.3, nugget=0.01), fill.all = T)
    
    
    tgh.vg <- variogram(tgh.i)
    tgh.quartfit = fit.lmc(tgh.vg, tgh.i, fit.ranges = F)
    
    # now do the co-kriging
    cokrig.tghquart <- predict(tgh.quartfit, test_data[,-3])
    F_test_list = list()
    for (j in 1:dim(test_data)[1]) {
      
      
      Fsim = c(cokrig.tghquart@data$tgh1.pred[j],cokrig.tghquart@data$tgh2.pred[j],cokrig.tghquart@data$tgh3.pred[j],
               cokrig.tghquart@data$tgh4.pred[j],cokrig.tghquart@data$tgh5.pred[j],
               cokrig.tghquart@data$tgh6.pred[j],cokrig.tghquart@data$tgh7.pred[j],cokrig.tghquart@data$tgh8.pred[j],cokrig.tghquart@data$tgh9.pred[j],
               cokrig.tghquart@data$tgh10.pred[j],cokrig.tghquart@data$tgh11.pred[j],cokrig.tghquart@data$tgh12.pred[j],
               cokrig.tghquart@data$tgh13.pred[j],cokrig.tghquart@data$tgh14.pred[j],cokrig.tghquart@data$tgh15.pred[j])
      
      
      ############ Extrapolate to 0 and 1  ############
      xseq= quartz
      if(Fsim[length(quartz)]==1)
      {Fsim=Fsim} else {Fsim = c(Fsim,1)
      xseq = c(xseq,(max(train_data$tgh)+xseq[length(xseq)])/2)}
      
      if(Fsim[1]==0)
      {Fsim=Fsim} else {Fsim= c(0,Fsim)
      xseq = c( (min(train_data$tgh)+xseq[1])/2 ,xseq)}
      F_test_list[[j]] = cbind(xseq, Fsim)
    }
    
    Fsimv_gh_list[[i]] = F_test_list
    tv_gh_list[[i]] = test_data$tgh
    data_tgh_list[[i]] = dat$tgh
  }
  parallel_sim_var[[k]]= list(Fsim =Fsimv_gh_list, tv= tv_gh_list, data_tgh= data_tgh_list )
  
}

#hist(parallel_sim_var[[3]]$data_tgh[[3]])

####################### Marginal distribution: Skew and ht (g=0.4, h=0.2) #################################
L=200; 
crpsF_y_var_a = matrix(0, nrow = sim_size, ncol = length(test_data_list[[1]][[1]]$tgh));  smoothF_gh = list()
pit_sim_var_a =  vector(mode = "list", length = sim_size) ;
rmse_sim_var_a=vector(); mad_sim_var_a=vector()

c=1
x_seq = seq(min(unlist(data_list[[1]])), max(unlist(data_list[[1]])), length.out = 10000) #skew and ht
x_mid = x_seq[-length(x_seq)] + diff(x_seq)/2
delta_x = diff(x_seq)[1]  


for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim_var[[c]]$Fsim[[i]], 
                    function(a) smoothglobal(Fx=a, globalF = parallel_sim_var[[c]]$data_tgh[[i]]) )
  median_pred_sim_var_a = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
  
  quantile_pred_sim_l = matrix(0, nrow=L, ncol = length(median_pred_sim_var_a))
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_sim_l[l,] = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  
  tv_gh = parallel_sim_var[[c]]$tv[[i]]
  
  crpsF_y_iter= vector()
  for (j  in 1:length(tv_gh)) {
    
    pit_sim_var_a[[i]][j] =   mean(quantile_pred_sim_l[,j]<=tv_gh[j])
    f_x = function(x,y) { {mean(quantile_pred_sim_l[,j]<=x) - ifelse( y <= x,1,0)}^2 }
    fx_mid = sapply(x_mid, function(x) f_x(x, tv_gh[j]))
    crpsF_y_iter[j] = sum(delta_x*fx_mid, na.rm = T)
  }
  crpsF_y_var_a[i,] =  crpsF_y_iter
}

hist(unlist(pit_sim_var_a), freq=F, xlab = "Proability Integral Transform",
     main=paste("Variogram (g=0,h=0) "), breaks=10)
abline(h=1)
mean(apply(crpsF_y_var_a, 2,median))
#mean(rmse_sim_var_a); mean(mad_sim_var_a)



#### Marginal distribution: Skewed (g=0.5, h=0)
L=200; 
crpsF_y_var_b = matrix(0, nrow = sim_size, ncol = length(test_data_list[[1]][[1]]$tgh));  smoothF_gh = list()
pit_sim_var_b =  vector(mode = "list", length = sim_size) ;
rmse_sim_var_b=vector(); mad_sim_var_b=vector()

c=2
x_seq = seq(min(unlist(data_list[[2]])), max(unlist(data_list[[2]])), length.out = 10000) #skew
x_mid = x_seq[-length(x_seq)] + diff(x_seq)/2
delta_x = diff(x_seq)[1]  

for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim_var[[c]]$Fsim[[i]], function(a) smoothglobal(Fx=a, globalF = parallel_sim_var[[c]]$data_tgh[[i]]) )
  median_pred_sim_var_b = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
  
  quantile_pred_sim_l = matrix(0, nrow=L, ncol = length(median_pred_sim_var_b))
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_sim_l[l,] = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  
  tv_gh = parallel_sim_var[[c]]$tv[[i]]
  
  crpsF_y_iter= vector()
  for (j  in 1:length(tv_gh)) {
    
    pit_sim_var_b[[i]][j] =   mean(quantile_pred_sim_l[,j]<=tv_gh[j])
    f_x = function(x,y) { {mean(quantile_pred_sim_l[,j]<=x) - ifelse( y <= x,1,0)}^2 }
    fx_mid = sapply(x_mid, function(x) f_x(x, tv_gh[j]))
    crpsF_y_iter[j] = sum(delta_x*fx_mid, na.rm = T)
  }
  crpsF_y_var_b[i,] =  crpsF_y_iter
}

hist(unlist(pit_sim_var_b), freq=F, xlab = "Proability Integral Transform",
     main=paste("Variogram (g=0.5,h=0) "), breaks=10)
abline(h=1)
mean(apply(crpsF_y_var_b, 2, median))
#mean(rmse_sim_var_b); mean(mad_sim_var_b)


#### Marginal distribution: Heavy tailed (g=0, h=0.5)
L=200; 
crpsF_y_var_c = matrix(0, nrow = sim_size, ncol = length(test_data_list[[1]][[1]]$tgh));  smoothF_gh = list()
pit_sim_var_c =  vector(mode = "list", length = sim_size) ;
rmse_sim_var_c=vector(); mad_sim_var_c=vector()

c=3
x_seq = seq(min(unlist(data_list[[3]])), max(unlist(data_list[[3]])), length.out = 10000) #ht
x_mid = x_seq[-length(x_seq)] + diff(x_seq)/2
delta_x = diff(x_seq)[1]  

for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim_var[[c]]$Fsim[[i]], function(a) smoothglobal(Fx=a, globalF = parallel_sim_var[[c]]$data_tgh[[i]]) )
  median_pred_sim_var_c = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
  
  quantile_pred_sim_l = matrix(0, nrow=L, ncol = length(median_pred_sim_var_c))
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_sim_l[l,] = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  
  tv_gh = parallel_sim_var[[c]]$tv[[i]]
  
  crpsF_y_iter= vector()
  for (j  in 1:length(tv_gh)) {
    
    pit_sim_var_c[[i]][j] =   mean(quantile_pred_sim_l[,j]<=tv_gh[j])
    f_x = function(x,y) { {mean(quantile_pred_sim_l[,j]<=x) - ifelse( y <= x,1,0)}^2 }
    fx_mid = sapply(x_mid, function(x) f_x(x, tv_gh[j]))
    crpsF_y_iter[j] = sum(delta_x*fx_mid, na.rm = T)
  }
  crpsF_y_var_c[i,] =  crpsF_y_iter
}

hist(unlist(pit_sim_var_c), freq=F, xlab = "Proability Integral Transform",
     main=paste("Variogram (g=0,h=0.5) "), breaks=10)
abline(h=1)
mean(apply(crpsF_y_var_c, 2, median))
#mean(rmse_sim_var_c); mean(mad_sim_var_c)











##### average length of ci and coverage probability #################
alpha = 0.1
### gaussian kriging
c=1
length_ci_g_a = vector(); cp_g_a = vector()

for (i in 1:sim_size) {
  tv = parallel_g[[c]][[i]]$tv
  predicted = parallel_g[[c]][[i]]$predicted
  var = parallel_g[[c]][[i]]$var
  
  ll = qnorm(alpha/2, mean = predicted, sd = sqrt(var))
  ul = qnorm(1-alpha/2, mean = predicted, sd =sqrt(var))
  length_ci_g_a[i] = mean(ul-ll)
  cp_g_a[i] = sum(tv > ll & tv < ul)/length(tv)
}
mean(cp_g_a)
mean(length_ci_g_a)

c=2
length_ci_g_b = vector(); cp_g_b = vector()

for (i in 1:sim_size) {
  tv = parallel_g[[c]][[i]]$tv
  predicted = parallel_g[[c]][[i]]$predicted
  var = parallel_g[[c]][[i]]$var
  
  ll = qnorm(alpha/2, mean = predicted, sd = sqrt(var))
  ul = qnorm(1-alpha/2, mean = predicted, sd =sqrt(var))
  length_ci_g_b[i] = mean(ul-ll)
  cp_g_b[i] = sum(tv > ll & tv < ul)/length(tv)
}
mean(cp_g_b)
mean(length_ci_g_b)

c=3
length_ci_g_c = vector(); cp_g_c = vector()

for (i in 1:sim_size) {
  tv = parallel_g[[c]][[i]]$tv
  predicted = parallel_g[[c]][[i]]$predicted
  var = parallel_g[[c]][[i]]$var
  
  ll = qnorm(alpha/2, mean = predicted, sd = sqrt(var))
  ul = qnorm(1-alpha/2, mean = predicted, sd =sqrt(var))
  length_ci_g_c[i] = mean(ul-ll)
  cp_g_c[i] = sum(tv > ll & tv < ul)/length(tv)
}
mean(cp_g_c)
mean(length_ci_g_c)


## variogram
############## average length ci and coverage proabability
c=1
length_ci_var_a = vector(); cp_var_a= vector()
for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim_var[[c]]$Fsim[[i]], 
                              function(a) smoothglobal(Fx=a, globalF = parallel_sim_var[[c]]$data_tgh[[i]]) )
  
  ll = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= alpha/2)])}))
  ul = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= 1- alpha/2)])}))
  
  tv_gh = parallel_sim_var[[c]]$tv[[i]]
  length_ci_var_a[i] = mean(ul-ll)
  cp_var_a[i] = sum(tv_gh > ll & tv_gh < ul)/length(tv_gh)
  
}
mean(length_ci_var_a)
mean(cp_var_a)

c=2
length_ci_var_b = vector(); cp_var_b= vector()
for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim_var[[c]]$Fsim[[i]], 
                              function(a) smoothglobal(Fx=a, globalF = parallel_sim_var[[c]]$data_tgh[[i]]) )
  
  ll = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= alpha/2)])}))
  ul = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= 1- alpha/2)])}))
  
  tv_gh = parallel_sim_var[[c]]$tv[[i]]
  length_ci_var_b[i] = mean(ul-ll)
  cp_var_b[i] = sum(tv_gh > ll & tv_gh < ul)/length(tv_gh)
  
}
mean(length_ci_var_b)
mean(cp_var_b)

c=3
length_ci_var_c = vector(); cp_var_c= vector()
for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim_var[[c]]$Fsim[[i]], 
                              function(a) smoothglobal(Fx=a, globalF = parallel_sim_var[[c]]$data_tgh[[i]]) )
  
  ll = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= alpha/2)])}))
  ul = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= 1- alpha/2)])}))
  
  tv_gh = parallel_sim_var[[c]]$tv[[i]]
  length_ci_var_c[i] = mean(ul-ll)
  cp_var_c[i] = sum(tv_gh > ll & tv_gh < ul)/length(tv_gh)
  
}
mean(length_ci_var_c)
mean(cp_var_c)


## Gaussian copula
c=1
length_ci_a = vector(); cp_a= vector()
for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim[[c]][[i]]$Fsim, function(a) smoothglobal(Fx=a, globalF = parallel_sim[[c]][[i]]$data_tgh) )
  
  
  ll = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= alpha/2)])}))
  ul = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= 1- alpha/2)])}))
  
  tv_gh = parallel_sim[[c]][[i]]$tv
  length_ci_a[i] = mean(ul-ll)
  cp_a[i] = sum(tv_gh > ll & tv_gh < ul)/length(tv_gh)
  
}
mean(length_ci_a)
mean(cp_a)

c=2
length_ci_b = vector(); cp_b= vector()
for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim[[c]][[i]]$Fsim, function(a) smoothglobal(Fx=a, globalF = parallel_sim[[c]][[i]]$data_tgh) )
  
  
  ll = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= alpha/2)])}))
  ul = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= 1- alpha/2)])}))
  
  tv_gh = parallel_sim[[c]][[i]]$tv
  length_ci_b[i] = mean(ul-ll)
  cp_b[i] = sum(tv_gh > ll & tv_gh < ul)/length(tv_gh)
  
}
mean(length_ci_b)
mean(cp_b)


c=3
length_ci_c = vector(); cp_c= vector()
for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim[[c]][[i]]$Fsim, function(a) smoothglobal(Fx=a, globalF = parallel_sim[[c]][[i]]$data_tgh) )
  
  
  ll = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= alpha/2)])}))
  ul = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= 1- alpha/2)])}))
  
  tv_gh = parallel_sim[[c]][[i]]$tv
  length_ci_c[i] = mean(ul-ll)
  cp_c[i] = sum(tv_gh > ll & tv_gh < ul)/length(tv_gh)
  
}
mean(length_ci_c)
mean(cp_c)


### t copula

c=1
length_ci_t_a = vector(); cp_t_a= vector()
for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim_t[[c]][[i]]$Fsim, 
                              function(a) smoothglobal(Fx=a, globalF = parallel_sim_t[[c]][[i]]$data_tgh) )
  
  ll = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= alpha/2)])}))
  ul = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= 1- alpha/2)])}))
  
  tv_gh = parallel_sim_t[[c]][[i]]$tv
  length_ci_t_a[i] = mean(ul-ll)
  cp_t_a[i] = sum(tv_gh > ll & tv_gh < ul)/length(tv_gh)
  
}
mean(length_ci_t_a)
mean(cp_t_a)

c=2
length_ci_t_b = vector(); cp_t_b= vector()
for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim_t[[c]][[i]]$Fsim, 
                              function(a) smoothglobal(Fx=a, globalF = parallel_sim_t[[c]][[i]]$data_tgh) )
  
  ll = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= alpha/2)])}))
  ul = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= 1- alpha/2)])}))
  
  tv_gh = parallel_sim_t[[c]][[i]]$tv
  length_ci_t_b[i] = mean(ul-ll)
  cp_t_b[i] = sum(tv_gh > ll & tv_gh < ul)/length(tv_gh)
  
}
mean(length_ci_t_b)
mean(cp_t_b)

c=3
length_ci_t_c = vector(); cp_t_c= vector()
for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim_t[[c]][[i]]$Fsim, 
                              function(a) smoothglobal(Fx=a, globalF = parallel_sim_t[[c]][[i]]$data_tgh) )
  
  ll = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= alpha/2)])}))
  ul = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= 1- alpha/2)])}))
  
  tv_gh = parallel_sim_t[[c]][[i]]$tv
  length_ci_t_c[i] = mean(ul-ll)
  cp_t_c[i] = sum(tv_gh > ll & tv_gh < ul)/length(tv_gh)
  
}
mean(length_ci_t_c)
mean(cp_t_c)







### shortest ci
## the length is short but it decreases the coverage
alpha = 0.1
length_ci_b = vector(); cp_b= vector()
for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim[[c]][[i]]$Fsim, function(a) smoothglobal(Fx=a, globalF = parallel_sim[[c]][[i]]$data_tgh) )
  
  lowlim = uplim = length_gamma = matrix(0, nrow = 100, ncol = 20) 
  
  gamma = seq(0,alpha, length.out=100)
  for (k in 1:length(gamma)) {
    lowlim[k,] =  as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>= gamma[k])])}))
    uplim[k,] =  as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=1-gamma[k])])}))
    length_gamma[k,] = uplim[k,]-lowlim[k,]
  }
  #length_gamma[,2]
  lim_index = apply(length_gamma,2, function(x) which.min(x) )
  
  ll = ul = vector()
  for (j in 1:20) {
    ll[j] =  lowlim[lim_index[j],j]
    ul[j] = uplim[lim_index[j],j]
  }
  tv_gh = parallel_sim[[c]][[i]]$tv
  length_ci_b[i] = mean(ul-ll)
  cp_b[i] = sum(tv_gh > ll & tv_gh < ul)/length(tv_gh)
  
}
mean(length_ci_b)
mean(cp_b)


####################### Figures for manuscript ############################

#  PIT

## Skewed
hist(unlist(pit_g_b), freq=F, main="Gaussian Kriging (Skewed)", breaks=10, ylim=c(0,2), 
     xlab = "Proability Integral Transform", cex.main=1.5, cex.lab=1.4)
abline(h=1)  
hist(unlist(pit_sim_var_b), freq=F, xlab = "Proability Integral Transform",ylim=c(0,2),
     main=paste("Variogram (Skewed) "), breaks=10, cex.main=1.5, cex.lab=1.4)
abline(h=1)
hist(unlist(pit_sim_b), freq=F, xlab = "Proability Integral Transform",ylim=c(0,2),
     main=paste("Gaussian Copula (Skewed) "), breaks = 10, cex.main=1.5, cex.lab=1.4)
abline(h=1)
hist(unlist(pit_sim_t_b), freq=F, xlab = "Proability Integral Transform", ylim=c(0,2),
     main=paste("t-Copula (Skewed)"), breaks=10, cex.main=1.5, cex.lab=1.4)
abline(h=1)


## heavy tailed
hist(unlist(pit_g_c), freq=F, main="Gaussian Kriging (Heavy-tailed)", breaks=10, ylim=c(0,2), 
     xlab = "Proability Integral Transform", cex.main=1.5, cex.lab=1.4)
abline(h=1)  
hist(unlist(pit_sim_var_c), freq=F, xlab = "Proability Integral Transform",  ylim=c(0,2),
     main=paste("Variogram (Heavy-tailed) "), breaks=10, cex.main=1.5, cex.lab=1.4)
abline(h=1)
hist(unlist(pit_sim_c), freq=F, xlab = "Proability Integral Transform", ylim=c(0,2),
     main=paste("Gaussian Copula (Heavy-tailed) "), breaks = 10, cex.main=1.5, cex.lab=1.4)
abline(h=1)
hist(unlist(pit_sim_t_c), freq=F, xlab = "Proability Integral Transform", ylim=c(0,2),
     main=paste("t-Copula (Heavy-tailed)"), breaks=10, cex.main=1.5, cex.lab=1.4)
abline(h=1)

## skewed and heavy tailed
hist(unlist(pit_g_a), freq=F, main="Gaussian Kriging (Skewed & Heavy-tailed)", breaks=10, ylim=c(0,2),
     xlab = "Proability Integral Transform" , cex.main=1.5, cex.lab=1.4)
abline(h=1)  
hist(unlist(pit_sim_var_a), freq=F, xlab = "Proability Integral Transform",ylim=c(0,2),
     main=paste("Variogram (Skewed & Heavy-tailed ) "), breaks=10, cex.main=1.5, cex.lab=1.4)
abline(h=1)
hist(unlist(pit_sim_a), freq=F, xlab = "Proability Integral Transform",ylim=c(0,2),
     main=paste("Gaussian Copula (Skewed & Heavy-tailed) "), breaks = 10, cex.main=1.5, cex.lab=1.4)
abline(h=1)
hist(unlist(pit_sim_t_a), freq=F, xlab = "Proability Integral Transform", ylim=c(0,2),
     main=paste("t-Copula (Skewed & Heavy-tailed)"), breaks=10, cex.main=1.5, cex.lab=1.4)
abline(h=1)



## ggplot

## Skewed
ggplot() + geom_histogram(aes(x=unlist(pit_g_b), y=..density..),
                          breaks=hist(unlist(pit_g_b), plot=F, breaks = 10)$breaks, col="black")+ ylim(0,2)+
  xlab("Proability Integral Transform")+ ylab("Density")+
  ggtitle("GK (skewed)")+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.8)),
        plot.title = element_text(size=rel(2), hjust = 0.5))+
  geom_hline(yintercept = 1)

ggplot() + geom_histogram(aes(x=unlist(pit_sim_var_b), y=..density..),breaks=hist(unlist(pit_sim_var_b), plot=F, breaks = 10)$breaks, col="black")+ ylim(0,2)+
  xlab("Proability Integral Transform")+ ylab("Density")+
  ggtitle("VMIK (skewed)")+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.8)),
        plot.title = element_text(size=rel(2), hjust = 0.5))+
  geom_hline(yintercept = 1)

ggplot() + geom_histogram(aes(x=unlist(pit_sim_b), y=..density..),breaks=hist(unlist(pit_sim_b), plot=F, breaks = 10)$breaks, col="black")+ ylim(0,2)+
  xlab("Proability Integral Transform")+ ylab("Density")+
  ggtitle(expression(paste("CMIK"[G], " (skewed)")))+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.8)),
        plot.title = element_text(size=rel(2), hjust = 0.5))+
  geom_hline(yintercept = 1)

ggplot() + geom_histogram(aes(x=unlist(pit_sim_t_b), y=..density..), breaks=hist(unlist(pit_sim_t_b), plot=F, breaks = 10)$breaks,  col="black")+ ylim(0,2)+
  xlab("Proability Integral Transform")+ ylab("Density")+
  ggtitle(expression(paste("CMIK"[t], " (skewed)")))+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.8)),
        plot.title = element_text(size=rel(2), hjust = 0.5))+
  geom_hline(yintercept = 1)




## heavy tailed
ggplot() + geom_histogram(aes(x=unlist(pit_g_c), y=..density..),
                          breaks=hist(unlist(pit_g_c), plot=F, breaks = 10)$breaks, col="black")+ ylim(0,2)+
  xlab("Proability Integral Transform")+ ylab("Density")+
  ggtitle("GK (heavy-tailed)")+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.8)),
        plot.title = element_text(size=rel(2), hjust = 0.5))+
  geom_hline(yintercept = 1)

ggplot() + geom_histogram(aes(x=unlist(pit_sim_var_c), y=..density..),
                          breaks=hist(unlist(pit_sim_var_c), plot=F, breaks = 10)$breaks, col="black")+ ylim(0,2)+
  xlab("Proability Integral Transform")+ ylab("Density")+
  ggtitle("VMIK (heavy-tailed)")+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.8)),
        plot.title = element_text(size=rel(2), hjust = 0.5))+
  geom_hline(yintercept = 1)

ggplot() + geom_histogram(aes(x=unlist(pit_sim_c), y=..density..),
                          breaks=hist(unlist(pit_sim_c), plot=F, breaks = 10)$breaks, col="black")+ ylim(0,2)+
  xlab("Proability Integral Transform")+ ylab("Density")+
  ggtitle(expression(paste("CMIK"[G], " (heavy-tailed)")))+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.8)),
        plot.title = element_text(size=rel(2), hjust = 0.5))+
  geom_hline(yintercept = 1)

ggplot() + geom_histogram(aes(x=unlist(pit_sim_t_c), y=..density..),
                          breaks=hist(unlist(pit_sim_t_c), plot=F, breaks = 10)$breaks,  col="black")+ ylim(0,2)+
  xlab("Proability Integral Transform")+ ylab("Density")+
  ggtitle(expression(paste("CMIK"[t], " (heavy-tailed)")))+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.8)),
        plot.title = element_text(size=rel(2), hjust = 0.5))+
  geom_hline(yintercept = 1)


## skewed and heavy tailed
ggplot() + geom_histogram(aes(x=unlist(pit_g_a), y=..density..),
                          breaks=hist(unlist(pit_g_a), plot=F, breaks = 10)$breaks, col="black")+ ylim(0,2)+
  xlab("Proability Integral Transform")+ ylab("Density")+
  ggtitle("GK (skewed & heavy-tailed)")+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.8)),
        plot.title = element_text(size=rel(2), hjust = 0.5))+
  geom_hline(yintercept = 1)

ggplot() + geom_histogram(aes(x=unlist(pit_sim_var_a), y=..density..),
                          breaks=hist(unlist(pit_sim_var_a), plot=F, breaks = 10)$breaks, col="black")+ ylim(0,2)+
  xlab("Proability Integral Transform")+ ylab("Density")+
  ggtitle("VMIK (skewed & heavy-tailed)")+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.8)),
        plot.title = element_text(size=rel(2), hjust = 0.5))+
  geom_hline(yintercept = 1)

ggplot() + geom_histogram(aes(x=unlist(pit_sim_a), y=..density..),
                          breaks=hist(unlist(pit_sim_a), plot=F, breaks = 10)$breaks, col="black")+ ylim(0,2)+
  xlab("Proability Integral Transform")+ ylab("Density")+
  ggtitle(expression(paste("CMIK"[G], " (skewed & heavy-tailed)")))+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.8)),
        plot.title = element_text(size=rel(2), hjust = 0.5))+
  geom_hline(yintercept = 1)

ggplot() + geom_histogram(aes(x=unlist(pit_sim_t_a), y=..density..), 
                          breaks=hist(unlist(pit_sim_t_a), plot=F, breaks = 10)$breaks,  col="black")+ ylim(0,2)+
  xlab("Proability Integral Transform")+ ylab("Density")+
  ggtitle(expression(paste("CMIK"[t], " (skewed & heavy-tailed)")))+
  theme(axis.text = element_text(size=rel(1.2)),
        axis.title = element_text(size=rel(1.8)),
        plot.title = element_text(size=rel(2), hjust = 0.5))+
  geom_hline(yintercept = 1)






# 2) CRPS

boxplot(crps_g_b,crpsF_y_var_b,crpsF_y_b, crpsF_y_t_b , main="CRPS (Skewed)",
        names =c("Gaussian prediction", "Variogram", "Gaussian copula", "t-copula" ))

boxplot(crps_g_c,crpsF_y_var_c,crpsF_y_c, crpsF_y_t_c , main="CRPS (heavy tailed)",
        names =c("Gaussian prediction", "Variogram", "Gaussian copula", "t-copula" ))

boxplot(crps_g_a,crpsF_y_var_a,crpsF_y_a, crpsF_y_t_a , main="CRPS (Skewed & heavy tailed)",
        names =c("Gaussian prediction", "Variogram", "Gaussian copula", "t-copula" ))



## Numerical assessment

## CRPS
data.frame(Skewed= c(mean(apply(crps_g_b,2,mean)),mean(apply(crpsF_y_var_b,2, mean)),mean(apply(crpsF_y_b,2,mean)), mean(apply(crpsF_y_t_b,2,mean))),
           ht = c(mean(apply(crps_g_c,2,mean)),mean(apply(crpsF_y_var_c,2,mean)),mean(apply(crpsF_y_c,2,mean)), mean(apply(crpsF_y_t_c,2,mean))),
           s_ht= c(mean(apply(crps_g_a,2,mean)),mean(apply(crpsF_y_var_a,2,mean)),mean(apply(crpsF_y_a,2,mean)), mean(apply(crpsF_y_t_a,2,mean)))
           )

##standard error
data.frame(Skewed= c(mean(apply(crps_g_b,2,sd)/sqrt(sim_size)),mean(apply(crpsF_y_var_b,2, sd)/sqrt(sim_size)),mean(apply(crpsF_y_b,2,sd)/sqrt(sim_size)), mean(apply(crpsF_y_t_b,2,sd)/sqrt(sim_size))),
           ht = c(mean(apply(crps_g_c,2,sd)/sqrt(sim_size)),mean(apply(crpsF_y_var_c,2,sd)/sqrt(sim_size)),mean(apply(crpsF_y_c,2,sd)/sqrt(sim_size)), mean(apply(crpsF_y_t_c,2,sd)/sqrt(sim_size))),
           s_ht= c(mean(apply(crps_g_a,2,sd)/sqrt(sim_size)),mean(apply(crpsF_y_var_a,2,sd)/sqrt(sim_size)),mean(apply(crpsF_y_a,2,sd)/sqrt(sim_size)), mean(apply(crpsF_y_t_a,2,sd)/sqrt(sim_size)))
)

##mCRPS
data.frame(Skewed= c(mean(apply(crps_g_b,2,median)),mean(apply(crpsF_y_var_b,2, median)),mean(apply(crpsF_y_b,2,median)), mean(apply(crpsF_y_t_b,2,median))),
           ht = c(mean(apply(crps_g_c,2,median)),mean(apply(crpsF_y_var_c,2,median)),mean(apply(crpsF_y_c,2,median)), mean(apply(crpsF_y_t_c,2,median))),
           s_ht= c(mean(apply(crps_g_a,2,median)),mean(apply(crpsF_y_var_a,2,median)),mean(apply(crpsF_y_a,2,median)), mean(apply(crpsF_y_t_a,2,median)))
)

## average length ci and coverage ###
#length
data.frame(Skewed= c(mean(length_ci_g_b),mean(length_ci_var_b),mean(length_ci_b),mean(length_ci_t_b)),
           ht = c(mean(length_ci_g_c),mean(length_ci_var_c),mean(length_ci_c),mean(length_ci_t_c)),
           s_ht= c(mean(length_ci_g_a),mean(length_ci_var_a),mean(length_ci_a),mean(length_ci_t_a))
)
# coverage
data.frame(Skewed= c(mean(cp_g_b),mean(cp_var_b),mean(cp_b),mean(cp_t_b)),
           ht = c(mean(cp_g_c),mean(cp_var_c),mean(cp_c),mean(cp_t_c)),
           s_ht= c(mean(cp_g_a),mean(cp_var_a),mean(cp_a),mean(cp_t_a))
)



