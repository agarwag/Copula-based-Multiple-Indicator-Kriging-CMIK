
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
sim_size = 500

g_values = c(0, 0,   0,   0.5, 0.75, 0.5, 0.3 )
h_values = c(0, 0.3, 0.5, 0,   0,    0.3, 0.4)
train_data_list =  vector(mode = "list", length = length(g_values))
test_data_list =  vector(mode = "list", length = length(g_values))
data_list = vector(mode = "list", length = length(g_values))

n=100
dt = sort(sample(n, size = n*.8, replace = F))

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




##################### Gaussian prediction for Tukey g and h process ####################################
## assuming that the distribution is gaussian
library(fields)
library(geoR)
library(scoringRules)


gaussian_pred = list(); gaussian_var = list(); tv = list(); parallel_g = list()
for (c in 1:length(g_values)) {
  
  for(i in 1:sim_size){
    
    train_data = train_data_list[[c]][[i]]
    test_data = test_data_list[[c]][[i]]
    
    train_loc = train_data[,1:2]
    test_loc = test_data[,1:2]
    
    if(c==1 | c==4 | c==5){
      ml <- likfit(data=train_data$tgh,coords=train_loc,
                   fix.nugget=F,cov.model="exponential",
                   ini.cov.pars  = c(0.98, 0.23), nugget = 0.01)
    }
    else if(c==2 | c==3 |c==6|c==7){
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
    
    gaussian_pred[[i]] = list(tv = test_data$tgh, predicted = pred$predict, var = pred$krige.var)
    
  }
  parallel_g[[c]] = gaussian_pred
}





### Gaussian marginal distribution (g=0.4,h=0.2)

rmse_g= matrix(0, nrow = sim_size, ncol = 20); mad_g = matrix(0, nrow = sim_size, ncol = 20)
rmse_g_c = matrix(0,nrow = length(g_values), ncol = 20); mad_g_c = matrix(0,nrow = length(g_values), ncol = 20)
sd_g_c = matrix(0,nrow = length(g_values), ncol = 20)

for (c in 1:length(g_values)) {
  
  for (i in 1:sim_size) {
    tv = parallel_g[[c]][[i]]$tv
    predicted = parallel_g[[c]][[i]]$predicted
    var = parallel_g[[c]][[i]]$var
    
    # rmse_g[i,] = sqrt(mean((tv-predicted)^2)) ## RMSE
    mad_g[i,] = abs(tv-predicted) ## MAD
  }
  #rmse_g_c[c] = mean(rmse_g);
  sd_g_c[c,] = 1.2533*apply(mad_g, 2, sd)/sqrt(sim_size)
  mad_g_c[c,] = apply(mad_g, 2, median)
}

rmse_g_c
mad_g_c




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
  
  if(c==1){ qseq = seq(0.05,0.95, length.out = 15) ## for gaussian
  } else if(c==2| c==3) {qseq = c(0.025,seq(0.05,0.95, length.out = 13),0.975)  ## for heavy tailed
  }  else if(c==4){ qseq = c(seq(0.05,0.75, length.out = 10) , seq(0.8,0.975, length.out = 5)) 
  }  else if(c==5){ qseq = c(seq(0.05,0.735, length.out = 10) , seq(0.8,0.975, length.out = 5)) ## for skewed
  } else if(c==6|c==7){qseq = c(seq(0.025,0.7, length.out = 10) , seq(0.8,0.975, length.out = 5))}
  
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
smoothF_gh = list()

rmse_sim= matrix(0, nrow = sim_size, ncol = 20); mad_sim = matrix(0, nrow = sim_size, ncol = 20)
rmse_c = matrix(0,nrow = length(g_values), ncol = 20); mad_c = matrix(0,nrow = length(g_values), ncol = 20)
sd_c = matrix(0,nrow = length(g_values), ncol = 20)

######### Measures for all marginal distributions
for(c in 1:length(g_values)){
  
for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim[[c]][[i]]$Fsim, function(a) smoothglobal(Fx=a, globalF = parallel_sim[[c]][[i]]$data_tgh) )
  median_pred_sim = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
  
  quantile_pred_sim_l = matrix(0, nrow=L, ncol = length(median_pred_sim))
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_sim_l[l,] = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  
  mean_pred_sim = apply(quantile_pred_sim_l, 2, mean)
  tv_gh = parallel_sim[[c]][[i]]$tv
 # rmse_sim[i,] = sqrt(mean((tv_gh-median_pred_sim)^2)) ## RMSE
  mad_sim[i,] = abs(tv_gh-median_pred_sim) ## MAD
  
}
#rmse_c[c,] = mean(rmse_sim)
sd_c[c,] = 1.2533*apply(mad_sim, 2, sd)/sqrt(sim_size)
mad_c[c,] = apply(mad_sim, 2, median)

}
rmse_c
mad_c



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
  #u = pnorm(Y) nu= par[2] #
  nu= 4
  logl <- dmvt(qt(u, df = nu), df=nu, sigma=R.m,log = TRUE) - sum(dt(qt(u, df=nu), df=nu, log=TRUE)) # log-density for X2(s) (denominator)
  
  #logl = log(dCopula(u, myCop))
  return(-logl)
}



library(doParallel)
registerDoParallel(cores=50)

start.time= Sys.time()
parallel_sim_t = list(); dat_gh= list()
for (c in 1:length(g_values)) {
  
  if(c==1){ qseq = seq(0.05,0.95, length.out = 15) ## for gaussian
  } else if(c==2| c==3) {qseq = c(0.025,seq(0.05,0.95, length.out = 13),0.975)  ## for heavy tailed
  }  else if(c==4){ qseq = c(seq(0.05,0.75, length.out = 10) , seq(0.8,0.975, length.out = 5)) 
  }  else if(c==5){ qseq = c(seq(0.05,0.735, length.out = 10) , seq(0.8,0.975, length.out = 5)) ## for skewed
  } else if(c==6|c==7){qseq = c(seq(0.025,0.7, length.out = 10) , seq(0.8,0.975, length.out = 5))}
  
  parallel_sim_t[[c]] = foreach(s=1:sim_size)  %dopar% {
    
    
    dat = data_list[[c]][[s]]
    train_data = train_data_list[[c]][[s]]
    test_data = test_data_list[[c]][[s]]
    
    Y = train_data$tgh
    coord = train_data[,1:2]
    
    coordinates(train_data) = ~x+y
    
    theta.hat_vec <- optim(par = c(0.25),logL.tcopula.Exp,Y=Y,coord=coord, method="L-BFGS-B",
                           lower=0.01,upper=10^5)$par
    
    theta.hat = theta.hat_vec[1]
    nu = 4 #round(theta.hat_vec[2],0)  #4
    
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



#### rmse/mad

L=200; 
smoothF_gh = list()

rmse_t= matrix(0, nrow = sim_size, ncol = 20); mad_t = matrix(0, nrow = sim_size, ncol = 20)
rmse_t_c = matrix(0,nrow = length(g_values), ncol = 20); mad_t_c = matrix(0,nrow = length(g_values), ncol = 20)
sd_t_c = matrix(0,nrow = length(g_values), ncol = 20)

for (c in 1:length(g_values)) {
  
for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim_t[[c]][[i]]$Fsim, function(a) smoothglobal(Fx=a, globalF = parallel_sim_t[[c]][[i]]$data_tgh) )
  median_pred_sim_t = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
  
  quantile_pred_sim_l = matrix(0, nrow=L, ncol = length(median_pred_sim_t))
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_sim_l[l,] = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  
  mean_pred_sim_t = apply(quantile_pred_sim_l, 2, mean)
  tv_gh = parallel_sim_t[[c]][[i]]$tv
 # rmse_t[i] = sqrt(mean((tv_gh-median_pred_sim_t)^2)) ## RMSE
  mad_t[i,] = abs(tv_gh-median_pred_sim_t) ## MAD
  
}
#rmse_t_c[c] = mean(rmse_t)
sd_t_c[c,] = 1.2533*apply(mad_t, 2, sd)/sqrt(sim_size)
mad_t_c[c,] = apply(mad_t, 2, median)
}
rmse_t_c
mad_t_c


###################################### Variogram ###################################
Fsimv_gh_list = list(); tv_gh_list = list(); data_tgh_list= list()
parallel_sim_var= list()

library(gstat)
for (k in 1:length(g_values)) {
  
  if(c==1){ qseq = seq(0.05,0.95, length.out = 15) ## for gaussian
  } else if(c==2| c==3) {qseq = c(0.025,seq(0.05,0.95, length.out = 13),0.975)  ## for heavy tailed
  }  else if(c==4){ qseq = c(seq(0.05,0.75, length.out = 10) , seq(0.8,0.975, length.out = 5)) 
  }  else if(c==5){ qseq = c(seq(0.05,0.735, length.out = 10) , seq(0.8,0.975, length.out = 5)) ## for skewed
  } else if(c==6|c==7){qseq = c(seq(0.025,0.7, length.out = 10) , seq(0.8,0.975, length.out = 5))}
  
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


####################### Marginal distribution: Skew and ht (g=0.4, h=0.2) #################################
L=200; 

rmse_var= matrix(0, nrow = sim_size, ncol = 20); mad_var = matrix(0, nrow = sim_size, ncol = 20)
rmse_var_c = matrix(0,nrow = length(g_values), ncol = 20); mad_var_c = matrix(0,nrow = length(g_values), ncol = 20)
sd_var_c = matrix(0,nrow = length(g_values), ncol = 20)

for (c in 1:length(g_values)) {
  
for (i in 1:sim_size) {
  smoothF_gh[[i]] =   lapply (parallel_sim_var[[c]]$Fsim[[i]], 
                              function(a) smoothglobal(Fx=a, globalF = parallel_sim_var[[c]]$data_tgh[[i]]) )
  median_pred_sim_var = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
  
  quantile_pred_sim_l = matrix(0, nrow=L, ncol = length(median_pred_sim_var))
  for (l in 1:L) {
    pl=l/(L+1)
    quantile_pred_sim_l[l,] = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
  }
  
  mean_pred_sim_var = apply(quantile_pred_sim_l, 2, mean)
  tv_gh = parallel_sim_var[[c]]$tv[[i]]
  #rmse_var[i] = sqrt(mean((tv_gh-median_pred_sim_var)^2)) ## RMSE
  mad_var[i,] = abs(tv_gh-median_pred_sim_var) ## MAD
  
}
#rmse_var_c[c] = mean(rmse_var)
sd_var_c[c,] = 1.2533*apply(mad_var, 2, sd)/sqrt(sim_size)
mad_var_c[c,] = apply(mad_var, 2, median)
}

rmse_var_c
mad_var_c

cbind(g_values, h_values,rmse_g_c,
rmse_var_c,
rmse_c,
rmse_t_c)


cbind(g_values, h_values, apply(mad_g_c, 1, mean),
apply(mad_var_c, 1, mean),
apply(mad_c, 1, mean),
apply(mad_t_c, 1, mean))

































###### theta  = 0.5

### Simulate data from Tukey g and h

theta = 0.375

g_values = c(0, 0,   0,   0.5, 0.75, 0.5, 0.3 )
h_values = c(0, 0.3, 0.5, 0,   0,    0.3, 0.4)
train_data_list =  vector(mode = "list", length = length(g_values))
test_data_list =  vector(mode = "list", length = length(g_values))
data_list = vector(mode = "list", length = length(g_values))

for (c in 1:length(g_values)) {
  
  for(i in 1:sim_size){
    
    g= g_values[c]
    h= h_values[c]
    
    dat = dat2 = Dtgh1(nx=10, theta=theta, g, h)
    dat = dat2[[1]]
    
    dt = sort(sample(nrow(dat), nrow(dat)*.8, replace = F))
    
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
  
  if(c==1){ qseq = seq(0.05,0.95, length.out = 15) ## for gaussian
  } else if(c==2| c==3) {qseq = c(0.025,seq(0.05,0.95, length.out = 13),0.975)  ## for heavy tailed
  }  else if(c==4){ qseq = c(seq(0.05,0.75, length.out = 10) , seq(0.8,0.975, length.out = 5)) 
  }  else if(c==5){ qseq = c(seq(0.05,0.735, length.out = 10) , seq(0.8,0.975, length.out = 5)) ## for skewed
  } else if(c==6|c==7){qseq = c(seq(0.025,0.7, length.out = 10) , seq(0.8,0.975, length.out = 5))}
  
  parallel_sim[[c]] = foreach(s=1:sim_size)  %dopar% {
    
    
    dat = data_list[[c]][[s]]
    train_data = train_data_list[[c]][[s]]
    test_data = test_data_list[[c]][[s]]
    
    Y = train_data$tgh
    coord = train_data[,1:2]
    
    coordinates(train_data) = ~x+y
    
    theta.hat <- optim(par=0.5,logL.Gcopula.Exp,Y=Y,coord=coord, method="L-BFGS-B",
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

smoothF_gh= list()
rmse_sim= matrix(0, nrow = sim_size, ncol = 20); mad_sim = matrix(0, nrow = sim_size, ncol = 20)
rmse_c0.5 = matrix(0,nrow = length(g_values), ncol = 20); mad_c0.5 = matrix(0,nrow = length(g_values), ncol = 20)
sd_c0.5 = matrix(0,nrow = length(g_values), ncol = 20)

######### Measures for all marginal distributions
for(c in 1:length(g_values)){
  
  for (i in 1:sim_size) {
    smoothF_gh[[i]] =   lapply (parallel_sim[[c]][[i]]$Fsim, function(a) smoothglobal(Fx=a, globalF = parallel_sim[[c]][[i]]$data_tgh) )
    median_pred_sim = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
    
    quantile_pred_sim_l = matrix(0, nrow=L, ncol = length(median_pred_sim))
    for (l in 1:L) {
      pl=l/(L+1)
      quantile_pred_sim_l[l,] = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
    }
    
    mean_pred_sim = apply(quantile_pred_sim_l, 2, mean)
    tv_gh = parallel_sim[[c]][[i]]$tv
    # rmse_sim[i,] = sqrt(mean((tv_gh-median_pred_sim)^2)) ## RMSE
    mad_sim[i,] = abs(tv_gh-median_pred_sim) ## MAD
    
  }
  #rmse_c0.5[c,] = mean(rmse_sim)
  sd_c0.5[c,] = 1.2533*apply(mad_sim, 2, sd)/sqrt(sim_size)
  mad_c0.5[c,] = apply(mad_sim, 2, median)
  
}
rmse_c0.5
mad_c0.5


##################### Gaussian prediction for Tukey g and h process ####################################
## assuming that the distribution is gaussian
library(fields)
library(geoR)
library(scoringRules)


gaussian_pred = list(); gaussian_var = list(); tv = list(); parallel_g = list()
for (c in 1:length(g_values)) {
  
  for(i in 1:sim_size){
    
    train_data = train_data_list[[c]][[i]]
    test_data = test_data_list[[c]][[i]]
    
    train_loc = train_data[,1:2]
    test_loc = test_data[,1:2]
    
    if(c==1 | c==4 | c==5){
      ml <- likfit(data=train_data$tgh,coords=train_loc,
                   fix.nugget=F,cov.model="exponential",
                   ini.cov.pars  = c(0.98, 0.23), nugget = 0.01)
    }
    else if(c==2 | c==3 |c==6|c==7){
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
    
    gaussian_pred[[i]] = list(tv = test_data$tgh, predicted = pred$predict, var = pred$krige.var)
    
  }
  parallel_g[[c]] = gaussian_pred
}





### Gaussian marginal distribution (g=0.4,h=0.2)

rmse_g= matrix(0, nrow = sim_size, ncol = 20); mad_g = matrix(0, nrow = sim_size, ncol = 20)
rmse_g_c0.5 = matrix(0,nrow = length(g_values), ncol = 20); mad_g_c0.5 = matrix(0,nrow = length(g_values), ncol = 20)
sd_g_c0.5 = matrix(0,nrow = length(g_values), ncol = 20)

for (c in 1:length(g_values)) {
  
  for (i in 1:sim_size) {
    tv = parallel_g[[c]][[i]]$tv
    predicted = parallel_g[[c]][[i]]$predicted
    var = parallel_g[[c]][[i]]$var
    
    # rmse_g[i,] = sqrt(mean((tv-predicted)^2)) ## RMSE
    mad_g[i,] = abs(tv-predicted) ## MAD
  }
  #rmse_g_c0.5[c] = mean(rmse_g);
  sd_g_c0.5[c,] = 1.2533*apply(mad_g, 2, sd)/sqrt(sim_size)
  mad_g_c0.5[c,] = apply(mad_g, 2, median)
}

rmse_g_c0.5
mad_g_c0.5



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
  #u = pnorm(Y) nu= par[2] #
  nu= 4
  logl <- dmvt(qt(u, df = nu), df=nu, sigma=R.m,log = TRUE) - sum(dt(qt(u, df=nu), df=nu, log=TRUE)) # log-density for X2(s) (denominator)
  
  #logl = log(dCopula(u, myCop))
  return(-logl)
}



library(doParallel)
registerDoParallel(cores=50)

start.time= Sys.time()
parallel_sim_t = list(); dat_gh= list()
for (c in 1:length(g_values)) {
  
  if(c==1){ qseq = seq(0.05,0.95, length.out = 15) ## for gaussian
  } else if(c==2| c==3) {qseq = c(0.025,seq(0.05,0.95, length.out = 13),0.975)  ## for heavy tailed
  }  else if(c==4){ qseq = c(seq(0.05,0.75, length.out = 10) , seq(0.8,0.975, length.out = 5)) 
  }  else if(c==5){ qseq = c(seq(0.05,0.735, length.out = 10) , seq(0.8,0.975, length.out = 5)) ## for skewed
  } else if(c==6|c==7){qseq = c(seq(0.025,0.7, length.out = 10) , seq(0.8,0.975, length.out = 5))}
  
  parallel_sim_t[[c]] = foreach(s=1:sim_size)  %dopar% {
    
    
    dat = data_list[[c]][[s]]
    train_data = train_data_list[[c]][[s]]
    test_data = test_data_list[[c]][[s]]
    
    Y = train_data$tgh
    coord = train_data[,1:2]
    
    coordinates(train_data) = ~x+y
    
    theta.hat_vec <- optim(par = c(0.5),logL.tcopula.Exp,Y=Y,coord=coord, method="L-BFGS-B",
                           lower=0.01,upper=10^5)$par
    
    theta.hat = theta.hat_vec[1]
    nu = 4 #round(theta.hat_vec[2],0)  #4
    
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



#### rmse/mad

L=200; 
smoothF_gh = list()

rmse_t= matrix(0, nrow = sim_size, ncol = 20); mad_t = matrix(0, nrow = sim_size, ncol = 20)
rmse_t_c0.5 = matrix(0,nrow = length(g_values), ncol = 20); mad_t_c0.5 = matrix(0,nrow = length(g_values), ncol = 20)
sd_t_c0.5 = matrix(0,nrow = length(g_values), ncol = 20)

for (c in 1:length(g_values)) {
  
  for (i in 1:sim_size) {
    smoothF_gh[[i]] =   lapply (parallel_sim_t[[c]][[i]]$Fsim, function(a) smoothglobal(Fx=a, globalF = parallel_sim_t[[c]][[i]]$data_tgh) )
    median_pred_sim_t = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
    
    quantile_pred_sim_l = matrix(0, nrow=L, ncol = length(median_pred_sim_t))
    for (l in 1:L) {
      pl=l/(L+1)
      quantile_pred_sim_l[l,] = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
    }
    
    mean_pred_sim_t = apply(quantile_pred_sim_l, 2, mean)
    tv_gh = parallel_sim_t[[c]][[i]]$tv
    # rmse_t[i] = sqrt(mean((tv_gh-median_pred_sim_t)^2)) ## RMSE
    mad_t[i,] = abs(tv_gh-median_pred_sim_t) ## MAD
    
  }
  #rmse_t_c0.5[c] = mean(rmse_t)
  sd_t_c0.5[c,] = 1.2533*apply(mad_t, 2, sd)/sqrt(sim_size)
  mad_t_c0.5[c,] = apply(mad_t, 2, median)
}
rmse_t_c0.5
mad_t_c0.5



###################################### Variogram ###################################
Fsimv_gh_list = list(); tv_gh_list = list(); data_tgh_list= list()
parallel_sim_var= list()

library(gstat)
for (k in 1:length(g_values)) {
  
  if(c==1){ qseq = seq(0.05,0.95, length.out = 15) ## for gaussian
  } else if(c==2| c==3) {qseq = c(0.025,seq(0.05,0.95, length.out = 13),0.975)  ## for heavy tailed
  }  else if(c==4){ qseq = c(seq(0.05,0.75, length.out = 10) , seq(0.8,0.975, length.out = 5)) 
  }  else if(c==5){ qseq = c(seq(0.05,0.735, length.out = 10) , seq(0.8,0.975, length.out = 5)) ## for skewed
  } else if(c==6|c==7){qseq = c(seq(0.025,0.7, length.out = 10) , seq(0.8,0.975, length.out = 5))}
  
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


####################### Marginal distribution: Skew and ht (g=0.4, h=0.2) #################################
L=200; 

rmse_var= matrix(0, nrow = sim_size, ncol = 20); mad_var = matrix(0, nrow = sim_size, ncol = 20)
rmse_var_c0.5 = matrix(0,nrow = length(g_values), ncol = 20); mad_var_c0.5 = matrix(0,nrow = length(g_values), ncol = 20)
sd_var_c0.5 = matrix(0,nrow = length(g_values), ncol = 20)

for (c in 1:length(g_values)) {
  
  for (i in 1:sim_size) {
    smoothF_gh[[i]] =   lapply (parallel_sim_var[[c]]$Fsim[[i]], 
                                function(a) smoothglobal(Fx=a, globalF = parallel_sim_var[[c]]$data_tgh[[i]]) )
    median_pred_sim_var = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
    
    quantile_pred_sim_l = matrix(0, nrow=L, ncol = length(median_pred_sim_var))
    for (l in 1:L) {
      pl=l/(L+1)
      quantile_pred_sim_l[l,] = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
    }
    
    mean_pred_sim_var = apply(quantile_pred_sim_l, 2, mean)
    tv_gh = parallel_sim_var[[c]]$tv[[i]]
    #rmse_var[i] = sqrt(mean((tv_gh-median_pred_sim_var)^2)) ## RMSE
    mad_var[i,] = abs(tv_gh-median_pred_sim_var) ## MAD
    
  }
  #rmse_var_c0.5[c] = mean(rmse_var)
  sd_var_c0.5[c,] = 1.2533*apply(mad_var, 2, sd)/sqrt(sim_size)
  mad_var_c0.5[c,] = apply(mad_var, 2, median)
}

rmse_var_c0.5
mad_var_c0.5

cbind(g_values, h_values,rmse_g_c0.5,
      rmse_var_c0.5,
      rmse_c0.5,
      rmse_t_c0.5)


cbind(g_values, h_values,mad_g_c0.5,
      mad_var_c0.5,
      mad_c0.5,
      mad_t_c0.5)



























###### theta  = 0.75

### Simulate data from Tukey g and h
theta = 0.5

g_values = c(0, 0,   0,   0.5, 0.75, 0.5, 0.3 )
h_values = c(0, 0.3, 0.5, 0,   0,    0.3, 0.4)
train_data_list =  vector(mode = "list", length = length(g_values))
test_data_list =  vector(mode = "list", length = length(g_values))
data_list = vector(mode = "list", length = length(g_values))

for (c in 1:length(g_values)) {
  
  for(i in 1:sim_size){
    
    g= g_values[c]
    h= h_values[c]
    
    dat = dat2 = Dtgh1(nx=10, theta=theta, g, h)
    dat = dat2[[1]]
    
    dt = sort(sample(nrow(dat), nrow(dat)*.8, replace = F))
    
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
  
  if(c==1){ qseq = seq(0.05,0.95, length.out = 15) ## for gaussian
  } else if(c==2| c==3) {qseq = c(0.025,seq(0.05,0.95, length.out = 13),0.975)  ## for heavy tailed
  }  else if(c==4){ qseq = c(seq(0.05,0.75, length.out = 10) , seq(0.8,0.975, length.out = 5)) 
  }  else if(c==5){ qseq = c(seq(0.05,0.735, length.out = 10) , seq(0.8,0.975, length.out = 5)) ## for skewed
  } else if(c==6|c==7){qseq = c(seq(0.025,0.7, length.out = 10) , seq(0.8,0.975, length.out = 5))}
  
  parallel_sim[[c]] = foreach(s=1:sim_size)  %dopar% {
    
    
    dat = data_list[[c]][[s]]
    train_data = train_data_list[[c]][[s]]
    test_data = test_data_list[[c]][[s]]
    
    Y = train_data$tgh
    coord = train_data[,1:2]
    
    coordinates(train_data) = ~x+y
    
    theta.hat <- optim(par=0.75,logL.Gcopula.Exp,Y=Y,coord=coord, method="L-BFGS-B",
                       lower=0.01,upper=1)$par
    
    ################################# Mutiple thresholds ###########################################
    th = quantile(Y, qseq) 
    
    plot(ecdf(dat$tgh), main="Global ECDF of RR conc.")
    abline(v=th, lty=2)
    
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

smoothF_gh= list()
rmse_sim= matrix(0, nrow = sim_size, ncol = 20); mad_sim = matrix(0, nrow = sim_size, ncol = 20)
rmse_c0.75 = matrix(0,nrow = length(g_values), ncol = 20); mad_c0.75 = matrix(0,nrow = length(g_values), ncol = 20)
sd_c0.75 = matrix(0,nrow = length(g_values), ncol = 20)


######### Measures for all marginal distributions
for(c in 1:length(g_values)){
  
  for (i in 1:sim_size) {
    smoothF_gh[[i]] =   lapply (parallel_sim[[c]][[i]]$Fsim, function(a) smoothglobal(Fx=a, globalF = parallel_sim[[c]][[i]]$data_tgh) )
    median_pred_sim = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
    
    quantile_pred_sim_l = matrix(0, nrow=L, ncol = length(median_pred_sim))
    for (l in 1:L) {
      pl=l/(L+1)
      quantile_pred_sim_l[l,] = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
    }
    
    mean_pred_sim = apply(quantile_pred_sim_l, 2, mean)
    tv_gh = parallel_sim[[c]][[i]]$tv
    # rmse_sim[i,] = sqrt(mean((tv_gh-median_pred_sim)^2)) ## RMSE
    mad_sim[i,] = abs(tv_gh-median_pred_sim) ## MAD
    
  }
  #rmse_c0.5[c,] = mean(rmse_sim)
  sd_c0.75[c,] = 1.2533*apply(mad_sim, 2, sd)/sqrt(sim_size)
  mad_c0.75[c,] = apply(mad_sim, 2, median)
  
}
rmse_c0.75
mad_c0.75

##################### Gaussian prediction for Tukey g and h process ####################################
## assuming that the distribution is gaussian
library(fields)
library(geoR)
library(scoringRules)


gaussian_pred = list(); gaussian_var = list(); tv = list(); parallel_g = list()
for (c in 1:length(g_values)) {
  
  for(i in 1:sim_size){
    
    train_data = train_data_list[[c]][[i]]
    test_data = test_data_list[[c]][[i]]
    
    train_loc = train_data[,1:2]
    test_loc = test_data[,1:2]
    
    if(c==1 | c==4 | c==5){
      ml <- likfit(data=train_data$tgh,coords=train_loc,
                   fix.nugget=F,cov.model="exponential",
                   ini.cov.pars  = c(0.98, 0.23), nugget = 0.01)
    }
    else if(c==2 | c==3 |c==6|c==7){
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
    
    gaussian_pred[[i]] = list(tv = test_data$tgh, predicted = pred$predict, var = pred$krige.var)
    
  }
  parallel_g[[c]] = gaussian_pred
}





### Gaussian marginal distribution (g=0.4,h=0.2)

rmse_g= matrix(0, nrow = sim_size, ncol = 20); mad_g = matrix(0, nrow = sim_size, ncol = 20)
rmse_g_c0.75 = matrix(0,nrow = length(g_values), ncol = 20); mad_g_c0.75 = matrix(0,nrow = length(g_values), ncol = 20)
sd_g_c0.75 = matrix(0,nrow = length(g_values), ncol = 20)

for (c in 1:length(g_values)) {
  
  for (i in 1:sim_size) {
    tv = parallel_g[[c]][[i]]$tv
    predicted = parallel_g[[c]][[i]]$predicted
    var = parallel_g[[c]][[i]]$var
    
    # rmse_g[i,] = sqrt(mean((tv-predicted)^2)) ## RMSE
    mad_g[i,] = abs(tv-predicted) ## MAD
  }
  #rmse_g_c0.5[c] = mean(rmse_g);
  sd_g_c0.75[c,] = 1.2533*apply(mad_g, 2, sd)/sqrt(sim_size)
  mad_g_c0.75[c,] = apply(mad_g, 2, median)
}

rmse_g_c0.75
mad_g_c0.75



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
  #u = pnorm(Y) nu= par[2] #
  nu= 4
  logl <- dmvt(qt(u, df = nu), df=nu, sigma=R.m,log = TRUE) - sum(dt(qt(u, df=nu), df=nu, log=TRUE)) # log-density for X2(s) (denominator)
  
  #logl = log(dCopula(u, myCop))
  return(-logl)
}



library(doParallel)
registerDoParallel(cores=50)

start.time= Sys.time()
parallel_sim_t = list(); dat_gh= list()
for (c in 1:length(g_values)) {
  
  if(c==1){ qseq = seq(0.05,0.95, length.out = 15) ## for gaussian
  } else if(c==2| c==3) {qseq = c(0.025,seq(0.05,0.95, length.out = 13),0.975)  ## for heavy tailed
  }  else if(c==4){ qseq = c(seq(0.05,0.75, length.out = 10) , seq(0.8,0.975, length.out = 5)) 
  }  else if(c==5){ qseq = c(seq(0.05,0.735, length.out = 10) , seq(0.8,0.975, length.out = 5)) ## for skewed
  } else if(c==6|c==7){qseq = c(seq(0.025,0.7, length.out = 10) , seq(0.8,0.975, length.out = 5))}
  
  parallel_sim_t[[c]] = foreach(s=1:sim_size)  %dopar% {
    
    
    dat = data_list[[c]][[s]]
    train_data = train_data_list[[c]][[s]]
    test_data = test_data_list[[c]][[s]]
    
    Y = train_data$tgh
    coord = train_data[,1:2]
    
    coordinates(train_data) = ~x+y
    
    theta.hat_vec <- optim(par = c(0.75),logL.tcopula.Exp,Y=Y,coord=coord, method="L-BFGS-B",
                           lower=0.01,upper=10^5)$par
    
    theta.hat = 0.75#theta.hat_vec[1]
    nu = 4 #round(theta.hat_vec[2],0)  #4
    
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



#### rmse/mad

L=200; 
smoothF_gh = list()

rmse_t= matrix(0, nrow = sim_size, ncol = 20); mad_t = matrix(0, nrow = sim_size, ncol = 20)
rmse_t_c0.75 = matrix(0,nrow = length(g_values), ncol = 20); mad_t_c0.75 = matrix(0,nrow = length(g_values), ncol = 20)
sd_t_c0.75 = matrix(0,nrow = length(g_values), ncol = 20)

for (c in 1:length(g_values)) {
  
  for (i in 1:sim_size) {
    smoothF_gh[[i]] =   lapply (parallel_sim_t[[c]][[i]]$Fsim, function(a) smoothglobal(Fx=a, globalF = parallel_sim_t[[c]][[i]]$data_tgh) )
    median_pred_sim_t = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
    
    quantile_pred_sim_l = matrix(0, nrow=L, ncol = length(median_pred_sim_t))
    for (l in 1:L) {
      pl=l/(L+1)
      quantile_pred_sim_l[l,] = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
    }
    
    mean_pred_sim_t = apply(quantile_pred_sim_l, 2, mean)
    tv_gh = parallel_sim_t[[c]][[i]]$tv
    # rmse_t[i] = sqrt(mean((tv_gh-median_pred_sim_t)^2)) ## RMSE
    mad_t[i,] = abs(tv_gh-median_pred_sim_t) ## MAD
    
  }
  #rmse_t_c0.75[c] = mean(rmse_t)
  sd_t_c0.75[c,] = 1.2533*apply(mad_t, 2, sd)/sqrt(sim_size)
  mad_t_c0.75[c,] = apply(mad_t, 2, median)
}
rmse_t_c0.75
mad_t_c0.75


###################################### Variogram ###################################
Fsimv_gh_list = list(); tv_gh_list = list(); data_tgh_list= list()
parallel_sim_var= list()

library(gstat)
for (k in 1:length(g_values)) {
  
  if(c==1){ qseq = seq(0.05,0.95, length.out = 15) ## for gaussian
  } else if(c==2| c==3) {qseq = c(0.025,seq(0.05,0.95, length.out = 13),0.975)  ## for heavy tailed
  }  else if(c==4){ qseq = c(seq(0.05,0.75, length.out = 10) , seq(0.8,0.975, length.out = 5)) 
  }  else if(c==5){ qseq = c(seq(0.05,0.735, length.out = 10) , seq(0.8,0.975, length.out = 5)) ## for skewed
  } else if(c==6|c==7){qseq = c(seq(0.025,0.7, length.out = 10) , seq(0.8,0.975, length.out = 5))}
  
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


####################### Marginal distribution: Skew and ht (g=0.4, h=0.2) #################################
L=200; 

rmse_var= matrix(0, nrow = sim_size, ncol = 20); mad_var = matrix(0, nrow = sim_size, ncol = 20)
rmse_var_c0.75 = matrix(0,nrow = length(g_values), ncol = 20); mad_var_c0.75 = matrix(0,nrow = length(g_values), ncol = 20)
sd_var_c0.75 = matrix(0,nrow = length(g_values), ncol = 20)

for (c in 1:length(g_values)) {
  
  for (i in 1:sim_size) {
    smoothF_gh[[i]] =   lapply (parallel_sim_var[[c]]$Fsim[[i]], 
                                function(a) smoothglobal(Fx=a, globalF = parallel_sim_var[[c]]$data_tgh[[i]]) )
    median_pred_sim_var = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=0.5)])}))
    
    quantile_pred_sim_l = matrix(0, nrow=L, ncol = length(median_pred_sim_var))
    for (l in 1:L) {
      pl=l/(L+1)
      quantile_pred_sim_l[l,] = as.numeric(lapply(smoothF_gh[[i]], function(a) {min(a[,1][which(a[,2]>=pl)])}) )
    }
    
    mean_pred_sim_var = apply(quantile_pred_sim_l, 2, mean)
    tv_gh = parallel_sim_var[[c]]$tv[[i]]
    #rmse_var[i] = sqrt(mean((tv_gh-median_pred_sim_var)^2)) ## RMSE
    mad_var[i,] = abs(tv_gh-median_pred_sim_var) ## MAD
    
  }
  #rmse_var_c0.75[c] = mean(rmse_var)
  sd_var_c0.75[c,] = 1.2533*apply(mad_var, 2, sd)/sqrt(sim_size)
  mad_var_c0.75[c,] = apply(mad_var, 2, median)
}

#rmse_var_c0.75
#mad_var_c0.75

#cbind(g_values, h_values,rmse_g_c0.75,
#      rmse_var_c0.75,
#      rmse_c0.75,
#      rmse_t_c0.75)




### MAD
cbind(g_values, h_values, apply(mad_g_c,1, mean),
      apply(mad_var_c, 1, mean),
      apply(mad_c,1,mean),
      apply(mad_t_c,1,mean))


cbind(g_values, h_values, apply(mad_g_c0.5,1,mean),
      apply(mad_var_c0.5,1,mean),
      apply(mad_c0.5,1,mean),
      apply(mad_t_c0.5,1,mean))

cbind(g_values, h_values, apply(mad_g_c0.75,1,mean),
     apply(mad_var_c0.75,1,mean),
      apply(mad_c0.75,1,mean),
      apply(mad_t_c0.75,1,mean))


## standard error of MAD 
cbind(g_values, h_values, apply(sd_g_c,1, mean),
      apply(sd_var_c, 1, mean),
      apply(sd_c,1,mean),
      apply(sd_t_c,1,mean))

cbind(g_values, h_values, apply(sd_g_c0.5,1, mean),
      apply(sd_var_c0.5, 1, mean),
      apply(sd_c0.5,1,mean),
      apply(sd_t_c0.5,1,mean))

cbind(g_values, h_values, apply(sd_g_c0.75,1, mean),
      apply(sd_var_c0.75, 1, mean),
      apply(sd_c0.75,1,mean),
      apply(sd_t_c0.75,1,mean))


## sd of MAD over locations
cbind(g_values, h_values, apply(mad_g_c,1, sd),
      apply(mad_var_c, 1, sd),
      apply(mad_c,1,sd),
      apply(mad_t_c,1,sd))

cbind(g_values, h_values, apply(mad_g_c0.5,1, sd),
      apply(mad_var_c0.5, 1, sd),
      apply(mad_c0.5,1,sd),
      apply(mad_t_c0.5,1,sd))

cbind(g_values, h_values, apply(mad_g_c0.75,1, sd),
      apply(mad_var_c0.75, 1, sd),
      apply(mad_c0.75,1,sd),
      apply(mad_t_c0.75,1,sd))

## range
#cbind(g_values, h_values, c(diff(apply(mad_g_c,1, range))),
#      c(diff(apply(mad_var_c, 1, range))),
#      c(diff(apply(mad_c,1,range))),
#      c(diff(apply(mad_t_c,1,range))))

