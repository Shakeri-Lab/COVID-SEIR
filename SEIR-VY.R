library(foreach)
library(doParallel)
numCores <- commandArgs(trailingOnly=TRUE)[1]
numCores <- as.numeric(numCores) - 1
registerDoParallel(cores=numCores)
print(paste('number of cores is ',numCores, ' for YV'))

library(doRNG)
registerDoRNG(625904618)
library(tidyverse)
library(pomp)
options(stringsAsFactors=FALSE)
stopifnot(packageVersion("pomp")>="3.0")


run_level <- 3
covid_Np <-          switch(run_level,100, 1e3, 2e4)
covid_Nmif <-        switch(run_level, 10, 100, 100)
covid_Nreps_eval <-  switch(run_level,  2,  10,  10)
covid_Nreps_local <- switch(run_level, 10,  10,  10)
covid_Nreps_global <-switch(run_level, 10,  20,  20)
covid_Nsim <-        switch(run_level, 50, 100, 100) 


set.seed(1350254336)
setwd("/home/mf4yc/SEIRR")

cache_address = "/project/shakeri-lab/cache/cache_abm"


#########################################################################
#--------------------------|  rproc  |----------------------------------#
#########################################################################
rproc <- Csnippet("
                  double beta, foi, dw, births, mu_SE;
                  
                  //we only consider those that participate in the epidemic:
                  double pop = S + E + I + R;
                  
                  // transmission rate
                  beta = b0;
                  
                  // expected force of infection. iota: imported infections
                  // alpha mixing parameter, = 1:homogeneous mixing
                  foi = beta*pow(I+iota, alpha)/pop;
                  
                  // white noise (extrademographic stochasticity)
                  dw = rgammawn(sigmaSE,dt);
                  
                  mu_SE = foi*dw/dt;  // stochastic force of infection
                  
                  // Poisson births: fraction of leak into S from N
                  births = rpois(br*dt);
                  
                  
                  // State Updates:
                  double dN_SE  = rbinom(S , 1-exp(-mu_SE  *dt));
                  double dN_EI  = rbinom(E , 1-exp(-mu_EI  *dt));
                  double dN_IR  = rbinom(I , 1-exp(-mu_IR *dt));
                  S += births - dN_SE;
                  E += dN_SE  - dN_EI;
                  I += dN_EI  - dN_IR;
                  R += dN_IR;
                  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
                  ")

#########################################################################
#--------------------------|  rinit  |----------------------------------#
#########################################################################
rinit <- Csnippet("
                  double m = eta*N;
                  S = nearbyint(m*S_0);
                  E = nearbyint(m*E_0);
                  I = nearbyint(m*I_0);
                  R = nearbyint(m*R_0);
                  W = 0;
                  ")

#########################################################################
#--------------------------|  dmeas  |----------------------------------#
#########################################################################
dmeas <- Csnippet("
                  // Model for Viral Load
                  double shed_cases = E + I;
                  double mu_V = rho_V*shed_cases;
                  //double std_V = sqrt(mu_V*(1+od_V));
                  double lik_V = dnorm(V, mu_V, sd_V, 1);
                  
                  // Model for Case Counts
                  double mu_Y = rho_Y*I;
                  double std_Y = sqrt(mu_Y*(1+od_Y));
                  double lik_Y;
                  
                  if (Y > 0.0) {
                  lik_Y = pnorm(Y+0.5,mu_Y,std_Y,1,1)
                  - pnorm(Y-0.5,mu_Y,std_Y,1,1);
                  } else {
                  lik_Y = pnorm(Y+0.5,mu_Y,std_Y,1,1);
                  }
                  
                  // Combined likelihood
                  lik = lik_V + lik_Y;
                  //lik = lik_V;
                  //lik = lik_Y;
                  lik = (give_log) ? lik : exp(lik);
                  
                  ")

#########################################################################
#--------------------------|  rmeas  |----------------------------------#
#########################################################################
rmeas <- Csnippet("
                  // Viral Load
                  double shed_cases = E + I;
                  double mu_V = rho_V*shed_cases;
                  //double std_V = sqrt(mu_V*(1+od_V));
                  V = rnorm(mu_V, sd_V);
                  
                  // Case Counts
                  double mu_Y = rho_Y*I;
                  double std_Y = sqrt(mu_Y*(1+od_Y));
                  Y = rnorm(mu_Y, std_Y);
                  if (Y > 0.0) {
                  Y = nearbyint(Y);
                  } else {
                  Y = 0.0;
                  }
                  
                  ")

#########################################################################
#-------------------------|  Load Data  |-------------------------------#
#########################################################################
NewHaven = read_csv("/home/mf4yc/SEIRR/Data/abm.csv")

#########################################################################
#-------------------------|  Parameters  |------------------------------#
#########################################################################
parameters = c(
  "b0", "alpha", "iota",      
  "sigmaSE",                  
  "br",                       
  "mu_EI", "mu_IR",           
  "N",                        
  "eta",                      
  "rho_V", "sd_V",            
  "rho_Y", "od_Y",            
  "S_0","E_0","I_0", "R_0")


par_trans = parameter_trans(
  log = c(
    "b0", "alpha", "iota","sigmaSE", "br",
    "rho_V", "sd_V", "od_Y"),
  logit = c("mu_EI", "mu_IR","eta", "rho_Y"),
  barycentric=c("S_0","E_0","I_0", "R_0")
)
states = c("S", "E", "I", "R", "W")

#########################################################################
#-------------------------|  Covariates  |------------------------------#
#########################################################################

sdm_covar <- covariate_table(
  t=      NewHaven[["day"]],
  sdmm=   NewHaven[["sdm"]],
  event=  NewHaven[["events"]],
  order=  "constant",
  times=  "t"
)

# shifting case counts by the assumed reporting delay
rep_del = 5
NewHaven %>% mutate_at(c("Y_1"), 
                       tibble::lst("Y_1"=lead), 
                       n=rep_del) %>%
  mutate_at(c("Y_2"), 
            tibble::lst("Y_2"=lead), 
            n=rep_del)%>%
  mutate_at(c("Y_3"), 
            tibble::lst("Y_3"=lead), 
            n=rep_del)%>%
  mutate_at(c("Y"), 
            tibble::lst("Y"=lead), 
            n=rep_del)-> NewHaven_c

# focusing on the first peak for now
NewHaven <- NewHaven_c[1:70,]

#########################################################################
#-------------------------|  pomp Model  |------------------------------#
#########################################################################

covidSEIRsR = NewHaven %>%
  select(-logV) %>%
  rename(V = V,
         Y = Y
  ) %>%
  pomp(
    times = "day", # column name of data that corresponds to time
    t0 = 0,        # starting time
    # rprocess = discrete_time(rproc, delta.t=1), # daily
    rprocess = euler(rproc, delta.t=1/6), # every four
    rinit = rinit,
    rmeasure = rmeas,
    dmeasure = dmeas,
    accumvars= c("W"),
    partrans = par_trans,
    statenames = states,
    paramnames = parameters,
    covar=sdm_covar
  )


#########################################################################
#-------------------------|  Simulations  |-----------------------------#
#########################################################################

params_guess = c(
  b0=0.013, alpha=1.62, iota=5,
  sigmaSE=0.8,
  br=2,
  mu_EI=.16, mu_IR=0.13, # state transition
  rho_V=150, sd_V=1000,                    # measurement V
  rho_Y=.14, od_Y=0,                   # measurement Y
  eta=.05, N=50000,                   # initial value parameters
  S_0=.95, E_0=.04, I_0=.01, R_0=.0)


# y = covidSEIRsR %>%
#   simulate(params=params_guess, nsim=250, format="data.frame")
# 
# y_avg = y %>% group_by(day) %>% summarize_at(vars(S:R, V, Y), mean)
# 
# 
# observed = NewHaven %>%
#   mutate(actual.cases = Y / params_guess['rho_Y']) %>%
#   select(day, V = V, Y = actual.cases) %>%
#   pivot_longer(c(V, Y))
# 
# y %>% pivot_longer(c(V, Y)) %>%
#   ggplot(aes(x = day, y = value)) +
#   geom_line(aes(color = factor(.id))) +
#   geom_line(data = y_avg %>% pivot_longer(c(V, Y)),
#             size=2, color="blue") +
#   geom_line(data = observed, color="black", size=2) +
#   scale_color_brewer(type = 'qual', palette = 3) +
#   guides(color = FALSE) +
#   facet_wrap(~name, scales="free_y")

#########################################################################
#----------------------|  Particle Filtering  |-------------------------#
#########################################################################

# tic <- Sys.time()
# 
# L_pf=0
# pf=0
# foreach(i=1:10,.combine=c) %dopar% {
#   library(pomp)
#   covidSEIRsR %>% pfilter(params=params_guess,Np=5000)
# } -> pf
# 
# pf %>% logLik() %>% logmeanexp(se=TRUE) -> L_pf
# L_pf
# toc <- Sys.time()


#########################################################################
#-----------------------------| local |---------------------------------#
#########################################################################

registerDoRNG(482947940)
bake(file=paste(cache_address,"/local_search_simple_VY.rds", sep = ''),{
  foreach(i=1:covid_Nreps_local,.combine=c, .errorhandling="remove") %dopar% {
    library(pomp)
    library(tidyverse)
    covidSEIRsR %>%
      mif2(
        params=params_guess,
        Np=covid_Np, Nmif=covid_Nmif,
        cooling.fraction.50=0.5,
        rw.sd=rw.sd(b0=0.02, alpha=0.02, iota=0.02,
                    sigmaSE=0.02, br=0.01, mu_EI=0.00, mu_IR=0.00, 
                    S_0=ivp(0.00), E_0=ivp(0.00), I_0=ivp(0.00), R_0=ivp(0.00), 
                    eta=ivp(0.00), rho_V=0.01, rho_Y=0.0, sd_V=0.5, od_Y=0.0)
      ) %>%
      mif2(cooling.fraction.50=0.3) %>%
      mif2(cooling.fraction.50=0.1)
  } -> mifs_local
  attr(mifs_local,"ncpu") <- getDoParWorkers()
  mifs_local
}) -> mifs_local
t_loc <- attr(mifs_local,"system.time")
ncpu_loc <- attr(mifs_local,"ncpu")

# #plotting the parallel tasks for the local search}
# mifs_local %>%
#   traces() %>%
#   melt() %>%
#   ggplot(aes(x=iteration, y=value, group=L1, color=factor(L1)))+
#   geom_line()+
#   guides(color=FALSE)+
#   facet_wrap(~variable,scales="free_y")


#parallel runs for calculating exact likelihood}
registerDoRNG(900242057)
bake(file=paste(cache_address,"/lik_local_simple_VY.rds", sep = ''),{
  foreach(mf=mifs_local,.combine=rbind, .errorhandling="remove") %dopar% {
    library(pomp)
    library(tidyverse)
    evals <- replicate(covid_Nreps_local, logLik(pfilter(mf,Np=covid_Np)))
    ll <- logmeanexp(evals,se=TRUE)
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- getDoParWorkers()
  results
}) -> results
t_local <- attr(results,"system.time")
ncpu_local <- attr(results,"ncpu")
# pairs(~loglik+b0+alpha+iota+br+rho_V+sd_V,data=results,pch=16)

#########################################################################
#---------------------------| Global |----------------------------------#
#########################################################################

set.seed(${seed})
#initial guesses, creating parameter ranges
runif_design(
  lower=c(b0=0.0001,  rho_V=100,  sd_V=500,
          alpha=0.5,  iota=1.0, sigmaSE=0.4, br=1.0),
  upper=c(b0=1.0000,  rho_V=300, sd_V=3000,
          alpha=1.5,  iota=20., sigmaSE=2, br=20.),
  nseq=$G
) -> guesses

mf1 <- mifs_local[[1]]

fixed_params <- c(N=50000, mu_EI=.16, mu_IR=0.13, rho_Y=0.14, S_0=.95, E_0=.04, I_0=.01, R_0=.0, eta=0.05, od_Y=0)

#run global search}
bake(file=paste(cache_address,"/global_search_simple_VY_$m.rds", sep = ''),{
  registerDoRNG(1270401374)
  foreach(guess=iter(guesses,"row"), .combine=rbind, .errorhandling="remove") %dopar% {
    library(pomp)
    library(tidyverse)
    mf1 %>%
      mif2(params=c(unlist(guess),fixed_params),
           Np=covid_Np, Nmif=covid_Nmif,
           cooling.fraction.50=0.5,
           rw.sd=rw.sd(b0=0.02, alpha=0.02, iota=0.02,
                       sigmaSE=0.02, br=0.01, mu_EI=0.00, mu_IR=0.00, 
                       S_0=ivp(0.00), E_0=ivp(0.00), I_0=ivp(0.00), R_0=ivp(0.00), 
                       eta=ivp(0.00), rho_V=0.01, rho_Y=0.0, sd_V=0.5, od_Y=0.0)) %>%
      mif2(cooling.fraction.50=0.3) %>%
      mif2() %>%
      mif2(cooling.fraction.50=0.1) %>%
      mif2() -> mf
    replicate(
      covid_Nreps_global,
      mf %>% pfilter(Np=covid_Np) %>% logLik()
    ) %>%
      logmeanexp(se=TRUE) -> ll
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- getDoParWorkers()
  results
}) %>%
  filter(is.finite(loglik)) -> results
t_global <- attr(results,"system.time")
t_gncpu_global <- attr(results,"ncpu")



# write.table(results, "./results/covid_params_simple.csv", sep = ",",
#             col.names = names(results),
#             append = T)

write.table(results, "./results/covid_params_simple_VY.csv", sep = ",",
            col.names = !file.exists("./results/covid_params_simple.csv"),
            append = T)
