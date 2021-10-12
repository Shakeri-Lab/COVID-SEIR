library(foreach)
library(doParallel)
library(latex2exp)
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
covid_Np <-          switch(run_level,100, 1e3, 5e4)
covid_Nmif <-        switch(run_level, 10, 100, 100)
covid_Nreps_eval <-  switch(run_level,  2,  10,  10)
covid_Nreps_local <- switch(run_level, 10,  10,  10)
covid_Nreps_global <-switch(run_level, 10,  20,  20)
covid_Nsim <-        switch(run_level, 50, 100, 100) 


set.seed(1350254336)
setwd("/home/mf4yc/SEIRR")

cache_address = "/project/shakeri-lab/cache/cache_abm"

NewHaven = read_csv("/home/mf4yc/SEIRR/Data/abm.csv")
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
NewHaven_cal <- NewHaven_c[1:70,]
NewHaven_pred <- NewHaven_c[71:100,]


####################################################################
#----------------------| Reading Results |-------------------------#
####################################################################

read_csv("/home/mf4yc/SEIRR/results/covid_params_simple_VY.csv") %>%
  filter(is.finite(loglik)) %>%
  select(-X1) %>%
  arrange(-loglik) %>%
  unique() -> res_VY

read_csv("/home/mf4yc/SEIRR/results/covid_params_simple_V.csv") %>%
  filter(is.finite(loglik)) %>%
  select(-X1) %>%
  arrange(-loglik) %>%
  unique() -> res_V

read_csv("/home/mf4yc/SEIRR/results/covid_params_simple_Y.csv") %>%
  filter(is.finite(loglik)) %>%
  select(-X1) %>%
  arrange(-loglik) %>%
  unique() -> res_Y

res <- bind_rows("V & Y" = res_VY, "V" = res_V, "Y" = res_Y, .id = "type")

####################################################################
#----------------------| Poorman Profiles |------------------------#
####################################################################

res %>%
  filter(type=="Y") %>%
  filter(loglik > max(loglik)-3) %>%
  group_by(cut=round(sigmaSE,2)) %>%
  filter(rank(-loglik)<2) %>%
  ggplot(aes(x=sigmaSE, y=loglik, color=type)) +
  geom_point()+
  labs(title = "Y")


####################################################################
#--------------------| Comparing Simulations |---------------------#
####################################################################

bake(file=paste(cache_address,"/local_search_simple_VY.rds", sep = '')) -> mifs_local
mifs_local[[1]] -> SEIR_VY

bake(file=paste(cache_address,"/local_search_simple_V.rds", sep = '')) -> mifs_local
mifs_local[[1]] -> SEIR_V

bake(file=paste(cache_address,"/local_search_simple_Y.rds", sep = '')) -> mifs_local
mifs_local[[1]] -> SEIR_Y


#----------------------------| Getting MLE parameters |-----------------------------#

res %>%
  filter(type=="V & Y") %>%
  filter(loglik == max(loglik)) %>%
  select(-loglik.se, -loglik, -type) -> params_opt_VY

res %>%
  filter(type=="V") %>%
  filter(loglik == max(loglik)) %>%
  select(-loglik.se, -loglik, -type) -> params_opt_V

res %>%
  filter(type=="Y") %>%
  filter(rank(-loglik)==2) %>%
  select(-loglik.se, -loglik, -type) -> params_opt_Y

#-------------------------------------------------------------------------#

params_sim = params_opt_Y
y_Y = SEIR_Y %>%
  simulate(params=params_sim, nsim=1000, format="data.frame")

params_sim = params_opt_V
y_V = SEIR_V %>%
  simulate(params=params_sim, nsim=1000, format="data.frame")

params_sim = params_opt_VY
y_VY = SEIR_VY %>%
  simulate(params=params_sim, nsim=1000, format="data.frame")

y <- bind_rows("VY" = y_VY, "V" = y_V, "Y" = y_Y, .id = "type")

y %>% group_by(day, type) %>% summarize_at(vars(S:R, V, Y), mean) -> y_avg
y %>% group_by(day, type) %>% summarize_at(vars(V, Y), sd) -> y_sd

observed = NewHaven_cal %>%
  mutate(actual.cases = Y / 0.14) %>%
  select(day, V = V, Y = actual.cases) 

# observed = NewHaven_cal %>%
#   mutate(actual.cases = Y / 0.14) %>%
#   select(day, V = V, Y = actual.cases) %>%
#   pivot_longer(c(V, Y))

# y %>% pivot_longer(c(V, Y)) %>%
#   ggplot(aes(x = day, y = value)) +
#   geom_line(aes(color = factor(type))) +
#   geom_line(data = y_avg %>% pivot_longer(c(V, Y)),
#             size=2, color="blue") +
#   geom_line(data = observed, color="black", size=2) +
#   scale_color_brewer(type = 'qual', palette = 3) +
#   guides(color = FALSE) +
#   facet_wrap(~name, scales="free_y")


y_avg$V_low <- y_avg$V - 0.5*y_sd$V
y_avg$V_up  <- y_avg$V + 0.5*y_sd$V
y_avg$Y_low <- y_avg$Y - 0.5*y_sd$Y
y_avg$Y_up  <- y_avg$Y + 0.5*y_sd$Y

observed$V_low <- observed$V
observed$V_up  <- observed$V
observed$Y_low <- observed$Y
observed$Y_up  <- observed$Y
observed$type  <- "data"

y_avg <- y_avg %>%
  select(day, type, V, V_low, V_up, Y, Y_low, Y_up)

y_avg <- rbind(y_avg, observed)


cbPalette <- c("#000000", "#CC0000", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#E69F00")

y_avg %>% 
  group_by(type) %>%
  ggplot(aes(x = day, y = V)) +
  geom_line(aes(color = type), size=1.5) + 
  geom_ribbon(aes(ymin=V_low, ymax=V_up, fill = type), alpha = 0.3)+
  # geom_line(data = observed, aes(x=day, y=V), 
  #           color="black", size=1.5, show.legend = c("data")) +
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette)+
  labs(y="Viral Load") +
  theme_bw() +
  theme(plot.title = element_text(size=16, face="bold", family="Times New Roman", hjust = 0.5),
        axis.title.x = element_text(size = 22, family = "Times New Roman"),
        axis.title.y = element_text(size = 22, family = "Times New Roman"),
        axis.text = element_text( size = 21))


y_avg %>% 
  group_by(type) %>%
  ggplot(aes(x = day, y = Y)) +
  geom_line(aes(color = type), size=1.5) +
  geom_ribbon(aes(ymin=Y_low, ymax=Y_up, fill = type), alpha = 0.3)+
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette)+
  # geom_line(data = observed %>% pivot_wider(names_from = name), 
  #           color="black", size=1.5) +
  labs(y="Reported Cases") +
  theme_bw() +
  theme(plot.title = element_text(size=16, face="bold", family="Times New Roman", hjust = 0.5),
        axis.title.x = element_text(size = 22, family = "Times New Roman"),
        axis.title.y = element_text(size = 22, family = "Times New Roman"),
        axis.text = element_text( size = 21))


x = rpois(1000, 30)
hist(x, 50)


############################################################################
#-----------------------------| Forecast |---------------------------------#
############################################################################


forecast_pomp <- function(h, model, ranges){
  
  sobol_design(
    lower=ranges[,"min"],
    upper=ranges[,"max"],
    nseq=20
  ) -> params
  
  library(foreach)
  library(doParallel)
  library(iterators)
  library(doRNG)
  
  registerDoParallel()
  registerDoRNG(887851050L)
  
  ## ----forecasts2b----------------------------------------------------------
  foreach(p=iter(params,by="row"),
          .inorder=FALSE,
          .combine=bind_rows
  ) %dopar% {
    
    library(pomp)
    
    ## ----forecasts2c----------------------------------------------------------
    M1 <- model
    
    M1 %>% pfilter(params=p, Np=1000, save.states=TRUE) -> pf
    
    ## ----forecasts2d----------------------------------------------------------
    pf %>%
      saved.states() %>% ## latent state for each particle
      tail(1) %>%        ## last timepoint only
      melt() %>%         ## reshape and rename the state variables
      spread(variable,value) %>%
      group_by(rep) %>%
      summarize(S_0=S/(S+E+I+R), E_0=E/(S+E+I+R), I_0=I/(S+E+I+R), R_0=R/(S+E+I+R)) %>%
      gather(variable,value,-rep) %>%
      spread(rep,value) %>%
      column_to_rownames("variable") %>%
      as.matrix() -> x
    
    ## ----forecasts2e1----------------------------------------------------------
    pp <- parmat(unlist(p),ncol(x))
    
    ## ----forecasts2e2----------------------------------------------------------
    M1 %>%
      simulate(params=pp,format="data.frame") %>%
      select(.id,day,Y, V) %>%
      mutate(
        period="calibration",
        loglik=logLik(pf)
      ) -> calib
    
    ## ----forecasts2f----------------------------------------------------------
    M2 <- M1
    time(M2) <- max(time(M1))+seq_len(h)
    timezero(M2) <- max(time(M1))
    
    ## ----forecasts2g----------------------------------------------------------
    pp[rownames(x),] <- x
    
    M2 %>%
      simulate(params=pp,format="data.frame") %>%
      select(.id,day,Y,V) %>%
      mutate(
        period="projection",
        loglik=logLik(pf)
      ) -> proj
    
    ## ----forecasts2h----------------------------------------------------------
    bind_rows(calib,proj) -> sims
    return(sims)
  }
}


simq_f <- function(sims){
  sims %>%
    mutate(weight=exp(loglik-mean(loglik))) %>%
    arrange(day,.id) -> sims
  
  ## ----forecasts2k----------------------------------------------------------
  sims %>%
    filter(day==max(day)) %>%
    summarize(ess=sum(weight)^2/sum(weight^2))
  
  ## ----forecasts2l----------------------------------------------------------
  sims %>%
    group_by(day,period) %>%
    summarize(
      p=c(0.025,0.5,0.975),
      q=quantile(Y,weights=weight,probs=p),
      label=c("lower","median","upper")
    ) %>%
    select(-p) %>%
    spread(label,q) %>%
    ungroup() %>%
    mutate(date=day) -> simq
  
  return(simq)
}

##############################

res_VY %>%
  select(-loglik.se) %>%
  filter(loglik>max(loglik)-0.5*qchisq(df=1, p=0.95)) %>%
  gather(parameters,value) %>%
  group_by(parameters) %>%
  summarize(min=min(value),max=max(value)) %>%
  ungroup() %>%
  filter(parameters!="loglik") %>%
  column_to_rownames("parameters") %>%
  as.matrix() -> ranges_VY

res_V %>%
  select(-loglik.se) %>%
  filter(loglik>max(loglik)-0.5*qchisq(df=1, p=0.95)) %>%
  gather(parameters,value) %>%
  group_by(parameters) %>%
  summarize(min=min(value),max=max(value)) %>%
  ungroup() %>%
  filter(parameters!="loglik") %>%
  column_to_rownames("parameters") %>%
  as.matrix() -> ranges_V

res_Y %>%
  select(-loglik.se) %>%
  filter(loglik>max(loglik)-0.5*qchisq(df=1, p=0.95)) %>%
  gather(parameters,value) %>%
  group_by(parameters) %>%
  summarize(min=min(value),max=max(value)) %>%
  ungroup() %>%
  filter(parameters!="loglik") %>%
  column_to_rownames("parameters") %>%
  as.matrix() -> ranges_Y

sims_VY <- forecast_pomp(30, SEIR_VY, ranges_VY)
sims_V <- forecast_pomp(30, SEIR_V, ranges_V)
sims_Y <- forecast_pomp(30, SEIR_Y, ranges_Y)

simq_VY <- simq_f(sims_VY)
simq_V <- simq_f(sims_V)
simq_Y <- simq_f(sims_Y)

simq <- bind_rows("VY" = simq_VY, "V" = simq_V, "Y" = simq_Y, .id = "type") %>%
  filter(period=='projection') %>%
  select(-date)

#-----------------------| ARIMA |------------------------#
library(forecast)
cases <- NewHaven_c[1:70,]$Y/0.14
viral <- NewHaven_c[1:70,]$V

cases %>% Arima(order=c(4,1,0)) -> arima_case

viral %>% Arima(order=c(1,1,1)) -> arima_viral

case_f_arima <- forecast(arima_case, 30)

data_frame(type='arima', day=71:100, period='projection', lower=case_f_arima$lower[,2], median=case_f_arima$mean, 
           upper=case_f_arima$upper[,2]) -> arima_case_f

arima_case_f$lower[3:30] = 0

simq <- bind_rows(simq, arima_case_f)


observed = NewHaven_c[1:100,] %>%
  mutate(actual.cases = Y / 0.14) %>%
  select(day, V = V, Y = actual.cases)

observed %>% select(-V) %>%
  mutate(median=Y, lower=Y, upper=Y, period="projection", type="data") %>%
  select(-Y) -> observed

simq <- rbind(simq, observed)
#--------------------------------------------------------#

cbPalette <- c("#F0E442", "#000000","#CC0000", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#E69F00")

simq %>%
  # filter(type!='Y') %>%
  # filter(type!='VY') %>%
  ggplot(aes(x=day))+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=type),alpha=0.3,color=NA)+
  geom_line(aes(y=median,color=type), size=1.5)+
  # geom_point(data=observed, mapping=aes(x=day,y=Y),color="black")+
  # geom_line(data=observed, mapping=aes(x=day,y=Y))+
  labs(y="cases") +
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette)+
    labs(y="Case Count") +
  coord_cartesian(ylim = c(0, 500), xlim = c(0,100)) +
  theme_bw()+
  theme(plot.title = element_text(size=16, face="bold", family="Times New Roman", hjust = 0.5),
        axis.title.x = element_text(size = 22, family = "Times New Roman"),
        axis.title.y = element_text(size = 22, family = "Times New Roman"),
        axis.text = element_text( size = 21)) -> p

p.zoom <- ggplot(simq, aes(x=day)) +
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=type),alpha=0.3,color=NA) +
  geom_line(aes(y=median,color=type), size=1.5) +
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  coord_cartesian(xlim=c(65,100), ylim=c(0,50)) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        legend.position = "none")

p +  annotation_custom(ggplotGrob(p.zoom), xmin = 0, xmax = 65, 
                       ymin = 200, ymax = 500)


simq %>%
  filter(type!='Y') %>%
  filter(type!='VY') %>%
  ggplot(aes(x=day))+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=type),alpha=0.4)+
  geom_line(aes(y=median,color=type))+
  # geom_point(data=observed, mapping=aes(x=day,y=Y),color="black")+
  # geom_line(data=observed, mapping=aes(x=day,y=Y))+
  #ylim(0,100)+
  labs(y="cases") +
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette)+
  labs(y="Case Count") +
  coord_cartesian(ylim = c(0, 55), xlim = c(0,100)) +
  theme_bw()+
  theme(plot.title = element_text(size=16, face="bold", family="Times New Roman", hjust = 0.5),
        axis.title.x = element_text(size = 15, family = "Times New Roman"),
        axis.title.y = element_text(size = 15, family = "Times New Roman"),
        axis.text = element_text( size = 12))














