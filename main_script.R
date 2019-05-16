# install.packages("deSolve", dep=TRUE)
# install.packages("gridExtra", dep = TRUE)
rm(list=ls())
# dev.off()

#Load packages
library(deSolve)
## Warning: package 'deSolve' was built under R version 3.4.4
library(gridExtra)
library(ggplot2)
library(reshape)

#Source functions
source("fun/TB.Basic.R")
source("fun/scale_up.R")

#-----------------------------------------------------------------------
#1.a) Given lecture and previous statement complete the table of parameters 

#1.b) looking at the resulting parameters could you infer what is the average
#duration of the infectious period? 
#-----------------------------------------------------------------------

# Model Parameters
T.lfx    <- 72              # life expectancy
phi      <- 0.1*(1/T.lfx)   # 10% of slow latent progress over a lifetime
selfcure <- 0.5*(1/3)       # TB cured spontaneusly
mu       <- 1/T.lfx         # Background mortality
mu_tb    <- 0.5*(1/3)       # TB mortality
beta     <- 6.5             # Transmission rate per capita
fast     <- 0.1             # fraction fast progressing to active TB
imm      <- 0.5             # Infectiousness decline (partial immunity)
relapse  <- 0.005           # Relapse rate
q        <- 0               # immunized

t.intervention <- 400      # years of simulation
t.scale        <- 3        # Scaling up time of interventions

# Initial Conditions
N <- 1                  # Total population equal to unity 
I0 <-1e-6               # TB seed at time 0

#Prepare tu run
params <- c(phi=phi, selfcure=selfcure, mu=mu, mu_tb=mu_tb
            ,beta=beta, fast=fast, imm=imm, relapse=relapse)   # running parameters


times  = seq(0, t.intervention, by=1)          # time scale

# Initial conditions
xstart <- c(U = N-I0,
            L = 0,
            I = I0,  
            R = 0,
            Incidence=0, 
            Irecent=0 , 
            Iremote=0)               

#Run the model
out0 <- as.data.frame(ode(y = xstart, times = times, 
                         func = TB.Basic, parms = params))   #

# Model output
N       <- out0$U+out0$L+out0$I+out0$R  
rate.inc<- 1e5*(diff(out0$Incidence)/N[1:length(N)-1])
fr.remo <- diff(out0$Iremote)/diff(out0$Incidence)
time    <- out0$time[1:length(out0$time)-1]

# Order in dataframe for plotting
dat<-data.frame(Years=time+(2019-400), incidence=rate.inc)

#Create plot
p<- ggplot(data=dat, mapping = aes(x=Years, y=incidence))
p + 
  geom_line(col="blue",  size=1.2) +
  ggtitle ('TB Incidence') +
  theme_bw() + ylab('Rate per 100,000 pop')


#-----------------------------------------------------------------------
#2) Can you modify manually the transmission rate per capita to achieve an
# incidence rate as that of Bangladesh in 2017?  
#-----------------------------------------------------------------------
Inc.Bangladesh<- 221

dot<-data.frame(Data="Bangladesh",Years=2017, incidence=Inc.Bangladesh)


p1<-p + 
  geom_line(col="blue",  size=1.2) +
  ggtitle ('TB Incidence') +
  theme_bw() + ylab('Rate per 100,000 pop') +
  geom_point(dot, mapping=aes(x=Years, y=incidence, col=Data), size=6, shape=18) 
  

# Pie of remote incidence in 2017
df <- data.frame(
  Source = c("Recent", "Remote"),
  value  = c(1-tail(fr.remo,1),tail(fr.remo,1))
)

mycols <- c("#0073C2FF", "#EFC000FF")
pie<- ggplot(df, aes(x="", y=value, fill=Source))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0)+
  scale_fill_manual(values = mycols) +
  theme_void()+
  ggtitle ('Baseline') 
  

grid.arrange(p1,pie, ncol=2, nrow =2)

########################################################################
## Simulation 0
# Project a baseline

# Prepare for new simulation
sfin       <- tail(out0,1)     
params_new <- params
params_old <- params
times_new  <- seq(t.intervention, t.intervention+25 , by=1)
t.interv   <- c(times_new[2], times_new[2]+t.scale)


# Starting conditions
xstart <- c(U = sfin$U, 
            L = sfin$L,  
            I = sfin$I, 
            R = sfin$R,
            Incidence= sfin$Incidence,
            Irecent=   sfin$Irecent,  
            Iremote=   sfin$Iremote) 


#Create function handle
fx<-TB.Basic

#Run the model
out <- as.data.frame(ode(y = xstart, times = times_new, 
                          func = fx, parms = params_new))  # ??
# Model output
N            <- out$U+out$L+out$I+out$R  
rate.inc    <- 1e5*(diff(out$Incidence)/N[1:length(N)-1])
fr.remo0     <- diff(out$Iremote)/diff(out$Incidence)
time         <- out$time[1:length(out$time)-1]
dat0         <- data.frame(Years=time+(2019-400), Incidence=rate.inc)
dat0$Sim     <- "Baseline"

########################################################################
## Simulation 1
# An Intervention simulating introduction of treatment

# Prepare for new simulation
sfin       <- tail(out0,1)     
params_new <- params
params_old <- params
times_new  <- seq(t.intervention, t.intervention+25 , by=1)
t.interv   <- c(times_new[2], times_new[2]+t.scale)

# Change parameters for intervention
T.cs <-   1   # Time delay (yrs) between developing symptoms and seeking for care
pDx  <-   0.5 # Probability of being diagnosed once sought care
pTx  <-   0.7 # probability of recieving correct Tx if diagnosed
T.rTx<-   0.5 # 6 months treatment duration

Tx <-pDx*pTx*(1/(T.cs+T.rTx))

params_new["selfcure"]<-selfcure + Tx

# Starting conditions
xstart <- c(U = sfin$U, 
            L = sfin$L,  
            I = sfin$I, 
            R = sfin$R,
            Incidence= sfin$Incidence,
            Irecent=   sfin$Irecent,  
            Iremote=   sfin$Iremote) 


#Create function handle
fx<-TB.Basic
scalefx<-function(t, state, parameters) scale_up(t, state, parameters,t.interv,params,fx)

#Run the model
out <- as.data.frame(ode(y = xstart, times = times_new, 
                          func = scalefx, parms = params_new))  # ??
# Model output
N            <- out$U+out$L+out$I+out$R  
rate.inc     <- 1e5*(diff(out$Incidence)/N[1:length(N)-1])
fr.remo1     <- diff(out$Iremote)/diff(out$Incidence)
time         <- out$time[1:length(out$time)-1]
dat1         <- data.frame(Years=time+(2019-400), Incidence=rate.inc)
dat1$Sim     <- "Treatment"

# Create plot
data  <-rbind(dat0, dat1)
remote<-fr.remo1
titl  <-"Treatment"

p<- ggplot(data=data, mapping = aes(x=Years, y=Incidence, col=Sim))

p1<-p + 
  geom_line(size=1.2) +
  ggtitle ('TB Incidence') +
  theme_bw() + ylab('Rate per 100,000 pop')+
  ylim(0,max(data$Incidence))


# Pie chart of remote vs recent incidence 
df <- data.frame(
  Source = c("Recent", "Remote"),
  value  = c(1-tail(remote,1),tail(remote,1))
)

mycols <- c("#0073C2FF", "#EFC000FF")
pie<- ggplot(df, aes(x="", y=value, fill=Source))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0)+
  scale_fill_manual(values = mycols) +
  theme_void()+
  ggtitle (titl) 


grid.arrange(p1,pie, ncol=2, nrow =2)


########################################################################
## Simulation 3
# An Intervention simulating transm reduction

# Prepare for new simulation
sfin       <- tail(out0,1)     
params_new <- params_new
params_old <- params
times_new  <- seq(t.intervention, t.intervention+25 , by=1)
t.interv   <- c(times_new[2], times_new[2]+t.scale)

# Change parameters for intervention
params_new["beta"]<-beta * 0

# Starting conditions
xstart <- c(U = sfin$U, 
            L = sfin$L,  
            I = sfin$I, 
            R = sfin$R,
            Incidence= sfin$Incidence,
            Irecent=   sfin$Irecent,  
            Iremote=   sfin$Iremote) 


#Create function handle
fx<-TB.Basic
scalefx<-function(t, state, parameters) scale_up(t, state, parameters,t.interv,params,fx)

#Run the model
out <- as.data.frame(ode(y = xstart, times = times_new, 
                         func = scalefx, parms = params_new))  # ??
# Model output
N            <- out$U+out$L+out$I+out$R  
rate.inc     <- 1e5*(diff(out$Incidence)/N[1:length(N)-1])
fr.remo2      <- diff(out$Iremote)/diff(out$Incidence)
time         <- out$time[1:length(out$time)-1]
dat2         <- data.frame(Years=time+(2019-400), Incidence=rate.inc)
dat2$Sim     <- "Transmission"

# Create plot
data  <-rbind(data, dat2)
remote<-fr.remo2
titl  <-"Transmission"

p<- ggplot(data=data, mapping = aes(x=Years, y=Incidence, col=Sim))

p1<-p + 
  geom_line(size=1.2) +
  ggtitle ('TB Incidence') +
  theme_bw() + ylab('Rate per 100,000 pop')+
  ylim(0,max(data$Incidence))


# Pie chart of remote vs recent incidence 
df <- data.frame(
  Source = c("Recent", "Remote"),
  value  = c(1-tail(remote,1),tail(remote,1))
)

mycols <- c("#0073C2FF", "#EFC000FF")
pie<- ggplot(df, aes(x="", y=value, fill=Source))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0)+
  scale_fill_manual(values = mycols) +
  theme_void()+
  ggtitle (titl) 


grid.arrange(p1,pie, ncol=2, nrow =2)


########################################################################
## Simulation 3
# An Intervention simulating LTBI treatment


# Prepare for new simulation
sfin       <- tail(out0,1)     
params_new <- params_new
params_old <- params
times_new  <- seq(t.intervention, t.intervention+25 , by=1)
t.interv   <- c(times_new[2], times_new[2]+t.scale)

# Change parameters for intervention
params_new["phi"] <-  0.01*(1/T.lfx) 

# Starting conditions
xstart <- c(U = sfin$U, 
            L = sfin$L,  
            I = sfin$I, 
            R = sfin$R,
            Incidence= sfin$Incidence,
            Irecent=   sfin$Irecent,  
            Iremote=   sfin$Iremote) 


#Create function handle
fx<-TB.Basic
scalefx<-function(t, state, parameters) scale_up(t, state, parameters,t.interv,params,fx)

#Run the model
out <- as.data.frame(ode(y = xstart, times = times_new, 
                         func = scalefx, parms = params_new))  # ??
# Model output
N            <- out$U+out$L+out$I+out$R  
rate.inc     <- 1e5*(diff(out$Incidence)/N[1:length(N)-1])
fr.remo3      <- diff(out$Iremote)/diff(out$Incidence)
time         <- out$time[1:length(out$time)-1]
dat3         <- data.frame(Years=time+(2019-400), Incidence=rate.inc)
dat3$Sim     <- "LTBI"

# Create plot
data  <-rbind(data, dat3)
remote<-fr.remo3
titl  <-"LTBI"

p<- ggplot(data=data, mapping = aes(x=Years, y=Incidence, col=Sim))

p1<-p + 
  geom_line(size=1.2) +
  ggtitle ('TB Incidence') +
  theme_bw() + ylab('Rate per 100,000 pop')+
  ylim(0,max(data$Incidence))


# Pie chart of remote vs recent incidence 
df <- data.frame(
  Source = c("Recent", "Remote"),
  value  = c(1-tail(remote,1),tail(remote,1))
)

mycols <- c("#0073C2FF", "#EFC000FF")
pie<- ggplot(df, aes(x="", y=value, fill=Source))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0)+
  scale_fill_manual(values = mycols) +
  theme_void()+
  ggtitle (titl) 


grid.arrange(p1,pie, ncol=2, nrow =2)

########################################################################
## Counterfactual

# Prepare for new simulation
# params_new <- params
# params_old <- params
# 
# times_new  <- seq(t.intervention, t.intervention+25 , by=1)
# 
# 
# out0 <- as.data.frame(ode(y = xstart, times = times_new, 
#                           func = TB.Basic, parms = params_new))  # ??
# 
# N <- out0$U+out0$L+out0$I+out0$R  
# rate.inc0<- 1e5*(diff(out0$Incidence)/N[1:length(N)-1])
# fr.remo0 <- diff(out0$Iremote)/diff(out0$Incidence)
# 
# ############################################
# ############################## Plot all
# 
# # Trajectories
# Years<-time2+1617
# allruns<-data.frame(Years , rate.inc0, rate.inc2, rate.inc3, rate.inc4)
# sz<-1.2
# p<- ggplot(data=allruns, mapping = aes(x=Years, y=value))
# 
# p1<-p + 
#   geom_line(aes(y=rate.inc0, colour="Baseline"),  size=sz)+
#   geom_line(aes(y=rate.inc2, colour="Treatment"),  size=sz)+
#   geom_line(aes(y=rate.inc3, colour="Transmission"),  size=sz)+
#   geom_line(aes(y=rate.inc4, colour="LTBI"), size=sz)+
#   scale_colour_manual(name="Intervention",
#                       breaks=c("Baseline","Treatment","Transmission","LTBI"), 
#                       values = c("Baseline"="indianred2", "Treatment"="yellow3", 
#                                  "Transmission"="springgreen4", "LTBI"="royalblue")) +
#   geom_hline(yintercept=rate.inc0[1]*0.1, linetype="dashed", color = "black") +
#   ylim(0,max(rate.inc0))+
#   ggtitle ('TB Interventions') +
#   theme_bw() + ylab('Rate per 100,000')
# 
# 
# # Bars
# remoteI<-data.frame(
#   Intervention=factor(c("Baseline","Treatment","Transmission","LTBI"),
#                       levels=c("Baseline","Treatment","Transmission","LTBI")),
#   remote=c(tail(fr.remo0, n=1), tail(fr.remo2, n=1),
#            tail(fr.remo3, n=1) , tail(fr.remo4, n=1)) )
# 
# remoteI$remote<-remoteI$remote*100
# 
# b<- ggplot(data=remoteI, mapping=aes(x=Intervention, y=remote, fill=Intervention))
# p2<- b +
#   geom_bar(stat="identity")+
#   scale_fill_manual(values = c("Baseline"="indianred2", 
#                                "Treatment"="yellow3", 
#                                "Transmission"="springgreen4", 
#                                "LTBI"="royalblue")) +
#   ylab('Fraction Remote') +
#   ggtitle ('Incidence from remote source') +
#   theme_bw()
#   
#   
# 
# 
# grid.arrange(p1,p2, ncol=1, nrow=2)
# 
# 

# plot(time2, rate.inc0, col='blue',type='l',ylim = c(0,max(rate.inc0)),
#      xlab ='Years', ylab = 'TB Incidence x 100K')
# lines(time2, rate.inc2, col='black')
# lines(time3, rate.inc3, col='red')
# lines(time4, rate.inc4, col='green')
# lines(c(400,430),c(rate.inc0[1]*0.1 , rate.inc0[1]*0.1 ),lty=3)



