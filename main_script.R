# install.packages("deSolve", dep=TRUE)
# install.packages("gridExtra", dep = TRUE)
# install.packages("reshape", dep = TRUE)

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
source("fun/get_intervention.R")

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
xstart <- data.frame(U = N-I0,
            L = 0,
            I = I0,  
            R = 0,
            Incidence=0, 
            Irecent=0 , 
            Iremote=0)               

fx<-TB.Basic

#run the model
out  <-get_intervention(xstart, params, NA,times, 
                        NA, NA,fx, "Initial", NA) 
#plot
out$lines



#-----------------------------------------------------------------------
#2) Can you modify manually the transmission rate per capita to achieve an
# incidence rate as that of Bangladesh in 2017?  
#-----------------------------------------------------------------------
Inc.country<- 221

params["beta"]<-6.7

#run the model
out0  <-get_intervention(xstart, params, NA,times, 
                        NA, NA,fx, "Initial", NA) 

#plot
dot<-data.frame(Data="Bangladesh",Years=2017, incidence=Inc.country)

p1<-out$lines +
  geom_point(dot, mapping=aes(x=Years, y=incidence, col=Data), size=6, shape=18) 

grid.arrange(p1,out$pie)


########################################################################
## Simulation 0
# Project a baseline

# Prepare for new simulation
sfin       <- tail(out0$out,1)     
params_new <- params
params_old <- params
times_new  <- seq(t.intervention, t.intervention+25 , by=1)
t.interv   <- c(times_new[2], times_new[2]+t.scale)
int_name   <- "Baseline"

fx_scale<-function(t, state, parameters) scale_up(t, state, parameters,t.interv,params,fx)


data0<-get_intervention(sfin, params_new, params_old,times_new,
                             NA, fx_scale,fx, "Baseline", NA) 

grid.arrange(data0$lines ,data0$pie )
 
########################################################################
## Simulation 1
# An Intervention simulating introduction of treatment

# Change parameters for intervention
T.cs <-   1   # Time delay (yrs) between developing symptoms and seeking for care
pDx  <-   0.5 # Probability of being diagnosed once sought care
pTx  <-   0.7 # probability of recieving correct Tx if diagnosed
T.rTx<-   0.5 # 6 months treatment duration

Tx <-pDx*pTx*(1/(T.cs+T.rTx))

params_new["selfcure"]<-selfcure + Tx

data1<-get_intervention(sfin, params_new, params_old,times_new,
                       t.interv, fx_scale,fx, "Treatment", data0$data) 

p1<-data1$lines +
  # EndTb
  geom_hline(yintercept=Inc.country*0.1, linetype="dashed", color = "black", size=1)+
  # Elimination
  geom_hline(yintercept=0.1/1e5, linetype="dashed", color = "red", size=1)


grid.arrange(p1 ,data1$pie )
########################################################################
## Simulation 3
# An Intervention simulating transm reduction

# Change parameters for intervention
params_new["beta"]<-beta * 0

data2<-get_intervention(sfin, params_new, params_old,times_new,
                        t.interv, fx_scale,fx, "Transmission", data1$data) 
p1<-data2$lines +
  # EndTb
  geom_hline(yintercept=Inc.country*0.1, linetype="dashed", color = "black", size=1)+
  # Elimination
  geom_hline(yintercept=0.1/1e5, linetype="dashed", color = "red", size=1)


grid.arrange(p1 ,data2$pie )


########################################################################
## Simulation 3
# An Intervention simulating LTBI treatment

# Change parameters for intervention
params_new["phi"] <-  0.1*(1/T.lfx) 

data3<-get_intervention(sfin, params_new, params_old,times_new,
                        t.interv, fx_scale,fx, "Transmission", data2$data) 

p1<-data3$lines +
  # EndTb
  geom_hline(yintercept=Inc.country*0.1, linetype="dashed", color = "black", size=1)+
  # Elimination
  geom_hline(yintercept=0.1/1e5, linetype="dashed", color = "red", size=1)


grid.arrange(p1 ,data3$pie )

