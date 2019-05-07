 # install.packages("deSolve", dep=TRUE)
 # install.packages("gridExtra", dep = TRUE)
rm(list=ls())
# dev.off()

#Load packages
library(deSolve)
## Warning: package 'deSolve' was built under R version 3.4.4
library(gridExtra)

#Source functions
source("fun/TB.Basic.R")


# Model Parameters
T.lat    <- 67              # latency
T.durat  <- 6.5             # duration of infectious period
T.lxpc   <- 67              # Life expec
T.tbout  <- 6.5             # mort Tb untreated

phi      <- 0.1*(1/T.lat)  # 10% of slow latent progress over a lifetime
selfcure <- 1/T.durat       # TB cured spontaneusly
mu       <- 1/T.lxpc        # Background mortality
mu_tb    <- 1/T.tbout       # TB mortality
beta     <- 3              # Transmission rate per capita
fast     <- 0.14            # fraction fast progressing to active TB
imm      <- 0.5             # Infectiousness decline (partial immunity)
relapse  <- 0.0032           # Relapse rate


t.intervention = 400      # ??

# Initial Conditions
N <- 1                  # ??
I0 <-1e-6               # ??

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
out <- as.data.frame(ode(y = xstart, times = times, 
                         func = TB.Basic, parms = params))   #


#Plot some results
plot(out$time, out$U, col='blue',type='l',ylim = c(0,max(out$U)),
     xlab ='Years', ylab = 'Number')
lines(out$time, out$L, col='black')
lines(out$time, out$I, col='red')
lines(out$time, out$R, col='green')

legend('topleft',c('U','L','I','R'),lwd=2,col=c('blue','black','red','green'))


N <- out$U+out$L+out$I+out$R  

rate.inc<- 1e5*(diff(out$Incidence)/N[1:length(N)-1])
time<-out$time[1:length(out$time)-1]
plot(time, rate.inc, col='blue',type='l',ylim = c(0,max(rate.inc)),
     xlab ='Years', ylab = 'TB Incidence x 100K')


########################################################################
## Simulation 2
# An Intervention simulating introduction of treatment

params2 <- params
params2["selfcure"]<-1/(T.durat*0.01)
        
times2  = seq(t.intervention, t.intervention+25 , by=1)
# Starting conditions
xstart2 <- c(U = tail(out$U,1), 
             L = tail(out$L,1),  
             I = tail(out$I,1), 
             R = tail(out$R,1),
             Incidence= tail(out$Incidence,1),
             Irecent= tail(out$Irecent,1),  
             Iremote= tail(out$Iremote,1)) 

out2 <- as.data.frame(ode(y = xstart2, times = times2, 
                          func = TB.Basic, parms = params2))  # ??

N <- out2$U+out2$L+out2$I+out2$R  

rate.inc2<- 1e5*(diff(out2$Incidence)/N[1:length(N)-1])
time2<-out2$time[1:length(out2$time)-1]
plot(time2, rate.inc2, col='blue',type='l',ylim = c(0,max(rate.inc2)),
     xlab ='Years', ylab = 'TB Incidence x 100K')


########################################################################
## Simulation 3
# An Intervention simulating transm reduction

params3 <- params2
params3["beta"] <-beta*0

times3  = seq(t.intervention, t.intervention+25 , by=1)
# Starting conditions
xstart3 <- c(U = tail(out2$U,1), 
             L = tail(out2$L,1),  
             I = tail(out2$I,1), 
             R = tail(out2$R,1),
             Incidence= tail(out2$Incidence,1),
             Irecent= tail(out2$Irecent,1),  
             Iremote= tail(out2$Iremote,1)) 

out3 <- as.data.frame(ode(y = xstart3, times = times3, 
                          func = TB.Basic, parms = params3))  # ??

N <- out3$U+out3$L+out3$I+out3$R  

rate.inc3<- 1e5*(diff(out3$Incidence)/N[1:length(N)-1])
time3<-out3$time[1:length(out3$time)-1]
plot(time3, rate.inc3, col='blue',type='l',ylim = c(0,max(rate.inc3)),
     xlab ='Years', ylab = 'TB Incidence x 100K')


########################################################################
## Simulation 4
# An Intervention simulating LTBI treatment

params4 <- params3
params4["phi"] <-  0.01*(1/T.lat) 

times4  = seq(t.intervention, t.intervention+25 , by=1)
# Starting conditions
xstart4 <- c(U = tail(out2$U,1), 
             L = tail(out2$L,1),  
             I = tail(out2$I,1), 
             R = tail(out2$R,1),
             Incidence= tail(out2$Incidence,1),
             Irecent= tail(out2$Irecent,1),  
             Iremote= tail(out2$Iremote,1)) 

out4 <- as.data.frame(ode(y = xstart4, times = times4, 
                          func = TB.Basic, parms = params4))  # ??

N <- out4$U+out4$L+out4$I+out4$R  

rate.inc4<- 1e5*(diff(out4$Incidence)/N[1:length(N)-1])
time4<-out4$time[1:length(out4$time)-1]
plot(time4, rate.inc4, col='blue',type='l',ylim = c(0,max(rate.inc4)),
     xlab ='Years', ylab = 'TB Incidence x 100K')

########################################################################
############################## Plot all
## Counterfactual

times0  = seq(t.intervention, t.intervention+25 , by=1)
# Starting conditions
xstart0 <- c(U = tail(out$U,1), 
             L = tail(out$L,1),  
             I = tail(out$I,1), 
             R = tail(out$R,1),
             Incidence= tail(out$Incidence,1),
             Irecent= tail(out$Irecent,1),  
             Iremote= tail(out$Iremote,1)) 

out0 <- as.data.frame(ode(y = xstart0, times = times0, 
                          func = TB.Basic, parms = params))  # ??

N <- out0$U+out0$L+out0$I+out0$R  
rate.inc0<- 1e5*(diff(out0$Incidence)/N[1:length(N)-1])



plot(time2, rate.inc0, col='blue',type='l',ylim = c(0,max(rate.inc0)),
     xlab ='Years', ylab = 'TB Incidence x 100K')
lines(time2, rate.inc2, col='black')
lines(time3, rate.inc3, col='red')
lines(time4, rate.inc4, col='green')




