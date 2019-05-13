# install.packages("deSolve", dep=TRUE)
# install.packages("gridExtra", dep = TRUE)
rm(list=ls())
# dev.off()

#Load packages
library(deSolve)
## Warning: package 'deSolve' was built under R version 3.4.4
library(gridExtra)
library(ggplot2)

#Source functions
source("fun/TB.Basic.R")
source("fun/scale_up.R")


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


t.intervention <- 400      # ??
t.scale        <- 3        # Scaling up time

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

# Prepare for new simulation
sfin       <- out     
params_new <- params
params_old <- params
times_new  <- seq(t.intervention, t.intervention+25 , by=1)
t.interv   <- c(times_new[2], times_new[2]+t.scale)

# Change parameters for intervention
params_new["selfcure"]<-1/(T.durat*0.5)

# Starting conditions
xstart <- c(U = tail(sfin$U,1), 
            L = tail(sfin$L,1),  
            I = tail(sfin$I,1), 
            R = tail(sfin$R,1),
            Incidence= tail(sfin$Incidence,1),
            Irecent= tail(sfin$Irecent,1),  
            Iremote= tail(sfin$Iremote,1)) 


#Create function handle
fx<-TB.Basic
scalefx<-function(t, state, parameters) scale_up(t, state, parameters,t.interv,params,fx)

out2 <- as.data.frame(ode(y = xstart, times = times_new, 
                          func = scalefx, parms = params_new))  # ??

N <- out2$U+out2$L+out2$I+out2$R  

rate.inc2<- 1e5*(diff(out2$Incidence)/N[1:length(N)-1])
fr.remo2 <- diff(out2$Iremote)/diff(out2$Incidence)
time2<-out2$time[1:length(out2$time)-1]
plot(time2, rate.inc2, col='blue',type='l',ylim = c(0,max(rate.inc2)),
     xlab ='Years', ylab = 'TB Incidence x 100K')


########################################################################
## Simulation 3
# An Intervention simulating transm reduction

# Prepare for new simulation
params_new <- params_new
params_old <- params


# Change parameters for intervention
params_new["beta"] <-beta*0


# Create function handle
fx<-TB.Basic
scalefx<-function(t, state, parameters) scale_up(t, state, parameters,t.interv,params_old,fx)


out3 <- as.data.frame(ode(y = xstart, times = times_new, 
                          func = scalefx, parms = params_new))  # ??

N <- out3$U+out3$L+out3$I+out3$R  

rate.inc3<- 1e5*(diff(out3$Incidence)/N[1:length(N)-1])
fr.remo3 <- diff(out3$Iremote)/diff(out3$Incidence)
time3<-out3$time[1:length(out3$time)-1]
plot(time3, rate.inc3, col='blue',type='l',ylim = c(0,max(rate.inc3)),
     xlab ='Years', ylab = 'TB Incidence x 100K')


########################################################################
## Simulation 4
# An Intervention simulating LTBI treatment

# Prepare for new simulation
params_new <- params_new
params_old <- params


# Change parameters for intervention
params_new["phi"] <-  0.01*(1/T.lat) 



# Create function handle
fx<-TB.Basic
scalefx<-function(t, state, parameters) scale_up(t, state, parameters,t.interv,params_old,fx)


out4 <- as.data.frame(ode(y = xstart, times = times_new, 
                          func = scalefx, parms = params_new))  # ??

N <- out4$U+out4$L+out4$I+out4$R  

rate.inc4<- 1e5*(diff(out4$Incidence)/N[1:length(N)-1])
fr.remo4 <- diff(out4$Iremote)/diff(out4$Incidence)
time4<-out4$time[1:length(out4$time)-1]
plot(time4, rate.inc4, col='blue',type='l',ylim = c(0,max(rate.inc4)),
     xlab ='Years', ylab = 'TB Incidence x 100K')

########################################################################
## Counterfactual

# Prepare for new simulation
params_new <- params
params_old <- params

times_new  <- seq(t.intervention, t.intervention+25 , by=1)


out0 <- as.data.frame(ode(y = xstart, times = times_new, 
                          func = TB.Basic, parms = params_new))  # ??

N <- out0$U+out0$L+out0$I+out0$R  
rate.inc0<- 1e5*(diff(out0$Incidence)/N[1:length(N)-1])
fr.remo0 <- diff(out0$Iremote)/diff(out0$Incidence)

############################################
############################## Plot all

# Trajectories
Years<-time2+1617
allruns<-data.frame(Years , rate.inc0, rate.inc2, rate.inc3, rate.inc4)
sz<-1.2
p<- ggplot(data=allruns, mapping = aes(x=Years, y=value))

p1<-p + 
  geom_line(aes(y=rate.inc0, colour="Baseline"),  size=sz)+
  geom_line(aes(y=rate.inc2, colour="Treatment"),  size=sz)+
  geom_line(aes(y=rate.inc3, colour="Transmission"),  size=sz)+
  geom_line(aes(y=rate.inc4, colour="LTBI"), size=sz)+
  scale_colour_manual(name="Intervention",
                      breaks=c("Baseline","Treatment","Transmission","LTBI"), 
                      values = c("Baseline"="indianred2", "Treatment"="yellow3", 
                                 "Transmission"="springgreen4", "LTBI"="royalblue")) +
  geom_hline(yintercept=rate.inc0[1]*0.1, linetype="dashed", color = "black") +
  ylim(0,max(rate.inc0))+
  ggtitle ('TB Interventions') +
  theme_bw() + ylab('Rate per 100,000')


# Bars
remoteI<-data.frame(
  Intervention=factor(c("Baseline","Treatment","Transmission","LTBI"),
                      levels=c("Baseline","Treatment","Transmission","LTBI")),
  remote=c(tail(fr.remo0, n=1), tail(fr.remo2, n=1),
           tail(fr.remo3, n=1) , tail(fr.remo4, n=1)) )

remoteI$remote<-remoteI$remote*100

b<- ggplot(data=remoteI, mapping=aes(x=Intervention, y=remote, fill=Intervention))
p2<- b +
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("Baseline"="indianred2", 
                               "Treatment"="yellow3", 
                               "Transmission"="springgreen4", 
                               "LTBI"="royalblue")) +
  ylab('Fraction Remote') +
  ggtitle ('Incidence from remote source') +
  theme_bw()
  
  


grid.arrange(p1,p2, ncol=1, nrow=2)



# plot(time2, rate.inc0, col='blue',type='l',ylim = c(0,max(rate.inc0)),
#      xlab ='Years', ylab = 'TB Incidence x 100K')
# lines(time2, rate.inc2, col='black')
# lines(time3, rate.inc3, col='red')
# lines(time4, rate.inc4, col='green')
# lines(c(400,430),c(rate.inc0[1]*0.1 , rate.inc0[1]*0.1 ),lty=3)



