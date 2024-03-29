
get_intervention<- function (sfin, params_new, params_old,times_new,
                             t.interv, fx_scale,fx_basic, int_name, data_stub) {
  
  # Starting conditions
  xstart <- c(U = sfin$U, 
              L = sfin$L,  
              I = sfin$I, 
              R = sfin$R,
              Incidence= sfin$Incidence,
              Irecent=   sfin$Irecent,  
              Iremote=   sfin$Iremote) 
  
  #Run the model
  if (is.na(t.interv))
  {
    fx<-fx_basic
  }  else {
    fx<-fx_scale
    }
  
  out <- as.data.frame(ode(y = xstart, times = times_new, 
                           func = fx, parms = params_new))  # 
  # Model output
  N            <- out$U+out$L+out$I+out$R  
  rate.inc    <- 1e5*(diff(out$Incidence)/N[1:length(N)-1])
  fr.remo     <- diff(out$Iremote)/diff(out$Incidence)
  time         <- out$time[1:length(out$time)-1]
  dat         <- data.frame(Years=time+(2019-400), Incidence=rate.inc)
  dat$Sim     <- int_name
  
  # Create plot
  if (is.na(data_stub))
  {
    data<-dat
  }
  else 
  {
    data  <-rbind(data_stub, dat)
  }
  
  
  remote<-fr.remo
  titl  <-int_name
  
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
  
  
 
  
  output<-list("out"=out, "lines"=p1, "pie"=pie, "data"=data)
  
  return(output)
  
}