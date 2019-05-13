scale_up<- function (t, state, parameters,t.interv,parameters_old,fx) {
  
  
  scale <- min((t-t.interv[1])/(t.interv[2]-t.interv[1]),1); 
  
  if (scale<0) 
    {
    scale=0
    }
  
  pars_scaled <- parameters;
  
  pars_scaled <- parameters_old + scale*(parameters-parameters_old)
  
  
  return(fx(t, state,pars_scaled))
   
}