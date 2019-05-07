
TB.Basic <- function (t, state, parameters) {
  
 
  with(as.list(c(state,parameters)),             
       {
         N      <- U + L + I + R                                      # Total population
         births <- I * mu_tb + N * mu                                 # Births for stable population
         lambda <-   beta * I/N                                       # foi
         
         dU <- births - U * (lambda+mu)                              # Uninfected state
         
         dL <-  U * lambda * (1-fast) + R * (lambda * (1-fast) * imm) - 
           L * (mu + phi)                                               # LTBI
         
         dI <-  U * (lambda * fast) + R * (lambda * fast * imm) + 
           L * phi + R * relapse - I * (mu + mu_tb + selfcure)          # Infectious
         
         dR <-  I * selfcure - R * (lambda * imm + relapse + mu)        # recovered
         
         dIncidence <- U * (lambda * fast) + R * (lambda * fast * imm) + L * phi + R * relapse 
         dIrecent   <- U * (lambda * fast) + R * (lambda * fast * imm)
         +  R * relapse                                               # 
         dIremote   <- L * phi +  R * relapse                                          # 
         
         # 
         dx <- c(dU, dL, dI, dR, dIncidence, dIrecent , dIremote)
         list(dx)
       }
  )
}