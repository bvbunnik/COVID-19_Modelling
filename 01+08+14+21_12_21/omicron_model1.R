library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("tidyverse"); library("RColorBrewer")
library("gridExtra"); library("cowplot")
rm(list=ls())

#Function for the generation time - a function of R0 and the doubling time
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

GenTime(20,1.1)

GT = 5
gamma = 0.2
R0_delta = 2.0
R0_omicron = 6.0

SIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
      dS = -vr*S - beta_SI1*S*I1S  - beta_SI2*S*I2S 
      dV = vr*S - (1-eff1)*beta_VI1*V*(I1S+I1V) - (1-eff2)*beta_VI2*V*(I2S+I2V)
      dI1V = (1-eff1)*beta_VI1*V*(I1S+I1V) - gamma*I1V 
      dI2V = (1-eff2)*beta_VI2*V*(I2S+I2V) - gamma*I2V
      dI1S = beta_SI1*S*I1S - gamma*I1S 
      dI2S = beta_SI2*S*I2S - gamma*I2S
      dR1 = gamma*I1S + gamma*I1V + gamma*I1R2 - beta_R1I2*R1*I2S
      dR2 = gamma*I2S + gamma*I2V + gamma*I2R1 - beta_R2I1*R2*I1S
      dI2R1 = beta_R1I2*R1*I2S - gamma*I2R1
      dI1R2 = beta_R2I1*R2*I1S - gamma*I1R2
      
      #H = hosp_frac_IV*I1V + hosp_frac_IV*I2V + hosp_frac_IS*I1S + 2*hosp_frac_IS*I2S + hosp_frac_IS*I2R1 + hosp_frac_IV*I1R2
      H = IHR1*I1V + IHR2*I2V + IHR3*I1S + IHR4*I2S + IHR5*I2R1 + IHR6*I1R2
      H_IHR2 = IHR2*I2V
      return(list(c(dS,dV,dI1V,dI2V,dI1S,dI2S,dR1,dR2,dI2R1,dI1R2),"H"=H, "H_IHR2" = H_IHR2))
    })
}

times <- seq(0,180,by = 1)


init = c(S=(0.48-1.82e-5), V=0.2, I1V=0.0,I2V=0,I1S=0.02,I2S=1.82e-5,R1=0.3,R2=0,I2R1=0,I1R2=0)
parms = list(gamma = gamma,
             beta_SI1 = R0_delta*gamma,
             beta_SI2 = 2*R0_delta*gamma,
             beta_VI1 = R0_delta*gamma,
             beta_VI2 = 2*R0_delta*gamma,
             beta_R1I2 = 2*R0_delta*gamma,
             beta_R2I1 = 0.0,
             vr = 0.0045,
             eff1 = 0.5,
             eff2 = 0.0,
             IHR1 = 0.0001,
             IHR2 = 0.001,
             IHR3 = 0.001,
             IHR4 = 2*0.001,
             IHR5 = 0.001,
             IHR6 = 0.0001)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

out_m = melt(out, id.vars = "time")
 
ggplot(out_m, aes(x=time, y=value)) + 
  geom_line(size=1.1, colour="darkred") + 
  facet_wrap(~variable, scales = "free_y") +
  ylab("Proportion of population")
  

library("fast")

para_models = fast_parameters(minimum = c(0.1*1.82e-5,2.0,1*0.001,1*0.0001,0.5*0.0045,0), maximum = c(10*1.82e-5,8.0,4*0.001,10*0.0001,2*0.0045,0.5), names = c("I2_0", "R0_2","IHR4","IHR2","vr","eff2"), factor=9)

total_H = vector(mode = "numeric", length=nrow(para_models))
for (i in 1:nrow(para_models)){
  init = c(S=(0.48-para_models[i,1]), V=0.2, I1V=0.0,I2V=0,I1S=0.02,I2S=para_models[i,1],R1=0.3,R2=0,I2R1=0,I1R2=0)
  
  parms = list(gamma = gamma,
               beta_SI1 = R0_delta*gamma,
               beta_SI2 = para_models[i,2]*gamma,
               beta_VI1 = R0_delta*gamma,
               beta_VI2 = para_models[i,2]*gamma,
               beta_R1I2 = para_models[i,2]*gamma,
               beta_R2I1 = 0.0,
               vr = para_models[i,5],
               eff1 = 0.5,
               eff2 = para_models[i,6],
               IHR1 = 0.0001,
               IHR2 = para_models[i,4],
               IHR3 = 0.001,
               IHR4 = para_models[i,3],
               IHR5 = 0.001,
               IHR6 = 0.0001)
  
  out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))
  #get total H
  total_H[i] = sum(out$H)
}

#Sensitivity analysis
sens<-sensitivity(x=total_H, numberf=6, make.plot=T, names = c("I2_0", "R0_2","IHR4","IHR2","vr","eff2"))

#Plot partial variances (Figure 2)
df.equilibrium <- data.frame(parameter=rbind("I2_0", "R0_2","IHR4","IHR2","vr","eff2"), value=sens)
windows(12,8)
p <- ggplot(df.equilibrium, aes(parameter, value))
p + geom_bar(stat="identity", fill=brewer.pal(3,"Set1")[2]) + 
  scale_x_discrete("Parameter", waiver(), c("I2_0", "R0_2","IHR4","IHR2","vr","eff2")) +
  ylab("Partial Variance") + theme_bw()+theme(axis.text=element_text(size = 16, colour = "black"), axis.title=element_text(size=20))



#Graph H for minimum - baseline - maximum of parameter range FAST analysis.
minimum_values = c(0.1*1.82e-5,2.0,1*0.001,1*0.0001,0.5*0.0045,0)
baseline_values = c(1.82e-5, 4.0, 0.001, 0.0001,0.0045,0.5)
maximum_values = c(10*1.82e-5,8.0,4*0.001,10*0.0001,2*0.0045,1.0)

times <- seq(0,180,by = 1)

#baseline:
init = c(S=(0.48-1.82e-5), V=0.2, I1V=0.0,I2V=0,I1S=0.02,I2S=1.82e-5,R1=0.3,R2=0,I2R1=0,I1R2=0)
parms = list(gamma = gamma,
             beta_SI1 = R0_delta*gamma,
             beta_SI2 = 2*R0_delta*gamma,
             beta_VI1 = R0_delta*gamma,
             beta_VI2 = 2*R0_delta*gamma,
             beta_R1I2 = 2*R0_delta*gamma,
             beta_R2I1 = 0.0,
             vr = 0.0045,
             eff1 = 0.5,
             eff2 = 0.0,
             IHR1 = 0.0001,
             IHR2 = 0.0001,
             IHR3 = 0.001,
             IHR4 = 2*0.001,
             IHR5 = 0.001,
             IHR6 = 0.0001)
rm(out)
out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

out_m = melt(out, id.vars = "time")

p_baseline <- ggplot(out_m %>% filter(variable=="H"), aes(x=time, y=value)) + 
  geom_line(size=1.1, colour="darkred") + 
  ylab("") + xlab("")

p_baseline
#"I2_0", "R0_2","IHR4","IHR2","vr","eff2"

#min_I2_0:
rm(init, parms, out, out_m)
init = c(S=(0.48-minimum_values[1]), V=0.2, I1V=0.0,I2V=0,I1S=0.02,I2S=minimum_values[1],R1=0.3,R2=0,I2R1=0,I1R2=0)
parms = list(gamma = gamma,
             beta_SI1 = R0_delta*gamma,
             beta_SI2 = 2*R0_delta*gamma,
             beta_VI1 = R0_delta*gamma,
             beta_VI2 = 2*R0_delta*gamma,
             beta_R1I2 = 2*R0_delta*gamma,
             beta_R2I1 = 0.0,
             vr = 0.0045,
             eff1 = 0.5,
             eff2 = 0.0,
             IHR1 = 0.0001,
             IHR2 = 0.0001,
             IHR3 = 0.001,
             IHR4 = 2*0.001,
             IHR5 = 0.001,
             IHR6 = 0.0001)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

out_m = melt(out, id.vars = "time")

p_min_I20 <- ggplot(out_m %>% filter(variable=="H"), aes(x=time, y=value)) + 
  geom_line(size=1.1, colour="darkred") + 
  ylab("Prop. Hosp.") + xlab("")
#max_I2_0:
rm(init, parms, out, out_m)
init = c(S=(0.48-maximum_values[1]), V=0.2, I1V=0.0,I2V=0,I1S=0.02,I2S=maximum_values[1],R1=0.3,R2=0,I2R1=0,I1R2=0)
parms = list(gamma = gamma,
             beta_SI1 = R0_delta*gamma,
             beta_SI2 = 2*R0_delta*gamma,
             beta_VI1 = R0_delta*gamma,
             beta_VI2 = 2*R0_delta*gamma,
             beta_R1I2 = 2*R0_delta*gamma,
             beta_R2I1 = 0.0,
             vr = 0.0045,
             eff1 = 0.5,
             eff2 = 0.0,
             IHR1 = 0.0001,
             IHR2 = 0.0001,
             IHR3 = 0.001,
             IHR4 = 2*0.001,
             IHR5 = 0.001,
             IHR6 = 0.0001)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

out_m = melt(out, id.vars = "time")

p_max_I20 <- ggplot(out_m %>% filter(variable=="H"), aes(x=time, y=value)) + 
  geom_line(size=1.1, colour="darkred") + 
  ylab("") + xlab("")


#p_min_I20 <- p_min_I20 + scale_y_continuous(limits = c(0, 5e-5))
#p_baseline <- p_baseline + scale_y_continuous(limits = c(0, 5e-5))
#plots combined:
P_I20 <- plot_grid(p_min_I20, p_baseline, p_max_I20, nrow=1)

P_I20
#"I2_0", "R0_2","IHR4","IHR2","vr","eff2"

#min_R2_0:
rm(init, parms, out, out_m)
init = c(S=(0.48-baseline_values[1]), V=0.2, I1V=0.0,I2V=0,I1S=0.02,I2S=baseline_values[1],R1=0.3,R2=0,I2R1=0,I1R2=0)
parms = list(gamma = gamma,
             beta_SI1 = R0_delta*gamma,
             beta_SI2 = minimum_values[2]*gamma,
             beta_VI1 = R0_delta*gamma,
             beta_VI2 = minimum_values[2]*gamma,
             beta_R1I2 = minimum_values[2]*gamma,
             beta_R2I1 = 0.0,
             vr = 0.0045,
             eff1 = 0.5,
             eff2 = 0.0,
             IHR1 = 0.0001,
             IHR2 = 0.0001,
             IHR3 = 0.001,
             IHR4 = 2*0.001,
             IHR5 = 0.001,
             IHR6 = 0.0001)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

out_m = melt(out, id.vars = "time")

p_min_R20 <- ggplot(out_m %>% filter(variable=="H"), aes(x=time, y=value)) + 
  geom_line(size=1.1, colour="darkred") + 
  ylab("Prop. Hosp") + xlab("")

#max_R2_0:
rm(init, parms, out, out_m)
init = c(S=(0.48-baseline_values[1]), V=0.2, I1V=0.0,I2V=0,I1S=0.02,I2S=baseline_values[1],R1=0.3,R2=0,I2R1=0,I1R2=0)
parms = list(gamma = gamma,
             beta_SI1 = R0_delta*gamma,
             beta_SI2 = maximum_values[2]*gamma,
             beta_VI1 = R0_delta*gamma,
             beta_VI2 = maximum_values[2]*gamma,
             beta_R1I2 = maximum_values[2]*gamma,
             beta_R2I1 = 0.0,
             vr = 0.0045,
             eff1 = 0.5,
             eff2 = 0.0,
             IHR1 = 0.0001,
             IHR2 = 0.0001,
             IHR3 = 0.001,
             IHR4 = 2*0.001,
             IHR5 = 0.001,
             IHR6 = 0.0001)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

out_m = melt(out, id.vars = "time")

p_max_R20 <- ggplot(out_m %>% filter(variable=="H"), aes(x=time, y=value)) + 
  geom_line(size=1.1, colour="darkred") + 
  ylab("") + xlab("")


#plots combined:
P_R20 <- plot_grid(p_min_R20, p_baseline, p_max_R20, nrow=1, align="t")

#"I2_0", "R0_2","IHR4","IHR2","vr","eff2"
#min_IHR4:
rm(init, parms, out, out_m)
init = c(S=(0.48-baseline_values[1]), V=0.2, I1V=0.0,I2V=0,I1S=0.02,I2S=baseline_values[1],R1=0.3,R2=0,I2R1=0,I1R2=0)
parms = list(gamma = gamma,
             beta_SI1 = R0_delta*gamma,
             beta_SI2 = baseline_values[2]*gamma,
             beta_VI1 = R0_delta*gamma,
             beta_VI2 = baseline_values[2]*gamma,
             beta_R1I2 = baseline_values[2]*gamma,
             beta_R2I1 = 0.0,
             vr = 0.0045,
             eff1 = 0.5,
             eff2 = 0.0,
             IHR1 = 0.0001,
             IHR2 = 0.0001,
             IHR3 = 0.001,
             IHR4 = minimum_values[3],
             IHR5 = 0.001,
             IHR6 = 0.0001)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

out_m = melt(out, id.vars = "time")

p_min_IHR4 <- ggplot(out_m %>% filter(variable=="H"), aes(x=time, y=value)) + 
  geom_line(size=1.1, colour="darkred") + 
  ylab("Prop. Hosp") + xlab("")

#max_IHR4:
rm(init, parms, out, out_m)
init = c(S=(0.48-baseline_values[1]), V=0.2, I1V=0.0,I2V=0,I1S=0.02,I2S=baseline_values[1],R1=0.3,R2=0,I2R1=0,I1R2=0)
parms = list(gamma = gamma,
             beta_SI1 = R0_delta*gamma,
             beta_SI2 = baseline_values[2]*gamma,
             beta_VI1 = R0_delta*gamma,
             beta_VI2 = baseline_values[2]*gamma,
             beta_R1I2 = baseline_values[2]*gamma,
             beta_R2I1 = 0.0,
             vr = 0.0045,
             eff1 = 0.5,
             eff2 = 0.0,
             IHR1 = 0.0001,
             IHR2 = 0.0001,
             IHR3 = 0.001,
             IHR4 = maximum_values[3],
             IHR5 = 0.001,
             IHR6 = 0.0001)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

out_m = melt(out, id.vars = "time")

p_max_IHR4 <- ggplot(out_m %>% filter(variable=="H"), aes(x=time, y=value)) + 
  geom_line(size=1.1, colour="darkred") + 
  ylab("") + xlab("")


#plots combined:
P_IHR4 <- plot_grid(p_min_IHR4, p_baseline, p_max_IHR4, nrow=1, align="t")

#"I2_0", "R0_2","IHR4","IHR2","vr","eff2"
#min_IHR2:
rm(init, parms, out, out_m)
init = c(S=(0.48-baseline_values[1]), V=0.2, I1V=0.0,I2V=0,I1S=0.02,I2S=baseline_values[1],R1=0.3,R2=0,I2R1=0,I1R2=0)
parms = list(gamma = gamma,
             beta_SI1 = R0_delta*gamma,
             beta_SI2 = baseline_values[2]*gamma,
             beta_VI1 = R0_delta*gamma,
             beta_VI2 = baseline_values[2]*gamma,
             beta_R1I2 = baseline_values[2]*gamma,
             beta_R2I1 = 0.0,
             vr = 0.0045,
             eff1 = 0.5,
             eff2 = 0.0,
             IHR1 = 0.0001,
             IHR2 = minimum_values[4],
             IHR3 = 0.001,
             IHR4 = baseline_values[3],
             IHR5 = 0.001,
             IHR6 = 0.0001)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

out_m = melt(out, id.vars = "time")

p_min_IHR2 <- ggplot(out_m %>% filter(variable=="H"), aes(x=time, y=value)) + 
  geom_line(size=1.1, colour="darkred") + 
  ylab("Prop. Hosp") + xlab("")

#max_IHR2:
rm(init, parms, out, out_m)
init = c(S=(0.48-baseline_values[1]), V=0.2, I1V=0.0,I2V=0,I1S=0.02,I2S=baseline_values[1],R1=0.3,R2=0,I2R1=0,I1R2=0)
parms = list(gamma = gamma,
             beta_SI1 = R0_delta*gamma,
             beta_SI2 = baseline_values[2]*gamma,
             beta_VI1 = R0_delta*gamma,
             beta_VI2 = baseline_values[2]*gamma,
             beta_R1I2 = baseline_values[2]*gamma,
             beta_R2I1 = 0.0,
             vr = 0.0045,
             eff1 = 0.5,
             eff2 = 0.0,
             IHR1 = 0.0001,
             IHR2 = maximum_values[4],
             IHR3 = 0.001,
             IHR4 = baseline_values[3],
             IHR5 = 0.001,
             IHR6 = 0.0001)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

out_m = melt(out, id.vars = "time")

p_max_IHR2 <- ggplot(out_m %>% filter(variable=="H"), aes(x=time, y=value)) + 
  geom_line(size=1.1, colour="darkred") + 
  ylab("") + xlab("")

#plots combined:
P_IHR2 <- plot_grid(p_min_IHR2, p_baseline,p_max_IHR2, nrow=1, align="t")


#"I2_0", "R0_2","IHR4","IHR2","vr","eff2"
#min_vr:
rm(init, parms, out, out_m)
init = c(S=(0.48-baseline_values[1]), V=0.2, I1V=0.0,I2V=0,I1S=0.02,I2S=baseline_values[1],R1=0.3,R2=0,I2R1=0,I1R2=0)
parms = list(gamma = gamma,
             beta_SI1 = R0_delta*gamma,
             beta_SI2 = baseline_values[2]*gamma,
             beta_VI1 = R0_delta*gamma,
             beta_VI2 = baseline_values[2]*gamma,
             beta_R1I2 = baseline_values[2]*gamma,
             beta_R2I1 = 0.0,
             vr = minimum_values[5],
             eff1 = 0.5,
             eff2 = 0.0,
             IHR1 = 0.0001,
             IHR2 = baseline_values[4],
             IHR3 = 0.001,
             IHR4 = baseline_values[3],
             IHR5 = 0.001,
             IHR6 = 0.0001)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

out_m = melt(out, id.vars = "time")

p_min_vr <- ggplot(out_m %>% filter(variable=="H"), aes(x=time, y=value)) + 
  geom_line(size=1.1, colour="darkred") + 
  ylab("Prop. Hosp") + xlab("")

#max_vr:
rm(init, parms, out, out_m)
init = c(S=(0.48-baseline_values[1]), V=0.2, I1V=0.0,I2V=0,I1S=0.02,I2S=baseline_values[1],R1=0.3,R2=0,I2R1=0,I1R2=0)
parms = list(gamma = gamma,
             beta_SI1 = R0_delta*gamma,
             beta_SI2 = baseline_values[2]*gamma,
             beta_VI1 = R0_delta*gamma,
             beta_VI2 = baseline_values[2]*gamma,
             beta_R1I2 = baseline_values[2]*gamma,
             beta_R2I1 = 0.0,
             vr = maximum_values[5],
             eff1 = 0.5,
             eff2 = 0.0,
             IHR1 = 0.0001,
             IHR2 = baseline_values[4],
             IHR3 = 0.001,
             IHR4 = baseline_values[3],
             IHR5 = 0.001,
             IHR6 = 0.0001)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

out_m = melt(out, id.vars = "time")

p_max_vr <- ggplot(out_m %>% filter(variable=="H"), aes(x=time, y=value)) + 
  geom_line(size=1.1, colour="darkred") + 
  ylab("") + xlab("")


#plots combined:
P_vr <- plot_grid(p_min_vr, p_baseline, p_max_vr, nrow=1, align="t")



#"I2_0", "R0_2","IHR4","IHR2","vr","eff2"
#min_eff2:
rm(init, parms, out, out_m)
init = c(S=(0.48-baseline_values[1]), V=0.2, I1V=0.0,I2V=0,I1S=0.02,I2S=baseline_values[1],R1=0.3,R2=0,I2R1=0,I1R2=0)
parms = list(gamma = gamma,
             beta_SI1 = R0_delta*gamma,
             beta_SI2 = baseline_values[2]*gamma,
             beta_VI1 = R0_delta*gamma,
             beta_VI2 = baseline_values[2]*gamma,
             beta_R1I2 = baseline_values[2]*gamma,
             beta_R2I1 = 0.0,
             vr = baseline_values[5],
             eff1 = 0.5,
             eff2 = minimum_values[6],
             IHR1 = 0.0001,
             IHR2 = baseline_values[4],
             IHR3 = 0.001,
             IHR4 = baseline_values[3],
             IHR5 = 0.001,
             IHR6 = 0.0001)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

out_m = melt(out, id.vars = "time")

p_min_eff2 <- ggplot(out_m %>% filter(variable=="H"), aes(x=time, y=value)) + 
  geom_line(size=1.1, colour="darkred") + 
  ylab("Prop. Hosp") + xlab("Time")

#max_eff2:
rm(init, parms, out, out_m)
init = c(S=(0.48-baseline_values[1]), V=0.2, I1V=0.0,I2V=0,I1S=0.02,I2S=baseline_values[1],R1=0.3,R2=0,I2R1=0,I1R2=0)
parms = list(gamma = gamma,
             beta_SI1 = R0_delta*gamma,
             beta_SI2 = baseline_values[2]*gamma,
             beta_VI1 = R0_delta*gamma,
             beta_VI2 = baseline_values[2]*gamma,
             beta_R1I2 = baseline_values[2]*gamma,
             beta_R2I1 = 0.0,
             vr = baseline_values[5],
             eff1 = 0.5,
             eff2 = maximum_values[6],
             IHR1 = 0.0001,
             IHR2 = baseline_values[4],
             IHR3 = 0.001,
             IHR4 = baseline_values[3],
             IHR5 = 0.001,
             IHR6 = 0.0001)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

out_m = melt(out, id.vars = "time")

p_max_eff2 <- ggplot(out_m %>% filter(variable=="H"), aes(x=time, y=value)) + 
  geom_line(size=1.1, colour="darkred") + 
  ylab("") + xlab("")


#plots combined:
P_eff2 <- plot_grid(p_min_eff2, p_baseline, p_max_eff2, nrow=1, align="t")


###############
#Combine all plots
###############
#"I2_0", "R0_2","IHR4","IHR2","vr","eff2"
p_all <- plot_grid(P_I20, P_R20, P_IHR4, P_IHR2, P_vr, P_eff2, nrow=6, align="l", axis="l")
p_all
p_all + ggtext()

###############################
# R0(x) variations
################################
#baseline

init = c(S=(0.48-baseline_values[1]), V=0.2, I1V=0.0,I2V=0,I1S=0.02,I2S=baseline_values[1],R1=0.3,R2=0,I2R1=0,I1R2=0)
parms = list(gamma = gamma,
             beta_SI1 = R0_delta*gamma,
             beta_SI2 = baseline_values[2]*gamma,
             beta_VI1 = R0_delta*gamma,
             beta_VI2 = baseline_values[2]*gamma,
             beta_R1I2 = baseline_values[2]*gamma,
             beta_R2I1 = 0.0,
             vr = 0.0045,
             eff1 = 0.5,
             eff2 = 0.5,
             IHR1 = 0.0001,
             IHR2 = 0.0001,
             IHR3 = 0.001,
             IHR4 = 2*0.001,
             IHR5 = 0.001,
             IHR6 = 0.0001)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

out_m = melt(out, id.vars = "time")

p_baseline_R20 <- ggplot(out_m %>% filter(variable=="H"), aes(x=time, y=value)) + 
  geom_line(size=1.1, colour="darkred") + 
  ylab("Proportion of population")


#min_R2_0:
init = c(S=(0.48-baseline_values[1]), V=0.2, I1V=0.0,I2V=0,I1S=0.02,I2S=baseline_values[1],R1=0.3,R2=0,I2R1=0,I1R2=0)
parms = list(gamma = gamma,
             beta_SI1 = 1.0*gamma,
             beta_SI2 = 2.0*gamma,
             beta_VI1 = 1.0*gamma,
             beta_VI2 = 2.0*gamma,
             beta_R1I2 = 2.0*gamma,
             beta_R2I1 = 0.0,
             vr = 0.0045,
             eff1 = 0.5,
             eff2 = 0.5,
             IHR1 = 0.0001,
             IHR2 = 0.0001,
             IHR3 = 0.001,
             IHR4 = 2*0.001,
             IHR5 = 0.001,
             IHR6 = 0.0001)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

out_m = melt(out, id.vars = "time")

p_min_R20 <- ggplot(out_m %>% filter(variable=="H"), aes(x=time, y=value)) + 
  geom_line(size=1.1, colour="darkred") + 
  ylab("Proportion of population")

#max_I2_0:
init = c(S=(0.48-baseline_values[1]), V=0.2, I1V=0.0,I2V=0,I1S=0.02,I2S=baseline_values[1],R1=0.3,R2=0,I2R1=0,I1R2=0)
parms = list(gamma = gamma,
             beta_SI1 = 4.0*gamma,
             beta_SI2 = 8.0*gamma,
             beta_VI1 = 4.0*gamma,
             beta_VI2 = 8.0*gamma,
             beta_R1I2 = 8.0*gamma,
             beta_R2I1 = 0.0,
             vr = 0.0045,
             eff1 = 0.5,
             eff2 = 0.5,
             IHR1 = 0.0001,
             IHR2 = 0.0001,
             IHR3 = 0.001,
             IHR4 = 2*0.001,
             IHR5 = 0.001,
             IHR6 = 0.0001)

out <- data.frame(ode(y = init, func = SIR, times = times, parms = parms))

out_m = melt(out, id.vars = "time")

p_max_R20 <- ggplot(out_m %>% filter(variable=="H"), aes(x=time, y=value)) + 
  geom_line(size=1.1, colour="darkred") + 
  ylab("Proportion of population")


#plots combined:
P_R20 <- plot_grid(p_min_R20, p_baseline, p_max_R20, nrow=1, align="t")
P_R20
