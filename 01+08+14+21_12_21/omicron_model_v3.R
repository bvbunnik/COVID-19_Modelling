library("deSolve"); library("ggpubr"); library("tidyverse"); library("RColorBrewer")
library("gridExtra"); library("grid");library("cowplot"); library("reshape2")
library("parallel"); library("doParallel")
rm(list=ls())

#Function for the generation time - a function of R0 and the doubling time
GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

GenTime(2,4*1.4)

GT = 5
Re_delta = 1.4
gamma = 1/GT
Re_omicron = 4*Re_delta



vaccination_function <- function(time, end_time, vacc_rate){
  if (time<end_time){
    return(vacc_rate)
  } else {
    return(0)
  }
}



VBIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    vr1 = vaccination_function(time, 60, vr)
    
    if ((V-vr1)<0){vr1=0}
    dV = -vr*V - beta_VI1*V*(I1V+I1B+I1R2+I1R2B) - beta_VI2*V*(I2V+I2B+I2R1+I2R1B)
    dB = vr*V - (1-eff1)*beta_VI1*B*(I1V+I1B+I1R2+I1R2B) - (1-eff2)*beta_VI2*B*(I2V+I2B+I2R1+I2R1B)
    dI1V = beta_VI1*V*(I1V+I1B+I1R2+I1R2B) - gamma*I1V
    dI2V = beta_VI2*V*(I2V+I2B+I2R1+I2R1B) - gamma*I2V
    dI1B = (1-eff1)*beta_VI1*B*(I1V+I1B+I1R2+I1R2B) - gamma*I1B
    dI2B = (1-eff2)*beta_VI2*B*(I2V+I2B+I2R1+I2R1B) - gamma*I2B
    dR1 = gamma*I1V - vr*R1 - (1-eff4)*beta_VI2*R1*(I2V+I2B+I2R1+I2R1B)
    dR2 = gamma*I2V - vr*R2 - (1-eff3)*beta_VI1*R2*(I1V+I1B+I1R2+I1R2B)
    dI2R1 = (1-eff4)*beta_VI2*R1*(I2V+I2B+I2R1+I2R1B) - gamma*I2R1
    dI1R2 = (1-eff3)*beta_VI1*R2*(I1V+I1B+I1R2+I1R2B) - gamma*I1R2
    dR1B = gamma*I1B + vr*R1 - (1-eff6)*beta_VI2*R1B*(I2V+I2B+I2R1+I2R1B)
    dR2B = gamma*I2B + vr*R2 - (1-eff5)*beta_VI1*R2B*(I1V+I1B+I1R2+I1R2B)
    dI1R2B = (1-eff5)*beta_VI1*R2B*(I1V+I1B+I1R2+I1R2B) - gamma*I1R2B
    dI2R1B = (1-eff6)*beta_VI2*R1B*(I2V+I2B+I2R1+I2R1B) - gamma*I2R1B
    dR1R2 = gamma*I2R1 + gamma*I1R2 - vr*R1R2
    dRall = vr*R1R2 + gamma*I1R2B + gamma*I2R1B
  
    dcumB = vr*V + vr*R1 + vr*R2 + vr*R1R2
    totalI = I1V+I2V+I1B+I2B+I2R1+I1R2+I1R2B+I2R1B
    totalI1 = I1V+I1B+I1R2+I1R2B
    totalI2 = I2V+I2B+I2R1+I2R1B
    
    H = IHR1*I1V + IHR2*I2V + he1*IHR1*I1B + he2*IHR2*I2B + he4*IHR2*I2R1 + he3*IHR1*I1R2 + he5*IHR1*I1R2B + he6*IHR2*I2R1B
    dcumH = H
    return(list(c(dV,dB,dI1V,dI2V,dI1B,dI2B,dR1,dR2,dI2R1,dI1R2,dR1B,dR2B,dI1R2B,dI2R1B,dR1R2,dRall,dcumH, dcumB), 
                "H"=H,"total_I"=totalI,"total_I1"=totalI1,"total_I2"=totalI2))
  })
}

times <- seq(0,180,by = 0.1)
  
  
init = c(V = (0.455 - 0.01130435 - 7.494118e-06),#- 1.028696e-05
           B = (0.245-0.008695652 - 5.764706e-06),#- 7.913043e-06
           I1V = 0.01130435, #0.02*0.455/(0.35+0.455),
           I2V =  7.494118e-06,#1.82E-05*0.455/(0.3+0.35+0.455),
           I1B = 0.008695652, #0.02*0.35/(0.35+0.455),
           I2B = 5.764706e-06, #1.82E-05*0.35/(0.3+0.35+0.455)
           R1 = 0.3-0.105-4.941176e-06,
           R2 = 0,
           I2R1 = 4.941176e-06, #1.82E-05*0.3/(0.3+0.35+0.455)
           I1R2 = 0,
           R1B = 0.105,
           R2B = 0,
           I1R2B = 0,
           I2R1B = 0,
           R1R2 = 0,
           Rall = 0,
           cumH = 0,
           cumB = 0.245+0.105
        )


parms = list(gamma = gamma,
               vr = 0.025,
               beta_VI1 = Re_delta*gamma,
               beta_VI2 = 0.4*Re_omicron*gamma,#Re_omicron*gamma,
               eff1 = 0,
               eff2 = 0,
               eff3 = 0,
               eff4 = 0,
               eff5 = 1-(1-0.0)*(1-0.0), #1-1(1-eff1)*(1-eff3)
               eff6 = 1-(1-0.0)*(1-0.0), #1-1(1-eff2)*(1-eff4)
               IHR1 = 0.0012,
               IHR2 = 1.0*0.0012,
               he1=0.01,
               he2=1*0.01,
               he3=1,
               he4=1,
               he5=1*0.01, #he1*he3
               he6=1*0.01) #he2*he4
  
out <- data.frame(ode(y = init, func = VBIR, times = times, parms = parms))

max(out$H)/out$H[1]
parms$beta_VI2/parms$beta_VI1

total_I2_t0 = out$total_I2[1]
target = 2*total_I2_t0

t_d = which(abs(out$total_I2 - target) == min(abs(out$total_I2 - target)))

out$time[t_d]

out_m = melt(out, id.vars = "time")
  
#plot all comps on free scales: 
#ggplot(out_m, aes(x=time, y=value)) + 
#    geom_line(size=1.0, colour="darkred") + 
#    facet_wrap(~variable, scales = "free_y") +
#    ylab("Proportion of population")

#make V,B,R1,R2 same scale:
max1=out_m %>% filter(variable %in% c("V","B","R1","R2")) %>% summarise(max1=max(value)) %>% as.numeric

#make V,B,R1,R2 same scale:
max2=out_m %>% filter(variable %in% c("I1V","I2V","I1B","I2B","I2R1", "I1R2","I1R2B","I2R1B")) %>% summarise(max2=max(value)) %>% as.numeric

#make plots for all comps:


myplots <- vector('list', ncol(out))
for (i in seq_along(out)) {
  message(i)
  myplots[[i]] <- local({
    i <- i
    p1 <- ggplot(out, aes(x = out[[1]], y=out[[i]])) +
      geom_line(colour="darkred", size=1.1) + xlab("") +ylab("") + ggtitle(colnames(out)[i]) + theme_bw() +
      theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
            axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
            legend.position="bottom")
  })
}

for(i in c(2,3,8,9)){
  myplots[[i]] <- myplots[[i]] + ylim(c(0,max1))
}

for(i in c(4,5,6,7,10,11,14,15 )){
  myplots[[i]] <- myplots[[i]] + ylim(c(0,max2))
}

myplots[[18]] <- ggplot(out, aes(x = out[[1]]+7, y=out[[18]])) +
  geom_line(colour="darkred", size=1.1) + xlab("") +ylab("") + ggtitle(colnames(out)[18]) + theme_bw() + 
  theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
        axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position="bottom")
myplots[[20]] <- ggplot(out, aes(x = out[[1]]+7, y=out[[20]])) +
  geom_line(colour="darkred", size=1.1) + xlab("") +ylab("") + ggtitle(colnames(out)[20]) + theme_bw() + 
  theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
        axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position="bottom")

grid <- plot_grid(myplots[[2]],myplots[[3]],myplots[[4]],myplots[[5]],myplots[[6]],
          myplots[[7]],myplots[[8]],myplots[[9]],myplots[[10]],myplots[[11]],
          myplots[[12]],myplots[[13]],myplots[[14]],myplots[[15]],myplots[[16]],
          myplots[[17]],myplots[[18]],myplots[[19]],myplots[[20]],myplots[[21]],
          myplots[[22]], myplots[[23]], nnrow = 5, axis="lbtr")

y.grob <- textGrob("Prop. of pop.", 
                   gp=gpar(fontface="bold", fontsize=15), rot=90)
x.grob <- textGrob("Time (days)", 
                   gp=gpar(fontface="bold", fontsize=15))

grid.arrange(arrangeGrob(grid, left = y.grob, bottom = x.grob))

ggsave(grid.arrange(arrangeGrob(grid, left = y.grob, bottom = x.grob)), filename = "baseline_2strain_Re_omicron0.4x.png", dpi = 300, 
       type = "cairo", width = 25, height = 20, units = "in", path = "figures/")



total_I2_t0 = out$total_I2[1]
target = 2*total_I2_t0

which(abs(out$total_I2 - target) == min(abs(out$total_I2 - target)))

#parameters to vary:
#I2(0), vr, beta_VI1, beta_VI2, eff1, eff2, eff3, eff4, IHR2, he1, he2, he3, he4
pars_min = c(20/5.5e6,0.025,0.07,2*0.07,0,0,0,0,0.1*0.0012,0.005,0.1*0.005,0.1,0.1)
pars_max = c(500/5.5e6,0.05,0.28,6*0.28,1,1,1,1,2*0.0012,0.02,2*0.02,1,1)
pars_names = c("I2(0)","vr","beta_VI1","beta_VI2","eff1","eff2","eff3","eff4","IHR2","he1","he2","he3","he4")
library("fast")

para_models = fast_parameters(minimum = pars_min, maximum = pars_max, names = pars_names, factor=9)
i=1
total_H = vector(mode = "numeric", length=nrow(para_models))

parallel::detectCores()
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(my.cluster)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()
#how many workers are available? (optional)
foreach::getDoParWorkers()

total_H = vector(mode = "numeric", length=nrow(para_models))
for(i in 1:nrow(para_models)) 
  {
  
  init1 = c(V = (0.455 - 0.01130435 - para_models[i,1]*0.455/(0.3+0.35+0.455)),#- 1.028696e-05
           B = (0.245-0.008695652 - 5.764706e-06),#- 7.913043e-06
           I1V = 0.01130435, #0.02*0.455/(0.35+0.455),
           I2V =  para_models[i,1]*0.455/(0.3+0.35+0.455),#1.82E-05*0.455/(0.3+0.35+0.455),
           I1B = 0.008695652, #0.02*0.35/(0.35+0.455),
           I2B = para_models[i,1]*0.35/(0.3+0.35+0.455), #1.82E-05*0.35/(0.3+0.35+0.455)
           R1 = 0.3-0.105-4.941176e-06,
           R2 = 0,
           I2R1 = para_models[i,1]*0.3/(0.3+0.35+0.455), #1.82E-05*0.3/(0.3+0.35+0.455)
           I1R2 = 0,
           R1B = 0.105,
           R2B = 0,
           I1R2B = 0,
           I2R1B = 0,
           R1R2 = 0,
           Rall = 0,
           cumH = 0,
           cumB = 0.245+0.105
  )
  #I2(0), vr, beta_VI1, beta_VI2, eff1, eff2, eff3, eff4, IHR2, he1, he2, he3, he4
  parms1 = list(gamma = gamma,
               vr = para_models[i,2],
               beta_VI1 = para_models[i,3],
               beta_VI2 = para_models[i,4],#Re_omicron*gamma,
               eff1 = para_models[i,5],
               eff2 = para_models[i,6],
               eff3 = para_models[i,7],
               eff4 = para_models[i,8],
               eff5 = 1-(1-para_models[i,5])*(1-para_models[i,7]), #1-1(1-eff1)*(1-eff3)
               eff6 = 1-(1-para_models[i,6])*(1-para_models[i,8]), #1-1(1-eff2)*(1-eff4)
               IHR1 = 0.0012,
               IHR2 = para_models[i,9],
               he1=para_models[i,10],
               he2=para_models[i,11],
               he3=para_models[i,12],
               he4=para_models[i,13],
               he5=para_models[i,10]*para_models[i,12], #he1*he3
               he6=para_models[i,11]*para_models[i,13]) #he2*he4
  
  out1 <- data.frame(ode(y = init1, func = VBIR, times = times, parms = parms1))
  #get total H
  total_H[i]=out1$cumH[1801]
  rm(init1, parms1, out1)
  }

parallel::stopCluster(cl = my.cluster)

write_csv(as.data.frame(total_H), "data/FAST_sens_totalH.csv")
#Sensitivity analysis
sens<-sensitivity(x=total_H, numberf=13, make.plot=T, names = pars_names)

#Plot partial variances (Figure 2)
df.equilibrium <- data.frame(parameter=rbind("I2(0)","vr","beta_VI1","beta_VI2","eff1","eff2","eff3","eff4","IHR2","he1","he2","he3","he4"), value=sens)
#windows(12,8)
p <- ggplot(df.equilibrium, aes(parameter, value))
p <- p + geom_bar(stat="identity", fill=brewer.pal(3,"Set1")[2]) + 
  scale_x_discrete("Parameter", waiver(), c("I2(0)","vr",expression(paste(beta[VI1])),expression(paste(beta[VI2])),"eff1","eff2","eff3","eff4","IHR2","he1","he2","he3","he4")) +
  ylab("Partial Variance") + theme_bw()+theme(axis.text=element_text(size = 16, colour = "black"), axis.title=element_text(size=20))
p

ggsave(p, filename = "FAST_sens.png", dpi = 300, type = "cairo", width = 16, height = 9, units = "in",
       path = "figures/")


# Sensitivity Analysis ----------------------------------------------------

#The list below are the explored values for each parameter (a vector for each parameter nested in a list)

sens_parms <- list("vr" = c(0.025, 0.05), 
                   "beta_VI1" = c(0.28, 0.07), 
                   "beta_VI2" = c(0.28*2, 0.28*4, 0.28*6), 
                   "eff1"= c(0, 0.5, 1), 
                   "eff2" = c(0, 0.5, 1), 
                   "eff3" = c(0, 0.5, 1), 
                   "eff4" = c(0, 0.5, 1), 
                   "IHR2" = c(0.0012*0.1, 0.0012, 0.0012*2),
                   "he1" = c(0.005, 0.01, 0.02), 
                   "he2" = c(0.01*0.1, 0.01, 0.01*2), 
                   "he3" = c(0.1, 0.5, 1), 
                   "he4" = c(0.1, 0.5, 1))


times <- seq(0, 180, by = 0.1)
gamma = 0.2
Reff_delta = 1.4
Reff_omicron = 1.4 * 4

init = c(V = (0.455 - 0.01130435 - 7.494118e-06),#- 1.028696e-05
         B = (0.245-0.008695652 - 5.764706e-06),#- 7.913043e-06
         I1V = 0.01130435, #0.02*0.455/(0.35+0.455),
         I2V =  7.494118e-06,#1.82E-05*0.455/(0.3+0.35+0.455),
         I1B = 0.008695652, #0.02*0.35/(0.35+0.455),
         I2B = 5.764706e-06, #1.82E-05*0.35/(0.3+0.35+0.455)
         R1 = 0.3-0.105-4.941176e-06,
         R2 = 0,
         I2R1 = 4.941176e-06, #1.82E-05*0.3/(0.3+0.35+0.455)
         I1R2 = 0,
         R1B = 0.105,
         R2B = 0,
         I1R2B = 0,
         I2R1B = 0,
         R1R2 = 0,
         Rall = 0,
         cumH = 0,
         cumB = 0.245+0.105
)

store_dat <- list()

for(i in 1:length(sens_parms)) {
  alter_vec <- sens_parms[[i]]
  parm_list <- list()
  for(z in 1:length(alter_vec)) {
    parms1 <- parms
    parms1[names(sens_parms)[i]] <- alter_vec[[z]]
    out <- data.frame(ode(y = init, func = VBIR, times = times, parms = parms1))
  
    
    parm_list[[z]] <- data.frame("time" = out["time"],
                                 "H" = out["H"],
                                 "value" = as.character(alter_vec[[z]]),
                                 "group" = names(sens_parms)[i])
  }
  store_dat[[i]] <- parm_list
}

#Plotting Sensitivity Analyses

plot_list <- list()

for(i in 1:length(store_dat)) {
  plotdat <- do.call(rbind.data.frame, store_dat[[i]])
  m_plotdat <- melt(plotdat, id.vars = c("time", "H"), measure.vars = "value")
  
  plot_list[[i]] <- ggplot(m_plotdat, aes(x = time, y= H, col = value)) + geom_line(size = 1.5) + theme_bw() + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, max(m_plotdat["H"])*1.2)) + 
    labs(title= paste0("Parameter = ",  names(sens_parms)[i]), 
         y="Incidence of Hospitalisations", x = "Time (days)", colour = "") + 
    theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
          axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
          legend.position="bottom")
}

p2 <- ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],
                plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]],
                plot_list[[9]], plot_list[[10]], plot_list[[11]], plot_list[[12]],
                ncol = 4, nrow=3)

p2
ggsave(p2, filename = "univariate_sens1.png", dpi = 300, type = "cairo", width = 20, height = 20, units = "in",
       path = "figures/")

#special case for I2V(0):
#lower value:
init = c(V = (0.455 - 0.01130435 - (20/100)*7.494118e-06),#- 1.028696e-05
         B = (0.245-0.008695652 - (20/100)*5.764706e-06),#- 7.913043e-06
         I1V = 0.01130435, #0.02*0.455/(0.35+0.455),
         I2V =  (20/100)*7.494118e-06,#1.82E-05*0.455/(0.3+0.35+0.455),
         I1B = 0.008695652, #0.02*0.35/(0.35+0.455),
         I2B = (20/100)*5.764706e-06, #1.82E-05*0.35/(0.3+0.35+0.455)
         R1 = 0.3-0.105-(20/100)*4.941176e-06,
         R2 = 0,
         I2R1 = (20/100)*4.941176e-06, #1.82E-05*0.3/(0.3+0.35+0.455)
         I1R2 = 0,
         R1B = 0.105,
         R2B = 0,
         I1R2B = 0,
         I2R1B = 0,
         R1R2 = 0,
         Rall = 0,
         cumH = 0,
         cumB = 0.245+0.105
)

parms1 <- parms
out <- data.frame(ode(y = init, func = VBIR, times = times, parms = parms1))


parm_list[[1]] <- data.frame("time" = out["time"],
                             "H" = out["H"],
                             "value" = "20",
                             "group" = "I2(0)")

init = c(V = (0.455 - 0.01130435 - (100/100)*7.494118e-06),#- 1.028696e-05
         B = (0.245-0.008695652 - (100/100)*5.764706e-06),#- 7.913043e-06
         I1V = 0.01130435, #0.02*0.455/(0.35+0.455),
         I2V =  (100/100)*7.494118e-06,#1.82E-05*0.455/(0.3+0.35+0.455),
         I1B = 0.008695652, #0.02*0.35/(0.35+0.455),
         I2B = (100/100)*5.764706e-06, #1.82E-05*0.35/(0.3+0.35+0.455)
         R1 = 0.3-0.105-(100/100)*4.941176e-06,
         R2 = 0,
         I2R1 = (100/100)*4.941176e-06, #1.82E-05*0.3/(0.3+0.35+0.455)
         I1R2 = 0,
         R1B = 0.105,
         R2B = 0,
         I1R2B = 0,
         I2R1B = 0,
         R1R2 = 0,
         Rall = 0,
         cumH = 0,
         cumB = 0.245+0.105
)
parms1 <- parms
out <- data.frame(ode(y = init, func = VBIR, times = times, parms = parms1))


parm_list[[2]] <- data.frame("time" = out["time"],
                             "H" = out["H"],
                             "value" = "100",
                             "group" = "I2(0)")

init = c(V = (0.455 - 0.01130435 - (500/100)*7.494118e-06),#- 1.028696e-05
         B = (0.245-0.008695652 - (500/100)*5.764706e-06),#- 7.913043e-06
         I1V = 0.01130435, #0.02*0.455/(0.35+0.455),
         I2V =  (500/100)*7.494118e-06,#1.82E-05*0.455/(0.3+0.35+0.455),
         I1B = 0.008695652, #0.02*0.35/(0.35+0.455),
         I2B = (500/100)*5.764706e-06, #1.82E-05*0.35/(0.3+0.35+0.455)
         R1 = 0.3-0.105-(500/100)*4.941176e-06,
         R2 = 0,
         I2R1 = (500/100)*4.941176e-06, #1.82E-05*0.3/(0.3+0.35+0.455)
         I1R2 = 0,
         R1B = 0.105,
         R2B = 0,
         I1R2B = 0,
         I2R1B = 0,
         R1R2 = 0,
         Rall = 0,
         cumH = 0,
         cumB = 0.245+0.105
)
parms1 <- parms
out <- data.frame(ode(y = init, func = VBIR, times = times, parms = parms1))


parm_list[[3]] <- data.frame("time" = out["time"],
                             "H" = out["H"],
                             "value" = "500",
                             "group" = "I2(0)")


store_dat[[1]] <- parm_list

plot_list = list()
for(i in 1:length(store_dat)) {
  plotdat <- do.call(rbind.data.frame, store_dat[[i]])
  m_plotdat <- melt(plotdat, id.vars = c("time", "H"), measure.vars = "value")
  
  plot_list[[i]] <- ggplot(m_plotdat, aes(x = time, y= H, col = value)) + geom_line(size = 1.5) + theme_bw() + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, max(m_plotdat["H"])*1.2)) + 
    labs(title="Parameter = I2(0)", y="Incidence of Hospitalisations", x = "Time (days)", colour = "") + 
    theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
          axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
          legend.position="bottom")
}
plot_list[[1]]


#vary all betas similar:

sens_parms <- list("beta_VI1" = c(1*Re_delta*gamma, 0.75*Re_delta*gamma, 0.5*Re_delta*gamma, 0.25*Re_delta*gamma), 
                   "beta_VI2" = c(1*Re_omicron*gamma, 0.75*Re_omicron*gamma, 0.5*Re_omicron*gamma, 0.25*Re_omicron*gamma))
values = c("1", "0.75", "0.5", "0.25")

times <- seq(0, 180, by = 0.1)
gamma = 0.2
Reff_delta = 1.4
Reff_omicron = 1.4 * 4

init = c(V = (0.455 - 0.01130435 - 7.494118e-06),#- 1.028696e-05
         B = (0.245-0.008695652 - 5.764706e-06),#- 7.913043e-06
         I1V = 0.01130435, #0.02*0.455/(0.35+0.455),
         I2V =  7.494118e-06,#1.82E-05*0.455/(0.3+0.35+0.455),
         I1B = 0.008695652, #0.02*0.35/(0.35+0.455),
         I2B = 5.764706e-06, #1.82E-05*0.35/(0.3+0.35+0.455)
         R1 = 0.3-0.105-4.941176e-06,
         R2 = 0,
         I2R1 = 4.941176e-06, #1.82E-05*0.3/(0.3+0.35+0.455)
         I1R2 = 0,
         R1B = 0.105,
         R2B = 0,
         I1R2B = 0,
         I2R1B = 0,
         R1R2 = 0,
         Rall = 0,
         cumH = 0,
         cumB = 0.245+0.105
)
i=1
store_dat <- list()
parm_list = list()
for(i in 1:4) {
  beta1 <- sens_parms[[1]][i]
  beta2 <- sens_parms[[2]][i]
  
  parms1 <- parms
  parms1["beta_VI1"] <- beta1
  parms1["beta_VI2"] <- beta2
  
  out <- data.frame(ode(y = init, func = VBIR, times = times, parms = parms1))
    
  parm_list[[i]] <- data.frame("time" = out["time"],
                                 "H" = out["H"],
                                 "value" = values[i])
}

plot_data <- rbind(parm_list[[1]], parm_list[[2]], parm_list[[3]], parm_list[[4]])

p <- ggplot(plot_data, aes(x=time, y=H, col=value)) + geom_line(size = 1.5) + theme_bw() + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, max(plot_data["H"])*1.2)) + 
  labs(y="Incidence of Hospitalisations", x = "Time (days)", colour = "") + 
  theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
        axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position="bottom")

ggsave(p, filename = "betas_plot.png", dpi = 300, type = "cairo", width = 9, height = 6, units = "in",
       path = "figures/")

#Sinlge optimsitic run:

parms = list(gamma = gamma,
             vr = 0.025,
             beta_VI1 = 0.28,
             beta_VI2 = 0.56,#Re_omicron*gamma,
             eff1 = 0.5,
             eff2 = 0.5,
             eff3 = 0.5,
             eff4 = 0.5,
             eff5 = 0.75, #1-1(1-eff1)*(1-eff3)
             eff6 = 0.75, #1-1(1-eff2)*(1-eff4)
             IHR1 = 0.0012,
             IHR2 = 0.0006,
             he1=0.01,
             he2=0.01,
             he3=0.1,
             he4=0.1,
             he5=0.001, #he1*he3
             he6=0.001) #he2*he4

out <- data.frame(ode(y = init, func = VBIR, times = times, parms = parms))

total_I2_t0 = out$total_I2[1]
target = 2*total_I2_t0

which(abs(out$total_I2 - target) == min(abs(out$total_I2 - target)))


out_m = melt(out, id.vars = "time")

#make V,B,R1,R2 same scale:
max1=out_m %>% filter(variable %in% c("V","B","R1","R2")) %>% summarise(max1=max(value)) %>% as.numeric

#make V,B,R1,R2 same scale:
max2=out_m %>% filter(variable %in% c("I1V","I2V","I1B","I2B","I2R1", "I1R2","I1R2B","I2R1B")) %>% summarise(max2=max(value)) %>% as.numeric

#make plots for all comps:


myplots <- vector('list', ncol(out))
for (i in seq_along(out)) {
  message(i)
  myplots[[i]] <- local({
    i <- i
    p1 <- ggplot(out, aes(x = out[[1]], y=out[[i]])) +
      geom_line(colour="darkred", size=1.1) + xlab("") +ylab("") + ggtitle(colnames(out)[i]) + theme_bw() +
      theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
            axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
            legend.position="bottom")
  })
}

for(i in c(2,3,8,9)){
  myplots[[i]] <- myplots[[i]] + ylim(c(0,max1))
}

for(i in c(4,5,6,7,10,11,14,15 )){
  myplots[[i]] <- myplots[[i]] + ylim(c(0,max2))
}

myplots[[18]] <- ggplot(out, aes(x = out[[1]]+7, y=out[[18]])) +
  geom_line(colour="darkred", size=1.1) + xlab("") +ylab("") + ggtitle(colnames(out)[18]) + theme_bw() + 
  theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
        axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position="bottom")
myplots[[20]] <- ggplot(out, aes(x = out[[1]]+7, y=out[[20]])) +
  geom_line(colour="darkred", size=1.1) + xlab("") +ylab("") + ggtitle(colnames(out)[20]) + theme_bw() + 
  theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
        axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position="bottom")

grid <- plot_grid(myplots[[2]],myplots[[3]],myplots[[4]],myplots[[5]],myplots[[6]],
                  myplots[[7]],myplots[[8]],myplots[[9]],myplots[[10]],myplots[[11]],
                  myplots[[12]],myplots[[13]],myplots[[14]],myplots[[15]],myplots[[16]],
                  myplots[[17]],myplots[[18]],myplots[[19]],myplots[[20]],myplots[[21]],
                  myplots[[22]], myplots[[23]], nnrow = 5, axis="lbtr")

y.grob <- textGrob("Prop. of pop.", 
                   gp=gpar(fontface="bold", fontsize=15), rot=90)
x.grob <- textGrob("Time (days)", 
                   gp=gpar(fontface="bold", fontsize=15))

grid.arrange(arrangeGrob(grid, left = y.grob, bottom = x.grob))

ggsave(grid.arrange(arrangeGrob(grid, left = y.grob, bottom = x.grob)), filename = "optmistic_parameters.png", dpi = 300, 
       type = "cairo", width = 25, height = 20, units = "in", path = "figures/")
