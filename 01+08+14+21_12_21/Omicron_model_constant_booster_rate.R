library("deSolve"); library("ggpubr"); library("tidyverse"); library("RColorBrewer")
library("gridExtra"); library("grid");library("cowplot"); library("reshape2")
library("parallel"); library("doParallel"); library("sigmoid")
rm(list=ls())


# Functions and init ------------------------------------------------------

GenTime <- function(T2, R0) {
  G = T2 * ((R0-1)/log(2))
  return(G)
}

vaccination_function <- function(time, end_time, vacc_rate){
  if (time<end_time){
    return(vacc_rate)
  } else {
    return(0)
  }
}

trajectory_plots <- function(data, plot=FALSE, save=TRUE, name="")
{
  warn<-options(warn=-1)
  out_m = melt(data, id.vars = "time")
  
  #make V,B,R1,R2 same scale:
  max1=out_m %>% filter(variable %in% c("V","B","R1","R2")) %>% summarise(max1=max(value)) %>% as.numeric
  
  #make V,B,R1,R2 same scale:
  max2=out_m %>% filter(variable %in% c("I1V","I2V","I1B","I2B","I2R1", "I1R2","I1R2B","I2R1B")) %>% summarise(max2=max(value)) %>% as.numeric
  
  #make plots for all comps:
  myplots <- vector('list', ncol(data))
  for (i in seq_along(data)) {
    myplots[[i]] <- local({
      i <- i
      p1 <- ggplot(data[,c(1,i)], aes(x = data[[1]], y=data[[i]])) +
        geom_line(colour="darkred", size=1.1) + xlab("") +ylab("") + ggtitle(colnames(data)[i]) + theme_bw() +
        theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
              axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
              legend.position="bottom")
    })
  }
  
  for(i in c(2,3,8,9,16,17)){
    myplots[[i]] <- myplots[[i]] + ylim(c(0,max1))
  }
  
  for(i in c(4,5,6,7,10,11,14,15 )){
    myplots[[i]] <- myplots[[i]] + ylim(c(0,max2))
  }
  
  myplots[[18]] <- ggplot(data, aes(x = data[[1]]+7, y=data[[18]])) +
    geom_line(colour="darkred", size=1.1) + xlab("") +ylab("") + ggtitle(colnames(data)[18]) + theme_bw() + 
    theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
          axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
          legend.position="bottom")
  myplots[[22]] <- ggplot(data, aes(x = data[[1]]+7, y=data[[22]])) +
    geom_line(colour="darkred", size=1.1) + xlab("") +ylab("") + ggtitle(colnames(data)[22]) + theme_bw() + 
    theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
          axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
          legend.position="bottom")
  myplots[[23]] <- ggplot(data, aes(x = data[[1]]+7, y=data[[23]])) +
    geom_line(colour="darkred", size=1.1) + xlab("") +ylab("") + ggtitle(colnames(data)[23]) + theme_bw() + 
    theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
          axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
          legend.position="bottom")
  myplots[[24]] <- ggplot(data, aes(x = data[[1]]+7, y=data[[24]])) +
    geom_line(colour="darkred", size=1.1) + xlab("") +ylab("") + ggtitle(colnames(data)[24]) + theme_bw() + 
    theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
          axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
          legend.position="bottom")
  
  grid <- plot_grid(myplots[[2]],myplots[[3]],myplots[[4]],myplots[[5]],myplots[[6]],
                    myplots[[7]],myplots[[8]],myplots[[9]],myplots[[10]],myplots[[11]],
                    myplots[[12]],myplots[[13]],myplots[[14]],myplots[[15]],myplots[[16]],
                    myplots[[17]],myplots[[18]],myplots[[19]],myplots[[20]],myplots[[21]],
                    myplots[[22]], myplots[[23]],myplots[[24]],myplots[[25]],myplots[[26]],myplots[[27]], nrow = 5, axis="lbtr")
  
  y.grob <- textGrob("Prop. of pop.", 
                     gp=gpar(fontface="bold", fontsize=15), rot=90)
  x.grob <- textGrob("Time (days)", 
                     gp=gpar(fontface="bold", fontsize=15))
  if(plot){grid.arrange(arrangeGrob(grid, left = y.grob, bottom = x.grob))}
  if(save){
    ggsave(grid.arrange(arrangeGrob(grid, left = y.grob, bottom = x.grob)), filename = name, dpi = 300, 
           type = "cairo", width = 25, height = 20, units = "in", path = "figures/")
  }
  options(warn)
}

GT = 5
Re_delta = 1.4
gamma = 1/GT
Re_omicron = 4*Re_delta

# Function will transition between 0 and 1 when V and vr are approximately equal
smooth.transition <- function(V, vr, tune = 0.01){
  sigmoid((V/vr - 1)/tune)
}

#we can boost in V, R1, R2 & R1R2.
#Divide to ratio over V, R1, R2 & R1R2

#if ((V+R1+R2+R1R2-vr1)<0){vr1=0} #todo check for each comp if not < 0
#message(vr1*V/(V+R1+R2+R1R2)," ", V)
# vr1V <- ifelse(vr1*V/(V+R1+R2+R1R2)<=V, vr1, 0)
# vr1R1 <- ifelse(vr1*R1/(V+R1+R2+R1R2)<=R1, vr1, 0)
# vr1R2 <- ifelse(vr1*R2/(V+R1+R2+R1R2)<=R2, vr1, 0)
# vr1R1R2 <- ifelse(vr1*R1R2/(V+R1+R2+R1R2)<=R1R2, vr1, 0)


# Main model --------------------------------------------------------------
VBIR_constant <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    vr1 = vaccination_function(time, 60, vr)
    tot_V = V+R1+R2+R1R2
    j <- smooth.transition(tot_V, vr1)
    
    #smooth the transition when switching of vaccination
    vr1 = j*vr1
    
    dV = -vr1*V/(V+R1+R2+R1R2) - beta_VI1*V*(I1V+I1B+I1R2+I1R2B) - beta_VI2*V*(I2V+I2B+I2R1+I2R1B)
    dB = vr1*V/(V+R1+R2+R1R2) - (1-eff1)*beta_VI1*B*(I1V+I1B+I1R2+I1R2B) - (1-eff2)*beta_VI2*B*(I2V+I2B+I2R1+I2R1B)
    dI1V = beta_VI1*V*(I1V+I1B+I1R2+I1R2B) - gamma*I1V
    dI2V = beta_VI2*V*(I2V+I2B+I2R1+I2R1B) - gamma*I2V
    dI1B = (1-eff1)*beta_VI1*B*(I1V+I1B+I1R2+I1R2B) - gamma*I1B
    dI2B = (1-eff2)*beta_VI2*B*(I2V+I2B+I2R1+I2R1B) - gamma*I2B
    dR1 = gamma*I1V - vr1*R1/(V+R1+R2+R1R2) - (1-eff4)*beta_VI2*R1*(I2V+I2B+I2R1+I2R1B)
    dR2 = gamma*I2V - vr1*R2/(V+R1+R2+R1R2) - (1-eff3)*beta_VI1*R2*(I1V+I1B+I1R2+I1R2B)
    dI1R2 = (1-eff3)*beta_VI1*R2*(I1V+I1B+I1R2+I1R2B) - gamma*I1R2
    dI2R1 = (1-eff4)*beta_VI2*R1*(I2V+I2B+I2R1+I2R1B) - gamma*I2R1
    dR1B = gamma*I1B + vr1*R1/(V+R1+R2+R1R2) - (1-eff6)*beta_VI2*R1B*(I2V+I2B+I2R1+I2R1B)
    dR2B = gamma*I2B + vr1*R2/(V+R1+R2+R1R2) - (1-eff5)*beta_VI1*R2B*(I1V+I1B+I1R2+I1R2B)
    dI1R2B = (1-eff5)*beta_VI1*R2B*(I1V+I1B+I1R2+I1R2B) - gamma*I1R2B
    dI2R1B = (1-eff6)*beta_VI2*R1B*(I2V+I2B+I2R1+I2R1B) - gamma*I2R1B
    dR1R2 = gamma*I2R1 + gamma*I1R2 - vr1*R1R2/(V+R1+R2+R1R2)
    dRall = vr1*R1R2/(V+R1+R2+R1R2) + gamma*I1R2B + gamma*I2R1B
    
    dcumB = vr1*V/(V+R1+R2+R1R2) + vr1*R1/(V+R1+R2+R1R2) + vr1*R2/(V+R1+R2+R1R2) + vr1*R1R2/(V+R1+R2+R1R2)
    
    totalI = I1V+I2V+I1B+I2B+I2R1+I1R2+I1R2B+I2R1B
    totalI1 = I1V+I1B+I1R2+I1R2B
    totalI2 = I2V+I2B+I2R1+I2R1B
    dcumI2 = beta_VI2*V*(I2V+I2B+I2R1+I2R1B)+(1-eff2)*beta_VI2*B*(I2V+I2B+I2R1+I2R1B)+
      (1-eff4)*beta_VI2*R1*(I2V+I2B+I2R1+I2R1B)+(1-eff6)*beta_VI2*R1B*(I2V+I2B+I2R1+I2R1B)
    dcumI1 = beta_VI1*V*(I1V+I1B+I1R2+I1R2B)+(1-eff1)*beta_VI1*B*(I1V+I1B+I1R2+I1R2B)+
      (1-eff3)*beta_VI1*R2*(I1V+I1B+I1R2+I1R2B)+(1-eff5)*beta_VI1*R2B*(I1V+I1B+I1R2+I1R2B)
    
    H = IHR1*I1V + IHR2*I2V + (1-he1)*IHR1*I1B + (1-he2)*IHR2*I2B + (1-he4)*IHR2*I2R1 + (1-he3)*IHR1*I1R2 + (1-he5)*IHR1*I1R2B + (1-he6)*IHR2*I2R1B
    H1 = IHR1*I1V + (1-he1)*IHR1*I1B + (1-he3)*IHR1*I1R2 + (1-he5)*IHR1*I1R2B
    H2 = IHR2*I2V + (1-he2)*IHR2*I2B + (1-he4)*IHR2*I2R1 + (1-he6)*IHR2*I2R1B
    
    dcumH = H
    return(list(c(dV,dB,dI1V,dI2V,dI1B,dI2B,dR1,dR2,dI2R1,dI1R2,dR1B,dR2B,dI1R2B,dI2R1B,dR1R2,dRall,dcumH, dcumB, dcumI1, dcumI2), 
                "H"=H, "H1"=H1, "H2"=H2, "total_I"=totalI,"total_I1"=totalI1,"total_I2"=totalI2))
  })
}



# Baseline simulation -----------------------------------------------------

times <- seq(0,180,by = 0.1)


init = c(V = (0.455 - 0.5*0.01130435 - 7.494118e-06 ),#- 1.028696e-05
         B = (0.245- 0.5*0.008695652 - 5.764706e-06),#- 7.913043e-06
         I1V = 0.5*0.01130435, #0.02*0.455/(0.35+0.455),
         I2V =  7.494118e-06,#1.82E-05*0.455/(0.3+0.35+0.455),
         I1B = 0.5*0.008695652, #0.02*0.35/(0.35+0.455),
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
         cumB = 0.245+0.105,
         cumI1 = 0.01130435 + 0.008695652,
         cumI2 = 7.494118e-06+5.764706e-06+4.941176e-06
)

parms = list(gamma = gamma,
                 vr = 0.02,
                 beta_VI1 = 0.8*Re_delta*gamma,
                 beta_VI2 = 0.4*Re_omicron*gamma,#Re_omicron*gamma,
                 eff1 = 0,
                 eff2 = 0,
                 eff3 = 0,
                 eff4 = 0,
                 #eff5 = 1-(1-parms$eff1)*(1-parms$eff3), #1-1(1-eff1)*(1-eff3)
                 #eff6 = 1-(1-parms$eff2)*(1-parms$eff4), #1-1(1-eff2)*(1-eff4)
                 IHR1 = 0.0024,
                 IHR2 = 0.0024,
                 he1=0.99,
                 he2=0.99,
                 he3=0,
                 he4=0
                 #he5=1-(1-parms$he1)*(1-parms$he3),
                 #he6=1-(1-parms$he2)*(1-parms$he4)
)

parms_dt2 = list(gamma = gamma,
             vr = 0.01,
             beta_VI1 = 0.8*Re_delta*gamma,
             beta_VI2 = 0.5*Re_omicron*gamma,#Re_omicron*gamma,
             eff1 = 0,
             eff2 = 0,
             eff3 = 0,
             eff4 = 0,
             #eff5 = 1-(1-parms$eff1)*(1-parms$eff3), #1-1(1-eff1)*(1-eff3)
             #eff6 = 1-(1-parms$eff2)*(1-parms$eff4), #1-1(1-eff2)*(1-eff4)
             IHR1 = 0.0024,
             IHR2 = 0.5*0.0012,
             he1=0.99,
             he2=0.99,
             he3=0,
             he4=0
             #he5=1-(1-parms$he1)*(1-parms$he3),
             #he6=1-(1-parms$he2)*(1-parms$he4)
)

parms$eff5 = 1-(1-parms$eff1)*(1-parms$eff3)
parms$eff6 = 1-(1-parms$eff2)*(1-parms$eff4)
parms$he5=1-(1-parms$he1)*(1-parms$he3)
parms$he6=1-(1-parms$he2)*(1-parms$he4)


baseline_parms = parms
out <- data.frame(ode(y = init, func = VBIR_constant, times = times, parms = parms))

# Checks for calibrations -------------------------------------------------
out$H[1]*4.4e6
out$H[181]*4.4e6

# Initial R delta:
0.8*Re_delta * ((0.455 - 0.5*0.01130435 - 7.494118e-06)+(0.245- 0.5*0.008695652 - 5.764706e-06)+ (0.3-0.105-4.941176e-06))

#doubling time omicron:
total_I2_t0 = out$total_I2[1]
target = 2*total_I2_t0
t_d = which(abs(out$total_I2 - target) == min(abs(out$total_I2 - target)))
out$time[t_d]
t_d1 = which(out$total_I2 - target>0)[1]
out$time[t_d1]

#initial cases:
(out$cumI1[21]-out$cumI1[11])*4.4e6
1-out$H[181]/out$H[1]

#% vaccinated after 60 days
out$cumB[601]
total_B = (out$B[601]+out$I1B[601]+out$I2B[601]+out$R1B[601]+out$R2B[601]+out$I1R2B[601]+out$I2R1B[601]+out$Rall[601])


trajectory_plots(out, plot=T, save=F, name = "doubling_time_2_days/trajectories_dt_2days_v1.png")

out$time[which(out$total_I2>out$total_I1)[1]]

# Vary betas & start time -------------------------------------------------
VBIR_constant_st <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    vr1 = vaccination_function(time, 60, vr)
    tot_V = V+R1+R2+R1R2
    j <- smooth.transition(tot_V, vr1)
    
    #smooth the transition when switching of vaccination
    vr1 = j*vr1
    if(time>=st){
      beta_VI1 = betas*beta_VI1
      beta_VI2 = betas*beta_VI2
      }
    
    dV = -vr1*V/(V+R1+R2+R1R2) - beta_VI1*V*(I1V+I1B+I1R2+I1R2B) - beta_VI2*V*(I2V+I2B+I2R1+I2R1B)
    dB = vr1*V/(V+R1+R2+R1R2) - (1-eff1)*beta_VI1*B*(I1V+I1B+I1R2+I1R2B) - (1-eff2)*beta_VI2*B*(I2V+I2B+I2R1+I2R1B)
    dI1V = beta_VI1*V*(I1V+I1B+I1R2+I1R2B) - gamma*I1V
    dI2V = beta_VI2*V*(I2V+I2B+I2R1+I2R1B) - gamma*I2V
    dI1B = (1-eff1)*beta_VI1*B*(I1V+I1B+I1R2+I1R2B) - gamma*I1B
    dI2B = (1-eff2)*beta_VI2*B*(I2V+I2B+I2R1+I2R1B) - gamma*I2B
    dR1 = gamma*I1V - vr1*R1/(V+R1+R2+R1R2) - (1-eff4)*beta_VI2*R1*(I2V+I2B+I2R1+I2R1B)
    dR2 = gamma*I2V - vr1*R2/(V+R1+R2+R1R2) - (1-eff3)*beta_VI1*R2*(I1V+I1B+I1R2+I1R2B)
    dI1R2 = (1-eff3)*beta_VI1*R2*(I1V+I1B+I1R2+I1R2B) - gamma*I1R2
    dI2R1 = (1-eff4)*beta_VI2*R1*(I2V+I2B+I2R1+I2R1B) - gamma*I2R1
    dR1B = gamma*I1B + vr1*R1/(V+R1+R2+R1R2) - (1-eff6)*beta_VI2*R1B*(I2V+I2B+I2R1+I2R1B)
    dR2B = gamma*I2B + vr1*R2/(V+R1+R2+R1R2) - (1-eff5)*beta_VI1*R2B*(I1V+I1B+I1R2+I1R2B)
    dI1R2B = (1-eff5)*beta_VI1*R2B*(I1V+I1B+I1R2+I1R2B) - gamma*I1R2B
    dI2R1B = (1-eff6)*beta_VI2*R1B*(I2V+I2B+I2R1+I2R1B) - gamma*I2R1B
    dR1R2 = gamma*I2R1 + gamma*I1R2 - vr1*R1R2/(V+R1+R2+R1R2)
    dRall = vr1*R1R2/(V+R1+R2+R1R2) + gamma*I1R2B + gamma*I2R1B
    
    dcumB = vr1*V/(V+R1+R2+R1R2) + vr1*R1/(V+R1+R2+R1R2) + vr1*R2/(V+R1+R2+R1R2) + vr1*R1R2/(V+R1+R2+R1R2)
    
    totalI = I1V+I2V+I1B+I2B+I2R1+I1R2+I1R2B+I2R1B
    totalI1 = I1V+I1B+I1R2+I1R2B
    totalI2 = I2V+I2B+I2R1+I2R1B
    dcumI2 = beta_VI2*V*(I2V+I2B+I2R1+I2R1B)+(1-eff2)*beta_VI2*B*(I2V+I2B+I2R1+I2R1B)+
      (1-eff4)*beta_VI2*R1*(I2V+I2B+I2R1+I2R1B)+(1-eff6)*beta_VI2*R1B*(I2V+I2B+I2R1+I2R1B)
    dcumI1 = beta_VI1*V*(I1V+I1B+I1R2+I1R2B)+(1-eff1)*beta_VI1*B*(I1V+I1B+I1R2+I1R2B)+
      (1-eff3)*beta_VI1*R2*(I1V+I1B+I1R2+I1R2B)+(1-eff5)*beta_VI1*R2B*(I1V+I1B+I1R2+I1R2B)
    
    H = IHR1*I1V + IHR2*I2V + (1-he1)*IHR1*I1B + (1-he2)*IHR2*I2B + (1-he4)*IHR2*I2R1 + (1-he3)*IHR1*I1R2 + (1-he5)*IHR1*I1R2B + (1-he6)*IHR2*I2R1B
    H1 = IHR1*I1V + (1-he1)*IHR1*I1B + (1-he3)*IHR1*I1R2 + (1-he5)*IHR1*I1R2B
    H2 = IHR2*I2V + (1-he2)*IHR2*I2B + (1-he4)*IHR2*I2R1 + (1-he6)*IHR2*I2R1B
    
    dcumH = H
    return(list(c(dV,dB,dI1V,dI2V,dI1B,dI2B,dR1,dR2,dI2R1,dI1R2,dR1B,dR2B,dI1R2B,dI2R1B,dR1R2,dRall,dcumH, dcumB, dcumI1, dcumI2), 
                "H"=H, "H1"=H1, "H2"=H2, "total_I"=totalI,"total_I1"=totalI1,"total_I2"=totalI2))
  })
}

n.cores <- 7 #parallel::detectCores() - 20
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

times=seq(0,180, by=1)

#eff=0;he=1;IHR=0.0024;betas=1;
#rm(eff,he,IHR,betas)
res_st_beta <- foreach (start_time = iter(seq(1,60,by=1)), .combine='rbind') %:% 
  foreach (betas = iter(seq(0,1.0, by=0.01)), .packages=c("deSolve", "sigmoid"), .combine='rbind') %dopar% {
    parms1 = list(gamma = 0.2,
                  vr = 0.01,
                  beta_VI1 = 0.8*1.4*0.2,
                  beta_VI2 = 0.4*5.6*0.2,
                  eff1 = 0,
                  eff2 = 0,
                  eff3 = 0,
                  eff4 = 0,
                  eff5 = 1.0-(1.0-0)*(1.0-0), 
                  eff6 = 1.0-(1.0-0)*(1.0-0), 
                  IHR1 = 0.0024,
                  IHR2 = 2*0.0012,
                  he1=0.99,
                  he2=0.99,
                  he3=0,
                  he4=0,
                  he5=1.0-(1.0-0.99)*(1.0-0),
                  he6=1.0-(1.0-0.99)*(1.0-0),
                  st = start_time,
                  beta=betas
    )
    
    output = data.frame(ode(y = init, func = VBIR_constant_st, times = times, parms = parms1))
    #get peak hospitalisations, from t=1 onwards.
    peak_H = max(output$H)
    peak_I = max(output$total_I)
    data.frame("start_time"=start_time, "beta"=betas,"peak_H"=peak_H,"peak_I"=peak_I)
  }

parallel::stopCluster(cl = my.cluster)

library(hrbrthemes); library(metR); library(shadowtext);

p1 <- ggplot(res_st_beta, aes(x=start_time, y=beta, fill=peak_H*4.4e6)) + 
  geom_raster() +
  stat_contour(aes(z=peak_H*4.4e6),breaks=c(100,200,300,400),colour="white") +
  geom_text_contour(aes(z = peak_H*4.4e6), breaks=c(100,200,300,400),
                    stroke = 0.2, label.placer = label_placer_flattest(), skip=0) +
  scale_fill_distiller(palette = "RdBu", name = "Peak hospsital\nadmissions per day") +
  xlab("Start time (days after t=0)") +
  ylab("Level of reduction in beta") +
  theme_bw(base_size = 14)

p1

ggsave(plot=p1, filename="figures/intervention_start_betas_pop_peak_H_v2.png",dpi = 300, 
       type = "cairo", width = 12.5, height = 9.38, units = "in")


labels <- data.frame(x=c(5,5,5), y=c(0.64,0.68, 0.76), label=c("1 in 20", "1 in 15", "1 in 10"))
p3 <- ggplot() + 
  geom_raster(data=res_st_beta, aes(x=start_time, y=beta, fill=peak_I*4.4e6)) +
  stat_contour(data=res_st_beta, aes(x=start_time, y=beta, z=peak_I*4.4e6), breaks=c(4.4e6/20, 4.4e6/15, 4.4e6/10), colour="white")+
  scale_fill_distiller(palette = "RdBu", name = "Peak number of\ninfections per day", breaks=c(4.4e6/20, 4.4e6/15, 4.4e6/10), labels=c("1 in 20","1 in 15","1 in 10")) +
  geom_shadowtext(data=labels, aes(x=x, y=y, label=label), bg.r=0.3, bg.colour="white", colour="black") +
  xlim(c(0,31))+
  xlab("Start time (days after t=0)") +
  ylab("Level of reduction in beta") +
  theme_bw(base_size = 14) + theme(legend.key.height= unit(2, 'cm'))
  
p3

ggsave(plot=p3, filename="figures/intervention_start_betas_pop_peak_I_v2.png",dpi = 300, 
       type = "cairo", width = 12.5, height = 9.38, units = "in")

# p2 <- ggplot(res_st_beta, aes(x=start_time, y=beta, fill=peak_H)) + geom_tile() +
#   scale_fill_distiller(palette = "RdYlGn", name = "Peak hospsital\nadmissions per day\n(pop. frac.)") +
#   xlab("Start time (days after t=0)") +
#   ylab("Level of reduction in beta") +
#   theme_bw(base_size = 14)
# 
# ggsave(plot=p2, filename="figures/intervention_start_betas_frac_peak_H.png",dpi = 300, 
#        type = "cairo", width = 12.5, height = 9.38, units = "in")

# p4 <- ggplot(res_st_beta, aes(x=start_time, y=beta, fill=peak_I)) + geom_tile() +
#   scale_fill_distiller(palette = "RdYlGn", name = "Peak infections per day\n (pop. frac.)") +
#   xlab("Start time (days after t=0)") +
#   ylab("Level of reduction in beta") +
#   theme_bw(base_size = 14)
# 
# p4
# 
# ggsave(plot=p4, filename="figures/intervention_start_betas_frac_peak_I.png",dpi = 300, 
#        type = "cairo", width = 12.5, height = 9.38, units = "in")


#################################################################
# Parallel version
#################################################################
n.cores <- 7 #parallel::detectCores() - 20
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

times=seq(0,180, by=1)

#eff=0;he=1;IHR=0.0024;betas=1;
#rm(eff,he,IHR,betas)
res1 <- foreach (eff = iter(c(0,0.1,0.2,0.3,0.4,0.5)), .combine='rbind') %:% 
  foreach (he = iter(c(0,0.2,0.4,0.6,0.8,1.0)), .combine='rbind') %:%
  foreach (IHR = iter(seq(0,0.0012, by = 0.00012)), .combine='rbind') %:% 
  foreach (betas = iter(seq(0,1.0, by=0.01)), .packages=c("deSolve"), .combine='rbind') %dopar% {
    parms1 = list(gamma = 0.2,
                  vr = 0.01,
                  beta_VI1 = betas*0.8*1.4*0.2,
                  beta_VI2 = betas*0.5*5.6*0.2,
                  eff1 = eff,
                  eff2 = eff,
                  eff3 = eff,
                  eff4 = eff,
                  eff5 = 1.0-(1.0-eff)*(1.0-eff), 
                  eff6 = 1.0-(1.0-eff)*(1.0-eff), 
                  IHR1 = 0.0024,
                  IHR2 = IHR,
                  he1=0.99,
                  he2=he,
                  he3=0,
                  he4=he,
                  he5=1.0-(1.0-he)*(1.0-he),
                  he6=1.0-(1.0-he)*(1.0-he)
    )
    
    output = data.frame(ode(y = init, func = VBIR_constant, times = times, parms = parms1))
    #get peak hospitalisations, from t=1 onwards.
    H0 = output$H[1]
    peak_H = max(output$H)
    t_peak_H = output$time[which.max(output$H)]
    data.frame("beta"=betas,"IHR2"=IHR,"he1"=he,"eff1"=eff,"H0"=H0,"peak_H"=peak_H,"t_peak_H"=t_peak_H)
  }

parallel::stopCluster(cl = my.cluster)

res1$he <- factor(res1$he1, levels=c("0", "0.2", "0.4", "0.6", "0.8", "1"))
res1$eff <- factor(res1$eff1, levels=rev(c("0","0.1","0.2","0.3","0.4","0.5")))

res1 <- res1 %>% mutate(col=case_when(peak_H <= H0 ~ -1,
                                      peak_H > H0 & peak_H <= 2*H0 ~ 0,
                                      peak_H > 2*H0 ~ 1,
                                      TRUE ~ -99))

panel_plot <- ggplot(res1, aes(x=beta, y=IHR2, fill=col)) + geom_raster() +
  scale_fill_gradientn(colours=c("green","#ffea00","red"), breaks=c(-1,0,1)) + 
  facet_grid(c("he", "eff"), labeller = "label_both") +
  theme_bw(base_size = 14) + theme(legend.position = "none", strip.background = element_rect(fill="white"), strip.text = element_text(size = 18))


panel_plot

ggsave(plot=panel_plot, filename="figures/doubling_time_2_days/panel_heatmap_peak_H_baseline.png",dpi = 300, 
       type = "cairo", width = 12.5, height = 9.38, units = "in")


#temp code#
n.cores <- 7 #parallel::detectCores() - 20
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

times=seq(0,180, by=1)

#eff=0;he=1;IHR=0.0024;betas=1;
#rm(eff,he,IHR,betas)
res2 <- foreach (eff = iter(c(0,0.1,0.2,0.3,0.4,0.5)), .combine='rbind') %:% 
  foreach (he = iter(c(0,0.2,0.4,0.6,0.8,1.0)), .combine='rbind') %:%
  foreach (IHR = iter(seq(0,2*0.00048, by = 0.000024)), .combine='rbind') %:% 
  foreach (betas = iter(seq(0,1.0, by=0.05)), .packages=c("deSolve", "sigmoid"), .combine='rbind') %dopar% {
    parms1 = list(gamma = 0.2,
                  vr = 0.01,
                  beta_VI1 = betas*0.8*1.4*0.2,
                  beta_VI2 = betas*0.4*5.6*0.2,
                  eff1 = eff,
                  eff2 = eff,
                  eff3 = eff,
                  eff4 = eff,
                  eff5 = 1.0-(1.0-eff)*(1.0-eff), 
                  eff6 = 1.0-(1.0-eff)*(1.0-eff), 
                  IHR1 = 0.0024,
                  IHR2 = IHR,
                  he1=0.99,
                  he2=he,
                  he3=0,
                  he4=he,
                  he5=1.0-(1.0-he)*(1.0-he),
                  he6=1.0-(1.0-he)*(1.0-he)
    )
    
    output = data.frame(ode(y = init, func = VBIR_constant, times = times, parms = parms1))
    #get peak hospitalisations, from t=1 onwards.
    H0 = output$H[1]
    peak_H = max(output$H)
    t_peak_H = output$time[which.max(output$H)]
    print(paste(betas, IHR, he, eff, H0, peak_H, t_peak_H, sep = ","))
    data.frame("beta"=betas,"IHR2"=IHR,"he1"=he,"eff1"=eff,"H0"=H0,"peak_H"=peak_H,"t_peak_H"=t_peak_H)
  }

parallel::stopCluster(cl = my.cluster)

res2$he <- factor(res2$he1, levels=c("0", "0.2", "0.4", "0.6", "0.8", "1"))
res2$eff <- factor(res2$eff1, levels=rev(c("0","0.1","0.2","0.3","0.4","0.5")))

res2 <- res2 %>% mutate(col=case_when(peak_H <= H0 ~ -1,
                                      peak_H > H0 & peak_H <= 2*H0 ~ 0,
                                      peak_H > 2*H0 ~ 1,
                                      TRUE ~ -99))

panel_plot <- ggplot(res2, aes(x=beta, y=IHR2, fill=col)) + geom_raster() +
  scale_fill_gradientn(colours=c("green","#ffea00","red"), breaks=c(-1,0,1)) + 
  facet_grid(c("he", "eff"), labeller = "label_both") +
  theme_bw(base_size = 14) + theme(legend.position = "none", strip.background = element_rect(fill="white"), strip.text = element_text(size = 18))


panel_plot

ggsave(plot=panel_plot, filename="figures/panel_heatmap_peak_H_IHR2=0.2IHR1.png",dpi = 300, 
       type = "cairo", width = 12.5, height = 9.38, units = "in")
###stop temp code


# Univariate Sensitivity Analysis ----------------------------------------------------


sens_parms <- list("vr" = c(0.5*0.01, 0.01, 2*0.01), 
                   "beta_VI1" = c(0.5*0.224, 0.224, 2*0.224), 
                   "beta_VI2" = c(0.5*0.448, 0.448, 2*0.448), 
                   "eff1"= c(0, 0.5, 1), 
                   "eff2" = c(0, 0.5, 1), 
                   "eff3" = c(0, 0.5, 1), 
                   "eff4" = c(0, 0.5, 1),
                   "IHR2" = c(0.0024*0.5, 0.0024, 0.0024*2),
                   "he1" = c(0, 0.5, 1), 
                   "he2" = c(0, 0.5, 1), 
                   "he3" = c(0, 0.5, 1), 
                   "he4" = c(0, 0.5, 1))


times <- seq(0, 180, by = 0.1)
gamma = 0.2
Reff_delta = 1.4
Reff_omicron = 1.4 * 4
init = c(V = (0.455 - 0.5*0.01130435 - 7.494118e-06),#- 1.028696e-05
         B = (0.245- 0.5*0.008695652 - 5.764706e-06),#- 7.913043e-06
         I1V = 0.5*0.01130435, #0.02*0.455/(0.35+0.455),
         I2V =  7.494118e-06,#1.82E-05*0.455/(0.3+0.35+0.455),
         I1B = 0.5*0.008695652, #0.02*0.35/(0.35+0.455),
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
         cumB = 0.245+0.105,
         cumI1 = 0.01130435 + 0.008695652,
         cumI2 = 7.494118e-06+5.764706e-06+4.941176e-06
)


store_dat <- list()

for(i in 1:length(sens_parms)) {
  alter_vec <- sens_parms[[i]]
  parm_list <- list()
  for(z in 1:length(alter_vec)) {
    parms1 <- parms
    parms1[names(sens_parms)[i]] <- alter_vec[[z]]
    parms1$eff5 = 1-(1-parms1$eff1)*(1-parms1$eff3)
    parms1$eff6 = 1-(1-parms1$eff2)*(1-parms1$eff4)
    parms1$he5=1-(1-parms1$he1)*(1-parms1$he3)
    parms1$he6=1-(1-parms1$he2)*(1-parms1$he4)
    out <- data.frame(ode(y = init, func = VBIR_constant, times = times, parms = parms1))
    
    
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
ggsave(p2, filename = "univariate_sens_const_vr.png", dpi = 300, type = "cairo", width = 20, height = 15, units = "in",
       path = "figures/")



# explore beta_VI2 --------------------------------------------------------

sens_parms <- list("beta_VI2" = c(0.3*0.448, 0.35*0.448, 0.4*0.448, 0.45*0.448, 0.5*0.448, 0.55*0.448, 0.6*0.448)) 


times <- seq(0, 180, by = 0.1)
gamma = 0.2
Reff_delta = 1.4
Reff_omicron = 1.4 * 4
init = c(V = (0.455 - 0.5*0.01130435 - 7.494118e-06),#- 1.028696e-05
         B = (0.245- 0.5*0.008695652 - 5.764706e-06),#- 7.913043e-06
         I1V = 0.5*0.01130435, #0.02*0.455/(0.35+0.455),
         I2V =  7.494118e-06,#1.82E-05*0.455/(0.3+0.35+0.455),
         I1B = 0.5*0.008695652, #0.02*0.35/(0.35+0.455),
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
         cumB = 0.245+0.105,
         cumI1 = 0.01130435 + 0.008695652,
         cumI2 = 7.494118e-06+5.764706e-06+4.941176e-06
)


store_dat <- list()

for(i in 1:length(sens_parms)) {
  alter_vec <- sens_parms[[i]]
  parm_list <- list()
  for(z in 1:length(alter_vec)) {
    parms1 <- baseline_parms
    parms1[names(sens_parms)[i]] <- alter_vec[[z]]
    parms1$eff5 = 1-(1-parms1$eff1)*(1-parms1$eff3)
    parms1$eff6 = 1-(1-parms1$eff2)*(1-parms1$eff4)
    parms1$he5=1-(1-parms1$he1)*(1-parms1$he3)
    parms1$he6=1-(1-parms1$he2)*(1-parms1$he4)
    out <- data.frame(ode(y = init, func = VBIR, times = times, parms = parms1))
    
    
    parm_list[[z]] <- data.frame("time" = out["time"],
                                 "H" = out["H"],
                                 "I2" = out$total_I2,
                                 "value" = as.character(alter_vec[[z]]),
                                 "group" = names(sens_parms)[i])
  }
  store_dat[[i]] <- parm_list
}

#Plotting Sensitivity Analyses

plot_list <- list()

for(i in 1:length(store_dat)) {
  plotdat <- do.call(rbind.data.frame, store_dat[[i]])
  m_plotdat <- melt(plotdat, id.vars = c("time", "H", "I2"), measure.vars = "value")
  
  plot_list[[i]] <- ggplot(m_plotdat, aes(x = time, y= H, col = value)) + geom_line(size = 1.5) + theme_bw() + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, max(m_plotdat["H"])*1.2)) + 
    labs(title= paste0("Parameter = ",  names(sens_parms)[i]), 
         y="Incidence of Hospitalisations", x = "Time (days)", colour = "") + 
    theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
          axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
          legend.position="bottom") + scale_color_hue(labels = c("0.3x", "0.35x","0.4x", "0.45x", "0.5x", "0.55x", "0.6x"))
  plot_list[[i+1]] <- ggplot(m_plotdat, aes(x = time, y=I2, col = value)) + geom_line(size = 1.5) + theme_bw() + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, max(m_plotdat["I2"])*1.2)) + 
    labs(title= paste0("Parameter = ",  names(sens_parms)[i]), 
         y="Total I2", x = "Time (days)", colour = "") + 
    theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
          axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
          legend.position="bottom")+ scale_color_hue(labels = c("0.3x", "0.35x","0.4x", "0.45x", "0.5x", "0.55x", "0.6x"))
  
}
plot_list[[1]]
p3<-ggarrange(plot_list[[2]], plot_list[[1]], nrow=1)

p3


ggsave(p3, filename = "exploration_beta_VI2_const_vr.png", dpi = 300, type = "cairo", width = 10, height = 5, units = "in",
       path = "figures/")



# Check for edge cases ----------------------------------------------------

#Low IHR2, he2&4=0

parms = list(gamma = gamma,
             vr = 0.02,
             beta_VI1 = 0.8*Re_delta*gamma,
             beta_VI2 = 0.4*Re_omicron*gamma,#Re_omicron*gamma,
             eff1 = 0,
             eff2 = 0,
             eff3 = 0,
             eff4 = 0,
             #eff5 = 1-(1-parms$eff1)*(1-parms$eff3), #1-1(1-eff1)*(1-eff3)
             #eff6 = 1-(1-parms$eff2)*(1-parms$eff4), #1-1(1-eff2)*(1-eff4)
             IHR1 = 2.0*0.0012,
             IHR2 = 0.2*0.0012,
             he1=0.99,
             he2=0.99,
             he3=0,
             he4=0.99
             #he5=1-(1-parms$he1)*(1-parms$he3),
             #he6=1-(1-parms$he2)*(1-parms$he4)
)
parms$eff5 = 1-(1-parms$eff1)*(1-parms$eff3)
parms$eff6 = 1-(1-parms$eff2)*(1-parms$eff4)
parms$he5=1-(1-parms$he1)*(1-parms$he3)
parms$he6=1-(1-parms$he2)*(1-parms$he4)

out <- data.frame(ode(y = init, func = VBIR, times = times, parms = parms))
trajectory_plots(out, plot=T, save=F)
max(out$total_I)/out$total_I[1]


betas=0.03
eff=0
he=0
IHR=0
times <- seq(0, 180, by = 1)
parms1 = list(gamma = 0.2,
              vr = 0.02,
              beta_VI1 = betas*0.8*1.4*0.2,
              beta_VI2 = betas*0.5*5.6*0.2,
              eff1 = eff,
              eff2 = eff,
              eff3 = eff,
              eff4 = eff,
              eff5 = 1.0-(1.0-eff)*(1.0-eff), 
              eff6 = 1.0-(1.0-eff)*(1.0-eff), 
              IHR1 = 0.0024,
              IHR2 = IHR,
              he1=0.99,
              he2=he,
              he3=0,
              he4=he,
              he5=1.0-(1.0-he)*(1.0-he),
              he6=1.0-(1.0-he)*(1.0-he)
)
out1 <- data.frame(ode(y = init, func = VBIR_constant, times = times, parms = parms1))
trajectory_plots(out1, plot=T, save=F)

library(sigmoid)

# Function will transition between 0 and 1 when h and Q are approximately equal
smooth.transition <- function(h, Q, tune = 0.01){
  sigmoid((h/Q - 1)/tune)
}

Q <- 1
h <- seq(0.001, 5, by = 0.001)
j <- smooth.transition(h, Q)

plot(h/Q, j, type = "l")
