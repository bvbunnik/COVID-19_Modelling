library("deSolve"); library("ggpubr"); library("tidyverse"); library("RColorBrewer")
library("gridExtra"); library("grid");library("cowplot"); library("reshape2")
library("parallel"); library("doParallel")
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



# Main model --------------------------------------------------------------
VBIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dV = -vr*V - beta_VI1*V*(I1V+I1B+I1R2+I1R2B) - beta_VI2*V*(I2V+I2B+I2R1+I2R1B)
    dB = vr*V - (1-eff1)*beta_VI1*B*(I1V+I1B+I1R2+I1R2B) - (1-eff2)*beta_VI2*B*(I2V+I2B+I2R1+I2R1B)
    dI1V = beta_VI1*V*(I1V+I1B+I1R2+I1R2B) - gamma*I1V
    dI2V = beta_VI2*V*(I2V+I2B+I2R1+I2R1B) - gamma*I2V
    dI1B = (1-eff1)*beta_VI1*B*(I1V+I1B+I1R2+I1R2B) - gamma*I1B
    dI2B = (1-eff2)*beta_VI2*B*(I2V+I2B+I2R1+I2R1B) - gamma*I2B
    dR1 = gamma*I1V - vr*R1 - (1-eff4)*beta_VI2*R1*(I2V+I2B+I2R1+I2R1B)
    dR2 = gamma*I2V - vr*R2 - (1-eff3)*beta_VI1*R2*(I1V+I1B+I1R2+I1R2B)
    dI1R2 = (1-eff3)*beta_VI1*R2*(I1V+I1B+I1R2+I1R2B) - gamma*I1R2
    dI2R1 = (1-eff4)*beta_VI2*R1*(I2V+I2B+I2R1+I2R1B) - gamma*I2R1
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
    dcumI2 = beta_VI2*V*(I2V+I2B+I2R1+I2R1B)+(1-eff2)*beta_VI2*B*(I2V+I2B+I2R1+I2R1B)+
      (1-eff4)*beta_VI2*R1*(I2V+I2B+I2R1+I2R1B)+(1-eff6)*beta_VI2*R1B*(I2V+I2B+I2R1+I2R1B)
    dcumI1 = beta_VI1*V*(I1V+I1B+I1R2+I1R2B)+(1-eff1)*beta_VI1*B*(I1V+I1B+I1R2+I1R2B)+
      (1-eff3)*beta_VI1*R2*(I1V+I1B+I1R2+I1R2B)+(1-eff5)*beta_VI1*R2B*(I1V+I1B+I1R2+I1R2B)
    
    H = IHR1*I1V + IHR2*I2V + (1-he1)*IHR1*I1B + (1-he2)*IHR2*I2B + (1-he4)*IHR2*I2R1 + (1-he3)*IHR1*I1R2 + (1-he5)*IHR1*I1R2B + (1-he6)*IHR2*I2R1B
    H1 = IHR1*I1V + (1-he1)*IHR1*I1B + (1-he3)*IHR1*I1R2 + (1-he5)*IHR1*I1R2B
    H2 = IHR2*I2V + (1-he2)*IHR2*I2B + (1-he4)*IHR2*I2R1 + (1-he6)*IHR2*I2R1B
    
    dcumH = H
    return(list(c(dV,dB,dI1V,dI2V,dI1B,dI2B,dR1,dR2,dI2R1,dI1R2,dR1B,dR2B,dI1R2B,dI2R1B,dR1R2,dRall,dcumH, dcumB, dcumI1, dcumI2), 
                "H"=H, "H1"=H1, "H2"=H2,"total_I"=totalI,"total_I1"=totalI1,"total_I2"=totalI2))
  })
}



# Baseline simulation -----------------------------------------------------

times <- seq(0,180,by = 1)


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


parms = list(gamma = gamma,
             vr = 0.03,
             beta_VI1 = 0.8*Re_delta*gamma,
             beta_VI2 = 0.4*Re_omicron*gamma,#Re_omicron*gamma,
             eff1 = 0,
             eff2 = 0,
             eff3 = 0,
             eff4 = 0,
             #eff5 = 1-(1-parms$eff1)*(1-parms$eff3), #1-1(1-eff1)*(1-eff3)
             #eff6 = 1-(1-parms$eff2)*(1-parms$eff4), #1-1(1-eff2)*(1-eff4)
             IHR1 = 2.0*0.0012,
             IHR2 = 2.0*0.0012,
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

out <- data.frame(ode(y = init, func = VBIR, times = times, parms = parms))
out_baseline = out

# Checks for calibrations -------------------------------------------------
out$H[18]*4.4e6

# Initial R delta:
0.8*Re_delta * ((0.455 - 0.5*0.01130435 - 7.494118e-06)+(0.245- 0.5*0.008695652 - 5.764706e-06)+ (0.3-0.105-4.941176e-06))

#doubling time omicron:
total_I2_t0 = out$total_I2[1]
target = 2*total_I2_t0
t_d = which(abs(out$total_I2 - target) == min(abs(out$total_I2 - target)))
out$time[t_d]

#initial cases:
(out$cumI1[2]-out$cumI1[1])*4.4e6
1-out$H[18]/out$H[1]

#% vaccinated after 60 days
out$cumB[60]
total_B = (out$B[60]+out$I1B[60]+out$I2B[60]+out$R1B[60]+out$R2B[60]+out$I1R2B[60]+out$I2R1B[60]+out$Rall[60])

#plot trajectories
trajectory_plots(out, name="baseline_trajactories_v1.png")


#create 5x5 panel plot, parameters to vary: 
# - all betas scaled 0-1,
# - IHR2 scaled 0 to IHR1
# - he2&he4, scaled by 0, 0.2, 0.4, 0.6, 0.8, 1.0, he6 to follow
# - eff1-4 scaled by 0, 0.2, 0.4, 0.6, 0.8, 1.0, eff5, eff6 to follow
res3 = data.frame(matrix(ncol=5,nrow=1, dimnames=list(NULL, c("beta", "IHR2", "he", "eff", "peak_H"))))
times=seq(0,180, by=1)

#eff=0;he=0.01;IHR=0.01;betas=1;
rm(eff,he,IHR,betas)
for (eff in c(0.5,0.6,0.7,0.8,0.9,1.0)){
  for (he in c(0.01, 0.11, 0.21, 0.31, 0.41, 0.51)){
    for (IHR in seq(0,0.01, by = 0.001)){
      for (betas in seq(0.5,1.0, by=0.05)){
        parms1 = list(gamma = 0.2,
                     vr = 0.025,
                     beta_VI1 = betas*0.4*1.4*0.2,
                     beta_VI2 = betas*0.4*5.6*0.2,
                     eff1 = eff,
                     eff2 = eff,
                     eff3 = eff,
                     eff4 = eff,
                     eff5 = 1.0-(1.0-eff)*(1.0-eff), 
                     eff6 = 1.0-(1.0-eff)*(1.0-eff), 
                     IHR1 = 0.01,
                     IHR2 = IHR,
                     he1=0.99,
                     he2=1-he,
                     he3=0,
                     he4=1-he,
                     he5=0.99,
                     he6=1.0-(1.0-he)*(1.0-he)
        )
        output = data.frame(ode(y = init, func = VBIR, times = times, parms = parms1))
        #get peak hospitalisations, from t=1 onwards.
        peak_H = 0
        if(output$time[which.max(output$H)]>1){
          peak_H = max(output$H)
        } else {
          peak_H = output$H[1]
        }
        res3 = rbind(res3, c(betas,IHR,he,eff,peak_H))
        rm(output,parms, peak_H)
      }
    }
  }
}

max(output$H)
output$total_I1[1]
max(output$H)>2*baseline_H0

res3 <- res3[-1,]

#names(res) <- c("beta", "IHR2", "he", "eff", "peak_H")
#write_csv(res, "data/sensitivity_panel.csv")

baseline_H0 = out_baseline$H[1]
new_baseline_H0 = output$H[1]

res3$he_f <- factor(res3$he, levels=rev(c("0.01", "0.11", "0.21", "0.31", "0.41", "0.51")))
res3$eff_f <- factor(res3$eff, levels=c("1", "0.8", "0.6", "0.4", "0.2","0"))

res3 <- res3 %>% mutate(col=case_when(peak_H<new_baseline_H0 ~ -1,
                              peak_H>=new_baseline_H0 & peak_H<=2*new_baseline_H0 ~ 0,
                              peak_H>2*new_baseline_H0 ~ 1,
                              TRUE ~ -99))
min(res2$col)
ggplot(res3, aes(x=beta, y=IHR2, fill=col)) + geom_raster() +
  scale_fill_gradientn(colours=c("green","yellow","red"), breaks=c(-1,0,1)) + facet_grid(c("he_f", "eff_f"), labeller = "label_both")


  
ggsave(filename="figures/panel_heatmap_peak_H_v6.png",dpi = 300, 
       type = "cairo", width = 12.5, height = 9.38, units = "in")


ggplot(output, aes(x=time, y=H)) + geom_line()



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
  foreach (IHR = iter(seq(0,2*0.0024, by = 0.00024)), .combine='rbind') %:% 
  foreach (betas = iter(seq(0,1.0, by=0.05)), .packages=c("deSolve"), .combine='rbind') %dopar% {
    parms1 = list(gamma = 0.2,
                  vr = 0.03,
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

    output = data.frame(ode(y = init, func = VBIR, times = times, parms = parms1))
    #get peak hospitalisations, from t=1 onwards.
    H0 = output$H[1]
    peak_H = max(output$H)
    t_peak_H = output$time[which.max(output$H)]
    data.frame("beta"=betas,"IHR2"=IHR,"he1"=he,"eff1"=eff,"H0"=H0,"peak_H"=peak_H,"t_peak_H"=t_peak_H)
  }

parallel::stopCluster(cl = my.cluster)

baseline_H0 = output$H[1]

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

ggsave(plot=panel_plot, filename="figures/panel_heatmap_peak_H_doubleIHR2.png",dpi = 300, 
       type = "cairo", width = 12.5, height = 9.38, units = "in")


# Investigative runs for beta=0 -------------------------------------------
times <- seq(0,180,by = 1)
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
parms_adj = list(gamma = gamma,
             vr = 0.03,
             beta_VI1 = 0.8*Re_delta*gamma,
             beta_VI2 = 0.4*Re_omicron*gamma,#Re_omicron*gamma,
             eff1 = 0,
             eff2 = 0,
             eff3 = 0,
             eff4 = 0,
             #eff5 = 1-(1-parms$eff1)*(1-parms$eff3), #1-1(1-eff1)*(1-eff3)
             #eff6 = 1-(1-parms$eff2)*(1-parms$eff4), #1-1(1-eff2)*(1-eff4)
             IHR1 = 2.0*0.0012,
             IHR2 = 2.0*0.0012,
             he1=0.99,
             he2=0.99,
             he3=0,
             he4=0
             #he5=1-(1-parms$he1)*(1-parms$he3),
             #he6=1-(1-parms$he2)*(1-parms$he4)
)

parms_adj$beta_VI1 = 0*parms_adj$beta_VI1
parms_adj$beta_VI2 = 0*parms_adj$beta_VI2
parms_adj$eff5 = 1-(1-parms_adj$eff1)*(1-parms_adj$eff3)
parms_adj$eff6 = 1-(1-parms_adj$eff2)*(1-parms_adj$eff4)
parms_adj$he2 = parms_adj$he4 = 0.8
parms_adj$he5=1-(1-parms_adj$he1)*(1-parms_adj$he3)
parms_adj$he6=1-(1-parms_adj$he2)*(1-parms_adj$he4)

output_adj = data.frame(ode(y = init, func = VBIR, times = times, parms = parms_adj))


out_m = melt(output_adj, id.vars = "time")


#make V,B,R1,R2 same scale:
max1=out_m %>% filter(variable %in% c("V","B","R1","R2")) %>% summarise(max1=max(value)) %>% as.numeric

#make V,B,R1,R2 same scale:
max2=out_m %>% filter(variable %in% c("I1V","I2V","I1B","I2B","I2R1", "I1R2","I1R2B","I2R1B")) %>% summarise(max2=max(value)) %>% as.numeric

#make plots for all comps:
myplots <- vector('list', ncol(output_adj))
for (i in seq_along(output_adj)) {
  message(i)
  myplots[[i]] <- local({
    i <- i
    p1 <- ggplot(output_adj[,c(1,i)], aes(x = output_adj[[1]], y=output_adj[[i]])) +
      geom_line(colour="darkred", size=1.1) + xlab("") +ylab("") + ggtitle(colnames(output_adj)[i]) + theme_bw() +
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

myplots[[18]] <- ggplot(output_adj, aes(x = output_adj[[1]]+7, y=output_adj[[18]])) +
  geom_line(colour="darkred", size=1.1) + xlab("") +ylab("") + ggtitle(colnames(output_adj)[18]) + theme_bw() + 
  theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
        axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position="bottom")
myplots[[22]] <- ggplot(output_adj, aes(x = output_adj[[1]]+7, y=output_adj[[21]])) +
  geom_line(colour="darkred", size=1.1) + xlab("") +ylab("") + ggtitle(colnames(output_adj)[21]) + theme_bw() + 
  theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
        axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position="bottom")
myplots[[23]] <- ggplot(output_adj, aes(x = output_adj[[1]]+7, y=output_adj[[22]])) +
  geom_line(colour="darkred", size=1.1) + xlab("") +ylab("") + ggtitle(colnames(output_adj)[22]) + theme_bw() + 
  theme(plot.title = element_text(size=18), legend.text=element_text(size=15), axis.text=element_text(size=15), 
        axis.title.y=element_text(size=15), axis.title.x= element_text(size=15), plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position="bottom")
myplots[[24]] <- ggplot(output_adj, aes(x = output_adj[[1]]+7, y=output_adj[[23]])) +
  geom_line(colour="darkred", size=1.1) + xlab("") +ylab("") + ggtitle(colnames(output_adj)[23]) + theme_bw() + 
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

pp <- grid.arrange(arrangeGrob(grid, left = y.grob, bottom = x.grob))

pp

ggsave(plot=pp, filename="check_beta_he=0.8_betas=0.png", dpi = 300, 
       type = "cairo", width = 25, height = 20, units = "in", path = "figures/")

output_adj$time[which.max(output_adj$H)]>14
which.max(output_adj$H)
output_adj$time[which.max(output_adj$H)]
output_adj$H[1]

# Univariate Sensitivity Analysis ----------------------------------------------------

# vr = 0.5x, baseline, 2x
# beta_VI1 = baseline, 2x, 4x 
# beta_VI2 = baseline, 2x, 4x 
# eff1= 0, 0.5, 1, 
# eff2 = 0, 0.5, 1, 
# eff3 = 0, 0.5, 1, 
# eff4 = 0, 0.5, 1, 
# eff5&eff6 to follow eff1,eff3 & eff2,eff4
# IHR2 = 0.1x, baseline, 2x),
# he1 = 0.1x, 0.5x, 1x, 
# he2 = 0.1x, 0.5x, 1x, 
# he3 = 0.1x, 0.5x, 1x, 
# he4 = 0.1x, 0.5x, 1x
# he5&he6 to follow he1,he3 & he2,h24

sens_parms <- list("vr" = c(0.015, 0.03, 0.06), 
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
ggsave(p2, filename = "univariate_sens_v1.png", dpi = 300, type = "cairo", width = 20, height = 15, units = "in",
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




ggsave(p3, filename = "exploration_beta_VI2.png", dpi = 300, type = "cairo", width = 10, height = 5, units = "in",
       path = "figures/")

