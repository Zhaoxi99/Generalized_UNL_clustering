library(dplyr)
library(lubridate)
library(ggplot2)
require(mcclust.ext)
library(UNL.est)
library(corrplot)
library(LSBP)

source("D:/PhD_study2_desktop/lddp_functions.R")

load("D:/PhD_study2_desktop/new_app/dde_app/dde.RData")

dde$GAD=dde$GAD/7
ggplot(data=dde, aes(x=DDE,y=GAD)) + geom_point(alpha=.5, cex=.5) + geom_smooth( method="loess", span = 1, col=1) + xlab("DDE (mg/L)") + ylab("Gestational age at delivery (in weeks)") + theme_bw() 

data=data.frame(scale(dde))

##conditional fit
fit_lm <- lm(GAD ~DDE,data = data)
summary(fit_lm)
X <- model.matrix(fit_lm);
k <- ncol(X)
mcmc <- list(nsave = 10000, nburn = 10000, nskip = 1)

#data dependent prior
prior <- list(m0 = fit_lm$coefficients,
              S0 = solve(t(X)%*%X)*(sigma(fit_lm)^2),
              nu = k + 2,
              Psi = 30*(solve(t(X)%*%X)*(sigma(fit_lm)^2)),
              a = 2,
              b = (sigma(fit_lm)^2)/2,
              aalpha = 2, balpha=2,
              L = 20)

# prior <- list(m0 = c(0,0),
#               S0 = solve(t(X)%*%X)*(sigma(fit_lm)^2),
#               nu = k + 2,
#               Psi = diag(2),
#               a = 1,
#               b = 1,
#               aalpha = 2, balpha=2,
#               L = 20)

set.seed(176)
fit_lddp <- lddp_moves(y = data$GAD,
                       X = X,
                       prior = prior, 
                       mcmc = mcmc,
                       standardise = FALSE)


psm_lddp <- comp.psm(fit_lddp$z)
output_vi_lddp <- minVI(psm_lddp, fit_lddp$z)
save.image("//csce.datastore.ed.ac.uk/csce/maths/groups/mdt/clustering_paper/application_data/dde/dde_fit_partitions.RData")

plotpsm(psm_lddp)

##is wasabi worthwhile?
library(WASABI)

ptm=proc.time()
set.seed(16789)
out_WASABI_2p <- WASABI(fit_lddp$z, psm = psm_lddp, L = 2,
                        method.init = "++", method = "salso")
proc.time()-ptm

ptm=proc.time()
set.seed(15359)
out_WASABI_3p <- WASABI(fit_lddp$z, psm = psm_lddp, L = 3,
                        method.init = "++", method = "salso")
proc.time()-ptm

ggsummary(out_WASABI_2p)
ggsummary(out_WASABI_3p)

table(out_WASABI_2p$particles[1,])
table(out_WASABI_3p$particles[1,])
plotpsm(psm_lddp)

table(out_WASABI_3p$particles[1,])
table(out_WASABI_3p$particles[2,])
table(output_vi_lddp$cl)

ggplot(dde)+geom_point(aes(x=DDE,y=GAD,col=as.factor(out_WASABI_2p$particles[1,])),alpha=.5, cex=.5)+
  xlab("DDE (mg/L)") + ylab("Gestational age at delivery (in weeks)") +labs(col="cluster")+
  geom_hline(yintercept = 37, color = "red", linewidth = 0.45)+
  geom_hline(yintercept = 42, color = "blue", linewidth = 0.45)+
  theme_bw()

ggplot(dde)+geom_point(aes(x=DDE,y=GAD,col=as.factor(out_WASABI_3p$particles[1,])),alpha=.5, cex=.5)+
  xlab("DDE (mg/L)") + ylab("Gestational age at delivery (in weeks)") +labs(col="cluster")+
  geom_hline(yintercept = 37, color = "red", linewidth = 0.45)+
  geom_hline(yintercept = 42, color = "blue", linewidth = 0.45)+
  theme_bw()


###
ind1 <- which(output_vi_lddp$cl == 1)
ind2 <- which(output_vi_lddp$cl == 2)
ind3 <- which(output_vi_lddp$cl == 3)

sum(dde$GAD[ind2]<37)/length(ind2)

sum(dde$GAD[ind3]>=42)/length(ind3)

sum(dde$GAD[ind1]>=37 & dde$GAD[ind1]<42)/length(ind1)
##regression in each cluster
fit_lm1 <- lm(GAD ~DDE,data = dde[ind1,])
summary(fit_lm1)

fit_lm2 <- lm(GAD ~DDE,data = dde[ind2,])
summary(fit_lm2)

fit_lm3 <- lm(GAD ~DDE,data = dde[ind3,])
summary(fit_lm3)

png(filename = "D:/PhD_study2_desktop/plots/application_plots/dde_plots/dde_cluster_partition.png",
    width = 609*5, height = 469*5,res = 72*5)
theme_set(theme_bw())
ggplot(dde)+geom_point(aes(x=DDE,y=GAD,col=as.factor(output_vi_lddp$cl)),alpha=.5, cex=.5)+
  xlab("DDE (mg/L)") + ylab("Gestational age at delivery (in weeks)") +labs(col="cluster")+
  geom_hline(yintercept = 37, color = "red", linewidth = 0.45,linetype = "dashed")+
  geom_hline(yintercept = 42, color = "blue", linewidth = 0.45,linetype = "dashed")+
  geom_abline(intercept = coef(fit_lm1)[1], slope = coef(fit_lm2)[2], 
              color = "red", linewidth = 0.45)+
  geom_abline(intercept = coef(fit_lm2)[1], slope = coef(fit_lm2)[2], 
              color = "#4DAF4A", linewidth = 0.45)+
  geom_abline(intercept = coef(fit_lm3)[1], slope = coef(fit_lm2)[2], 
              color = "blue", linewidth = 0.45)+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        plot.title.position = "plot",            # 可选：使标题相对于整个绘图区域定位
        plot.subtitle = element_text(hjust = 0.5),
        text =  element_text(size = 15),
        legend.text = element_text(size = 15)
  )
dev.off()


# ggplot(dde)+geom_point(aes(x=DDE,y=GAD,col=as.factor(out_WASABI$particles[1,])))+theme_bw()
# ggplot(dde)+geom_point(aes(x=DDE,y=GAD,col=as.factor(out_WASABI2$particles[3,])))+theme_bw()
# ggplot(dde)+geom_point(aes(x=DDE,y=GAD,col=as.factor(out_WASABI2$particles[2,])))+theme_bw()
# 
# plotpsm(psm_lddp)


##

x_c1 <- as.matrix(X[ind1,c(2)]);samples1=list(y=x_c1,y_cate=NULL)
x_c2 <- as.matrix(X[ind2,c(2)]);samples2=list(y=x_c2,y_cate=NULL)
x_c3 <- as.matrix(X[ind3,c(2)]);samples3=list(y=x_c3,y_cate=NULL)


##prior for group 1,2,3,4
L=10
set.seed(158)
prior1=prior_dpm(samples1, L=L, K=L%/%2, nstart = 5,categories = NULL)
set.seed(158)
prior2=prior_dpm(samples2, L=L, K=L%/%2, nstart = 5,categories = NULL)
set.seed(158)
prior3=prior_dpm(samples3, L=L, K=L%/%2, nstart = 5,categories = NULL)



ptm=proc.time()
set.seed(159)
res_1 <- dpm_MN_Mcate_shrink(y_all = samples1, prior = prior1$prior_shrink, mcmc = mcmc, standardise = FALSE) 
set.seed(159)
res_2 <- dpm_MN_Mcate_shrink(y_all = samples2, prior = prior2$prior_shrink, mcmc = mcmc, standardise = FALSE) 
set.seed(159)
res_3 <- dpm_MN_Mcate_shrink(y_all = samples3, prior = prior3$prior_shrink, mcmc = mcmc, standardise = FALSE) 
proc.time()-ptm
##calculate UNL
library(future)
library(future.apply)
nsave_list <- list()
for (i in 1:mcmc$nsave) {
  nsave_list[[as.character(i)]] <- i
}

plan(multisession, workers = detectCores()-2)
res_list=list(res_1,res_2,res_3)
#UNL for joint density
ptm=proc.time()
UNL_imp_lddp=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_list,
                           n_imp=5000,continuous_slct_index=c(1),
                           cate_slct_index=NULL,future.seed = TRUE)
proc.time()-ptm

tranform_list<-function(simulation_object){
  #simulation_object=UNL_imp
  nsave=length(names(simulation_object))
  unl=rep(0,nsave)
  eff_size=rep(0,nsave)
  for (i in 1:nsave){
    unl[i]=simulation_object[[i]]$unl
    eff_size[i]=simulation_object[[i]]$eff_size
  }
  return(list(unl=unl,eff_size=eff_size))
}

# png(filename = "D:/PhD_study2_desktop/plots/application_plots/dde_plots/dde_UNL_hist.png",
#     width = 609*5, height = 469*5,res = 72*5)
# hist(tranform_list(UNL_imp_lddp)$unl,xlim=c(1,3),col = rgb(0, 0, 1, 0.2),
#      freq = FALSE,xlab = "UNL",main = "3 clusters (conditional)")
# dev.off()

cols <- c(
  rgb(0, 0, 1, 0.2),
  rgb(0, 1, 0, 0.2),
  rgb(1, 0, 0, 0.2)
)

png(filename = "D:/PhD_study2_desktop/plots/application_plots/dde_plots/dde_UNL_hist.png",
    width = 609*5, height = 469*5,res = 72*5)

# 合并为 tidy 格式
df <- data.frame(
  unl = tranform_list(UNL_imp_lddp)$unl
)

# 绘图
ggplot(df, aes(x = unl)) +
  geom_histogram(aes(y = ..density..),        # freq = FALSE -> density
                 position = "identity",       # 叠加而非堆叠
                 alpha = 0.45,                # 透明度，便于比较
                 bins = 30,                   # 可调整或改为 binwidth = ...
                 color = "black",fill=rgb(0, 0, 1, 0.2)) +          # 柱子边框 +
  coord_cartesian(xlim = c(1, 3), ylim = c(0, 12)) +   # 保留你原来的 x/y 范围
  labs(x = "UNL", y = "Density") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    legend.position = c(0.98, 0.98),                  # 图内右上角（模拟 base::legend("topright")）
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank()
  )
dev.off()

summary(tranform_list(UNL_imp_lddp)$unl)

##conditional densities
DDE.points  <- round(quantile(data$DDE,c(0.1,0.6,0.9,0.99)),2)
X_pred_lddp<- cbind(1,DDE.points) 
sequenceGAD <- seq(from=min(data$GAD),to=max(data$GAD),length=100)

# Posterior density - Gibbs sampling
pred_density_lddp <- array(0,c(mcmc$nsave,length(sequenceGAD),4))

ptm=proc.time()
pred_density_lddp=density_lddp(fit_lddp = fit_lddp,X=X_pred_lddp,grid = sequenceGAD)/sd(dde$GAD)
proc.time()-ptm

estimate_lddp <- apply(pred_density_lddp,c(2,3),mean)
lower_lddp    <- apply(pred_density_lddp,c(2,3),function(x) quantile(x,0.025))
upper_lddp    <- apply(pred_density_lddp,c(2,3),function(x) quantile(x,0.975))

#fit lsbp
prior_lbp  <- prior_LSBP(2,2, 
                         b_mixing = rep(0,2), B_mixing=diag(1,2), 
                         b_kernel = rep(0,2), B_kernel=diag(1,2), 
                         a_tau = 1, b_tau= 1)
model_formula <- as.list(Formula::as.Formula(GAD ~ DDE | DDE))

ptm=proc.time()
set.seed(110) # The seed is setted so that the Gibbs sampler is reproducible.
fit_lbp   <- LSBP_Gibbs(Formula=model_formula, data=data, H=20, prior=prior_lbp, 
                        control=control_Gibbs(R=mcmc$nsave,burn_in=10000,method_init="random"), verbose=TRUE)
proc.time()-ptm

X1           <- cbind(1,DDE.points)             # Design matrix for the kernel
X2           <- cbind(1,DDE.points)

pred_density_lsbp <- array(0,c(mcmc$nsave,length(sequenceGAD),4))
for(r in 1:mcmc$nsave){      # Cycle over the iterations of the MCMC chain
  for(i in 1:100){  # Cycle over the GAD grid
    pred_density_lsbp[r,i,] <- c(LSBP_density(sequenceGAD[i],X1,X2,
                                              fit_lbp$param$beta_mixing[r,,],
                                              fit_lbp$param$beta_kernel[r,,],
                                              fit_lbp$param$tau[r,]))/sd(dde$GAD)
  }
}

estimate_lsbp <- apply(pred_density_lsbp,c(2,3),mean)
lower_lsbp    <- apply(pred_density_lsbp,c(2,3),function(x) quantile(x,0.025))
upper_lsbp    <- apply(pred_density_lsbp,c(2,3),function(x) quantile(x,0.975))

data_density_plot1 <- data.frame(
  estimate_lddp  = estimate_lddp[,1],
  lower_lddp       = lower_lddp[,1],
  upper_lddp       = upper_lddp[,1],
  estimate_lsbp  = estimate_lsbp[,1],
  lower_lsbp       = lower_lsbp[,1],
  upper_lsbp       = upper_lsbp[,1],
  sequenceGAD = sequenceGAD*sd(dde$GAD) + mean(dde$GAD))

data_density_plot2 <- data.frame(
  estimate_lddp  = estimate_lddp[,2],
  lower_lddp       = lower_lddp[,2],
  upper_lddp       = upper_lddp[,2],
  estimate_lsbp  = estimate_lsbp[,2],
  lower_lsbp       = lower_lsbp[,2],
  upper_lsbp       = upper_lsbp[,2],
  sequenceGAD = sequenceGAD*sd(dde$GAD) + mean(dde$GAD))

data_density_plot3 <- data.frame(
  estimate_lddp  = estimate_lddp[,3],
  lower_lddp       = lower_lddp[,3],
  upper_lddp       = upper_lddp[,3],
  estimate_lsbp  = estimate_lsbp[,3],
  lower_lsbp       = lower_lsbp[,3],
  upper_lsbp       = upper_lsbp[,3],
  sequenceGAD = sequenceGAD*sd(dde$GAD) + mean(dde$GAD))

data_density_plot4 <- data.frame(
  estimate_lddp  = estimate_lddp[,4],
  lower_lddp       = lower_lddp[,4],
  upper_lddp       = upper_lddp[,4],
  estimate_lsbp  = estimate_lsbp[,4],
  lower_lsbp       = lower_lsbp[,4],
  upper_lsbp       = upper_lsbp[,4],
  sequenceGAD = sequenceGAD*sd(dde$GAD) + mean(dde$GAD))

data_hist_plot1 <- data.frame(
  GAD = dde$GAD[which(dde$DDE < 20.505)])
data_hist_plot2 <- data.frame(
  GAD = dde$GAD[which(dde$DDE >= 20.505 & dde$DDE < 41.08)])
data_hist_plot3 <- data.frame(
  GAD = dde$GAD[which(dde$DDE >= 41.08 & dde$DDE < 79.6)])
data_hist_plot4 <- data.frame(
  GAD = dde$GAD[which(dde$DDE > 79.6)])

png(filename = "D:/PhD_study2_desktop/plots/application_plots/dde_plots/dde_LSBP_lddp_conditional1.png",
    width = 609*5, height = 469*5,res = 72*5)
binwidth=2*(quantile(data_hist_plot1$GAD,0.75)-quantile(data_hist_plot1$GAD,0.25))/(length(data_hist_plot1$GAD)^(1/3))
ggplot(data=data_density_plot1) + 
  geom_line(aes(x=sequenceGAD,y=estimate_lddp),col='red')+ 
  geom_line(aes(x=sequenceGAD,y=estimate_lsbp),col='blue')+ 
  geom_vline(xintercept = 37, color = "red", linewidth = 0.45,linetype = "dashed")+
  geom_vline(xintercept = 42, color = "blue", linewidth = 0.45,linetype = "dashed")+
  ylab("Density") + geom_ribbon(alpha=0.4,fill='red',aes(x=sequenceGAD,ymin=lower_lddp,ymax=upper_lddp))+
  geom_ribbon(alpha=0.4,fill='blue',aes(x=sequenceGAD,ymin=lower_lsbp,ymax=upper_lsbp))+
  xlab("Gestational age at delivery (in weeks)") +
  geom_histogram(data=data_hist_plot1,aes(x=GAD,y=after_stat(density)),alpha=0.2,binwidth = binwidth)+
  ggtitle("DDE = 12.57 (10th percentile)")+
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
    plot.title.position = "plot",            # 可选：使标题相对于整个绘图区域定位
    plot.subtitle = element_text(hjust = 0.5)
  )
dev.off()

png(filename = "D:/PhD_study2_desktop/plots/application_plots/dde_plots/dde_LSBP_lddp_conditional2.png",
    width = 609*5, height = 469*5,res = 72*5)
binwidth=2*(quantile(data_hist_plot2$GAD,0.75)-quantile(data_hist_plot2$GAD,0.25))/(length(data_hist_plot2$GAD)^(1/3))
ggplot(data=data_density_plot2) + 
  geom_line(aes(x=sequenceGAD,y=estimate_lddp),col='red')+ 
  geom_line(aes(x=sequenceGAD,y=estimate_lsbp),col='blue')+ 
  geom_vline(xintercept = 37, color = "red", linewidth = 0.45,linetype = "dashed")+
  geom_vline(xintercept = 42, color = "blue", linewidth = 0.45,linetype = "dashed")+
  ylab("Density") + geom_ribbon(alpha=0.4,fill='red',aes(x=sequenceGAD,ymin=lower_lddp,ymax=upper_lddp))+
  geom_ribbon(alpha=0.4,fill='blue',aes(x=sequenceGAD,ymin=lower_lsbp,ymax=upper_lsbp))+
  xlab("Gestational age at delivery") +
  geom_histogram(data=data_hist_plot2,aes(x=GAD,y=after_stat(density)),alpha=0.2,binwidth = binwidth)+
  ggtitle("DDE = 28.44 (60th percentile) (in weeks)")+
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        plot.title.position = "plot",            # 可选：使标题相对于整个绘图区域定位
        plot.subtitle = element_text(hjust = 0.5)
  )
dev.off()

png(filename = "D:/PhD_study2_desktop/plots/application_plots/dde_plots/dde_LSBP_lddp_conditional3.png",
    width = 609*5, height = 469*5,res = 72*5)
binwidth=2*(quantile(data_hist_plot3$GAD,0.75)-quantile(data_hist_plot3$GAD,0.25))/(length(data_hist_plot3$GAD)^(1/3))
ggplot(data=data_density_plot3) + 
  geom_line(aes(x=sequenceGAD,y=estimate_lddp),col='red')+ 
  geom_line(aes(x=sequenceGAD,y=estimate_lsbp),col='blue')+ 
  geom_vline(xintercept = 37, color = "red", linewidth = 0.45,linetype = "dashed")+
  geom_vline(xintercept = 42, color = "blue", linewidth = 0.45,linetype = "dashed")+
  ylab("Density") + geom_ribbon(alpha=0.4,fill='red',aes(x=sequenceGAD,ymin=lower_lddp,ymax=upper_lddp))+
  geom_ribbon(alpha=0.4,fill='blue',aes(x=sequenceGAD,ymin=lower_lsbp,ymax=upper_lsbp))+
  xlab("Gestational age at delivery (in weeks)") +
  geom_histogram(data=data_hist_plot3,aes(x=GAD,y=after_stat(density)),alpha=0.2,binwidth = binwidth)+
  ggtitle("DDE = 53.72 (90th percentile)")+
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        plot.title.position = "plot",            # 可选：使标题相对于整个绘图区域定位
        plot.subtitle = element_text(hjust = 0.5)
  )
dev.off()

png(filename = "D:/PhD_study2_desktop/plots/application_plots/dde_plots/dde_LSBP_lddp_conditional4.png",
    width = 609*5, height = 469*5,res = 72*5)
binwidth=2*(quantile(data_hist_plot4$GAD,0.75)-quantile(data_hist_plot4$GAD,0.25))/(length(data_hist_plot4$GAD)^(1/3))
ggplot(data=data_density_plot4) + 
  geom_line(aes(x=sequenceGAD,y=estimate_lddp),col='red')+ 
  geom_line(aes(x=sequenceGAD,y=estimate_lsbp),col='blue')+ 
  geom_vline(xintercept = 37, color = "red", linewidth = 0.45,linetype = "dashed")+
  geom_vline(xintercept = 42, color = "blue", linewidth = 0.45,linetype = "dashed")+
  ylab("Density") + geom_ribbon(alpha=0.4,fill='red',aes(x=sequenceGAD,ymin=lower_lddp,ymax=upper_lddp))+
  geom_ribbon(alpha=0.4,fill='blue',aes(x=sequenceGAD,ymin=lower_lsbp,ymax=upper_lsbp))+
  xlab("Gestational age at delivery (in weeks)") +
  geom_histogram(data=data_hist_plot4,aes(x=GAD,y=after_stat(density)),alpha=0.2,binwidth = binwidth)+
  ggtitle("DDE = 105.47 (99th percentile)")+
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        plot.title.position = "plot",            # 可选：使标题相对于整个绘图区域定位
        plot.subtitle = element_text(hjust = 0.5)
  )
dev.off()


##posterior predictive
ptm=proc.time()
set.seed(167)
y_post_lddp=predic_lddp(fit_lddp =fit_lddp,X=X)*sd(dde$GAD) + mean(dde$GAD)
proc.time()-ptm

dim(y_post_lddp)
df <- reshape2::melt(y_post_lddp)
df$x_value <- dde$DDE[df$Var1]

ggplot(data=dde, aes(x=DDE,y=GAD)) + coord_cartesian(xlim = c(0, 200), ylim = c(25, 45))+
  geom_point(alpha=.5, cex=.5) + xlab("DDE (mg/L)") + ylab("Gestational age at delivery (in weeks)") + theme_bw() 

ggplot(df, aes(x = x_value, y = value)) +
  coord_cartesian(xlim = c(0, 200), ylim = c(25, 45))+
  stat_bin2d(aes(fill = after_stat(density)))+
  scale_fill_gradient(low = "white", high = "blue") +
  xlab("DDE (mg/L)") + ylab("Gestational age at delivery (in weeks)")+
  theme_bw()

png(filename = "D:/PhD_study2_desktop/plots/application_plots/dde_plots/post_check_dde_kurtosis.png",
    width = 609*5, height = 469*5,res = 72*5)
data_kurtosis=data.frame(kurtosis=apply(y_post_lddp,FUN = moments::kurtosis,MARGIN = 1))
original_kurtosis=moments::kurtosis(dde$GAD)
theme_set(theme_bw())
ggplot(data_kurtosis, aes(x=kurtosis)) + 
  geom_histogram( fill="slategray4", color="slategray4")+
  geom_vline(xintercept = original_kurtosis, color="lightblue", linewidth=1)+
  xlab("Kurtosis of GAD")+ylab("Counts")+  theme(     plot.title = element_text(size = 20, face = "bold", hjust = 0.5),     axis.title = element_text(size = 20),     axis.text = element_text(size = 15)   )
dev.off()

png(filename = "D:/PhD_study2_desktop/plots/application_plots/dde_plots/post_check_dde_skewness.png",
    width = 609*5, height = 469*5,res = 72*5)
data_skewness=data.frame(skewness=apply(y_post_lddp,FUN = moments::skewness,MARGIN = 1))
original_skewness=moments::skewness(dde$GAD)
theme_set(theme_bw())
ggplot(data_skewness, aes(x=skewness)) + 
  geom_histogram(fill="slategray4", color="slategray4")+
  geom_vline(xintercept = original_skewness, color="lightblue", linewidth=1)+
  xlab("Skewness of GAD")+ylab("Counts")+  theme(     plot.title = element_text(size = 20, face = "bold", hjust = 0.5),     axis.title = element_text(size = 20),     axis.text = element_text(size = 15)   )
dev.off()

png(filename = "D:/PhD_study2_desktop/plots/application_plots/dde_plots/post_check_dde_min.png",
    width = 609*5, height = 469*5,res = 72*5)
data_min=data.frame(min=apply(y_post_lddp,FUN = min,MARGIN = 1))
original_min=min(dde$GAD)
theme_set(theme_bw())
ggplot(data_min, aes(x=min)) + 
  geom_histogram( fill="slategray4", color="slategray4")+
  geom_vline(xintercept = original_min, color="lightblue", linewidth=1)+
  xlab("Minmum value of GAD")+ylab("Counts")+  theme(     plot.title = element_text(size = 20, face = "bold", hjust = 0.5),     axis.title = element_text(size = 20),     axis.text = element_text(size = 15)   )
dev.off()

png(filename = "D:/PhD_study2_desktop/plots/application_plots/dde_plots/post_check_dde_max.png",
    width = 609*5, height = 469*5,res = 72*5)
data_max=data.frame(max=apply(y_post_lddp,FUN = max,MARGIN = 1))
original_max=max(dde$GAD)
theme_set(theme_bw())
ggplot(data_max, aes(x=max)) + 
  geom_histogram( fill="slategray4", color="slategray4")+
  geom_vline(xintercept = original_max, color="lightblue", linewidth=1)+
  xlim(45,50)+
  xlab("Maximum value of GAD")+ylab("Counts")+  theme(     plot.title = element_text(size = 20, face = "bold", hjust = 0.5),     axis.title = element_text(size = 20),     axis.text = element_text(size = 15)   )
dev.off()

##

thresholds=DDE.points
X1=cbind(1,data$DDE[data$DDE<thresholds[1]])
X2=cbind(1,data$DDE[data$DDE>thresholds[1]&data$DDE<thresholds[2]])
X3=cbind(1,data$DDE[data$DDE>thresholds[2]&data$DDE<thresholds[3]])
X4=cbind(1,data$DDE[data$DDE>thresholds[3]])


ptm=proc.time()
set.seed(167)
y_post_lddp_bin1=predic_lddp(fit_lddp =fit_lddp,X=X1)*sd(dde$GAD) + mean(dde$GAD)
y_post_lddp_bin2=predic_lddp(fit_lddp =fit_lddp,X=X2)*sd(dde$GAD) + mean(dde$GAD)
y_post_lddp_bin3=predic_lddp(fit_lddp =fit_lddp,X=X3)*sd(dde$GAD) + mean(dde$GAD)
y_post_lddp_bin4=predic_lddp(fit_lddp =fit_lddp,X=X4)*sd(dde$GAD) + mean(dde$GAD)
proc.time()-ptm

png(filename = "D:/PhD_study2_desktop/plots/application_plots/dde_plots/post_check_dde_dense1.png",
    width = 609*5, height = 469*5,res = 72*5)
plot(density(y_post_lddp_bin1[,1],n=2048),xlim=c(25,50),ylim=c(0,0.3),xlab = "Gestational age at delivery (in weeks)",
     main = "DDE level < 12.57 (10th percentile)",col="slategray4",cex.main=1.7,cex.lab=1.5)
for(r in 2:mcmc$nsave){
  lines(density(y_post_lddp_bin1[,r]),col="slategray4")
}
lines(density(dde$GAD[data$DDE<thresholds[1]],n=2048)$x,
      density(dde$GAD[data$DDE<thresholds[1]],n=2048)$y,col="lightblue",lwd=4)
dev.off()


png(filename = "D:/PhD_study2_desktop/plots/application_plots/dde_plots/post_check_dde_dense2.png",
    width = 609*5, height = 469*5,res = 72*5)
plot(density(y_post_lddp_bin2[,1],n=2048),xlim=c(25,50),ylim=c(0,0.3),xlab = "Gestational age at delivery (in weeks)",
     main = "12.57< DDE level < 28.44 (between the 10th and 60th percentiles)",col="slategray4",cex.main=1.7,cex.lab=1.5)
for(r in 2:mcmc$nsave){
  lines(density(y_post_lddp_bin2[,r]),col="slategray4")
}
lines(density(dde$GAD[data$DDE>thresholds[1]&data$DDE<thresholds[2]],n=2048)$x,
      density(dde$GAD[data$DDE>thresholds[1]&data$DDE<thresholds[2]],n=2048)$y,col="lightblue",lwd=4)
dev.off()


png(filename = "D:/PhD_study2_desktop/plots/application_plots/dde_plots/post_check_dde_dense3.png",
    width = 609*5, height = 469*5,res = 72*5)
plot(density(y_post_lddp_bin3[,1],n=2048),xlim=c(25,50),ylim=c(0,0.3),xlab = "Gestational age at delivery (in weeks)",
     main = "28.44< DDE level < 53.72 (between the 60th and 90th percentiles)",col="slategray4",cex.main=1.7,cex.lab=1.5)
for(r in 2:mcmc$nsave){
  lines(density(y_post_lddp_bin3[,r]),col="slategray4")
}
lines(density(dde$GAD[data$DDE>thresholds[2]&data$DDE<thresholds[3]],n=2048)$x,
      density(dde$GAD[data$DDE>thresholds[2]&data$DDE<thresholds[3]],n=2048)$y,col="lightblue",lwd=4)
dev.off()


png(filename = "D:/PhD_study2_desktop/plots/application_plots/dde_plots/post_check_dde_dense4.png",
    width = 609*5, height = 469*5,res = 72*5)
plot(density(y_post_lddp_bin4[,1],n=2048),xlim=c(25,50),ylim=c(0,0.3),xlab = "Gestational age at delivery (in weeks)",
     main = "DDE level >53.72 (90th percentile)",col="slategray4",cex.main=1.7,cex.lab=1.5)
for(r in 2:mcmc$nsave){
  lines(density(y_post_lddp_bin4[,r]),col="slategray4")
}
lines(density(dde$GAD[data$DDE>thresholds[3]],n=2048)$x,
      density(dde$GAD[data$DDE>thresholds[3]],n=2048)$y,col="lightblue",lwd=4)
dev.off()

save.image("D:/PhD_study2_desktop/new_app/dde_app/dde_test5.RData")
