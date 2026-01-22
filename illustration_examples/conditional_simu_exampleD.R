library(ggplot2)
require(mcclust.ext)
library(UNL.est)

source("D:/PhD_study2_desktop/lddp_functions.R")
set.seed(27149)
n=1000;p=20
sigma2_x=diag(4,p)
sigma2_x_non0_index1_raw=seq(2,2*floor(p/2),by=2)
sigma2_x_non0_index2_raw=seq(1,2*floor((p-1)/2)+1,by=2)

sigma2_x_non0_index1=combn(sigma2_x_non0_index1_raw, 2)
sigma2_x_non0_index2=combn(sigma2_x_non0_index2_raw, 2)

for(i in 1:dim(sigma2_x_non0_index1)[2]){
  sigma2_x[sigma2_x_non0_index1[1,i],sigma2_x_non0_index1[2,i]]=3
  sigma2_x[sigma2_x_non0_index1[2,i],sigma2_x_non0_index1[1,i]]=3
}

for(i in 1:dim(sigma2_x_non0_index2)[2]){
  sigma2_x[sigma2_x_non0_index2[1,i],sigma2_x_non0_index2[2,i]]=3
  sigma2_x[sigma2_x_non0_index2[2,i],sigma2_x_non0_index2[1,i]]=3
}

x=mvrnorm(n=n,mu=rep(4,p),Sigma=sigma2_x)
beta_10=0;beta_11=1;
beta_20=4.5;beta_21=0.1;
sigma1=1/4;sigma2=sqrt(1/8);mu1=4;mu2=6;tau1=2;tau2=2

p_weight=tau1*exp(-0.5*tau1^2*(x[,1]-mu1)^2)/(tau1*exp(-0.5*tau1^2*(x[,1]-mu1)^2)
                                              +tau2*exp(-0.5*tau2^2*(x[,1]-mu2)^2))
y1=rnorm(n,mean=beta_10+beta_11*x[,1],sd=sigma1)
y2=rnorm(n,mean=beta_20+beta_21*x[,1],sd=sigma2)

u_y=runif(n)
y1_index=u_y<p_weight
y2_index=u_y>=p_weight
y=rep(0,n)
y[y1_index]=y1[y1_index];y[y2_index]=y2[y2_index]
data=data.frame(y=y,x)

plot(x[,1],y)

theme_set(theme_bw())
ggplot(data)+geom_point(aes(x=X1,y=y,col=as.factor(u_y<p_weight)))+
  labs(x="x",y="y",col="cluster")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        plot.title.position = "plot",            # 可选：使标题相对于整个绘图区域定位
        plot.subtitle = element_text(hjust = 0.5),
        text =  element_text(size = 15),
        legend.text = element_text(size = 15)
  )

##conditional fit
fit_lm <- lm(y ~.,data = data)
X <- model.matrix(fit_lm);
k <- ncol(X)

#data dependent prior
prior <- list(m0=fit_lm$coefficients,
              S0 = solve(t(X)%*%X)*(sigma(fit_lm)^2), 
              nu = k + 2, 
              Psi = 30*(solve(t(X)%*%X)*(sigma(fit_lm)^2)),
              a = 2,
              b = (sigma(fit_lm)^2)/2,
              aalpha = 2, balpha=2,
              L = 20)

mcmc <- list(nsave = 10000, nburn = 10000, nskip = 1)


set.seed(179)
fit_lddp <- lddp_moves(y = y,
                       X = X,
                       prior = prior, 
                       mcmc = mcmc,
                       standardise = FALSE)


psm_lddp <- comp.psm(fit_lddp$z)
output_vi_lddp <- minVI(psm_lddp, fit_lddp$z)


png("D:/PhD_study2_desktop/plots/example_plots/example_D_lddp_fit.png",
    width = 609*5, height = 469*5,res = 72*5)
theme_set(theme_bw())
ggplot(data)+geom_point(aes(x=X1,y=y,col=as.factor(output_vi_lddp$cl)))+
  labs(x="x1",y="y",col="cluster")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        plot.title.position = "plot",            # 可选：使标题相对于整个绘图区域定位
        plot.subtitle = element_text(hjust = 0.5),
        text =  element_text(size = 15),
        legend.text = element_text(size = 15)
  )
dev.off()

ind1 <- which(output_vi_lddp$cl == 1)
ind2 <- which(output_vi_lddp$cl == 2)

x_c1 <- x[ind1,];samples1=list(y=x_c1,y_cate=NULL)
x_c2 <- x[ind2,];samples2=list(y=x_c2,y_cate=NULL)

##prior for group 1,2,3
L=10
set.seed(123)
prior1=prior_dpm(samples1, L=L, K=50, nstart = 5)
set.seed(123)
prior2=prior_dpm(samples2, L=L, K=50, nstart = 5)
#det(prior2$prior_full$L0)

ptm=proc.time()
set.seed(159)
res_1_shrink <- dpm_MN_Mcate_shrink(y_all = samples1, prior = prior1$prior_shrink, mcmc = mcmc, standardise = FALSE) 
set.seed(159)
res_2_shrink <- dpm_MN_Mcate_shrink(y_all = samples2, prior = prior2$prior_shrink, mcmc = mcmc, standardise = FALSE) 
set.seed(159)
res_1 <- dpm_MN_Mcate(y_all = samples1, prior = prior1$prior_full, mcmc = mcmc, standardise = FALSE) 
set.seed(159)
res_2 <- dpm_MN_Mcate(y_all = samples2, prior = prior2$prior_full, mcmc = mcmc, standardise = FALSE)  
proc.time()-ptm

##calculate UNL
library(future)
library(future.apply)
nsave_list <- list()
for (i in 1:mcmc$nsave) {
  nsave_list[[as.character(i)]] <- i
}

plan(multisession, workers = detectCores()-2)
res_list=list(res_1,res_2);res_list_shrink=list(res_1_shrink,res_2_shrink)
#UNL for joint density
ptm=proc.time()
set.seed(125)
UNL_imp_full_shrink=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_list_shrink,
                                  n_imp=5000,continuous_slct_index=1:20,
                                  cate_slct_index=NULL,future.seed = TRUE)
set.seed(125)
UNL_imp_pindex1_shrink=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_list_shrink,
                                     n_imp=5000,continuous_slct_index=sigma2_x_non0_index1_raw,
                                     cate_slct_index=NULL,future.seed = TRUE)
set.seed(125)
UNL_imp_pindex2_shrink=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_list_shrink,
                                     n_imp=5000,continuous_slct_index=sigma2_x_non0_index2_raw,
                                     cate_slct_index=NULL,future.seed = TRUE)
set.seed(125)
UNL_imp_p1_shrink=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_list_shrink,
                                n_imp=5000,continuous_slct_index=c(1),
                                cate_slct_index=NULL,future.seed = TRUE)

ptm=proc.time()
set.seed(125)
UNL_imp_full=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_list,
                           n_imp=5000,continuous_slct_index=1:20,
                           cate_slct_index=NULL,future.seed = TRUE)
set.seed(125)
UNL_imp_pindex1=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_list,
                              n_imp=5000,continuous_slct_index=sigma2_x_non0_index1_raw,
                              cate_slct_index=NULL,future.seed = TRUE)
set.seed(125)
UNL_imp_pindex2=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_list,
                              n_imp=5000,continuous_slct_index=sigma2_x_non0_index2_raw,
                              cate_slct_index=NULL,future.seed = TRUE)
set.seed(125)
UNL_imp_p1=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_list,
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

png("D:/PhD_study2_desktop/plots/example_plots/example_D_UNLhist_shrink.png",
    width = 609*5, height = 469*5,res = 72*5)
cols <- c(
  rgb(0, 0, 1, 0.2),
  rgb(0, 1, 0, 0.2),
  rgb(1, 0, 0, 0.2),
  rgb(0.5, 0, 0.5, 0.2)
)
m1    <- tranform_list(UNL_imp_full_shrink)$unl
m2    <- tranform_list(UNL_imp_pindex1_shrink)$unl
m3    <- tranform_list(UNL_imp_pindex2_shrink)$unl
m4    <- tranform_list(UNL_imp_p1_shrink)$unl

df <- data.frame(
  unl = c(m1, m2,m3,m4),
  group = factor(rep(
    c("all covariates",
      "even-indexed covariates","odd-indexed covariates",expression(x[1])),
    times = c(length(m1), length(m2),length(m3),length(m4))
  ),levels =c("all covariates",
              "even-indexed covariates","odd-indexed covariates",expression(x[1])))
)

mycols <- setNames(cols, c("all covariates",
                           "even-indexed covariates","odd-indexed covariates",expression(x[1])))

# 绘图（density：y = ..density..）
ggplot(df, aes(x = unl, fill = group)) +
  geom_histogram(aes(y = after_stat(density)),
                 position = "identity",    # 叠加显示
                 alpha = 0.45,             # 透明度，便于比较
                 bins = 30,          # 可根据需要调整 binwidth 或用 bins = 30
                 color = "black") +        # 柱子边框
  coord_cartesian(xlim = c(1, 2)) +   # 保留原来你给的 xlim/ylim
  scale_fill_manual(values = mycols,
                    name = NULL) +
  labs(x = "UNL", y = "Density") +
  theme_minimal() +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    legend.position = c(0.26, 0.98),                  # 图内右上角（模拟 base::legend("topright")）
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank()
  )
dev.off()


png("D:/PhD_study2_desktop/plots/example_plots/example_D_UNLhist.png",
    width = 609*5, height = 469*5,res = 72*5)
cols <- c(
  rgb(0, 0, 1, 0.2),
  rgb(0, 1, 0, 0.2),
  rgb(1, 0, 0, 0.2),
  rgb(0.5, 0, 0.5, 0.2)
)
m1    <- tranform_list(UNL_imp_full)$unl
m2    <- tranform_list(UNL_imp_pindex1)$unl
m3    <- tranform_list(UNL_imp_pindex2)$unl
m4    <- tranform_list(UNL_imp_p1)$unl

df <- data.frame(
  unl = c(m1, m2,m3,m4),
  group = factor(rep(
    c("all covariates",
      "even-indexed covariates","odd-indexed covariates","x_1"),
    times = c(length(m1), length(m2),length(m3),length(m4))
  ),levels =c("all covariates",
              "even-indexed covariates","odd-indexed covariates","x_1"))
)

mycols <- setNames(cols, c("all covariates",
                           "even-indexed covariates","odd-indexed covariates","x_1"))

# 绘图（density：y = ..density..）
ggplot(df, aes(x = unl, fill = group)) +
  geom_histogram(aes(y = after_stat(density)),
                 position = "identity",    # 叠加显示
                 alpha = 0.45,             # 透明度，便于比较
                 bins = 30,          # 可根据需要调整 binwidth 或用 bins = 30
                 color = "black") +        # 柱子边框
  coord_cartesian(xlim = c(1, 2)) +   # 保留原来你给的 xlim/ylim
  scale_fill_manual(values = mycols,
                    name = NULL) +
  labs(x = "UNL", y = "Density") +
  theme_minimal() +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    legend.position = c(0.38, 0.98),                  # 图内右上角（模拟 base::legend("topright")）
    legend.justification = c(1, 1),
    legend.text = element_text(size = 16),
    legend.background = element_blank(),
    legend.key = element_blank()
  )
dev.off()

thresholds=quantile(x[,1],probs = c(0.25,0.5,0.75))
X1=X[data$X1<thresholds[1],]
X2=X[data$X1>=thresholds[1]&data$X1<thresholds[2],]
X3=X[data$X1>=thresholds[2]&data$X1<thresholds[3],]
X4=X[data$X1>=thresholds[3],]

ptm=proc.time()
set.seed(167)
y_post_lddp_bin1=predic_lddp(fit_lddp =fit_lddp,X=X1)
set.seed(167)
y_post_lddp_bin2=predic_lddp(fit_lddp =fit_lddp,X=X2)
set.seed(167)
y_post_lddp_bin3=predic_lddp(fit_lddp =fit_lddp,X=X3)
set.seed(167)
y_post_lddp_bin4=predic_lddp(fit_lddp =fit_lddp,X=X4)
proc.time()-ptm

png("D:/PhD_study2_desktop/plots/example_plots/example_D_postcheck_dens1.png",
    width = 609*5, height = 469*5,res = 72*5)
par(cex.axis = 1.5,
    cex.lab  = 1.8,
    cex.main = 1.8,
    mar = c(5.1, 5.1, 4.1, 1.1))
plot(density(y_post_lddp_bin1[,1],n=2048),xlim=c(-2.5,6),ylim=c(0,0.8),xlab = "y",
     main = expression(paste(x[1],"<2.66 (25th quantile)"),col="slategray4"))
for(r in 2:mcmc$nsave){
  lines(density(y_post_lddp_bin1[,r]),col="slategray4")
}
lines(density(data$y[data$X1<thresholds[1]],n=2048)$x,
      density(data$y[data$X1<thresholds[1]],n=2048)$y,col="lightblue",lwd=4)
dev.off()

png("D:/PhD_study2_desktop/plots/example_plots/example_D_postcheck_dens2.png",
    width = 609*5, height = 469*5,res = 72*5)
par(cex.axis = 1.5,
    cex.lab  = 1.8,
    cex.main = 1.8,
    mar = c(5.1, 5.1, 4.1, 1.1))
plot(density(y_post_lddp_bin2[,1],n=2048),xlim=c(1,6),ylim=c(0,1.2),xlab = "y",
     main = expression(paste("2.66<",x[1],"<3.99 (between the 25th and 50th quantile)"),col="slategray4"))
for(r in 2:mcmc$nsave){
  lines(density(y_post_lddp_bin2[,r]),col="slategray4")
}
lines(density(data$y[data$X1>=thresholds[1]&data$X1<thresholds[2]],n=2048)$x,
      density(data$y[data$X1>=thresholds[1]&data$X1<thresholds[2]],n=2048)$y,col="lightblue",lwd=4)
dev.off()

png("D:/PhD_study2_desktop/plots/example_plots/example_D_postcheck_dens3.png",
    width = 609*5, height = 469*5,res = 72*5)
par(cex.axis = 1.5,
    cex.lab  = 1.8,
    cex.main = 1.8,
    mar = c(5.1, 5.1, 4.1, 1.1))
plot(density(y_post_lddp_bin3[,1],n=2048),xlim=c(2,8),ylim=c(0,1.2),xlab = "y",
     main = expression(paste("3.99<",x[1],"<5.43 (between the 50th and 75th quantile)"),col="slategray4"))
for(r in 2:mcmc$nsave){
  lines(density(y_post_lddp_bin3[,r]),col="slategray4")
}
lines(density(data$y[data$X1>=thresholds[2]&data$X1<thresholds[3]],n=2048)$x,
      density(data$y[data$X1>=thresholds[2]&data$X1<thresholds[3]],n=2048)$y,col="lightblue",lwd=4)
dev.off()

png("D:/PhD_study2_desktop/plots/example_plots/example_D_postcheck_dens4.png",
    width = 609*5, height = 469*5,res = 72*5)
par(cex.axis = 1.5,
    cex.lab  = 1.8,
    cex.main = 1.8,
    mar = c(5.1, 5.1, 4.1, 1.1))
plot(density(y_post_lddp_bin4[,1],n=2048),xlim=c(2,10),ylim=c(0,1.1),xlab = "y",
     main = expression(paste(x[1],">5.43 (75th quantile)"),col="slategray4"))
for(r in 2:mcmc$nsave){
  lines(density(y_post_lddp_bin4[,r]),col="slategray4")
}
lines(density(data$y[data$X1>=thresholds[3]],n=2048)$x,
      density(data$y[data$X1>=thresholds[3]],n=2048)$y,col="lightblue",lwd=4)
dev.off()

##conditioned on x2 kernel density plots
thresholds=quantile(x[,2],probs = c(0.25,0.5,0.75))
X1=X[data$X2<thresholds[1],]
X2=X[data$X2>=thresholds[1]&data$X2<thresholds[2],]
X3=X[data$X2>=thresholds[2]&data$X2<thresholds[3],]
X4=X[data$X2>=thresholds[3],]

ptm=proc.time()
set.seed(167)
y_post_lddp_bin1_x2=predic_lddp(fit_lddp =fit_lddp,X=X1)
set.seed(167)
y_post_lddp_bin2_x2=predic_lddp(fit_lddp =fit_lddp,X=X2)
set.seed(167)
y_post_lddp_bin3_x2=predic_lddp(fit_lddp =fit_lddp,X=X3)
set.seed(167)
y_post_lddp_bin4_x2=predic_lddp(fit_lddp =fit_lddp,X=X4)
proc.time()-ptm

png("D:/PhD_study2_desktop/plots/example_plots/example_D_postcheck_dens1_x2.png",
    width = 609*5, height = 469*5,res = 72*5)
par(cex.axis = 1.5,
    cex.lab  = 1.8,
    cex.main = 1.8,
    mar = c(5.1, 5.1, 4.1, 1.1))
plot(density(y_post_lddp_bin1_x2[,1],n=2048),xlim=c(-2.5,10),ylim=c(0,0.4),xlab = "y",
     main = expression(paste(x[2],"<2.57 (25th quantile)"),col="slategray4"))
for(r in 2:mcmc$nsave){
  lines(density(y_post_lddp_bin1_x2[,r]),col="slategray4")
}
lines(density(data$y[data$X2<thresholds[1]],n=2048)$x,
      density(data$y[data$X2<thresholds[1]],n=2048)$y,col="lightblue",lwd=4)
dev.off()

png("D:/PhD_study2_desktop/plots/example_plots/example_D_postcheck_dens2_x2.png",
    width = 609*5, height = 469*5,res = 72*5)
par(cex.axis = 1.5,
    cex.lab  = 1.8,
    cex.main = 1.8,
    mar = c(5.1, 5.1, 4.1, 1.1))
plot(density(y_post_lddp_bin2_x2[,1],n=2048),xlim=c(-2.5,10),ylim=c(0,0.4),xlab = "y",
     main = expression(paste("2.57<",x[2],"<3.98 (between the 25th and 50th quantile)"),col="slategray4"))
for(r in 2:mcmc$nsave){
  lines(density(y_post_lddp_bin2_x2[,r]),col="slategray4")
}
lines(density(data$y[data$X2>=thresholds[1]&data$X2<thresholds[2]],n=2048)$x,
      density(data$y[data$X2>=thresholds[1]&data$X2<thresholds[2]],n=2048)$y,col="lightblue",lwd=4)
dev.off()

png("D:/PhD_study2_desktop/plots/example_plots/example_D_postcheck_dens3_x2.png",
    width = 609*5, height = 469*5,res = 72*5)
par(cex.axis = 1.5,
    cex.lab  = 1.8,
    cex.main = 1.8,
    mar = c(5.1, 5.1, 4.1, 1.1))
plot(density(y_post_lddp_bin3_x2[,1],n=2048),xlim=c(-2.5,10),ylim=c(0,0.4),xlab = "y",
     main = expression(paste("3.98<",x[2],"<5.16 (between the 50th and 75th quantile)"),col="slategray4"))
for(r in 2:mcmc$nsave){
  lines(density(y_post_lddp_bin3_x2[,r]),col="slategray4")
}
lines(density(data$y[data$X2>=thresholds[2]&data$X2<thresholds[3]],n=2048)$x,
      density(data$y[data$X2>=thresholds[2]&data$X2<thresholds[3]],n=2048)$y,col="lightblue",lwd=4)
dev.off()

png("D:/PhD_study2_desktop/plots/example_plots/example_D_postcheck_dens4_x2.png",
    width = 609*5, height = 469*5,res = 72*5)
par(cex.axis = 1.5,
    cex.lab  = 1.8,
    cex.main = 1.8,
    mar = c(5.1, 5.1, 4.1, 1.1))
plot(density(y_post_lddp_bin4_x2[,1],n=2048),xlim=c(-2.5,10),ylim=c(0,0.4),xlab = "y",
     main = expression(paste(x[2],">5.16 (75th quantile)"),col="slategray4"))
for(r in 2:mcmc$nsave){
  lines(density(y_post_lddp_bin4_x2[,r]),col="slategray4")
}
lines(density(data$y[data$X2>=thresholds[3]],n=2048)$x,
      density(data$y[data$X2>=thresholds[3]],n=2048)$y,col="lightblue",lwd=4)
dev.off()


##conditioned on x3 kernel density plots
thresholds=quantile(x[,3],probs = c(0.25,0.5,0.75))
X1=X[data$X3<thresholds[1],]
X2=X[data$X3>=thresholds[1]&data$X3<thresholds[2],]
X3=X[data$X3>=thresholds[2]&data$X3<thresholds[3],]
X4=X[data$X3>=thresholds[3],]

ptm=proc.time()
set.seed(167)
y_post_lddp_bin1_x3=predic_lddp(fit_lddp =fit_lddp,X=X1)
set.seed(167)
y_post_lddp_bin2_x3=predic_lddp(fit_lddp =fit_lddp,X=X2)
set.seed(167)
y_post_lddp_bin3_x3=predic_lddp(fit_lddp =fit_lddp,X=X3)
set.seed(167)
y_post_lddp_bin4_x3=predic_lddp(fit_lddp =fit_lddp,X=X4)
proc.time()-ptm

png("D:/PhD_study2_desktop/plots/example_plots/example_D_postcheck_dens1_x3.png",
    width = 609*5, height = 469*5,res = 72*5)
par(cex.axis = 1.5,
    cex.lab  = 1.8,
    cex.main = 1.8,
    mar = c(5.1, 5.1, 4.1, 1.1))
plot(density(y_post_lddp_bin1_x3[,1],n=2048),xlim=c(-2.5,8),ylim=c(0,0.6),xlab = "y",
     main = expression(paste(x[3],"<2.67 (25th quantile)"),col="slategray4"))
for(r in 2:mcmc$nsave){
  lines(density(y_post_lddp_bin1_x3[,r]),col="slategray4")
}
lines(density(data$y[data$X3<thresholds[1]],n=2048)$x,
      density(data$y[data$X3<thresholds[1]],n=2048)$y,col="lightblue",lwd=4)
dev.off()

png("D:/PhD_study2_desktop/plots/example_plots/example_D_postcheck_dens2_x3.png",
    width = 609*5, height = 469*5,res = 72*5)
par(cex.axis = 1.5,
    cex.lab  = 1.8,
    cex.main = 1.8,
    mar = c(5.1, 5.1, 4.1, 1.1))
plot(density(y_post_lddp_bin2_x3[,1],n=2048),xlim=c(-2,8),ylim=c(0,0.6),xlab = "y",
     main = expression(paste("2.67<",x[3],"<4.03 (between the 25th and 50th quantile)"),col="slategray4"))
for(r in 2:mcmc$nsave){
  lines(density(y_post_lddp_bin2_x3[,r]),col="slategray4")
}
lines(density(data$y[data$X3>=thresholds[1]&data$X3<thresholds[2]],n=2048)$x,
      density(data$y[data$X3>=thresholds[1]&data$X3<thresholds[2]],n=2048)$y,col="lightblue",lwd=4)
dev.off()

png("D:/PhD_study2_desktop/plots/example_plots/example_D_postcheck_dens3_x3.png",
    width = 609*5, height = 469*5,res = 72*5)
par(cex.axis = 1.5,
    cex.lab  = 1.8,
    cex.main = 1.8,
    mar = c(5.1, 5.1, 4.1, 1.1))
plot(density(y_post_lddp_bin3_x3[,1],n=2048),xlim=c(0,8),ylim=c(0,0.6),xlab = "y",
     main = expression(paste("4.03<",x[3],"<5.47 (between the 50th and 75th quantile)"),col="slategray4"))
for(r in 2:mcmc$nsave){
  lines(density(y_post_lddp_bin3_x3[,r]),col="slategray4")
}
lines(density(data$y[data$X3>=thresholds[2]&data$X3<thresholds[3]],n=2048)$x,
      density(data$y[data$X3>=thresholds[2]&data$X3<thresholds[3]],n=2048)$y,col="lightblue",lwd=4)
dev.off()

png("D:/PhD_study2_desktop/plots/example_plots/example_D_postcheck_dens4_x3.png",
    width = 609*5, height = 469*5,res = 72*5)
par(cex.axis = 1.5,
    cex.lab  = 1.8,
    cex.main = 1.8,
    mar = c(5.1, 5.1, 4.1, 1.1))
plot(density(y_post_lddp_bin4_x3[,1],n=2048),xlim=c(2,10),ylim=c(0,0.8),xlab = "y",
     main = expression(paste(x[3],">5.47 (75th quantile)"),col="slategray4"))
for(r in 2:mcmc$nsave){
  lines(density(y_post_lddp_bin4_x3[,r]),col="slategray4")
}
lines(density(data$y[data$X3>=thresholds[3]],n=2048)$x,
      density(data$y[data$X3>=thresholds[3]],n=2048)$y,col="lightblue",lwd=4)
dev.off()

ptm=proc.time()
set.seed(167)
y_post_lddp=predic_lddp(fit_lddp =fit_lddp,X=X)
proc.time()-ptm


png(filename = "D:/PhD_study2_desktop/plots/example_plots/example_D_postcheck_kurtosis.png",
    width = 609*5, height = 469*5,res = 72*5)
data_kurtosis=data.frame(kurtosis=apply(y_post_lddp,FUN = moments::kurtosis,MARGIN = 1))
original_kurtosis=moments::kurtosis(data$y)
theme_set(theme_bw())
ggplot(data_kurtosis, aes(x=kurtosis)) + 
  geom_histogram( fill="slategray4", color="slategray4")+
  geom_vline(xintercept = original_kurtosis, color="lightblue", linewidth=1.6)+
  xlab("Kurtosis of y")+ylab("Counts")+  theme(     plot.title = element_text(size = 20, face = "bold", hjust = 0.5),     axis.title = element_text(size = 24),     axis.text = element_text(size = 15)   )
dev.off()

png(filename = "D:/PhD_study2_desktop/plots/example_plots/example_D_postcheck_skewness.png",
    width = 609*5, height = 469*5,res = 72*5)
data_skewness=data.frame(skewness=apply(y_post_lddp,FUN = moments::skewness,MARGIN = 1))
original_skewness=moments::skewness(data$y)
theme_set(theme_bw())
ggplot(data_skewness, aes(x=skewness)) + 
  geom_histogram(fill="slategray4", color="slategray4")+
  geom_vline(xintercept = original_skewness, color="lightblue", linewidth=1.6)+
  xlab("Skewness of y")+ylab("Counts")+  theme(     plot.title = element_text(size = 20, face = "bold", hjust = 0.5),     axis.title = element_text(size = 24),     axis.text = element_text(size = 15)   )
dev.off()

png(filename = "D:/PhD_study2_desktop/plots/example_plots/example_D_postcheck_min.png",
    width = 609*5, height = 469*5,res = 72*5)
data_min=data.frame(min=apply(y_post_lddp,FUN = min,MARGIN = 1))
original_min=min(data$y)
theme_set(theme_bw())
ggplot(data_min, aes(x=min)) + 
  geom_histogram( fill="slategray4", color="slategray4")+
  geom_vline(xintercept = original_min, color="lightblue", linewidth=1.6)+
  xlab("Minmum value of y")+ylab("Counts")+  theme(     plot.title = element_text(size = 20, face = "bold", hjust = 0.5),     axis.title = element_text(size = 24),     axis.text = element_text(size = 15)   )
dev.off()

png(filename = "D:/PhD_study2_desktop/plots/example_plots/example_D_postcheck_max.png",
    width = 609*5, height = 469*5,res = 72*5)
data_max=data.frame(max=apply(y_post_lddp,FUN = max,MARGIN = 1))
original_max=max(data$y)
theme_set(theme_bw())
ggplot(data_max, aes(x=max)) + 
  geom_histogram( fill="slategray4", color="slategray4")+
  geom_vline(xintercept = original_max, color="lightblue", linewidth=1.6)+
  xlab("Maximum value of y")+ylab("Counts")+  theme(     plot.title = element_text(size = 20, face = "bold", hjust = 0.5),     axis.title = element_text(size = 24),     axis.text = element_text(size = 15)   )
dev.off()

png(filename = "D:/PhD_study2_desktop/plots/example_plots/example_D_postcheck_sd.png",
    width = 609*5, height = 469*5,res = 72*5)
data_sd=data.frame(sd=apply(y_post_lddp,FUN = sd,MARGIN = 1))
original_sd=sd(data$y)
theme_set(theme_bw())
ggplot(data_sd, aes(x=sd)) + 
  geom_histogram( fill="slategray4", color="slategray4")+
  geom_vline(xintercept = original_sd, color="lightblue", linewidth=1.6)+
  xlab("Standard deviation of y")+ylab("Counts")+  theme(     plot.title = element_text(size = 20, face = "bold", hjust = 0.5),     axis.title = element_text(size = 24),     axis.text = element_text(size = 15)   )
dev.off()

##conditional density heatmap
ygrid=seq(-3, 8,length=100)
xgrid=seq(-3, 11,length=60)

q <- apply(x[,2:20], 2, quantile, probs = 0.25)
q_mat <- matrix(q, 
                nrow = length(xgrid), 
                ncol = length(q), 
                byrow = TRUE)
X_pred1=cbind(1,xgrid ,q_mat)

q <- apply(x[,2:20], 2, quantile, probs = 0.5)
q_mat <- matrix(q, 
                nrow = length(xgrid), 
                ncol = length(q), 
                byrow = TRUE)
X_pred2=cbind(1,xgrid ,q_mat)

q <- apply(x[,2:20], 2, quantile, probs = 0.75)
q_mat <- matrix(q, 
                nrow = length(xgrid), 
                ncol = length(q), 
                byrow = TRUE)
X_pred3=cbind(1,xgrid ,q_mat)

ptm=proc.time()
dense_q1=density_lddp(fit_lddp=fit_lddp,X=X_pred1,grid=ygrid)
dense_q2=density_lddp(fit_lddp=fit_lddp,X=X_pred2,grid=ygrid)
dense_q3=density_lddp(fit_lddp=fit_lddp,X=X_pred3,grid=ygrid)
proc.time()-ptm


p_weight_grid=tau1*exp(-0.5*tau1^2*(xgrid-mu1)^2)/(tau1*exp(-0.5*tau1^2*(xgrid-mu1)^2)
                                              +tau2*exp(-0.5*tau2^2*(xgrid-mu2)^2))

true_dense_grid=matrix(0,nrow = length(ygrid),ncol=length(xgrid))
for(i in 1:length(ygrid)){
  for(j in 1: length(xgrid)){
    true_dense_grid[i,j]=mixAK::dMVNmixture(x=ygrid[i],
                                            weight = c(p_weight_grid[j],1-p_weight_grid[j]),
                                            mean = c(beta_10+beta_11*xgrid[j],beta_20+beta_21*xgrid[j]),
                                            Sigma = c(sigma1^2,sigma2^2))
  }
}

mean_dense_q1=apply(dense_q1,FUN = mean,MARGIN = c(2,3))
mean_dense_q2=apply(dense_q2,FUN = mean,MARGIN = c(2,3))
mean_dense_q3=apply(dense_q3,FUN = mean,MARGIN = c(2,3))

library(dplyr)
heat_df1 <- expand.grid(y = ygrid, x = xgrid) %>%
  dplyr::mutate(density = as.vector(mean_dense_q1))

heat_df2 <- expand.grid(y = ygrid, x = xgrid) %>%
  dplyr::mutate(density = as.vector(mean_dense_q2))

heat_df3 <- expand.grid(y = ygrid, x = xgrid) %>%
  dplyr::mutate(density = as.vector(mean_dense_q3))

heat_dftrue <- expand.grid(y = ygrid, x = xgrid) %>%
  dplyr::mutate(density = as.vector(true_dense_grid))

png("D:/PhD_study2_desktop/plots/example_plots/example_D_heat_q1.png",
    width = 500*5, height = 334*5,res = 72*5)
theme_set(theme_bw())
ggplot(heat_df1, aes(x, y, fill = density)) +
  coord_cartesian(xlim = c(-3, 11), ylim = c(-3, 8))+
  geom_raster(interpolate = TRUE) + scale_fill_gradient(low = "white", high = "blue")+
  labs(x = "x1", y = "y",title="LDDP; conditioned on x2:x19 at the 25th quantile") +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    
    legend.key = element_blank(),panel.grid=element_blank()
  )
dev.off()

png("D:/PhD_study2_desktop/plots/example_plots/example_D_heat_q2.png",
    width = 500*5, height = 334*5,res = 72*5)
theme_set(theme_bw())
ggplot(heat_df2, aes(x, y, fill = density)) +
  coord_cartesian(xlim = c(-3, 11), ylim = c(-3, 8))+
  geom_raster(interpolate = TRUE) + scale_fill_gradient(low = "white", high = "blue")+
  labs(x = "x1", y = "y",title="LDDP; conditioned on x2:x19 at the 50th quantile") +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    
    legend.key = element_blank(),panel.grid=element_blank()
  )
dev.off()

png("D:/PhD_study2_desktop/plots/example_plots/example_D_heat_q3.png",
    width = 500*5, height = 334*5,res = 72*5)
theme_set(theme_bw())
ggplot(heat_df3, aes(x, y, fill = density)) +
  coord_cartesian(xlim = c(-3, 11), ylim = c(-3, 8))+
  geom_raster(interpolate = TRUE) + scale_fill_gradient(low = "white", high = "blue")+
  labs(x = "x1", y = "y",title="LDDP; conditioned on x2:x19 at the 75th quantile") +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    
    legend.key = element_blank(),panel.grid=element_blank()
  )
dev.off()

png("D:/PhD_study2_desktop/plots/example_plots/example_D_heat_true.png",
    width = 500*5, height = 334*5,res = 72*5)
theme_set(theme_bw())
ggplot(heat_dftrue, aes(x, y, fill = density)) +
  coord_cartesian(xlim = c(-3, 11), ylim = c(-3, 8))+
  geom_raster(interpolate = TRUE) + scale_fill_gradient(low = "white", high = "blue")+
  labs(x = "x1", y = "y",title="Truth") +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    
    legend.key = element_blank(),panel.grid=element_blank()
  )
dev.off()



