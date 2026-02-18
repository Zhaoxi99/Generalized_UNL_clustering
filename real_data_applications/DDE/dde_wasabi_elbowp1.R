library(dplyr)
library(lubridate)
library(ggplot2)
require(mcclust.ext)
library(UNL.est)
library(corrplot)
library(LSBP)
library(WASABI)

source("functions/lddp_functions.R")

load("dde.RData")

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


set.seed(176)
fit_lddp <- lddp_moves(y = data$GAD,
                       X = X,
                       prior = prior, 
                       mcmc = mcmc,
                       standardise = FALSE)


psm_lddp <- comp.psm(fit_lddp$z)
output_vi_lddp <- minVI(psm_lddp, fit_lddp$z)

plotpsm(psm_lddp)

ptm=proc.time()
set.seed(129)
out_elbow <- elbow(fit_lddp$z, L_max = 6, psm = psm_lddp, ncores = parallel::detectCores()-2,
                   multi.start = 4, mini.batch = 300,
                   method.init = "++", method = "salso")
proc.time()-ptm

plot(out_elbow$wass_vec, type = "b", ylab = "Wass distance", xlab = "Number of particles")

output_WASABI <- out_elbow$output_list[[4]]

ggsummary(output_WASABI)

table(output_WASABI$particles[1,])
table(output_WASABI$particles[2,])
table(output_vi_lddp$cl)


###
ind1 <- which(output_WASABI$particles[1,] == 1)
ind2 <- which(output_WASABI$particles[1,] == 2)

sum(dde$GAD[ind2]>=42)/length(ind2)

sum(dde$GAD[ind1]<42)/length(ind1)
##regression in each cluster
fit_lm1 <- lm(GAD ~DDE,data = dde[ind1,])
summary(fit_lm1)

fit_lm2 <- lm(GAD ~DDE,data = dde[ind2,])
summary(fit_lm2)




theme_set(theme_bw())
ggplot(dde)+geom_point(aes(x=DDE,y=GAD,col=as.factor(output_WASABI$particles[1,])),alpha=.5, cex=.5)+
  xlab("DDE (mg/L)") + ylab("Gestational age at delivery (in weeks)") +labs(col="cluster")+
  geom_hline(yintercept = 37, color = "red", linewidth = 0.45,linetype = "dashed")+
  geom_hline(yintercept = 42, color = "blue", linewidth = 0.45,linetype = "dashed")+
  geom_abline(intercept = coef(fit_lm1)[1], slope = coef(fit_lm2)[2], 
              color = "red", linewidth = 0.45)+
  geom_abline(intercept = coef(fit_lm2)[1], slope = coef(fit_lm2)[2], 
              color = "blue", linewidth = 0.45)+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        plot.title.position = "plot",            # 可选：使标题相对于整个绘图区域定位
        plot.subtitle = element_text(hjust = 0.5),
        text =  element_text(size = 15),
        legend.text = element_text(size = 15)
  )


##

x_c1 <- as.matrix(X[ind1,c(2)]);samples1=list(y=x_c1,y_cate=NULL)
x_c2 <- as.matrix(X[ind2,c(2)]);samples2=list(y=x_c2,y_cate=NULL)


##prior for group 1,2
L=10
set.seed(158)
prior1=prior_dpm(samples1, L=L, K=L%/%2, nstart = 5,categories = NULL)
set.seed(158)
prior2=prior_dpm(samples2, L=L, K=L%/%2, nstart = 5,categories = NULL)



ptm=proc.time()
set.seed(159)
res_1 <- dpm_MN_Mcate_shrink(y_all = samples1, prior = prior1$prior_shrink, mcmc = mcmc, standardise = FALSE) 
set.seed(159)
res_2 <- dpm_MN_Mcate_shrink(y_all = samples2, prior = prior2$prior_shrink, mcmc = mcmc, standardise = FALSE) 
proc.time()-ptm
##calculate UNL
library(future)
library(future.apply)
nsave_list <- list()
for (i in 1:mcmc$nsave) {
  nsave_list[[as.character(i)]] <- i
}

plan(multisession, workers = detectCores()-2)
res_list=list(res_1,res_2)
#UNL for joint density
ptm=proc.time()
set.seed(145)
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




cols <- c(
  rgb(0, 0, 1, 0.2),
  rgb(0, 1, 0, 0.2),
  rgb(1, 0, 0, 0.2)
)


df <- data.frame(
  unl = tranform_list(UNL_imp_lddp)$unl
)

ggplot(df, aes(x = unl)) +
  geom_histogram(aes(y = ..density..),        
                 position = "identity",       
                 alpha = 0.45,                
                 bins = 30,                   
                 color = "black",fill=rgb(0, 0, 1, 0.2)) +          
  coord_cartesian(xlim = c(1, 2), ylim = c(0, 12)) +   
  labs(x = "UNL", y = "Density") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    legend.position = c(0.98, 0.98),                 
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank()
  )

summary(tranform_list(UNL_imp_lddp)$unl)