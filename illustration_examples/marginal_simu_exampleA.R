library(ggplot2)
require(mcclust.ext)
library(UNL.est)

source("D:/PhD_study2_desktop/lddp_functions.R")
#example a
set.seed(111)
n=600
x=runif(n,min=-3,max=3)
eps=rnorm(n=n,mean=0,sd=0.1)
means=ifelse(x<=-1,2,ifelse(x>=1,-5,0))
y=means+eps
sign=ifelse(x<=-1,1,ifelse(x>=1,3,2))

data=data.frame(cbind(x,eps,means,y,sign))

data=data.frame(cbind(x,eps,means,y))
png("D:/PhD_study2_desktop/plots/example_plots/example_A1.png",
    width = 500*5, height = 334*5,res = 72*5)
ggplot(data)+geom_point(aes(x=x,y=y,col=as.factor(sign)))+labs(col="")+theme_bw()
dev.off()

# L=20
# prior <- list(m0 = 0, S0 = 10, a = 2, b = 0.5, aalpha = 2,  balpha = 2, L = 20)
# # prior <- list(m0 = apply(y,FUN = mean,MARGIN = 2,na.rm=TRUE), L0 = var(y),
# #               a = 2, b = 0.5, alpha = 1, L = 10)
# mcmc <- list(nsave = 10000, nburn = 2000, nskip = 1)
# 
# res <- dpm(y = y, prior = prior, mcmc = mcmc, standardise = TRUE)


L=10
set.seed(123)
prior=prior_dpm(list(y=as.matrix(y),y_cate=NULL), L=L, K=L%/%2, nstart = 5)
mcmc <- list(nsave = 10000, nburn = 10000, nskip = 1)
ptm=proc.time()
set.seed(159)
res <- dpm_MN_Mcate(y_all = list(y=as.matrix(y),y_cate=NULL), prior = prior$prior_full, mcmc = mcmc, standardise = FALSE)
proc.time()-ptm

psm_dpm <- comp.psm(res$z)
output_vi_dpm <- minVI(psm_dpm, res$z)

png("D:/PhD_study2_desktop/plots/example_plots/example_A_fit.png",
    width = 609*5, height = 469*5,res = 72*5)
theme_set(theme_bw())
ggplot(data)+geom_point(aes(x=x,y=y,col=as.character(output_vi_dpm$cl)))+labs(col="cluster")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        plot.title.position = "plot",            # 可选：使标题相对于整个绘图区域定位
        plot.subtitle = element_text(hjust = 0.5),
        text =  element_text(size = 20),
        legend.text = element_text(size = 20)
  )
dev.off()

ind1 <- which(output_vi_dpm$cl == 1)
ind2 <- which(output_vi_dpm$cl == 2)
ind3 <- which(output_vi_dpm$cl == 3)

x_c1 <- as.matrix(x[ind1]);samples1=list(y=x_c1,y_cate=NULL)
x_c2 <- as.matrix(x[ind2]);samples2=list(y=x_c2,y_cate=NULL)
x_c3 <- as.matrix(x[ind3]);samples3=list(y=x_c3,y_cate=NULL)

##prior for group 1,2,3
set.seed(123)
prior1=prior_dpm(samples1, L=L, K=L%/%2, nstart = 5)
set.seed(123)
prior2=prior_dpm(samples2, L=L, K=L%/%2, nstart = 5)
set.seed(123)
prior3=prior_dpm(samples3, L=L, K=L%/%2, nstart = 5)

ptm=proc.time()
set.seed(159)
res_1 <- dpm_MN_Mcate(y_all = samples1, prior = prior1$prior_full, mcmc = mcmc, standardise = FALSE) 
set.seed(159)
res_2 <- dpm_MN_Mcate(y_all = samples2, prior = prior2$prior_full, mcmc = mcmc, standardise = FALSE) 
set.seed(159)
res_3 <- dpm_MN_Mcate(y_all = samples3, prior = prior3$prior_full, mcmc = mcmc, standardise = FALSE) 
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
set.seed(123)
UNL_imp_full=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_list,
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

# png("D:/PhD_study2_desktop/plots/example_plots/example_A1_UNLhist.png",
#     width = 500*5, height = 334*5,res = 72*5)
# hist(tranform_list(UNL_imp_full)$unl,xlim=c(1,3),col = rgb(1, 0, 0, 0.2),freq = FALSE,xlab = "UNL",main = "3 clusters (marginal); 3 clusters (conditional)")
# hist(tranform_list(UNL_imp_full)$unl,xlim=c(1,2),col = rgb(0, 0, 1, 0.2), add = TRUE,freq = FALSE)
# dev.off()

png("D:/PhD_study2_desktop/plots/example_plots/example_A_UNLhist.png",
    width = 609*5, height = 469*5,res = 72*5)
cols <- c(
  rgb(0, 0, 1, 0.2),
  rgb(0, 1, 0, 0.2),
  rgb(1, 0, 0, 0.2)
)

df <- data.frame(
  unl = tranform_list(UNL_imp_full)$unl
)


# 绘图（density：y = ..density..）
ggplot(df, aes(x = unl, )) +
  geom_histogram(aes(y = after_stat(density)),
                 position = "identity",    # 叠加显示
                 alpha = 0.45,             # 透明度，便于比较
                 bins = 15,          # 可根据需要调整 binwidth 或用 bins = 30
                 color = "black",fill = cols[1]) +        # 柱子边框
  coord_cartesian(xlim = c(1, 3)) +   # 保留原来你给的 xlim/ylim
  labs(x = "UNL", y = "Density") +
  theme_minimal() +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 20),
    legend.position = c(0.17, 0.98),                  # 图内右上角（模拟 base::legend("topright")）
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank()
  )
dev.off()

