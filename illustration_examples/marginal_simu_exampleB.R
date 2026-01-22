library(ggplot2)
require(mcclust.ext)
library(UNL.est)

source("D:/PhD_study2_desktop/lddp_functions.R")
#example b
set.seed(111)

set.seed(151)
n=600
x1=runif(n,min=-2,max=2)
x2=runif(n,min=-2,max=2)
eps=rnorm(n=n,mean=0,sd=0.1)
means=ifelse(sin(x1*x2*pi/2)<=0,1,-1)
y=means+eps
data=data.frame(cbind(x1,x2,eps,means,y))

ggplot(data)+geom_point(aes(x=x1,y=x2,col=y))+theme_bw()

L=10
set.seed(173)
prior=prior_dpm(list(y=as.matrix(y),y_cate=NULL), L=L, K=L%/%2, nstart = 5)
mcmc <- list(nsave = 10000, nburn = 10000, nskip = 1)
ptm=proc.time()
set.seed(189)
res <- dpm_MN_Mcate(y_all = list(y=as.matrix(y),y_cate=NULL), prior = prior$prior_full, mcmc = mcmc, standardise = FALSE)
proc.time()-ptm

psm_dpm <- comp.psm(res$z)
output_vi_dpm <- minVI(psm_dpm, res$z)

ggplot() +
  geom_point(aes(x = x1, y = x2, color = as.factor(output_vi_dpm$cl))) +
  theme_bw() +
  labs( x = "x1", y = "x2", color = "Cluster") 

png("D:/PhD_study2_desktop/plots/example_plots/example_B_fit.png",
    width = 609*5, height = 469*5,res = 72*5)
theme_set(theme_bw())
ggplot(data)+geom_point(aes(x=x1,y=x2,col=as.character(output_vi_dpm$cl)))+labs(col="cluster")+
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


x_c1 <- as.matrix(data[ind1,c(1,2)]);samples1=list(y=x_c1,y_cate=NULL)
x_c2 <- as.matrix(data[ind2,c(1,2)]);samples2=list(y=x_c2,y_cate=NULL)



##prior for group 1,2,3
set.seed(123)
prior1=prior_dpm(samples1, L=L, K=L%/%2, nstart = 5)
set.seed(123)
prior2=prior_dpm(samples2, L=L, K=L%/%2, nstart = 5)


ptm=proc.time()
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
res_list=list(res_1,res_2)
#UNL for joint density
ptm=proc.time()
set.seed(125)
UNL_imp_full=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_list,
                           n_imp=5000,continuous_slct_index=c(1,2),
                           cate_slct_index=NULL,future.seed = TRUE)
#UNL for marginal density
set.seed(125)
UNL_imp_m1=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_list,
                         n_imp=5000,continuous_slct_index=c(1),
                         cate_slct_index=NULL,future.seed = TRUE)
set.seed(125)
UNL_imp_m2=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_list,
                         n_imp=5000,continuous_slct_index=c(2),
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




png("D:/PhD_study2_desktop/plots/example_plots/example_B_UNLhist.png",
    width = 609*5, height = 469*5,res = 72*5)
cols <- c(
  rgb(0, 0, 1, 0.2),
  rgb(0, 1, 0, 0.2),
  rgb(1, 0, 0, 0.2)
)
full  <- tranform_list(UNL_imp_full)$unl
m1 <- tranform_list(UNL_imp_m1)$unl
m2 <- tranform_list(UNL_imp_m2)$unl

# 合并为 tidy 格式
df <- data.frame(
  unl = c(full, m1, m2),
  group = factor(rep(
    c("x1 and x2",
      "x1",
      "x2"),
    times = c(length(full), length(m1), length(m2))
  ),level=c("x1 and x2",
            "x1",
            "x2"))
)
mycols <- setNames(cols[1:3], c("x1 and x2",
                                "x1",
                                "x2"))
# 绘图
ggplot(df, aes(x = unl, fill = group)) +
  geom_histogram(aes(y = ..density..),        # freq = FALSE -> density
                 position = "identity",       # 叠加而非堆叠
                 alpha = 0.45,                # 透明度，便于比较
                 bins = 90,                   # 可调整或改为 binwidth = ...
                 color = "black") +          # 柱子边框
  scale_fill_manual(values = mycols,
                    name = NULL) +
  coord_cartesian(xlim = c(1, 2), ylim = c(0, 18)) +   # 保留你原来的 x/y 范围
  labs(x = "UNL", y = "Density") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 20),
    legend.position = c(0.98, 0.98),                  # 图内右上角（模拟 base::legend("topright")）
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank()
  )
dev.off()