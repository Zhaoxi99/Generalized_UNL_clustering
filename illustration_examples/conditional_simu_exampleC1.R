library(ggplot2)
require(mcclust.ext)
library(UNL.est)

source("functions/lddp_functions.R")
#example c1
set.seed(111)
n=800
x=runif(n,min=-3,max=3)
eps=rnorm(n=n,mean=0,sd=0.4)
set.seed(121)
sign=sample(c(0,1,2,3),size=n,replace = TRUE)
y=rep(0,n)
y[sign==0]=2*x[sign==0]+eps[sign==0]
y[sign==1]=-2*x[sign==1]+eps[sign==1]
y[sign==2]=12*x[sign==2]+eps[sign==2]+80
y[sign==3]=-12*x[sign==3]+eps[sign==3]+80

data=data.frame(cbind(x,eps,y,sign))
x_cate=ifelse(sign==0|sign==1,1,2)
data=data.frame(cbind(x,eps,y,sign,x_cate))
data$x_cate=as.factor(x_cate);data$sign1=as.factor(sign+1)
#ggplot(data)+geom_point(aes(x=x,y=y,col=sign1))+labs(col="")+theme_bw()
cols <- c("1" = "#E41A1C",  # red
          "2" = "#377EB8",  # blue
          "3" = "#4DAF4A",  # green
          "4" = "#984EA3")  # purple

ggplot(data)+geom_point(aes(x=x,y=y,col=sign1))+
  scale_colour_manual(values = cols, drop = FALSE) +labs(col="")+theme_bw()+
  coord_cartesian(xlim = c(-3, 3), ylim = c(-65, 150))


##conditional fit
fit_lm <- lm(y ~x+x_cate,data = data)
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

ptm=proc.time()
set.seed(156)
fit_lddp <- lddp_moves(y = y,
                       X = X,
                       prior = prior, 
                       mcmc = mcmc,
                       standardise = FALSE)

proc.time()-ptm

psm_lddp <- comp.psm(fit_lddp$z)
output_vi_lddp <- minVI(psm_lddp, fit_lddp$z)
#output_vi_lddp$cl

ptm=proc.time()
y_post_lddp=predic_lddp(fit_lddp =fit_lddp,X=X)
proc.time()-ptm

ggplot() +
  geom_point(aes(x = x, y = y, color = as.factor(output_vi_lddp$cl))) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-65, 150))+
  theme_bw() +
  labs( x = "x", y = "y", color = "") 


ggplot() +
  geom_point(aes(x = x, y = y, color = as.factor(output_vi_lddp$cl))) +
  scale_colour_manual(values = cols, drop = FALSE)+
  coord_cartesian(xlim = c(-3, 3), ylim = c(-65, 150))+
  labs( x = expression(x^c), y = "y", color = "cluster")+
  theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        plot.title.position = "plot",            # 可选：使标题相对于整个绘图区域定位
        plot.subtitle = element_text(hjust = 0.5),
        text =  element_text(size = 20),
        legend.text = element_text(size = 20)
  )


ggplot() +
  geom_point(aes(x = x[x_cate==2], y = y[x_cate==2], color = as.factor(output_vi_lddp$cl[x_cate==2]))) +
  scale_colour_manual(values = cols, drop = FALSE)+
  coord_cartesian(xlim = c(-3, 3), ylim = c(-65, 150))+
  labs( x = "x", y = "y", color = "cluster") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 23, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 23),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )

#summary(y_post_lddp)



ind1 <- which(output_vi_lddp$cl == 1)
ind2 <- which(output_vi_lddp$cl == 2)
ind3 <- which(output_vi_lddp$cl == 3)
ind4 <- which(output_vi_lddp$cl == 4)


x_c1 <- as.matrix(x[ind1]);x_cate_c1 <- as.matrix(x_cate[ind1]);samples1=list(y=x_c1,y_cate=x_cate_c1)
x_c2 <- as.matrix(x[ind2]);x_cate_c2 <- as.matrix(x_cate[ind2]);samples2=list(y=x_c2,y_cate=x_cate_c2)
x_c3 <- as.matrix(x[ind3]);x_cate_c3 <- as.matrix(x_cate[ind3]);samples3=list(y=x_c3,y_cate=x_cate_c3)
x_c4 <- as.matrix(x[ind4]);x_cate_c4 <- as.matrix(x_cate[ind4]);samples4=list(y=x_c4,y_cate=x_cate_c4)


##prior for group 1,2,3,4
L=10
set.seed(124)
prior1=prior_dpm(samples1, L=L, K=L%/%2, nstart = 5,categories = 2)
set.seed(124)
prior2=prior_dpm(samples2, L=L, K=L%/%2, nstart = 5,categories = 2)
set.seed(124)
prior3=prior_dpm(samples3, L=L, K=L%/%2, nstart = 5,categories = 2)
set.seed(124)
prior4=prior_dpm(samples4, L=L, K=L%/%2, nstart = 5,categories = 2)



ptm=proc.time()
set.seed(159)
res_1 <- dpm_MN_Mcate(y_all = samples1, prior = prior1$prior_full, mcmc = mcmc, standardise = FALSE) 
set.seed(159)
res_2 <- dpm_MN_Mcate(y_all = samples2, prior = prior2$prior_full, mcmc = mcmc, standardise = FALSE) 
set.seed(159)
res_3 <- dpm_MN_Mcate(y_all = samples3, prior = prior3$prior_full, mcmc = mcmc, standardise = FALSE) 
set.seed(159)
res_4 <- dpm_MN_Mcate(y_all = samples4, prior = prior4$prior_full, mcmc = mcmc, standardise = FALSE) 
proc.time()-ptm
##calculate UNL
library(future)
library(future.apply)
nsave_list <- list()
for (i in 1:mcmc$nsave) {
  nsave_list[[as.character(i)]] <- i
}

plan(multisession, workers = detectCores()-2)
res_list=list(res_1,res_2,res_3,res_4)
#UNL for joint density
ptm=proc.time()
set.seed(123)
UNL_imp_lddp_full=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_list,
                                n_imp=5000,continuous_slct_index=c(1),
                                cate_slct_index=c(1),future.seed = TRUE)
set.seed(123)
UNL_imp_lddp_conti=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_list,
                                 n_imp=5000,continuous_slct_index=c(1),
                                 cate_slct_index=NULL,future.seed = TRUE)
set.seed(123)
UNL_imp_lddp_cate=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_list,
                                n_imp=5000,continuous_slct_index=NULL,
                                cate_slct_index=c(1),future.seed = TRUE)
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
full  <- tranform_list(UNL_imp_lddp_full)$unl
m1 <- tranform_list(UNL_imp_lddp_cate)$unl
m2 <- tranform_list(UNL_imp_lddp_conti)$unl

df <- data.frame(
  unl = c(full, m1, m2),
  group = factor(rep(
    c("x^c and x^d",
      "x^d",
      "x^c"),
    times = c(length(full), length(m1), length(m2))
  ),level=c("x^c and x^d",
            "x^d",
            "x^c"))
)
mycols <- setNames(cols[1:3], c("x^c and x^d",
                                "x^d",
                                "x^c"))
ggplot(df, aes(x = unl, fill = group)) +
  geom_histogram(aes(y = after_stat(density)),      
                 position = "identity",    
                 alpha = 0.45,              
                 #binwidth =0.02,
                 bins = 80,                  
                 color = "black") +     
  scale_fill_manual(values = mycols,
                    name = NULL) +
  coord_cartesian(xlim = c(1, 4), ylim = c(0, 10)) +  
  labs(x = "UNL", y = "Density") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 20),
    legend.position = c(0.98, 0.98),                  
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank()
  )


ygrid=seq(-65, 150,by=1)
xgrid=seq(-3, 3,by=0.1)
X_predcate1=cbind(1,xgrid,0)
X_predcate2=cbind(1,xgrid,1)
ptm=proc.time()
dense_cate1=density_lddp(fit_lddp=fit_lddp,X=X_predcate1,grid=ygrid)
dense_cate2=density_lddp(fit_lddp=fit_lddp,X=X_predcate2,grid=ygrid)
proc.time()-ptm

true_dense_cate1=matrix(0,nrow = length(ygrid),ncol=length(xgrid))
for(i in 1:length(ygrid)){
  for(j in 1: length(xgrid)){
    true_dense_cate1[i,j]=mixAK::dMVNmixture(x=ygrid[i],weight = c(0.5,0.5),mean = c(2*xgrid[j],-2*xgrid[j]),Sigma = c(0.4^2,0.4^2))
  }
}

true_dense_cate2=matrix(0,nrow = length(ygrid),ncol=length(xgrid))
for(i in 1:length(ygrid)){
  for(j in 1: length(xgrid)){
    true_dense_cate2[i,j]=mixAK::dMVNmixture(x=ygrid[i],weight = c(0.5,0.5),mean = c(-12*xgrid[j]+80,12*xgrid[j]+80),Sigma = c(0.4^2,0.4^2))
  }
}

med_dense_cate1=apply(dense_cate1,FUN = mean,MARGIN = c(2,3))
med_dense_cate2=apply(dense_cate2,FUN = mean,MARGIN = c(2,3))

library(dplyr)
heat_df1 <- expand.grid(y = ygrid, x = xgrid) %>%
  dplyr::mutate(density = as.vector(med_dense_cate1))

heat_df2 <- expand.grid(y = ygrid, x = xgrid) %>%
  dplyr::mutate(density = as.vector(med_dense_cate2))

heat_dftrue1 <- expand.grid(y = ygrid, x = xgrid) %>%
  dplyr::mutate(density = as.vector(true_dense_cate1))

heat_dftrue2 <- expand.grid(y = ygrid, x = xgrid) %>%
  dplyr::mutate(density = as.vector(true_dense_cate2))


theme_set(theme_bw())
ggplot(heat_df1, aes(x, y, fill = density)) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-65, 150))+
  geom_raster(interpolate = TRUE) + scale_fill_gradient(low = "white", high = "blue")+
  labs(x = expression(x^c), y = "y",title="LDDP") +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    
    legend.key = element_blank(),panel.grid=element_blank()
  )



theme_set(theme_bw())
ggplot(heat_df2, aes(x, y, fill = density)) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-65, 150))+
  geom_raster(interpolate = TRUE) + scale_fill_gradient(low = "white", high = "blue")+
  labs(x = expression(x^c), y = "y",title="LDDP") +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    
    legend.key = element_blank(),panel.grid=element_blank()
  )



theme_set(theme_bw())
ggplot(heat_dftrue1, aes(x, y, fill = density)) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-65, 150))+
  geom_raster(interpolate = TRUE) + scale_fill_gradient(low = "white", high = "blue")+
  labs(x = expression(x^c), y = "y",title="Truth") +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    
    legend.key = element_blank(),panel.grid=element_blank()
  )


theme_set(theme_bw())
ggplot(heat_dftrue2, aes(x, y, fill = density)) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-65, 150))+
  geom_raster(interpolate = TRUE) + scale_fill_gradient(low = "white", high = "blue")+
  labs(x = expression(x^c), y = "y",title="Truth") +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    
    legend.key = element_blank(),panel.grid=element_blank()
  )
