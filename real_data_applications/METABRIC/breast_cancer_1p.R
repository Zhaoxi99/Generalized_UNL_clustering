library(dplyr)
library(lubridate)
library(ggplot2)
require(mcclust.ext)
library(UNL.est)
library(corrplot)
library(WASABI)

source("D:/PhD_study2_desktop/lddp_functions.R")
#source("D:/PhD_study2_desktop/multi_dpm.R")
clinical_raw=readr::read_tsv("D:/PhD_study2_desktop/new_app/metabric/brca_metabric_clinical_data.tsv", 
                             na = c("", "NA"))  
gene_raw=readr::read_table("D:/PhD_study2_desktop/new_app/metabric/mRNA expression_with_BECN1_BRCA1.txt")
CNA_raw=readr::read_table("D:/PhD_study2_desktop/new_app/metabric/cna.txt") 
data_clinical_raw=clinical_raw%>%
  dplyr::select(`Sample ID`,`Lymph nodes examined positive`,
                `Neoplasm Histologic Grade`,`Tumor Size`,
                `Nottingham prognostic index`)
data_clinical=data_clinical_raw[complete.cases(data_clinical_raw),]
data_clinical$node_stage=ifelse(data_clinical$`Lymph nodes examined positive`==0,1,ifelse(data_clinical$`Lymph nodes examined positive`<=3,2,3))

data_raw=gene_raw%>%left_join(data_clinical,
                              by=c("SAMPLE_ID"="Sample ID"))
data=data_raw[complete.cases(data_raw),]

samples=list(y=as.matrix(data_clinical$`Tumor Size`),y_cate=cbind(data_clinical$node_stage,data_clinical$`Neoplasm Histologic Grade`))

##prior
L=10
set.seed(123)
prior=prior_dpm(samples, L=L, K=4, nstart = 5,categories = c(3,3))
mcmc <- list(nsave = 10000, nburn = 10000, nskip = 1)

ptm=proc.time()
set.seed(159)
res_dpm <- dpm_MN_Mcate(y_all = samples, prior = prior$prior_full, mcmc = mcmc, standardise = FALSE) 
proc.time()-ptm

psm_dpm <- comp.psm(res_dpm$z)
output_vi_dpm <- minVI(psm_dpm, res_dpm$z)

ptm=proc.time()
set.seed(129)
out_WASABI_2p <- WASABI(res_dpm$z, psm = psm_dpm, L = 2,
                        method.init = "++", method = "salso")
proc.time()-ptm

ptm=proc.time()
set.seed(1519)
out_WASABI_3p <- WASABI(res_dpm$z, psm = psm_dpm, L = 3,
                        method.init = "++", method = "salso")
proc.time()-ptm

ggsummary(out_WASABI_2p)
ggsummary(out_WASABI_3p)

table(out_WASABI_2p$particles[1,])
table(out_WASABI_2p$particles[2,])
table(output_vi_dpm$cl)

# summarize each group
summary(data_clinical$`Tumor Size`[output_vi_dpm$cl==1])
summary(data_clinical$`Tumor Size`[output_vi_dpm$cl==2])
summary(data_clinical$`Tumor Size`[output_vi_dpm$cl==3])

mean(data_clinical$`Tumor Size`[output_vi_dpm$cl==1])
mean(data_clinical$`Tumor Size`[output_vi_dpm$cl==2])
mean(data_clinical$`Tumor Size`[output_vi_dpm$cl==3])

# sd(data$`Tumor Size`[output_vi_dpm$cl==1])^2
# sd(data$`Tumor Size`[output_vi_dpm$cl==2])^2
# sd(data$`Tumor Size`[output_vi_dpm$cl==3])^2

table(data_clinical$`Neoplasm Histologic Grade`[output_vi_dpm$cl==1])/sum(output_vi_dpm$cl==1)
table(data_clinical$`Neoplasm Histologic Grade`[output_vi_dpm$cl==2])/sum(output_vi_dpm$cl==2)
table(data_clinical$`Neoplasm Histologic Grade`[output_vi_dpm$cl==3])/sum(output_vi_dpm$cl==3)

table(data_clinical$node_stage[output_vi_dpm$cl==1])/sum(output_vi_dpm$cl==1)
table(data_clinical$node_stage[output_vi_dpm$cl==2])/sum(output_vi_dpm$cl==2)
table(data_clinical$node_stage[output_vi_dpm$cl==3])/sum(output_vi_dpm$cl==3)

summary(data_clinical$`Nottingham prognostic index`[output_vi_dpm$cl==1])
summary(data_clinical$`Nottingham prognostic index`[output_vi_dpm$cl==2])
summary(data_clinical$`Nottingham prognostic index`[output_vi_dpm$cl==3])

data_plot=data_clinical%>%mutate(group=output_vi_dpm$cl)

png("D:/PhD_study2_desktop/plots/application_plots/metabric_plots/1p/nodestage_hist.png",
    width = 500*5, height = 334*5,res = 72*5)
df_prop <- data.frame(
  grade = as.character(data_clinical$node_stage),
  cluster = as.character(output_vi_dpm$cl)
) %>%
  filter(!is.na(grade) & grade != "") %>%       # 根据需要去掉缺失/空值
  group_by(cluster, grade) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# 2. 保证等级顺序和 cluster 顺序（可按需调整）
df_prop$grade <- factor(df_prop$grade, levels = c("1","2","3"))
df_prop$cluster <- factor(df_prop$cluster, levels = c("1","2","3"))  # 如果你有不同顺序，改这里

# 3. 颜色向量（如果你已有 cols，使用它；否则用示例）
cols <- c("#377EB8","#4DAF4A","#E41A1C")

mycols <- setNames(cols[1:3], levels(df_prop$cluster))

# 4. 绘图：以 grade 为 x，cluster 用颜色区分（每个 grade 里并列比较不同 cluster）
ggplot(df_prop, aes(x = grade, y = prop, fill = cluster)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8, color = "black") +
  scale_fill_manual(values = mycols, name = "Cluster") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "Lymph node stage", y = "Proportion", 
       title = "Proportion of lymph node stage by cluster") +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1)),
            position = position_dodge(width = 0.9), vjust = -0.25, size = 3) +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))
dev.off()

png("D:/PhD_study2_desktop/plots/application_plots/metabric_plots/1p/tumorgrade_hist.png",
    width = 500*5, height = 334*5,res = 72*5)
df_prop <- data.frame(
  grade = as.character(data_clinical$`Neoplasm Histologic Grade`),
  cluster = as.character(output_vi_dpm$cl)
) %>%
  filter(!is.na(grade) & grade != "") %>%       # 根据需要去掉缺失/空值
  group_by(cluster, grade) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# 2. 保证等级顺序和 cluster 顺序（可按需调整）
df_prop$grade <- factor(df_prop$grade, levels = c("1","2","3"))
df_prop$cluster <- factor(df_prop$cluster, levels = c("1","2","3"))  # 如果你有不同顺序，改这里

# 3. 颜色向量（如果你已有 cols，使用它；否则用示例）
cols <- c("#377EB8","#4DAF4A","#E41A1C")

mycols <- setNames(cols[1:3], levels(df_prop$cluster))

# 4. 绘图：以 grade 为 x，cluster 用颜色区分（每个 grade 里并列比较不同 cluster）
ggplot(df_prop, aes(x = grade, y = prop, fill = cluster)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8, color = "black") +
  scale_fill_manual(values = mycols, name = "Cluster") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "Neoplasm histologic grade", y = "Proportion", 
       title = "Proportion of tumor histologic grade by cluster") +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1)),
            position = position_dodge(width = 0.9), vjust = -0.25, size = 3) +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),     axis.title = element_text(size = 20),     axis.text = element_text(size = 15))
dev.off()

png("D:/PhD_study2_desktop/plots/application_plots/metabric_plots/1p/tumorsize_hist.png",
    width = 500*5, height = 334*5,res = 72*5)

cols <- c("#377EB8","#4DAF4A","#E41A1C")
ggplot(data_plot, aes(x = `Tumor Size`, fill = as.factor(group))) +
  geom_histogram(aes(y = ..density..),        # freq = FALSE -> density
                 position = "identity",       # 叠加而非堆叠
                 alpha = 0.45,                # 透明度，便于比较
                 bins = 30,                   # 可调整或改为 binwidth = ...
                 color = "black") +          # 柱子边框
  scale_fill_manual(values = cols) +labs(fill="cluster")+
  #scale_x_continuous(breaks = 1:7) +
  labs(x = "Tumor Size(mm)", y = "Density",title = "Histogram of tumor size by cluster") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),           # 标题居中
    #legend.position = c(0.02, 0.98),                  # 图内右上角（模拟 base::legend("topright")）
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank()
  )
dev.off()

png("D:/PhD_study2_desktop/plots/application_plots/metabric_plots/1p/npi_hist.png",
    width = 500*5, height = 334*5,res = 72*5)
cols <- c("#377EB8","#4DAF4A","#E41A1C")

ggplot(data_plot, aes(x = `Nottingham prognostic index`, fill = as.factor(group))) +
  geom_histogram(aes(y = ..density..),        # freq = FALSE -> density
                 position = "identity",       # 叠加而非堆叠
                 alpha = 0.45,                # 透明度，便于比较
                 bins = 30,                   # 可调整或改为 binwidth = ...
                 color = "black") +          # 柱子边框
  scale_fill_manual(values = cols) +labs(fill="cluster")+
  scale_x_continuous(breaks = 1:7) +
  labs(x = "NPI", y = "Density",title = "Histogram of NPI score by cluster") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    #legend.position = c(0.02, 0.98),                  # 图内右上角（模拟 base::legend("topright")）
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank()
  )
dev.off()


##est unl
ind1 <- which(data$SAMPLE_ID %in%(data_clinical$`Sample ID`[output_vi_dpm$cl == 1]))
ind2 <- which(data$SAMPLE_ID %in%(data_clinical$`Sample ID`[output_vi_dpm$cl == 2]))
ind3 <- which(data$SAMPLE_ID %in%(data_clinical$`Sample ID`[output_vi_dpm$cl == 3]))

x_c1 <- as.matrix(data[ind1,c(3,4,5)]);samples1=list(y=x_c1,y_cate=NULL)
x_c2 <- as.matrix(data[ind2,c(3,4,5)]);samples2=list(y=x_c2,y_cate=NULL)
x_c3 <- as.matrix(data[ind3,c(3,4,5)]);samples3=list(y=x_c3,y_cate=NULL)


##prior for group 1,2,3
L=10
set.seed(123)
prior1=prior_dpm(samples1, L=L, K=L%/%2, nstart = 5,categories = NULL)
set.seed(123)
prior2=prior_dpm(samples2, L=L, K=L%/%2, nstart = 5,categories = NULL)
set.seed(123)
prior3=prior_dpm(samples3, L=L, K=L%/%2, nstart = 5,categories = NULL)
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
                           n_imp=5000,continuous_slct_index=c(1,2,3),
                           cate_slct_index=NULL,future.seed = TRUE)
set.seed(123)
UNL_imp_metsr1=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_list,
                             n_imp=5000,continuous_slct_index=c(1,2),
                             cate_slct_index=NULL,future.seed = TRUE)
set.seed(123)
UNL_imp_metsr2=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_list,
                             n_imp=5000,continuous_slct_index=c(1,3),
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


sum(tranform_list(UNL_imp_metsr1)$unl>tranform_list(UNL_imp_metsr2)$unl)


cols <- c(
  rgb(0, 0, 1, 0.2),
  rgb(0, 1, 0, 0.2),
  rgb(1, 0, 0, 0.2)
)

png("D:/PhD_study2_desktop/plots/application_plots/metabric_plots/1p/met_esr_UNL.png",
    width = 500*5, height = 334*5,res = 72*5)
full  <- tranform_list(UNL_imp_full)$unl
m1    <- tranform_list(UNL_imp_metsr1)$unl
m2    <- tranform_list(UNL_imp_metsr2)$unl

df <- data.frame(
  unl = c(full, m1, m2),
  group = factor(rep(
    c("MET, ESR1 & ESR2",
      "MET & ESR1",
      "MET & ESR2"),
    times = c(length(full), length(m1), length(m2))
  ),levels =c("MET, ESR1 & ESR2",
              "MET & ESR1",
              "MET & ESR2"))
)

mycols <- setNames(cols[1:3], c("MET, ESR1 & ESR2",
                                "MET & ESR1",
                                "MET & ESR2"))

# 绘图（density：y = ..density..）
ggplot(df, aes(x = unl, fill = group)) +
  geom_histogram(aes(y = after_stat(density)),
                 position = "identity",    # 叠加显示
                 alpha = 0.45,             # 透明度，便于比较
                 bins = 30,          # 可根据需要调整 binwidth 或用 bins = 30
                 color = "black") +        # 柱子边框
  coord_cartesian(xlim = c(1, 3), ylim = c(0, 8)) +   # 保留原来你给的 xlim/ylim
  scale_fill_manual(values = mycols,
                    name = NULL) +
  labs(x = "UNL", y = "Density") +
  theme_minimal() +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    legend.position = c(0.98, 0.98),                  # 图内右上角（模拟 base::legend("topright")）
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank()
  )
dev.off()




###BECN1 BRCA1 examination
xb_c1 <- as.matrix(data[ind1,c(7,8)]);samplesb1=list(y=xb_c1,y_cate=NULL)
xb_c2 <- as.matrix(data[ind2,c(7,8)]);samplesb2=list(y=xb_c2,y_cate=NULL)
xb_c3 <- as.matrix(data[ind3,c(7,8)]);samplesb3=list(y=xb_c3,y_cate=NULL)


##prior for group 1,2,3
L=10
set.seed(123)
priorb1=prior_dpm(samplesb1, L=L, K=L%/%2, nstart = 5,categories = NULL)
set.seed(123)
priorb2=prior_dpm(samplesb2, L=L, K=L%/%2, nstart = 5,categories = NULL)
set.seed(123)
priorb3=prior_dpm(samplesb3, L=L, K=L%/%2, nstart = 5,categories = NULL)
ptm=proc.time()
set.seed(159)
res_b1 <- dpm_MN_Mcate(y_all = samplesb1, prior = priorb1$prior_full, mcmc = mcmc, standardise = FALSE) 
set.seed(159)
res_b2 <- dpm_MN_Mcate(y_all = samplesb2, prior = priorb2$prior_full, mcmc = mcmc, standardise = FALSE) 
set.seed(159)
res_b3 <- dpm_MN_Mcate(y_all = samplesb3, prior = priorb3$prior_full, mcmc = mcmc, standardise = FALSE) 
proc.time()-ptm

##calculate UNL
library(future)
library(future.apply)
nsave_list <- list()
for (i in 1:mcmc$nsave) {
  nsave_list[[as.character(i)]] <- i
}

plan(multisession, workers = detectCores()-2)
res_listb=list(res_b1,res_b2,res_b3)
#UNL for joint density
ptm=proc.time()
set.seed(123)
UNL_imp_b_all=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_listb,
                            n_imp=5000,continuous_slct_index=c(1,2),
                            cate_slct_index=NULL,future.seed = TRUE)
set.seed(123)
UNL_imp_becn1=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_listb,
                            n_imp=5000,continuous_slct_index=c(1),
                            cate_slct_index=NULL,future.seed = TRUE)
set.seed(123)
UNL_imp_brca1=future_lapply(X=nsave_list,FUN = imp_unl3,res_list=res_listb,
                            n_imp=5000,continuous_slct_index=c(2),
                            cate_slct_index=NULL,future.seed = TRUE)
proc.time()-ptm


sum(tranform_list(UNL_imp_becn1)$unl>tranform_list(UNL_imp_brca1)$unl)


cols <- c(
  rgb(0, 0, 1, 0.2),
  rgb(0, 1, 0, 0.2),
  rgb(1, 0, 0, 0.2)
)

png("D:/PhD_study2_desktop/plots/application_plots/metabric_plots/1p/BECN_BRCA_UNL.png",
    width = 500*5, height = 334*5,res = 72*5)
full  <- tranform_list(UNL_imp_b_all)$unl
becn1 <- tranform_list(UNL_imp_becn1)$unl
brca1 <- tranform_list(UNL_imp_brca1)$unl

# 合并为 tidy 格式
df <- data.frame(
  unl = c(full, becn1, brca1),
  group = factor(rep(
    c("BECN1 & BRCA1",
      "BECN1",
      "BRCA1"),
    times = c(length(full), length(becn1), length(brca1))
  ),level=c("BECN1 & BRCA1",
            "BECN1",
            "BRCA1"))
)
mycols <- setNames(cols[1:3], c("BECN1 & BRCA1",
                                "BECN1",
                                "BRCA1"))
# 绘图
ggplot(df, aes(x = unl, fill = group)) +
  geom_histogram(aes(y = ..density..),        # freq = FALSE -> density
                 position = "identity",       # 叠加而非堆叠
                 alpha = 0.45,                # 透明度，便于比较
                 bins = 30,                   # 可调整或改为 binwidth = ...
                 color = "black") +          # 柱子边框
  scale_fill_manual(values = mycols,
                    name = NULL) +
  coord_cartesian(xlim = c(1, 3), ylim = c(0, 8)) +   # 保留你原来的 x/y 范围
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

save.image("D:/PhD_study2_desktop/new_app/metabric/data/breastcancer_1p.RData")

