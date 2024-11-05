library(FactoMineR)
library(factoextra)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(dplyr)
library(tidyr)
library(ropls)
library(mixOmics)
library(rstatix)
########tidyr################################################PCA analysis
################################################
load("sample_A/B/C/D.RData")

t_sample_A/B/C/D <- t(sample_A/B/C/D)
t_sample_A/B/C/D <- data.frame(t_sample_A/B/C/D)

t_sample_A/B/C/D$group <- ifelse(rownames(t_sample_A/B/C/D) %in% paste0("N",1:12),"Nomal",
                            ifelse(rownames(t_sample_A/B/C/D) %in% paste0("P",1:8),"Poorly_differentiated","signet_ring"))  #"Poorly_differentiated","signet_ring",Nomal

t_sample_A/B/C/D.pca=PCA(t_sample_A/B/C/D[,-ncol(t_sample_A/B/C/D)],
                    graph = TRUE)

pca_score1_8 <- t_sample_A/B/C/D.pca$ind$coord %>%
  data.frame()%>%
  dplyr::mutate(group=t_sample_A/B/C/D$group)

P_PCA <- ggplot(pca_score1_8, aes(Dim.1,Dim.2,color = group),color=Sampletype)+
  geom_point(size=4)+
  stat_ellipse(type="t",linetype=2)+
  labs(x=paste("Dim 1(",round(get_eig(t_sample_A/B/C/D.pca)[1,2],2),"%)",sep = ""),
       y=paste("Dim 2(",round(get_eig(t_sample_A/B/C/D.pca)[2,2],2),"%)",sep = ""))+
  labs(title = "PCA")+
  theme_bw()+
  scale_color_manual(values = c('#1B9E77','#D95F02',"#7570B3"))+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)

ggsave("PCA.png", plot = P_PCA, width = 8, height = 6, dpi = 300)

################################################PLS-DA analysis
################################################
df_plsda1_8<- opls(x=t_sample_A/B/C/D[,-ncol(t_sample_A/B/C/D)], 
                   y=t_sample_A/B/C/D$group,
                   orthoI = 0)

sample.score1_8 = df_plsda1_8@scoreMN %>%  
  as.data.frame() %>%
  mutate(group = t_sample_A/B/C/D$group)

ggplot(sample.score1_8, aes(p1, p2, color = group),color=Sampletype)+
  geom_point(size=3)+
  stat_ellipse(type="t",linetype=2)+
  labs(x=paste("X-variate 1(",df_plsda1_8@modelDF[1,"R2X"]*100,"%)",sep = ""),
       y=paste("X-variate 2(",df_plsda1_8@modelDF[2,"R2X"]*100,"%)",sep = ""))+
  labs(title = "PLS-DA")+
  geom_label_repel(aes(label = rownames(pca_score1_8)))+
  theme_bw()+
  scale_color_manual(values = c('#1B9E77','#D95F02',"#7570B3"))+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend(title = NULL),shape=guide_legend(title = NULL))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)

################################################Volcano plot
#########In differentially differentiated regions, pairwise comparisons will be conducted, followed by a univariate analysis, and a volcano plot will be created. 
#########Normal tissue regions are labeled with an "N" prefix, poorly differentiated regions with a "P" prefix, and signet ring cell carcinoma regions with an "S" prefix. 
#########For example, in sample A, we will perform a univariate analysis between the signet ring cell carcinoma region (S) and the poorly differentiated carcinoma region (P) to identify upregulated and downregulated metabolites. 

signet1 <- sample_A[,-c(12:18)]
signet1 <- mutate(signet1,m_z=rownames(signet1),.before = 1)
rownames(signet1) <- c(1:10053)  #row number

group_signet1 <- names(signet1)[-1]%>% data.frame()
group_signet1 <- dplyr::mutate(group_signet1,class=rep(c("Nomal","signet_ring"),c(11,4)))  #choice for signet_ring,Poorly_differentiated，"Nomal
names(group_signet1)[1] <- "sample"

df_signet1 <- signet1 %>%
  as_tibble() %>%
  pivot_longer(-1,names_to = "sample",values_to = "value") %>%  
  left_join(group_signet1,by=c("sample" = "sample"))

FC_signet1 = df_signet1 %>%
  group_by(m_z,class) %>%  
  summarise(mean = mean(value,na.rm=T)) %>%               
  pivot_wider(names_from = class,values_from = mean) %>%  
  summarise(FC = signet_ring/Nomal) #choice for signet_ring,Poorly_differentiated，"Nomal

P_signet1 = df_signet1 %>%
  group_by(m_z) %>%
  t_test(value ~ class,var.equal=T) 

FDR_signet1 = P_signet1 %>%
  select(1,last_col()) %>%
  mutate(FDR = p.adjust(.$p,method = "BH"))  

Data_signet1 <- df_signet1 %>%
  as_tibble() %>%
  left_join(FC_signet1) %>%
  left_join(FDR_signet1)

Data_signet1 <- distinct(Data_signet1,Data_signet1$m_z,.keep_all = TRUE)

FC = 1.5
PValue = 0.05 

Data_signet1$sig = ifelse(Data_signet1$p < PValue & abs(log(Data_signet1$FC)) >= log(FC),
                          ifelse(log(Data_signet1$FC)> log(FC) ,'Up','Down'),
                          'Stable')
table(data$sig)
sum(is.na(data_signet1$sig))

volcano_signet1 <- ggplot(Data_signet1, aes(log(Data_signet1$FC),-1*log10(Data_signet1$p))) +
  geom_hline(yintercept =-log10(0.05), linetype = 'dashed', size = 0.5) +
  geom_vline(xintercept =c(-log(1.5),log(1.5)), linetype = 'dashed', size = 0.5) +
  geom_point(aes(color = sig),
             size = 1, 
             alpha = 0.5) +
  labs(title="volcanoplot",
       x="log(FC)", 
       y="-log[10](PValue)",
       colour="Group") +
  scale_color_manual(values = c("#0709F7","black","#FB0F15")) +
  theme_bw(base_size = 12) +
  theme(legend.title = element_text(),
        legend.text = element_text(color = 'black',size = 12, family = 'Arial', face = 'plain'),
        plot.title = element_text(size=15,hjust = 0.5),
        panel.background = element_blank())

ggsave("sampleA_volcano.png", plot = volcano_signet1, width = 10, height = 8, dpi = 300)



