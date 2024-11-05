library(clusterProfiler)
library(org.Hs.eg.db)
##################################P/N
P_uniprot <- mapIds(org.Hs.eg.db,
                    keys = dot_P[dot_P$group=="enzyme",]$id,
                    column = "SYMBOL",
                    keytype = "UNIPROT",
                    multiVals = "first")

length(combined_P)

P_uniprot <- na.omit(unlist(P_uniprot, use.names = FALSE))

combined_P <- union(dot_P[dot_P$group=="gene",]$id, P_uniprot)

combined_P <- unique(combined_P)

diffgene<-bitr(combined_P,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")

kk <- clusterProfiler::enrichKEGG(gene = diffgene$ENTREZID,
                                  keyType = "kegg",
                                  organism= "human", 
                                  qvalueCutoff = 0.05 , 
                                  pvalueCutoff=0.05,
                                  pAdjustMethod = "BH")

ek.rt <- data.frame(kk@result)
ek.rt <- separate(data=ek.rt, col=GeneRatio, into = c("GR1", "GR2"), sep = "/") 
ek.rt <- separate(data=ek.rt, col=BgRatio, into = c("BR1", "BR2"), sep = "/") 
ek.rt <- mutate(ek.rt, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) 
ek.rt <- arrange(ek.rt,desc(enrichment_factor))

ek.rt10 <- ek.rt %>% filter(row_number() >= 1,row_number() <= 20)
ek.rt10$Description <- factor(ek.rt10$Description, levels =rev(ek.rt10$Description) )

plot_P <- ggplot(ek.rt10, aes(y=Description, x=enrichment_factor, fill=p.adjust)) +
  geom_bar(stat = "identity", width=0.7) +
  scale_fill_gradient(low = "#2E6E9E", high = "#B62024") +
  labs(title = "KEGG Pathways Enrichment", x = "enrichment factor", y = "Pathways") +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14), 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )


##################################S/N
S_uniprot <- mapIds(org.Hs.eg.db,
                    keys = dot_S[dot_S$group=="enzyme",]$id,
                    column = "SYMBOL",
                    keytype = "UNIPROT",
                    multiVals = "first")

length(combined_P)

S_uniprot <- na.omit(unlist(S_uniprot, use.names = FALSE))

combined_S <- union(dot_S[dot_S$group=="gene",]$id, S_uniprot)

combined_S <- unique(combined_S)

diffgene_S<-bitr(combined_S,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")

kk_S <- clusterProfiler::enrichKEGG(gene = diffgene_S$ENTREZID,
                                    keyType = "kegg",
                                    organism= "human", 
                                    qvalueCutoff = 0.05 , 
                                    pvalueCutoff=0.05,
                                    pAdjustMethod = "BH")

ek.rt_S <- data.frame(kk_S@result)
ek.rt_S <- separate(data=ek.rt_S, col=GeneRatio, into = c("GR1", "GR2"), sep = "/") 
ek.rt_S <- separate(data=ek.rt_S, col=BgRatio, into = c("BR1", "BR2"), sep = "/") 
ek.rt_S <- mutate(ek.rt_S, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) 
ek.rt_S <- arrange(ek.rt_S,desc(enrichment_factor))

ek.rt10_S <- ek.rt_S %>% filter(row_number() >= 1,row_number() <= 20)
ek.rt10_S$Description <- factor(ek.rt10_S$Description, levels =rev(ek.rt10_S$Description) )

plot_S <- ggplot(ek.rt10_S, aes(y=Description, x=enrichment_factor, fill=p.adjust)) +
  geom_bar(stat = "identity", width=0.7) +
  scale_fill_gradient(low = "#2E6E9E", high = "#B62024") +
  labs(title = "KEGG Pathways Enrichment", x = "enrichment factor", y = "Pathways") +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),  
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )

##################################S/P
SP_uniprot <- mapIds(org.Hs.eg.db,
                     keys = dot_S[dot_S_P$group=="enzyme",]$id,
                     column = "SYMBOL",
                     keytype = "UNIPROT",
                     multiVals = "first")

length(combined_P)

SP_uniprot <- na.omit(unlist(SP_uniprot, use.names = FALSE))

combined_SP <- union(dot_S_P[dot_S_P$group=="gene",]$id, SP_uniprot)

combined_SP <- unique(combined_SP)

diffgene_SP<-bitr(combined_SP,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")

kk_SP <- clusterProfiler::enrichKEGG(gene = diffgene_SP$ENTREZID,
                                     keyType = "kegg",
                                     organism= "human", 
                                     qvalueCutoff = 0.05 , 
                                     pvalueCutoff=0.05,
                                     pAdjustMethod = "BH")

ek.rt_SP <- data.frame(kk_SP@result)
ek.rt_SP <- separate(data=ek.rt_SP, col=GeneRatio, into = c("GR1", "GR2"), sep = "/")
ek.rt_SP <- separate(data=ek.rt_SP, col=BgRatio, into = c("BR1", "BR2"), sep = "/") 
ek.rt_SP <- mutate(ek.rt_SP, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) 
ek.rt_SP <- arrange(ek.rt_SP,desc(enrichment_factor))

ek.rt10_SP <- ek.rt_SP %>% filter(row_number() >= 1,row_number() <= 20) %>% head(20)

ek.rt10_SP$Description <- factor(ek.rt10_SP$Description, levels =rev(ek.rt10_SP$Description) )


plot_SP <- ggplot(ek.rt10_SP, aes(y=ek.rt10_SP$Description, x=ek.rt10_SP$enrichment_factor, fill=p.adjust)) +
  geom_bar(stat = "identity", width=0.7) +
  scale_fill_gradient(low = "#2E6E9E", high = "#B62024") +
  labs(title = "KEGG Pathways Enrichment", x = "enrichment factor", y = "Pathways") +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14), 
    axis.text.y = element_text(size = 14), 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )


