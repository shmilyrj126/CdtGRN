library(openxlsx)
cor_gene = read.csv("C:\\Users\\hsc\\Desktop\\yixieyuanma\\p\\p3.csv",header = T,sep = ",")
genes_ind = read.xlsx("C:\\Users\\hsc\\Desktop\\yixieyuanma\\p\\40a1.xlsx",sheet = 1,colNames = T)
rownames(cor_gene) = cor_gene[,1]
cor_gene = cor_gene[,-1]

highcor_ind = which(cor_gene>=0.8,arr.ind = T)
highcor_ind = data.frame(highcor_ind)
n = nrow(highcor_ind)
genesname = rownames(cor_gene) 

newind = data.frame()
for (i in (1:n)){
  newind[i,1] = genesname[highcor_ind[i,1]]
  newind[i,2] = genesname[highcor_ind[i,2]]
}

newgenes = data.frame()
for (i in (1:n)){
  newgenes[i,1] = genes_ind$Gene.Symbol[genes_ind$Probe.Set.ID ==newind[i,1]]
  newgenes[i,2] = genes_ind$Gene.Symbol[genes_ind$Probe.Set.ID ==newind[i,2]]
}

write.csv(newgenes,"C:\\Users\\hsc\\Desktop\\a3.csv",row.names = F)
