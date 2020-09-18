library(openxlsx)
cor_gene = read.csv("C:\\Users\\hsc\\Desktop\\11_2501\\p1.csv",header = T,sep = ",")
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

write.csv(newind,"C:\\Users\\hsc\\Desktop\\w4.csv",row.names = F)