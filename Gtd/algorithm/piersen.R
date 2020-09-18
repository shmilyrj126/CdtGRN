
library(openxlsx)

data2 = read.xlsx("C:\\Users\\hsc\\Desktop\\x40.xlsx",sheet = 1,colNames = T)
str(data2)

##corr

data = as.matrix(data2[1:11,2:22])

data = t(data)

p = ncol(data)

cormatrix = matrix(0,p,p)

for (i in 1:(p-1)){
  for(j in (i+1):p){
    cormatrix[i,j] = cor(as.numeric(data[,i]),as.numeric(data[,j]),method = "pearson")
  }
}

#cormatrix = cormatrix+t(cormatrix)

#diag(cormatrix) = rep(1,8867)

cormatrix = as.data.frame(cormatrix)
colnames(cormatrix) = data2[1:11,1]
rownames(cormatrix) = data2[1:11,1]
write.csv(cormatrix,"C:\\R\\p40.csv",row.names = TRUE)

#count
m=0
for (i in 1:(p-1)){
  for(j in (i+1):p){
    if(cormatrix[i,j]<(-0.8)|cormatrix[i,j]>(0.8))
      m=m+1
  }
}

m
  



