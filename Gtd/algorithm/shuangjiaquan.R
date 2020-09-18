library(openxlsx)
##duquwenjian
data3 = read.xlsx("C:\\Users\\hsc\\Desktop\\s2.xlsx",sheet = 1,colNames = T)
str(data3)

##huansunjuzhen
data=as.matrix(data3[1:5,2:8])

data=t(data)

p = ncol(data)

cormatrix = matrix(0,p,p)

#jisuanxishu
fun=function(x,y){
  a=sum(x*y)/(sqrt(sum(x*x))*sqrt(sum(y*y)))
  return(a)
}

#xunhuanbianli
for(i in 1:(p-1)){
  for(j in (i+1):p){
    cormatrix[i,j]=fun(as.numeric(data[,i]),as.numeric(data[,j]))
  }
}
#yuanlai+zhuanzhi
#cormatrix = cormatrix+t(cormatrix)
#buchongduijiao
#diag(cormatrix) = rep(1,58)
#zhunahuaweishujukuang
cormatrix = as.data.frame(cormatrix)
colnames(cormatrix) = data3[1:5,1]
rownames(cormatrix) = data3[1:5,1]
write.csv(cormatrix,"C:\\R\\cormatrix2.csv",row.names = TRUE)