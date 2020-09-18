#install.packages("openxlsx")
library(openxlsx)
data1 = read.xlsx("D:\\R\\Rfile\\_paper\\liver.xlsx",sheet = 1,colNames = F)
data1 = as.data.frame(data1)

time = c(4:13,16,18,19)

APOL1 =as.numeric(data1[2,2:14])
ATM = as.numeric(data1[3,2:14])
BAX = as.numeric(data1[4,2:14])
BCL2 = as.numeric(data1[5,2:14])
Casein = as.numeric(data1[6,2:14])
CDKN2A = as.numeric(data1[7,2:14])
CHEK1 = as.numeric(data1[8,2:14])
CHEK2 = as.numeric(data1[9,2:14])
E2F1 = as.numeric(data1[10,2:14])
GBP = as.numeric(data1[11,2:14])
MAM2 = as.numeric(data1[12,2:14])
MYC = as.numeric(data1[13,2:14])
ras = as.numeric(data1[14,2:14])
SMAD4 = as.numeric(data1[15,2:14])
TP53 = as.numeric(data1[16,2:14])

##range(APOL1)

plot(x = time, y = APOL1,type = "o",col = "black",xaxt="n",yaxt="n",xlim = c(4,21),ylim=c(0,150),xlab = "time/weeks",ylab = "Gene Expression")
axis(1,c(4:19))
axis(2,seq(0,150,10))

lines(x = time, y = ATM,type = "o",col = "red")
lines(x = time, y = BAX,type = "o",col = "green")
lines(x = time, y = BCL2,type = "o",col = "blue")
lines(x = time, y = Casein,type = "o",col = "yellow")
lines(x = time, y = CDKN2A,type = "o",col = "springgreen2")
lines(x = time, y = CHEK1,type = "o",col = "skyblue")
lines(x = time, y = CHEK2,type = "o",col = "orange")
lines(x = time, y = E2F1,type = "o",col = "salmon2")
lines(x = time, y = GBP,type = "o",col = "purple")
lines(x = time, y = MAM2,type = "o",col = "pink")
lines(x = time, y = MYC,type = "o",col = "lightblue")
lines(x = time, y = ras,type = "o",col = "brown")
lines(x = time, y = SMAD4,type = "o",col = "deeppink")
lines(x = time, y = TP53,type = "o",col = "tomato3")

legend("topright",data1[2:16,1],
       col=c("black","red","green","blue","yellow","springgreen2","skyblue",
             "orange","salmon2","purple", "pink", "lightblue","brown","deeppink","tomato3"), 
       cex = 0.5,text.width =0.5,
       text.col=c("black","red","green","blue","yellow","springgreen2","skyblue",
                  "orange","salmon2","purple", "pink", "lightblue","brown","deeppink","tomato3"),
       lty = c(1,1,1),bty=n,x.intersp = 0.5, y.intersp = 0.5)

##corr
data = as.matrix(data1[2:16,2:14])
data = t(data)
p = ncol(data)
cormatrix = matrix(0,p,p)
for (i in 1:(p-1)){
  for(j in (i+1):p){
    cormatrix[i,j] = cor(as.numeric(data[,i]),as.numeric(data[,j]))
  }
}
cormatrix = cormatrix+t(cormatrix)
diag(cormatrix) = rep(1,15)

cormatrix = as.data.frame(cormatrix)
colnames(cormatrix) = data1[2:16,1]
rownames(cormatrix) = data1[2:16,1]
write.csv(cormatrix,"D:\\R\\Rfile\\_paper\\cormatrix.csv",row.names = TRUE)
