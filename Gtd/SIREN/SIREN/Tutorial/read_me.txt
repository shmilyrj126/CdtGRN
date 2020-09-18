# follow the rules :-)

source("SIREN.R")

exp=as.matrix(read.table("Expression_Format.txt",sep="\t"))
w=as.matrix(read.table("Weighting_Matrix.txt",sep="\t"))
net=as.matrix(read.table("Network_Format.txt",sep="\t"))

Result=SIREN(exp,w,net)
write.table(Result,file='Result.txt',sep="\t")