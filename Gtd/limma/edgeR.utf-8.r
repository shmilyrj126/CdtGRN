#Bioconductor 安装 edgeR
source('http://bioconductor.org/biocLite.R')
biocLite('edgeR')

#如果连接失败可尝试 source('https://bioconductor.org/biocLite.R')

#加载 edgeR 的同时会自动加载 limma
library(edgeR)

#读取数据，设置分组
targets <- as.matrix(read.delim('gene.txt', sep = '\t', row.names = 1))
group <- rep(c('C1', 'C2'), each = 5)

#（1）构建 DGEList 对象
dgelist <- DGEList(counts = targets, group = group)

#（2）过滤 low count 数据

#例如，直接根据选定的 count 值过滤
keep <- rowSums(dgelist$counts) >= 50
dgelist <- dgelist[keep, ,keep.lib.sizes = FALSE]

#再例如，CPM 标准化（推荐）
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, ,keep.lib.sizes = FALSE]

#（3）标准化，TMM 标准化
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')

plotMDS(dgelist_norm, col = rep(c('red', 'blue'), each = 5), dim = c(1, 2))	#样本无监督聚类

#（4）估算离散值
design <- model.matrix(~group)	#构建分组矩阵
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)	#估算离散值

plotBCV(dge)	#作图查看

#（5）差异分析
#negative binomial generalized log-linear model 拟合
fit <- glmFit(dge, design, robust = TRUE)	#拟合模型
lrt <- glmLRT(fit)	#统计检验

topTags(lrt)
write.csv(topTags(lrt, n = nrow(dgelist$counts)), 'glmLRT.csv', quote = FALSE)	#输出主要结果

dge_de <- decideTestsDGE(lrt, adjust.method = 'fdr', p.value = 0.05)	#查看默认方法获得的差异基因
summary(dge_de)

plotMD(lrt, status = dge_de, values = c(1, -1), col = c('blue', 'red'))	#作图观测
abline(h = c(-1, 1), col = 'gray', lty = 2)

#quasi-likelihood negative binomial generalized log-linear model 拟合
fit <- glmQLFit(dge, design, robust = TRUE)	#拟合模型
lrt <- glmQLFTest(fit)	#统计检验

topTags(lrt)
write.csv(topTags(lrt, n = nrow(dgelist$counts)), 'glmQLFTest.csv', quote = FALSE)	#输出主要结果

dge_de <- decideTestsDGE(lrt, adjust.method = 'fdr', p.value = 0.05)	#查看默认方法获得的差异基因
summary(dge_de)

plotMD(lrt, status = dge_de, values = c(1, -1), col = c('blue', 'red'))	#作图观测
abline(h = c(-1, 1), col = 'gray', lty = 2)

#配对检验
dge_et <- exactTest(dge)	#检验

topTags(dge_et)
write.csv(topTags(dge_et, n = nrow(dgelist$counts)), 'exactTest.csv', quote = FALSE)	#输出主要结果

dge_de <- decideTestsDGE(dge_et, adjust.method = 'fdr', p.value = 0.05)	#查看默认方法获得的差异基因
summary(dge_de)

detags <- rownames(dge)[as.logical(dge_de)]
plotSmear(dge_et, de.tags = detags, cex = 0.5)	#作图观测
abline(h = c(-1, 1), col = 'gray', lty = 2)

#voom 线性建模（limma）
limma_voom <- voom(dgelist, design, plot = TRUE)

fit <- lmFit(limma_voom, design)	#拟合
fit <- eBayes(fit)

topTable(fit, coef = 2)
write.csv(topTable(fit, coef = 2, number = nrow(dgelist$counts)), 'limma_voom.csv', quote = FALSE)	#输出主要结果
