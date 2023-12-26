library(infercnaLight)
library(pheatmap)
#useGenome('hg19') #USE hg39 !!!!
useGenome('hg38')
retrieveGenome()

cna = infercna(m = mgh125, refCells = NULL, n = 5000, noise = 0.1, isLog = TRUE, verbose = FALSE)
cnaM = cna[, !colnames(cna) %in% unlist(refCells)]

plot(cnaM[1,])
#行固定
#転置
pheatmap(t(cnaM), cluster_cols = FALSE)
#match(rownames(cna)[1:100],hg19$symbol)

range(unlogtpm(cnaM))
