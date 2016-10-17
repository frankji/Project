load('~/Project/hongyu/rotation/kang/gene_symbol.rda')
load('~/Project/hongyu/rotation/kang/data_all.rda')
load('~/Project/hongyu/rotation/kang/time_p.rda')
load('~/Project/hongyu/rotation/kang/brain_regions.rda')
load('gene_exp.rda')
TADA <- read.table('~/Project/hongyu/Project/lympho/TADAs.csv', sep = '\t', as.is = TRUE, header = TRUE)
TADA <- TADA[, !apply(TADA, 2, function(x) all(is.na(x)))]
tada <- TADA$tadaFdrAscSscExome
names(tada) <- TADA$RefSeqGeneName

library('WGCNA')
selected_gene <- names(which(!is.na(tada[gene_symbol])))
tada <- tada[selected_gene]
test <- unlist(lapply(gene_exp, function(x) as.vector(x)))
test <- as.data.frame(matrix(test, ncol=208, byrow = TRUE))
test <- collapseRows(datET = test, rowGroup = gene_symbol, rowID = rownames(test), connectivityBasedCollapsing = TRUE)
test <- test$datETcollapsed
test <- test[selected_gene, ]
test <- as.data.frame(test)
test[, 'tada'] <- tada
library(ranger)
RF <- ranger(tada ~. , data = test, write.forest = TRUE, importance = 'impurity')

