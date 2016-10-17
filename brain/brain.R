genes <- scan('~/Project/hongyu/GPL6883-11606.csv', skip = 27, sep = '\n', what = character())
labs <- as.character(sapply(genes, function(x) strsplit(x, split = '\t')[[1]][1]))
genes <- as.character(sapply(genes, function(x) strsplit(x, split = '\t')[[1]][6]))
names(genes) <- labs
genes <- genes[-1]
brain <- read.table('GSE28521_series_matrix.csv', skip = 59, as.is = TRUE, sep = ',')
brain <- brain[, !apply(brain, 2, function(x) all(is.na(x)))]
brain <- brain[!apply(brain, 1, function(x) any(is.na(x))), ]
ID <- brain[, 1]
brain <- brain[, -1]

info <- scan('GSE28521_series_matrix.csv', nlines = 59, what = character(), sep = '\n')

sample.name <- strsplit(info[27], split = ',')[[1]][2:80]
sample.loc <- strsplit(info[36], split = ',')[[1]][2:80]
Frontal <- grep('Frontal', sample.loc)
rownames(brain) <- as.character(genes[ID])
colnames(brain) <- sample.name
samples <- strsplit(info[35], split = ',')[[1]][2:80]

cases <- grep('cases', samples)
controls <- grep('controls', samples)

TADA <- read.table('../lympho/TADAs.csv', sep = '\t', as.is = TRUE, header = TRUE)
TADA <- TADA[, !apply(TADA, 2, function(x) all(is.na(x)))]
tada <- TADA$tadaFdrAscSscExome
names(tada) <- TADA$RefSeqGeneName


library('WGCNA')
genes <- as.character(genes[ID])
datET <- collapseRows(datET = brain, rowGroup = genes, rowID = rownames(brain), connectivityBasedCollapsing = TRUE)
brain.cases <- datET$datETcollapsed[, cases]
brain.controls <- datET$datETcollapsed[, controls]

selected.genes <- rownames(datET$datETcollapsed)
selected.genes <- selected.genes[!is.na(tada[selected.genes])]
tada <- tada[selected.genes]

brain.cases <- brain.cases[selected.genes, ]
brain.controls <- brain.controls[selected.genes, ]

cases.sim <- cor(t(brain.cases))
diag(cases.sim) <- 0
cases.sim <- atanh(cases.sim)
controls.sim <- cor(t(brain.controls))
diag(controls.sim) <- 0
controls.sim <- atanh(controls.sim)
controls.sim1 <- controls.sim * sqrt(37)
controls.sim2 <- 2*pnorm(-abs(controls.sim1))
controls.sim2 <- p.adjust(controls.sim2, method = 'BH')
controls.sim2 <- matrix(controls.sim2, ncol = length(selected.genes))
rownames(controls.sim2) <- selected.genes
colnames(controls.sim2) <- selected.genes
z.sim <- abs(cases.sim - controls.sim)/(sqrt(1/(36)+1/(37)))

case.control <- 2*pnorm(-abs(z.sim))
case.control.adjusted<- p.adjust(case.control, method = 'BH')
case.control <- case.control.adjusted < 0.05
case.control <- matrix(case.control, ncol = length(selected.genes))
rownames(case.control) <- selected.genes
colnames(case.control) <- selected.genes
test <- apply(case.control, 1, function(x) sum(x))
q <- quantile(as.numeric(test), 0.90)
impacted <- names(which(test > q))
impacted <- impacted[order(test[impacted], decreasing = TRUE)]

# Get who is enriched in the interaction within the impacted.
enrich <- c()
for(gene in selected.genes){
  x <- sum(controls.sim2[impacted, gene] < 0.01)
  m <- sum(controls.sim2[gene, ] < 0.01)
  k <- length(impacted)
  n <- length(selected.genes) - m
  enrich <- c(enrich, phyper(x, m, n, k, lower.tail = FALSE))
}
enrich <- p.adjust(enrich, method = 'bonferroni')
names(enrich) <- selected.genes

test1 <- apply(controls.sim2, 1, function(x) sum((x<0.05) & (tada<0.05)))
names(test1) <- selected.genes
random.tada <- rep(FALSE, length(selected.genes))
names(random.tada) <- selected.genes
random.tada['DYRK1A'] <- TRUE
test2 <- apply(controls.sim2, 1, function(x) sum((x<0.01) & (random.tada)))

library(ggplot2)
test <- as.data.frame(cbind(test, test1))

vio <- ggplot(data = test)
vio + geom_violin(aes(x = factor(test1), y  = test))
mat.2nd <- matrix(FALSE, nrow = length(selected.genes), ncol = length(selected.genes))

for(i in 1:nrow(mat.2nd)){
  mat.2nd[i, controls.edges[i, ]] <- test[controls.edges[i, ]]
}
test2 <- apply(mat.2nd, 1, sum)
names(test2) <- selected.genes
top30 <- names(sort(tada)[1:50])