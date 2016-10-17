library('WGCNA')
source("http://bioconductor.org/biocLite.R")
library(org.Hs.eg.db)
library(annotate)
getid <- function(x) mapIds(org.Hs.eg.db, keys = x, column = "UNIGENE", keytype = 'SYMBOL')
last <- function(x){x[[length(x)]]}
#################### Read DAWN genes information ##########################
dawn <- read.csv(file = 'DAWN.csv', as.is = TRUE)
rasd <- dawn$Gene[dawn$rASD == 'yes']
rasd.names <- getid(rasd)
###########################################################################

#################### Read Case-control expression data ####################
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
#rownames(brain) <- as.character(genes[ID])
colnames(brain) <- sample.name
samples <- strsplit(info[35], split = ',')[[1]][2:80]

cases <- grep('cases', samples)
controls <- grep('controls', samples)

genes <- as.character(genes[ID])
datET <- collapseRows(datET = brain, rowGroup = genes, rowID = rownames(brain), connectivityBasedCollapsing = TRUE)
brain.cases <- datET$datETcollapsed[, cases]
brain.controls <- datET$datETcollapsed[, controls]
#controls.names <- getid(rownames(brain.controls))
#rasd.controls <- as.logical(sapply(rasd.names, function(x) x %in% controls.names))
#rasd.controls <- rasd.controls & !is.na(rasd.names)
rasd.controls <- as.logical(sapply(rasd, function(x) x %in% rownames(brain.controls)))
rasd.controls <- rasd[rasd.controls]
brain.controls.rsd <- brain.controls[rasd.controls, ]
brsd.cor <- cor(t(brain.controls.rsd))
brsd.p <- 2*pnorm(-abs(atanh(brsd.cor) * sqrt(ncol(brain.controls) - 3)))
brsd.q <- p.adjust(brsd.p, method = 'BH')
brsd.n.mat <- matrix(as.numeric(brsd.q <= 0.01), nrow = length(rasd.controls))
brsd.n <- as.numeric(brsd.n.mat)
tmp1 <- rep(1:length(rasd.controls), each = length(rasd.controls))
tmp2 <- rep(1:length(rasd.controls), length(rasd.controls))
brsd.out <- as.data.frame(cbind(tmp2, tmp1, brsd.n), stringsAsFactors = FALSE)
colnames(brsd.out) <- c('A', 'B', 'P')
brsd.out <- brsd.out[brsd.out$A > brsd.out$B & brsd.out$P > 0, ]
write.table(brsd.out, file = 'brsd.csv', sep = ',', row.names = FALSE)
###########################################################################
###################Read Kang's data########################################
load('~/Project/hongyu/Project/DL/gene_exp.rda')
load('~/Project/hongyu/rotation/kang/gene_symbol.rda')
kang <- t(sapply(gene_exp, function(x) as.numeric(x[, c('DFC', 'STC', 'CBC')])))
datET <- collapseRows(datET = kang[!is.na(gene_symbol), ], 
                      rowGroup = gene_symbol[!is.na(gene_symbol)], 
                      rowID = 1:sum(!is.na(gene_symbol)), 
                      connectivityBasedCollapsing = TRUE)
kang.controls <- datET$datETcollapsed
kang.controls.rsd <- kang.controls[rasd.controls, ]
krsd.cor <- cor(t(kang.controls.rsd))
krsd.p <- 2*pnorm(-abs(atanh(krsd.cor) * sqrt(ncol(kang.controls) - 3)))
krsd.q <- p.adjust(krsd.p, method = 'BH')
krsd.n.mat <- matrix(as.numeric(krsd.q <= 0.01), nrow = length(rasd.controls))
krsd.n <- as.numeric(krsd.n.mat)
tmp1 <- rep(1:length(rasd.controls), each = length(rasd.controls))
tmp2 <- rep(1:length(rasd.controls), length(rasd.controls))
krsd.out <- as.data.frame(cbind(tmp2, tmp1, krsd.n), stringsAsFactors = FALSE)
colnames(krsd.out) <- c('A', 'B', 'P')
krsd.out <- krsd.out[krsd.out$A > krsd.out$B & krsd.out$P > 0, ]
write.table(krsd.out, file = 'krsd.csv', sep = ',', row.names = FALSE)

#### Read original kang's data ##############################
load('~/Project/hongyu/20160627/kang/data_all.rda')
load('~/Project/hongyu/20160627/kang/brain_regions.rda')
load('~/Project/hongyu/20160627/kang/time_p.rda')
load('~/Project/hongyu/20160627/kang/gene_symbol.rda')
names(data_all) <- brain_regions
names(time_p) <- brain_regions
gene_exp <- data_all[c('DFC', 'STC', 'CBC')]
times <- time_p[c('DFC', 'STC', 'CBC')]
gene_exp <- mapply(function(x, y) x[, y >= 10] , gene_exp, times)
kang <- do.call(cbind, gene_exp)
BA <- gsub('.+_([A-Z]+)_.+', '\\1', colnames(kang))
datET <- collapseRows(datET = kang[!is.na(gene_symbol), ], 
                      rowGroup = gene_symbol[!is.na(gene_symbol)], 
                      rowID = rownames(kang)[!is.na(gene_symbol)], 
                      connectivityBasedCollapsing = TRUE)
kang <- datET$datETcollapsed
shape <- unlist(lapply(times, function(x) x[x >= 10]))
kang.pca <- prcomp(t(kang), scale. = TRUE)
library('ggplot2')
PCA.plot <- ggplot(data = as.data.frame(kang.pca$x))
PCA.plot <- PCA.plot + geom_point(aes(x = PC1, y = PC2, colour = BA), size = 0.0, alpha = 0)
PCA.plot <- PCA.plot + geom_text(aes(x = PC1, y = PC2, colour = BA, label = shape), size = 5, show.legend = FALSE)
PCA.plot <- PCA.plot + scale_color_manual('Location', values = rainbow(3))
PCA.plot <- PCA.plot + guides(color = guide_legend(override.aes = list(shape = 16, size = 3, alpha = 1)))
PCA.plot <- PCA.plot + theme_classic() + theme(
  text = element_text(size = 25),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))


