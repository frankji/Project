# Read Data
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
sample.loc <- gsub(".+: (\\w+)", "\\1", sample.loc)
colnames(brain) <- sample.name
samples <- strsplit(info[35], split = ',')[[1]][2:80]
samples <- gsub(".+: (\\w+)", "\\1", samples)

# PCA of expression
library('ggplot2')
brain.pca <- prcomp(t(brain), scale. = TRUE)
PCA.plot <- ggplot(data = as.data.frame(brain.pca$x))
PCA.plot <- PCA.plot + geom_point(aes(x = PC1, y = PC2, colour = sample.loc, shape = samples), size = 3)
PCA.plot <- PCA.plot + scale_color_manual('Location', values = rainbow(3))
PCA.plot <- PCA.plot + theme_classic() + theme(
  text = element_text(size = 25),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

cerebellum.shape <- samples[sample.loc == 'Cerebellum']
cerebellum.pca <- prcomp(t(brain[, sample.loc == 'Cerebellum']), scale. = TRUE)
PCA.plot <- ggplot(data = as.data.frame(cerebellum.pca$x))
PCA.plot <- PCA.plot + geom_point(aes(x = PC1, y = PC2, color = cerebellum.shape), size = 3)
PCA.plot <- PCA.plot + scale_color_manual('Location', values = rainbow(2))
PCA.plot <- PCA.plot + theme_classic() + theme(
  text = element_text(size = 25),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

cortex.groups <- samples[grep('cortex', sample.loc)]
cortex.color <- sample.loc[grep('cortex', sample.loc)]
cortex.pca <- prcomp(t(brain[, grep('cortex', sample.loc)]), scale. = TRUE)
PCA.plot <- ggplot(data = as.data.frame(cortex.pca$x))
PCA.plot <- PCA.plot + geom_point(aes(x = PC1, y = PC2, colour = cortex.color, shape = cortex.groups), size = 3)
PCA.plot <- PCA.plot + scale_color_manual('Location', values = rainbow(2))
PCA.plot <- PCA.plot + theme_classic() + theme(
  text = element_text(size = 25),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

### Read Kang's Data####
load('~/Project/hongyu/Project/DL/gene_exp.rda')
load('~/Project/hongyu/rotation/kang/gene_symbol.rda')
kang <- t(sapply(gene_exp, function(x) as.numeric(x[8:13, c('DFC', 'STC', 'CBC')])))
datET <- collapseRows(datET = kang[!is.na(gene_symbol), ], 
                      rowGroup = gene_symbol[!is.na(gene_symbol)], 
                      rowID = 1:sum(!is.na(gene_symbol)), 
                      connectivityBasedCollapsing = TRUE)
kang <- datET$datETcollapsed
colnames(kang) <- rep(c('DFC', 'STC', 'CBC'), each = 6)
shape <- as.character(rep(10:15, 3))
kang.pca <- prcomp(t(kang), scale. = TRUE)
PCA.plot <- ggplot(data = as.data.frame(kang.pca$x))
PCA.plot <- PCA.plot + geom_point(aes(x = PC1, y = PC2, colour = colnames(kang)), size = 0.0, alpha = 0)
PCA.plot <- PCA.plot + geom_text(aes(x = PC1, y = PC2, colour = colnames(kang), label = shape), size = 5, show.legend = FALSE)
PCA.plot <- PCA.plot + scale_color_manual('Location', values = rainbow(3))
PCA.plot <- PCA.plot + guides(color = guide_legend(override.aes = list(shape = 16, size = 3, alpha = 1)))
PCA.plot <- PCA.plot + theme_classic() + theme(
  text = element_text(size = 25),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))


