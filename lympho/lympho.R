genes <- scan('~/Project/hongyu/GPL6883-11606.csv', skip = 27, sep = '\n', what = character())
labs <- as.character(sapply(genes, function(x) strsplit(x, split = '\t')[[1]][1]))
genes <- as.character(sapply(genes, function(x) strsplit(x, split = '\t')[[1]][6]))
names(genes) <- labs
genes <- genes[-1]
lympho <- read.csv('~/Project/hongyu/GSE37772_series_matrix.csv', skip = 71, as.is = TRUE)
lympho <- lympho[-nrow(lympho), ]
ID <- lympho[, 1]
genes <- as.character(genes[ID])

TADA <- read.table('TADAs.csv', sep = '\t', as.is = TRUE, header = TRUE)
TADA <- TADA[, !apply(TADA, 2, function(x) all(is.na(x)))]
tada <- TADA$tadaFdrAscSscExome
names(tada) <- TADA$RefSeqGeneName

lympho <- lympho[, 2:440]
lympho.text <- scan(file = '~/Project/hongyu/GSE37772_series_matrix.csv', nlines = 71, sep = '\n', what = character())
samples <- strsplit(lympho.text[36], split = ',')[[1]] 
samples <- samples[2:440]
cases <- grep('Autism', samples)
controls <- grep('Control', samples)


library('WGCNA')
datET <- collapseRows(datET = lympho, rowGroup = genes, rowID = rownames(lympho), connectivityBasedCollapsing = TRUE)
lympho.cases <- datET$datETcollapsed[, cases]
lympho.controls <- datET$datETcollapsed[, controls]
p.values <- rep(NA, 18630)

for(i in 1:18630){
  p.values[i] <- t.test(lympho.cases[i, ], lympho.controls[i, ])$p.value
}

q <- p.adjust(p.values, method = 'BH')
names(q) <- rownames(datET$datETcollapsed)


selected.genes <- rownames(datET$datETcollapsed)
selected.genes <- selected.genes[!is.na(tada[selected.genes])]
q.values <- q[selected.genes]
tada.values <- tada[selected.genes]
plot(tada.values~q.values)

lympho.cases <- lympho.cases[selected.genes, ]
lympho.controls <- lympho.controls[selected.genes, ]

cases.sim <- adjacency(t(lympho.cases), type = 'unsigned', power = 1)
controls.sim <- adjacency(t(lympho.controls), type = 'unsigned', power = 1)

cases.edges <- cases.sim > 0.75
controls.edges <- controls.sim > 0.75

case.control <- matrix(NA, nrow = length(selected.genes), ncol = length(selected.genes)) 
for (i in 1:length(selected.genes)) {
 case.control[i, ] <- xor(cases.edges[i, ], controls.edges[i, ])
}

rownames(case.control) <- selected.genes
top30 <- names(sort(tada.values)[1:30])
change.30 <- apply(case.control[top30, ], 1, sum)
test <- apply(case.control, 1, sum)