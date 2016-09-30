lympho <- read.csv('~/Project/hongyu/GSE37772_series_matrix.csv', skip = 71, as.is = TRUE)
ID <- lympho[, 1]
lympho <- lympho[-1, 1:440]
rownames(llympho) <- ID
