# Reading target file
targets <- readTargets("Targets.txt")
# Read the data into R
RG <- read.maimages(targets$FileName, source="agilent")

###Task
#MA plot for a specific sample (ex sample2)
# Density plot
plotMA3by2(RG)

MA <- normalizeWithinArrays(RG, method="loess", bc.method = "normexp", offset = 50)
plotMA3by2(MA)
#### Task
# Density plot after normalizzation

MA.q <- normalizeBetweenArrays(MA, method = "quantile")
# Density plot after quantile normailsation

design <- modelMatrix(targets, ref = "normal")
fit <- lmFit(MA.q, design)
eb <- eBayes(fit)
tb <- topTable (eb, adjust.method="BH")
tt <- topTable(eb, n=nrow(MA.q)) # or tt <- topTable(eb, n=30)

topGenes <- tt[tt[, "P.Value"]<0.001,]
select <- rownames(topGenes)

select <- as.numeric(select)
selectedProbes <- MA.q[select,]


# To extract the original values of the probes from the array

select <-rownames(tt)[tt$P.Value<=0.05 & (tt$logFC>=0.5 | tt$logFC<=-0.5)]

# this will extract all the probes having p value less than 0.05 and log ratio values of >= 0.5 and <= 0.5

# convert them from character vectors to numeric vectors

select <- as.numeric(select)

# 
selectedProbes <- MA.q[select,]

# we extract the original values of all significant genes and try to cluster them for representation.

write.table(selectedProbes, "outputSelected.txt", row.names=T, col.names=T, sep = "\t") # or write.table(tt, "outputSelected.txt", row.names=T, col.names=T, sep = "\t")