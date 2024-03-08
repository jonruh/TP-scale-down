library("DESeq2")
library(dplyr)
library(EnhancedVolcano)

metadata<-read.csv(file = "../../data/RNAseq_data/metadata.tsv", sep= "\t", row.names = 1,check.names = FALSE)
countdata<-read.csv(file = "../../data/RNAseq_data/FeatureCounts.tsv", sep= "\t", row.names = 1,check.names = FALSE)

dds <- DESeqDataSetFromMatrix(countData = countdata,
                             colData = metadata,
                             design= ~ condition)

dds <- DESeq(dds)

resultsNames(dds) # lists the coefficients

res <- results(dds, contrast=c("condition","DDB35_scale-down","DDB35_control"))
write.csv(res,"../../data/RNAseq_data/DESeq_res_DDB35_sd_vs_ctrl.csv")

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.0
)

res <- results(dds, contrast=c("condition","HMP3071_scale-down","HMP3071_control"))
write.csv(res,"../../data/RNAseq_data/DESeq_res_HMP3071_sd_vs_ctrl.csv")

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.0
)

### Strain

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = metadata,
                              ~ process_condition + strain_ID)

dds <- DESeq(dds)

# res <- results(dds, contrast=c("strain_ID","HMP3071","DDB35"))
res <- results(dds, name="strain_ID_HMP3071_vs_DDB35")
write.csv(res,"../../data/RNAseq_data/DESeq_res_HMP3071_vs_DDB35.csv")

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.0
)