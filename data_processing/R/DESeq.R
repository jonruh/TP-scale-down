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

res <- results(dds, name="condition_HMP3071_scale.down_vs_HMP3071_control")
write.csv(res,"../../data/RNAseq_data/DESeq_res_HMP3071_sd_vs_ctrl.csv")

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-5,
                FCcutoff = 1.0
)

res <- results(dds, contrast=c("condition","SDT714_scale.down","SDT714_control"))
write.csv(res,"../../data/RNAseq_data/DESeq_res_SDT714_sd_vs_ctrl.csv")

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-5,
                FCcutoff = 1.0
)

res <- results(dds, contrast=c("condition","SDT714_control","HMP3071_control"))
write.csv(res,"../../data/RNAseq_data/DESeq_res_SDT714_ctrl_vs_HMP3071_ctrl.csv")

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-5,
                FCcutoff = 1.0
)

res <- results(dds, contrast=c("condition","SDT714_scale.down","HMP3071_scale.down"))
write.csv(res,"../../data/RNAseq_data/DESeq_res_SDT714_sd_vs_HMP3071_sd.csv")

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-5,
                FCcutoff = 1.0
)