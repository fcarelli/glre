# qPCR enrichment plot

qpcr_enrichment = matrix(c(0.157490131236859, 0.143587294374629, 0.154963462492373, 0.148307829363560), ncol=2)

pdf(file = "plots/qPCR_T05F1.2_fold_change.pdf", width = 4, height = 6, useDingbats = FALSE)

par(mfrow=c(1,1), mar=c(4,4,3,2))
k = barplot(colMeans(qpcr_enrichment), ylim=c(0, 1), las = 1, ylab = "fold-change mut vs control",
        names.arg=c("primer set 1", "primer set 2"), main = "qPCR T05F1.2")
stripchart(list(qpcr_enrichment[,1], qpcr_enrichment[,2]), method = "jitter", vertical = T, add = T, at = k[,1], pch = 19)

dev.off()


HIM17_qpcr_enrichment = read.table("HIM17_qPCR/HIM17_pct_input.txt")

pdf(file = "plots/qPCR_HIM17_pct_input.pdf", width = 7, height = 6, useDingbats = FALSE)

par(mfrow=c(1,1), mar=c(4,4,3,2))
k = barplot(t(matrix(c(rowMeans(HIM17_qpcr_enrichment[,c(1,2)]),
                       rowMeans(HIM17_qpcr_enrichment[,c(3,4)]),
                       rowMeans(HIM17_qpcr_enrichment[,c(5,6)])), ncol = 3, byrow = F)), ylim = c(0, 1.5),  las = 1, ylab = "% input", beside=T, 
            names.arg=c("+ ctrl", "- ctrl 1", "- ctrl 2", "endogenous", "construct"), main = "ChIP-qPCR HIM-17")
stripchart(list(HIM17_qpcr_enrichment[1, c(1,2)],
                HIM17_qpcr_enrichment[2, c(1,2)],
                HIM17_qpcr_enrichment[3, c(1,2)],
                HIM17_qpcr_enrichment[4, c(1,2)],
                HIM17_qpcr_enrichment[5, c(1,2)]), method = "jitter", vertical = T, add = T, at = k[1,], pch = 19)
stripchart(list(HIM17_qpcr_enrichment[1, c(3,4)],
                HIM17_qpcr_enrichment[2, c(3,4)],
                HIM17_qpcr_enrichment[3, c(3,4)],
                HIM17_qpcr_enrichment[4, c(3,4)],
                HIM17_qpcr_enrichment[5, c(3,4)]), method = "jitter", vertical = T, add = T, at = k[2,], pch = 19)
stripchart(list(HIM17_qpcr_enrichment[1, c(5,6)],
                HIM17_qpcr_enrichment[2, c(5,6)],
                HIM17_qpcr_enrichment[3, c(5,6)],
                HIM17_qpcr_enrichment[4, c(5,6)],
                HIM17_qpcr_enrichment[5, c(5,6)]), method = "jitter", vertical = T, add = T, at = k[3,], pch = 19)


dev.off()
