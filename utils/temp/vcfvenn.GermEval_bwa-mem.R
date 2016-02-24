require(VennDiagram)
pdf("vcfvenn.GermEval.pdf")
colortable <- c("orange", "red", "blue")
nametable <- c("Freebayes", "GATK (HC)", "Illumina Plantinum")
venn.plot <- draw.triple.venn(
            102782,92595,57944,30508+52517,2298+52517,743+52517,52517,
            c("Freebayes (Parrallel)", "GATK (HC)", "Illumina Plantinum v0.7"),
            fill = c("blue", "red", "green"),
            lty = "blank",
            cex = 2,
            cat.cex = 2,
            cat.col = c("blue", "red", "green"))
grid.draw(venn.plot);

dev.off()

