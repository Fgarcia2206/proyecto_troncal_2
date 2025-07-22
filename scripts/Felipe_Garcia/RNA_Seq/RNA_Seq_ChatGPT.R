if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c("org.Hs.eg.db", "biomaRt", "goseq", "AnnotationDbi"))
install.packages("tidyverse")

if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) {
  BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
}
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot")






# Librerías
library(tidyverse)
library(org.Hs.eg.db)
library(biomaRt)
library(goseq)
library(AnnotationDbi)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(dplyr)
library(clusterProfiler)
library(enrichplot)

# Ruta de resultados
output_dir <- "C:/Users/fgarc/OneDrive/Escritorio/Doctorado/Ramos/1° Semestre/Troncal/proyecto_troncal2/results/RNA_Seq/RNA_Seq_3"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Leer archivo
rna_data <- read.csv("C:/Users/fgarc/OneDrive/Escritorio/Doctorado/Ramos/1° Semestre/Troncal/proyecto_troncal2/data/RNA-Seq-expression-Norilsk2019.csv")

# Eliminar duplicados
rna_data <- rna_data[!duplicated(rna_data$Gene), ]

# Verificar estructura
head(rna_data)



##########VOLCANO PLOT##########

# Crear variable de estado de significancia
rna_data <- rna_data %>%
  mutate(Significant = ifelse(log_2.fold.change >= 1 & Adjusted.p.value <= 0.05, "Sobreexpresado",
                              ifelse(log_2.fold.change <= -1 & Adjusted.p.value <= 0.05, "Subexpresado", "No significativo")))


# Volcano plot con fondo blanco
volcano_plot <- ggplot(rna_data, aes(x = log_2.fold.change, y = -log10(Adjusted.p.value), color = Significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Sobreexpresado" = "red", "Subexpresado" = "blue", "No significativo" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_bw() +  # Fondo blanco clásico
  theme(
    panel.grid.major = element_blank(),   # Quitar líneas mayores
    panel.grid.minor = element_blank(),   # Quitar líneas menores
    panel.border = element_rect(color = "black"),  # Bordes negros
    panel.background = element_rect(fill = "white"), # Fondo blanco
    plot.background = element_rect(fill = "white")   # Fondo blanco del plot completo
  ) +
  labs(title = "Volcano Plot - Fiebre de Norilsk", x = "Log2 Fold Change", y = "-Log10(P-valor ajustado)")

# Guardar el gráfico
ggsave(filename = file.path(output_dir, "volcano_plot_norilsk.png"), plot = volcano_plot, width = 10, height = 6)



###############################################################################
###############################################################################























####Paso 3: Filtrado de genes sobreexpresados###

# Filtrar genes sobreexpresados
up_genes <- rna_data %>%
  filter(log_2.fold.change >= 1, Adjusted.p.value <= 0.05) %>%
  pull(Gene)

# Guardar lista
write.table(up_genes, file.path(output_dir, "upregulated_genes.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)









#### Paso 4: Anotación GO y KEGG con clusterProfiler ####

# Convertir los IDs ENSEMBL a ENTREZID usando bitr para mejor manejo
ensembl2entrez <- bitr(up_genes,
                       fromType = "ENSEMBL",
                       toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)

# Filtrar IDs válidos
up_entrez_ids <- unique(na.omit(ensembl2entrez$ENTREZID))
cat("Genes sobreexpresados mapeados a ENTREZ:", length(up_entrez_ids), "\n")

# Enriquecimiento GO (BP, MF, CC)
ego_all <- enrichGO(gene          = up_entrez_ids,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = "ENTREZID",
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.2,
                    readable      = TRUE)

# Simplificar GO para reducir redundancias
ego_simplified <- simplify(ego_all, cutoff = 0.7, by = "p.adjust", select_fun = min)

# Guardar resultados GO
write.csv(as.data.frame(ego_simplified),
          file = file.path(output_dir, "GO_enrichment_upregulated_simplified.csv"),
          row.names = FALSE)

# Visualizar enriquecimiento GO solo si hay resultados
if (nrow(as.data.frame(ego_simplified)) > 0) {
  go_plot <- dotplot(ego_simplified, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ .)
  ggsave(file.path(output_dir, "GO_dotplot_upregulated.png"), go_plot, width = 10, height = 8)
} else {
  message("No se encontraron términos GO enriquecidos significativos.")
}

# Enriquecimiento KEGG
up_entrez_chr <- as.character(up_entrez_ids)

ekegg <- enrichKEGG(gene         = up_entrez_chr,
                    organism     = "hsa",
                    pvalueCutoff = 0.1)

# Convertir IDs KEGG a nombres legibles
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Guardar resultados KEGG
write.csv(as.data.frame(ekegg),
          file = file.path(output_dir, "KEGG_enrichment_upregulated.csv"),
          row.names = FALSE)

# Visualizar enriquecimiento KEGG solo si hay resultados
if (nrow(as.data.frame(ekegg)) > 0) {
  kegg_plot <- dotplot(ekegg)
  ggsave(file.path(output_dir, "KEGG_dotplot_upregulated.png"), kegg_plot, width = 10, height = 6)
} else {
  message("No se encontraron términos KEGG enriquecidos significativos.")
}



















































