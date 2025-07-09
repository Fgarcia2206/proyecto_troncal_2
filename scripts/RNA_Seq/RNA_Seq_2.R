# Instalar paquetes si es necesario
if (!requireNamespace("ggplot2")) install.packages("ggplot2")
if (!requireNamespace("biomaRt")) install.packages("biomaRt")
if (!requireNamespace("clusterProfiler")) install.packages("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("pathview")) BiocManager::install("pathview")

# Cargar librerías
library(ggplot2)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)

# Ruta de salida (reemplaza "\" por "/" en R)
output_path <- "C:/Users/fgarc/OneDrive/Escritorio/Doctorado/Ramos/1° Semestre/Troncal/proyecto_troncal2/results/RNA_Seq"

# Crear carpeta si no existe
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# Leer datos del RNA-seq
data <- read.csv("C:/Users/fgarc/OneDrive/Escritorio/Doctorado/Ramos/1° Semestre/Troncal/proyecto_troncal2/data/RNA-Seq-expression-Norilsk2019.csv")

# Crear columna de significancia
data$significance <- "No significativo"
data$significance[data$log2FoldChange > 1 & data$padj < 0.05] <- "Sobreexpresado"
data$significance[data$log2FoldChange < -1 & data$padj < 0.05] <- "Subexpresado"

# --- Volcano Plot ---
pdf(file.path(output_path, "volcano_plot.pdf"))
ggplot(data, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("gray", "red", "blue")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "Volcano Plot RNA-seq Norilsk",
       x = "Log2 Fold Change",
       y = "-Log10 adjusted p-value (padj)")
dev.off()

# --- Genes sobreexpresados ---
genes_up <- data$gene_id[data$log2FoldChange > 1 & data$padj < 0.05]

# --- Anotación con biomaRt ---
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"),
                   filters = "ensembl_gene_id",
                   values = genes_up,
                   mart = ensembl)

# Filtrar genes con ID Entrez
gene_info <- gene_info[!is.na(gene_info$entrezgene_id), ]
write.csv(gene_info, file = file.path(output_path, "genes_sobreexpresados.csv"), row.names = FALSE)

# --- GO Enrichment ---
ego <- enrichGO(gene         = gene_info$entrezgene_id,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable     = TRUE)

# Guardar resultados GO
write.csv(as.data.frame(ego), file = file.path(output_path, "GO_enrichment_BP.csv"))

# Gráfico GO
pdf(file.path(output_path, "GO_BP_barplot.pdf"))
barplot(ego, showCategory = 10, title = "Top GO términos (BP) - Sobreexpresados")
dev.off()

# --- KEGG Enrichment ---
# KEGG usa código organismo (hsa = Homo sapiens)
ekegg <- enrichKEGG(gene         = gene_info$entrezgene_id,
                    organism     = 'hsa',
                    keyType      = 'kegg',
                    pvalueCutoff = 0.05)

# Hacer IDs legibles
ekegg_readable <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Guardar resultados KEGG
write.csv(as.data.frame(ekegg_readable), file = file.path(output_path, "KEGG_enrichment.csv"))

# Gráfico KEGG
pdf(file.path(output_path, "KEGG_barplot.pdf"))
barplot(ekegg_readable, showCategory = 10, title = "Top KEGG Pathways - Sobreexpresados")
dev.off()
