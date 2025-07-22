# 1. Instalación y carga de paquetes -------------------------------------------------

# Instalar paquetes de Bioconductor si no están instalados
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("ggplot2", "ggrepel", "clusterProfiler", "org.Hs.eg.db", "goseq", "biomaRt"))
install.packages("dplyr")

# Cargar los paquetes
library(ggplot2)
library(ggrepel)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(goseq)
library(biomaRt)

# 2. Configuración de rutas ----------------------------------------------------------

# Definir la ruta base donde se guardarán los resultados
base_path <- "C:/Users/fgarc/OneDrive/Escritorio/Doctorado/Ramos/1° Semestre/Troncal/proyecto_troncal2/results/RNA_Seq/RNA_Seq_DeepSeek"

# Crear el directorio si no existe
if (!dir.exists(base_path)) {
  dir.create(base_path, recursive = TRUE)
}

# 3. Cargar y preparar los datos ----------------------------------------------------

# Leer los datos (ajusta la ruta según donde tengas el archivo)
rna_seq <- read.csv("C:/Users/fgarc/OneDrive/Escritorio/Doctorado/Ramos/1° Semestre/Troncal/proyecto_troncal2/data/RNA-Seq-expression-Norilsk2019.csv", stringsAsFactors = FALSE)

# Eliminar duplicados manteniendo el registro con el menor p-value
rna_seq <- rna_seq %>%
  group_by(Gene) %>%
  arrange(Adjusted.p.value) %>%
  slice(1) %>%
  ungroup()

# Convertir a data.frame
rna_seq <- as.data.frame(rna_seq)

# Renombrar columnas para facilitar el trabajo
colnames(rna_seq) <- c("Gene", "Species", "log2FC", "padj")

# Verificar estructura
str(rna_seq)

# 4. Volcano plot -------------------------------------------------------------------

# Definir umbrales (ajustables según necesidad)
fc_threshold <- 2  # |log2FC| > 2
p_threshold <- 0.05  # padj < 0.05

# Crear columnas para el volcano plot
rna_seq$diffexpressed <- "NO"
rna_seq$diffexpressed[rna_seq$log2FC > fc_threshold & rna_seq$padj < p_threshold] <- "UP"
rna_seq$diffexpressed[rna_seq$log2FC < -fc_threshold & rna_seq$padj < p_threshold] <- "DOWN"

# Seleccionar algunos genes para etiquetar (los más significativos)
rna_seq$delabel <- NA
rna_seq$delabel[rna_seq$diffexpressed != "NO"] <- rna_seq$Gene[rna_seq$diffexpressed != "NO"]

# Volcano plot con ajuste del eje Y
volcano_plot <- ggplot(data = rna_seq, aes(x = log2FC, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_point(aes(color = diffexpressed), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey50", "UP" = "red"),
                     name = "Expression",
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), col = "black", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(p_threshold), col = "black", linetype = "dashed", linewidth = 0.5) +
  labs(title = "Volcano Plot - Expresión diferencial en Fiebre de Norilsk",
       subtitle = "Genes diferencialmente expresados (|log2FC| > 2, padj < 0.05)",
       x = "log2 Fold Change",
       y = "-log10(Adjusted p-value)") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  geom_text_repel(
    data = subset(rna_seq, abs(log2FC) > 4 & padj < 1e-50),
    aes(label = Gene),
    size = 3,
    box.padding = 0.5,
    max.overlaps = 20,
    segment.color = 'grey50',
    segment.size = 0.2,
    min.segment.length = 0.1
  ) +
  scale_y_continuous(
    limits = c(0, 400),  # Límite superior a 500
    expand = expansion(mult = c(0, 0.05))
  )

# Mostrar el plot
print(volcano_plot)

# Guardar el plot
ggsave(file.path(base_path, "volcano_plot_Norilsk.png"), 
       plot = volcano_plot, width = 10, height = 8, dpi = 300, bg = "white")






# 5. Análisis de genes sobreexpresados (Versión corregida) -------------------------

# Filtrar genes significativamente sobreexpresados
up_genes <- rna_seq[rna_seq$diffexpressed == "UP", "Gene"]

# Solución 1: Usar el archivo de anotación más reciente
if (!require("biomaRt")) install.packages("biomaRt")
library(biomaRt)

# Intentar conectar al servidor más reciente de Ensembl
tryCatch({
  mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
}, error = function(e) {
  # Si falla, intentar con el archivo local
  message("Error al conectar con Ensembl, usando versión local...")
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "http://grch37.ensembl.org")
})

# Convertir IDs ENSEMBL a Entrez ID
gene_ids <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                  filters = "ensembl_gene_id",
                  values = up_genes,
                  mart = mart)

# Solución alternativa si falla biomaRt (usando org.Hs.eg.db)
if (nrow(gene_ids) == 0) {
  if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
  library(org.Hs.eg.db)
  
  gene_ids <- mapIds(org.Hs.eg.db,
                     keys = up_genes,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")
  gene_ids <- data.frame(ensembl_gene_id = names(gene_ids),
                         entrezgene_id = unname(gene_ids))
}

# Continuar con el análisis...
pwf <- nullp(gene_vector, bias.data = NULL, plot.fit = FALSE)
go_results <- goseq(pwf, gene2cat = NULL, method = "Wallenius")

# Usar goseq para análisis de enriquecimiento GO
pwf <- nullp(gene_vector, bias.data = NULL, plot.fit = FALSE)
go_results <- goseq(pwf, gene2cat = NULL, method = "Wallenius")

# Filtrar resultados significativos (FDR < 0.05)
significant_go <- go_results[go_results$over_represented_pvalue < 0.05,]

# Anotar los términos GO
significant_go$term <- Term(GOTERM)[significant_go$category]
significant_go$ontology <- Ontology(GOTERM)[significant_go$category]

# Ver los resultados
head(significant_go[order(significant_go$over_represented_pvalue),])

# Guardar resultados en la ruta especificada
write.csv(significant_go, 
          file.path(base_path, "GO_enrichment_results.csv"), 
          row.names = FALSE)

# Análisis de rutas KEGG
kegg_results <- goseq(pwf, gene2cat = NULL, test.cats = "KEGG", method = "Wallenius")
significant_kegg <- kegg_results[kegg_results$over_represented_pvalue < 0.05,]

# Ver resultados KEGG
head(significant_kegg[order(significant_kegg$over_represented_pvalue),])

# Guardar resultados KEGG en la ruta especificada
write.csv(significant_kegg, 
          file.path(base_path, "KEGG_pathway_results.csv"), 
          row.names = FALSE)

# 6. Visualización de resultados GO -------------------------------------------------

# Filtrar los 10 términos GO más significativos por ontología
top_go_bp <- significant_go %>%
  filter(ontology == "BP") %>%
  arrange(over_represented_pvalue) %>%
  head(10)

top_go_mf <- significant_go %>%
  filter(ontology == "MF") %>%
  arrange(over_represented_pvalue) %>%
  head(10)

top_go_cc <- significant_go %>%
  filter(ontology == "CC") %>%
  arrange(over_represented_pvalue) %>%
  head(10)

# Función para crear plots GO
create_go_plot <- function(data, title) {
  ggplot(data, aes(x = -log10(over_represented_pvalue), y = reorder(term, -log10(over_represented_pvalue)))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = title,
         x = "-log10(p-value)",
         y = "GO Term") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
}

# Crear plots
go_bp_plot <- create_go_plot(top_go_bp, "Top 10 Biological Processes")
go_mf_plot <- create_go_plot(top_go_mf, "Top 10 Molecular Functions")
go_cc_plot <- create_go_plot(top_go_cc, "Top 10 Cellular Components")

# Mostrar plots
print(go_bp_plot)
print(go_mf_plot)
print(go_cc_plot)

# Guardar plots en la ruta especificada
ggsave(file.path(base_path, "GO_BP_plot.png"), go_bp_plot, width = 10, height = 6)
ggsave(file.path(base_path, "GO_MF_plot.png"), go_mf_plot, width = 10, height = 6)
ggsave(file.path(base_path, "GO_CC_plot.png"), go_cc_plot, width = 10, height = 6)

# 7. Mensaje final ------------------------------------------------------------------

cat("Análisis completado. Todos los resultados se han guardado en:\n", base_path)