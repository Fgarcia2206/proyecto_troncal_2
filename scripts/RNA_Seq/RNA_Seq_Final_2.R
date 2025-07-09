### INSTALACIÓN DE PAQUETES ####
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c("org.Hs.eg.db", "biomaRt", "goseq", "AnnotationDbi"))
install.packages("tidyverse")
install.packages("ggrepel")

if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) {
  BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
}
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot")


### CARGA DE LIBRERÍAS ####
library(tidyverse)
library(org.Hs.eg.db)
library(biomaRt)
library(goseq)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(forcats)

### DEFINIR RUTA DE SALIDA ####
output_dir <- "C:/Users/fgarc/OneDrive/Escritorio/Doctorado/Ramos/1° Semestre/Troncal/proyecto_troncal2/results/RNA_Seq/RNA_Seq_Final_2"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


### LECTURA Y PREPROCESAMIENTO DE DATOS ####
rna_data <- read.csv("C:/Users/fgarc/OneDrive/Escritorio/Doctorado/Ramos/1° Semestre/Troncal/proyecto_troncal2/data/RNA-Seq-expression-Norilsk2019.csv")
rna_data <- rna_data[!duplicated(rna_data$Gene), ]  # Eliminar duplicados


################################################################################
### VOLCANO PLOT SIN ETIQUETAS + TABLA TOP 10 UP / DOWN ###
################################################################################

fc_threshold <- 2
p_threshold <- 0.05

# Clasificar genes según expresión diferencial
rna_data <- rna_data %>%
  mutate(diffexpressed = case_when(
    log_2.fold.change >= fc_threshold & Adjusted.p.value <= p_threshold ~ "UP",
    log_2.fold.change <= -fc_threshold & Adjusted.p.value <= p_threshold ~ "DOWN",
    TRUE ~ "NO"
  )) %>%
  arrange(Adjusted.p.value)

# Seleccionar los top 10 genes UP y DOWN
top_up <- rna_data %>%
  filter(diffexpressed == "UP") %>%
  arrange(Adjusted.p.value) %>%
  slice(1:10)

top_down <- rna_data %>%
  filter(diffexpressed == "DOWN") %>%
  arrange(Adjusted.p.value) %>%
  slice(1:10)

# Guardar tabla combinada
top_genes_table <- bind_rows(
  top_up %>% mutate(Regulacion = "UP"),
  top_down %>% mutate(Regulacion = "DOWN")
)

# Guardar como archivo CSV
write.csv(top_genes_table, file.path(output_dir, "top_genes_up_down.csv"), row.names = FALSE)

# Crear Volcano Plot SIN etiquetas
volcano_plot <- ggplot(rna_data, aes(x = log_2.fold.change, y = -log10(Adjusted.p.value), col = diffexpressed)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey50", "UP" = "red")) +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot - Fiebre de Norilsk",
       subtitle = "Genes diferencialmente expresados (|log2FC| > 2, padj < 0.05)",
       x = expression(Log[2]~"Fold Change"),
       y = expression(-Log[10]~"Valor p ajustado")) +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        legend.position = "bottom")

# Guardar imagen
ggsave(file.path(output_dir, "volcano_plot_norilsk_clean.png"), plot = volcano_plot, width = 10, height = 6)


################################################################################
### ANÁLISIS DE ENRIQUECIMIENTO GO (genes sobreexpresados - top 10 términos) ###
################################################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(forcats)

# Enriquecimiento GO para genes sobreexpresados (upregulated)
ego_all <- enrichGO(
  gene          = up_entrez_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# Simplificar términos GO redundantes
ego_simplified <- simplify(ego_all, cutoff = 0.7, by = "p.adjust", select_fun = min)

# Guardar todos los términos simplificados en CSV
write.csv(as.data.frame(ego_simplified),
          file = file.path(output_dir, "GO_enrichment_upregulated_simplified.csv"),
          row.names = FALSE)

# Extraer dataframe con resultados
ego_df <- as.data.frame(ego_simplified)

if (nrow(ego_df) > 0) {
  
  # Función para graficar top 10 y guardar gráfico por ontología
  plot_go_ontology <- function(df, ont, output_dir) {
    df_sub <- df %>%
      filter(ONTOLOGY == ont) %>%
      arrange(p.adjust) %>%
      slice_head(n = 10) %>%
      mutate(Description = factor(Description, levels = rev(unique(Description))))
    
    if (nrow(df_sub) == 0) {
      message(paste("No hay términos GO significativos para ontología", ont))
      return(NULL)
    }
    
    p <- ggplot(df_sub, aes(x = Count, y = Description)) +
      geom_point(aes(size = Count, color = -log10(p.adjust))) +
      scale_color_gradient(low = "blue", high = "red", name = expression(-log[10](p.adjust))) +
      scale_size_continuous(name = "Número de genes") +
      theme_minimal(base_size = 14) +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
      ) +
      labs(
        title = paste("Top 10 términos GO:", ont),
        x = "Número de genes"
      )
    
    # Guardar gráfico
    ggsave(
      filename = file.path(output_dir, paste0("GO_dotplot_upregulated_top10_", ont, ".png")),
      plot = p,
      width = 9, height = 6
    )
  }
  
  # Generar gráficos por cada ontología
  plot_go_ontology(ego_df, "BP", output_dir)
  plot_go_ontology(ego_df, "CC", output_dir)
  plot_go_ontology(ego_df, "MF", output_dir)
  
} else {
  message("No se encontraron términos GO enriquecidos significativos para genes sobreexpresados.")
}



################################################################################
### ANÁLISIS DE ENRIQUECIMIENTO KEGG ####
################################################################################

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