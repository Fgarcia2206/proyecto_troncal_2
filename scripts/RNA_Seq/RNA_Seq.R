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


### DEFINIR RUTA DE SALIDA ####
output_dir <- "C:/Users/fgarc/OneDrive/Escritorio/Doctorado/Ramos/1° Semestre/Troncal/proyecto_troncal2/results/RNA_Seq/RNA_Seq_Final"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


### LECTURA Y PREPROCESAMIENTO DE DATOS ####
rna_data <- read.csv("C:/Users/fgarc/OneDrive/Escritorio/Doctorado/Ramos/1° Semestre/Troncal/proyecto_troncal2/data/RNA-Seq-expression-Norilsk2019.csv")

# Eliminar genes duplicados
rna_data <- rna_data[!duplicated(rna_data$Gene), ]



################################################################################
### VOLCANO PLOT MEJORADO (versión actualizada) ####
################################################################################

# Parámetros (ajustados para coincidir con tu preferencia)
fc_threshold <- 2  # |log2FC| > 2 (antes era 1)
p_threshold <- 0.05

# Clasificar según significancia (adaptado a tus nombres de columnas)
rna_data <- rna_data %>%
  mutate(diffexpressed = case_when(
    log_2.fold.change >= fc_threshold & Adjusted.p.value <= p_threshold ~ "UP",
    log_2.fold.change <= -fc_threshold & Adjusted.p.value <= p_threshold ~ "DOWN",
    TRUE ~ "NO"
  )) %>%
  arrange(Adjusted.p.value)

# Seleccionar genes para etiquetar (los más significativos)
rna_data$delabel <- NA
rna_data$delabel[rna_data$diffexpressed != "NO"] <- rna_data$Gene[rna_data$diffexpressed != "NO"]

# Volcano plot con tu estilo preferido
volcano_plot <- ggplot(rna_data, aes(x = log_2.fold.change, y = -log10(Adjusted.p.value), 
                                     col = diffexpressed, label = delabel)) +
  geom_point(aes(color = diffexpressed), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey50", "UP" = "red"),
                     name = "Expression",
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
             col = "black", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(p_threshold), 
             col = "black", linetype = "dashed", linewidth = 0.5) +
  labs(title = "Volcano Plot - Fiebre de Norilsk",
       subtitle = "Genes diferencialmente expresados (|log2FC| > 2, padj < 0.05)",
       x = expression(Log[2]~"Fold Change"),
       y = expression(-Log[10]~"Valor p ajustado")) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  geom_text_repel(
    data = subset(rna_data, abs(log_2.fold.change) > 4 & Adjusted.p.value < 1e-50),
    aes(label = Gene),
    size = 3,
    box.padding = 0.5,
    max.overlaps = 20,
    segment.color = 'grey50',
    segment.size = 0.2,
    min.segment.length = 0.1
  ) +
  scale_y_continuous(
    limits = c(0, 400),  # Límite superior ajustado a 400 como en tu ejemplo
    expand = expansion(mult = c(0, 0.05))
  )

# Guardar el plot
ggsave(filename = file.path(output_dir, "volcano_plot_norilsk_improved.png"), 
       plot = volcano_plot, width = 10, height = 8, dpi = 300, bg = "white")

# Mostrar el plot
print(volcano_plot)



#################################################################################
### ANÁLISIS DE ENRIQUECIMIENTO FUNCIONAL - VERSIÓN ROBUSTA ####
#################################################################################

# 1. Preparación de genes sobreexpresados --------------------------------------
up_genes <- rna_data %>% 
  filter(diffexpressed == "UP") %>% 
  pull(Gene)

cat("\nPreparación de genes para análisis funcional:\n")
cat("- Total de genes sobreexpresados:", length(up_genes), "\n")

# 2. Mapeo ENSEMBL a ENTREZID con diagnóstico mejorado -------------------------
cat("\nIniciando mapeo ENSEMBL a ENTREZID...\n")

ensembl2entrez <- tryCatch({
  # Mapeo incluyendo SYMBOL para diagnóstico
  result <- bitr(up_genes, 
                 fromType = "ENSEMBL", 
                 toType = c("ENTREZID", "SYMBOL"),
                 OrgDb = org.Hs.eg.db)
  
  # Verificación de resultados
  if(nrow(result) == 0) stop("El mapeo devolvió 0 filas")
  
  # Reporte detallado
  unmapped <- setdiff(up_genes, result$ENSEMBL)
  mapping_rate <- round(nrow(result)/length(up_genes)*100, 1)
  
  cat("\nResultados del mapeo:\n")
  cat("- Genes de entrada:", length(up_genes), "\n")
  cat("- Genes mapeados:", nrow(result), "\n")
  cat("- Porcentaje de éxito:", mapping_rate, "%\n")
  
  if(length(unmapped) > 0){
    cat("\nGenes no mapeados (primeros 10):\n")
    print(head(unmapped, 10))
    
    # Guardar lista de no mapeados
    write.table(data.frame(ENSEMBL=unmapped),
                file.path(output_dir, "unmapped_genes.txt"),
                row.names = FALSE, quote = FALSE)
  }
  
  # Filtrado riguroso
  result <- result[!is.na(result$ENTREZID) & result$ENTREZID != "", ]
  cat("\nENTREZIDs válidos después de filtrar NAs:", nrow(result), "\n")
  
  # Verificar duplicados
  if(any(duplicated(result$ENSEMBL))){
    dup_genes <- result$ENSEMBL[duplicated(result$ENSEMBL)]
    cat("\nAdvertencia: Hay", length(dup_genes), "genes ENSEMBL con múltiples ENTREZIDs\n")
    print(result[result$ENSEMBL %in% dup_genes, ])
  }
  
  result
}, error = function(e) {
  cat("\nError en el mapeo ENSEMBL→ENTREZID:\n", e$message, "\n")
  return(NULL)
})

# Verificación crítica
if(is.null(ensembl2entrez) || nrow(ensembl2entrez) == 0){
  stop("Fallo crítico: No hay genes mapeados para continuar con el análisis")
}

# 3. Preparación de IDs ENTREZ para análisis -----------------------------------
up_entrez_ids <- unique(ensembl2entrez$ENTREZID)
up_entrez_chr <- as.character(up_entrez_ids)  # Conversión explícita a character

cat("\nResumen de genes para análisis funcional:\n")
cat("- ENTREZIDs únicos:", length(up_entrez_chr), "\n")
cat("- Ejemplo de IDs (primeros 5):", paste(head(up_entrez_chr, 5), collapse = ", "), "\n")

# 4. Análisis de Enriquecimiento GO -------------------------------------------
cat("\nIniciando análisis de enriquecimiento GO...\n")

ego_all <- enrichGO(gene          = up_entrez_ids,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = "ENTREZID",
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.2,
                    readable      = TRUE)

# Procesamiento de resultados GO
if(!is.null(ego_all) && nrow(ego_all) > 0){
  ego_simplified <- simplify(ego_all, 
                             cutoff = 0.7, 
                             by = "p.adjust", 
                             select_fun = min)
  
  # Guardado de resultados
  write.csv(as.data.frame(ego_simplified),
            file.path(output_dir, "GO_enrichment_results.csv"),
            row.names = FALSE)
  
  # Visualización mejorada
  if(nrow(ego_simplified) > 0){
    # Dotplot por ontologías
    go_dotplot <- dotplot(ego_simplified, 
                          split = "ONTOLOGY",
                          showCategory = 10,
                          font.size = 10) + 
      facet_grid(ONTOLOGY ~ ., scales = "free") +
      scale_color_gradient(low = "red", high = "blue") +
      labs(title = "Enriquecimiento de Términos GO") +
      theme_minimal(base_size = 12) +
      theme(axis.text.y = element_text(size = 9))
    
    ggsave(file.path(output_dir, "GO_dotplot_improved.png"), 
           go_dotplot, width = 12, height = 10, dpi = 300)
    
    # Red de genes y términos GO (solo si no hay demasiados genes)
    if(length(up_entrez_ids) <= 500){
      go_cnet <- cnetplot(ego_simplified,
                          categorySize = "pvalue",
                          showCategory = 5,
                          circular = FALSE,
                          colorEdge = TRUE,
                          node_label = "all") +
        labs(title = "Red de Genes y Términos GO") +
        theme_minimal()
      
      ggsave(file.path(output_dir, "GO_cnetplot_improved.png"), 
             go_cnet, width = 14, height = 10, dpi = 300)
    }
  } else {
    message("No se encontraron términos GO enriquecidos significativos después de simplificar.")
  }
} else {
  message("El análisis GO no produjo resultados significativos.")
}

# 5. Análisis de Enriquecimiento KEGG -----------------------------------------
cat("\nIniciando análisis de enriquecimiento KEGG...\n")

ekegg <- enrichKEGG(gene         = up_entrez_chr,
                    organism     = "hsa",
                    pvalueCutoff = 0.1,
                    minGSSize    = 10,
                    maxGSSize    = 500)

# Procesamiento de resultados KEGG
if(!is.null(ekegg) && nrow(ekegg) > 0){
  # Conversión a nombres legibles
  ekegg <- setReadable(ekegg, 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "ENTREZID")
  
  # Guardado de resultados
  write.csv(as.data.frame(ekegg),
            file.path(output_dir, "KEGG_enrichment_results.csv"),
            row.names = FALSE)
  
  # Visualización mejorada
  kegg_plot <- dotplot(ekegg,
                       showCategory = 15,
                       font.size = 12,
                       color = "p.adjust") +
    scale_color_gradient(low = "red", high = "blue") +
    labs(title = "Rutas KEGG Enriquecidas (p.adj < 0.1)") +
    theme_minimal()
  
  ggsave(file.path(output_dir, "KEGG_dotplot.png"), 
         kegg_plot, width = 12, height = 8, dpi = 300)
  
  # Enlace al pathway más significativo
  if(nrow(ekegg) > 0){
    top_pathway <- ekegg$ID[1]
    cat("\nPathway KEGG más significativo:", top_pathway, "\n")
    cat("Enlace directo: https://www.kegg.jp/kegg-bin/show_pathway?", top_pathway, "\n", sep="")
  }
} else {
  message("\nNo se encontraron rutas KEGG enriquecidas significativas.")
  
  # Diagnóstico avanzado
  cat("\nDiagnóstico para KEGG:\n")
  cat("- Genes enviados a KEGG:", length(up_entrez_chr), "\n")
  
  # Verificar conexión con KEGG
  if(requireNamespace("KEGGREST", quietly = TRUE)){
    # Verificar si los genes existen en KEGG
    kegg_genes <- unique(KEGGREST::keggLink("hsa", "gene"))
    matched_genes <- sum(up_entrez_chr %in% gsub("hsa:", "", kegg_genes))
    cat("- Genes encontrados en KEGG:", matched_genes, "/", length(up_entrez_chr), "\n")
    
    if(matched_genes == 0){
      cat("\nPosibles causas:\n")
      cat("1. Los IDs ENTREZ no coinciden con la base de datos KEGG\n")
      cat("2. Problema de conexión con el servidor KEGG\n")
      cat("3. Los genes no están asociados a rutas conocidas\n")
      
      # Ejemplo de verificación manual
      cat("\nEjemplo de verificación manual (usando KEGGREST):\n")
      if(length(up_entrez_chr) > 0){
        test_gene <- up_entrez_chr[1]
        cat("- Gene:", test_gene, "\n")
        try({
          pathways <- KEGGREST::keggGet(paste0("hsa:", test_gene))
          if(length(pathways) > 0){
            cat("- Encontrado en:", paste(names(pathways[[1]]$PATHWAY), collapse = ", "), "\n")
          } else {
            cat("- No encontrado en KEGG\n")
          }
        }, silent = TRUE)
      }
    }
  } else {
    cat("\nInstale KEGGREST para diagnóstico avanzado:\n")
    cat('BiocManager::install("KEGGREST")\n')
  }
}

#################################################################################
### RESUMEN FINAL DE RESULTADOS ####
#################################################################################

cat("\n\nRESUMEN FINAL DEL ANÁLISIS\n")
cat("=================================\n")
cat("Genes analizados:", nrow(rna_data), "\n")
cat("Genes sobreexpresados:", sum(rna_data$diffexpressed == "UP"), "\n")
cat("Genes subexpresados:", sum(rna_data$diffexpressed == "DOWN"), "\n")

cat("\nResultados de mapeo:\n")
cat("- Genes ENSEMBL de entrada:", length(up_genes), "\n")
cat("- Genes mapeados a ENTREZID:", length(up_entrez_ids), "\n")
cat("- Porcentaje de éxito:", round(length(up_entrez_ids)/length(up_genes)*100, 1), "%\n")

if(exists("ego_simplified")){
  cat("\nResultados GO:\n")
  if(nrow(ego_simplified) > 0){
    cat("- Términos GO significativos:", nrow(ego_simplified), "\n")
    cat("- Términos BP:", sum(ego_simplified$ONTOLOGY == "BP"), "\n")
    cat("- Términos MF:", sum(ego_simplified$ONTOLOGY == "MF"), "\n")
    cat("- Términos CC:", sum(ego_simplified$ONTOLOGY == "CC"), "\n")
  } else {
    cat("- No se encontraron términos GO significativos\n")
  }
}

if(exists("ekegg")){
  cat("\nResultados KEGG:\n")
  if(nrow(as.data.frame(ekegg)) > 0){
    cat("- Rutas KEGG significativas:", nrow(as.data.frame(ekegg)), "\n")
    cat("- Pathway más significativo:", ekegg$Description[1], "(p.adjust =", ekegg$p.adjust[1], ")\n")
  } else {
    cat("- No se encontraron rutas KEGG significativas\n")
  }
}

cat("\nArchivos guardados en:", output_dir, "\n")
cat("=================================\n")