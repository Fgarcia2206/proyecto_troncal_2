library(dplyr)

# Verificar duplicados
duplicated_genes <- rna_seq$Gene[duplicated(rna_seq$Gene)]
print(paste("Número de genes duplicados:", length(duplicated_genes)))

# Opción A: Conservar la entrada con menor p-valor ajustado (usando dplyr)
rna_seq_dedup <- rna_seq %>%
  arrange(padj) %>%  # Ordenar por significancia
  distinct(Gene, .keep_all = TRUE)  # Conservar solo la primera aparición

# Opción B: Alternativa con R base (si prefieres no usar dplyr)
rna_seq_dedup <- rna_seq[order(rna_seq$padj), ]
rna_seq_dedup <- rna_seq_dedup[!duplicated(rna_seq_dedup$Gene), ]





# Verificar que se eliminaron los duplicados
print(paste("Genes después de eliminar duplicados:", nrow(rna_seq_dedup)))

# Ver los genes que eran duplicados
print("Ejemplos de genes duplicados originales:")
print(head(duplicated_genes))