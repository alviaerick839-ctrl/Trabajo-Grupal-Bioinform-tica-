Pasos del analisis
# 02_Scripts/analysis_deseq2.R

Descripción: Script principal para realizar el Análisis Diferencial de Expresión (ADE) 
usando el paquete DESeq2 en R.

# -----------------------------------------------------------
# 1. Cargar datos y librerías 
# -----------------------------------------------------------
install.packages(DESeq2)
install.packages(ggplot2)
print("Librerías cargadas exitosamente.")

count_matrix <- read.table("../01_Data_Raw/gene_counts.tsv", header = TRUE, row.names = 1)
metadata <- read.table("../01_Data_Raw/samples_metadata.txt", header = TRUE, row.names = 1)

# -----------------------------------------------------------
# 2. Preprocesamiento y Construcción del Objeto DESeqDataSet (Simulado)
# -----------------------------------------------------------
 dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
colData = metadata, 
design = ~ Condición)
dds <- dds[rowSums(counts(dds)) > 1, ] # Filtrado de bajos conteos

# -----------------------------------------------------------
# 3. Ejecutar el Análisis (Simulado)
# -----------------------------------------------------------
dds <- DESeq(dds)
res <- results(dds, contrast=c("Condición", "Tumor", "Control"))

# 4. Guardar Resultados Clave
write.csv(as.data.frame(res), file="../03_Results/DE_results_table.csv")
saveRDS(dds, file="../03_Results/dds_object.rds")
