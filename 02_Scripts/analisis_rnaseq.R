#############################################
# Script básico para analizar datos RNA-seq
# Cargar datos, normalizar y graficar
#############################################

# 1. Cargar librerías
library(tidyverse)

# 2. Cargar archivo (cámbialo por tu ruta)
datos <- read.table("conteos.txt", header = TRUE, sep = "\t")

# 3. Convertir GeneID a fila
rownames(datos) <- datos$GeneID
datos <- datos[ , -1]

# 4. Exploración inicial
cat("Dimensiones:", dim(datos), "\n")
summary(datos)

# 5. Normalización CPM sencilla
cpm <- sweep(datos, 2, colSums(datos) / 1e6, "/")

# 6. Graficar distribución de muestras
boxplot(log2(cpm + 1),
        main = "Distribución de expresión (log2 CPM)",
        xlab = "Muestras", ylab = "Expresión")

# 7. Calcular genes más expresados
top_genes <- head(sort(rowMeans(cpm), decreasing = TRUE), 10)
print("Top 10 genes más expresados:")
print(top_genes)

# 8. Guardar matriz CPM
write.table(cpm, "cpm_resultados.txt", sep="\t", quote=FALSE)

#############################################
# Fin del script
#############################################
