#Read data

metabo_data <- read.csv("C:/Users/irene/OneDrive/Documentos/metaboData/Datasets/2018-MetabotypingPaper/DataValues_S013.csv")
metadata_metabo <- read.csv("C:/Users/irene/OneDrive/Documentos/metaboData/Datasets/2018-MetabotypingPaper/DataInfo_S013.csv")


#Create Summarized Experimetn

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")

library(SummarizedExperiment)


col= colnames(metabo_data)
row= rownames(metabo_data)

se <- SummarizedExperiment(assays = list(counts = metabo_data),
                           colData = col,
                           rowData = row,
                           metadata = metadata_metabo)



save(se, file = "SummarizedExperiment.rda")


#Exploración general

dim(se)
head(colnames(se))
head(rownames(se))

metabolic_data <- assay(se, "counts")


summary_result <- summary(metabolic_data [1:6])
head(summary_result)


#Representación de variables

opt <- par(mfrow = c(3, 3))

for (i in 6:10) {
  hist(metabolic_data[, i], main = names(metabolic_data)[i])
}

par(opt)

#Género
freq_gender <- table(metabolic_data$GENDER)
freq_gender

barplot(freq_gender, xlab= "Género", ylab= "Frecuencia")

#Cirugía

freq_ciru <- table(metabolic_data$SURGERY)

barplot(freq_ciru, xlab= "Tipo de cirugía", ylab= "Frecuencia")

#Edad

hist(metabolic_data$AGE, xlab= "Edad", ylab= "Frecuencia")
max(metabolic_data$AGE)
min(metabolic_data$AGE)

#Quito las primeras columnas ya que no son sobre datos metabólicos
subset_data <- metabolic_data[, 10:696]

#Quito columnas con NAs, no numéricas y con solo ceros


#Quitar columnas no numéricas
numeric_data <- subset_data[, sapply(subset_data, is.numeric)]

#Quitar columnas con solo ceros
numeric_data <- numeric_data[, colSums(numeric_data != 0, na.rm = TRUE) > 0]

#Quitar columnas con NAs
numeric_data <- numeric_data[ , colSums(is.na(numeric_data))==0]

#Hacemos PCA

library(ggfortify)

pca_result <- prcomp(numeric_data, scale. = TRUE)

#Representamos PCA
autoplot(pca_result)

autoplot(pca_result, data = metabolic_data , colour = "SURGERY")

autoplot(pca_result, data = metabolic_data , colour = "AGE")

autoplot(pca_result, data = metabolic_data , colour = "GENDER")


library(ggplot2)

pca_plot <- autoplot(pca_result, data = metabolic_data, colour = "SUBJECTS") +
  geom_text(aes(label = SUBJECTS), size = 3, vjust = -0.5) +
  labs(title = "PCA Plot con etiquetas de sujetos")

pca_plot


library(factoextra)

var <- get_pca_var(pca_result)

head(var$contrib, 4)


#Contribución se las variables a PC1
fviz_contrib(pca_result, choice = "var", axes = 1, top = 10)

#Heatmap

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")


matrix <- as.matrix(numeric_data)
heatmap(matrix, scale = "column")

