library(GEOquery)
# obtenemos el expressionSet
gse848 <- getGEO("GSE848",GSEMatrix=TRUE)
# obtenemos los datos de expresion
datExpr = exprs (gse848[[1]])
# cambiamos las columnas
targets <- pData(gse848[[1]])
colnames(datExpr) <- targets$title
# nos quedamos con los primeros 16
datExpr <- datExpr[,c(3:18)]
library("MultiBaC")
# Para ordenar los nombres de tu matriz en R por la última letra (A o B), puedes utilizar el siguiente código:
# Vector de nombres
nombres <- c("E2 8h A", "E2 8h B", "E2 48h A", "E2 48h B", "E2+ICI 8h A", "E2+ICI 8h B", 
             "E2+ICI 48h A", "E2+ICI 48h B", "E2+Ral 8h A", "E2+Ral 8h B", "E2+Ral 48h A", "E2+Ral 48h B",
             "E2+TOT 8h A", "E2+TOT 8h B", "E2+TOT 48h A", "E2+TOT 48h B")

# Ordenar por la última letra (A o B)
nombres_ordenados <- nombres[order(sub(".* ", "", nombres))]
# Renombrar en el formato A.E2_8h
nombres_renombrados <- gsub("([AB])$", "\\1", nombres_ordenados) # Coloca A. o B.
nombres_renombrados <- gsub(" ", "_", nombres_renombrados)        # Reemplaza espacio con _
nombres_renombrados <- gsub("\\+", "", nombres_renombrados)       # Elimina el '+'
# Reordenar la matriz datExpr usando los nombres de columna ordenados
colnames(datExpr) <- nombres_renombrados
# separamos en batch A y B
A.datexpr <- datExpr[,c(1:8)]
B.datexpr <- datExpr[,c(9:16)]
# Creamos MultiAssayExperiment
data_mbac<- createMbac(inputOmics = list(A.datexpr, B.datexpr), 
                       batchFactor = c("A", "B"),
                       experimentalDesign = list("A" =  c("E2_8h, E2_48h, E2ICI_8h, E2ICI_48h, E2Ral_8h, E2Ral_48h, E2TOT_8h, E2TOT_48h"),
                                                 "B" = c("E2_8h, E2_48h, E2ICI_8h, E2ICI_48h, E2Ral_8h, E2Ral_48h, E2TOT_8h, E2TOT_48h")),
                       omicNames = "RNA")
# to remove batch effects or unwanted noise from a single omic data matrix in the MultiBaC package is the ARSyNbac function
par(mfrow = c(1,2))
arsyn_4 <- ARSyNbac(data_mbac, modelName = "datexpr", beta = 0.5, 
                    batchEstimation = FALSE, filterNoise = TRUE,
                    Interaction = FALSE,
                    Variability = 0.95)
plot(arsyn_4, typeP="pca.cor", bty = "L",
     pch = custom_pch, cex = 3, col.per.group = custom_col,
     legend.text = c("Color: Batch", names(data_mbac$ListOfBatches),
                     "Fill: Cond.", unique(colData(data_mbac$ListOfBatches$A)$tfactor)),
     args.legend = list("x" = "topright",
                        "pch" = c(NA, 15, 15, 15, 
                                  NA, 15, 0),
                        "col" = c(NA, "brown2", "dodgerblue", "forestgreen",
                                  NA, 1, 1),
                        "bty" = "n",
                        "cex" = 1.2))