#==============================================================================#
#                                   PROJET "FRUITS"                            #
#                                ANALYSES MULTIVARIEES                         #
#==============================================================================#

# Se rendre dans le répertoire contenant le jeu de données
#=========================================================

setwd("/home/richard/Bureau/Bio-statistiques/multivariées/")

#==============================================================================#
#                           1. ANALYSE EXPLORATOIRE                            #
#==============================================================================#

# Importer le jeu de données
#===========================

tfruits <- read.table("teneurs_fruits.csv",
                 sep = ",",
                 dec = ",",
                 header = T
                 #row.names = T,
                 #check.names = F,
                 #fill = TRUE
                 )

# Virer les 0 après les décimales
#================================

#tfruits <- as.numeric(format(tfruits, drop0trailing = T))

# Formater le nom des colonnes
#=============================

#colClean <- function(x){ colnames(x) <- gsub("..", ".", colnames(x)); x }
colnames(tfruits)[colnames(tfruits) == 'X'] <- "especes"
colnames(tfruits) <- gsub(x=colnames(tfruits), pattern = "\\.$", replacement = "")
colnames(tfruits) <- gsub(x=colnames(tfruits), pattern = "\\.\\.", replacement = "\\.")

# Convertir les colonnes qui sont en "g" en "mg"
#===============================================

# Repérer les colonnes en "g" au lieu de "mg" (2 façons)
col_grammes <- grep(".g", names(tfruits), value = T)
col_grammes <- grep(".g", names(tfruits), value = T, fixed = T)

for (i in 1:length(col_grammes)) {
  # Multiplier par 1000 pour convertir en mg
  tfruits[col_grammes[i]] <- get(col_grammes[i], tfruits) * 1000
  # changer nom de la colonne
  colnames(tfruits)[match(colnames(tfruits[col_grammes[i]]), colnames(tfruits))] <- sub(".g", ".mg", col_grammes[i])
}

# Supprimer les colonnes "joules.kJ" et Calories.kcal
#====================================================

tfruits$joules.kJ <- NULL
tfruits$Calories.kcal <- NULL

#tfruits <- subset(tfruits, select = -c(especes, Calories.kcal))

# Voir matrice et description
#============================

View(tfruits)

dim(tfruits) #72 fruits

# Statistiques descriptives par variable
#=======================================

summary(tfruits)

# Distribution de chacune des variables
#======================================

library(ggplot2)

distribution <- function(df, figure = TRUE) {
  for (i in 1:ncol(df)) {
    if (is.numeric(df[,i]) == TRUE) {
      unite <- sub("(^.*)\\.(.*$)", "\\2", colnames(df[i]))
      variable <- sub("(^.*)\\.(.*$)", "\\1", colnames(df[i]))
      p <- ggplot(data = df, aes(as.character(x=especes), y = df[,i])) +
        geom_col(color="steelblue", fill="steelblue") +
        ggtitle(paste0("Distribution des quantités de ", variable, " au sein de ", nrow(df)," fruits")) +
        labs(x="Espèces", y = unite) +
        theme(axis.text.x = element_text(angle=45, hjust = 1, size = 6))
      print(p)
    }
    else {
      #cat(paste0("la variable ", colnames(df[i]), " n'est pas numérique\n"))
      message("la variable ", colnames(df[i]), " n'est pas numérique !!!\n")
    }
  }
  }
      
distribution(tfruits)

# Matrice de corrélation
#=======================

library(dplyr)

#num_tfruit <- select_if(tfruits, is.numeric)
num_tfruits <- tfruits %>% select_if(is.numeric)
rownames(num_tfruits) <- tfruits$especes

corr_tfruit <- cor(num_tfruits)

library(ggcorrplot)

#ggcorrplot(corr_tfruit, method = "square")
# option hc.order = TRUE ?
ggcorrplot(corr_tfruit,
           method = "circle",
           outline.col = "white",
           ggtheme = ggplot2::theme_light(),
           lab = T,
           insig = "blank",
           colors = c("blue", "white", "red")
           ) +
  ggtitle("Matrice de corrélation") +
  theme(
    panel.background = element_rect(fill = "lightblue",
                                    colour = "lightblue",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "white"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "white")
  )

# Scatter plot matrix
#====================

library("GGally")

ggpairs(num_tfruits, columns = 1:ncol(num_tfruits), title = "Scatter plot matrix, distribution et coefficient de Pearson")

#ggscatmat(tfruits, columns = 1:ncol(num_tfruits), alpha = 0.1) +
#  ggtitle("Scatter plot matrix, distribution et coefficient de Pearson")

#==============================================================================#
#                                   2. ACP                                     #
#==============================================================================#

#test <- prcomp(num_tfruits, center=T,scale.=T)
#biplot(prcomp(num_tfruits, center=T,scale.=T))


library(FactoMineR)
library(factoextra)

# Calculer l’ACP
#===============

res.pca <- PCA(num_tfruits, scale.unit = TRUE, graph = FALSE)

# Voir une interprétation automatique des résultats
#==================================================
library(FactoInvestigate)
Investigate(res.pca)

# Accéder aux résultats de l'ACP
#===============================

# Valeurs propres

res.pca$eig

# Résultats des variables

res.var <- res.pca$var
res.var$coord          # Coordonnées
res.var$contrib        # Contributions aux axes
res.var$cos2           # Qualité de représentation

# Résultats des individus

res.ind <- res.pca$ind
res.ind$coord          # Coordonnées
res.ind$contrib        # Contributions aux axes
res.ind$cos2           # Qualité de représentation


# Visualisation des valeurs propres (Eigenvalues)
#================================================

eig.val <- res.pca$eig

barplot(eig.val[, 2], 
        names.arg = 1:nrow(eig.val), 
        main = "Variances Expliquées par Composante Principale (%)",
        xlab = "Composantes Principales",
        ylab = "Pourcnetage des variances",
        col ="steelblue")

# Add connected line segments to the plot
lines(x = 1:nrow(eig.val), eig.val[, 2], 
      type = "b", pch = 19, col = "red")

# Avec Factoextra
eig.val <- get_eigenvalue(res.pca)

fviz_eig(res.pca,
         addlabels = TRUE,
         ylim = c(0, max(res.pca$eig[,2])),
         linecolor = "red") +
  labs(title = "Décomposition de l'inertie totale des axes",
       subtitle = "Histogramme des valeurs propres"
       ) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
        )

# Graphiques des variables
#========================

# Cercle des corrélations
#------------------------

plot(res.pca, choix = "var", autoLab = "yes")

# Avec Factoextra
fviz_pca_var(res.pca, col.var = "black", repel = TRUE) +
  ggtitle("ACP - variables") +
  theme(plot.title = element_text(hjust = 0.5))

# Qualité de représentation des variables
#----------------------------------------

var <- get_pca_var(res.pca)

#cos2 des variables sur toutes les dimensions

#library("corrplot")
#corrplot(var$cos2, is.corr=FALSE)

ggcorrplot(var$cos2,
           method = "circle",
           outline.col = "white",
           ggtheme = ggplot2::theme_light(),
           lab = T,
           insig = "blank",
           colors = c("blue", "white", "red")
           ) +
  ggtitle("Matrice de corrélation") +
  theme(
    panel.background = element_rect(fill = "lightblue",
                                    colour = "lightblue",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "white"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "white"),
    plot.title = element_text(hjust = 0.5)
    )

# Cos2 total des variables sur Dim.1 et Dim.2
fviz_cos2(res.pca, choice = "var", axes = 1:2)

# Cos2 total des variables sur Dim.3 et Dim.4
fviz_cos2(res.pca, choice = "var", axes = 3:4)

# Cos2 total des variables sur les 4 1ères dimensions
fviz_cos2(res.pca, choice = "var", axes = 1:4)

# Avec Factoextra : colorer en fonction du cos2 = qualité de représentation
# La transparence détermine la qualité
fviz_pca_var(res.pca, alpha.var = "cos2",
             repel = TRUE,
             title = "cercle des corrélations"
             ) +
  theme(plot.title = element_text(hjust = 0.5))
# La couleur détermine la qualité
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("deeppink2", "lightseagreen", "deepskyblue4"),
             repel = TRUE,
             title = "Cercle des corrélations"
             ) +
  theme(plot.title = element_text(hjust = 0.5))
# Couleur et transparence determinent la qualité
fviz_pca_var(res.pca, col.var = "cos2", alpha.var = "cos2",
             gradient.cols = c("deeppink2", "lightseagreen", "deepskyblue4"),
             repel = TRUE,
             title = "cercle des corrélations"
             ) +
  theme(plot.title = element_text(hjust = 0.5))

# Contributions des variables aux axes principaux
#------------------------------------------------

corrplot(var$contrib, is.corr=FALSE)

ggcorrplot(var$contrib,
           method = "circle",
           outline.col = "white",
           ggtheme = ggplot2::theme_light(),
           lab = T,
           insig = "blank",
           colors = c("blue", "white", "red")
           ) +
  ggtitle("Matrice de corrélation") +
  theme(
    panel.background = element_rect(fill = "lightblue",
                                    colour = "lightblue",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    colour = "white"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    colour = "white"),
    plot.title = element_text(hjust = 0.5)
  )

#La ligne en pointillé rouge, sur les graphiques ci-dessous indique la contribution moyenne attendue. Si la contribution des variables était uniforme, la valeur attendue serait 1/length(variables) = 1/7
# Contributions des variables à la Composante Principale 1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 7)
# Contributions des variables à la Composante Principale 2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 7)
# Contribution totale aux Composantes Principales 1 et 2
fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 7)
# Contributions des variables à la Composante Principale 3
fviz_contrib(res.pca, choice = "var", axes = 3, top = 7)
# Contributions des variables à la Composante Principale 4
fviz_contrib(res.pca, choice = "var", axes = 4, top = 7)
# Contribution totale aux Composantes Principales 3 et 4
fviz_contrib(res.pca, choice = "var", axes = 3:4, top = 7)
# Contribution totale aux 4 1ères Composantes Principales
fviz_contrib(res.pca, choice = "var", axes = 1:4, top = 7)

# Mise en évidence des variables les plus importantes (ou contributives)
fviz_pca_var(res.pca, col.var = "contrib",
             alpha.var = "contrib",
             gradient.cols = c("deeppink2", "lightseagreen", "deepskyblue4"),
             repel = TRUE,
             title = "Cercle des corrélations"
             ) +
  theme(plot.title = element_text(hjust = 0.5)
        )


# Description des dimensions
#===========================

# Identifier les variables les plus corrélées avec une dimension donnée

res.desc <- dimdesc(res.pca, axes = c(1:4), proba = 0.05)
# Description de la dimension 1
res.desc$Dim.1
# Description de la dimension 2
res.desc$Dim.2
# Description de la dimension 3
res.desc$Dim.3
# Description de la dimension 4
res.desc$Dim.4


# Graphique des individus
#========================

ind <- get_pca_ind(res.pca)

# Les 20 individus les plus représentatifs de PC1 et PC2 sont labellisés
plot.PCA(res.pca, choix = "ind", autoLab = "yes", select = "contrib 20")

# Avec Factoextra

# Cas des Composantes 1 et 2
#---------------------------
# Coloré selon le cosinus carré
fviz_pca_ind(res.pca,
             axes = c(1,2),
             col.ind = "cos2",
             gradient.cols = c("deeppink2", "lightseagreen", "deepskyblue4"),
             label = "none",
             repel = TRUE # Évite le chevauchement de texte
             ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_set(theme_gray())

# Coloré selon la contribution
fviz_pca_ind(res.pca,
             axes = c(1,2),
             col.ind = "contrib",
             gradient.cols = c("deeppink2", "lightseagreen", "deepskyblue4"),
             label = "none",
             repel = TRUE # Évite le chevauchement de texte
             ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_set(theme_gray())

# Visualiser les 20 1ers individus représentatifs pour PC1 et PC2
fviz_pca_ind(res.pca,
             axes = c(1,2),
             col.ind = "contrib",
             gradient.cols = c("deeppink2", "lightseagreen", "deepskyblue4"),
             select.ind = list(contrib = 20),  
             repel = TRUE # Évite le chevauchement de texte
             ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_set(theme_gray())

# Taille et couleur des points en fonction du cos2
fviz_pca_ind(res.pca,
             axes = c(1,2),
             col.ind = "cos2",
             pointsize = "cos2",
             gradient.cols = c("deeppink2", "lightseagreen", "deepskyblue4"),
             label = "none",
             repel = TRUE
             ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_set(theme_gray())

# Cas des Composantes 3 et 4
#---------------------------
# Coloré selon le cosinus carré
fviz_pca_ind(res.pca,
             axes = c(3,4),
             col.ind = "cos2",
             gradient.cols = c("deeppink2", "lightseagreen", "deepskyblue4"),
             label = "none",
             repel = TRUE # Évite le chevauchement de texte
) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_set(theme_gray())

# Coloré selon la contribution
fviz_pca_ind(res.pca,
             axes = c(3,4),
             col.ind = "contrib",
             gradient.cols = c("deeppink2", "lightseagreen", "deepskyblue4"),
             label = "none",
             repel = TRUE # Évite le chevauchement de texte
) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_set(theme_gray())

# Visualiser les 20 1ers individus représentatifs pour PC1 et PC2
fviz_pca_ind(res.pca,
             axes = c(3,4),
             col.ind = "contrib",
             gradient.cols = c("deeppink2", "lightseagreen", "deepskyblue4"),
             select.ind = list(contrib = 20),  
             repel = TRUE # Évite le chevauchement de texte
) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_set(theme_gray())

# Taille et couleur des points en fonction du cos2
fviz_pca_ind(res.pca,
             axes = c(3,4),
             col.ind = "cos2",
             pointsize = "cos2",
             gradient.cols = c("deeppink2", "lightseagreen", "deepskyblue4"),
             label = "none",
             repel = TRUE
) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_set(theme_gray())


# Bar plot de la qualité de représentation (cos2) des individus
fviz_cos2(res.pca, choice = "ind", axes = 1) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, size = 6))
fviz_cos2(res.pca, choice = "ind", axes = 2) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, size = 6))
fviz_cos2(res.pca, choice = "ind", axes = 3) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, size = 6))
fviz_cos2(res.pca, choice = "ind", axes = 4) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, size = 6))


# Contribution des 20 premièrs individus sur PC1
fviz_contrib(res.pca,
             choice = "ind",
             axes = 1,
             top = 20)

# Contribution des 20 premièrs individus sur PC2
fviz_contrib(res.pca,
             choice = "ind",
             axes = 2,
             top = 20)

# Contribution des 20 premièrs individus sur PC3
fviz_contrib(res.pca,
             choice = "ind",
             axes = 3,
             top = 20)

# Contribution des 20 premièrs individus sur PC4
fviz_contrib(res.pca,
             choice = "ind",
             axes = 4,
             top = 20)

# Contribution des 20 premiers individus sur PC1 et PC2
fviz_contrib(res.pca, choice = "ind", axes = 1:2, top = 20)

# Contribution des 20 premiers individus sur PC3 et PC4
fviz_contrib(res.pca, choice = "ind", axes = 3:4, top = 20)

# # Contribution des 20 premiers individus sur PC1, PC2, PC3 et PC4
fviz_contrib(res.pca, choice = "ind", axes = 1:4, top = 20)


# Biplot des individus et des variables
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "deeppink2",
                col.ind = "deepskyblue3",
                label = "var"
                ) +
  ggtitle("Biplot des individus et des variables") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_set(theme_gray())

#==============================================================================#
#                                   3. Clustering                              #
#==============================================================================#

# 1 / Classification hierarchique ascendante (CAH)
#=================================================

# Calcul des distances & graphiques
#----------------------------------

num_tfruits.dist <- dist(num_tfruits, method = "euclidean")
num_fruits.hc.complete <- hclust(num_tfruits.dist, method = "complete")
num_fruits.hc.ward <- hclust(num_tfruits.dist, method = "ward.D")
plot(num_fruits.hc.complete, cex = 0.6)
plot(num_fruits.hc.ward, cex = 0.6)

# Avec les variables centrées & réduites (mais ça change rien)
num_tfruits.scale <- scale(num_tfruits, center = TRUE, scale = TRUE)
num_tfruits.scale.dist <- dist(num_tfruits.scale, method = "euclidean")
num_fruits.scale.hc.complete <- hclust(num_tfruits.dist, method = "complete")
num_fruits.scale.hc.ward <- hclust(num_tfruits.dist, method = "ward.D")
plot(num_fruits.scale.hc.complete, cex = 0.6)
rect.hclust(num_fruits.scale.hc.complete, k = 3)
plot(num_fruits.scale.hc.ward, cex = 0.6)
rect.hclust(num_fruits.scale.hc.ward, k = 3)

# Graphiques
#-----------


#HCPC(num_tfruits, graph = TRUE)

# Avec FactoMineR/Factoextra
#---------------------------

res.cah <- HCPC(num_tfruits, graph = FALSE)

# Dendogramme
fviz_dend(res.cah,
          cex = 0.7,
          palette = c("deeppink2", "lightseagreen", "deepskyblue4"),
          rect = TRUE, rect_fill = TRUE,
          rect_border = c("deeppink2", "lightseagreen", "deepskyblue4"),
          labels_track_height = 0.8
          )

# Factor map
fviz_cluster(res.cah,
             axes = c(1,2),
             repel = TRUE,
             show.clust.cent = TRUE,
             palette = c("deeppink2", "lightseagreen", "deepskyblue4"),
             ggtheme = theme_minimal(),
             main = "Factor map"
)

plot(res.cah, choice = "3D.map", ind.names = F)

fviz_cluster(res.cah,
             axes = c(3,4),
             repel = TRUE,
             show.clust.cent = TRUE,
             palette = c("deeppink2", "lightseagreen", "deepskyblue4"),
             ggtheme = theme_minimal(),
             main = "Factor map"
)

# Taille des clusters
ggplot(data = res.cah$data.clust, aes(x=clust)) +
  geom_bar(color="steelblue", fill="steelblue") +
  ggtitle("Taille des clusters")

# 2 / Classification Hiérarchique sur Composantes Principales (HCPC)
#===================================================================

res.hcpc <- HCPC(res.pca, graph = FALSE)

fviz_dend(res.hcpc, 
          cex = 0.7,
          palette = c("deeppink2", "lightseagreen", "deepskyblue4"),               # Palette de couleur ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE,
          rect_border = c("deeppink2", "lightseagreen", "deepskyblue4"),           # Couleur du rectangle
          labels_track_height = 0.8
)

fviz_cluster(res.hcpc,
             repel = TRUE,
             show.clust.cent = TRUE,
             palette = c("deeppink2", "lightseagreen", "deepskyblue4"),
             ggtheme = theme_minimal(),
             main = "Factor map"
)

# Taille des clusters
ggplot(data = res.hcpc$data.clust, aes(x=clust)) +
  geom_bar(color="steelblue", fill="steelblue") +
  ggtitle("Taille des clusters")

#plot.HCPC(res.hcpc, choice = "3D.map", ind.names = F)


# 3 / Partitionnement en K-moyennes (k-means)
#============================================

# Calcul du nombre optimal de clusters
#-------------------------------------

# Silhouette method
fviz_nbclust(num_tfruits, hcut, method = "silhouette") +
  labs(subtitle = "Silhouette method")

# Elbow method
fviz_nbclust(num_tfruits, hcut, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method")

# Gap statistic
fviz_nbclust(num_tfruits, hcut, method = "gap_stat", nboot = 500) +
  labs(subtitle = "Gap statistic method")

# Avec NbClust
library(NbClust)
NbClust(data = num_tfruits, distance = "euclidean",
        min.nc = 2, max.nc = 10, method = "complete")

NbClust(data = num_tfruits, distance = "euclidean",
        min.nc = 2, max.nc = 10, method = "ward.D")

NbClust <- NbClust(data = num_tfruits, distance = "euclidean",
        min.nc = 2, max.nc = 10, method = "kmeans")

NbClust_Best.nc <- as.data.frame(NbClust$Best.nc[1,])
ggplot(data=NbClust_Best.nc, aes(x=as.character(NbClust_Best.nc[,1]))) +
  geom_bar(stat="count", fill="steelblue") +
  labs(title="Fréquences des nombres de clusters calculés selon NbClust",
       x="Nombre de clusters k",
       y = "Fréquence parmi tous les indices") +
  theme(plot.title = element_text(hjust = 0.5))


# Affectation des kmeans
#-----------------------

num_fruits.kmeans <- kmeans(num_tfruits, centers = 3, nstart = 20)


#test <- PCA(num_tfruits, ncp = 3, graph = FALSE)
#fviz_dend(test, 
#          cex = 0.7,
#          palette = c("deeppink2", "lightseagreen", "deepskyblue4", "gold"),               # Palette de couleur ?ggpubr::ggpar
#          rect = TRUE, rect_fill = TRUE, # Rectangle autour des groupes
#          rect_border = c("deeppink2", "lightseagreen", "deepskyblue4", "gold"),           # Couleur du rectangle
#          labels_track_height = 0.8      # Augment l'espace pour le texte
#          )

# Visualiser la relation entre les clusters et chacune des variables utilisées
pairs(num_tfruits, col=c(1:3)[num_fruits.kmeans$cluster])

fviz_cluster(num_fruits.kmeans,
             axes = c(1,2),
             num_tfruits,
             ellipse.type = "norm"
)

fviz_cluster(num_fruits.kmeans,
             axes = c(3,4),
             num_tfruits,
             ellipse.type = "norm"
)

# Taille des clusters
num_fruits.kmeans$size


#==============================================================================#
#                              4. Modélisation                                 #
#==============================================================================#

library(pls)
library(MASS)


# Créer un nouveau le jeu de données
#===================================

tfruits2 <- read.table("teneurs_fruits.csv",
                      sep = ",",
                      dec = ",",
                      header = T
                      #row.names = T,
                      #check.names = F,
                      #fill = TRUE
)

#colClean <- function(x){ colnames(x) <- gsub("..", ".", colnames(x)); x }
colnames(tfruits2)[colnames(tfruits2) == 'X'] <- "especes"
colnames(tfruits2) <- gsub(x=colnames(tfruits2), pattern = "\\.$", replacement = "")
colnames(tfruits2) <- gsub(x=colnames(tfruits2), pattern = "\\.\\.", replacement = "\\.")
# Convertir les colonnes qui sont en "g" en "mg"
# Repérer les colonnes en "g" au lieu de "mg" (2 façons)
col_grammes <- grep(".g", names(tfruits2), value = T)

for (i in 1:length(col_grammes)) {
  # Multiplier par 1000 pour convertir en mg
  tfruits2[col_grammes[i]] <- get(col_grammes[i], tfruits2) * 1000
  # changer nom de la colonne
  colnames(tfruits2)[match(colnames(tfruits2[col_grammes[i]]), colnames(tfruits2))] <- sub(".g", ".mg", col_grammes[i])
}
# Supprimer les colonnes "joules.kJ" et Calories.kcal
tfruits2$joules.kJ <- NULL

YC <- tfruits2$Calories.kcal
XC <- tfruits
null=lm(YC ~ -1, data = data.frame(XC)) 
# définition du modèle complet
full=lm(YC ~ -1+., data = data.frame(XC)) 
# choix du meilleur modèle par SW en utilisant le critère d'AIC pour optimiser le modèle
fit.sw=stepAIC(null, scope=list(lower=null, upper=full), direction="forward",k=log(n),trace=0) 
pander(summary(fit.sw))

# PCR
fit.pcr <- pcr(YC~.,data=data.frame(XC), validation = "CV")
plot(fit.pcr$validation$PRESS[1:10], type= "b")
fit.pcr.2 <- pcr(YC~.,data=data.frame(XC), ncomp=2)
summary(fit.pcr)
yfit <- predict(fit.pcr,ncomp = 2)

# Plot Y-Yfit
plot(YC,yfit,pch=20,main="Y-Yfit")
