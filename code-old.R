# Partie 1 : Impact de la structure de départ des simulations 
#            sur la distance d50 et la forme de la PR2


# 1. Représenter la distribution de la distance d50 pour 
#    les structures extraites des deux simulations. 
#    Commenter les résultats obtenus.


# Lire le fichier csv, et en extraire les données correspondant à chaque conformation
d50 <- read.csv("d50_3ebz_1hsi.csv", sep=' ',dec='.')
d50 <- head(d50, nrow(d50) - 6)

d50_1hsi <- d50[grep("1hsi", d50$structure), ]
d50_3ebz <- d50[grep("3ebz", d50$structure), ]

# Représenter les d50 de la structure 1HSI
hist(d50_1hsi$d50, main = "Distribution of d50 values for 1hsi",
     xlab = "d50 (Å) ", breaks = 20, col = "cadetblue")

# Représenter les d50 de la structure 3EBZ
hist(d50_3ebz$d50, main = "Distribution of d50 values for 3ebz",
     xlab = "d50 (Å)", breaks = 20, col = "coral")

# Représenter les deux sur le même histogramme avec des couleurs transparentes pouvant se chevaucher
hist(d50_1hsi$d50, main = "Distribution of d50 values for 1hsi and 3ebz",
     xlab = "d50 (Å) ", breaks = 20, col = adjustcolor( "cadetblue", alpha.f = 0.55))
hist(d50_3ebz$d50, breaks = 20, col = adjustcolor("coral", alpha.f = 0.55), add = T)
legend("topright", c("1HSI", "3EBZ"), fill=c(adjustcolor( "cadetblue", alpha.f = 0.55),
                                             adjustcolor( "coral", alpha.f = 0.55)))

# 2. Comparaison de la distance d50 moyenne dans les deux simulations. 
#    Quel est l’impact de la structure de départ sur la distance d50 de la PR2 ?

# a. la variable aléatoire : la distance d50 (quantitative comtinue)
# b. les paramètres dans la population et dans les échantillons :
#    Les paramètres calculés ci-dessous à partir des échantillons 
#    sont des estimateurs des paramètres de la population

mean_d50_1hsi <- mean(d50_1hsi$d50)
mean_d50_3ebz <- mean(d50_3ebz$d50)

median_d50_1hsi <- median(d50_1hsi$d50)
median_d50_3ebz <- median(d50_3ebz$d50)

sd_d50_1hsi <- sd(d50_1hsi$d50)
sd_d50_3ebz <- sd(d50_3ebz$d50)

var_d50_1hsi <- var(d50_1hsi$d50)
var_d50_3ebz <- var(d50_3ebz$d50)

# c. le test que nous allons réaliser : un test de comparaison 
#    entre deux moyennes, soit le test t pour les échantillons indépendants 
#     - H0: les moyennes sont égalles -  µ1 = µ2
#     - H1: les moyennes sont differents - µ1 != µ2
# d. les conditions de validité du test réalisé :
#     - les deux échantillons - 1hsi et 3ebz - sont indépendants (supposons)
#     - la taille des deux échantillons leur est suffisamment grande (1001)
#     - la variance des deux groupes est égale (on va verifier avec le test
#       de Fligner-Killeen car notres echantillons ne sont pas vraiment 
#       distribuées normalement)

fligner.test(d50_1hsi$d50, d50_3ebz$d50, sigma.x=sd_d50_1hsi, y, sigma.y=sd_d50_3ebz)

#       le test de Fligner-Killeen a retourné une p-value de 0.4684 > 0.05,
#       donc au niveau de signification de 0,05 il n'y a pas suffisamment de 
#       preuves pour rejeter H0. Par conséquent, on peut supposer que les 
#       variances des deux groupes sont égales.
# e. d’énoncer la règle de décision des tests réalisés

t.test(d50_1hsi$d50, d50_3ebz$d50, var.equal = TRUE)

# f. de conclure à la question biologique

#    Les résultats montrent que la différence moyenne de la distance d50 entre 
#    les deux structures est statistiquement significative au niveau de 0.05, 
#    avec une p-valeur < 2.2e-16. On peut donc rejeter H0 et conclure que les 
#    moyennes sont differents. L'intervalle de confiance à 95 % pour la 
#    différence moyenne se situe entre -2.463885 et -2.338162. Les moyennes des 
#    échantillons "1hsi" et "3ebz" sont respectivement de 5.045917 et 7.446941.

#    Cela indique que la structure de départ a un impact significatif sur la 
#    distance d50 de la PR2, étant donné que la distance d50 moyenne est 
#    significativement différente entre les deux structures.


# 3. Etudier et comparer l’évolution de la distance d50 au cours du temps pour 
#    les deux simulations.


# Ajout du temps dans les data frames
library(stringr)
time <- as.numeric(str_match(str_match(d50_1hsi$structure, "[0-9]{1,3}.pdb"), "[0-9]+"))
d50_1hsi <- cbind(d50_1hsi, time)

time <- as.numeric(str_match(str_match(d50_3ebz$structure, "[0-9]{1,3}.pdb"), "[0-9]+"))
d50_3ebz <- cbind(d50_3ebz, time)

# Représentation des données
library(ggplot2)
library(gridExtra)

# Représentation des deux données sur le même graphique
ggplot() +
  geom_line(data = d50_1hsi_summary, alpha = 0.4, aes(x = time, y = mean_d50, color = "cadetblue"), size = 0.01) +
  geom_line(data = d50_3ebz_summary, alpha = 0.4, aes(x = time, y = mean_d50, color = "coral"), size = 0.01) +
  geom_point(data=d50_1hsi, aes(time, d50), color = "cadetblue", size=1) +
  geom_smooth(data=d50_1hsi, aes(time, d50,, color = "1hsi")) +
  geom_point(data=d50_3ebz, aes(time, d50), color = "coral", size=1) +
  geom_smooth(data=d50_3ebz, aes(time, d50, color = "3ebz")) +
  xlab("Time (ns)") + ylab("D50 (Å)") +
  ggtitle("Evolution of D50 in 1hsi and 3ebz structure over time") +
  scale_color_manual(values = c("1hsi" = "darkblue", "3ebz" = "red"))


# Partie 2 : Impact de la conformation de départ sur la structure de la PR2


## 4. La première étape va consister à supprimer les distances dSL−Cα qui 
##    ont une variance nulle dans les 2002 structures. Pourquoi supprime-t-on 
##    ces distances ? Combien de distances dSL−Cα avez-vous supprimées ?

SLCalpha_1hsi <- read.csv("dist_SLCalpha_superlig_1HSI.csv", sep=' ',dec='.')
SLCalpha_3ebz <- read.csv("dist_SLCalpha_superlig_3EBZ.csv", sep=' ',dec='.')
SLCalpha <- rbind(SLCalpha_1hsi, SLCalpha_3ebz)
SLCalpha <- data.frame(SLCalpha, row.names = NULL)

library(reshape2)
SLCalpha <- dcast(SLCalpha, structure ~ CA, value.var = "dSLCA")
SLCalpha_var <- apply(SLCalpha, 2, var)
SLCalpha_var

##    Si la distance a une variance nulle, cela signifie que la distance est  
##    constante t a la même valeur quel que soit le type de structures ou 
##    la forme de la structure. Par conséquent, ces distances sont supprimées 
##    car ils n'apportent aucune information.

##    On a supprimé 0 distances dSL−Cα.

## 5. Réaliser une analyse en composante principale des 2002 structures à 
##    partir de distances dSL−Cα avec une variance non nulle. Analyser les  
##    résultats obtenus.

SLCalpha_cor_matrix <- cor(sapply(SLCalpha[-1,-1], as.numeric))
pca <- prcomp(SLCalpha_cor_matrix)

high_corr_pairs <- which(SLCalpha_cor_matrix > 0.85 & SLCalpha_cor_matrix < 1, arr.ind = TRUE)

colors <- rep("cadetblue", ncol(SLCalpha)-1)
colors[high_corr_pairs[,1]] <- "coral"
biplot(pca, cex = 0.02, arrow.len = 0.01, col = colors, expand = 1)
segments(0, 0, pca$rotation[high_corr_pairs[,1],1], pca$rotation[high_corr_pairs[,1],2], col = "coral")
  

library(factoextra)
library(FactoMineR)
pca <- PCA(SLCalpha_cor_matrix)
fviz_pca_var(pca, col.var = "contrib", repel = TRUE)
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 50))


## 6. Calculer et représenter la matrice de corrélation des distances dSL−Cα 
##    calculées avec les 2002 structures. Pour la représentation, vous pouvez 
##    utiliser la fonction corrplot() de la librairie corrplot().Commenter les 
##    résultats obtenus.

install.packages("corrplot")
library(corrplot)

corrplot(SLCalpha_cor_matrix, method = "color", type = "upper",  tl.cex = 0.2, tl.col = "black")

## 7. Supprimer les distances dSL−Cα ayant un coefficient de corrélation 
##    supérieur à 0.85 en utilisant la fonction findCorrelation(cutoff = 0.85) 
##    du package caret. Pourquoi supprimer ces distances dSL−Cα ? Combien de 
##    distances dSL−Cα avez-vous supprimées ?
  
install.packages("caret")
library(caret)

high_corr_indices <- findCorrelation(SLCalpha_cor_matrix, cutoff = 0.85)
low_corr_dist <- SLCalpha[, -high_corr_indices]

cat("Number of distances removed:", ncol(SLCalpha) - ncol(low_corr_dist))

## 8. A l’aide de PyMOL, réaliser une figure qui localise sur la structure les 
##    carbones Cα impliqués dans les distances dSL−Cα sélectionnées. Commenter 
##    la localisation de ces résidus sur la structure de la PR2.

colnames(SLCalpha[,-high_corr_indices])

## 9. A partir des distances dSL−Cα sélectionnées, calculer une nouvelle ACP des 
##    2002 structures. Commenter les résultats obtenus.

SLCalpha_low_cor_matrix <- cor(sapply(low_corr_dist[-1,-1], as.numeric))
pca_2 <- prcomp(SLCalpha_low_cor_matrix)
biplot(pca_2, cex = 0.2, arrow.len = 0.02, col = "cadetblue", expand = 1)

library(factoextra)
library(FactoMineR)
pca_2 <- PCA(SLCalpha_low_cor_matrix)
fviz_pca_var(pca_2, col.var = "contrib", repel = TRUE)
fviz_eig(pca_2, addlabels = TRUE, ylim = c(0, 50))

## 10. A partir des distances dSL−Cα sélectionnées, calculer une classification 
##     hiérarchique des 2002 structures. Commenter la classification obtenue.

hc <- hclust(dist(low_corr_dist), method = "ward.D2")
plot(hc, hang = -1, cex = 0.001)

## 11. A partir de la classification, extraire 6 groupes en utilisant la 
##     fonction cutree(). Commenter les différents groupes obtenus. Etudier et 
##     commenter la distribution des deux types de structures (structures 
##     extraites de la simulation lancée avec 1HSI ou avec 3EBZ) dans les 
##     différents groupes.

groups <- cutree(hc, k = 6)
low_corr_dist_with_groups <- data.frame(low_corr_dist, groups)
str <- rep(c("1hsi", "3bzez"), each = 1001)
low_corr_dist_with_groups <- data.frame(low_corr_dist_with_groups, str)
table(groups)

table_df <- table(low_corr_dist_with_groups$groups, low_corr_dist_with_groups$str)
my_colors <- c("cadetblue", "coral")
ggplot(data = as.data.frame(table_df), aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(x = "Group", y = "Count", fill = "Type") +
  ggtitle("Distribution of Types between Six Groups")
