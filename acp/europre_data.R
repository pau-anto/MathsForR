# ============================================================
# Analyse en Composantes Principales (ACP) - Données Europe
# ============================================================

# ── Librairies ───────────────────────────────────────────────
library(ggplot2)
library(readr)
library(corrplot)
library(FactoMineR)
library(ggrepel)


# ── Chargement des données ───────────────────────────────────
europe <- read.table(
  "C:/Users/Julien ANTOGNELLI/Documents/ESGI/Maths/MathsForR/data/Donnees Europe.txt",
  header = TRUE,   # 1ère ligne = noms des colonnes
  row.names = 1    # 1ère colonne = noms des lignes (at, be, cy...)
)

# Conversion requise par FactoMineR
europe <- as.data.frame(europe)


# ============================================================
# PREMIÈRE ACP
# ============================================================

# ── Introduction au jeu de données ──────────────────────────
head(europe)


# ── Graphes de couples de variables ─────────────────────────

# Couple chom / depsoc
ggplot(europe, aes(x = chom, y = depsoc, label = rownames(europe))) +
  geom_point() +
  geom_text_repel(size = 3) +
  labs(x = "chom", y = "depsoc")

# Couple pib / devel
ggplot(europe, aes(x = pib, y = devel, label = rownames(europe))) +
  geom_point() +
  geom_text_repel(size = 3) +
  labs(x = "PIB", y = "devel")

# Couple benev / pvap
ggplot(europe, aes(x = benev, y = pvap, label = rownames(europe))) +
  geom_point() +
  geom_text_repel(size = 3) +
  labs(x = "benev", y = "pvap")

# Interprétation :
# (chom, depsoc) -> nuage dispersé, pas de tendance
# (PIB, devel)   -> tendance croissante nette, mais un point isolé avec un très fort PIB
# (benev, pvap)  -> tendance négative modérée


# ── Matrice des corrélations ─────────────────────────────────

# Matrice brute
cor_mat <- round(cor(europe, method = "pearson"), 2)
cor_mat
write.csv(cor_mat, file = "./graphs_&_tables/corr_matrix.csv", row.names = TRUE)

# Visualisation
corrplot(cor_mat,
         method      = "color",
         type        = "full",
         diag        = TRUE,
         col         = colorRampPalette(c("#FFD54F", "#FFFDE7", "#B71C1C"))(200),
         tl.col      = "black",
         tl.cex      = 0.8,
         tl.srt      = 35,
         addCoef.col = "black",
         number.cex  = 0.65,
         mar         = c(1, 1, 2, 1))


# ── ACP normée ───────────────────────────────────────────────
# scale.unit = TRUE  : centre et réduit les données (ACP normée)
# ncp = 4            : on conserve les 4 premières composantes principales
res.pca <- PCA(europe, scale.unit = TRUE, ncp = 4, graph = FALSE)


# ── Variances des composantes principales ───────────────────
variances_cp <- data.frame(
  Axe      = paste0("Dim.", 1:4),
  Variance = round(res.pca$eig[1:4, 1], 2)
)
variances_cp
write.csv(variances_cp, file = "./graphs_&_tables/variances_cp.csv", row.names = FALSE)


# ── Tableau de corrélation des variables ─────────────────────
corr_variables <- round(res.pca$var$cor[, 1:4], 2)
corr_variables
write.csv(corr_variables, file = "./graphs_&_tables/corr_variables.csv", row.names = TRUE)

# Coordonnées des variables = corrélations variable / composante
# Toutes les variables ont une forte corrélation avec Dim1 -> effet de taille


# ── Coordonnées des individus sur les axes ───────────────────
coord_individus <- round(res.pca$ind$coord[, 1:4], 2)
coord_individus
write.csv(coord_individus, file = "./graphs_&_tables/coord_individus.csv", row.names = TRUE)


# ── Contribution des individus aux axes (en %) ───────────────
contrib_individus <- round(res.pca$ind$contrib[, 1:4], 2)
contrib_individus
write.csv(contrib_individus, file = "./graphs_&_tables/contrib_individus.csv", row.names = TRUE)


# ── Graphe des valeurs propres (scree plot) ──────────────────
png("./graphs_&_tables/valeurs_propres.png", width = 800, height = 600, res = 120)

vp <- res.pca$eig[, 1]

plot(1:length(vp), vp, type = "b", pch = 19,
     xlab = "Axes principaux",
     ylab = "Valeurs propres")

abline(h = 1, col = "red", lty = 2)   # Critère de Kaiser

legend("topright",
       legend = expression("Critère de Kaiser (" * lambda == 1 * ")"),
       col    = "red",
       lty    = 2,
       lwd    = 1,
       bty    = "n")

dev.off()

# Selon la règle de Kaiser : 2 composantes à retenir (Dim1 ~8.5 et Dim2 ~1.2)


# ── Table complète valeurs propres / % variance / cumul ──────
res_filtre <- round(res.pca$eig, 2)
res_filtre

# Filtrage sur valeur propre > 1
res_filtre <- round(res.pca$eig[res.pca$eig[, 1] > 1, ], 2)
res_filtre
# 2 axes retenus, expliquant 73.91 % de l'inertie


# ── Cercle de corrélation des variables ──────────────────────
png("./graphs_&_tables/cercle_correlation.png", width = 800, height = 800, res = 120)
plot(res.pca, choix = "var", axes = c(1, 2), title = " ")
dev.off()

plot(res.pca, choix = "var", axes = c(1, 2), title = "")


# ── Tableau des coordonnées des variables sur le cercle ──────
corr_vars <- round(res.pca$var$coord[, 1:2], 2)

axe_dominant <- apply(corr_vars, 1, function(x) {
  if (abs(x[1]) > abs(x[2])) {
    paste0("Dim 1 (", ifelse(x[1] >= 0, "+", "-"), ")")
  } else if (abs(x[2]) > abs(x[1])) {
    paste0("Dim 2 (", ifelse(x[2] >= 0, "+", "-"), ")")
  } else {
    "Dim 1 et 2"
  }
})

tab_axes <- data.frame(
  Variable      = rownames(corr_vars),
  `Axe dominant` = axe_dominant,
  `r (Dim 1)`   = corr_vars[, 1],
  `r (Dim 2)`   = corr_vars[, 2],
  row.names     = NULL,
  check.names   = FALSE
)

tab_axes
write.csv(tab_axes, file = "./graphs_&_tables/axes_variables.csv", row.names = FALSE)

# Toutes les variables sont dominées par Dim1 -> effet de taille important


# ============================================================
# ACP SANS `lu`
# ============================================================

# ── Retrait de lu et nouvelle ACP ────────────────────────────
no_lu <- europe[rownames(europe) != "lu", ]

res.pca_lu <- PCA(no_lu, scale.unit = TRUE, ncp = 4, graph = FALSE)


# ── Variances des CP (sans lu) ───────────────────────────────
variances_cp_lu <- data.frame(
  Axe      = paste0("Dim.", 1:4),
  Variance = round(res.pca_lu$eig[1:4, 1], 2)
)
variances_cp_lu
write.csv(variances_cp_lu, file = "./graphs_&_tables/variances_cp_lu.csv", row.names = FALSE)


# ── Corrélations des variables (sans lu) ─────────────────────
corr_variables_lu <- round(res.pca_lu$var$cor[, 1:4], 2)
corr_variables_lu
write.csv(corr_variables_lu, file = "./graphs_&_tables/corr_variables_lu.csv", row.names = TRUE)


# ── Coordonnées de `lu` en individu supplémentaire ───────────
# predict() projette lu dans l'espace de l'ACP construite sans lui
lu_sup   <- predict(res.pca_lu, newdata = europe["lu", ])
lu_coord <- as.data.frame(round(lu_sup$coord[, 1:4], 2))
lu_coord


# ── Coordonnées des individus actifs (sans lu) ───────────────
coord <- round(res.pca_lu$ind$coord[, 1:4], 2)
coord
write.csv(coord, file = "./graphs_&_tables/coord_no_lu.csv", row.names = TRUE)


# ── Contributions aux axes (en %) ────────────────────────────
contrib <- as.data.frame(round(res.pca_lu$ind$contrib[, 1:4], 2))
contrib
write.csv(contrib, file = "./graphs_&_tables/contrib_no_lu.csv", row.names = TRUE)


# ── Qualités de représentation (cos²) ────────────────────────
cos2 <- as.data.frame(round(res.pca_lu$ind$cos2[, 1:4] * 100, 2))
cos2
write.csv(cos2, file = "./graphs_&_tables/cos2_no_lu.csv", row.names = TRUE)


# ── Scree plot sans lu ───────────────────────────────────────
png("./graphs_&_tables/scree_plot_lu.png", width = 800, height = 600, res = 120)

vp_lu <- res.pca_lu$eig[, 1]

plot(1:length(vp_lu), vp_lu, type = "b", pch = 19,
     xlab = "Axes principaux",
     ylab = "Valeurs propres")

abline(h = 1, col = "red", lty = 2)

legend("topright",
       legend = expression("Critère de Kaiser (" * lambda == 1 * ")"),
       col    = "red",
       lty    = 2,
       lwd    = 1,
       bty    = "n")

dev.off()


# ── Table complète valeurs propres (sans lu) ─────────────────
res_filtre_lu <- round(res.pca_lu$eig, 2)
res_filtre_lu

# Conclusion :
# Avec lu    : (3.36 + 1.46 + 1.25) / 8 ≈ 75.9 % de variance expliquée sur 3 axes
# Sans lu    : (4.56 + 1.40 + 1.13) / 9 ≈ 78.8 % de variance expliquée sur 3 axes
# → La qualité globale s'améliore légèrement après retrait de lu.
#   Dim.1 passe de 3.36 à 4.56 : lu perturbait le premier axe via son PIB très atypique.
#   Le nombre de composantes à retenir reste identique (3).