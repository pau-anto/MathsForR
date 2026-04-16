# ============================================================
# PROJET T1 - PARTIE C : ACP SANS LUXEMBOURG
# ESGI 4ABD - Mathématiques pour le Big Data
# ============================================================

library(FactoMineR)

# --- Chargement des données ---
data <- read.table("Donnees Europe V2.txt", header = TRUE)
rownames(data) <- c("at","be","cy","cz","de","dk","ee","es","fi","fr",
                    "gr","hu","ie","it","lt","lu","lv","mt","nl","pl",
                    "pt","se","si","sk","uk")

# Séparation de lu (individu supplémentaire)
lu_data    <- data["lu", ]
data_actif <- data[rownames(data) != "lu", ]
n <- nrow(data_actif)
p <- ncol(data_actif)



# C.1 : CENTRAGE-RÉDUCTION ET TABLEAUX

# ACP normée sur les individus actifs (lu mis en supplémentaire)
# On remet lu temporairement dans le jeu pour qu'il soit projeté
data_avec_lu <- rbind(data_actif, lu_data)
lu_index     <- nrow(data_avec_lu)  # dernière ligne = lu

res.pca <- PCA(data_avec_lu, scale.unit = TRUE, ncp = 4,
               ind.sup = lu_index, graph = FALSE)

lambda      <- res.pca$eig[, 1]       # valeurs propres
pays_actifs <- rownames(data_actif)

cat("=== Variances des 4 premières composantes principales ===\n")
variances_cp <- data.frame(
  Axe          = paste0("F", 1:4),
  Valeur_propre = round(lambda[1:4], 4),
  Part_inertie  = round(res.pca$eig[1:4, 2], 2),
  Inertie_cum   = round(res.pca$eig[1:4, 3], 2)
)
print(variances_cp)

cat("\n=== Tableau des corrélations variables / composantes principales ===\n")
corr_var <- round(res.pca$var$cor[, 1:4], 4)
print(corr_var)

cat("\n=== Coordonnées des individus actifs (4 premiers axes) ===\n")
coords_ind <- round(res.pca$ind$coord[, 1:4], 4)
print(coords_ind)

cat("\n=== Coordonnées de lu (individu supplémentaire) ===\n")
print(round(res.pca$ind.sup$coord[, 1:4], 4))

cat("\n=== Contributions des individus aux axes (%) — seuil = 1/n =",
    round(100/n, 2), "% ===\n")
contrib <- round(res.pca$ind$contrib[, 1:4], 2)
print(contrib)

cat("\n=== Qualités de représentation cos² * 100 ===\n")
cos2 <- round(100 * res.pca$ind$cos2[, 1:4], 2)
print(cos2)





# C.2 : NOMBRE DE COMPOSANTES ET QUALITÉ GLOBALE

n_retain <- sum(lambda > 1)
cat("\n=== C.2 : Critère de Kaiser ===\n")
cat("Axes retenus (λ > 1) :", n_retain, "\n")
cat("Inertie expliquée    :", round(res.pca$eig[n_retain, 3], 2), "%\n")

png("eboulis_sans_lu.png", width = 800, height = 500)
barplot(lambda, names.arg = paste0("F", seq_along(lambda)),
        col  = ifelse(lambda > 1, "steelblue", "lightgray"),
        main = "Éboulis des valeurs propres (sans lu)",
        xlab = "Composante principale", ylab = "Valeur propre")
abline(h = 1, col = "red", lty = 2, lwd = 2)
legend("topright", c("λ > 1 (retenu)", "λ ≤ 1", "Seuil Kaiser"),
       fill = c("steelblue","lightgray", NA),
       lty  = c(NA, NA, 2), col = c(NA, NA, "red"), lwd = c(NA, NA, 2))
dev.off()





# C.3 : CERCLE DE CORRÉLATION

png("cercle_corr_sans_lu.png", width = 700, height = 700)
plot(res.pca, choix = "var", axes = c(1, 2),
     title = "Cercle de corrélation F1-F2 (sans lu)")
dev.off()





# C.4 : VARIABLES QUI DÉTERMINENT LES AXES RETENUS
# Critère : |r(variable, axe)| > 1/sqrt(p) ≈ 0.35 
# ou contribution > 100/p %

cat("\n=== C.4 : Variables déterminantes par axe (|r| > 1/sqrt(p)) ===\n")
seuil_r <- round(1 / sqrt(p), 3)
cat("Seuil retenu : |r| >", seuil_r, "\n\n")

for (k in 1:n_retain) {
  r_k  <- res.pca$var$cor[, k]
  vars_det <- names(r_k)[abs(r_k) > seuil_r]
  cat(sprintf("F%d : %s\n", k, paste(vars_det, collapse = ", ")))
  cat(sprintf("  (corrélations : %s)\n\n",
              paste(sprintf("%s=%.2f", vars_det, r_k[vars_det]), collapse = ", ")))
}

# Tableau synthétique des corrélations |r| > seuil pour les axes retenus
tab_det <- data.frame(
  Variable = rownames(corr_var)
)
for (k in 1:n_retain) {
  tab_det[[paste0("r_F", k)]] <- round(res.pca$var$cor[, k], 2)
}
cat("Tableau complet corrélations (axes retenus) :\n")
print(tab_det)





# C.5 : PLAN FACTORIEL F1-F2

new_2004 <- c("cy","cz","ee","hu","lt","lv","mt","pl","si","sk")
pays_col <- ifelse(pays_actifs %in% new_2004, "red", "steelblue")

pct <- round(res.pca$eig[1:3, 2], 1)

png("plan_factoriel_12_sans_lu.png", width = 900, height = 700)
plot(res.pca, choix = "ind", axes = c(1, 2),
     title = "Plan factoriel F1-F2 (sans lu)\n* = individu supplémentaire")
dev.off()





# C.6 : INTERPRÉTATION DES AXES RETENUS

cat("\n=== C.6 : Interprétation des axes retenus ===\n")

# Rappel des corrélations utiles (|r| > seuil)
corr_axes <- res.pca$var$cor[, 1:n_retain]

cat("\n--- Axe F1 ---\n")
cat("Variables positivement corrélées (pays développés) :\n")
cat("  depsoc (+0.70), devel (+0.83), pib (+0.65)\n")
cat("  => pays avec dépenses sociales élevées, fort IDH, PIB élevé\n")
cat("Variables négativement corrélées (pays moins développés) :\n")
cat("  chom (-0.70), pvap (-0.69), trvpv (-0.68)\n")
cat("  => pays avec fort chômage et forte pauvreté\n")
cat("Interprétation : F1 est un AXE DE DÉVELOPPEMENT SOCIO-ÉCONOMIQUE.\n")
cat("  Côté positif : anciens membres UE (dk, se, fi, at, be)\n")
cat("  Côté négatif : nouveaux membres 2004 (lt, lv, pl, sk, ee)\n\n")

cat("--- Axe F2 ---\n")
r2 <- sort(corr_axes[, 2], decreasing = TRUE)
cat("Variables positivement corrélées :\n")
cat("  depedu (+0.74), pvav (+0.67)\n")
cat("  => pays avec fortes dépenses d'éducation ET fort taux de pauvreté avant transferts\n")
cat("Variables négativement corrélées :\n")
cat("  pib (-0.56)\n")
cat("  => PIB faible relativement\n")
cat("Interprétation : F2 oppose des pays avec un système de REDISTRIBUTION FORT\n")
cat("  (dk, se : dépenses éducation + pvav élevé mais pvap réduit après transferts)\n")
cat("  aux pays à PIB élevé mais dépenses éducation modestes (nl, at).\n\n")

if (n_retain >= 3) {
  cat("--- Axe F3 ---\n")
  cat("Variables positivement corrélées :\n")
  cat("  pvap (+0.61), pvav (+0.61), trvpv (+0.50)\n")
  cat("  => pauvreté résiduelle élevée y compris parmi les travailleurs\n")
  cat("Variables négativement corrélées :\n")
  cat("  (aucune très forte au seuil retenu)\n")
  cat("Interprétation : F3 capture une PAUVRETÉ STRUCTURELLE non expliquée par F1.\n")
  cat("  ie (Irlande) se distingue avec fort pvap malgré PIB correct.\n\n")
}





# C.7 : INDIVIDUS LES MOINS BIEN REPRÉSENTÉS

cos2_actif <- res.pca$ind$cos2
cos2_cum   <- 100 * rowSums(cos2_actif[, 1:n_retain])

cat("=== C.7 : Individus les moins bien représentés (plan F1-F", n_retain,") ===\n")
ord <- order(cos2_cum)
df_mal_rep <- data.frame(
  pays       = pays_actifs[ord[1:5]],
  cos2_F1    = round(100 * cos2_actif[ord[1:5], 1], 1),
  cos2_F2    = round(100 * cos2_actif[ord[1:5], 2], 1),
  cos2_F3    = round(100 * cos2_actif[ord[1:5], 3], 1),
  cum_axes   = round(cos2_cum[ord[1:5]], 1)
)
print(df_mal_rep)




# C.8 : POSITION DE LU (individu supplémentaire)

lu_coords <- res.pca$ind.sup$coord[1, ]
cat("\n=== C.8 : Coordonnées de lu ===\n")
cat("F1:", round(lu_coords[1], 3),
    " F2:", round(lu_coords[2], 3),
    " F3:", round(lu_coords[3], 3), "\n")
cat("lu est à l'extrême gauche sur F2 (très négatif) :\n")
cat("  PIB très élevé (271.8) et depsoc faibles => atypique sur cet axe.\n")
cat("  Sur F1, lu est du côté positif (pays développé) mais pas extrême.\n")
cat("  Sa position excentrée sur F2 explique pourquoi il perturbait l'ACP initiale.\n")

