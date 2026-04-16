# ============================================================
# PROJET : Décomposition SVD et Pseudo-inverse
# Mathématiques pour le Big Data — ESGI 4ABD
# ============================================================

# ============================================================
# PARTIE 1 — Nos outils de calcul
# ============================================================

# Création de notre propre fonction SVD
ma_svd <- function(A) {
  # On commence par calculer la matrice de corrélation
  AtA <- t(A) %*% A
  eigen_res <- eigen(AtA)
  
  # On trie les résultats par ordre d'importance (valeurs propres décroissantes)
  ordre <- order(eigen_res$values, decreasing = TRUE)
  valeurs_propres <- eigen_res$values[ordre]
  V <- eigen_res$vectors[, ordre]
  
  # On identifie le rang en ignorant les valeurs proches de zéro
  r <- sum(valeurs_propres > 1e-10)
  sigma <- sqrt(valeurs_propres[1:r])
  
  # Construction de la matrice U à partir de A et V
  Vr <- V[, 1:r, drop = FALSE]
  U_r_brut <- A %*% Vr
  Ur <- apply(U_r_brut, 2, function(col) col / sqrt(sum(col^2)))
  
  # Préparation de la matrice Sigma (aux bonnes dimensions)
  m <- nrow(A); n <- ncol(A)
  Sigma <- matrix(0, nrow = m, ncol = n)
  for (i in 1:r) Sigma[i, i] <- sigma[i]
  
  # Si nécessaire, on complète Ur pour qu'elle soit une base orthonormée complète
  if (ncol(Ur) < m) {
    Ur <- qr.Q(qr(cbind(Ur, matrix(rnorm(m * (m - r)), m, m - r))))
  }
  
  return(list(U = Ur, Sigma = Sigma, V = V, sigma = sigma, r = r))
}

# Calcul de la pseudo-inverse de Moore-Penrose
pseudo_inverse <- function(A, tol = 1e-10) {
  res <- ma_svd(A)
  r   <- res$r
  Ur  <- res$U[, 1:r, drop = FALSE]
  Vr  <- res$V[, 1:r, drop = FALSE]
  
  # On inverse uniquement les valeurs singulières non nulles
  D_inv <- diag(1 / res$sigma[1:r], nrow = r)
  return(Vr %*% D_inv %*% t(Ur))
}

# Fonction pour résoudre et analyser un système linéaire
resoudre_systeme <- function(A, b) {
  cat("--- Analyse des dimensions ---\n")
  cat("La matrice A possède", nrow(A), "lignes et", ncol(A), "colonnes.\n")
  
  rang_A    <- qr(A)$rank
  A_aug     <- cbind(A, b)
  rang_aug  <- qr(A_aug)$rank
  
  cat("\nRésultats de l'analyse :\n")
  cat("- Rang de la matrice A :", rang_A, "\n")
  cat("- Rang de la matrice augmentée :", rang_aug, "\n")
  cat("- Nombre d'inconnues à trouver :", ncol(A), "\n")
  
  cat("\n--- Verdict du système ---\n")
  if (rang_A == rang_aug && rang_A == ncol(A)) {
    cat("C'est un cas idéal : une seule solution exacte existe.\n")
  } else if (rang_A == rang_aug && rang_A < ncol(A)) {
    cat("Il y a une infinité de solutions. On choisit celle avec la plus petite norme.\n")
  } else {
    cat("Le système n'a pas de solution exacte. On cherche le meilleur compromis possible.\n")
  }
  
  # Calcul de la solution via la pseudo-inverse
  A_plus <- pseudo_inverse(A)
  x_sol  <- A_plus %*% b
  residu <- sqrt(sum((A %*% x_sol - b)^2))
  
  cat("\n--- Solution trouvée (x) ---\n")
  print(round(x_sol, 2))
  
  cat("\nPetit check de vérification (A * x) :\n")
  print(round(A %*% x_sol, 2))
  
  cat("\nValeur cible (b) :\n")
  print(b)
  
  cat("\nÉcart constaté (résidu) :", round(residu, 6), "\n")
  if (residu < 1e-8) {
    cat("Le résidu est quasi nul : on a trouvé une solution exacte.\n")
  } else {
    cat("On a une solution approchée (optimale au sens des moindres carrés).\n")
  }
  
  return(invisible(list(x = x_sol, residu = residu, 
                        rang_A = rang_A, rang_aug = rang_aug)))
}

# ============================================================
# PARTIE 2 — Étude du premier cas (Section 3.2)
# ============================================================

A1 <- matrix(c(-18, 2,-14, -2,
               13,19, 11, 21,
               -4,-4,-12,  4,
               4,12,  8,  8), 
             nrow = 4, byrow = FALSE)

B1 <- matrix(c( 6, 2, 0,-1,
                -8, 7,-1,-2,
                -4,-5,-8, 4,
                5,-6, 2, 4,
                -4, 4, 2,-8), 
             nrow = 4, byrow = FALSE)

cat("\n****************************************\n")
cat("  ANALYSE SVD : MATRICE A1\n")
cat("****************************************\n")
svd_A1 <- ma_svd(A1)
cat("Valeurs singulières (sigma) :", round(svd_A1$sigma, 2), "\n")
cat("Rang détecté :", svd_A1$r, "\n")

# On vérifie si on arrive à reconstruire la matrice de base
reconstitution_A1 <- svd_A1$U %*% svd_A1$Sigma %*% t(svd_A1$V)
cat("Précision de la reconstruction (écart max) :", 
    round(max(abs(reconstitution_A1 - A1)), 8), "\n")

# Test rapide avec la fonction standard de R pour comparer
svd_native_A1 <- svd(A1)
cat("Valeurs singulières de référence (R natif) :", round(svd_native_A1$d, 2), "\n")

cat("\n****************************************\n")
cat("  ANALYSE SVD : MATRICE B1\n")
cat("****************************************\n")
svd_B1 <- ma_svd(B1)
cat("Valeurs singulières (sigma) :", round(svd_B1$sigma, 2), "\n")
cat("Rang détecté :", svd_B1$r, "\n")

reconstitution_B1 <- svd_B1$U %*% svd_B1$Sigma %*% t(svd_B1$V)
cat("Précision de la reconstruction (écart max) :", 
    round(max(abs(reconstitution_B1 - B1)), 8), "\n")

svd_native_B1 <- svd(B1)
cat("Valeurs singulières de référence (R natif) :", round(svd_native_B1$d, 2), "\n")

# ============================================================
# PARTIE 3 — Étude du second cas (Section 3.3)
# ============================================================

A2 <- matrix(c(-3,-1, 0, 0,
               -3,-1, 0, 0,
               -6,-1,-1, 1,
               6, 1, 1, 1,
               1,-2,-1,-1), 
             nrow = 4, byrow = FALSE)

B2 <- matrix(c( 4,-5, 2, 6,
                0, 0, 0, 0,
                -1, 3,-1,-3,
                -2, 5,-2,-6,
                0, 0, 0, 0), 
             nrow = 4, byrow = FALSE)

b <- c(6, -1, -4, 6)

cat("\n****************************************\n")
cat("  CALCUL DES PSEUDO-INVERSES\n")
cat("****************************************\n")

A2_plus <- pseudo_inverse(A2)
cat("Voici la pseudo-inverse de A2 :\n")
print(round(A2_plus, 2))

B2_plus <- pseudo_inverse(B2)
cat("\nEt celle de B2 :\n")
print(round(B2_plus, 2))

cat("\n****************************************\n")
cat("  RÉSOLUTION DU SYSTÈME A2 * x = b\n")
cat("****************************************\n")
res_A <- resoudre_systeme(A2, b)

cat("\n****************************************\n")
cat("  RÉSOLUTION DU SYSTÈME B2 * x = b\n")
cat("****************************************\n")
res_B <- resoudre_systeme(B2, b)