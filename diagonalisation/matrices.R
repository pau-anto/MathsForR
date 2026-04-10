

#SECTION 3.1 : Implémentation manuelle de la SVD

ma_svd <- function(A) {
  # Étape 1 : calculer AᵀA
  AtA <- t(A) %*% A
  
  # Étape 2 : valeurs propres et vecteurs propres de AᵀA
  eigen_res <- eigen(AtA)
  
  # Les valeurs singulières = racines carrées des valeurs propres
  ordre <- order(eigen_res$values, decreasing = TRUE)
  valeurs_propres <- eigen_res$values[ordre]
  V <- eigen_res$vectors[, ordre]  # matrice V
  
  # On garde seulement les valeurs propres > 0 (numériquement)
  r <- sum(valeurs_propres > 1e-10)
  sigma <- sqrt(valeurs_propres[1:r])
  
  # Étape 3 : construire U à partir des r premiers vecteurs de V
  Vr <- V[, 1:r, drop = FALSE]
  U_r_brut <- A %*% Vr
  
  # On normalise chaque colonne pour avoir des vecteurs unitaires
  Ur <- apply(U_r_brut, 2, function(col) col / sqrt(sum(col^2)))
  
  # Étape 4 : construire Sigma (matrice rectangulaire m x n)
  m <- nrow(A); n <- ncol(A)
  Sigma <- matrix(0, nrow = m, ncol = n)
  diag(Sigma)[1:r] <- sigma
  
  # Compléter U si besoin (vecteurs orthogonaux arbitraires)
  # On utilise qr() pour ça
  if (ncol(Ur) < m) {
    Ur <- qr.Q(qr(cbind(Ur, matrix(rnorm(m * (m - r)), m, m - r))))
  }
  
  return(list(U = Ur, Sigma = Sigma, V = V, sigma = sigma, r = r))
}

# --- SECTION 3.2 : Calcul SVD pour les matrices du cas 1 ---

A1 <- matrix(c(-18,2,-14,-2, 13,19,11,21, -4,-4,-12,4, 4,12,8,8),
             nrow = 4, byrow = FALSE)
B1 <- matrix(c(6,2,0,-1, -8,7,-1,-2, -4,-5,-8,4, 5,-6,2,4, -4,4,2,-8),
             nrow = 4, byrow = FALSE)

cat("=== SVD de A (cas 1) ===\n")
svd_A1 <- ma_svd(A1)
cat("Valeurs singulières :", round(svd_A1$sigma, 2), "\n")

# Vérification : U * Sigma * t(V) doit redonner A
reconstitution <- svd_A1$U %*% svd_A1$Sigma %*% t(svd_A1$V)
cat("Erreur de reconstruction :", round(max(abs(reconstitution - A1)), 6), "\n")

# Comparaison avec la fonction native R
svd_native <- svd(A1)
cat("Valeurs singulières (R natif) :", round(svd_native$d, 2), "\n")

# SECTION 3.3 : Pseudo-inverse et résolution de systèmes

pseudo_inverse <- function(A, tol = 1e-10) {
  res <- ma_svd(A)
  r   <- res$r
  Ur  <- res$U[, 1:r, drop = FALSE]
  Vr  <- res$V[, 1:r, drop = FALSE]
  D_inv <- diag(1 / res$sigma[1:r], nrow = r)
  
  # A⁺ = Vr * D⁻¹ * Urᵀ
  return(Vr %*% D_inv %*% t(Ur))
}

A2 <- matrix(c(-3,-1,0,0, -3,-1,0,0, -6,-1,-1,1, 6,1,1,1, 1,-2,-1,-1),
             nrow = 4, byrow = FALSE)
B2 <- matrix(c(4,-5,2,6, 0,0,0,0, -1,3,-1,-3, -2,5,-2,-6, 0,0,0,0),
             nrow = 4, byrow = FALSE)

b  <- c(6, -1, -4, 6)

# Résolution Ax = b
A2_plus <- pseudo_inverse(A2)
x_A <- A2_plus %*% b
cat("\n=== Résolution A2·x = b ===\n")
cat("Solution x :", round(x_A, 2), "\n")
cat("Vérification A·x :", round(A2 %*% x_A, 2), "\n")
cat("Résidu ||Ax - b|| :", round(sqrt(sum((A2 %*% x_A - b)^2)), 6), "\n")

# Résolution Bx = b
B2_plus <- pseudo_inverse(B2)
x_B <- B2_plus %*% b
cat("\n=== Résolution B2·x = b ===\n")
cat("Solution x :", round(x_B, 2), "\n")
cat("Vérification B·x :", round(B2 %*% x_B, 2), "\n")
cat("Résidu ||Bx - b|| :", round(sqrt(sum((B2 %*% x_B - b)^2)), 6), "\n")