library(MASS)   # pour fractions()

# 0.  FONCTIONS UTILITAIRES

# Vérification de diagonalisabilité :
#   pour chaque valeur propre λ, dim(Eλ) = m(λ)  <=>  mult. géom. = mult. alg.
est_diagonalisable <- function(A, tol = 1e-6) {
  n    <- nrow(A)
  vals <- round(Re(eigen(A)$values), 6)
  for (lam in unique(vals)) {
    m_alg <- sum(vals == lam)
    m_geo <- n - qr(A - lam * diag(n), tol = tol)$rank
    if (m_geo < m_alg) return(FALSE)
  }
  TRUE
}

# Bloc de Jordan J_m(λ) de taille m
jordan_bloc <- function(lambda, m) {
  J <- lambda * diag(m)
  if (m > 1) for (i in seq_len(m - 1)) J[i, i + 1] <- 1
  J
}

# Puissance d'un bloc de Jordan : (J_m(λ))^n
# Formule : (J^n)_{i,j} = C(n, j-i) · λ^{n-(j-i)}  pour j >= i
jordan_bloc_puissance <- function(lambda, m, n) {
  J_n <- matrix(0, m, m)
  for (i in seq_len(m))
    for (j in i:m) {
      k <- j - i
      if (n >= k) J_n[i, j] <- choose(n, k) * lambda^(n - k)
    }
  J_n
}

# Sous-espace propre E_lambda = Ker(A - lambda*I) via eigen() :
#   les vecteurs propres de M associés à la valeur propre 0 sont exactement
#   les vecteurs du noyau de M. On filtre en vérifiant M·v ≈ 0.
noyau_via_eigen <- function(M, tol = 1e-8) {
  e   <- eigen(M)
  idx <- which(abs(Re(e$values)) < tol)
  if (length(idx) == 0) return(matrix(0, nrow(M), 0))
  candidats <- Re(e$vectors[, idx, drop = FALSE])
  # Conserver uniquement les vrais vecteurs du noyau (résidu M·v ≈ 0)
  ok <- apply(candidats, 2, function(v) max(abs(M %*% v)) < sqrt(tol))
  if (!any(ok)) return(matrix(0, nrow(M), 0))
  candidats[, ok, drop = FALSE]
}


# DÉCOMPOSITION DE JORDAN  A = P · J · P^{-1}
#
#      Pour chaque valeur propre λ :
#        1. Vecteurs propres v1 via eigen()  [noyau de A - λI]
#        2. Vecteurs propres généralisés : résoudre (A - λI)·v_{k+1} = v_k
#           avec qr.solve(), jusqu'à ce qu'aucune extension ne soit possible.
#        Chaque chaîne de longueur d  <=>  un bloc de Jordan J_d(λ) dans J
decomposer_jordan <- function(A, tol = 1e-8) {
  n    <- nrow(A)
  vals <- round(Re(eigen(A)$values), 6)
  
  P_mat   <- matrix(0, n, 0)
  J_blocs <- list()
  
  for (lam in unique(vals)) {
    m_alg <- sum(vals == lam)
    M     <- A - lam * diag(n)
    
    # Vecteurs propres : noyau de M via eigen()
    V_base <- noyau_via_eigen(M, tol)
    m_geo  <- ncol(V_base)
    
    for (j in seq_len(m_geo)) {
      
      # --- Démarrer la chaîne avec le vecteur propre v1 ---
      chain  <- list(V_base[, j, drop = FALSE])
      v_prev <- chain[[1]]
      
      # --- Étendre la chaîne : résoudre (A - λI)·v_next = v_prev ---
      while (length(chain) < m_alg) {
        
        v_next <- tryCatch(
          qr.solve(M, v_prev, tol = tol),
          error = function(e) NULL
        )
        
        # Arrêter si le système n'a pas de solution valide
        if (is.null(v_next)) break
        if (max(abs(M %*% v_next - v_prev)) > sqrt(tol)) break
        
        # Arrêter si v_next est linéairement dépendant des vecteurs déjà dans P
        P_actuel <- cbind(P_mat, do.call(cbind, chain))
        if (ncol(P_actuel) > 0) {
          Q_P    <- qr.Q(qr(P_actuel))
          v_comp <- v_next - Q_P %*% (t(Q_P) %*% v_next)
          if (sqrt(sum(v_comp^2)) < sqrt(tol)) break
        }
        
        chain  <- c(chain, list(v_next))
        v_prev <- v_next
      }
      
      P_mat   <- cbind(P_mat, do.call(cbind, chain))
      J_blocs <- c(J_blocs, list(jordan_bloc(lam, length(chain))))
    }
  }
  
  # Assembler J en matrice bloc-diagonale
  sz    <- sum(sapply(J_blocs, nrow))
  J_mat <- matrix(0, sz, sz)
  idx   <- 1
  for (bl in J_blocs) {
    s <- nrow(bl)
    J_mat[idx:(idx+s-1), idx:(idx+s-1)] <- bl
    idx <- idx + s
  }
  
  list(P = Re(P_mat), J = J_mat, blocs = J_blocs)
}

# Calcule A^n = P · J^n · P^{-1}  (J^n calculé bloc par bloc)
puissance_par_jordan <- function(decomp, n) {
  P <- decomp$P; blocs <- decomp$blocs; N <- nrow(P)
  J_n <- matrix(0, N, N)
  idx <- 1
  for (bl in blocs) {
    s <- nrow(bl)
    J_n[idx:(idx+s-1), idx:(idx+s-1)] <- jordan_bloc_puissance(bl[1,1], s, n)
    idx <- idx + s
  }
  Re(P %*% J_n %*% solve(P))
}





cat(" QUESTION 1 : Système spécifique (diagonalisation explicite)\n")

A  <- matrix(c(1, 0, 2,
               0,-1, 0,
               2, 0, 1), nrow = 3, byrow = TRUE)
X0 <- c(1, 2, 3)

cat("Matrice de transition A :\n"); print(A)
cat("Conditions initiales X0 :", X0, "\n\n")

# Polynôme caractéristique : det(A - λI) = (1+λ)²(λ-3)
# Valeurs propres : λ = -1 (mult. alg. 2) et λ = 3 (mult. alg. 1)
eig <- eigen(A)
cat("Valeurs propres :", round(Re(eig$values), 2), "\n")

# Sous-espaces propres :
#   E_{-1} = Ker(A+I) = vect{(1,0,-1), (0,1,0)}  dim=2=m(-1) => diagonalisable
#   E_3    = Ker(A-3I) = vect{(1,0,1)}             dim=1=m(3)  => diagonalisable
cat("Diagonalisable ? (dim(Eλ) = m(λ) pour tout λ) :",
    est_diagonalisable(A), "\n\n")

# Matrice de passage P et matrice diagonale D
P     <- matrix(c( 1, 0, 1,
                   0, 1, 0,
                   -1, 0, 1), nrow = 3, byrow = TRUE)
D_val <- c(-1, -1, 3)
D_mat <- diag(D_val)
P_inv <- solve(P)

cat("Matrice de passage P :\n"); print(P)
cat("\nP^{-1} :\n"); print(fractions(P_inv))
cat("\nVérification  A·P = P·D  (erreur max) :",
    max(abs(A %*% P - P %*% D_mat)), "\n\n")

# X_n = P · D^n · P^{-1} · X0  avec  D^n = diag((-1)^n, (-1)^n, 3^n)
solution_q1 <- function(n) {
  An <- P %*% diag(D_val^n) %*% P_inv
  round(Re(An %*% X0), 2)
}

cat("Formules analytiques :\n")
cat("  u_n = 2·3^n  -  (-1)^n\n")
cat("  v_n = 2·(-1)^n\n")
cat("  w_n = 2·3^n  +  (-1)^n\n\n")

cat("Vérification numérique (n = 0 à 4) :\n")
for (n in 0:4) {
  Xn <- solution_q1(n)
  cat(sprintf("  n=%d : u=%.2f,  v=%.2f,  w=%.2f\n", n, Xn[1], Xn[2], Xn[3]))
}

cat("\nValeurs pour n = 20 :\n")
cat(sprintf("  u_20 = 2·3^20 - 1  = %s\n", format(2*3^20 - 1, big.mark = " ")))
cat(sprintf("  v_20 = 2\n"))
cat(sprintf("  w_20 = 2·3^20 + 1  = %s\n", format(2*3^20 + 1, big.mark = " ")))





cat(" QUESTION 2 : Méthode générale par diagonalisation\n")

resoudre_diagonalisation <- function(A, X0, n_values = 0:5, verbose = TRUE) {
  k <- nrow(A)
  
  # Test de diagonalisabilité (méthode du cours : dim(Eλ) = m(λ))
  if (!est_diagonalisable(A))
    stop("A n'est pas diagonalisable. Utilisez resoudre_jordanisation() (Q3).")
  
  eig     <- eigen(A)
  lambdas <- Re(eig$values)
  P       <- Re(eig$vectors)
  P_inv   <- solve(P)
  
  if (verbose) {
    cat("Matrice A :\n"); print(A)
    cat("\nValeurs propres :", round(lambdas, 2), "\n")
    cat("Matrice de passage P :\n"); print(round(P, 6))
    cat("Vérification A·P = P·D (erreur max) :",
        round(max(abs(A %*% P - P %*% diag(lambdas))), 10), "\n\n")
  }
  
  res <- lapply(n_values, function(n) {
    An <- Re(P %*% diag(lambdas^n) %*% P_inv)
    c(n = n, round(An %*% X0, 2))
  })
  
  df <- as.data.frame(do.call(rbind, res))
  colnames(df) <- c("n", paste0("x", seq_len(k)))
  df
}

cat("--- Application au système de la Question 1 ---\n")
res_q1 <- resoudre_diagonalisation(A, X0, n_values = 0:5, verbose = FALSE)
print(res_q1)

cat("\n--- Exemple 2×2 : x_{n+1}=3x_n+y_n,  y_{n+1}=x_n+3y_n,  X0=(1,0) ---\n")
A2   <- matrix(c(3, 1, 1, 3), nrow = 2, byrow = TRUE)
X0_2 <- c(1, 0)
res_q2 <- resoudre_diagonalisation(A2, X0_2, n_values = 0:6)
print(res_q2)





cat(" QUESTION 3 : Cas non diagonalisable — Jordanisation\n")

# --- 3a. Illustration de la formule sur un bloc de Jordan ---
cat("Illustration : J_3(2)^4  (formule du cours)\n")
print(jordan_bloc_puissance(2, 3, 4))
cat("Résultat attendu : [16,32,24 | 0,16,32 | 0,0,16]\n\n")

# --- 3b. Fonction générale par jordanisation ---
resoudre_jordanisation <- function(A, X0, n_values = 0:5, verbose = TRUE) {
  decomp <- decomposer_jordan(A)
  
  if (verbose) {
    cat("Matrice A :\n"); print(A)
    cat("\nForme de Jordan J  (A = P·J·P^{-1}) :\n")
    print(round(decomp$J, 6))
    cat("\nMatrice de passage P :\n")
    print(round(decomp$P, 6))
    erreur <- max(abs(A - decomp$P %*% decomp$J %*% solve(decomp$P)))
    cat("Vérification A = P·J·P^{-1} (erreur max) :", round(erreur, 10), "\n\n")
  }
  
  res <- lapply(n_values, function(n) {
    An <- puissance_par_jordan(decomp, n)
    c(n = n, round(An %*% X0, 2))
  })
  
  df <- as.data.frame(do.call(rbind, res))
  colnames(df) <- c("n", paste0("x", seq_len(length(X0))))
  df
}

# --- 3c. Exemple 1 : bloc de Jordan 2×2 ---
cat("=== Exemple 1 : A = J_2(2)  (2×2, non diagonalisable) ===\n")
A_nd2  <- matrix(c(2, 1, 0, 2), nrow = 2, byrow = TRUE)
X0_nd2 <- c(1, 1)
cat("Diagonalisable ?", est_diagonalisable(A_nd2), "\n\n")

res_nd2 <- resoudre_jordanisation(A_nd2, X0_nd2, n_values = 0:5)
print(res_nd2)
cat("\nFormules analytiques pour X0=(1,1) :\n")
cat("  x1_n = 2^n + n·2^(n-1)\n")
cat("  x2_n = 2^n\n\n")

# --- 3d. Exemple 2 : bloc de Jordan 3×3 ---
cat("=== Exemple 2 : A = J_3(3)  (3×3, non diagonalisable) ===\n")
A_nd3  <- matrix(c(3, 1, 0,
                   0, 3, 1,
                   0, 0, 3), nrow = 3, byrow = TRUE)
X0_nd3 <- c(1, 0, 1)
cat("Diagonalisable ?", est_diagonalisable(A_nd3), "\n\n")

res_nd3 <- resoudre_jordanisation(A_nd3, X0_nd3, n_values = 0:4)
print(res_nd3)
cat("\nFormules analytiques pour X0=(1,0,1) :\n")
cat("  x1_n = 3^n + C(n,2)·3^(n-2)\n")
cat("  x2_n = n·3^(n-1)\n")
cat("  x3_n = 3^n\n\n")

# --- 3e. Vérification croisée : système Q1 via jordanisation ---
cat("=== Vérification : système Q1 résolu par jordanisation ===\n")
res_q1_jordan <- resoudre_jordanisation(A, X0, n_values = 0:5, verbose = FALSE)
print(res_q1_jordan)
cat("\nÉcart max diagonalisation vs jordanisation :",
    round(max(abs(res_q1[,-1] - res_q1_jordan[,-1])), 10), "\n")





cat(" RÉCAPITULATIF — QUESTION 1\n")
cat("  Valeurs propres : λ₁=λ₂=-1  (E_{-1}=vect{(1,0,-1),(0,1,0)})\n")
cat("                    λ₃=3       (E_3=vect{(1,0,1)})\n\n")
cat("  X_n = P · diag((-1)^n, (-1)^n, 3^n) · P^{-1} · X0\n\n")
cat("  u_n = 2·3^n  -  (-1)^n\n")
cat("  v_n = 2·(-1)^n\n")
cat("  w_n = 2·3^n  +  (-1)^n\n\n")
cat(sprintf("  u_20 = 2·3^20 - 1 = %s\n", format(2*3^20 - 1, big.mark = " ")))
cat(sprintf("  v_20 = 2\n"))
cat(sprintf("  w_20 = 2·3^20 + 1 = %s\n", format(2*3^20 + 1, big.mark = " ")))