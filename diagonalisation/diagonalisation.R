#   RÉSOLUTION DE SYSTÈMES DE SUITES RÉCURRENTES
#   par Diagonalisation et Jordanisation

library(MASS)
library(Matrix)

cat("  QUESTION 1 : Résolution du système spécifique\n")

A <- matrix(c(1, 0, 2,
              0,-1, 0,
              2, 0, 1), nrow=3, byrow=TRUE)

X0 <- c(1, 2, 3)

cat("Matrice de transition A :\n")
print(A)
cat("\nConditions initiales X0 : ", paste(X0, collapse=" "))

#Valeurs propres et vecteurs propres
eig <- eigen(A)
cat("\nValeurs propres :\n")
print(eig$values)
cat("\nVecteurs propres (colonnes de P) :\n")
print(round(Re(eig$vectors), 8))

#Matrice de passage P
P <- matrix(c( 1, 0, 1,
               0, 1, 0,
               -1, 0, 1), nrow=3, byrow=TRUE)

D_vals <- c(-1, -1, 3)

cat("\nMatrice de passage P :\n")
print(P)

# Vérification : A*P == P*D
D_mat <- diag(D_vals)
cat("\nVérification A*P - P*D (doit être ~0) :\n")
print(round(A %*% P - P %*% D_mat, 10))

#Inverse de P
P_inv <- solve(P)
cat("\nInverse de P :\n")
print(fractions(P_inv))

#Fonction générale : X_n = P * D^n * P^{-1} * X0
solution_q1 <- function(n) {
  Dn <- diag(D_vals^n)
  An <- P %*% Dn %*% P_inv
  round(An %*% X0)
}

# Vérification pour n=0,1,2,3
cat("\nVérification numérique :\n")
for (n in 0:3) {
  Xn <- solution_q1(n)
  print(paste0("n=", n, " : u=", Xn[1], ", v=", Xn[2], ", w=", Xn[3]))
}

# Valeurs pour n=20
cat("\nCalcul pour n=20\n")
val_3_20 <- 3^20
print(paste0("u_20 = 2*3^20 - 1 = ", format(2*val_3_20 - 1, big.mark=" ")))
print(paste0("v_20 = 2"))
print(paste0("w_20 = 2*3^20 + 1 = ", format(2*val_3_20 + 1, big.mark=" ")))
print(paste0("3^20 = ", format(val_3_20, big.mark=" ")))


cat("\n  QUESTION 2 : Méthode générale par Diagonalisation\n")

resoudre_diagonalisation <- function(A, X0, n_values = 0:5, verbose = TRUE) {
  
  k <- nrow(A)
  eig <- eigen(A)
  lambdas <- eig$values
  P <- eig$vectors
  
  if (abs(det(P)) < 1e-10) {
    stop("La matrice A n'est pas diagonalisable ! Utilisez la méthode de Jordanisation (Question 3).")
  }
  
  P_inv <- solve(P)
  
  if (verbose) {
    cat("Matrice A :"); print(A)
    cat("\nValeurs propres :"); print(lambdas)
    cat("\nMatrice de passage P :"); print(round(Re(P), 6))
    cat("\nInverse P^{-1} :"); print(round(Re(P_inv), 6))
  }
  
  resultats <- matrix(NA, nrow=length(n_values), ncol=k+1)
  colnames(resultats) <- c("n", paste0("x", 1:k))
  
  for (i in seq_along(n_values)) {
    n <- n_values[i]
    Dn <- diag(lambdas^n)
    An <- Re(P %*% Dn %*% P_inv)
    Xn <- round(An %*% X0, 8)
    resultats[i,] <- c(n, Xn)
  }
  
  return(as.data.frame(resultats))
}

cat("Application au système de la Question 1 :\n")
res <- resoudre_diagonalisation(A, X0, n_values = 0:5)
print(res)

cat("\n Exemple 2x2\n")
cat("Système : x_{n+1} = 3*x_n + y_n,  y_{n+1} = x_n + 3*y_n\n")
cat("Conditions initiales : x0=1, y0=0\n")

A2 <- matrix(c(3, 1, 1, 3), nrow=2, byrow=TRUE)
X0_2 <- c(1, 0)
res2 <- resoudre_diagonalisation(A2, X0_2, n_values = 0:6)
print(res2)


cat("\n  QUESTION 3 : Cas non diagonalisable — Jordanisation\n")

jordan_block_power <- function(lambda, m, n) {
  J_n <- matrix(0, nrow=m, ncol=m)
  for (i in 1:m) {
    for (j in i:m) {
      k_idx <- j - i
      if (n >= k_idx) {
        binom_coeff <- choose(n, k_idx)
        J_n[i, j] <- binom_coeff * lambda^(n - k_idx)
      }
    }
  }
  return(J_n)
}

cat("Illustration : bloc de Jordan J_3(2)^n pour n=4 :\n")
print(jordan_block_power(lambda=2, m=3, n=4))

cat("\nVérification manuelle : J_3(2)^4 doit donner :\n")
cat("  [ 16,  32,  24 ]\n")
cat("  [  0,  16,  32 ]\n")
cat("  [  0,   0,  16 ]\n")

resoudre_jordanisation <- function(A, X0, n_values = 0:5) {
  
  eig <- eigen(A)
  P <- eig$vectors
  lambdas <- eig$values
  
  if (qr(P)$rank == nrow(P)) {
    P_inv <- solve(P)
    mat_power_diag <- function(n) {
      Dn <- diag(lambdas^n)
      return(Re(P %*% Dn %*% P_inv))
    }
    
  } else {
    cat("Matrice non diagonalisable → utilisation méthode générale (Jordan implicite)\n")
    mat_power_diag <- function(n) {
      if (n == 0) return(diag(nrow(A)))
      result <- diag(nrow(A))
      base <- A
      exp <- n
      while (exp > 0) {
        if (exp %% 2 == 1) result <- result %*% base
        base <- base %*% base
        exp <- exp %/% 2
      }
      return(result)
    }
  }
  
  res <- data.frame()
  for (n in n_values) {
    An <- mat_power_diag(n)
    Xn <- Re(An %*% X0)
    res <- rbind(res, c(n, Xn))
  }
  colnames(res) <- c("n", paste0("x", 1:length(X0)))
  return(res)
}

cat("Exemple : matrice non diagonalisable\n")
cat("A = [[2,1],[0,2]] (bloc de Jordan 2x2, valeur propre lambda=2)\n")

A_nd <- matrix(c(2, 1, 0, 2), nrow=2, byrow=TRUE)
X0_nd <- c(1, 1)

cat("\nCalcul de J_2(2)^n manuellement :\n")
for (n in 0:5) {
  Jn <- jordan_block_power(lambda=2, m=2, n=n)
  Xn <- round(Jn %*% X0_nd)
  print(paste0("n=", n, " : x1=", Xn[1], ",  x2=", Xn[2],
               "   [x1 = 2^n + n*2^(n-1), x2 = 2^n]"))
}

cat("\nFormules analytiques :\n")
cat("  x1_n = 2^n + n*2^(n-1)\n")
cat("  x2_n = 2^n\n")

cat("\nExemple 3x3 non diagonalisable\n")
cat("A = [[3,1,0],[0,3,1],[0,0,3]] (bloc de Jordan 3x3)\n")

A_nd3 <- matrix(c(3,1,0, 0,3,1, 0,0,3), nrow=3, byrow=TRUE)
X0_nd3 <- c(1, 0, 1)

cat("Calcul via J_3(3)^n :\n")
for (n in 0:4) {
  Jn <- jordan_block_power(lambda=3, m=3, n=n)
  Xn <- round(Jn %*% X0_nd3)
  print(paste0("n=", n, " : x=(", Xn[1], ", ", Xn[2], ", ", Xn[3], ")"))
}

cat("\nFormules analytiques (avec X0=(1,0,1)) :\n")
cat("  x1_n = 3^n + C(n,2)*3^(n-2)\n")
cat("  x2_n = n*3^(n-1)\n")
cat("  x3_n = 3^n\n")


cat("\n  RÉCAPITULATIF DES RÉSULTATS — QUESTION 1\n")
cat("  u_n = 2*3^n - (-1)^n\n")
cat("  v_n = 2*(-1)^n\n")
cat("  w_n = 2*3^n + (-1)^n\n")
print(paste0("  u_20 = 2*3^20 - 1 = ", format(2*3^20 - 1, big.mark=" ")))
print(paste0("  v_20 = 2"))
print(paste0("  w_20 = 2*3^20 + 1 = ", format(2*3^20 + 1, big.mark=" ")))