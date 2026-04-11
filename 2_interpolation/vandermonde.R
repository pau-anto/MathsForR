
# Exercice 2 - Question 1 & 2 : 

generer_nuage <- function(n, a, b, c, d) {
  if (b - a < n) stop("Intervalle [a,b] trop petit.")
  x <- sort(sample(a:b, n + 1))
  y <- sample(c:d, n + 1, replace = TRUE)
  plot(x, y, pch = 19, col = "steelblue", cex = 1.5,
       main = paste0("Nuage de points Delta (n = ", n, ")"),
       xlab = "x", ylab = "y",
       xlim = c(a-1, b+1), ylim = c(c-2, d+2))
  grid()
  return(data.frame(x = x, y = y))
}

vandermonde <- function(x) {
  n <- length(x)
  V <- matrix(0, nrow = n, ncol = n)
  for (j in 1:n) V[, j] <- x^(j - 1)
  return(V)
}

coeff_vandermonde <- function(x, y) {
  round(solve(vandermonde(x), y), 2)
}

eval_poly <- function(a, x_val) {
  sapply(x_val, function(xv) sum(a * xv^(0:(length(a)-1))))
}

plot_vandermonde <- function(x, y, titre = "Interpolation de Vandermonde") {
  V <- vandermonde(x)
  cat("Matrice de Vandermonde (n =", length(x)-1, ")\n")
  print(round(V, 2))
  cat("Conditionnement :", round(kappa(V), 2), "\n")

  a <- tryCatch(coeff_vandermonde(x, y), error = function(e) {
    cat("Matrice numeriquement singuliere\n")
    return(NULL)
  })

  if (!is.null(a)) {
    cat("Coefficients [a0, ..., an] :\n")
    print(a)
    x_dense <- seq(min(x), max(x), length.out = 500)
    y_poly  <- eval_poly(a, x_dense)
    plot(x_dense, y_poly, type = "l", col = "firebrick", lwd = 2,
         main = titre, xlab = "x", ylab = "y",
         ylim = range(c(y, y_poly)))
    points(x, y, pch = 19, col = "steelblue", cex = 1.5)
    legend("topleft", legend = c("Points", "Polynome Vandermonde"),
           col = c("steelblue","firebrick"), pch = c(19,NA), lty = c(NA,1), lwd = 2)
    grid()
  }
}

# Application n = 9, 19, 29
set.seed(2026)
x9  <- sort(sample(0:60, 10)) ; y9  <- sample(-50:100, 10, replace = TRUE)
x19 <- sort(sample(0:60, 20)) ; y19 <- sample(-50:100, 20, replace = TRUE)
x29 <- sort(sample(0:60, 30)) ; y29 <- sample(-50:100, 30, replace = TRUE)

plot_vandermonde(x9,  y9,  "Vandermonde - n = 9")
plot_vandermonde(x19, y19, "Vandermonde - n = 19")
plot_vandermonde(x29, y29, "Vandermonde - n = 29")

# Application au nuage impose
x_data <- 0:9
y_data <- c(2, 6, 12, 20, 30, 42, 56, 72, 90, 110)
plot_vandermonde(x_data, y_data, "Vandermonde - nuage impose")
