
# Exercice 2 - Question 3 : Interpolation de Newton

differences_divisees <- function(x, y) {
  n  <- length(x)
  dd <- matrix(0, nrow = n, ncol = n)
  dd[, 1] <- y
  for (j in 2:n)
    for (i in j:n)
      dd[i,j] <- (dd[i,j-1] - dd[i-1,j-1]) / (x[i] - x[i-j+1])
  return(dd)
}

coeff_newton <- function(x, y) {
  round(diag(differences_divisees(x, y)), 2)
}

eval_newton <- function(c_n, x_noeuds, x_val) {
  sapply(x_val, function(xv) {
    res <- c_n[length(c_n)]
    for (k in (length(c_n)-1):1) res <- res * (xv - x_noeuds[k]) + c_n[k]
    res
  })
}

plot_newton <- function(x, y, titre = "Interpolation de Newton") {
  c_n <- coeff_newton(x, y)
  cat("Newton | n =", length(x)-1, "\n")
  cat("Coefficients :\n") ; print(c_n)

  x_dense <- seq(min(x), max(x), length.out = 500)
  y_poly  <- eval_newton(c_n, x, x_dense)
  plot(x_dense, y_poly, type = "l", col = "darkgreen", lwd = 2,
       main = titre, xlab = "x", ylab = "y",
       ylim = range(c(y, y_poly)))
  points(x, y, pch = 19, col = "steelblue", cex = 1.5)
  legend("topleft", legend = c("Points", "Polynome Newton"),
         col = c("steelblue","darkgreen"), pch = c(19,NA), lty = c(NA,1), lwd = 2)
  grid()
}

# Application n = 9, 19, 29
set.seed(2026)
x9  <- sort(sample(0:60, 10)) ; y9  <- sample(-50:100, 10, replace = TRUE)
x19 <- sort(sample(0:60, 20)) ; y19 <- sample(-50:100, 20, replace = TRUE)
x29 <- sort(sample(0:60, 30)) ; y29 <- sample(-50:100, 30, replace = TRUE)

plot_newton(x9,  y9,  "Newton - n = 9")
plot_newton(x19, y19, "Newton - n = 19")
plot_newton(x29, y29, "Newton - n = 29")

# Application au nuage impose + previsions
x_data <- 0:9
y_data <- c(2, 6, 12, 20, 30, 42, 56, 72, 90, 110)
c_n    <- coeff_newton(x_data, y_data)

x_dense <- seq(0, 20, length.out = 500)
y_poly  <- eval_newton(c_n, x_data, x_dense)
plot(x_dense, y_poly, type = "l", col = "darkgreen", lwd = 2,
     main = "Newton - nuage impose + extrapolation",
     xlab = "x", ylab = "y", ylim = range(c(y_data, y_poly)))
points(x_data, y_data, pch = 19, col = "steelblue", cex = 1.5)
abline(v = 9, lty = 3, col = "gray50")
grid()

cat("Previsions :\n")
for (xp in c(10, 15, 20))
  cat(sprintf("P(%d) = %.2f\n", xp, round(eval_newton(c_n, x_data, xp), 2)))
