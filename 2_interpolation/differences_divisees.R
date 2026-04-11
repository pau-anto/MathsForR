
# Exercice 2 - Questions 4 & 5 : Differences divisees + polynome canonique

differences_divisees <- function(x, y) {
  n  <- length(x)
  dd <- matrix(0, nrow = n, ncol = n)
  dd[, 1] <- y
  for (j in 2:n)
    for (i in j:n)
      dd[i,j] <- (dd[i,j-1] - dd[i-1,j-1]) / (x[i] - x[i-j+1])
  return(dd)
}

coeff_newton <- function(x, y) round(diag(differences_divisees(x, y)), 2)

eval_newton <- function(c_n, x_noeuds, x_val) {
  sapply(x_val, function(xv) {
    res <- c_n[length(c_n)]
    for (k in (length(c_n)-1):1) res <- res * (xv - x_noeuds[k]) + c_n[k]
    res
  })
}

afficher_dd <- function(x, y) {
  dd  <- differences_divisees(x, y)
  n   <- length(x)
  tab <- as.data.frame(round(dd, 2))
  colnames(tab) <- paste0("Ordre_", 0:(n-1))
  rownames(tab) <- paste0("x", 0:(n-1), "=", x)
  cat("Tableau des differences divisees :\n")
  print(tab)
  cat("Coefficients de Newton (diagonale) :\n")
  print(round(diag(dd), 2))
}

newton_vers_canonique <- function(c_n, x_noeuds) {
  n   <- length(c_n)
  res <- rep(0, n)
  res[1] <- c_n[n]
  for (k in (n-1):1) {
    res_new <- rep(0, n)
    for (j in 1:(n-1)) res_new[j+1] <- res_new[j+1] + res[j]
    for (j in 1:n)     res_new[j]   <- res_new[j] - x_noeuds[k] * res[j]
    res_new[1] <- res_new[1] + c_n[k]
    res <- res_new
  }
  return(round(res, 2))
}

# Question 4 : exemple sur petit nuage
x_ex <- c(0, 1, 3, 4, 9)
y_ex <- c(47, 38, 36, 13, -4)
afficher_dd(x_ex, y_ex)

# Question 5 : nuage impose
x_data <- 0:9
y_data <- c(2, 6, 12, 20, 30, 42, 56, 72, 90, 110)
print(data.frame(x = x_data, y = y_data))
afficher_dd(x_data, y_data)

c_n   <- coeff_newton(x_data, y_data)
a_can <- newton_vers_canonique(c_n, x_data)
cat("Coefficients base canonique [a0, ..., a9] :\n") ; print(a_can)
cat("P(x) = x^2 + 3x + 2\n")

x_dense <- seq(0, 20, length.out = 500)
y_poly  <- eval_newton(c_n, x_data, x_dense)
plot(x_dense, y_poly, type = "l", col = "purple", lwd = 2,
     main = "Q5 - Polynome d'interpolation + previsions",
     xlab = "x", ylab = "P(x)", ylim = range(c(y_data, y_poly)))
points(x_data, y_data, pch = 19, col = "steelblue", cex = 1.5)
points(c(10,15,20), eval_newton(c_n, x_data, c(10,15,20)),
       pch = 17, col = "red", cex = 2)
abline(v = 9, lty = 3, col = "gray50")
legend("topleft", legend = c("Points", "P(x)", "Previsions"),
       col = c("steelblue","purple","red"),
       pch = c(19,NA,17), lty = c(NA,1,NA), lwd = 2)
grid()

cat("Previsions :\n")
for (xp in c(10, 15, 20))
  cat(sprintf("P(%d) = %.2f\n", xp, round(eval_newton(c_n, x_data, xp), 2)))
