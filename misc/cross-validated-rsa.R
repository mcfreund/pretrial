## cross-validated euclidean ----

contrast_matrix <- function(n, condition.names) {
  # n <- 10
  # condition.names <- letters[1:n]
  
  if (n < 2) stop("you need more than 1 condition")
  
  W <- matrix(0, nrow = n^2, ncol = n)
  
  if (missing(condition.names)) {
    dimnames(W) <- list(contrast = NULL, condition = NULL)
  } else {
    dimnames(W) <- list(
      contrast = paste0(rep(condition.names, each = n), "_", rep(condition.names, n)),
      condition = condition.names
    )
  }
  
  
  for (condition.i in seq_len(n)) {
    # condition.i = 1
    
    row.beg <- (condition.i - 1) * n + 1
    row.end <- (condition.i - 1) * n + n
    W.i <- W[row.beg:row.end, ]  ## square matrix
    
    W.i[, condition.i] <- 1
    diag(W.i) <- diag(W.i) - 1
    
    W[row.beg:row.end, ] <- W.i
    
  }
  
  W
  
}


W <- contrast_matrix(10, LETTERS[1:10])
B <- matrix(rnorm(10 * 100), ncol = 100, nrow = 10)
E <- matrix(rnorm(100 * 1E4), ncol = 100)
S <- solve(cov(E))

## standard diagonal of matrix product method:

matprod <- diag(W %*% B %*% t(B) %*% t(W))
matprod.w <- diag(W %*% B %*% S %*% t(B) %*% t(W))

## hamard product method:
# https://stackoverflow.com/questions/42569698/how-to-just-calculate-the-diagonal-of-a-matrix-product-in-r
hamard <- colSums(t(W %*% B) * t(B) %*% t(W))
hamard.w <- colSums(t(W %*% B %*% S) * t(B) %*% t(W))

all.equal(matprod, hamard)
all.equal(matprod.w, hamard.w)


## cross-validated correlation ----

B2 <- matrix(rnorm(10 * 100), ncol = 100, nrow = 10)

B1 <- sweep(B, 2, colMeans(B))  ## center
B2 <- sweep(B2, 2, colMeans(B))
V <- B1 %*% t(B2)  ## inner product matrix
## divide each entry by sqrt(product of corresponding diagonal entries):
## NB: if product of diagonal entries is negative, NaN will be produced with warnings
R <- cov2cor(V)

## now for prewhitening...

Vw <- B1 %*%  S %*% t(B2)  ## inner product matrix
Rw <- cov2cor(Vw)
