
## vanilla RSA ----

## models for full run-wise similarity structure

## cross-validated RSA ----

## models for condition x condition similarity structure
dnames <- dimnames(vanil.2tr1knot)$.row
m <- diag(length(dnames))
dimnames(m) <- dimnames(vanil.2tr1knot)[c(".row", ".col")]
m[] <- 0
lt <- lower.tri(m)

is.incon <- grepl("InCon", dnames)
is.congr <- grepl("PC50Con|biasCon", dnames)
is.pc50  <- grepl("PC50", dnames)
is.bias  <- grepl("bias", dnames)
is.run1  <- grepl("run1", dnames)
is.run2  <- grepl("run2", dnames)

m.incon <- m
m.congr <- m
m.pc50 <- m
m.bias <- m
m.bias.any <- m
m.run <- m
m.ond <- m
m.offd <- m

m.incon[is.incon, is.incon] <- 1
m.congr[is.congr, is.congr] <- 1
m.pc50[is.pc50, is.pc50] <- 1
m.bias.any[is.bias, ] <- 1
m.bias.any[, is.bias] <- 1
m.bias[is.bias, is.bias] <- 1
m.run[is.run1, is.run1] <- 1
m.run[is.run2, is.run2] <- 1
diag(m.ond[is.run1, is.run2]) <- 1
diag(m.ond[is.run2, is.run1]) <- 1
m.offd <- 1 - m.ond
diag(m.offd) <- 0

mlist <- list(m.incon, m.congr, m.pc50, m.bias, m.bias.any, m.run, m.ond, m.offd)
names(mlist) <- c("incon", "congr", "pc50", "bias", "bias.any", "run", "ond", "offd")

rm(m.incon, m.congr, m.pc50, m.bias, m.bias.any, m.run, m.ond, m.offd)
