library(mikeutils)

n <- 100

a <- runif(n)
# b <- runif(n)
b <- a + runif(n)
m <- (a + b) / 2
mu <- scale2unit(m)
ac <- a - m
bc <- b - m
ssq(ac) == ssq(bc)
# acs <- scale2unit(ac)
# bcs <- scale2unit(bc)
acs <- scale2unit(ac) * 1 / ssq(ac)
bcs <- scale2unit(bc) * 1 / ssq(bc)
all.equal(cosinesim(acs, bcs), -1)
all.equal(c(acs %*% ac), 1)
all.equal(c(bcs %*% bc), 1)

## project new point p onto discriminant
pa <- a * 10000
pb <- b * 10000
# cosinesim(a, b)
# cosinesim(a, p)
# cosinesim(b, p)

## WRONG
all(acs %*% pa > 0, acs %*% pb < 0)

## RIGHT
all(
  acs %*% (pa - c(mu %*% pa %*% mu)) > 0,
  acs %*% (pb - c(mu %*% pb %*% mu)) < 0,
  bcs %*% (pa - c(mu %*% pa %*% mu)) < 0,
  bcs %*% (pb - c(mu %*% pb %*% mu)) > 0
)

## adaptively scaling the mean vector preserves length information of the 'test' dataset.
## all length info can easily be ignored, however, by scaling all input patterns to length 1.
## this procedure would give the same result.

