library(tcor)

# Thresholded correlation (from the tcor examples, but adding a check for
# equivalence with the cor result)

set.seed(1)
s <- svd(matrix(rnorm(100 * 1000), nrow=100))
A <- s$u %*% (1 /( 1:100) * t(s$v)) 
C <- cor(A)
C <- C * upper.tri(C)
y <- which(C >= 0.98, arr.ind=TRUE)
x <- tcor(A, t=0.98)$indices[, 1:2]
# order x and y conformably for comparison
swap <- x[, 1] > x[, 2]
x2   <- x[, 2]
x[swap, 2] <- x[swap, 1]
x[swap, 1] <- x2[swap]
x <- x[order(x[, 2]), ]
y <- y[order(y[, 2]), ]
stopifnot(all.equal(x, y, check.attributes=FALSE))

# Thresholded distance (from the tdist example plus a comparison)
x  <- matrix(rnorm(100 * 20), nrow=100)
td <- tdist(x, 10, rank=TRUE)
d  <- dist(x)
stopifnot(all.equal(td$indices[1:10, "val"], sort(d)[1:10]))

# non-rank version of tdist
td2 <- tdist(x, t=td$indices[10, "val"])
stopifnot(all.equal(td$indices[1:10, "val"], td2$indices[, "val"]))
