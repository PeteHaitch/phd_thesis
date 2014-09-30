#### Trying to understand Eckhart et al. (2006) ####
# Peter Hickey
# 30/09/2014
# There are attempts to understand and derive the "percentage identical
# methylation" measurement that is shown in Figure 3a as evidence of
# co-methylation.

#### TODOs ####
# Tidy up into a .Rmd file

n1 <- 100
n2 <- 100
prob <- c(0.5, 0.5)
x <- sample(c(0, 1), size = n1, replace = TRUE, prob = prob)
y <- sample(c(0, 1), size = n2, replace = TRUE, prob = prob)
# x <- c(rep(0, n1), rep(1, n2))
# y <- c(rep(1, n1), rep(0, n2))

z <- matrix(rep(NA, n1 * n2), ncol = n2)

for (i in seq_len(n1)) {
  for (j in seq_len(n2)) {
    z[i, j]  <- x[i] == y[j]
  }
}

sum(z) / length(z)

p <- seq(0, 1, by = 0.01)
pp <- expand.grid(p, p)
ppp <- (1 - pp$Var1) * (1 - pp$Var2) + pp$Var1 * pp$Var2


f <- function(p, c) {
  ifelse(pmax(1 - 1/ p, 1 - 1 / (1 - p)) <= c,
         2 * p * (1 - c) * (p - 1) + 1,
         NA)
}

f(seq(0, 1, 0.01), 0)
f(seq(0, 1, 0.01), 0.5)
f(seq(0, 1, 0.01), 1)
f(seq(0, 1, 0.01), -0.5)
f(seq(0, 1, 0.01), -1)

g <- function(p1, p2, c) {
  # TODO: What is the sufficient condition on the value of c?
  2 * p1 * p2 + 2 * c * p1 * (1 - p1) - (p1 + p2) + 1
}

g(pp$Var1, pp$Var2, 0)
g(p, p, 0)
stopifnot(all.equal(f(p, 0), g(p, p, 0)))
g(pp$Var1, pp$Var2, 0.5)
g(p, p, 0.5)
stopifnot(all.equal(f(p, 0.5), g(p, p, 0.5)))
g(pp$Var1, pp$Var2, 1)
stopifnot(all.equal(f(p, 1), g(p, p, 1)))
g(pp$Var1, pp$Var2, -0.5)
g(p, p, -0.5)
stopifnot(all.equal(f(p, -0.5), g(p, p, -0.5)))
g(pp$Var1, pp$Var2, -1)
stopifnot(all.equal(f(p, -1), g(p, p, -1)))
# TODO: General case but with p1, p2 > 0.5

# TODO: Plot density of PIM
