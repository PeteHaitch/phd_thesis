#### Exploring EM to estimate within-fragment co-methylation ####
# Peter Hickey
# 29/09/2014

#### Some example data ####
n <- 1000
x <- data.frame(l1 = sample(c(0L, 1L, NA_integer_), size = n,
                            prob = c(0.65, 0.3, 0.1), replace = TRUE),
                l2 = sample(c(0L, 1L, NA_integer_), size = n,
                            prob = c(0.65, 0.3, 0.1), replace = TRUE),
                l3 = sample(c(0L, 1L, NA_integer_), size = n,
                            prob = c(0.65, 0.3, 0.1), replace = TRUE)
                )
xx <- xtabs(~ ., data = x, exclude = NULL, na.action = na.pass)

ll <- loglm(~ l1 + l2 + l3, xx)
