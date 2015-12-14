context("Standardize Covariates")

# 
var.names <- c("a", "b", "c")
x <- expand.grid(site = 1:4, rep = 1:4)
a_mean <- 1
b_mean <- 5
c_mean <- -2.2
a_sd <- 0.1
b_sd <- 1
c_sd <- 0.2
df_stds <- data.frame(var.names, c(a_mean, b_mean, c_mean), c(a_sd, b_sd, c_sd))
names(df_stds) <- c("var.names", "means", "stdevs")
x$a <- rnorm(nrow(x), a_mean, a_sd)
x$b <- rnorm(nrow(x), b_mean, b_sd)
x$c <- rnorm(nrow(x), c_mean, c_sd)

x_std <- stdCovs(x=x, y=df_stds, var.names = var.names)

# test that all variable names are present in x - change to expect error if not all true
test_that("all variables are present in dataframe", {
  expect_true(all(var.names %in% names(x)))
})

# test handling of factors

# test that column means approximate 0

# test that column sds approximate 1

