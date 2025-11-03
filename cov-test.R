# testing how R computes the sample cross-covariance
a <- matrix(1:15, nrow=5, ncol=3)
b <- matrix(25:39, 5, 3)

c <- cov(a, b)
d <- matrix(0, 3, 3)
for (i in 1:nrow(a)) {
  d = d + 0.25 * (a[i,] - colMeans(a)) %*% t(b[i,] - colMeans(b))
}

print(c) # sample covariance
print(d) # sample covariance computed by hand


