## 00_functions

# as per ## https://handbook-5-1.cochrane.org/chapter_7/table_7_7_a_formulae_for_combining_groups.htm

CombMean <- function(n1, m1, n2, m2){
  top <- n1*m1 + n2*m2
  top/(n1+n2)
}

CombSd <- function(n1, m1, n2, m2, s1, s2){
  top = (n1-1)*s1^2 + (n2-1)*s2^2 + (n1*n2/(n1+n2)) * (m1^2 + m2^2 - 2*m1*m2)
  bottom = n1 + n2 -1
  res = (top/bottom)^0.5
  res
}


## test this, it is correct
a <- rnorm(10, sample(1:100, 1), sample(5:15, 1))
b <- rnorm(10, sample(1:100, 1), sample(5:15, 1))
ab <- c(a, b)
means <- map_dbl(list(a, b, ab), mean)
sds <- map_dbl(list(a, b, ab), sd)
mean_res <- CombMean(n1 = 10, m1 = means[1], n2 = 10, m2 = means[2])
sd_res <- CombSd(n1 = 10, m1 = means[1], n2 = 10, m2 = means[2], s1 = sds[1], s2 = sds[2])

means
sds
mean_res
sd_res

print(round(means[3],5) == round(mean_res,5))
print(round(sds[3],5) == round(sd_res,5))
