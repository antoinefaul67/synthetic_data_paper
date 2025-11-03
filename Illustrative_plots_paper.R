library(rvinecopulib)
library(QRM)
library(GMCM)

## Gaussian copula with different marginal distributions

set.seed(031125)
gaussian_copula_samples = rbicop(n=1000,family="gaussian", parameters = c(0.8))

plot(gaussian_copula_samples)

var1 = qgamma(gaussian_copula_samples[,1], shape=1)

var2 = qnorm(gaussian_copula_samples[,2], mean = 5, sd =1)

bivariate_data = cbind(var1,var2)


p <- ggplot(bivariate_data, aes(var1, var2)) + geom_point() + theme_classic() + xlab("Var1")+
  ylab("Var2")+theme(axis.text=element_text(size=16,face="bold"), axis.title =element_text(size=16,face="bold") )

ggExtra::ggMarginal(p, theme_gray(), type = "histogram",
                    fill = "steelblue", col = "darkblue")

## Gaussian vs Student vs GMCM copulas

set.seed(03112025)

gaussian_copula_samples = rbicop(n=1000,family="gaussian", parameters = c(0.8))


S <- equicorr(d = 2, rho = 0.8)
student_copula_samples <- rcopula.t(1000, df = 4, Sigma = S)


m <- 2 #
theta <- rtheta(m =m, d = 2, method = c("old"))

gmcm_copula_samples <- SimulateGMCMData(n = 10000, d = 2, theta=theta, m=m)$u


par(mfrow = c(1, 3))

# Plot Gaussian Copula Samples
plot(gaussian_copula_samples, main = "Samples from Gaussian Copula",
     xlab = "U1", ylab = "U2", col = rgb(0.2, 0.4, 0.6, 0.7), pch = 16)

# Plot Student Copula Samples
plot(student_copula_samples, main = "Samples from Student Copula",
     xlab = "U1", ylab = "U2", col = rgb(0.2, 0.6, 0.4, 0.7), pch = 16)

# Plot GMCM Samples
plot(gmcm_copula_samples, main = "Samples from GMCM",
     xlab = "U1", ylab = "U2", col = rgb(0.6, 0.2, 0.4, 0.7), pch = 16)

# Reset the plotting layout
par(mfrow = c(1, 1))

