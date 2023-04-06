require(mvtnorm)
require(rethinking)

S <- rlkjcorr(1, 4, 1)
pm <- sin(S*(pi/2))
# Generate random correlation matrix
R <- matrix(0, ncol = 4, nrow = 4)
for(i in 1:(n-1)) {
  for(j in (i+1):n) {
    R[i,j] <- runif(1, 0.3, 1)
    R[j,i] <- R[i,j]
  }
}
diag(R) <- 1

#pear <- sin(0.5*pi*S)
S <- rlkjcorr(1, 4, 1)

# so I could generate a correlation matrix
# then generate correlations Xrel-Y
# then simulate with a copula to have such population
# could also compute 
difs_sam <- array(NA, dim=c(4,4,1000))
difs_pop <- array(NA, dim=c(4,4,1000))

for (i in 1:1000){
  AB <- MASS::mvrnorm(1000,c(0,0,0,0),S)
  U <- pnorm(AB)
  x <- qgamma(U[,1],shape=0.2)
  y <- qnorm(U[,2],0,1)
  b1 <- qbinom(U[,3],1,0.5)
  b2 <- qlogis(U[,4,drop=F],0,1)
  cts <- genOrd(b2, t(as.matrix(params$popProbsRel[[1]][1,])),1000)
  ABC <- cbind(x, y, b1, b2)
  ABC2 <- cbind(z, y, b1, b2)
  difs_sam[,,i] <- cor(AB, method="kendall") - cor(ABC,method = "kendall")
  difs_pop[,,i] <- 2/pi*asin(S) - cor(ABC,method = "kendall")
}

apply(difs_sam, c(1,2), range)
apply(difs_pop,c(1,2),range)

# check what are the differences with the package mvdc
# have exact R^2 in y each time
# or generate with assuming a distribution on Y (so add correlations to pop matrix)
# but then R^2 is also not preserved, since only non-parametric correlations are preserved
# which means in any case the population R^2 is not there
# but that's actually fine: I don't want linearity, that's why I'm doing all of this
# I just want some dependence between X and Y
# and if I know pop Kendall's correlation, that's fine

# for some reason R does not work: try tomorrow with only scaling
# and see if Y has enough of importance, else just scale them and add 0 as threshold
# and that's it

pm <- sin(S*(pi/2))

cor(AB)
S

#cor(AB)[1,2]
U <- pnorm(AB)
cor(U, method = "kendall")
#plot(1:1000,U[,1])
#hist(U[,1])
x <- qgamma(U[,1],shape=0.2)
y <- qnorm(U[,2],0,1)
b1 <- qbinom(U[,3],1,0.5)
b2 <- qlogis(U[,4,drop=F],0,1)
cts <- genOrd(b2, t(as.matrix(params$popProbsRel[[1]][1,])),1000)
ABC <- cbind(x, y, b1, cts)

ABC <- cbind(x, y, b1, b2)
ABC2 <- cbind(z, y, b1, b2)
cor(ABC,method = "kendall")
cor(ABC2,method = "kendall")
#cor(x,y,method = "kendall")
Sp <- rlkjcorr(1, 4, 0.1) # pop mat kendall
Sk <- 2/pi*asin(Sp)
difs <- replicate(200, matrix(NA, nrow=5, ncol=5))
S <- rlkjcorr(1, 5, 0.1) # pop mat kendall
diag(S) <- 1
cors <- matrix(NA, nrow=200, ncol=2)
for (i in 1:200){
  AB <- rmvnorm(mean=c(0,0,0,0,0),sig=S,n=1000)
  cors[i, 1] <- cor(AB[,1:2])[1,2]
  U <- pnorm(AB)
  x <- qgamma(U[,1],2)
  y <- qbeta(U[,2],1,2)
  cors[i, 2] <- cor(x, y)[1]
  b1 <- qbinom(U[,3],1,0.5)
  b2 <- qbinom(U[,4],1,0.5)
  z <- ifelse(U[,5] <= 0.2, 0, 
              ifelse(U[,1] <= 0.5, 1, 
                     ifelse(U[,1] < 0.6, 2,
                            ifelse(U[,1] < 0.8, 3, 4))))
  ABC2 <- cbind(x, y, b1, b2, z)
  #difs[,,i] <- cor(U, method = "kendall") - cor(ABC2, method = "kendall")
}

plot(cors[,1], cors[,2])
# transform uniform probabilities to MV probit probabilities to then 
# generate R????
apply(difs, c(1,2), max)
cor(b1, b2, method = "kendall")
cor(b1, x, method = "kendall")
cor(x, b2, method = "kendall")
cor(b1, y, method = "kendall")
cor(y, b2, method = "spearman")
#plot(x,z)
#cor(x, z)
#cor(x, y)
#S

obj <- BiCop(family = 2, par = 0.4, par2 = 6)
obj
## see the object's content or a summary
str(obj)
summary(obj)
## a selection of functions that can be used with BiCop objects
simdata <- BiCopSim(300, obj)
simdata


library(copula)
marginals <- c("norm", "norm", "binom", "binom")
corr_mat <- matrix(c(1.00, 0.50, -0.30, 0.20,
                     0.50, 1.00, -0.40, 0.10,
                     -0.30, -0.40, 1.00, -0.25,
                     0.20, 0.10, -0.25, 1.00), nrow=4, ncol=4)
Sigma <- matrix(c(
  1, 0.8, 0.4, 0.4,
  0.8, 1, 0.6, 0.4,
  0.4, 0.6, 1, 0.8,
  0.4, 0.4, 0.8, 1
), nrow = 4)
S <- rlkjcorr(1, 4, 1) # this is a matrix of kendall taus
vec_parasm <- sin(pi/2*S) # this is avector of copula parameters
# but it's correct for bivariate copulas only
gaussian_copula <- normalCopula(param=vec_parasm[lower.tri(vec_parasm)], dim=4, dispstr = "un")
correlated_data <- rCopula(1000, gaussian_copula)
cor(correlated_data, method = "kendall")

# Set the desired correlation matrix
Sigma <- matrix(c(
  1, 0.8, 0.4, 0.4,
  0.8, 1, 0.6, 0.4,
  0.4, 0.6, 1, 0.8,
  0.4, 0.4, 0.8, 1
), nrow = 4)

# Use Cholesky decomposition to get the lower triangular matrix L such that Sigma = L %*% t(L)
L <- chol(S)

# Set the number of observations
n <- 1000

# Generate the standard normal variables
Z <- matrix(rnorm(n * ncol(Sigma)), nrow = n)

# Generate the correlated variables
correlated_vars <- (t(L %*% t(Z)))

# Extract the binary and continuous variables from the correlated variables
binary_vars <- (correlated_vars[, 1:2] > 0) + 0
continuous_vars <- correlated_vars[, 3:4]

cor(cbind(binary_vars, continuous_vars))
cor(correlated_vars)







library(copula)

# Define the number of binary and continuous variables
num_binary_vars <- 2
num_continuous_vars <- 2

# Define the marginal distributions for the binary and continuous variables
marginals <- c("norm", "norm", "binom", "logis")
params <- list(list(mean = 0, sd = 1),
  list(mean = 0, sd = 1),
  list(size = 1, prob = 0.5),
  list(location = 0, scale = 1)
)

cors <- replicate(1000, matrix(NA, nrow=4, ncol=4))
# Define the copula and the correlation matrix
for (i in 1:1000){
  copula <- normalCopula(param = c(0.8, 0.4, 0.4,
                                   0.6, 0.4, 0.8), dim = 4,
                         dispstr = "un") # params should correspond to the correlations
  # Generate the joint distribution using the copula and the marginal distributions
  joint <- mvdc(copula, marginals, params)
  simulated_data <- rMvdc(1000, joint)
  # population kendall correlation matrix - the observed matrix kendall
  cprop <- matrix(cumsum(baseprobs), ncol=4)
  quant <- t(apply(cprop, 1, stats::qlogis))
  
  iLogisZ <- simulated_data[,4]
  n_obs_i <- length(iLogisZ)
  matlp <- matrix(rep(quant, n),
                  ncol = ncol(cprop),
                  byrow = TRUE
  )
  
  locateGrp <- (iLogisZ > cbind(-Inf, matlp))
  simulated_data[,4] <- apply(locateGrp, 1, sum)
    
  cors[,,i] <- 2/pi*asin(Sigma) - cor(simulated_data)
}
apply(cors, c(1,2), min)
apply(cors, c(1,2), max)
apply(cors, c(1,2), mean)


copula <- normalCopula(param = c(0.8, 0.4, 0.4,
                                 0.6, 0.4, 0.8), dim = 4,
                       dispstr = "un") # params should correspond to the correlations
# find a function to transform kendall correlation
Sigma <- matrix(c(
  1, 0.8, 0.4, 0.4,
  0.8, 1, 0.6, 0.4,
  0.4, 0.6, 1, 0.8,
  0.4, 0.4, 0.8, 1
), nrow = num_binary_vars + num_continuous_vars)


# Generate the joint distribution using the copula and the marginal distributions
joint <- mvdc(copula, marginals, params)
simulated_data <- rMvdc(1000, joint)
cor(simulated_data, method = "kendall")

# Extract the binary and continuous variables from the simulated data
binary_vars <- simulated_data[, 1:num_binary_vars]
continuous_vars <- simulated_data[, (num_binary_vars + 1):(num_binary_vars + num_continuous_vars)]

cor(simulated_data, method = "kendall")
2/pi*asin(Sigma)



# Set the number of variables
n <- 4

# Set the desired correlation matrix
corr_matrix <- matrix(c(
  1.00, 0.80, 0.50, 0.20,
  0.80, 1.00, 0.30, -0.10,
  0.50, 0.30, 1.00, -0.50,
  0.20, -0.10, -0.50, 1.00
), nrow = n)

# Generate random data with the specified correlation matrix
set.seed(123) # for reproducibility
data <- matrix(NA, nrow = 1000, ncol = n)
for (i in 1:n) {
  for (j in 1:n) {
    if (i < j) {
      # Generate the correlation coefficient
      rho <- corr_matrix[i, j]
      
      # Generate random data for the two variables
      if (i <= 2 & j <= 2) {
        x <- rnorm(1000)
        y <- rnorm(1000)
        data[, i] <- ifelse(x > 0, 1, 0)
        data[, j] <- ifelse(y > rho*x, 1, 0)
      } else if (i > 2 & j > 2) {
        x <- rnorm(1000)
        y <- rnorm(1000)
        data[, i] <- cut(x, quantile(x, probs = seq(0, 1, length.out = 5)), include.lowest = TRUE, labels = FALSE)
        data[, j] <- cut(y, quantile(y, probs = seq(0, 1, length.out = 5)), include.lowest = TRUE, labels = FALSE)
      } else {
        x <- rnorm(1000)
        y <- rnorm(1000)
        data[, i] <- cut(x, quantile(x, probs = seq(0, 1, length.out = 4)), include.lowest = TRUE, labels = FALSE)
        data[, j] <- ifelse(y > rho*(data[, i] - 1/2), 1, 0)
      }
    }
  }
}

# Check the correlations
cor(data, use = "pairwise.complete.obs", method = "pearson")
2/pi*asin(corr_matrix)


# Function to generate random correlation matrix
random_corr_matrix <- function(n) {
  A <- matrix(runif(n*n), n, n)
  R <- t(A) %*% A
  diag(R) <- 1
  R <- R/sqrt(diag(R) %*% t(diag(R)))
  return(R)
}

# Set seed for reproducibility
# Function to generate random correlation matrix
random_corr_matrix <- function(n) {
  A <- matrix(runif(n*n), n, n)
  R <- t(A) %*% A
  diag(R) <- 1
  R <- R/sqrt(diag(R) %*% t(diag(R)))
  return(R)
}

# Set seed for reproducibility
set.seed(123)

# Number of variables
n <- 5

# Generate random correlation matrix
R <- random_corr_matrix(n)

# Generate random data with specified correlation matrix
data <- matrix(0, nrow = 1000, ncol = n)

for(i in 1:n) {
  for(j in 1:n) {
    if(i < j) {
      data[, i] <- rnorm(1000)
      data[, j] <- rnorm(1000)
      data[, i] <- scale(data[, i])
      data[, j] <- scale(data[, j])
      data[, i] <- qnorm(pnorm(data[, i]) * R[i, j] + pnorm(rnorm(1000)) * sqrt(1 - R[i, j]))
      data[, j] <- qnorm(pnorm(data[, j]) * R[i, j] + pnorm(rnorm(1000)) * sqrt(1 - R[i, j]))
    }
  }
}

# Check the correlation matrix of the generated data
cor(data)


# Define the correlation matrix
cor_mat <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.3, 0.5, 0.3, 1), nrow = 3)

# Create the normal copula
library(copula)
norm_cop <- normalCopula(0.5)

# Set the sample size
n <- 1000

# Generate uniform random numbers
u <- matrix(runif(3*n), ncol = 3)

# Transform uniform random numbers to Gaussian
g <- qnorm(u)

# Generate correlated continuous variables
x <- rnorm(n, mean = 0, sd = 1)
y <- rnorm(n, mean = 0, sd = 1)
z <- rnorm(n, mean = 0, sd = 1)

# Generate correlated binary variables
b1 <- as.numeric(rnorm(n, mean = 0.5 + 0.2*x + 0.2*y + 0.1*z, sd = 1) > 0)
b2 <- as.numeric(rnorm(n, mean = 0.5 + 0.2*x + 0.1*y + 0.2*z, sd = 1) > 0)

# Combine continuous and binary variables into a data frame
df <- data.frame(x, y, z, b1, b2)

# Apply the normal copula to the data
df_g <- apply(df, 2, function(x) qnorm(pnorm(x)))
df_cop <- norm_cop(df_g)

# Transform Gaussian random numbers to uniform
u_cop <- apply(df_cop, 2, function(x) pnorm(x))
u_cop_b <- u_cop[, 4:5]

# Generate correlated binary variables using uniform random numbers from the copula
b1_cop <- as.numeric(runif(n) < u_cop_b[, 1])
b2_cop <- as.numeric(runif(n) < u_cop_b[, 2])

# Combine continuous and binary variables into a data frame
df_final <- data.frame(x = df[, 1], y = df[, 2], z = df[, 3], b1 = b1_cop, b2 = b2_cop)

# Check the correlation matrix of the final data frame
cor(df_final)


# Function to generate random correlation matrix with lower bound
# n: size of correlation matrix
# lower_bound: lower bound for off-diagonal elements
generate_cor_matrix <- function(n, lower_bound) {
  # Initialize correlation matrix with all elements drawn from a standard normal distribution
  R <- matrix(rnorm(n*n), n, n)
  # Set diagonal elements to 1
  diag(R) <- 1
  # Set off-diagonal elements lower than lower_bound to lower_bound
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (R[i, j] < lower_bound) {
        R[i, j] <- lower_bound
        R[j, i] <- lower_bound
      }
    }
  }
  # Check if R is positive definite
  while (min(eigen(matrix(R), symmetric = T)) < 0) {
    # If not, generate a new R
    R <- matrix(rnorm(n*n), n, n)
    diag(R) <- 1
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        if (R[i, j] < lower_bound) {
          R[i, j] <- lower_bound
          R[j, i] <- lower_bound
        }
      }
    }
  }
  return(R)
}

genCorData <- function(n, mu, sigma, corMatrix = NULL, rho, corstr = "ind",
                       cnames = NULL, idname = "id") {
  nvars <- length(mu)
  
  if (!is.null(cnames)) {
    nnames <- trimws(unlist(strsplit(cnames, split = ",")))
    
    if (length(nnames) != nvars) {
      stop("Invalid number of variable names")
    }
  }
  
  #corMatrix <- .buildCorMat(nvars, corMatrix, corstr, rho)
  
  if (length(sigma) == 1) {
    varMatrix <- (sigma^2) * corMatrix
  } else if (length(sigma) > 0) {
    D <- diag(sigma)
    
    if (length(diag(corMatrix)) != length(sigma)) {
      stop("Improper number of standard deviations")
    }
    
    varMatrix <- (D %*% corMatrix) %*% D
  }
  
  dt <- data.table(mvnfast::rmvn(n = n, mu = mu, sigma = varMatrix))
  
  
  if (!is.null(cnames)) setnames(dt, nnames)
  
  dtid <- data.table(1:nrow(dt))
  setnames(dtid, idname)
  
  dt <- cbind(dtid, dt)
  setkeyv(dt, idname)
  
  return(dt[])
}

.genQuantU <- function(nvars, n, rho, corstr, corMatrix, idname = "id") {
  
  # "Declare" vars to avoid R CMD warning
  seqid <- NULL
  period <- NULL
  Unew <- NULL
  Y <- NULL
  
  mu <- rep(0, nvars)
  if (is.null(corMatrix)) {
    dt <- genCorData(n, mu, sigma = 1, rho = rho, corstr = corstr, idname = idname)
  } else {
    dt <- genCorData(n, mu, sigma = 1, corMatrix = corMatrix, idname = idname)
  }
  
  dtM <- melt(dt, id.vars = idname, variable.factor = TRUE, value.name = "Y", variable.name = "seq")
  
  dtM[, period := as.integer(seq) - 1]
  setkeyv(dtM, idname)
  dtM[, seqid := .I]
  dtM[, Unew := stats::pnorm(Y)]
  
  return(dtM[, -"Y"])
}

m <- rlkjcorr(1, 2)
Ux <- U[,1:2]
m <- cor(Ux)

zs <- .genQuantU(2, 1000, corMatrix = m)
zs[, logisZ := stats::qlogis(p = zs$Unew)]
cprop <- t(apply(baseprobs, 1, cumsum))
quant <- t(apply(cprop, 1, stats::qlogis))

mycat <- matrix(NA, nrow=1000, ncol=2)

for (i in 1:nCats) {
  iLogisZ <- zs[period == i - 1, logisZ]
  n_obs_i <- length(iLogisZ)
  matlp <- matrix(rep(quant[i, ], n),
                  ncol = ncol(cprop),
                  byrow = TRUE
  )
  
  z <- rep(0, n_obs_i)
  
  npVar_mat <- matrix(rep(0, n_obs_i))
  
  #npAdj <- NULL
  
  npmat <- npVar_mat %*% c(0,0,0,0) # npAdj is #npVAR X
  matlp <- matlp - npmat - z
  
  locateGrp <- (iLogisZ > cbind(-Inf, matlp))
  mycat[,i] <- apply(locateGrp, 1, sum)
  
}
 cor(mycat, method = "kendall")

 