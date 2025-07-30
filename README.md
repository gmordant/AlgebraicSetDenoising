
# Denoising for latent variables on algebraic sets

This project addresses the issue of denoising when the latent variables are supported on an algebraic set.





## Functions

The main functions included in the script is:


### `multiexponentsL(dim,deg)`

Generate the exponents for the polynomials of degree at most deg in dimension dim.

### `Vandermonde(X, U, sigma)`

Construct the debiased moment matrix for the given parameters. 
X is the dataset, U is the ouput of multiexponents and sigma is the denoising level required. 



## Example

An example is provided within the script that involves generating data points from a circle with added Gaussian noise, computing the relevant matrices, and solving the quadratics to visualize the circle.

### Example Code Usage

```R
set.seed(123)

# Number of points
n <- 600
dim<- 2
deg <- 2

# Radius of the circle
radius <- sqrt(2)
sigma<-0.4

# Generate angles evenly spaced around the circle
angles <- seq(0, 2*pi, length.out = n)
x1 <- radius * cos(angles)
y1 <- radius * sin(angles)

x2 <- x1 + sigma * rnorm(n)
y2 <- y1 + sigma * rnorm(n)

# Create data matrix, compute moment matrix, and extract kernel
X <- cbind(x2, y2)

# Without correction
M <- Vandermonde(X, U=multiexponentsL(dim, deg), sigma=0)
temp <- svd(M)
q <- temp$v[,6]
q <- q / q[1]

S <- SolveOrder2(q)

# With correction
svd_result <- svd(Vandermonde(X, U=multiexponentsL(dim, deg), sigma=sigma))
q1 <- svd_result$v[,6]
q1 <- q1 / q1[1]

S1 <- SolveOrder2(q1)

SA <- SemiAlg(q1, lambda=n^{-1/2}*log(n))

idx <- sample(1:nrow(SA), 2000)

# Plotting the results
plot(x2, y2, col = "black", pch = 19, cex=0.4, xlab = "", ylab = "", xlim= c(-2.5,2.5), ylim = c(-2.5,2.5))
points(x1, y1, col = "green", pch = 19, cex=0.4)
points(S[,1], S[,2], col = "red", pch = 19, cex = 0.4)
points(S1[,1], S1[,2], col = "blue", pch = 19, cex = 0.4)
```

## License

This project is licensed under the MIT License.

