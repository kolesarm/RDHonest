## Tolerance for comparing results that pass through numerical optimization
## (optimal bandwidths, and anything computed at them). Near an optimum the
## objective is flat, so rounding differences of order .Machine$double.eps in
## the objective move the optimizer's answer by about sqrt(.Machine$double.eps),
## roughly 1.5e-8 relative. Such differences arise between mathematically
## equivalent computations, and between BLAS/LAPACK implementations (e.g. Intel
## MKL). This leaves a wide margin above that while still catching any real
## error.
tol_opt <- 1e-6
