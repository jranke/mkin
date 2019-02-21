logistic.solution <- function(t, parent.0, kmax, k0, r)
{
	parent = parent.0 * (kmax / (kmax - k0 + k0 * exp (r * t))) ^(kmax/r)
}
