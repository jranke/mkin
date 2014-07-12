IORE.solution <- function(t, parent.0, k.iore, N)
{
	parent = (parent.0^(1 - N) - (1 - N) * k.iore * t)^(1/(1 - N))
}
