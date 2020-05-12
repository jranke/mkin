mkinpredict <-
function (x, odeparms, odeini, outtimes = seq(0, 120, by = 0.10000000000000001), 
    solution_type = "deSolve", use_compiled = "auto", method.ode = "lsoda", 
    atol = 1e-08, rtol = 1e-10, map_output = TRUE, ...) 
{
    UseMethod("mkinpredict", x)
}
