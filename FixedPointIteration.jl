"""
Fixed-point iteration solver.

    * func works for arbitrary dimensional input.
    * The metric function needs to be supplied and it determines the error 
      between iterations.
    * x₀ is the input guess. It also determines the output dimension.
    * args are variadic input parameters passed to func
    * The tolerance and max iterations may need tuning.
"""
function FixedPointIteration( func, metric, x₀, args...; tolerance=1e-8, maxiter=5000 )
    old_point = copy(x₀) 
    new_point = zeros( size(x₀) )
    error = 1.
    iteration = 0
    while error >= tolerance && iteration < maxiter
        iteration += 1
        new_point = func(old_point, args...)
        error = metric( new_point, old_point )
        old_point = copy(new_point)
        # @show(iteration, error)
    end
    println("$iteration iterations completed with an error of $error.")
    return new_point
end