const error_notice_flag = -2.

"""
Fixed-point iteration solver.

    * func works for arbitrary dimensional input.
    * The metric function needs to be supplied and it determines the error 
      between iterations.
    * x₀ is the input guess. It also determines the output dimension.
    * args are variadic input parameters passed to func
    * The tolerance and max iterations may need tuning.
"""
function FixedPointIteration( func, metric, x₀, args...; tolerance=1e-10, maxiter=100 )
    old_point = copy(x₀) 
    new_point = zeros( size(x₀) )
    error = 1.
    all_errors = error_notice_flag * ones(maxiter+1)
    all_errors[1] = error
    iteration = 1
    error_decreasing = true
    while error_decreasing && error >= tolerance && iteration < maxiter
        iteration += 1
        new_point = func(old_point, args...)
        error = metric( new_point, old_point )
        old_point = copy(new_point)
        # @show(iteration, error)
        all_errors[iteration] = error
        # error_decreasing = all_errors[iteration] <= all_errors[iteration-1]
    end
    display(all_errors)
    println("$iteration iterations completed with an error of $error.")
    return new_point, all_errors
end