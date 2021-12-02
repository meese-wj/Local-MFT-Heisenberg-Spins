const error_notice_flag = -2.

"""
Fixed-point iteration solver.

    * func works for arbitrary dimensional input.
    * The metric function needs to be supplied and it determines the error 
      between iterations.
    * x₀ is the input guess. It also determines the output dimension.
    * args are variadic input parameters passed to func
    * The tolerance and max iterations may need tuning.
    * The state_function provides a way for measuring a scalar state
      during the fixed-point procedure. If the state_function takes in 
      parameters, then include it as a λ-function.
"""
function FixedPointIteration( func, metric, x₀, args...; 
                              tolerance=1e-10, maxiter=100, state_function=nothing )
    old_point = copy(x₀) 
    new_point = zeros( size(x₀) )
    error = 1.
    all_errors = error_notice_flag * ones(maxiter)
    all_errors[1] = error

    all_states = nothing
    if state_function !== nothing
        all_states = zeros( maxiter )
        all_states[1] = state_function( x₀ )
    end
    
    iteration = 1
    error_decreasing = true
    while error_decreasing && error >= tolerance && iteration < maxiter
        iteration += 1
        new_point = func(old_point, args...)
        error = metric( new_point, old_point )
        old_point = copy(new_point)
        all_errors[iteration] = error
        # error_decreasing = all_errors[iteration] <= all_errors[iteration-1]
        if state_function !== nothing 
            all_states[iteration] = state_function( new_point )
        end
    end
    display(all_errors)
    println("$iteration iterations completed with an error of $error.")
    if state_function === nothing 
        return new_point, all_errors, nothing
    else
        return new_point, all_errors, all_states
    end
end