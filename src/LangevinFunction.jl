"""
Compute the Langevin function as 
    coth βx - 1/x 
but be careful about small x. Right now 
this is accurate up to order O(β⁹x⁹).
"""
function LangevinFunction( x, β )
    arg = β * x
    if arg < 0.25
        output = (((-1 * arg^7. / 4725) + 2 * arg^5. / 945) - arg^3. / 45) + arg / 3
        return arg
    end
    return coth(arg) - 1/arg
end

"""
Compute the Langevin function at zero
temperature β = infinity
"""
function LangevinFunction( x )
    return 1.
end