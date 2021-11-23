"""
Compute the Langevin function as 
    coth βx - 1/x 
but be careful about small x. Right now 
this is accurate up to order O(beta¹¹x⁹).
"""
function LangevinFunction( x, β )
    arg = β * x
    if arg < 0.01
        output = (((-1 * arg^7. / 4725) + 2 * arg^5. / 945) - arg^3. / 45) + arg / 3
        return β * arg
    end
    return coth(arg) - 1/x
end