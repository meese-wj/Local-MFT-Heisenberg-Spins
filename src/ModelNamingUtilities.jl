"""
Determine whether the model is antiferromagnetic (AFM)
or ferromagnetic (FM).
"""
function model_name( coupling::Float64, latt_params::LatticeParameters ) 
    name = ""
    if coupling > 0.
        name = "AFM"
    elseif coupling < 0.
        name = "FM"
    end
    return "$(name)_Lx-$(latt_params.Lx)_Ly-$(latt_params.Ly)"
end      