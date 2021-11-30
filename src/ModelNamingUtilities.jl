include("HeisenbergModel.jl")
include("Heisenberg_J1-J2_Model.jl")

"""
Determine whether the model is antiferromagnetic (AFM)
or ferromagnetic (FM).
"""
function model_name( model_params::Any )
    return
end

"""
Name for the Heisenberg model with nearest neighbor 
interactions only.
"""
function model_name( model_params::ModelParameters )
    name = ""
    if model_params.Jex > 0.
        name = "AFM"
    elseif model_params.Jex < 0.
        name = "FM"
    end
    return name
end

"""
Name for the Heisenberg J1-J2 model.
"""
function model_name( model_params::J1_J2_ModelParameters )
    name = ""
    if model_params.J2_ex > 0.
        name = "AFM-J2"
    elseif model_params.J2_ex < 0.
        name = "FM-J2"
    end
    return name
end

"""
Combine the model_name and lattice parameters 
into a final string name.
"""
function model_name( model_params, latt_params::LatticeParameters ) 
    return "$(model_name(model_params))_Lx-$(latt_params.Lx)_Ly-$(latt_params.Ly)"
end      