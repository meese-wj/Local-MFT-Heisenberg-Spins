struct Spin3
    S₁::Float64
    S₂::Float64
    S₃::Float64
end

import Base.+
import Base.*

+(A::Spin3, B::Spin3) = Spin3( A.S₁ + B.S₁, A.S₂ + B.S₂, A.S₃ + B.S₃ )
*(A::Spin3, λ::Real) = Spin3( λ * A.S₁, λ * A.S₂, λ * A.S₃ )
*(λ::Real, A::Spin3) = A * λ
-(A::Spin3, B::Spin3) = A + (-1. * B)
⋅(A::Spin3, B::Spin3) = A.S₁ * B.S₁ + A.S₂ * B.S₂ + A.S₃ * B.S₃
abs2(A::Spin3) = A ⋅ A
abs(A::Spin3) = sqrt( abs2(A) )
unit_spin3(A::Spin3) = A * 1./abs(A)