struct Spin3
    S₁::Float64
    S₂::Float64
    S₃::Float64
end

using Base

Base.:+(A::Spin3, B::Spin3) = Spin3( A.S₁ + B.S₁, A.S₂ + B.S₂, A.S₃ + B.S₃ )
Base.:*(A::Spin3, λ::Real) = Spin3( λ * A.S₁, λ * A.S₂, λ * A.S₃ )
Base.:*(λ::Real, A::Spin3) = A * λ
Base.:-(A::Spin3, B::Spin3) = A + (-1. * B)
⋅(A::Spin3, B::Spin3) = A.S₁ * B.S₁ + A.S₂ * B.S₂ + A.S₃ * B.S₃
Base.abs2(A::Spin3) = A ⋅ A
Base.abs(A::Spin3) = sqrt( abs2(A) )
Base.copy(A::Spin3) = Spin3( A.S₁, A.S₂, A.S₃ )
unit_spin3(A::Spin3) = A * (1. / abs(A))

"""
Calculate the projection of 𝐀 in the direction of 𝐁
"""
proj( A::Spin3, B::Spin3 ) = ( ( A ⋅ B ) / abs2( B ) ) * B 
# Projection in the x-direction
projx( A::Spin3 ) = Spin3( A.S₁, 0., 0. )
# Projection in the y-direction
projy( A::Spin3 ) = Spin3( 0., A.S₂, 0. )
# Projection in the z-direction
projz( A::Spin3 ) = Spin3( 0., 0., A.S₃ )