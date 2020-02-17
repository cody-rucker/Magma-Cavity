# given parametrized curves c₁, c₂, c₃, c₄ where c₂, c₄ parametrize the bottom and top curves
# while c₁, c₃ parametrize left and right curves

#=
                              c₄ = (ρ(ξ), θ(ξ))
                         ————————————————————————————–
                        |                            |
                        |                            |
                        |                            |
     c₁ = (ρ(η), θ(η))  |                            |  c₃ = (ρ(η), θ(η))
                        |                            |
                        |                            |
                        |                            |
                        ——————————————————————————————
                              c₂ = (ρ(ξ), θ(ξ))
=#

function TFI(x₂, x₄, x₁, x₃, ξ, η)

    m, n = length(η), length(ξ)
 

    x = zeros(length(η),length(ξ))



    for i = 1:n
        for j = 1:m
            x[j,i] = (1-η[j])*x₂[i] + η[j]*x₄[i] + (1-ξ[i])*x₁[j] + ξ[i]*x₃[j]-
            (ξ[i]*η[j]*x₄[length(ξ)] + ξ[i]*(1-η[j])*x₂[length(ξ)]
               + η[j]*(1-ξ[i])*x₄[1] + (1-ξ[i])*(1-η[j])*x₂[1])
        end
    end

    return x
end




function LaplaceSmoothing(ρ, θ, m, n)

    α = zeros(m, n)
    β = zeros(m, n)
    γ = zeros(m, n)


    newρ = copy(ρ)
    newθ = copy(θ)

    for i = 2:m-1
        for j = 2:n-1

            α[i, j] = (1/4)*((ρ[i,j+1]-ρ[i,j-1])^2 + (θ[i,j+1] - θ[i,j-1])^2)

            β[i, j] =  (1/16)*((ρ[i+1,j]-ρ[i-1,j])*(ρ[i,j+1]-ρ[i,j-1])
                               +(θ[i+1,j] - θ[i-1,j])*(θ[i,j+1] - θ[i,j-1]))

            γ[i, j] = (1/4)*((ρ[i+1,j]-ρ[i-1,j])^2 + (θ[i+1,j] - θ[i-1,j])^2)
        end
    end

    for i = 2:m-1
        for j = 2:n-1

            newρ[i,j]=((-0.5)/(α[i,j]+γ[i,j]+10^-9))*(2*β[i,j]*(ρ[i+1,j+1]-ρ[i-1,j+1]-ρ[i+1,j-1] + ρ[i-1,j-1])
                                                      -α[i,j]*(ρ[i+1,j]+ρ[i-1,j])-γ[i,j]*(ρ[i,j+1]+ρ[i,j-1]))

            newθ[i,j]=((-0.5)/(α[i,j]+γ[i,j]+10^-9))*(2*β[i,j]*(θ[i+1,j+1]-θ[i-1,j+1]-θ[i+1,j-1] + θ[i-1,j-1])
                                                      -α[i,j]*(θ[i+1,j]+θ[i-1,j])-γ[i,j]*(θ[i,j+1]+θ[i,j-1]))
        end
    end
    return newρ, newθ

end

