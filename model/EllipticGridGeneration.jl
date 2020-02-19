

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




function PoissonSmoothing(ρ, θ, m, n)

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


