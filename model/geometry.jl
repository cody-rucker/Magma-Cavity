#using Pkg
#Pkg.add("Plots")
#Pkg.add("PyPlot")
#include("tfi.jl")
include("TFI.jl")
include("EllipticGridGeneration.jl")
using Plots

r = 1
D = 5
L = 5
ϕ = atan(D/L)

h = 0.1
ξ = 0:h: 1
η = 0:h: 1

n = length(ξ)
m = length(η)



R = sqrt(D^2 + L^2)

ρ₄ = zeros(size(ξ))
θ₄ = zeros(size(ξ))

ρ₂ = zeros(size(ξ))
θ₂ = zeros(size(ξ))

θ₁ = zeros(size(η))
ρ₁ = zeros(size(η))

ρ₃ = zeros(size(η))
θ₃ = zeros(size(η))

for i in 1:length(ξ)
    ρ₄[i] = D/sin((1-ξ[i])*pi/2 + ϕ*ξ[i])
    θ₄[i] = (1-ξ[i])*pi/2 + ϕ*ξ[i]

    ρ₂[i] = r
    θ₂[i] =(1- ξ[i])*pi/2 + ϕ*ξ[i] 
end


for i in 1:length(η)
    ρ₁[i] = D*η[i] + r*(1-η[i])
    θ₁[i] = pi/2
   
 

    ρ₃[i] = R*η[i] + r*(1-η[i])
    θ₃[i] = ϕ
end


#######################################################################
ρ₄ᶜ = zeros(size(ξ))
θ₄ᶜ = zeros(size(ξ))

ρ₂ᶜ = zeros(size(ξ))
θ₂ᶜ = zeros(size(ξ))

θ₁ᶜ = zeros(size(η))
ρ₁ᶜ = zeros(size(η))

ρ₃ᶜ = zeros(size(η))
θ₃ᶜ = zeros(size(η))

for i in 1:length(ξ)
    ρ₄ᶜ[i] = -D/sin(-ξ[i]*pi/2 - ϕ*(1-ξ[i]))
    θ₄ᶜ[i] = (-ξ[i])*pi/2 - ϕ*(1-ξ[i])

  #  ρ₂ᶜ[i] = r
  #  θ₂ᶜ[i] =(- ξ[i])*pi/2 - ϕ*(1-ξ[i]) 
end


for i in 1:length(η)
    ρ₁ᶜ[i] = R*η[i] + r*(1-η[i])
    θ₁ᶜ[i] = -ϕ
    
    

    ρ₃ᶜ[i] = D*η[i] + r*(1-η[i])
    θ₃ᶜ[i] = -pi/2
end

#####################################################################################################


ρ₄ᴮ = zeros(size(ξ))
θ₄ᴮ = zeros(size(ξ))

ρ₂ᴮ = zeros(size(ξ))
θ₂ᴮ = zeros(size(ξ))

θ₁ᴮ = zeros(size(η))
ρ₁ᴮ = zeros(size(η))

ρ₃ᴮ = zeros(size(η))
θ₃ᴮ = zeros(size(η))

for i in 1:length(ξ)
    ρ₄ᴮ[i] = L/cos(ϕ*(1-2*ξ[i]))
    θ₄ᴮ[i] = ϕ*(1-2*ξ[i])

#    ρ₂ᴮ[i] = r
#    θ₂ᴮ[i] =ϕ*(1-2*ξ[i]) 
end


for i in 1:length(η)
    ρ₁ᴮ[i] = R*η[i] + r*(1-η[i])
    θ₁ᴮ[i] = ϕ
    
    

    ρ₃ᴮ[i] = R*η[i] + r*(1-η[i])
    θ₃ᴮ[i] = -ϕ
end



######################################################################################################
ρᴬ = zeros(m, n)
θᴬ = zeros(m, n)

ρᴬ = TFI(ρ₂, ρ₄, ρ₁, ρ₃, ξ, η)
θᴬ = TFI(θ₂, θ₄, θ₁, θ₃, ξ, η)


######################################################################################################

ρᶜ = zeros(length(η),length(ξ))
θᶜ = zeros(length(η),length(ξ))

ρᶜ = TFI(ρ₃, ρ₄ᶜ, ρ₁ᶜ, ρ₃ᶜ, ξ, η) 
θᶜ = TFI(θ₃, θ₄ᶜ, θ₁ᶜ, θ₃ᶜ, ξ, η)
#ρₙ, θₙ = tfi(ρ₃, θ₃, ρ₄ᶜ, θ₄ᶜ, ρ₁ᶜ, θ₁ᶜ, ρ₃ᶜ, θ₃ᶜ)

##################################################################################

ρᴮ = zeros(length(η),length(ξ))
θᴮ = zeros(length(η),length(ξ)) 

ρᴮ = TFI(ρ₃ᶜ, ρ₄ᴮ, ρ₁ᴮ, ρ₃ᴮ,ξ, η)
θᴮ = TFI(θ₃ᶜ, θ₄ᴮ, θ₁ᴮ, θ₃ᴮ, ξ, η)
#ρₙ, θₙ = tfi(ρ₃, θ₃, ρ₄ᶜ, θ₄ᶜ, ρ₁ᶜ, θ₁ᶜ, ρ₃ᶜ, θ₃ᶜ)
###################################################################



for t = 1:1000

    # Conditions for block A, upper block

    ρᴬ[:, :], θᴬ[:, :] = LaplaceSmoothing(ρᴬ, θᴬ, m, n)

    ρᴬ[2:n-1, 1] = ρᴬ[2:n-1, 2]
    ρᴬ[m-1, :] = ρᴬ[m, :]
    θᴬ[m-1, :] = θᴬ[m,:]
    θᴬ[2, :] = θᴬ[1, :]

    # COnditions for block B, middle block.

    ρᴮ[:, :], θᴮ[:, :] = LaplaceSmoothing(ρᴮ, θᴮ, m, n)

    θᴮ[2, :] = θᴮ[1, :]

    ρᴮ[m-1, :] = ρᴮ[m, :]
    θᴮ[m-1, :] = θᴮ[m,:]


    # conditions for block C, lower block

    ρᶜ[:, :], θᶜ[:, :] = LaplaceSmoothing(ρᶜ, θᶜ, m, n)
    ρᶜ[2:n-1, m] = ρᶜ[2:n-1, m-1]

    ρᶜ[m-1, :] = ρᶜ[m, :]
    θᶜ[m-1, :] = θᶜ[m, :]

    θᶜ[1, :] = θᶜ[2, :]




#pyplot()
    pl = plot( proj=:polar, legend= :false, showaxis= :false, m=2, grid=:true)
    plot!(size=(700,700))

    
    plot!(θᴬ, ρᴬ, proj=:polar,linecolor= :black, lw=1)
    plot!(θᴬ', ρᴬ', proj=:polar,linecolor= :black, lw=1)
    plot!(θᴮ, ρᴮ, proj=:polar,linecolor= :black, lw=1)
    plot!(θᴮ', ρᴮ', proj=:polar,linecolor= :black, lw=1)
    plot!(θᶜ, ρᶜ, proj=:polar,linecolor= :black, lw=1)
    plot!(θᶜ', ρᶜ', proj=:polar,linecolor= :black, lw=1, xlim=[0,1])

    display(pl)


end

#savefig("tfi_elliptic_grid.png")











