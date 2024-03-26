# Vinícius Müller
# 25/03/2024

# Solving the diffusion equation with FTCS method
# Source at x=0: \rho(0,t)=1
# Drain at x=L: \rho(L,t)=0
# Initial condition: \rho(x>0,t=0)=0

using Plots
Plots.default(show = true)
gr() # backend for Plots

# Method
function ftcs(f,ti,tf,dt=0.1)
    t = ti # initializes time
    while t < tf # loop
        t+=dt # increase time
        f[2:L-1] += k*(f[3:L] + f[1:L-2] - 2f[2:L-1]) # ftcs discretization
        # Extreme values are treated as boundary conditions
    end
    return f
end

# Parameters
L = 100 # xmax
x = [i for i in 0:L] # space interval
println(x)
tf_list = [0, 25, 50, 75, 100] # time interval
k = 0.5 # k = D dt/dx^2 (k<=0.5 for FTCS), D: diffusion coefficient
ρ0 = 1 # initial condition

# Initialization
ρ = zeros(L+1) # initializes \rho as zeros
ρ[1] = ρ0 # sets initial condition [1, 0, ... , 0]
println(ρ)
plot!(x,ρ,label="t=0",dpi=1000)

# Integration
for i in 2:lastindex(tf_list)
    tf = tf_list[i]
    ti = tf_list[i-1]
    global ρ = ftcs(ρ,ti,tf)
    println(ρ)
    plot!(x,ρ,label="t=$tf",dpi=1000)
end

# Plot
plot!(legend=:bottomright)
title!("Diffusion")
xlabel!("x")
ylabel!("ρ")
readline() # keep plot open
savefig("Plots/ftcs_diffusion.png")