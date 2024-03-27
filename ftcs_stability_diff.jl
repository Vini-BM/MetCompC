# Vinícius Müller
# 26/03/2024

# Testing the stability of FTCS method
# Diffusion equation
# Source at x=0: \rho(0,t)=1
# Drain at x=L: \rho(L,t)=0
# Initial condition: \rho(x>0,t=0)=0

using Plots
Plots.default(show = true)
gr() # backend for Plots

# Method
function ftcs(f,tf,k,dt=0.1)
    t = 0 # initializes time
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
tf = 20 # time interval
k_list = [.3, .5, .7, .9] # parameter to be tested
plot_list = []
ρ0 = 1 # initial condition

# Initialization
ρ0 = zeros(L+1) # initializes \rho as zeros
ρ0[1] = 1 # sets initial condition [1, 0, ... , 0]
println(ρ0)

# Integration
for i=1:length(k_list)
    k = k_list[i]
    ρ = ftcs(ρ0,tf,k)
    println(ρ)
    push!(plot_list,plot(x,ρ,title="k=$k",legend=false,color=i,dpi=1000))
end

# Plot
plot!(plot_list...,plot_title="FTCS stability | t=$tf")
xlabel!("x")
ylabel!("ρ")
readline() # keep plot open
savefig("Plots/ftcs_stability.png")
