# Vinícius Müller
# 26/03/2024

# Testing the stability of explicit FTCS method
# Diffusion equation
# Source at x=0: f(0,t)=1
# Drain at x=L: f(L,t)=0
# Initial condition: f(x>0,t=0)=0

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

# Initialization
f0 = zeros(L+1) # initializes f as zeros
f0[1] = 1 # sets initial condition [1, 0, ... , 0]

# Integration
for i=1:length(k_list)
    k = k_list[i]
    f = ftcs(f0,tf,k)
    push!(plot_list,plot(x,f,title="k=$k",legend=false,color=i,dpi=1000))
end

# Plot
plot!(plot_list...,plot_title="Diffusion (Explicit FTCS stability) | t=$tf")
xlabel!("x")
ylabel!("f")
readline() # keep plot open
savefig("Plots/ftcsStability.png")
