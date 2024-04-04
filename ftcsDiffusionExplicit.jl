# Vinícius Müller
# 25/03/2024

# Solving the diffusion equation with explicit FTCS method
# Source at x=0: f(0,t)=1
# Drain at x=L: f(L,t)=0
# Initial condition: f(x>0,t=0)=0

using Plots
Plots.default(show = true)
gr() # backend for Plots

# Method
function ftcs(f,ti,tf,dt=0.1)
    t = ti # initializes time
    while t < tf # loop
        t+=dt # increase time
        f[2:end-1] += k*(f[3:end] + f[1:end-2] - 2f[2:end-1]) # ftcs discretization
        # Extreme values are treated as boundary conditions
    end
    #return f -> not necessary because f is global
end

# Parameters
L = 100 # xmax
x = [i for i in 1:L] # space interval
tf_list = [25*i for i in 0:6] # time interval
k = 0.5 # k = D dt/dx^2 (k<=0.5 for FTCS), D: diffusion coefficient
f0 = 1 # initial condition

# Initialization
f = zeros(L) # initializes f as zeros
f[1] = f0 # sets initial condition [1, 0, ... , 0]
plot!(x,f,label="t=0",dpi=1000)

# Integration
for i in 2:lastindex(tf_list)
    tf = tf_list[i]
    ti = tf_list[i-1]
    ftcs(f,ti,tf)
    plot!(f,label="t=$tf",dpi=1000)
end

# Plot
plot!(legend=:bottomright)
title!("Diffusion (Explicit)")
xlabel!("x")
ylabel!("f")
readline() # keep plot open
savefig("Plots/ftcsDiffusionExplicit.png")
