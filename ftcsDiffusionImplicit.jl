# Vinícius Müller
# 28/03/2024

# Solving the diffusion equation with implicit FTCS method (backward difference)
# Source at x=0: f(0,t)=1
# Drain at x=L: f(L,t)=0
# Initial condition: f(x>0,t=0)=0

using LinearAlgebra
using Plots
Plots.default(show = true)
gr() # backend for Plots

# Method
function ftcsImplicit(f,IA,ti,tf,bc0=1,bcL=0,dt=0.1)
    t = ti # initializes time
    while t < tf # loop
        t+=dt # increase time
        f = IA*f # ftcs discretization
        f[1] = bc0 # boundary condition
        f[end] = bcL # boundary condition
        # Extreme values are treated as boundary conditions
    end
    return f #-> not necessary because f is global
end

# Parameters
L = 100 # xmax
x = [i for i in 0:L] # space interval
tf_list = [25*i for i in 0:6] # time interval
k = 0.5 # k = D dt/dx^2 (k<=0.5 for FTCS), D: diffusion coefficient
f0 = 1 # initial condition

# Initialization
subd = [-k for i in 0:L-1]
diag = [1+2k for i in 0:L]
A = Tridiagonal(subd,diag,subd) # matrix
IA = inv(A) # inverse matrix
f = zeros(L+1) # initializes f as zeros

f[1] = f0 # sets initial condition [1, 0, ... , 0]
plot(x,f,label="t=0",dpi=1000) # plot initial condition

# Integration
for i in 2:lastindex(tf_list)
    tf = tf_list[i]
    ti = tf_list[i-1]
    global f = ftcsImplicit(f,IA,ti,tf) # didn't work without method returning f
    plot!(x,f,dpi=1000,label="t=$tf")
end


# Plot
plot!(legend=:bottomright)
title!("Diffusion (Implicit)")
xlabel!("x")
ylabel!("f")
readline() # keep plot open
savefig("Plots/ftcsDiffusionImplicit.png")
