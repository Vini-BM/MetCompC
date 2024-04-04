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
function ftcsImplicit(f,IA,k,ti,tf,bc0=1,bcL=0,dt=0.1)
    t = ti # initializes time
    while t < tf # loop
        t += dt # increase time
        f[1] += k*bc0 # boundary condition
        f[end] += k*bcL # boundary condition
        f = IA*f # ftcs discretization
        # Extreme values are treated as boundary conditions
        # It is necessary to add k*BC to 2nd and L-1th elements
    end
    return f
end

# Parameters
L = 100 # xmax
k = 0.5 # k = D dt/dx^2 (k<=0.5 for FTCS), D: diffusion coefficient
bc0 = 1 # boundary condition
bcL = 0 # boundary condition


# Initialization

## space, time, f
x = [i for i in 1:L] # space interval
tf_list = [25*i for i in 0:6] # time interval
f = zeros(L-2) # initializes f as zeros from 2nd element to L-1th element --> length(x) - 2

## Matrix
subd = [-k for i in 1:L-3] # subdiagonal = upperdiagonal --> L-3 elements
diag = [1+2k for i in 1:L-2] # diagonal --> L-2 elements
A = Tridiagonal(subd,diag,subd) # matrix
IA = inv(A) # inverse matrix

## Plot
plot(x,vcat(bc0,f,0),label="t=0",dpi=1000) # plot initial condition


# Integration
for i in 2:lastindex(tf_list)
    tf = tf_list[i]
    ti = tf_list[i-1]
    global f = ftcsImplicit(f,IA,k,ti,tf)
    full_f = vcat(bc0,f,bcL) # add BC to f
    plot!(x,full_f,dpi=1000,label="t=$tf")
end


# Plot
plot!(legend=:bottomright)
title!("Diffusion (Implicit)")
xlabel!("x")
ylabel!("f")
readline() # keep plot open
savefig("Plots/ftcsDiffusionImplicit.png")
