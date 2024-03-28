# Vinícius Müller
# 28/03/2024

# Solving the diffusion equation with implicit FTCS method (backward difference)
# Source at x=0: \rho(0,t)=1
# Drain at x=L: \rho(L,t)=0
# Initial condition: \rho(x>0,t=0)=0

using LinearAlgebra
using Plots
Plots.default(show = true)
gr() # backend for Plots

# Method
function ftcsImplicit(f,IA,ti,tf,dt=0.1)
    t = ti # initializes time
    while t < tf # loop
        t+=dt # increase time
        f*=IA # ftcs discretization
        # Extreme values are treated as boundary conditions
    end
    #return f -> not necessary because f is global
end

# Parameters
L = 100 # xmax
x = [i for i in 0:L] # space interval
tf_list = [0, 25, 50, 75, 100] # time interval
k = 0.5 # k = D dt/dx^2 (k<=0.5 for FTCS), D: diffusion coefficient
f0 = 1 # initial condition

# Initialization
subd = [-k for i in 0:L-1]
diag = [1+2k for i in 0:L]
A = Tridiagonal(subd,diag,subd) # matrix
IA = inv(A) # inverse matrix
f = zeros(L+1) # initializes \rho as zeros
f[1] = f0 # sets initial condition [1, 0, ... , 0]
plot(x,f,label="t=0",dpi=1000)


# Integration
#for i in 2:lastindex(tf_list)
#    tf = tf_list[i]
#    ti = tf_list[i-1]
#    #@time ftcs(ρ,ti,tf)
#    plot!(ρ,label="t=$tf",dpi=1000)
#end
ftcsImplicit(f,IA,0,500)
plot(x,f,dpi=1000)

# Plot
plot!(legend=:bottomright)
title!("Diffusion")
xlabel!("x")
ylabel!("f")
readline() # keep plot open
#savefig("Plots/ftcs_diffusion.png")
