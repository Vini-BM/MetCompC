# Vinícius Müller
# 04/04/2024

# Solving the diffusion equation with Crank-Nicholson method
# Radiation defects on material: f_t = Df_{xx} - Sf + s(x)
# Defect source: s(x) = s_0 exp[-(x-L/2)²/2\sigma] -> gaussian
# Drain S
# Initial condition: f(0,x) = 0
# Boundary conditions: f(t,0) = f(t,L) = 0

using LinearAlgebra
using Plots
Plots.default(show = true)

# Parameters
L = 100 # xmax
dt = .1 # timestep
tf = 50 # integration time
S = .2 # amplitude of drain
k = .5

# Initialization
x = collect(1:1:L) # x array
f = zeros(L) # initial condition
plot(x,f,label="t=0",dpi=1000) # plot initial condition

# Gaussian source
s0 = 2 # amplitude
sigma = 2 # variance
center = L/2 # center of gaussian
s = zeros(L) # initialize array
@. s = s0*exp(-(x-center)^2 / 2sigma^2) # broadcasting

# Matrix

## n+1
diagPlus = [2+2k+S*dt for i in 1:L]
updiagPlus = [-k for i in 1:L-1]
A = Tridiagonal(updiagPlus,diagPlus,updiagPlus)
A[1,L] = -k
A[L,1] = -k
IA = inv(A)

## n
diag = [2-2k-S*dt for i in 1:L]
updiag = [k for i in 1:L-1]
B = Tridiagonal(updiag,diag,updiag)
B[1,L] = k
B[L,1] = k

# Method
function cn(f,IA,B,source,ti,tf,dt)
    t = ti
    while t < tf
        t += dt
        f = IA*(B*f+source*dt)
    end
    return f
end

# Integration
f = cn(f,IA,B,s,0,tf,dt)

# Plot
plot!(x,f,label="t=$tf",dpi=1000)
title!("Radiation Defects (Crank-Nicholson)")
xlabel!("x")
ylabel!("f")
readline() # keep plot open
#savefig("Plots/cnRadDefects.png")
