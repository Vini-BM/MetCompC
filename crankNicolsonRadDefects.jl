# Vinícius Müller
# 06/04/2024

# Solving the diffusion equation with Crank-Nicolson method
# Radiation defects on material: f_t = Df_{xx} - Sf + s(x)
# Defect source: s(x) = s_0 exp[-(x-L/2)²/2\sigma] -> gaussian
# Drain S
# Initial condition: f(0,x) = 0
# Periodic boundary conditions: f(t,0) = f(t,L) -> enables stationary state

using LinearAlgebra
using Plots
using Distributions
using Printf
Plots.default(show = true)

# Parameters
xmax = 50 # xmax
dt = .1 # timestep
S = .7 # amplitude of drain
k = .5 # D dt / (dx)^2
tf_list = [.5*i for i in 0:10] # time interval

# Initialization
dx = .1
x = collect(1:.1:xmax) # x array
L = length(x) # for resolution
f = zeros(L) # initial condition
plot(x,f,label=@sprintf("t = %.2f", 0),dpi=1000) # plot initial condition

# Gaussian source
s0 = 2 # amplitude
sigma = rand(Uniform(0,xmax/3)) # variance
center = xmax/3 # center of gaussian
s = zeros(L) # initialize array
@. s = s0*exp(-(x-center)^2 / 2sigma^2) # broadcasting
println("Mean = $center")
println("Variance = $sigma")
println("Amplitude = $s0")

# Matrix

## n+1
diagPlus = [2+2k+S*dt for i in 1:L] # matrix diagonal
updiagPlus = [-k for i in 1:L-1] # upper diagonal = lower diagonal
A = Tridiagonal(updiagPlus,diagPlus,updiagPlus)
A = convert(Matrix,A) # change elements
A[1,L] = -k # PBC
A[L,1] = -k # PBC
IA = inv(A)

## n
diag = [2-2k-S*dt for i in 1:L]
updiag = [k for i in 1:L-1]
B = Tridiagonal(updiag,diag,updiag)
B = convert(Matrix,B)
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
for i in 2:length(tf_list)
    ti = tf_list[i-1]
    tf = tf_list[i]
    label = @sprintf("t = %.2f",tf)
    global f = cn(f,IA,B,s,ti,tf,dt)
    plot!(x,f,label=label,dpi=1000)
end

# Plot
title!("Radiation Defects (Crank-Nicolson)")
xlabel!("x")
ylabel!("f")
readline() # keep plot open
savefig("Plots/cnRadDefects.png")
