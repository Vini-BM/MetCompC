# Vinícius Müller
# 09/04/2024

using Plots
using LinearAlgebra
using Distributions
Plots.default(show=true)

# Method
function lax(f,k,ti,tf,dt)
    # k = vdt/2dx
    t = ti
    while t<tf
        t += dt
        # Inner loop
        f[2:end-1] = .5(f[3:end] + f[1:end-2]) - k*(f[3:end]-f[1:end-2])
        # Periodic Boundary Conditions
        bc1 = f[1]
        bcL = f[end]
        f[1] = .5(f[2]+bcL) - k*(f[2]-bcL)
        f[end] = .5(bc1+f[end-1]) - k*(bc1-f[end-1])
    end
    return f
end

# Parameters
L = 100 # xmax
x = collect(1:1:L) # space interval
klist = [.45,.5,.55]
ti = 0
tf = 5
dt = 0.1

# Initialization

f = zeros(L) # initialize f as zeros
sigma = rand(Uniform(0,L/4)) # standard deviation
@. f = exp(-(x-L/2)^2 / 2sigma^2) # gaussian package
println("Standard deviation = $sigma")
plot(x,f,label="t=$ti",ls=:dash)

# Loop
for k in klist
    g = lax(f,k,ti,tf,dt)
    println("ok")
    plot!(x,g,label="k=$k")
    println("plot ok")
end

# Plot
println("final")
title!("Drift Equation (Lax Method) | t=$tf")
xlabel!("x")
ylabel!("f")
readline() # keep plot open
