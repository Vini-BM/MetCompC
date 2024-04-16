# Vinícius Müller
# 10/04/2024

using Plots
#susing PlotlyJS
#using LinearAlgebra
Plots.default(show=true)


# Method: Simple Over-Relaxation
function overRelax(f,alpha,tol,bcL,bcR,lambda1,lambda2)
    L = size(f)[1] # length of x array
    x = collect(1:1:L) # x array
    t = 0 # timestep
    diff = 0 # loop control
    while abs(diff) < tol # condition based on tolerance
    #while t < 10
        lastf = deepcopy(f) # copy of f
        t += 1
        f[2:end-1,2:end-1] = -alpha*lastf[2:end-1,2:end-1] + (1+alpha)/4 * (lastf[3:end,2:end-1] + f[1:end-2,2:end-1] + lastf[2:end-1,3:end] + f[2:end-1,1:end-2]) # update f
        @. f[L,1:end] = bcR # right boundary condition
        @. f[1,1:end] = bcL # left boundary condition
        @. f[1:end,1] = exp(-lambda1*x) # bottom boundary condition
        @. f[1:end,end] = exp(-lambda2*x) # top boundary condition
        diff = maximum(abs.(f)) - maximum(abs.(lastf))
        println(diff)
    end
    println(diff)
    return f, x, t
end

# Initialization

L = 50
l1 = .1
l2 = .3
bcL = 10
bcR = 0
tol = .001
alpha = .4
#f = rand(Float64,(L,L))
f = ones(L,L)

# Loop
g, x, itertime = overRelax(f,alpha,tol,bcL,bcR,l1,l2)
println("Iteration time = $itertime")
#println(g)

# Plot
heatmap(x,x,g)
readline()
