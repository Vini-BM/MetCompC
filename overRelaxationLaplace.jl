# @ Vinícius Müller
# 10-04-2024

# Implements the "Over-Relaxation Method" for the Laplace equation
# Last updated: 18-04-2024
## Shows gif

using Plots
Plots.default(show=true)

# Method: Simple Over-Relaxation
function overRelax(f,alpha,tol)
    L = size(f)[1] # length of x array
    t = 0 # timestep
    iterate = true # loop control
    while iterate # condition based on tolerance
        lastf = deepcopy(f) # copy of f
        t += 1 # update time
        for i in 2:L-1, j in 2:L-1 # loop through space
            f[i,j] = -alpha*lastf[i,j] + (1+alpha) * (f[i-1,j] + lastf[i+1,j] + f[i,j-1] + f[i,j+1]) / 4 # update f
        end
        # Fixed boundary conditions
        global diff = maximum(abs.(f-lastf))
        if abs(diff) < tol
            iterate = false
            println("Difference = $diff")
        end
        heatmap(f, title="Thermal equilibrium ($t iterations)", xlabel="x", ylabel="y", colorbar_title="Temperature")
    end
    return f, t
end

# Initialization
## Parameters
L = 100 # size
l1 = .1 # left exponent
l2 = .5 # right exponent
bcT = 0 # top BC
bcB = 1 # bottom BC
tol = .00001 # tolerance
alpha = .7 # method coefficient
x = collect(1:1:L) # x array
## Function
f = 10*rand(Float64,(L,L)) # random initial condition
@. f[end,1:end] = bcT # top
@. f[1,1:end] = bcB # bottom
@. f[1:end,1] = exp(-l1*x) # left
@. f[1:end,end] = exp(-l2*x) # right

# Loop
g, itertime = overRelax(f,alpha,tol)
println("Iteration time = $itertime")
readline()
