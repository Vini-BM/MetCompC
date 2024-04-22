# @ Vinícius Müller
# 21-04-2024

# Implements the "Over-Relaxation Method" for the Poisson equation and creates an animation
# The animation is made simply by replotting the heatmap every iteration

using Plots
using LaTeXStrings # necessary to display LaTeX in plots
Plots.default(show=true)

# Method: Simple Over-Relaxation
function overRelax(f,source,alpha,tol)
    t = 0 # timestep
    diff = 1 # initialize loop control
    while diff > tol # condition based on tolerance
        lastf = deepcopy(f) # copy of f
        t += 1 # update time
        for i in 2:L-1, j in 2:L-1 # loop through space
            f[i,j] = -alpha*lastf[i,j] + (1+alpha) * (f[i-1,j] + lastf[i+1,j] + f[i,j-1] + f[i,j+1] + source[i,j]) / 4 # update f
        end
        # Fixed boundary conditions
        diff = maximum(abs.(f-lastf)) # difference between current and previous values
        heatmap(f, title=L"Over-Relaxation | $\alpha$ = %$alpha | %$t iterations", xlabel=L"$x$", ylabel=L"$y$", colorbar_title="Temperature")
    end
    return f, t, diff
end

# Initialization

## Parameters

L = 100 # size
l1 = .1 # left exponent
l2 = .5 # right exponent
bcT = 0 # top BC
bcB = 1 # bottom BC
tol = .000001 # tolerance
alpha = .935 # method coefficient
x = collect(1:1:L) # x array
density = 1 # source term
source = zeros(L,L) # initializes source array
source[Int(L/2),Int(L/2)] = density # sets source at center
# In order for the method to work, the source must be treated as an array; here it is simply [0,...0,1,0,...,0]

## Function

### Initialize vector
f = 20*rand(Float64,(L,L)) # random initial condition with values (0,20) for a more interesting animation
### Set source
f[Int(L/2),Int(L/2)] = density # sets source term at center
### Boundary conditions
@. f[end,1:end] = bcT # top
@. f[1,1:end] = bcB # bottom
@. f[1:end,1] = exp(-l1*x) # left
@. f[1:end,end] = exp(-l2*x) # right

# Loop

g, itertime, diff = overRelax(f,source,alpha,tol)
println("Over-Relaxation")
println("Iteration time = $itertime")
println("Difference = $diff")
readline() # keep plot open
