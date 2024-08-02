using Random
using Plots

"""
Solves the equations of motion for a particle in a double potential well with additive noise.
Uses Stokes' approximation, in which the particle's inertia is neglected (mass ~ 0).

The potential is given by
    V(x) = ax^4 - bx^2

The stochastic equation of motion is then given by
    dX = (\alpha X - X^3)dt + \beta dW(t)

"""

function integrate(x0,tmax,dt,alpha,beta)
    # Initialize values
    t = 0.0
    x = x0
    tlist = [t]
    xlist = [x]

    # Loop
    while t <= tmax
        t += dt
        r = randn()
        x += alpha*x*dt - x^3*dt + beta*r*sqrt(dt)
        push!(tlist,t)
        push!(xlist,x)
    end
    return tlist, xlist
end


# Parameters
tmax = 10 # integration time
dt = .01 # integration  timestep
x0 = 6. # initial position
paramslist = [[2,1],[2,3],[2,5]] # alpha (potential parameter) and beta (noise)

# Loop on parameters
for params in paramslist
    alpha, beta = params
    # Integration
    tlist, xlist = integrate(x0,tmax,dt,alpha,beta)
    global p = plot!(tlist,xlist,label="alpha=$alpha, beta=$beta")
end

# Plot
p = plot!(xlabel="t",ylabel="x",title="Particle in double potential well")
display(p)
readline()
