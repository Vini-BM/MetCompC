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
tmax = 5 # integration time
dt = .01 # integration  timestep
alpha = 2 # potential parameter
beta = 1 # noise
x0 = 6. # initial position

tlist, xlist = integrate(x0,tmax,dt,alpha,beta)
plot(tlist,xlist,xlabel="t",ylabel="x",title="Particle in double potential well",display=true)
readline()
