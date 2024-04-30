using Distributions
using Plots

function randomSetup(X,Y,N,temp)
    N = Int(N)
    dt = 0.1
    x = X*rand(N)
    y = Y*rand(N)
    d = Normal(0,sqrt(temp)) # mean=0, var=temp -> d is a distribution 'object', not a random variable
    vx = rand(d,N) # generates a N-vector with random values distributed according to d
    vy = rand(d,N)
    Km = sum(vx.^2 + vy.^2)/(2*N) # mean kinetic energy -> using Kb = m = 1, Km = temp
    println("Mean Kinetic Energy: $Km")
    return x, y, vx, vy, temp, dt
end

