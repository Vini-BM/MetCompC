using Random
using Plots
using Statistics

function ornst_uhl(x0,tmax,dt,k,beta)
    sqrt_dt = sqrt(dt)
    num_steps = Int(ceil(tmax/dt))
    x = zeros(num_steps)
    x[1] = x0
    dW = sqrt_dt*randn(num_steps)
    t = dt.*collect(0:num_steps-1)
    #println(length(t))
    #println(length(x))
    for i in 1:num_steps-1
        x[i+1] = x[i] -k*x[i]*dt + beta*dW[i]
    end
    return t, x
end

function mean_ou(N::Int,x0,tmax,dt,k,beta)
    tlist = []
    xlist = []
    for i in N
        t, x = ornst_uhl(x0,tmax,dt,k,beta)
        push!(tlist,t)
        push!(xlist,x)
    end
    
end

x0 = 1
dt = 0.1
beta = 2.0
k = 1
tmax = 10

t,x = ornst_uhl(x0,tmax,dt,k,beta)
plot(t,x,display=true)
readline()
