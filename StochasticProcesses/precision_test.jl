using Random
using Plots



function integrate(tmax,dt_big,m=10)
    # Parameters
    dt_small = dt_big/m
    intervals = Int(round(tmax/dt_big)) # intervals for time integration
    # at each interval we use m steps of dt_small
    sqrt_dt_small = sqrt(dt_small)
    
    # Initialization
    W_small = 0 # counter for Wiener increment with small dt
    W_big = 0 # counter for Wiener increment with big dt
    W = zeros(intervals+1) # one for each interval + one for initial condition
    int_small = zeros(intervals+1)
    int_big = zeros(intervals+1)
    t = collect(0:dt_big:tmax) # array with time
    integral_counter = 0

    # Loop
    for j in 1:intervals
        dW_big = 0 # Wiener differential
        for k in 1:m # small interval
            dW_small = sqrt_dtp * (randn(0)) # gaussian noise
            dW_big += dW_small # update differential with big dt 
            # update integral inside interval
            integral_counter += W_small^2 * dW_small + W_small*dt_small 
            W_small += dW_small # update noise
        end
        int_small[j+1] = integral_counter # update integral with small dt
        int_big[j+1] = int_big[j] + W_big^2 * dW_big + W_big * dt_big # update integral with big dt
        W_big += dW_big # update noise
        W[j] = W_big
    end
    return t, W, int_small, int_big
end

# Numerical
m = 10
t, W, int_small, int_big = integrate(tmax,dt_big,m)

# Analytic
int_ito = W.^3 / 3

# Plot
plot!(t,int_ito,label="Analytic")
plot!(t,int_small,label="Small dt")
plot!(t,int_big,label="Big dt")
readline()
