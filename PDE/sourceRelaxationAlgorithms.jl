# @ Vinícius Müller
# 20-4-2024

# Compares four different methods to solve the Poisson equation and plots the final result for each one

using Plots

# Method: Simple Over-Relaxation
function overRelax(f,source,tol,alpha)
    L = size(f)[1] # length of x array
    t = 0
    diff = 1
    while diff > tol # loop
        lastf = copy(f) # copy of f
        t += 1 # update time
        for i in 2:L-1, j in 2:L-1 # loop through space
            f[i,j] = -alpha*lastf[i,j] + (1+alpha) * (f[i-1,j] + lastf[i+1,j] + f[i,j-1] + f[i,j+1] + source[i,j]) / 4 # update f
            # Fixed boundary conditions
        end
        diff = maximum(abs.(f-lastf))
    end
    return f, t, diff
end

# Method: Gauss-Siedel
function gaussSiedel(f,source,tol)
    L = size(f)[1] # length of x array
    t = 0
    diff = 1
    while diff > tol # loop
        lastf = copy(f) # copy of f
        t += 1 # update time
        for i in 2:L-1, j in 2:L-1 # loop through space
            f[i,j] = (lastf[i-1,j] + f[i+1,j] + f[i,j-1] + lastf[i,j+1] + source[i,j]) / 4 # update f
            # Fixed boundary conditions
        end
        diff = maximum(abs.(f-lastf))
    end
    return f, t, diff
end

# Method: Jacobi
function jacobi(f,source,tol)
    L = size(f)[1] # length of x array
    t = 0
    diff = 1
    while diff > tol # loop
        lastf = copy(f) # copy of f
        t += 1 # update time
        for i in 2:L-1, j in 2:L-1 # loop through space
            f[i,j] = (f[i-1,j] + f[i+1,j] + f[i,j-1] + f[i,j+1] + source[i,j]) / 4 # update f
            # Fixed boundary conditions
        end
        diff = maximum(abs.(f-lastf))
    end
    return f, t, diff
end

# Method: FTCS
function ftcs(f,source,tol,k)
    L = size(f)[1] # length of x array
    t = 0
    diff = 1
    while diff > tol # loop
        lastf = copy(f) # copy of f
        t += 1 # update time
        for i in 2:L-1, j in 2:L-1 # loop through space
            f[i,j] += k*(f[i-1,j] + f[i+1,j] + f[i,j-1] + f[i,j+1] - 4*f[i,j] + source[i,j]) / 4 # update f
            # Fixed boundary conditions
        end
        diff = maximum(abs.(f-lastf))
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
density = 1 # source value
source = zeros(L,L) # source term
source[Int(L/2),Int(L/2)] = density
tol = .00001 # tolerance
k = 1/2 # FTCS coefficient
alpha = .935 # Over-Relaxation coefficient
x = collect(1:1:L) # x array
plotlist = [] # list of plots
methods = ["FTCS", "Jacobi", "Gauss-Siedel", "Over-Relaxation"]


## Function
f = rand(Float64,(L,L)) # random initial condition
f[Int(L/2),Int(L/2)] = density # sets source value at center
@. f[end,1:end] = bcT # top
@. f[1,1:end] = bcB # bottom
@. f[1:end,1] = exp(-l1*x) # left
@. f[1:end,end] = exp(-l2*x) # right

# Loop

## FTCS
f1 = copy(f)
f1, t1, diff1 = ftcs(f1,source,tol,k)
println("$(methods[1]): iterations = $t1, difference = $diff1")
plot1 = heatmap(f1, title="$(methods[1]) ($t1 iterations)", xlabel="x", ylabel="y", colorbar_title="Temperature")
push!(plotlist,plot1)

## Jacobi
f2 = copy(f)
f2, t2, diff2 = jacobi(f2,source,tol)
println("$(methods[2]): iterations = $t2, difference = $diff2")
plot2 = heatmap(f2, title="$(methods[2]) ($t2 iterations)", xlabel="x", ylabel="y", colorbar_title="Temperature")
push!(plotlist,plot2)

## Gauss-Siedel
f3 = copy(f)
f3, t3, diff3 = gaussSiedel(f3,source,tol)
println("$(methods[3]): iterations = $t3, difference = $diff3")
plot3 = heatmap(f3, title="$(methods[3]) ($t3 iterations)", xlabel="x", ylabel="y", colorbar_title="Temperature")
push!(plotlist,plot3)

## Over-Relaxation
f4 = copy(f)
f4, t4, diff4 = overRelax(f4,source,tol,alpha)
println("$(methods[4]): iterations = $t4, difference = $diff4")
plot4 = heatmap(f4, title="$(methods[4]) ($t4 iterations)", xlabel="x", ylabel="y", colorbar_title="Temperature")
push!(plotlist,plot4)

# Plot
p = plot(plotlist...,plot_title="Relaxation Algorithms",layout=(2,2),size=(1000,800))
#savefig("Plots/relaxationAlgorithms.png")
display(p)
readline()
