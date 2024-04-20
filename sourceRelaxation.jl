# @ Vinícius Müller
# 20-4-2024

using Plots
Plots.default(show=true)

# Method: Simple Over-Relaxation
function overRelax(f,source,alpha)
    L = size(f)[1] # length of x array
    lastf = copy(f) # copy of f
    for i in 2:L-1, j in 2:L-1 # loop through space
        f[i,j] = -alpha*lastf[i,j] + (1+alpha) * (f[i-1,j] + lastf[i+1,j] + f[i,j-1] + f[i,j+1] + source[i,j]) / 4 # update f
        # Fixed boundary conditions
    end
    diff = maximum(abs.(f-lastf))
    return f, diff
end

# Method: Gauss-Siedel
function gaussSiedel(f,source)
    L = size(f)[1] # length of x array
    lastf = copy(f) # copy of f
    for i in 2:L-1, j in 2:L-1 # loop through space
        f[i,j] = (lastf[i-1,j] + f[i+1,j] + f[i,j-1] + lastf[i,j+1] - source[i,j]) / 4 # update f
        # Fixed boundary conditions
    end
    diff = maximum(abs.(f-lastf))
    return f, diff
end

# Method: Jacobi
function jacobi(f,source)
    L = size(f)[1] # length of x array
    lastf = copy(f)
    for i in 2:L-1, j in 2:L-1 # loop through space
        f[i,j] = (f[i-1,j] + f[i+1,j] + f[i,j-1] + f[i,j+1] - source[i,j]) / 4 # update f
        # Fixed boundary conditions
    end
    diff = maximum(abs.(f-lastf))
    return f, diff
end

# Method: FTCS
function ftcs(f,source,k)
    L = size(f)[1] # length of x array
    lastf = copy(f)
    for i in 2:L-1, j in 2:L-1 # loop through space
        f[i,j] += k*(f[i-1,j] + f[i+1,j] + f[i,j-1] + f[i,j+1] - 4*f[i,j] - source[i,j]) / 4 # update f
        # Fixed boundary conditions
    end
    diff = maximum(abs.(f-lastf))
    return f, diff
end

# Loop

function loop(f,source,tol)
    plotlist = Vector{Any}(undef,4)
    difflist = ones(4)
    f_list = [copy(f) for i in 1:4]
    methods = ["FTCS", "Jacobi", "Gauss-Siedel", "Over-Relaxation"]
    t = 0 # timestep
    iterate = true
    anim = @animate while iterate # begins animation
        t += 1 # update time
        if difflist[1] > tol
            f1, diff1 = ftcs(f_list[1],source,1/2) # FTCS
            f_list[1] = f1
            difflist[1] = diff1
            plotlist[1] = heatmap(f1, title="$(methods[1]) ($t iterations)", xlabel="x", ylabel="y", colorbar_title="Temperature")
        end
        if difflist[2] > tol
            f2, diff2 = jacobi(f_list[2],source) # Jacobi
            f_list[2] = f2
            difflist[2] = diff2
            plotlist[2] = heatmap(f1, title="$(methods[2]) ($t iterations)", xlabel="x", ylabel="y", colorbar_title="Temperature")
        end
        if difflist[3] > tol
            f3, diff3 = gaussSiedel(f_list[3],source) # Gauss-Siedel
            f_list[3] = f3
            difflist[3] = diff3
            plotlist[3] = heatmap(f1, title="$(methods[3]) ($t iterations)", xlabel="x", ylabel="y", colorbar_title="Temperature")
        end
        if difflist[4] > tol
            f4, diff4 = overRelax(f_list[4],source,.935) # Over-Relaxation
            f_list[4] = f4
            difflist[4] = diff4
            plotlist[4] = heatmap(f1, title="$(methods[4]) ($t iterations)", xlabel="x", ylabel="y", colorbar_title="Temperature")
        end
        if count(x->x>tol,difflist) == 0
            iterate = false
            break
        end
    end
    return f_list, plotlist, anim
end


# Initialization

## Parameters
L = 100 # size
l1 = .1 # left exponent
l2 = .5 # right exponent
bcT = 0 # top BC
bcB = 1 # bottom BC
source = zeros(L,L) # source term
source[Int(L/2),Int(L/2)] = 1
tol = .00001 # tolerance
#alpha = .935 # method coefficient
x = collect(1:1:L) # x array

## Function
#f = 10*rand(Float64,(L,L)) # random initial condition
f = copy(source) # initial condition
@. f[end,1:end] = bcT # top
@. f[1,1:end] = bcB # bottom
@. f[1:end,1] = exp(-l1*x) # left
@. f[1:end,end] = exp(-l2*x) # right
#f[Int(L/2),Int(L/2)] = source

# Loop

g, plots, anim = loop(f,source,tol)
gif(anim)
readline()
