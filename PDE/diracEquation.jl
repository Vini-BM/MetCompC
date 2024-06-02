# Vinícius Müller
# 17-04-2024

# Dirac Equation in (1+1) dimensions

# x direction

using Plots
Plots.default(show=true)

function func(f, g, v, vsc, tf, dt, dx)
    t = 0
    ymin = 0
    ymax = 1
    pf = abs2.(f)
    pfint = sum(pf)*dx
    pg = abs2.(g)
    plot([x,x,x],[pf,pg,real(v)],ylim=(ymin,ymax),title="t=$t | $pfint")
    while t<tf
        oldf = deepcopy(f)
        oldg = deepcopy(g)
        f[2:end-1] = .5*(oldf[3:end]+oldf[1:end-2]) - dt/(2*dx)*(oldg[3:end]-oldg[1:end-2]) + im*(v[2:end-1]-(1 .+vsc[2:end-1]))*dt.*oldf[2:end-1]
        g[2:end-1] = .5*(oldg[3:end]+oldg[1:end-2]) - dt/(2*dx)*(oldf[3:end]-oldf[1:end-2]) + im*(v[2:end-1]+(1 .+vsc[2:end-1]))*dt.*oldg[2:end-1]
        #f[1] = 0
        #f[end] = 0
        #g[1] = 0
        #g[end] = 0
        t += dt
        pf = abs2.(f)
        pg = abs2.(g)
        pfint = sum(pf)*dx
        plot([x,x,x],[pf,pg,real(v)],ylim=(ymin,ymax),title="t=$t | $pfint")
        #plot([x,x],[abs2.(f),abs.(g)],ylim=(0,ymax))
    end
    return f, g
end

# Parameters
xmax = 200
dx = 1
x = collect(1:dx:xmax)
L = length(x)
dt = .01
#dxlist = [0.5,1,1.5]
t = 0

# Gaussian source
sigma = 6 # std
center = xmax/3 # center of gaussian
f = zeros(L) # initialize array
@. f = exp(-(x-center)^2 / 2sigma^2) # broadcasting
f = complex(f)
#println(f)
# g
g = zeros(L)
g = complex(g)
# v
#vsc = x .^(-1)
vsc = zeros(L)
vsc = complex(vsc)
v = zeros(L)
v[55] = 50
v[65] = 50
v = complex(v)
# Plot
#plot(x,real(f))
#plot!(x,g)


tf = 50
f, g = func(f,g,v,vsc,tf,dt,dx)
#for dx in dxlist
#f0 = copy(f)
#g0 = copy(g)
#ff, gf = func(f0,g0,tf,dt,dx)
#end
#plot!(x,real(f))
readline()




