# @ Vinícius Müller -- adpated from class
# 02-05-2024
# Calculates the execution time and memory allocation for a Molecular Dynamics program

using Plots
using Distributions
using Random
using LinearAlgebra
using Printf

function setup(sqN,distx,disty,T)

    # Define parameters
    N = sqN^2 # number of particles
    X = distx*sqN # number of sites on x axis
    Y = disty*sqN # number of sites on y axis
    ## if distx=disty, the system is a square and the lattice parameter is a^2 = distx^2 = X^2/N
    ## X*Y = cell_area*N

    # Initialize x and y
    x = Array{Float64,1}(undef,N)
    y = Array{Float64,1}(undef,N)
    for i = 1:N
        j=(i-1)%sqN # site indices on x axis --> division remainder
        l=div((i-1),sqN) # site indices on y axis --> integer division
        x[i]=distx*(j+1/2) # x coordinate = x index + center of square
        y[i]=disty*(l+1/2) # same for y
    end

    # Initialize velocities
    d = Normal(0,sqrt(T)) # normal distribution for velocities with mean=0 and variance=T
    Random.seed!(123456) # sets seed
    vx = rand(d,N) # x velocities
    vy = rand(d,N) # y velocities
    Tm = sum(vx .^2 + vy .^2)/(2N) # temperature = mean kinetic energy
    #println("Temperature = $Tm")

    # Initialize colors for plot
    colors = distinguishable_colors(N)
    return N,X,Y,x,y,vx,vy,Tm,colors
end

"""
function plotFrame(N,X,Y,x,y,Tm,t,colors,V,K,E)
     title=@sprintf("T=%.2f t=%.2f E=%.3f V=%.2f K=%.2f",Tm,t,E,V,K)
     global p=scatter(x,y,title=title,c=colors,ms=7,legend=false)
    [annotate!(x, y+.1, Plots.text(string(i),12)) for (i,x,y) in zip(1:N,x,y)]
     plot!(p,xlimits=(0,X),ylimits=(0,Y))
     return p
end
"""

function positionDiff(N,x,y)
     diff_x = x .- x' # matrix with delta x for each pair of particles
     diff_y = y .- y' # same for y
     diff_r = sqrt.(diff_x .^2 + diff_y .^2) # r^2 = x^2 + y^2
     diff_r += Diagonal(ones(N)) # distances will be used on division --> cannot divide by zero
     return diff_x,diff_y,diff_r
end

function boundary(X,Y,x,y,vx,vy)
     zx = findall(t->(t<0)||(t>X),x) # indices of particles outside of box on x axis
     vx[zx] *= -1 # reverse velocity
     zy = findall(t->(t<0)||(t>Y),y) # same for y
     vy[zy] *= -1
    return vx,vy
end

function forceMod(r,k0,r0)
    f = -k0 .*(r .- r0) # restoring force --> harmonic potential
return f
end

function force(N,x,y,k0,r0,r_lim)
    diff_x, diff_y, diff_r = positionDiff(N,x,y)
    f_mod = forceMod.(diff_r,k0,r0) # calculates force modulus for each pair of particles
    f_mod[diagind(f_mod)] .= 0.0 # sets 'autoforce' to zero
    z = findall(t->(t > r_lim),diff_r) # indices of pairs of particles farthest than r_lim
    f_mod[z] .= 0.0 # sets force between distant particles to zero
    fx_ind = f_mod .* diff_x # force between every pair of particles on x axis
    fy_ind = f_mod .* diff_y # same for y
    fx = sum(fx_ind,dims=2) # resultant force for each particle on x axis
    fy = sum(fy_ind,dims=2) # same for y
    return fx,fy
end

function potential(r,k0,r0,r_lim)
    Vr =  0.5 .*k0 .*(r .^2 .- r_lim .^2) .- k0 .*r0 .*(r .- r_lim) # truncated harmonic potential
    return Vr
end

function kineticEnergy(vx,vy)
        K = sum(vx .^2+vy .^2)/2
    return K
end

function potentialEnergy(N,x,y,k0,r0,r_lim)
    diff_x, diff_y, diff_r = positionDiff(N,x,y)
    Vmatrix = potential(diff_r,k0,r0,r_lim) # potential energy calculated from potential on distances
    Vmatrix[diagind(Vmatrix)] .= 0.0 # sets 'autoenergy' to zero
    z = findall(t->(t > r_lim),diff_r) # indices of pairs of particles farthest from r_lim
    Vmatrix[z] .= 0.0 # sets energy between distance particles to zero
    V = 0.25*sum(Vmatrix) # symmetry: divide by 2; two sums: divide by 4
    return V
end

function move(N,x,y,vx,vy,fxi,fyi,k0,r0,r_lim,dt)
    ## Makes one timestep
    ## Velocity-Verlet
    mdt2 = dt^2/2 # quadratic term on accelerated movement
    h = dt/2 # half step for velocities
    x += vx*dt + mdt2*fxi # update x position
    y += vy*dt + mdt2*fyi # update y position
    fxf, fyf = force(N,x,y,k0,r0,r_lim) # calculates new forces
    vx += h*(fxi+fxf) # updates x velocities --> mean over old and new forces
    vy += h*(fyi+fyf) # updates y velocities
    return x,y,vx,vy,fxf,fyf
end

function evol!(tmax,N,X,Y,x,y,vx,vy,fx,fy,k0,r0,r_lim,dt)
    ## Evolves the system until tmax time
    t = 0 # initialize time
    icount = 0 # initialize counter
    while t < tmax
        t += dt # updates time
        icount += 1 # updates counter
        x, y, vx, vy, fx, fy = move(N,x,y,vx,vy,fx,fy,k0,r0,r_lim,dt) # makes step
        vx, vy = boundary(X,Y,x,y,vx,vy)
        """
        # Print thermodynamic quantities every printstep steps
        printstep = 10
        if mod(icount,printstep) == 0
            V,K,TE = Total_Energy(x,y,vx,vy)
            p=plot_frame(x,y,Tm,t,N,X,Y,dy,c,V,K,TE)
            display(p)
            sleep(.1)
        end
        """
    end
end

# Program

## Parameters
k0 = 2000. # harmonic potential constant
r0 = 0.5 # harmonic potential center
r_lim = 1.2*r0 # separation distance for force
dx = r0 # lattice parameter for x
dy = r0 # lattice parameter for y
T = 1. # temperature
dt = 0.001 # timestep
tmax = 30 # integration time
sqNlist = collect(2:1:12) #list for √N

## Precompile
N,X,Y,x,y,vx,vy,Tm,colors = setup(2,dx,dy,T)
fx,fy = force(N,x,y,k0,r0,r_lim)
evol!(tmax,N,X,Y,x,y,vx,vy,fx,fy,k0,r0,r_lim,dt)

## Loop through N
timelist = []
memorylist = []
Nlist = sqNlist .^2
for sqN in sqNlist
    local N,X,Y,x,y,vx,vy,Tm,colors = setup(sqN,dx,dy,T)
    local fx,fy = force(N,x,y,k0,r0,r_lim)
    stats = @timed evol!(tmax,N,X,Y,x,y,vx,vy,fx,fy,k0,r0,r_lim,dt)
    time = stats.time
    push!(timelist,time)
    memory = stats.bytes
    push!(memorylist,memory)
    println("----- N: $N")
    println("Time: $time s")
    println("Allocated memory: $memory bytes")
end

## Linear fit
### Log
timelog = log10.(timelist)
memorylog = log10.(memorylist)
Nlog = log10.(Nlist)
### Angular coefficients
nlen = length(Nlist)
alpha_t = 0 # angular coefficient for time
alpha_m = 0 # angular coefficient for memory
ilist = []
for i in 2:nlen-1
    push!(ilist,i)
    global alpha_t += (timelog[i+1]-timelog[i])/(Nlog[i+1]-Nlog[i])
    global alpha_m += (memorylog[i+1]-memorylog[i])/(Nlog[i+1]-Nlog[i])
end
alpha_t /= length(ilist) # mean
alpha_m /= length(ilist) # mean
println("***** Fit *****")
println("Time: $alpha_t")
println("Allocated memory: $alpha_m")
### Linear coefficients --> just for plot
c_t = timelist[end-1]/Nlist[end-1]^alpha_t # linear coefficient for time
c_m = 1300000+(timelist[end]-timelist[1])/(Nlist[end]-Nlist[1])^alpha_m # linear coefficient for memory
### Fit
exptime = c_t * (Nlist .^ alpha_t) # power law function for time
expmem = c_m * (Nlist .^ alpha_m) # power law function for memory

## Plot
plotlist = []
### Titles and labels
timetitle = @sprintf("Execution time for tmax=%.1f | α=%.2f",tmax,alpha_t)
memorytitle = @sprintf("Allocated memory for tmax=%.1f | β=%.2f",tmax,alpha_m)
### Append subplots
#### Time
timeplot = plot(Nlist,timelist,xaxis=:log10,yaxis=:log10,xlabel="N",ylabel="Time (s)",title=timetitle,label="Data",margin=5Plots.mm)
plot!(Nlist,exptime,xaxis=:log10,yaxis=:log10,ls=:dash,label="Fit: t ∝ N^α")
push!(plotlist,timeplot)
#### Memory
memoryplot = plot(Nlist,memorylist,xaxis=:log10,yaxis=:log10,xlabel="N",ylabel="Allocated memory (bytes)",title=memorytitle,label="Data",margin=5Plots.mm)
plot!(Nlist,expmem,xaxis=:log10,yaxis=:log10,ls=:dash,label="Fit: M ∝ N^β")
push!(plotlist,memoryplot)
### Display
p = plot(plotlist...,size=(1200,600))
display(p)
readline()
savefig("dinamol_time_memory.png")