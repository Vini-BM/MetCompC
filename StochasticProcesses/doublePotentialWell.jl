using Random
using Plots

# Parameters
tmax = 10
dt = .1
alpha = 2
beta = .1

# Initialize values
t = 0.0
x = 6.0
tlist = [t]
xlist = [x]

# Loop
while t <= tmax
    global t += dt
    r = randn()
    global x += alpha*x*dt - x^3*dt + beta*r*sqrt(dt)
    push!(tlist,t)
    push!(xlist,x)
end
println(tlist)
plot!(tlist,xlist,display=true)
readline()
