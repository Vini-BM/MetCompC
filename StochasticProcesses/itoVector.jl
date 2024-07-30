using Random
using Plots

function func()
tmax = 10
dt = 0.01
w0 = 1
beta = 0.05
b2 = beta^2/2
sqrtdt = sqrt(dt)

num_floor = floor(Int,tmax/dt)
num = num_floor + 1
xi = zeros(num)
yi = zeros(num)
xe = zeros(num)
ye = zeros(num)
xe[1] = 1.
xi[1] = 1.
dW = sqrtdt*randn(num_floor)
W = 0

for j in 1:num_floor
    W += dW[j]
    phase = w0*dt*j + beta*W

    # Ito
    xij = xi[j]
    yij = yi[j]
    xi[j+1] = xij - (w0*yij+b2*xij)*dt - yij*beta*dW[j]
    yi[j+1] = yij + (w0*xij-b2*yij)*dt + xij*beta*dW[j]

    # Exact
    xe[j+1] = cos(phase)
    ye[j+1] = sin(phase)
end
return xi, yi, xe, ye
end

#tlist = dt*collect(0:num_floor)

xi, yi, xe, ye = func()

plot(xe,ye,label="Analytic")
plot!(xi,yi,label="Ito")
#println(xi)