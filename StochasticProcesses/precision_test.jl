using Random
using Plots

# Parameters
tmax = 10
dt = 0.1

dtg = dt
dtp = dt/10
intervals_g = Int(round(tmax/dtg))
sqrt_dtp = sqrt(dtp)

# Initialization
Wp = 0
Wg = 0
W = zeros(intervals_g+1)
intp = zeros(intervals_g+1)
intg = zeros(intervals_g+1)
t = collect(0:dtg:tmax)
integral_counter = 0

# Loop
for j in 1:intervals_g
    dWg = 0
    for k in 1:10
        dWp = sqrt_dtp * (randn(0))
        dWg += dWp
        global integral_counter += Wp^2 * dWp + Wp*dtp
        global Wp += dWp
    end
    global intp[j+1] = integral_counter
    global intg[j+1] = intg[j] + Wg^2 * dWg + Wg * dtg
    global Wg += dWg
    global W[j] = Wg
end

# Analytic
int_ito = W.^3 / 3

println(length(t))
println(length(int_ito))
println(length(intp))
# Plot
plot!(t,int_ito,label="Analytic")
plot!(t,intp,label="Small dt")
plot!(t,intg,label="Big dt")