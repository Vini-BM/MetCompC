module Setup

using Plots
using Printf
using LaTeXStrings
using Distributions
using Random
using LinearAlgebra

function normalSetup(sqN,dx,dy,T,seed)
    """
    Initializes the system using a normal distribution with mean=0 and variance=temperature.
    The system is a box with dimensions X and Y.
    To each particle is assigned a unit cell.
    Each unit cell has dimensions dx and dy.
    The system is then composed of N unit cells: S = X*Y = N*dx*dy.
    ---------------------------------------------------------------
    Arguments:
        sqN: square root of N
        dx: x-dimension of unit cell
        dy: y-dimension of unit cell
        T: temperature of ensemble
        seed: normal distribution seed
    """

    # Set seed
    Random.seed!(seed)

    # Define parameters
    N = sqN^2 # number of particles
    X = dx*sqN # number of sites on x axis
    Y = dy*sqN # number of sites on y axis
    size = [X,Y] # return array with size

    # Initialize positions
    x = Array{Float64,1}(undef,N)
    y = Array{Float64,1}(undef,N)
    for i = 1:N
        j=(i-1)%sqN # site indices on x axis --> division remainder
        l=div((i-1),sqN) # site indices on y axis --> integer division
        x[i]=dx*(j+1/2) # x coordinate = x index + center of square
        y[i]=dy*(l+1/2) # same for y
    end
    positions = hcat(x,y) # return matrix with positions

    # Initialize velocities
    d = Normal(0,sqrt(T)) # normal distribution for velocities with mean=0 and variance=T
    vx = rand(d,N) # x velocities
    vy = rand(d,N) # y velocities
    velocities = hcat(vx,vy) # return matrix with velocities

    # Recalculate temperature
    Km = sum(vx .^2 + vy .^2)/(2*N) # temperature = mean kinetic energy

    # Initialize colors for plot
    colors = distinguishable_colors(N) # returns array with N random colors --> for plotting purposes
    
    return N,size,positions,velocities,Km,colors
end

end # module