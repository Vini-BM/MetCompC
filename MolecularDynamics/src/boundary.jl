module Boundary

using Plots
using Printf
using LaTeXStrings
using Distributions
using Random
using LinearAlgebra

function reflect(size,positions,velocities)
    """
    Implements reflection of velocities on box walls for fixed boundary conditions.
    -------------------------------------------------------------------------------
    Arguments:
        size: array with dimensions X and Y
        positions: array with x and y coordinates
        velocities: array with x and y velocities
    """

    # Read parameters:
    ## Dimensions
    X = size[1]
    Y = size[2]
    ## Positions
    x = positions[1]
    y = positions[2]
    ## Velocities
    vx = velocities[1]
    vy = velocities[2]

    # Reflection on x
    zx = findall(t->(t<0)||(t>X),x) # find indices of particles outside of box on x axis
    vx[zx] *= -1 # reverse velocity

    # Reflection on y
    zy = findall(t->(t<0)||(t>Y),y) # same for y
    vy[zy] *= -1

    velocities = [vx,vy]
    return velocities
end