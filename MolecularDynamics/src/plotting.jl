module Plotting

using Plots
using Printf
using LaTeXStrings
using Distributions
using Random
using LinearAlgebra

function frame(size,positions,time,colors,Km,K,V,potential_label)
    """
    Plots the particles' positions. For a 2D system only, but may be updated for a 3D one.
    -------------------------------------------------------------------------------------
    Arguments:
        size: array with dimensions
        positions: array with particles' positions
        time: current time counter
        colors: list of particles' colors
        Km: mean kinetic energy = temperature
        K: kinetic energy
        V: potential energy
        potential_label: string with potential name
    """

    # Reading arguments
    X = size[1]
    Y = size[2]
    x = positions[:,1]
    y = positions[:,2]
    N = length(x)

    # Energy
    E = K+V # total energy

    # Strings
    varstring = @sprintf("t=%.3f | E=%.2f | V=%.2f | K=%.2f | T=%.2f",time,E,V,K,Km)
    title = potential_label * " | " * varstring # plot title

    # Scatter
    global p=scatter(x,y,title=title,c=colors,ms=5,legend=false) # scatter particles

    # Annotations
    [annotate!(x, y+.1, Plots.text(string(i),12)) for (i,x,y) in zip(1:N,x,y)] # annotate particle labels

    # Limits
    plot!(p,xlimits=(0,X),ylimits=(0,Y))
    return p
end

end # module