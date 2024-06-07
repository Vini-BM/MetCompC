module LJ

using Plots
using Printf
using LaTeXStrings
using Distributions
using Random
using LinearAlgebra

include("general.jl")
include("plotting.jl")

#############################################
########## LENNARD-JONES POTENTIAL ##########
#############################################

# Potential

function lennardJones(r,epsilon,sigma)
    """
    Calculates the Lennard-Jones potential for a pair of particles
    --------------------------------------------------------------
    Arguments:
        r: distance between particles
        epsilon: depth of the potential minimum
        sigma: distance at which the potential is zero
    """
    Vlj = 4*epsilon .* ((sigma ./ r) .^ 12 - (sigma ./ r) .^ 6)
    return Vlj
end

# Force

function forceModLJ(r,epsilon,sigma)
    """
    Returns the force modulus for a pair of particles subject to the Lennard-Jones potential
    --------------------------------------------------------------------------------------------
    Arguments:
        r: distance between particles
        epsilon: depth of the potential minimum
        sigma: distance at which the potential is zero
    """
    force = (24*epsilon ./ r) .* (2*(sigma ./ r) .^ 12 - (sigma ./ r) .^ 6)
    return force
end

function forceLJ(positions,epsilon,sigma)
    """
    Calculates the force between a pair of particles subject to the LJ potential
    ----------------------------------------------------------------------------
    Arguments:
        positions: matrix with x and y coordinates
        epsilon: depth of the potential minimum
        sigma: distance at which the potential is zero
    """

    # Position difference
    diff_x, diff_y, diff_r = General.positionDiff(positions)

    # Force modulus
    fmod = forceModLJ.(diff_r,epsilon,sigma) # calculates force modulus for each pair of particles
    fmod[diagind(fmod)] .= 0.0 # sets 'autoforce' to zero

    # Individual forces on x and y
    fx_ind = fmod .* diff_x # force between every pair of particles on x axis
    fy_ind = fmod .* diff_y # same for y

    # Resultant force
    fx = sum(fx_ind,dims=2) # resultant force for each particle on x axis
    fy = sum(fy_ind,dims=2) # same for y
    forces = hcat(fx,fy)

    return forces
end

# Potential energy
function potentialEnergyLJ(positions,epsilon,sigma)
    """
    Returns the potential energy for the system subject to the LJ potential.
    ---------------------------------------------------------------------------------------
    Arguments:
        positions: matrix with x and y coordinates
        epsilon: depth of the potential minimum
        sigma: distance at which the potential is zero
    """

    # Position difference
    diff_x, diff_y, diff_r = General.positionDiff(positions)

    # Calculate potential
    Vmatrix = lennardJones(diff_r,epsilon,sigma) # potential energy calculated from LJ on distances
    Vmatrix[diagind(Vmatrix)] .= 0.0 # sets 'autoenergy' to zero
    V = 0.25*sum(Vmatrix) # symmetry: divide by 2; two sums: divide by 4
    return V
end

####################
##### MOVEMENT #####
####################

# LJ
function moveLJ(positions,velocities,forces,epsilon,sigma,dt)
    """
    Integrates the equations of motion for the THP using Velocity-Verlet.
    Makes one timestep.
    --------------------------------------------------------------------
    Arguments:
        positions: matrix with x and y coordinates
        velocities: matrix with x and y velocities
        forces: matrix with initial x and y forces
        epsilon: depth of the potential minimum
        sigma: distance at which the potential is zero
        dt: timestep
    """
    # Update positions
    term = dt^2/2 # quadratic term on accelerated movement
    positions += velocities*dt + forces*term # uniformly accelerated movement
    # Update velocities
    newforces = forceLJ(positions,epsilon,sigma) # calculates new forces
    velocities += dt*(forces+newforces)/2 # mean over old and new forces
    return positions, velocities, newforces
end

# LJ time evolution
function evolveLJ!(size,positions,velocities,forces,epsilon,sigma,colors,tmax,dt,printstep)
    """
    Performs the time evolution of the system subject to the truncated harmonic potential.
    Does not return values; instead, plots the system configuration at each step.
    --------------------------------------------------------------------------------------
    Arguments:
        size: array with X and Y dimensions
        positions: matrix with x and y coordinates
        velocities: matrix with x and y velocities
        forces: matrix with x and y forces
        epsilon: depth of the potential minimum
        sigma: distance at which the potential is zero
        colors: array with particle colors for plotting
        tmax: final integration time
        dt: integration timestep
        printstep: timestep for plotting and calculating thermodynamic quantities
    """
    # Reading arguments
    N = length(positions[:,1]) # number of particles
    # Initialization
    t = 0 # time
    counter = 0 # counter for printing
    # Loop
    while t < tmax
        t += dt # updates time
        counter += 1 # updates counter
        positions, velocities, forces = moveLJ(positions,velocities,forces,epsilon,sigma,dt) # makes step
        velocities = General.reflect(size,positions,velocities)
        # Plotting
        label = "LJ"
        if counter%printstep == 0
            #println(pos[1,1])
            V = potentialEnergyLJ(positions,epsilon,sigma) # potential energy
            K = General.kineticEnergy(velocities) # kinetic energy
            temp = K/N # temperature = mean kinetic energy
            p = Plotting.frame(size,positions,t,colors,temp,K,V,label) # plot
            display(p) # show plot
        end
    end
    #return positions, velocities, forces
end

end # module