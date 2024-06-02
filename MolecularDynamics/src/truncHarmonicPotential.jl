module THP

using Plots
using Printf
using LaTeXStrings
using Distributions
using Random
using LinearAlgebra

include("general.jl")
include("plotting.jl")

##################################################
########## TRUNCATED HARMONIC POTENTIAL ##########
##################################################

# Potential
function truncHarmPotential(r,k0,r0,rlim)
    """
    Calculates the truncated harmonic potencial for the system.
    -----------------------------------------------------------
    Arguments:
        r: matrix with distances between particles
        k0: coupling constant
        r0: oscillation center
        rlim: maximum distance for which the potential is nonzero
    """
    # Potential
    Vr =  0.5 .*k0 .*(r .^2 .- rlim .^2) .- k0 .*r0 .*(r .- rlim) # potential calculated with broadcasting
    return Vr
end

# Force modulus
function forceModTHP(r,k0,r0)
    """
    Returns the force modulus for a pair of particles at distance r from each other.
    -------------------------------------------------------------------------------
    Arguments:
        r: distance
        k0: coupling constant
        r0: oscillation center
    """
    # Force
    force = -k0 * (r .- r0) # restoring force
    return force
end

# Forces
function forceTHP(positions,k0,r0,rlim)
    """
    Calculates the restoring force from the truncated harmonic potential for each pair of particles.
    -----------------------------------------------------------------------------------------------
    Arguments:
        positions: array with x and y coordinates
        k0: coupling constant
        r0: oscillation center
        rlim: maximum distance for which the potential is nonzero
    """

    # Position difference
    diff_x, diff_y, diff_r = General.positionDiff(positions)

    # Force modulus
    fmod = forceModTHP.(diff_r,k0,r0) # calculates force modulus for each pair of particles
    fmod[diagind(fmod)] .= 0.0 # sets 'autoforce' to zero
    z = findall(t->(t > rlim),diff_r) # indices of pairs of particles farthest than rlim
    fmod[z] .= 0.0 # sets force between distant particles to zero
    fx_ind = fmod .* diff_x # force between every pair of particles on x axis
    fy_ind = fmod .* diff_y # same for y
    fx = sum(fx_ind,dims=2) # resultant force for each particle on x axis
    fy = sum(fy_ind,dims=2) # same for y
    forces = hcat(fx,fy)
    return forces
end

# Potential energy
function potentialEnergyTHP(positions,k0,r0,rlim)
    """
    Returns the potential energy for the system subject to the truncated harmonic potential.
    ---------------------------------------------------------------------------------------
    Arguments:
        positions: array with x and y coordinates
        k0: coupling constant
        r0: oscillation center
        rlim: maximum distance for which the potential is nonzero
    """

    # Position difference
    diff_x, diff_y, diff_r = General.positionDiff(positions)

    # Calculate potential
    Vmatrix = truncHarmPotential(diff_r,k0,r0,rlim) # potential energy calculated from THP on distances
    Vmatrix[diagind(Vmatrix)] .= 0.0 # sets 'autoenergy' to zero
    z = findall(t->(t>rlim),diff_r) # indices of pairs of particles farthest from r_lim
    Vmatrix[z] .= 0.0 # sets energy between distant particles to zero
    V = 0.25*sum(Vmatrix) # symmetry: divide by 2; two sums: divide by 4
    return V
end

####################
##### MOVEMENT #####
####################

# THP
function moveTHP(positions,velocities,forces,k0,r0,rlim,dt)
    """
    Integrates the equations of motion for the THP using Velocity-Verlet.
    Makes one timestep.
    --------------------------------------------------------------------
    Arguments:
        positions: array with x and y coordinates
        velocities: array with x and y velocities
        forces: array with initial x and y forces
        k0: coupling constant
        r0: oscillation center
        rlim: maximum distance for which the potential is nonzero
        dt: timestep
    """
    # Update positions
    term = dt^2/2 # quadratic term on accelerated movement
    positions += velocities*dt + forces*term # uniformly accelerated movement
    # Update velocities
    newforces = forceTHP(positions,k0,r0,rlim) # calculates new forces
    velocities += dt*(forces+newforces)/2 # mean over old and new forces
    return positions, velocities, newforces
end

# THP time evolution
function evolveTHP!(size,positions,velocities,forces,k0,r0,rlim,colors,tmax,dt,printstep)
    """
    Performs the time evolution of the system subject to the truncated harmonic potential.
    Does not return values; instead, plots the system configuration at each step.
    --------------------------------------------------------------------------------------
    Arguments:
        size: array with X and Y dimensions
        positions: array with x and y coordinates
        velocities: array with x and y velocities
        forces: array with x and y forces
        k0: coupling constant
        r0: oscillation center
        rlim: maximum distance for which the potential is nonzero
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
        positions, velocities, forces = moveTHP(positions,velocities,forces,k0,r0,rlim,dt) # makes step
        velocities = General.reflect(size,positions,velocities)
        # Plotting
        label = "THP"
        if counter%printstep == 0
            #println(pos[1,1])
            V = potentialEnergyTHP(positions,k0,r0,rlim) # potential energy
            K = General.kineticEnergy(velocities) # kinetic energy
            temp = K/N # temperature = mean kinetic energy
            p = Plotting.frame(size,positions,t,colors,temp,K,V,label) # plot
            display(p) # show plot
        end
    end
    #return positions, velocities, forces
end

end # module