module General

using LinearAlgebra

################################
########### BOUNDARY ###########
################################

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
    x = positions[:,1]
    y = positions[:,2]
    ## Velocities
    vx = velocities[:,1]
    vy = velocities[:,2]

    # Reflection on x
    zx = findall(t->(t<0)||(t>X),x) # find indices of particles outside of box on x axis
    vx[zx] *= -1 # reverse velocity

    # Reflection on y
    zy = findall(t->(t<0)||(t>Y),y) # same for y
    vy[zy] *= -1

    velocities = hcat(vx,vy)
    return velocities
end

#################################
###### POSITION DIFFERENCE ######
#################################

function positionDiff(positions)
    """
    Calculates the difference in position for each pair of particles.
    Works by performing element-by-element difference in array with its transposed,
    thus creating a matrix containing the differences in position.
    The diagonal is irrelevant but set to zero by definition;
    further, this matrix will be used in divisions,
    so the diagonal is set to one and the 'autoforces' will be set to zero.
    -----------------------------------------------------------------------
    Arguments:
        positions: array with x and y coordinates
    """

    # Reading arguments
    x = positions[:,1] # x positions
    y = positions[:,2] # y positions
    N = length(x) # number of particles

    # Difference in position
    diff_x = x .- x' # matrix with delta x for each pair of particles
    diff_y = y .- y' # same for y
    diff_r = sqrt.(diff_x .^2 + diff_y .^2) # difference in total distance: r^2 = x^2 + y^2
    diff_r += Diagonal(ones(N)) # distances will be used on division --> cannot divide by zero

    differences = [diff_x,diff_y,diff_r] # return array
    return differences
end

####################################
########## KINETIC ENERGY ##########
####################################

# Kinetic energy
function kineticEnergy(velocities)
    """
    Returns the kinetic energy for the system.
    -----------------------------------------
    Arguments:
        velocities: array with x and y velocities
    """

    # Reading arguments
    vx = velocities[:,1]
    vy = velocities[:,2]

    # Kinetic energy
    K = sum(vx .^2 + vy .^2)/2
    return K
end

end # module