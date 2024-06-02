using MolecularDynamics

# Parameters
k0 = 200. # harmonic potential constant
r0 = 0.5 # harmonic potential center
rlim = 1.2*r0 # separation distance for force
dx = r0 # lattice parameter for x
dy = r0 # lattice parameter for y
T = 1.2 # temperature
dt = 0.001 # timestep
tmax = 50 # integration time
printstep = 10
sqN = 10

seed = 2352528
N, size, pos, vel, Km, colors = Setup.normalSetup(sqN,dx,dy,T,seed)
#println(pos)
#println(length(pos))
#println(length(pos[1]))
#println(length(pos[2]))
#println(pos)
forces = THP.forceTHP(pos,k0,r0,rlim)
#println(forces[1])
THP.evolveTHP!(size,pos,vel,forces,k0,r0,rlim,colors,tmax,dt,printstep)
#pos, vel, forces = THP.moveTHP(pos,vel,forces,k0,r0,rlim,dt)
#println(forces[1])
#println(pos)