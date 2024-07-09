#This version uses structures for the system state (xv) and for the forces (f)
#The forces are calculated using Newton's third law, so just half are calculated
using Plots
using Random,Distributions
using LinearAlgebra
using Printf
using Setfield

struct Coord
    x::Vector{Float64}
    y::Vector{Float64}
    vx::Vector{Float64}
    vy::Vector{Float64}
end

struct forc
    fx::Vector{Float64}
    fy::Vector{Float64}
end

struct DynPar
    r0 :: Float64
    k0 :: Float64
    ra :: Float64
    dtt :: Float64
    tmax :: Float64
    X :: Float64
    Y :: Float64
end

struct boxes
    nc::Integer
    nf::Integer
    Lx::Integer
    Ly::Integer
    v1::Vector{Integer}
    v2::Vector{Integer}
    v3::Vector{Integer}
    v4::Vector{Integer}
    ncel::Matrix{Integer}
    lista::Matrix{Integer}
    num::Vector{Integer}
end
function inicia()
    sqN=10
    N=sqN^2
    Temp=0.001
    Random.seed!(1234)
    x,y,vx,vy,X,Y=init_pos_vel(sqN,Temp)
    xv = Coord(x,y,vx,vy)
    N=length(xv.x)
    X,Y,box=init_boxes(X,Y,N)
    dynpar=DynPar(0.5,50.,0.6,0.001,10.,X,Y) #r0,k0,ra,dt,tmax,X,Y
    Tm=sum(vx .^2+vy .^2)/(2N)
    println("Temperatura=$Tm")
    return xv,dynpar,box
end

function init_pos_vel(sqN,Temp)
    N = sqN^2
    x = Array{Float64,1}(undef,N)
    y = Array{Float64,1}(undef,N)
    d=Normal(0,sqrt(Temp))
    ddx = 0.5  #initial distance between particles along X 
    ddy = 0.5  #initial distance between particles along Y
    X = ddx*(sqN+1) #define X size
    Y = ddy*(sqN+1) #define Y size
    for i = 1:N
        j=(i-1)%sqN;l=div((i-1),sqN)
        if (l%2==0)
            x[i]=j*ddx+ddx/2+.5
        else
            x[i]=j*ddx+.5
        end

        if (j%2==0)
            y[i]=l*ddy+ddy/2+.5
        else
            y[i]=l*ddy+.5
        end
    end
    vx=rand(d,N)
    vy=rand(d,N)
    #vx .= 0.0
    #vy .= 0.0
    return x,y,vx,vy,X,Y
end

function init_boxes(X,Y,N)
    l = 2 # define box side
    X = ceil(Int,X/l)*l #redefine X to integer multiple of l 
    Y = ceil(Int,Y/l)*l #redefine Y to integer multiple of l
    Lx = Int(X/l) #number of boxes in X direction
    Ly = Int(Y/l)  #number of boxes in Y direction
    nc = Int(Lx*Ly);nf=nc+1
    #matriz identifying the number of the cells
    ncel=[(i-1)*Lx+j for i=1:Lx, j=1:Ly]
    #defining the neighbor lists
    v1 = Array{Integer,1}(undef,nc-1)
    v2 = Array{Integer,1}(undef,nc-1)
    v3 = Array{Integer,1}(undef,nc-1)
    v4 = Array{Integer,1}(undef,nc-1)
    #neighbors excluding borders
    for j=2:Ly-1
        for i=2:Lx-1
            c=ncel[i,j]
            v1[c]=ncel[i,j+1]   #right neighbor
            v2[c]=ncel[i+1,j-1]#left below neighbor
            v3[c]=ncel[i+1,j]  #below neighbor
            v4[c]=ncel[i+1,j+1]#right below neighbor
        end
    end
    #neighbors of the first column
    [(c=ncel[i,1];v1[c]=ncel[i,2];v2[c]=nf;v3[c]=ncel[i+1,1];v4[c]=ncel[i+1,2]) for i=1:Ly-1] 
    #neighbors of the last column
    [(c=ncel[i,Lx];v1[c]=nf;v2[c]=ncel[i+1,Lx-1];v3[c]=ncel[i+1,1];v4[c]=nf) for i=1:Ly-1] 
    #neighbors of the last line
    [(c=ncel[Ly,j];v1[c]=ncel[Ly,j+1];v2[c]=nf;v3[c]=nf;v4[c]=nf) for j=1:Ly-1]
    lista=zeros(nf,ceil(Integer,2*N/nc))
    num=zeros(nc)
    box = boxes(nc,nf,Lx,Ly,v1,v2,v3,v4,ncel,lista,num)
    return X,Y,box
end

function build_lists(xv,box)
    N=length(xv.x)
    for n=1:N
        i=floor(Integer,xv.y[n])+1  #particle n cell row
        j=floor(Integer,xv.x[n])+1  #particle n cell column
        c=box.ncel[i,j]             #particle box
        box.num[c]+=1               #increases counting in box c
        lista[c,box.num[c]]=n       #particle n is in box c
    end
    return xv,box
end

function plot_frame(xv,dynpar,t,V,K,TE)
    x=xv.x;y=xv.y
    N=length(x)
    dy = 0.1
    c=distinguishable_colors(N)
    Tm=sum(xv.vx .^2+xv.vy .^2)/(2N)
    XX = dynpar.X;    YY = dynpar.Y
    title=@sprintf("T=%.2f t=%.2f TE=%.3f V=%.2f K=%.2f",Tm,t,TE,V,K)
    global p=scatter(xv.x,xv.y,title=title,c=c,ms=7,legend=false)
    [annotate!(x,y+dy, Plots.text(string(i), 12)) for (i,x,y) in zip(1:N,x,y)]
    plot!(p,xlimits=(0,XX),ylimits=(0,YY))
    return p
end

function contorno(xv,dynpar)
    zx=findall(t->(t<0)||(t>dynpar.X),xv.x)
    xv.vx[zx] .= -xv.vx[zx]
    zy=findall(t->(t<0)||(t>dynpar.Y),xv.y)
    xv.vy[zy] .= -xv.vy[zy]
    return xv
end

function force_cells_intra(xv,dynpar,f,box)
    f.fx .= 0.0 ; f.fy = 0.0
    #in box forces
    for c=1:box.nc
        n=num[c]  #particles in box c
        if nc > 1
            fpx=zeros(n,n)
            fpy=zeros(n,n)
            listc=box.lista[c,1:n]
            for k=1:n-1
                #p=[k+1:n;]
                for p in k+1:n
                    nk=listc[k]    #identity of particle k
                    np=listc[p]    #identity of particle p
                    xij=xv.x[nk]-xv.x[np]
                    yij=xv.y[nk]-xv.y[np]
                    rij=sqrt(xij^2+yij^2)
                    if rij < dynpar.ra
                        fmod = -dynpar.k0*(rij-dynpar.r0)
                        f.fx[k]+=fmod*xij/rij
                        f.fy[k]+=fmod*yij/rij
                        f.fx[p]-=fmod*xij/rij
                        f.fy[p]-=fmod*yij/rij
                    end
		        end
	        end
        end
	end
    return f
end

function force_cells_inter(xv,dynpar,f,box)
# Now we calculate the forces among particles of c and the neighboring boxes
    for c=1:nc-1
        n=num[c]  #particles in cell c
        listv=[]
        if nc > 0
            listc=box.lista[c,1:n]
            n1=num[v1[c]]
            n2=num[v2[c]]
            n3=num[v3[c]]
            n4=num[v4[c]]
            nv=n1+n2+n3+n4
            if nv > 0
                if n1>0 push!(listv,box.lista[v1[c],1:n1])
                end
                if n2>0 push!(listv,box.lista[v1[c],1:n2])
                end
                if n3>0 push!(listv,box.lista[v1[c],1:n3])
                end
                if n4>0 push!(listv,box.lista[v1[c],1:n4])	    
                end   
            end
	    end

        for i in listc
	        for j in listv
	        	xij=xv.x[i]-xv.x[j]
                	yij=xv.y[i]-xv.y[j]
                	rij=sqrt(xij^2+yij^2)
                	if rij<dynpar.ra
                       fmod = -dynpar.k0*(rij-dynpar.r0)
                       fx = fmod*xij/rij
                       fy = fmod*yij/rij
                       f.fx[i] += fx
                       f.fy[i] += fy
                       f.fx[j] -= fx #usando a terceira lei de Newton permite 
                       f.fy[j] -= fy #que iniciemos o laço interno em i+1
	                end
            end
        end
     end
     return f
end

function force(xv,dynpar,f)
    f.fx .= 0.0 ; f.fy .= 0.0
    N = length(f.fx)
    for i=1:N
        for j=i+1:N
            xij=xv.x[i]-xv.x[j]
            yij=xv.y[i]-xv.y[j]
            rij=sqrt(xij^2+yij^2)
            if rij<dynpar.ra
                fmod = -dynpar.k0*(rij-dynpar.r0)
                fx = fmod*xij/rij
                fy = fmod*yij/rij
                f.fx[i] += fx
                f.fy[i] += fy
                f.fx[j] -= fx #usando a terceira lei de Newton permite 
                f.fy[j] -= fy #que iniciemos o laço interno em i+1       
            end
        end
    end
    return f
end

function PotentialEnergy(xv,dynpar)
    N=length(xv.x)
    V = 0.0
    for i=1:N
        for j=i+1:N
            xij=xv.x[i]-xv.x[j]
            yij=xv.y[i]-xv.y[j]
            rij=sqrt(xij^2+yij^2)
            if rij<dynpar.ra
                V+=0.5*dynpar.k0*(rij^2-dynpar.ra^2)-dynpar.k0*dynpar.r0*(rij-dynpar.ra)
            end
        end
    end
    return V
end

function Kinetic(xv)
    K = sum(xv.vx .^2+xv.vy .^2)/2
    return K
end

function Total_Energy(xv,dynpar)
    V = PotentialEnergy(xv,dynpar)
    K = Kinetic(xv)
    TE = V + K
    return V,K,TE
end

function move(xv,dynpar,f,box)
    fxi=f.fx;fyi=f.fy
    dt=dynpar.dtt
    mdt2=dt^2/2; h=dt/2
    xv.x .+= xv.vx*dt + mdt2*fxi 
    xv.y .+= xv.vy*dt + mdt2*fyi 
    #f=force(xv,dynpar,f)
    f = forces_cell_intra(xv,dynpar,f,box)
    f = forces_cell_inter(xv,dynpar,f,box)
    xv.vx .+= h*(fxi+f.fx) -0.95*xv.vx*dt
    xv.vy .+= h*(fyi+f.fy) -0.95*xv.vy*dt
    return xv,f
end


function evol(xv,f,dynpar,box)
#    f=force(xv,dynpar,f)
    xv, box = build_lists(xv,box)
    f=forces_cell_intra(xv,dynpar,f,box)
    f=forces_cell_inter(xv,dynpar,f,box)
    t=0;icount=1
    dt = dynpar.dtt
    while t < dynpar.tmax
        t +=dt;icount += 1
        xv,f=move(xv,dynpar,f,box)
        xv=contorno(xv,dynpar)
        xv, box = build_lists(xv,box)
        # if mod(icount,1000) == 0
        #     #println(t)
        #     V,K,TE = Total_Energy(xv,dynpar)
        #     p=plot_frame(xv,dynpar,t,V,K,TE) 
        #     display(p)
        #     #sleep(.1)
        # end
    end
end
#**********main***********
xv,dynpar,box=inicia()
N=length(xv.x)
println("N=$N tmax=$(dynpar.tmax)")
f=forc(zeros(N),zeros(N))
@time evol(xv,f,dynpar,box)
@time evol(xv,f,dynpar,box)
readline()
