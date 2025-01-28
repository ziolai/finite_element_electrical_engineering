using LinearAlgebra
using StaticArrays 
using SparseArrays

const Point2D = SVector{2,Float64};

@enum bcType begin
           dirichlet    = 1
           neumann      = 2     
           periodic     = 3
           antiperiodic = 4 
       end

struct bcPatch
    bctype::bcType     # type of boundary condition being used 
    bcvalue::Float64   # value of boundary condition being used
end

""" 
Generates a one-dimensional uniform mesh between point 0 and 1

Input N(1) is number of elements. Output is the one-dimensional mesh. 

[`mesh(::NTuple{1,Int})`](@ref)
"""
function fdmmesh(N::NTuple{1,Int})
    m = tuple(range(0,1,length=N[1]+1))  
end

"""
Generates a one-dimensional uniform mesh between point a and b

Input N(1) is number of elements. Output is the one-dimensional mesh. 

[`mesh(::NTuple{1,Int},::Real,::Real)`](@ref)
"""
function fdmmesh(N::NTuple{1,Int},a::Real,b::Real)
    if (a>=b)
        ErrorException("Wrong input")
    else 
        m = tuple(range(a,b,length=N+1)) 
    end 
end 

"""
Generates a two-dimensional mesh between the points (0,0) and (1,1) 

Input N(1) and N(2) are the number of elements in x and y direction. Output is the two-dimensional mesh. 

[`mesh(::NTuple{2,Int})`](@ref)
"""
function fdmmesh(N::NTuple{2,Int})
    m = (range(0,1,length=N[1]+1),range(0,1,length=N[2]+1))
end 

""" 
Generates a two-dimensional mesh between the points a and b 

Input N(1) and N(2) are the number of elements in x and y direction. Output is the two-dimensional mesh. 

[`mesh(::NTuple{2,Int},::Point2D,::Point2D)`](@ref)
"""
function fdmmesh(N::NTuple{2,Int}, a::Point2D,b::Point2D)
    if ((a.x>=b.x)||(a.y>=b.y))
        ErrorException("Wrong input")
    else         
        m = (range(a.x,b.x,length=N[1]+1),range(a.y,b.y,length=N[2]+1))
    end 
end

"""
Generates a vector on a one-dimensional uniform mesh by nodal evaluation a function
"""
function loadvec(mesh::NTuple{1,AbstractRange},sourceFct::Function)  
    f = sourceFct.(mesh[1])
    return f      
end

""" 
Generates a vector on a two-dimensional uniform mesh by nodal evaluation a function
Remark: we deliberately store the vector as a 2D array to facilitate the setting on boundary conditions 
"""
function loadvec(mesh::NTuple{2,AbstractRange},sourceFct::Function)  
    f = vec(sourceFct.(mesh[1], mesh[2]')) 
    return f      
end

"""
Generates the one-dimensional diffusion or graph Laplacian matrix 

Input r is the one-dimensional mesh. 
Output A is the one-dimensional diffusion matrix with stencil (1/h2)*[-1 2 -1]. 
Remark: we deliberately convert from Tridiagonal to sparse to allow for periodic boundary conditions at later stage.
""" 
function stiffmat(mesh::NTuple{1,AbstractRange})
    xm = mesh[1]
    Np1 = length(xm)
    N = Np1-1; h = abs(xm[end] - xm[1])/N; h2 = h*h
    e = ones(Np1) 
    A = spdiagm(-1 => -e[2:end], 0 => 2*e, 1 => -e[2:end])
    A = (1/h2)*A
    return A     
end

"""
Generates the one-dimensional convection-diffusion matrix 

Input mesh is the one-dimensional mesh. 
Output A is the one-dimensional diffusion matrix with stencil (1/h)*[0 -1 1]. 
""" 
function convectionmat(mesh::NTuple{1,AbstractRange})
    xm = mesh[1]
    Np1 = length(xm)
    N = Np1-1; h = abs(xm[end] - xm[1])/N
    e = ones(Np1) 
    A = spdiagm(0 => -e, 1 => e[2:end])
    A = (1/h)*A
    return A     
end

""" 
Generates the two-dimensional diffusion or graph Laplacian matrix 
""" 
function stiffmat(mesh::NTuple{2,AbstractRange})
    xm = mesh[1]; ym = mesh[2]  
    Nxp1 = length(xm); Nyp1 = length(ym);   
    #..contribution in x-direction 
    A1Dxx = stiffmat(tuple(xm))
    A2Dxx = kron(A1Dxx,I(Nyp1))
    #..contribution in y-direction
    A1Dyy = stiffmat(tuple(ym))
    A2Dyy = kron(I(Nxp1),A1Dyy)
    #..sum of x-direction and y-direction contribution 
    A = A2Dxx + A2Dyy
    return A 
end

""" 
Generates boundary condition entries for matrix A and right-hand side vector f 
Remark: after treatment of boundary conditions, the matrix is no longer symmetric. 
"""
function applybc!(mesh::NTuple{1,AbstractRange},bc::NTuple{2,bcPatch},A,f)
 
    #..retrieve mesh information
    xm = mesh[1];
    Np1 = length(xm)
    N = Np1-1; h = abs(xm[end] - xm[1])/N; h2 = h*h
    
    #..treat coupled boundary conditions 
    patch = bc[1]
    if (patch.bctype==bcType(3))
        A = sparse(A)
        A[1,end-1] = -1/h2
        A[end,2]   = -1/h2 
        return A,f 
    end 
            
    #..treatment of left bc 
    patch = bc[1]           
    if (patch.bctype==bcType(1)) 
        A[1,1] = 1.; A[1,2] = 0. 
        f[1] = patch.bcvalue
    elseif (patch.bctype==bcType(2))
        A[1,1] = 1/h ;  A[1,2] = -1/h
        f[1] = patch.bcvalue 
    else
         error(" handle_bc!:: Error:: Left boundary conditions not set")           
    end 

    #..treatment of right bc 
    patch = bc[2]           
    if (patch.bctype==bcType(1)) 
        A[end,end] = 1.; A[end,end-1] = 0. 
        f[end] = patch.bcvalue
    elseif (patch.bctype==bcType(2))
        A[end,end] = 1/h ;  A[end,end-1] = -1/h 
        f[end] = patch.bcvalue 
    else
         error(" handle_bc!:: Error:: Right boundary conditions not set")           
    end 
        
    return A,f  

end

function applybc!(mesh::NTuple{2,AbstractRange},bc::NTuple{4,bcPatch},A,f)

    # requires chechking what is x and y
    Nxp1 = length(mesh[1])
    Nyp1 = length(mesh[2])
    linear = LinearIndices((1:Nxp1, 1:Nyp1)) 
    
    #..treatment of bottom bc from left to right 
    patch = bc[1]  
    for i in 1:Nxp1 
        I = linear[i,1]
        f[I] = patch.bcvalue 
        A[I,:] .= 0.; A[I,I] = 1.
    end

    #..treatment of right bc from bottom to top 
    patch = bc[2]           
    for j = 1:Nyp1  
        I = linear[Nxp1,j]
        f[I] = patch.bcvalue 
        A[I,:] .= 0.; A[I,I] = 1.       
    end
        
    #..treatment of top bc from left to right
    patch = bc[3]           
    for i in 1:Nxp1
        I = linear[i,Nyp1]
        f[I] = patch.bcvalue 
        A[I,:] .= 0.; A[I,I] = 1.       
    end
        
    #..treatment of left bc from bottom to top
    patch = bc[4]           
    for j = 1:Nyp1
        I = linear[1,j]
        f[I] = patch.bcvalue 
        A[I,:] .= 0.; A[I,I] = 1.        
    end
    
    return A,f 
        
end 

""" 
Generates discrete solution vector u - rename function later 
"""
function fdmsolve(A,f)
    u = A\f 
    return u 
end
;