"""
    struct Element 

Holds information for a single linear tetrahedral element. 
"""
struct Elem3DLin
  p1::Point3D                      # coordinates first node 
  p2::Point3D                      # coordinates second node 
  p3::Point3D                      # coordinates third node  
  p4::Point3D                      # coordinates fourth node      
  node_elem_tag::SVector{12,Int64} # node tags for element 
  Emat::MMatrix{4,4,Float64, 16}   # matrix of basis function coefficients 
  vol::Float64                     # area of the element 
end; 

"""
    struct Mesh 

Holds information for the entire mesh as an array of linear triangular elements. 
"""
struct Mesh 
  nnodes::Int64                # number of nodes 
  nelems::Int64                # number of elements
  dofPerElem::Int64            # number of dofs per element   
  Elements::Array{Elem3DLin,1} # list of Elements 
end; 

"""
    eval_volume(p1,p2,p3,p4)

Evaluate volume of tetrahedron with vertices p1,p2,p3 and p4  
"""
function evalVol(p1,p2,p3,p4)
    v1 = p2 - p1
    v2 = p3 - p1
    v3 = p4 - p1
    return abs(dot(cross(v1, v2), v3)) / 6.0
end; 

"""
    gen_basis(p1,p2,p3,p4)

Generate basis functions of tetrahedron with vertices p1,p2,p3 and p4 
"""
function genBasis(p1,p2,p3,p4) 
    Xmat = SMatrix{4,4,Float64, 16}(
        p1[1], p2[1], p3[1], p4[1], 
        p1[2], p2[2], p3[2], p4[2],
        p1[3], p2[3], p3[3], p4[3],
        1., 1., 1., 1.) 
    rhs  = SMatrix{4,4,Float64, 16}(1I) 
    Emat = SMatrix{4,4,Float64, 16}(Xmat\rhs);
    return Emat
end; 

"""
    gen_basis(p1,p2,p3,p4)

Evaluate at x the basis functions of element   
"""
function evalBasis(x,element)
    return Transpose(element.Emat)*push(x,1.)
end;  

"""
Generate local 12x1 vector for a tetrahedron.
"""
function genLocVec(element)
        
    base = SVector{12,Float64}(
        1.0, 1.0, 1.0, 1.0,
        0., 0., 0., 0.,
        0., 0., 0., 0.,
    )
    
    return (element.vol / 4) * base
end;

"""
Generate local 12x12 consistent mass matrix for a tetrahedron.
"""
function genLocMassMat(element)
        
    # Generate the base template matrix for P1 tets
    # (Diagonal elements = 2, Off-diagonal = 1)
    base = SMatrix{4, 4, Float64, 16}(
        2.0, 1.0, 1.0, 1.0,
        1.0, 2.0, 1.0, 1.0,
        1.0, 1.0, 2.0, 1.0,
        1.0, 1.0, 1.0, 2.0
    )
    
    base = (element.vol / 20.0) * base
    
    Z = zeros(SMatrix{4,4,Float64})
    
    ext_base = hcat(
        vcat(base, Z, Z),
        vcat(Z, base, Z),
        vcat(Z, Z, base)
    )
    
    return ext_base 
end;

# define kernel as function of destination and source coordinates 
kernel(rd, rs) = 1/(norm(rd - rs)+1e-6); 

"""
    integrate_tet(f, element; rtol=1e-8, atol=1e-8)

Integrates a 3D function `f(p)` over a tetrahedron defined by vertices 
`p1`, `p2`, `p3`, and `p4`. `f` should accept a 3-element vector/tuple.
"""
function integrate_tet(f, element; rtol=1e-8, atol=1e-8)
    p1 = element.p1; p2 = element.p2; 
    p3 = element.p3; p4 = element.p4; 
    # Compute the geometric Jacobian matrix and its determinant
    # This maps the reference simplex to the physical tetrahedron
    J_geom_mat = hcat(p2 .- p1, p3 .- p1, p4 .- p1)
    jac_geom = abs(det(J_geom_mat))
    
    # Define the transformed integrand over the unit cube [0, 1]^3
    function unit_cube_integrand(t)
        x, y, z = t[1], t[2], t[3]
        
        # Duffy transformation: unit cube -> reference simplex coordinates
        u = x
        v = y * (1.0 - x)
        w = z * (1.0 - x) * (1.0 - y)
        
        # Map reference coordinates to the physical tetrahedron space
        p = p1 .+ u .* (p2 .- p1) .+ v .* (p3 .- p1) .+ w .* (p4 .- p1)
        
        # Internal transformation Jacobian scaling factor
        jac_ref = (1.0 - x)^2 * (1.0 - y)
        
        return f(p) * jac_geom * jac_ref
    end
    
    # Integrate over the hypercube [0,0,0] to [1,1,1]
    val, err = hcubature(unit_cube_integrand, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], rtol=rtol, atol=atol)
    return val
end; 

"""
Integrates a scalar function `f(x::Vector)` over a 3D triangle 
defined by vertices p1, p2, and p3.
"""
function integrate_tri(f, p1, p2, p3; rtol=1e-8, atol=1e-12)
    # Calculate the cross product of the two spanning vectors
    v1 = p2 - p1
    v2 = p3 - p1
    cross_prod = cross(v1, v2)
    
    # 2 * Area of the triangle is the norm of the cross product
    twice_area = norm(cross_prod)
    
    # The transformed integrand mapping [0,1]² to the triangle
    # uv[1] is 'u', uv[2] is 'v'
    integrand = function(uv)
        u, v = uv[1], uv[2]
        
        # Parametric coordinate on the 3D triangle plane
        point = p1 + u * (p2 - p1) + u * v * (p3 - p2)
        
        # f(point) multiplied by the Jacobian factor: 2 * Area * u
        return f(point) * twice_area * u
    end
    
    # Integrate over the unit square [0, 1] × [0, 1]
    val, err = hcubature(integrand, [0.0, 0.0], [1.0, 1.0]; rtol=rtol, atol=atol)
    return val
end; 