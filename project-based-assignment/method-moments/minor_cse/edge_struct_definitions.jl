Base.@kwdef struct Edge3D 
    r1::Point3D
    r2::Point3D
    len::Float64 = norm(r2-r1)
    tau::Point3D = normalize(r2-r1)
    e::Point3D = tau/len 
end 

function (edge3D::Edge3D)(r)
    return [dot(edge3D.r2-r,edge3D.e), dot(r-edge3D.r1,edge3D.e)] 
end

################################################################

Base.@kwdef struct EdgePair
    e1::Edge3D   # first edge in pair 
    e2::Edge3D   # second edge in pair  
end 

function distanceEdgePair(edgePair::EdgePair)
    e1 = edgePair.e1 
    e2 = edgePair.e2  
    r11_hat = e2.r1 + dot(e1.r1-e2.r2,e2.tau)*e2.tau 
    dist = norm(e1.r1-r11_hat) 
    return dist
end 

function directionEdgePair(edgePair::EdgePair)
    e1 = edgePair.e1 
    e2 = edgePair.e2  
    dir = norm(cross(e1.tau,e2.tau)) 
    return dir 
end 

# return the intersection of two edges assumed to exist  
function intersectionEdgePair(edgePair::EdgePair)
    e1 = edgePair.e1 
    e2 = edgePair.e2    
    rdiff = e2.r1-e1.r1
    num1 = dot(e1.tau,e2.tau)*dot(rdiff,e2.tau)-dot(rdiff,e1.tau)
    num2 = -dot(e1.tau,e2.tau)*dot(rdiff,e1.tau)+dot(rdiff,e2.tau)
    den  = dot(e1.tau,e2.tau)^2-1
    return [e1.r1 + (num1/den)*e1.tau, e2.r1 + (num2/den)*e2.tau]
end

# abstract type to hold a pair of edges 
abstract type AbstractEdgePair end 

# intersecting pair of edges 
Base.@kwdef struct InterEdgePair <: AbstractEdgePair
    edgePair::EdgePair
    svec::Point3D = intersectionEdgePair(edgePair)[1] # point of intersection 
end 

# crossing pair of edges 
Base.@kwdef struct CrossEdgePair <: AbstractEdgePair
    edgePair::EdgePair
    e1 = edgePair.e1                                # first edge 
    e2 = edgePair.e2                                # second edge 
    nw::Point3D = normalize(cross(e1.tau,e2.tau))   # normal on the plane
    h::Float64 = abs(dot(nw,e1.r1-e2.r1))           # height between the edges 
    r1hat::Point3D = e1.r1+dot(e2.r1-e1.r1,nw)*nw   # projection of start of first edge onto plane 
    r2hat::Point3D = e1.r2+dot(e2.r2-e1.r2,nw)*nw   # projection of end of first edge onto plane
    e1hat::Edge3D = Edge3D(r1 = r1hat,r2 = r2hat)   # projected first edge 
    edgePairhat::EdgePair = EdgePair(e1=e1hat, e2=e2)    # projected edge pair 
    svec::Point3D = intersectionEdgePair(edgePairhat)[1] # point of intersection 
    R0::Float64 = norm(cross(e2.tau,(e2.r1-r2hat))) # distance 
end

# parallel pair of edges
Base.@kwdef struct ParalEdgePair <: AbstractEdgePair 
    edgePair::EdgePair
end

# coinciding or inline pair of edges    
Base.@kwdef struct CoinEdgePair <: AbstractEdgePair 
    edgePair::EdgePair
end

AllEdgePair = Union{InterEdgePair, CrossEdgePair, ParalEdgePair, CoinEdgePair}

function classifyEdgePair(edgePair::EdgePair)::AllEdgePair 
    edgePair = abstractEdgePair.EdgePair 
    dist = distEdgePair(edgePair) 
    dir = dirEdgePair(edgePair)
    small = 1e-11; 
    if ((dist<small)&&(dir<small))
        abstractEdgePair = InterEdgePair(edgePair = edgePair)
    elseif ((dist>=small)&&(dir<small))
        abstractEdgePair = CrossEdgePair(edgePair = edgePair)
    elseif ((dist<small)&&(dir>=small))
        abstractEdgePair = ParalEdgePair(edgePair = edgePair)
    else
        abstractEdgePair = CoinEdgePair(edgePair = edgePair)
    end
    return abstractEdgePair 
end 

################################################################