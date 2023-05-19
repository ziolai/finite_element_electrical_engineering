"""
    gausslegendre_points(order::Int64) -> Array{Float64,1}
    
Returns coordinates of the Gauss-Legendre quadrature points on the default interval [-1,1] 
for exact integration of polynomials up to the given order.
"""
function gausslegendre_points(order::Int64)
    if order <= 1
        return [0]
    elseif order <= 3
        val = 1/sqrt(3)
        return [-val, val]
    elseif order <= 5
        val = sqrt(3/5)
        return  [-val, 0, val]
    elseif order <= 7
        val = 2/7 * sqrt(6/5)
        val1 = sqrt(3/7 - val)
        val2 = sqrt(3/7 + val)
        return  [-val2, -val1, val1, val2]
    elseif order <= 9
        val = 2 * sqrt(10/7)
        val1 = 1/3 * sqrt(5 - val)
        val2 = 1/3 * sqrt(5 + val)
        return  [-val2, -val1, 0, val1, val2]
    else
        throw(ErrorException("Order $order not possible. " *
                                "Highest possible order for 1D is 9."))
    end
end

"""
    gausslegendre_weights(order::Int64) -> Array{Float64,1}
    
Returns weights of the Gauss-Legendre quadrature points on the default interval [-1,1] 
for exact integration of polynomials up to the given order.
"""
function gausslegendre_weights(order::Int64)
    if order <= 1
        return [2]
    elseif order <= 3
        return [1, 1]
    elseif order <= 5
        val = 5/9
        return  [val, 8/9, val]
    elseif order <= 7
        val = sqrt(30)
        val1 = (18 + val) / 36
        val2 = (18 - val) / 36
        return  [val2, val1, val1, val2]
    elseif order <= 9
        val = 13 * sqrt(70)
        val1 = (322 + val) / 900
        val2 = (322 - val) / 900
        return  [val2, val1, 128/225, val1, val2]
    else
        throw(ErrorException("Order $order not possible. " *
                                "Highest possible order for 1D is 9."))
    end
end

"""
    compute_coordinates_triangle(order::Int64) -> Array{Array{Float64,1},1}
    
Returns coordinates of the Gauss-Legendre quadrature points
on the 2-dimensional FEM reference element
for exact integration of polynomials up to the given order.
"""
function compute_coordinates_triangle(order::Int64)
    order1 = order + 1
    p1 = gausslegendre_points(order1)
    n = length(p1)

    c = Array{Array{Float64,1},1}(undef, n*n)

    r = 1
    for i = 1:n
        for j = 1:n
            c[r] = [(1 + p1[i]) / 2 , (1 - p1[i]) * (1 + p1[j]) / 4]
            r = r + 1
        end
    end
    return c
end

"""
    compute_weights_triangle(order::Int64) -> Array{Float64,1}
    
Returns weights of the Gauss-Legendre quadrature points
on the 2-dimensional FEM reference element
for exact integration of polynomials up to the given order.
"""
function compute_weights_triangle(order::Int64)
    order1 = order + 1
    p1 = gausslegendre_points(order1)
    w1 = gausslegendre_weights(order1)
    n = length(p1)
    
    w = Array{Float64,1}(undef, n*n)

    r = 1
    for i = 1:n
        for j = 1:n
            w[r] = (1 - p1[i]) / 8 * w1[i] * w1[j]
            r = r + 1
        end
    end
    return w
end


