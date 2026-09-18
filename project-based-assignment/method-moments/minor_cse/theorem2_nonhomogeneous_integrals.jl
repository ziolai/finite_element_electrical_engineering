using StaticArrays
using HCubature # provides adaptive numerical integration 

const Point3D = SVector{3,Float64};

function f1(h,eta)
    return 1/eta^2*(sqrt(h^2+eta^2)-h)
end 

function f2(h,eta)
    return sqrt(h^2+eta^2)/(2*eta^2)-h^2/(2*eta^3)*asinh(eta/h)
end 

function f3(h,eta)
    return (2*h^3+(h^2+eta^2)^(3/2)-3*h^2*sqrt(h^2+eta^2))/(3*eta^4)
end

function f4(h,eta)
    return ((2*eta^3-3*eta*h^2)*sqrt(h^2+eta^2)+3*h^4*asinh(eta/h))/(8*eta^5)
end 

function f5(h,eta)
    return (3*(h^2+eta^2)^(2.5)-10*h^2*(h^2+eta^2)^(1.5)+15*h^4*sqrt(h^2+eta^2)-8*h^5)/(15*eta^6)
end 

function I8integrand(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    arg = sqrt(R0^2+s^2)
    return f1(h,arg) 
end

function I8(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    t0 = sqrt(h^2+R0^2)
    return I14(s,h,R0) + asinh(s/t0) - h/R0*atan(s/R0)
end

function I9(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    t0 = sqrt(h^2+R0^2)
    return sqrt(t0^2+s^2)/h-asinh(h/sqrt(R0^2+s^2))-.5*log(R0^2+s^2) 
end

function I9integrand(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    arg = sqrt(R0^2+s^2)
    return s/h*f1(h,arg) 
end

function I10(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    t0 = sqrt(h^2+R0^2)
    return t0^2/(2*R0^2)*asinh(s/t0)-(h^2*s)/(2*R0^2*sqrt(R0^2+s^2))*asinh(sqrt(R0^2+s^2)/h)
end

function I10integrand(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    arg = sqrt(R0^2+s^2)
    return f2(h,arg) 
end

function I11(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    t0 = sqrt(h^2+R0^2)
    return sqrt(t0^2+s^2)/(2*h)+h*asinh(sqrt(R0^2+s^2)/h)/(2*sqrt(R0^2+s^2))
end

function I11integrand(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    arg = sqrt(R0^2+s^2)
    return s/h*f2(h,arg) 
end

function I12(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    return I14(s,h,R0)/3+I15(s,h,R0)/3+I18(s,h,R0)+2*I19(s,h,R0)/3 
end

function I12integrand(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    arg = sqrt(R0^2+s^2)
    return f3(h,arg) 
end

function I13(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    t0 = sqrt(h^2+R0^2)
    return ((t0^2+s^2)^(3/2)/h-h^2)/(3*(R0^2+s^2))
end

function I13integrand(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    arg = sqrt(R0^2+s^2)
    return s/h*f3(h,arg) 
end

function I14(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    t0 = sqrt(h^2+R0^2)
    argplus = sqrt(t0^2+s^2)+h
    argmin  = sqrt(t0^2+s^2)-h
    return sign(s)*(I17plus(argmin,h,R0) - I17min(argplus,h,R0)) 
end

function I15(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    t0 = sqrt(h^2+R0^2)
    return asinh(s/t0)
end

function I16(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    t0 = sqrt(h^2+R0^2)
    return h/R0*atan(s/R0)
end

function I17plus(x,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    t0 = sqrt(h^2+R0^2)
    arg = -R0^2/(x*t0)+h/t0
    return h/(2*R0)*asin(arg) 
end 

function I17min(x,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    t0 = sqrt(h^2+R0^2)
    arg  = -R0^2/(x*t0)-h/t0
    return h/(2*R0)*asin(arg)  
end

function I18(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    t0 = sqrt(h^2+R0^2)
    return h^3/(3*R0^3)*(R0*s/(R0^2+s^2)+atan(s/R0)) 
end

function I19(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    t0 = sqrt(h^2+R0^2)
    value = -(h^2*s*sqrt(t0^2+s^2))/(2*R0^2*(R0^2+s^2))
    value += h/(2*R0)*(1-h^2/R0^2)*atan(h*s/(R0*sqrt(t0^2+s^2)))
    return value-I14(s,h,R0)
end

function I14hat(s,h,R0) 
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    t0 = sqrt(h^2+R0^2)
    t1 = 1/8*(1+h^2/R0^2)*I14(s,h,R0)
    t2 = 1/4*(1-h^4/R0^4)*I15(s,h,R0)
    t3 = 1/4*I19(s,h,R0)
    t4 = h^4/(8*R0^4)*(R0^2*s/(R0^2+s^2)^(1.5)+2*s/sqrt(R0^2+s^2))*asinh(sqrt(R0^2+s^2)/h)
    return t1+t2+t3+t4 
end

function I14hatintegrand(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    arg = sqrt(R0^2+s^2)
    return f4(h,arg) 
end

function I15hat(s,h,R0) 
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    t0 = sqrt(h^2+R0^2)
    return sqrt(t0^2+s^2)/(4*h)+h*sqrt(t0^2+s^2)/(8*(R0^2+s^2))-h^3*asinh(sqrt(R0^2+s^2)/h)/(8*(R0^2+s^2)^(1.5))
end

function I15hatintegrand(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    arg = sqrt(R0^2+s^2)
    return s/h*f4(h,arg) 
end

function I16hat(s,h,R0) 
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    t0 = sqrt(h^2+R0^2)
    t1 = 2*h^2/(15*R0^2)*s*(t0^2+s^2)^(1.5)/(R0^2+s^2)^2 
    t2 = 1/5*asinh(s/t0)
    t3 = (h^4/(5*R0^4)-h^2/(5*R0^2))*s*sqrt(t0^2+s^2)/(R0^2+s^2)
    t4 = (h^5/(5*R0^5)+h/(15*R0))*atan(h*s/(R0*sqrt(t0^2+s^2)))
    t5 = -1/15*I14(s,h,R0)
    t6 = -h^5/(15*R0^5)*(2*R0^3*s/(R0^2+s^2)^2+3*R0*s/(R0^2+s^2)+3*atan(s/R0))
    return t1+t2+t3+t4+t5+t6 
end

function I16hatintegrand(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    arg = sqrt(R0^2+s^2)
    return f5(h,arg) 
end

function I17hat(s,h,R0) 
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    t0 = sqrt(h^2+R0^2)
    return (2*h)/15*(h^3-(t0^2+s^2)^(1.5))/(R0^2+s^2)^2+h/5*sqrt(t0^2+s^2)/(R0^2+s^2)+1/(5*h)*sqrt(t0^2+s^2)
end

function I17hatintegrand(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    arg = sqrt(R0^2+s^2)
    return s/h*f5(h,arg) 
end

function I3(a::Point3D,el::Edge3D)
    cutoff = 1e-12 
    num  = norm(el.r2-a) + norm(el.r1-a) + norm(el.r2-el.r1)
    den  = norm(el.r2-a) + norm(el.r1-a) - norm(el.r2-el.r1)
    if (den<cutoff) den = cutoff end 
    return log(num/den)  
end 

function I3integrand(s,R0)
    arg = sqrt(R0^2+s^2)
    return 1/arg
end

function R(a::Point3D,e::Edge3D,h) 
    R0 = norm(cross(e.r1-a,e.tau))
    s1 = dot(e.r1-a,e.tau)
    s2 = dot(e.r2-a,e.tau)
    return I8(s2,h,R0) - I8(s1,h,R0)   
end

function Rintegrand(s,h,R0)
    cutoff = 1e-12; if (R0<cutoff) R0 = cutoff end
    arg = sqrt(R0^2+s^2)
    return f1(h,arg)
end

function Rphi1(a::Point3D,el::Edge3D)
    t1 = dot(el.r2-a,el.tau)/norm(el.r2-el.r1)*I3(a,el)
    t2 = (norm(el.r1-a) - norm(el.r2-a)) /norm(el.r2-el.r1)
    return t1+t2 
end;

function Rphi1integrand(lmbd,a::Point3D,el::Edge3D)
    r  = el1.r1+lmbd*el1.tau
    return dot(el.r2-r,el.e)/norm(r-a)
end;

function myI8delta(s1,s2,h,R0)   
    val = hcubature(s -> f1(h,sqrt(s[1]^2+R0^2)), (s1,), (s2,))[1]
    return val 
end 

function myI9delta(s1,s2,h,R0)   
    val = hcubature(s -> (s[1]/h)*f1(h,sqrt(s[1]^2+R0^2)), (s1,), (s2,))[1]
    return val 
end 

function myI10delta(s1,s2,h,R0)   
    val = hcubature(s -> f2(h,sqrt(s[1]^2+R0^2)), (s1,), (s2,))[1]
    return val 
end 

function myI11delta(s1,s2,h,R0)   
    val = hcubature(s -> (s[1]/h)*f2(h,sqrt(s[1]^2+R0^2)), (s1,), (s2,))[1]
    return val 
end 

function myI12delta(s1,s2,h,R0)   
    val = hcubature(s -> f3(h,sqrt(s[1]^2+R0^2)), (s1,), (s2,))[1]
    return val 
end 

function myI13delta(s1,s2,h,R0)   
    val = hcubature(s -> (s[1]/h)*f3(h,sqrt(s[1]^2+R0^2)), (s1,), (s2,))[1]
    return val 
end 