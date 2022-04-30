module SIRProject

function SIS(du, u, p, t)
	s, I = u
	β, γ = p
	du[1] = ds = -β*s*I/(s+I) + γ*I
	du[2] = dI = β*s*I/(s+I) - γ*I
end

#="""
Its governed by the following set of equations 
``\frac{dS}{dt} = \frac{-βIs}{s+I+r}``
``\frac{dI}{dt} = \frac{βIs}{s+I+r} - γI``
``\frac{dI}{dt} = γI``=#
"""
``I`` = infected individuals who can pass on the disease to others
``s`` = individuals who're yet to be infected
``r`` = individuals who've been infected but can not transmit the disease.
``\beta`` and ``\gamma`` are positive constants representing the infection rate and the recovery rate. 

"""
function SIR(du, u, p, t)
	s, I, r = u
	β, γ = p
	du[1] = ds = -β*s*I/(s+I+r)
	du[2] = dI = β*s*I/(s+I+r) - γ*I
	du[3] = dr = γ*I
end

function SIRD(du, u, p, t)
	s, I, r, d = u
	β, γ, μ = p
	du[1] = ds = -β*s*I/(s+I+r-d)
	du[2] = dI = β*s*I/(s+I+r-d) - γ*I - μ*I
	du[3] = dr = γ*I
	du[4] = dd = μ*I
end

function SIRV(du, u, p, t)
	s, I, r, V = u
	β, γ, v = p
	du[1] = ds = -β*s*I/(s+I+r+V) - v*s
	du[2] = dI = β*s*I/(s+I+r+V) - γ*I
	du[3] = dr = γ*I
	du[4] = dv = v*s
end

function SEIR(du, u, p, t)
	s, e, I, r = u
	β, γ, a, μ = p
	du[1] = ds = μ*(I+r+e) -β*s*I/(s+I+r+e)
	du[2] = de = β*I*s/(s+I+r+e) - (μ+a)*e
	du[3] = dI = a*e - (γ + μ)*I
	du[4] = dr = γ*I - μ*r
end

export SIS, SIR, SIRV, SIRD, SEIR

end
