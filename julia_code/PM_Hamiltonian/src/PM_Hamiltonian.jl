using DifferentialEquations
#using LinearAlgebra
using DelimitedFiles

Base.@ccallable function julia_main()::Cint
	try
		real_main()
	catch
		Base.invokelatest(Base.display_error, Base.catch_stack())
		return 1
	end
	return 0
end

function real_main()::Cint
	file = ARGS[1]
	arr = readdlm(file, ' ', Float64, '\n')
	n = size(arr)[1] - 1
	G, M, L, T = arr[1,1], arr[1,2], arr[1,3], arr[1,4]
	q0 = arr[2:end,2:4]
	p0 = arr[2:end,5:7]
	c0 = arr[2:end,1]
	append!(c0,G)
	C = 1.000000000000; #Speed of light
    tspan = (0.0, 1 * 12000.0); # The amount of time for which the simulation runs

	function H(p, q, params)
		paramsCopy = copy(params)
		G = paramsCopy[end]
		m = deleteat!(paramsCopy, length(params))

		return sum(mline(a, m, p) for a in 1:n) + (-0.5 * G * sum((mline(a, m, p) * mline(b, m, p)) / norm(q[a,:] - q[b,:]) * (1 + dot(p[a,:],p[a,:])/(mline(a, m, p)*mline(a, m, p)) + dot(p[b,:],p[b,:])/(mline(b, m, p)*mline(b, m, p))) for a in 1:n, b in 1:n if a != b)) + (0.25 * G * sum(1 / norm(q[a,:] - q[b,:]) * (7 * dot(p[a,:], p[b,:]) + dot(p[a,:], nhat(q, a, b)) * dot(p[b,:], nhat(q, a, b))) for a in 1:n, b in 1:n if a != b)) + (-0.25 * G * sum(1 / (norm(q[a,:] - q[b,:]) * mline(a, m, p) * mline(b, m, p) * (y(m, q, p, b, a) + 1) * (y(m, q, p, b, a) + 1) * y(m, q, p, b, a)) * (2 * (2 * dot(p[a,:], p[b,:]) * dot(p[a,:], p[b,:]) * dot(p[b,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) - 2 * dot(p[a,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) * dot(p[a,:], p[b,:]) * dot(p[b,:], p[b,:]) + dot(p[a,:], nhat(q, b, a)) * dot(p[a,:], nhat(q, b, a)) * dot(p[b,:], p[b,:]) * dot(p[b,:], p[b,:]) - dot(p[a,:], p[b,:]) * dot(p[a,:], p[b,:]) * dot(p[b,:], p[b,:])) * 1 / (mline(b, m, p) * mline(b, m, p)) + 2 * (-dot(p[a,:], p[a,:]) * dot(p[b,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) + dot(p[a,:], nhat(q, b, a)) * dot(p[a,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) + 2 * dot(p[a,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) * dot(p[a,:], p[b,:]) + dot(p[a,:], p[b,:]) * dot(p[a,:], p[b,:]) - dot(p[a,:], nhat(q, b, a)) * dot(p[a,:], nhat(q, b, a)) * dot(p[b,:], p[b,:])) + (-3 * dot(p[a,:], p[a,:]) * dot(p[b,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) + dot(p[a,:], nhat(q, b, a)) * dot(p[a,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) + 8 * dot(p[a,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) * dot(p[a,:], p[b,:]) + dot(p[a,:], p[a,:]) * dot(p[b,:], p[b,:]) - 3 * dot(p[a,:], nhat(q, b, a)) * dot(p[a,:], nhat(q, b, a)) * dot(p[b,:], p[b,:])) * y(m, q, p, b, a)) for a in 1:n, b in 1:n if a != b))
	end

	prob = HamiltonianProblem(H, p0, q0, tspan, c0)
	sol = solve(prob, Vern9(), #=callback=cb_collision,=# reltol = 1.0e-9, abstol = 1.0e-9, saveat = 10);

end

function dot(q, p)
	return (q[1] * p[1]) + (q[2] * p[2]) + (q[3] * p[3])
end

function mline(a, m, p)
	return sqrt(m[a]*m[a] + dot(p[a,:], p[a,:]))
end

function norm(q)
	return sqrt(dot(q, q))
end

function nhat(q, a, b)
	rab = q[a,:] - q[b,:]
	normal = norm(rab)
	return rab / normal
end

function nhat2(q, a, b)
	return normalize(q[a,:] - q[b,:])
end

function y(m, q, p, a, b)
	return (1 / mline(a, m, p)) * sqrt(m[a]*m[a] + dot(nhat(q,a,b), p[a,:])*dot(nhat(q,a,b), p[a,:]))
end

function pmnbody(p, q, params)
	paramsCopy = copy(params)
	G = paramsCopy[end]
	m = deleteat!(paramsCopy, length(params))

	n = length(m)

	w1 = sum(mline(a, m, p) for a in 1:n)

	w2 = -0.5 * G * sum((mline(a, m, p) * mline(b, m, p)) / norm(q[a,:] - q[b,:]) * (1 + dot(p[a,:],p[a,:])/(mline(a, m, p)*mline(a, m, p)) + dot(p[b,:],p[b,:])/(mline(b, m, p)*mline(b, m, p))) for a in 1:n, b in 1:n if a != b)

	w3 = 0.25 * G * sum(1 / norm(q[a,:] - q[b,:]) * (7 * dot(p[a,:], p[b,:]) + dot(p[a,:], nhat(q, a, b)) * dot(p[b,:], nhat(q, a, b))) for a in 1:n, b in 1:n if a != b)

	w4 = -0.25 * G * sum(1 / (norm(q[a,:] - q[b,:]) * mline(a, m, p) * mline(b, m, p) * (y(m, q, p, b, a) + 1) * (y(m, q, p, b, a) + 1) * y(m, q, p, b, a)) * (2 * (2 * dot(p[a,:], p[b,:]) * dot(p[a,:], p[b,:]) * dot(p[b,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) - 2 * dot(p[a,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) * dot(p[a,:], p[b,:]) * dot(p[b,:], p[b,:]) + dot(p[a,:], nhat(q, b, a)) * dot(p[a,:], nhat(q, b, a)) * dot(p[b,:], p[b,:]) * dot(p[b,:], p[b,:]) - dot(p[a,:], p[b,:]) * dot(p[a,:], p[b,:]) * dot(p[b,:], p[b,:])) * 1 / (mline(b, m, p) * mline(b, m, p)) + 2 * (-dot(p[a,:], p[a,:]) * dot(p[b,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) + dot(p[a,:], nhat(q, b, a)) * dot(p[a,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) + 2 * dot(p[a,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) * dot(p[a,:], p[b,:]) + dot(p[a,:], p[b,:]) * dot(p[a,:], p[b,:]) - dot(p[a,:], nhat(q, b, a)) * dot(p[a,:], nhat(q, b, a)) * dot(p[b,:], p[b,:])) + (-3 * dot(p[a,:], p[a,:]) * dot(p[b,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) + dot(p[a,:], nhat(q, b, a)) * dot(p[a,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) + 8 * dot(p[a,:], nhat(q, b, a)) * dot(p[b,:], nhat(q, b, a)) * dot(p[a,:], p[b,:]) + dot(p[a,:], p[a,:]) * dot(p[b,:], p[b,:]) - 3 * dot(p[a,:], nhat(q, b, a)) * dot(p[a,:], nhat(q, b, a)) * dot(p[b,:], p[b,:])) * y(m, q, p, b, a)) for a in 1:n, b in 1:n if a != b)

	return w1 + w2 + w3 + w4
end

julia_main()
