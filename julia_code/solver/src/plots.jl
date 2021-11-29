using Plots
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
	arr = readdlm(file, ',', Float64, '\n')
	nbody = (size(arr)[1] - 1) / 9

	plot(arr[3,:], arr[1,:], arr[2,:], aspect_ratio = :equal)

	for i in 2:nbody
		ifactor = convert(Int64, (i - 1) * 6)
		plt = plot!(arr[ifactor + 3,:], arr[ifactor + 1,:], arr[ifactor + 2,:], aspect_ratio = :equal)
	end

	savefig("nOrbits.png")

	hamilArr = arr[size(arr)[1],:]
	deleteat!(hamilArr, 1)
	originalHamiltonian = copy(hamilArr[1])
	lenHamil = size(hamilArr)
	originalArr = ones(lenHamil) * originalHamiltonian
	hamilVariance = hamilArr - originalArr

	plot(hamilVariance)
	savefig("hamiltonian_variance.png")

	return 0
end

julia_main()
