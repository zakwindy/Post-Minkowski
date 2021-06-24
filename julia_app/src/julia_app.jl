module julia_app

using DifferentialEquations
using Plots
using LSODA

function julia_main()::Cint
    @show typeof(ARGS)
    @show ARGS[1]
    return 0
end

@userplot TwoBodyPlot
@recipe function f(tb::TwoBodyPlot)
    m, x, y, z, xarray, yarray, zarray, negx, posx, negy, posy, negz, posz, r = tb.args
    n = r
    xlims --> (negx, posx)
    ylims --> (negy, posy)
    zlims --> (negz, posz)
    append!(xarray, x)
    append!(yarray, y)
    append!(zarray, z)
    if size(xarray)[1] > n
        deleteat!(xarray, 1)
        deleteat!(yarray, 1)
        deleteat!(zarray, 1)
    end
    linewidth --> range(0, m * 5, length = n)
    #seriesalpha --> range(0, 1, length = n)
    aspect_ratio --> 1
    label --> false
    xarray, yarray, zarray
end

CSI = 3.00e8;
GSI = 6.647e-11;
MSUN = 1.989e30;

end # module
