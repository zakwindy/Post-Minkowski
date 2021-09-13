function write_equations()::Cint
    ndata, equations, newton = ARGS[1], ARGS[2], ARGS[3]

    ndata = open(ARGS[1], "r")

    nbody_string = readline(ndata)
    nbody = parse(Int64, nbody_string)
    ndim_string = readline(ndata)
    ndim = parse(Int64, ndim_string)

    close(ndata)

    equations = open(ARGS[2], "r")
    newton = open(ARGS[3], "r")
    output = open("FormattedEquations.jl", "w+")

    write(output, "module julia_app\n")

    write(output, "using DifferentialEquations\n")
    write(output, "using DelimitedFiles\n\n")

    write(output, "Base.@ccallable function julia_main()::Cint\n")
    write(output, "\ttry\n")
    write(output, "\t\treal_main()\n")
    write(output, "\tcatch\n")
    write(output, "\t\tBase.invokelatest(Base.display_error, Base.catch_stack())\n")
    write(output, "\t\treturn 1\n")
    write(output, "\tend\n")
    write(output, "\treturn 0\n")
    write(output, "end\n\n")

    write(output, "function real_main()::Cint\n")
    write(output, "\tG = 1.0\n\tC = 1.0\n\tmSun = 1.0\n")
    write(output, "\tfile = ARGS[1]\n")
    write(output, "\tarr = readdlm(file, ' ', Float64, '\\n')\n")
    write(output, "\tnbody = size(arr)[1]\n")
    write(output, "\tdata = arr[1,2:end]\n")
    write(output, "\tfor i in 2:nbody\n")
    write(output, "\t\tappend!(data, arr[i,2:end])\n")
    write(output, "\tend\n")
    write(output, "\tc0 = arr[1:end,1]\n")
    write(output, "\tappend!(c0,G)\n")
    write(output, "\ttspan = (0.0, 1 * 12000.0); # The amount of time for which the simulation runs\n")
    write(output, "\th01 = [0.0]\n\t#TAYDEN WUZ HERE\n")
    write(output, "\tu0 = collect(Base.Iterators.flatten([data, h01]));\n")
    write(output, "\tprob = ODEProblem(FormattedEquations.PM, u0, tspan, c0);\n")
    write(output, "\tsol = DifferentialEquations.solve(prob, Vern9(), #=callback=cb_collision,=# reltol = 1.0e-9, abstol = 1.0e-9, saveat = 10);\n\n")

    write(output, "\topen(\"PMdata.csv\", \"w+\") do io\n")
    write(output, "\t\twritedlm(io, sol, ',')\n")
    write(output, "\tend\n\n")

    write(output, "\treturn 0\n")
    write(output, "end\n\n")

    write(output, "CSI = 3.00e8;\n")
    write(output, "GSI = 6.647e-11;\n")
    write(output, "MSUN = 1.989e30;\n\n")

    ## This section writes out the Post-Minkowskian n-body function

    write(output, "function PM(du, u, p, t)\n")

    for i in 1:nbody
        u_index = (i - 1) * (2 * ndim)
        write(output, "\tqx", string(i), " = u[", string(u_index + 1), ']', '\n')
        write(output, "\tqy", string(i), " = u[", string(u_index + 2), ']', '\n')
        write(output, "\tqz", string(i), " = u[", string(u_index + 3), ']', '\n')
        write(output, "\tpx", string(i), " = u[", string(u_index + 4), ']', '\n')
        write(output, "\tpy", string(i), " = u[", string(u_index + 5), ']', '\n')
        write(output, "\tpz", string(i), " = u[", string(u_index + 6), ']', '\n')
    end

    write(output, '\n')

    for i in 1:nbody
        write(output, "\tm", string(i), " = p[", string(i), "]\n")
    end

    write(output, "\tG = p[", string(nbody + 1), "]\n\n")

    lastline = ""
    ibody = 0.0
    idim = 1
    while ! eof(equations)
        line = readline(equations)

        if line[end] == '\\'
            line = chop(line)
            lastline = string(lastline, line)
        elseif (lastline != "") || (line[1] == 'd')
            lastline = string(lastline, line)
            lastline = fix_dot(lastline)

            index = findfirst('=', lastline)
            lastline = chop(lastline, head = index - 1)
            if ibody < nbody
                lastline = string("du[", idim, ']', lastline)
                idim += 1
                if mod(idim - 1, ndim) == 0
                    idim += ndim
                    ibody += 0.5
                    if ibody == (nbody / 2.0)
                        idim = ndim + 1
                    end
                end
            else
                lastline = string("u[", (2 * nbody * ndim) + 1, ']', lastline)
            end

            write(output, '\t', lastline, "\n")
            lastline = ""
        elseif line[1] == 'o'
            line = fix_dot(line)
            write(output, '\t', line, "\n")
        end
    end

    write(output, "end\n\n")

    close(equations)

    ## This section writes out the Newtonian n-body Hamiltonian equations

    write(output, "function newton(du, u, p, t)\n")

    for i in 1:nbody
        u_index = (i - 1) * (2 * ndim)
        write(output, "\tqx", string(i), " = u[", string(u_index + 1), ']', '\n')
        write(output, "\tqy", string(i), " = u[", string(u_index + 2), ']', '\n')
        write(output, "\tqz", string(i), " = u[", string(u_index + 3), ']', '\n')
        write(output, "\tpx", string(i), " = u[", string(u_index + 4), ']', '\n')
        write(output, "\tpy", string(i), " = u[", string(u_index + 5), ']', '\n')
        write(output, "\tpz", string(i), " = u[", string(u_index + 6), ']', '\n')
    end

    write(output, '\n')

    for i in 1:nbody
        write(output, "\tm", string(i), " = p[", string(i), "]\n")
    end

    write(output, "\tG = p[", string(nbody + 1), "]\n\n")

    lastline = ""
    ibody = 0.0
    idim = 1
    while ! eof(newton)
        line = readline(newton)

        if line[end] == '\\'
            line = chop(line)
            lastline = string(lastline, line)
        elseif (lastline != "") || (line[1] == 'd')
            lastline = string(lastline, line)
            lastline = fix_dot(lastline)

            index = findfirst('=', lastline)
            lastline = chop(lastline, head = index - 1)
            if ibody < nbody
                lastline = string("du[", idim, ']', lastline)
                idim += 1
                if mod(idim - 1, ndim) == 0
                    idim += ndim
                    ibody += 0.5
                    if ibody == (nbody / 2.0)
                        idim = ndim + 1
                    end
                end
            else
                lastline = string("u[", (2 * nbody * ndim) + 1, ']', lastline)
            end

            write(output, '\t', lastline, "\n")
            lastline = ""
        elseif line[1] == 'o'
            line = fix_dot(line)
            write(output, '\t', line, "\n")
        end
    end

    write(output, "end\n\n")

    close(newton)

    write(output, "end # module\n")

    close(output)

    return 0
end

function fix_dot(str)::String
    str = replace(str, "./" => ".0/")
    str = replace(str, ".)" => ".0)")
    str = replace(str, ".*" => ".0*")
    str = replace(str, ".+" => ".0+")
    str = replace(str, ".-" => ".0-")
    str = replace(str, ',' => '^')
    str = replace(str, "pow" => "")
    return str
end

write_equations()
