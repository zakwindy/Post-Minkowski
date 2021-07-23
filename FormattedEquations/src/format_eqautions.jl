function write_equations()::Cint
    ndata, equations = ARGS[1], ARGS[2]

    ndata = open(ARGS[1], "r")

    nbody_string = readline(ndata)
    nbody = parse(Int64, nbody_string)
    ndim_string = readline(ndata)
    ndim = parse(Int64, ndim_string)

    close(ndata)

    equations = open(ARGS[2], "r")
    output = open("FormattedEquations.jl", "w+")

    write(output, "module FormattedEquations\n")
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
        elseif lastline != ""
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

    write(output, "end\nend\n")

    close(equations)
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
