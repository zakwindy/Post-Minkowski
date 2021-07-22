function write_equations()::Cint
    ndata, equations = ARGS[1], ARGS[2]

    ndata = open(ARGS[1], "r")

    nbody_string = readline(ndata)
    nbody = parse(Int64, nbody_string)
    ndim_string = readline(ndata)
    ndim = parse(Int64, ndim_string)

    close(ndata)


    equations = open(ARGS[2], "r")
    output = open("formatted_equations.jl", "w+")

    lastline = ""
    while ! eof(equations)
        line = readline(equations)

        if line[end] == '\\'
            line = chop(line)
            lastline = string(lastline, line)
        elseif lastline != ""
            lastline = string(lastline, line)
            write(output, lastline, "\n")
            lastline = ""
        elseif line[1] == 'o'
            write(output, line, "\n")
        end

    end

    close(equations)
    close(output)

    return 0
end

write_equations()
