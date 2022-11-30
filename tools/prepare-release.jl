using Dates
using OrderedCollections
using TextWrap
using Pkg.Types: destructure, read_project, write_project
using YAML


function update(func, tgt)
    src = tempname(cleanup=true)
    func(src)
    if !success(`cmp --quiet $src $tgt`)
        mv(src, tgt, force=true)
    end
end


function main(version)

    date = Dates.today()
    year = Dates.year(date)

    yaml = YAML.load_file("CITATION.cff", dicttype=OrderedDict{String,Any})
    yaml["version"] = version
    yaml["date-released"] = date
    update("CITATION.cff") do fname
        YAML.write_file(fname, yaml)
    end

    names = []
    authors = []
    for author in yaml["authors"]
        name = author["given-names"] * " " * author["family-names"]
        email = "<" * author["email"] * ">"
        push!(names, name)
        push!(authors, name * " " * email)
    end

    proj = destructure(read_project("Project.toml"))
    proj["version"] = version
    proj["authors"] = authors

    update("Project.toml") do fname
        write_project(proj, fname)
    end

    update("LICENSE.txt") do fname
        open(fname, "w") do io
            print_wrapped(io,
                replace(read("LICENSE.txt", String),
                    r"(?<=Copyright (c) ).*?(?=\n\n)"s
                        => "2020-$year $(join(names, ", "))"),
                width = 78,
                break_long_words = false,
                replace_whitespace = false)
        end
    end

end


if length(ARGS) != 1
    println("usage: julia $PROGRAM_FILE ‹version›")
    exit(1)
else
    main(ARGS[1])
end

