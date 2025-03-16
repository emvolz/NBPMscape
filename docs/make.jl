using NBPMscape
using Documenter

DocMeta.setdocmeta!(NBPMscape, :DocTestSetup, :(using NBPMscape); recursive=true)

makedocs(;
    modules=[NBPMscape],
    authors="Erik Volz <erik.volz@gmail.com>",
    repo="https://github.com/emvolz/NBPMscape.jl/blob/{commit}{path}#{line}",
    sitename="NBPMscape.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://emvolz.github.io/NBPMscape.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API Reference" => "api.md",
        "Examples" => "examples.md"
    ],
)

deploydocs(;
    repo="github.com/emvolz/NBPMscape.jl",
    devbranch="main",
)
