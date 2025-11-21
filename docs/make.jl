using ProlateSpheroidal
using Documenter

DocMeta.setdocmeta!(ProlateSpheroidal, :DocTestSetup, :(using ProlateSpheroidal); recursive=true)

makedocs(;
    modules=[ProlateSpheroidal],
    authors="Xuanzhao Gao <xgao@flatironinstitute.org> and contributors",
    sitename="ProlateSpheroidal.jl",
    format=Documenter.HTML(;
        canonical="https://ArrogantGao.github.io/ProlateSpheroidal.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ArrogantGao/ProlateSpheroidal.jl",
    devbranch="main",
)
