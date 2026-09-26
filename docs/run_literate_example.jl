using Literate

length(ARGS) == 2 || error("usage: julia run_literate_example.jl INPUT.jl OUTPUT_DIR")

input_file = abspath(ARGS[1])
output_dir = abspath(ARGS[2])
repo_root = normpath(joinpath(@__DIR__, ".."))

# Keep GR headless in isolated example processes as well as in the parent docs build.
get!(ENV, "GKSwstype", "100")

extra_literate_config = if isempty(get(ENV, "CI", ""))
    Dict("repo_root_path" => repo_root, "repo_root_url" => "file://" * repo_root)
else
    Dict()
end

function postprocess(content)
    # With execute=true, Literate would normally leave `#hide` markers for
    # Documenter to consume. These pages contain pre-executed ordinary Markdown
    # fences, so remove the hidden source lines here instead.
    return replace(content, r"(?m)^[^\n]*#\s*hide[ \t]*(?:\n|$)" => "")
end

function preprocess(content)
    # Configure Plots once inside the isolated Literate sandbox. The `#hide`
    # markers keep this worker-only setup out of the generated page.
    occursin(r"(?m)^using .*Plots", content) || return content
    return """
    using Plots #hide
    gr() #hide
    default(; fmt = :png) #hide
    nothing #hide

    $content
    """
end

Literate.markdown(
    input_file,
    output_dir;
    flavor = Literate.DefaultFlavor(),
    config = extra_literate_config,
    execute = true,
    preprocess = preprocess,
    postprocess = postprocess,
)
