using BenchmarkTools

const DEFAULT_ALERT_RATIO = 1.30

function _usage()
    return "usage: julia --project=benchmarks benchmarks/compare_benchmarks.jl BASELINE.json TARGET.json [ALERT_RATIO]"
end

length(ARGS) in (2, 3) || error(_usage())

baseline_path, target_path = ARGS[1:2]
alert_ratio = length(ARGS) == 3 ? parse(Float64, ARGS[3]) : DEFAULT_ALERT_RATIO
alert_ratio > 1 || error("alert ratio must be greater than one")

for path in (baseline_path, target_path)
    isfile(path) || error("benchmark result file not found: $path")
end

function _leaf_map(path)
    results = BenchmarkTools.load(path)[1]
    return Dict(
        Tuple(string.(keypath)) => estimate for
        (keypath, estimate) in BenchmarkTools.leaves(results)
    )
end

_ratio(current, baseline) =
    iszero(baseline) ? (iszero(current) ? 1.0 : Inf) : current / baseline
_key_string(keypath) = join(keypath, " / ")
_ratio_string(value) = isfinite(value) ? string(round(value; digits = 3), "x") : "inf"

baseline = _leaf_map(baseline_path)
target = _leaf_map(target_path)

baseline_keys = Set(keys(baseline))
target_keys = Set(keys(target))
common_keys = sort!(collect(intersect(baseline_keys, target_keys)); by = _key_string)
target_only = sort!(collect(setdiff(target_keys, baseline_keys)); by = _key_string)
baseline_only = sort!(collect(setdiff(baseline_keys, target_keys)); by = _key_string)

regressions = Tuple[]
for key in common_keys
    time_ratio = _ratio(target[key].time, baseline[key].time)
    time_ratio > alert_ratio && push!(regressions, key)
end

function _write_report(io)
    println(io, "## Paired benchmark comparison")
    println(io)
    println(
        io,
        "Base and target were measured in the same workflow job. A time ratio above ",
        alert_ratio,
        "x fails the gate, matching the existing 30% alert threshold.",
    )
    println(io)
    println(io, "| Benchmark | Base time | Target time | Time | Memory | Allocs |")
    println(io, "|---|---:|---:|---:|---:|---:|")
    for key in common_keys
        base = baseline[key]
        head = target[key]
        println(
            io,
            "| `",
            _key_string(key),
            "` | ",
            round(base.time; digits = 1),
            " ns | ",
            round(head.time; digits = 1),
            " ns | ",
            _ratio_string(_ratio(head.time, base.time)),
            " | ",
            _ratio_string(_ratio(head.memory, base.memory)),
            " | ",
            _ratio_string(_ratio(head.allocs, base.allocs)),
            " |",
        )
    end

    if !isempty(target_only)
        println(io)
        println(io, "Target-only benchmarks (reported but not gated):")
        for key in target_only
            println(io, "- `", _key_string(key), "`")
        end
    end

    if !isempty(baseline_only)
        println(io)
        println(io, "Base-only benchmarks (reported but not gated):")
        for key in baseline_only
            println(io, "- `", _key_string(key), "`")
        end
    end

    println(io)
    if isempty(regressions)
        println(
            io,
            "**Result: pass.** No common benchmark exceeded the paired regression threshold.",
        )
    else
        println(io, "**Result: fail.** Regressions above ", alert_ratio, "x:")
        for key in regressions
            ratio = _ratio(target[key].time, baseline[key].time)
            println(io, "- `", _key_string(key), "`: ", _ratio_string(ratio))
        end
    end
end

_write_report(stdout)

summary_path = get(ENV, "GITHUB_STEP_SUMMARY", "")
if !isempty(summary_path)
    open(summary_path, "a") do io
        _write_report(io)
    end
end

isempty(regressions) || exit(1)
