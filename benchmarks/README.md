# Benchmarks

`runbenchmarks.jl` executes the package benchmark suite and writes minimum-time
`BenchmarkTools` estimates to `benchmarks_output.json`.

For pull requests, CI benchmarks both GitHub's synthetic merge commit and the
exact current base SHA in the same workflow job. `compare_benchmarks.jl` compares
benchmark leaves present in both results and fails when target time exceeds base
time by more than 30%, matching the existing benchmark alert threshold. New or
removed benchmark leaves are reported but are not gated because there is no
paired measurement for them.

For `main`, the workflow continues to publish accepted benchmark results through
`github-action-benchmark` so the historical charts remain available. Historical
results are not used as the pull-request regression gate because hosted-runner
and dependency drift can otherwise dominate the comparison.
