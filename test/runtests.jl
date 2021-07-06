include("tests.jl")

runtests(slow_tests = !("fast" in ARGS))
