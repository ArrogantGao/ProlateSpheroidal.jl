# Generate reference data using the C++ PSWF library.
# Run from the package root: julia --project test/fixtures/generate_reference.jl

function pswf_lib_path()
    base_dir = joinpath(homedir(), "codes", "PSWF", "build")
    if Sys.iswindows()
        return joinpath(base_dir, "PSWF.dll")
    elseif Sys.isapple()
        return joinpath(base_dir, "libPSWF.dylib")
    else
        return joinpath(base_dir, "libPSWF.so")
    end
end

lib_path = pswf_lib_path()
if !isfile(lib_path)
    error("Reference PSWF library not found at $(lib_path)")
end

c_prolc180(eps) = ccall((:prolc180, lib_path), Float64, (Float64,), eps)
c_prolc180_der3(eps) = ccall((:prolc180_der3, lib_path), Float64, (Float64,), eps)
c_prolate0_eval(c, x) = ccall((:prolate0_eval, lib_path), Float64, (Float64, Float64), c, x)
c_prolate0_eval_derivative(c, x) = ccall((:prolate0_eval_derivative, lib_path), Float64, (Float64, Float64), c, x)
c_prolate0_int_eval(c, r) = ccall((:prolate0_int_eval, lib_path), Float64, (Float64, Float64), c, r)
c_prolate0_lambda(c) = ccall((:prolate0_lambda, lib_path), Float64, (Float64,), c)

c_values = [1.0, 5.0, 20.0, 50.0]
xs_eval = [-0.9, -0.5, 0.0, 0.25, 0.75, 1.0]
xs_der = [-0.8, -0.3, 0.1, 0.6, 0.9]
rs_int = [0.1, 0.5, 1.0]
eps_vals = [1e-6, 1e-10, 1e-16]
xs_from_eps = [-0.5, 0.0, 0.5]

reference = Dict{Symbol, Any}()
reference[:prolc180] = [(eps = eps, value = c_prolc180(eps)) for eps in eps_vals]
reference[:prolc180_der3] = [(eps = eps, value = c_prolc180_der3(eps)) for eps in eps_vals]
reference[:prolate0_eval] = [(c = c, x = x, value = c_prolate0_eval(c, x)) for c in c_values for x in xs_eval]
reference[:prolate0_eval_derivative] = [(c = c, x = x, value = c_prolate0_eval_derivative(c, x)) for c in c_values for x in xs_der]
reference[:prolate0_int_eval] = [(c = c, r = r, value = c_prolate0_int_eval(c, r)) for c in c_values for r in rs_int]
reference[:prolate0_lambda] = [(c = c, value = c_prolate0_lambda(c)) for c in c_values]
reference[:prolate0_eval_from_eps] = [(eps = eps, c = c_prolc180(eps), x = x,
    value = c_prolate0_eval(c_prolc180(eps), x)) for eps in eps_vals for x in xs_from_eps]

out_path = joinpath(@__DIR__, "pswf_reference.jl")
open(out_path, "w") do io
    println(io, "# Auto-generated reference data from the C++ PSWF library.")
    println(io, "# Regenerate with: julia --project test/fixtures/generate_reference.jl")
    println(io, "const PSWF_REFERENCE = ")
    show(io, reference)
    println(io)
end

println("Wrote reference data to $(out_path)")
