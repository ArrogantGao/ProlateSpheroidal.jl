using ProlateSpheroidal
using Plots

# Generate plots for multiple eps values.
const eps_values = [1e-6, 1e-10, 1e-16]
xs = range(-1.0, 1.0; length = 201)

plt_eval = plot(title = "psi0(x) for c = prolc180(eps)", xlabel = "x", ylabel = "psi0(x)")
plt_der = plot(title = "psi0'(x) for c = prolc180(eps)", xlabel = "x", ylabel = "d/dx psi0(x)")

for eps in eps_values
    c = prolc180(eps)
    params = Prolate0Params(c)
    vals = [params(x) for x in xs]
    ders = [prolate0_eval_derivative(params, x) for x in xs]
    label = "eps=$(eps), c=$(round(c; digits=3))"
    plot!(plt_eval, xs, vals; label = label, yscale = :log10, dpi = 300)
    plot!(plt_der, xs, ders; label = label, dpi = 300)
end

mkpath(joinpath(@__DIR__, "figs"))

savefig(plt_eval, joinpath(@__DIR__, "figs", "psi0_eval.png"))
savefig(plt_der, joinpath(@__DIR__, "figs", "psi0_derivative.png"))
