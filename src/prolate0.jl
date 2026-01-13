function prosinin(c::Float64, ts::AbstractVector{Float64}, whts::AbstractVector{Float64},
                  fs::AbstractVector{Float64}, x::Float64)
    rint = 0.0
    derrint = 0.0
    n = length(ts)
    for i in 1:n
        diff = x - ts[i]
        sin_term = sin(c * diff)
        cos_term = cos(c * diff)
        rint += whts[i] * fs[i] * sin_term / diff
        derrint += whts[i] * fs[i] / (diff * diff) * (c * diff * cos_term - sin_term)
    end
    return rint, derrint
end

function prolcoef(rlam::Float64, k::Int, c::Float64)
    d = k * (k - 1)
    d = d / (2 * k + 1) / (2 * k - 1)
    uk = d

    d = (k + 1) * (k + 1)
    d = d / (2 * k + 3)
    d2 = k * k
    d2 = d2 / (2 * k - 1)
    vk = (d + d2) / (2 * k + 1)

    d = (k + 1) * (k + 2)
    d = d / (2 * k + 1) / (2 * k + 3)
    wk = d

    alpha = -c * c * uk
    beta = rlam - k * (k + 1) - c * c * vk
    gamma = -c * c * wk

    return alpha, beta, gamma, uk, vk, wk
end

function prolmatr!(as::AbstractVector{Float64}, bs::AbstractVector{Float64}, cs::AbstractVector{Float64}, n::Int,
                   c::Float64, rlam::Float64, ifsymm::Int, ifodd::Int)
    done = 1.0
    half = done / 2.0
    k = 0

    if ifodd > 0
        for k0 in 1:2:(n + 2)
            k += 1
            alpha, beta, gamma, _, _, _ = prolcoef(rlam, k0, c)
            as[k] = alpha
            bs[k] = beta
            cs[k] = gamma

            if ifsymm != 0
                if k0 > 1
                    as[k] = as[k] / sqrt(k0 - 2 + half) * sqrt(k0 + half)
                end
                cs[k] = cs[k] * sqrt(k0 + half) / sqrt(k0 + half + 2)
            end
        end
    else
        for k0 in 0:2:(n + 2)
            k += 1
            alpha, beta, gamma, _, _, _ = prolcoef(rlam, k0, c)
            as[k] = alpha
            bs[k] = beta
            cs[k] = gamma

            if ifsymm != 0
                if k0 != 0
                    as[k] = as[k] / sqrt(k0 - 2 + half) * sqrt(k0 + half)
                end
                cs[k] = cs[k] * sqrt(k0 + half) / sqrt(k0 + half + 2)
            end
        end
    end
end

function prolql1!(d::AbstractVector{Float64}, e::AbstractVector{Float64})
    n = length(d)
    if n == 1
        return 0
    end

    for i in 1:(n - 1)
        e[i] = e[i + 1]
    end
    e[n] = 0.0

    for l in 1:n
        j = 0
        while true
            m = l
            while m < n
                tst1 = abs(d[m]) + abs(d[m + 1])
                tst2 = tst1 + abs(e[m])
                if tst2 == tst1
                    break
                end
                m += 1
            end

            if m == l
                break
            end
            if j == 30
                return l
            end
            j += 1

            g = (d[l + 1] - d[l]) / (2.0 * e[l])
            r = sqrt(g * g + 1.0)
            g = d[m] - d[l] + e[l] / (g + copysign(r, g))
            s = 1.0
            c = 1.0
            p = 0.0

            for i in (m - 1):-1:l
                f = s * e[i]
                b = c * e[i]
                r = sqrt(f * f + g * g)
                e[i + 1] = r
                if r == 0.0
                    d[i + 1] -= p
                    e[m] = 0.0
                    break
                end
                s = f / r
                c = g / r
                g = d[i + 1] - p
                r = (d[i] - g) * s + 2.0 * c * b
                p = s * r
                d[i + 1] = g + p
                g = c * r - b
            end

            if r == 0.0
                break
            end
            d[l] -= p
            e[l] = g
            e[m] = 0.0
        end

        if l == 1
            continue
        end
        for i in l:-1:2
            if d[i] >= d[i - 1]
                break
            end
            d[i], d[i - 1] = d[i - 1], d[i]
        end
    end

    return 0
end

function prolfact!(a::AbstractVector{Float64}, b::AbstractVector{Float64}, c::AbstractVector{Float64},
                   u::AbstractVector{Float64}, v::AbstractVector{Float64}, w::AbstractVector{Float64})
    n = length(a)
    for i in 1:(n - 1)
        d = c[i + 1] / a[i]
        a[i + 1] -= b[i] * d
        u[i] = d
    end

    for i in (n - 1):-1:1
        d = b[i] / a[i + 1]
        v[i + 1] = d
    end

    done = 1.0
    for i in 1:n
        w[i] = done / a[i]
    end
end

function prolsolv!(u::AbstractVector{Float64}, v::AbstractVector{Float64}, w::AbstractVector{Float64},
                   rhs::AbstractVector{Float64})
    n = length(rhs)
    for i in 1:(n - 1)
        rhs[i + 1] -= u[i] * rhs[i]
    end

    for i in (n - 1):-1:1
        rhs[i] -= rhs[i + 1] * v[i + 1]
    end

    for i in 1:n
        rhs[i] *= w[i]
    end
end

function prolfun0!(n::Int, c::Float64, as::AbstractVector{Float64}, bs::AbstractVector{Float64},
                   cs::AbstractVector{Float64}, xk::AbstractVector{Float64},
                   u::AbstractVector{Float64}, v::AbstractVector{Float64}, w::AbstractVector{Float64},
                   eps::Float64)
    delta = 1.0e-8
    ifsymm = 1
    numit = 4
    rlam = 0.0
    ifodd = -1

    prolmatr!(as, bs, cs, n, c, rlam, ifsymm, ifodd)

    halfn = n ÷ 2
    bs_half = @view bs[1:halfn]
    as_half = @view as[1:halfn]

    ierr = prolql1!(bs_half, as_half)
    if ierr != 0
        return 2048, 0, 0.0
    end

    rkhi = -bs_half[halfn]
    rlam = -bs_half[halfn] + delta

    fill!(xk, 1.0)

    prolmatr!(as, bs, cs, n, c, rlam, ifsymm, ifodd)

    u_half = @view u[1:halfn]
    v_half = @view v[1:halfn]
    w_half = @view w[1:halfn]
    cs_half = @view cs[1:halfn]
    prolfact!(bs_half, cs_half, as_half, u_half, v_half, w_half)

    for _ in 1:numit
        xk_half = @view xk[1:halfn]
        prolsolv!(u_half, v_half, w_half, xk_half)

        d = 0.0
        for j in 1:halfn
            d += xk_half[j] * xk_half[j]
        end
        d = sqrt(d)
        for j in 1:halfn
            xk_half[j] /= d
        end

        for j in 1:halfn
            as_half[j] = xk_half[j]
        end
    end

    nterms = 0
    half = 0.5
    xk_half = @view xk[1:halfn]
    for i in 1:halfn
        if abs(xk_half[i]) > eps
            nterms = i
        end
        xk_half[i] *= sqrt((i - 1) * 2 + half)
        cs_half[i] = xk_half[i]
    end

    j = 1
    for i in 1:(nterms + 1)
        xk[j] = cs[i]
        xk[j + 1] = 0.0
        j += 2
    end

    nterms *= 2
    return 0, nterms, rkhi
end

function prolps0i!(c::Float64, w::AbstractVector{Float64}, lenw::Int)
    ns = [48, 64, 80, 92, 106, 120, 130, 144, 156, 168, 178, 190, 202, 214, 224, 236, 248, 258, 268, 280]
    eps = 1.0e-16

    n = trunc(Int, c * 3)
    n ÷= 2

    i = floor(Int, c / 10)
    if i <= 19
        n = ns[i + 1]
    end

    ixk = 1
    lxk = n + 2

    ias = ixk + lxk
    las = n + 2

    ibs = ias + las
    lbs = n + 2

    ics = ibs + lbs
    lcs = n + 2

    iu = ics + lcs
    lu = n + 2

    iv = iu + lu
    lv = n + 2

    iw = iv + lv
    lw = n + 2

    ltot = iw + lw
    if ltot >= lenw
        return 512, 0, ltot, 0.0
    end

    as = @view w[ias:(ias + las - 1)]
    bs = @view w[ibs:(ibs + lbs - 1)]
    cs = @view w[ics:(ics + lcs - 1)]
    xk = @view w[ixk:(ixk + lxk - 1)]
    u = @view w[iu:(iu + lu - 1)]
    v = @view w[iv:(iv + lv - 1)]
    ww = @view w[iw:(iw + lw - 1)]

    ierr, nterms, rkhi = prolfun0!(n, c, as, bs, cs, xk, u, v, ww, eps)
    if ierr != 0
        return ierr, nterms, ltot, rkhi
    end

    return 0, nterms, ltot, rkhi
end

function prol0ini!(c::Float64, w::Vector{Float64}, lenw::Int)
    ier = 0
    thresh = 45.0
    iw = 11
    w[1] = iw + 0.1
    w[9] = thresh

    subw = @view w[iw:lenw]
    ier, nterms, ltot, rkhi = prolps0i!(c, subw, lenw - iw + 1)
    if ier != 0
        return ier, 0.0, rkhi, 0, ltot
    end

    if c >= thresh
        w[8] = c
        w[5] = nterms + 0.1
        keep = nterms + 3
        return 0, 0.0, rkhi, keep, ltot
    end

    ngauss = nterms * 2
    lw = nterms + 2
    its = iw + lw
    lts = ngauss + 2
    iwhts = its + lts
    lwhts = ngauss + 2
    ifs = iwhts + lwhts
    lfs = ngauss + 2

    keep = ifs + lfs
    if keep > ltot
        ltot = keep
    end
    if keep >= lenw
        return 1024, 0.0, rkhi, keep, ltot
    end

    w[2] = its + 0.1
    w[3] = iwhts + 0.1
    w[4] = ifs + 0.1

    itype = 1
    xs, ws, _, _ = legeexps(itype, ngauss)
    for i in 0:(ngauss - 1)
        w[its + i] = xs[i + 1]
        w[iwhts + i] = ws[i + 1]
    end

    for i in 0:(ngauss - 1)
        val = legeexev(w[its + i - 1], @view(w[iw:(iw + nterms - 1)]), nterms - 1)
        w[ifs + i - 1] = val
    end

    x0 = 0.0
    f0 = legeexev(x0, @view(w[iw:(iw + nterms - 1)]), nterms - 1)
    rlam, _ = prosinin(c, @view(w[its:(its + ngauss - 1)]), @view(w[iwhts:(iwhts + ngauss - 1)]),
                       @view(w[ifs:(ifs + ngauss - 1)]), x0)

    rlam = rlam / f0
    rlam20 = rlam

    w[5] = nterms + 0.1
    w[6] = ngauss + 0.1
    w[7] = rlam
    w[8] = c

    return 0, rlam20, rkhi, keep, ltot
end

function prol0eva(x::Float64, w::Vector{Float64})
    iw = Int(trunc(w[1]))
    its = Int(trunc(w[2]))
    iwhts = Int(trunc(w[3]))
    ifs = Int(trunc(w[4]))

    nterms = Int(trunc(w[5]))
    ngauss = Int(trunc(w[6]))
    rlam = w[7]
    c = w[8]
    thresh = w[9]

    if abs(x) > 1
        if c < thresh - 1.0e-10
            return 0.0, 0.0
        end

        psi0, derpsi0 = prosinin(c, @view(w[its:(its + ngauss - 1)]), @view(w[iwhts:(iwhts + ngauss - 1)]),
                                @view(w[ifs:(ifs + ngauss - 1)]), x)
        psi0 /= rlam
        derpsi0 /= rlam
        return psi0, derpsi0
    end

    psi0, derpsi0 = legeFDER(x, @view(w[iw:(iw + nterms - 2)]), nterms - 2)
    psi0 = sqrt(2.0) * psi0
    derpsi0 = sqrt(2.0) * derpsi0
    return psi0, derpsi0
end

function prol0int0r(w::Vector{Float64}, r::Float64, xs::Vector{Float64}, ws::Vector{Float64})
    val = 0.0
    for i in eachindex(xs)
        xs_r = (xs[i] + 1) * r / 2
        fval, _ = prol0eva(xs_r, w)
        val += ws[i] * r / 2 * fval
    end
    return val
end

struct Prolate0Params
    c::Float64
    lenw::Int
    keep::Int
    ltot::Int
    workarray::Vector{Float64}
    rlam20::Float64
    rkhi::Float64
    quad_nodes::Vector{Float64}
    quad_weights::Vector{Float64}
    int_nodes::Vector{Float64}
    int_weights::Vector{Float64}
end

function Prolate0Params(c::Real; lenw::Int = 10000, quad_npts::Int = 200, int_npts::Int = 200)
    c_val = Float64(c)
    workarray = zeros(Float64, lenw)
    ier, rlam20, rkhi, keep, ltot = prol0ini!(c_val, workarray, lenw)
    if ier != 0
        error("Unable to init Prolate0Params")
    end
    quad_nodes, quad_weights = gaussian_quadrature(quad_npts)
    int_nodes, int_weights, _, _ = legeexps(1, int_npts)
    return Prolate0Params(c_val, lenw, keep, ltot, workarray, rlam20, rkhi,
                          quad_nodes, quad_weights, int_nodes, int_weights)
end

function prolate0_eval_derivative(params::Prolate0Params, x::Real)
    _, der = prol0eva(Float64(x), params.workarray)
    return der
end

function prolate0_eval(params::Prolate0Params, x::Real)
    val, _ = prol0eva(Float64(x), params.workarray)
    return val
end

(params::Prolate0Params)(x::Real) = prolate0_eval(params, x)

function prolate0_int_eval(params::Prolate0Params, r::Real)
    return prol0int0r(params.workarray, Float64(r), params.int_nodes, params.int_weights)
end

function prolate0_lambda(params::Prolate0Params)
    lambda = 0.0
    for i in eachindex(params.quad_nodes)
        x = params.quad_nodes[i]
        lambda += params.quad_weights[i] * prolate0_eval(params, x) * cos(params.c * x * 0.5)
    end
    lambda /= prolate0_eval(params, 0.5)
    return lambda
end

function prolate0_eval_derivative(c::Real, x::Real)
    return prolate0_eval_derivative(Prolate0Params(c), x)
end

function prolate0_eval(c::Real, x::Real)
    return prolate0_eval(Prolate0Params(c), x)
end

function prolate0_int_eval(c::Real, r::Real)
    return prolate0_int_eval(Prolate0Params(c), r)
end

function prolate0_lambda(c::Real)
    return prolate0_lambda(Prolate0Params(c))
end
