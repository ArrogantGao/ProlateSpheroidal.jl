function legepol(x::Float64, n::Int)
    if n == 0
        return 1.0, 0.0
    elseif n == 1
        return x, 1.0
    end

    pk = 1.0
    pkp1 = x
    for k in 1:(n - 1)
        pkm1 = pk
        pk = pkp1
        pkp1 = ((2 * k + 1) * x * pk - k * pkm1) / (k + 1)
    end
    pol = pkp1
    der = n * (x * pkp1 - pk) / (x * x - 1)
    return pol, der
end

function legetayl(pol::Float64, der::Float64, x::Float64, h::Float64, n::Int, k::Int)
    done = 1.0
    q0 = pol
    q1 = der * h
    q2 = (2 * x * der - n * (n + done) * pol) / (1 - x * x)
    q2 = q2 * h * h / 2

    sum_val = q0 + q1 + q2
    sumder = q1 / h + q2 * 2 / h

    if k <= 2
        return sum_val, sumder
    end

    qi = q1
    qip1 = q2
    for i in 1:(k - 2)
        d = 2 * x * (i + 1) * (i + 1) / h * qip1 - (n * (n + done) - i * (i + 1)) * qi
        d = d / (i + 1) / (i + 2) * h * h / (1 - x * x)
        qip2 = d

        sum_val += qip2
        sumder += d * (i + 2) / h

        qi = qip1
        qip1 = qip2
    end

    return sum_val, sumder
end

function legerts(itype::Int, n::Int)
    k = 30
    d = 1.0
    d2 = d + 1.0e-24
    if d2 != d
        k = 54
    end

    ts = zeros(Float64, n)
    whts = zeros(Float64, n)

    half = n ÷ 2
    ifodd = n - 2 * half
    h = (atan(1.0) * 4.0) / (2.0 * n)

    ii = 0
    for i in 1:n
        if i < (n ÷ 2 + 1)
            continue
        end
        ii += 1
        t = (2.0 * i - 1.0) * h
        ts[ii] = -cos(t)
    end

    pol, der = legepol(0.0, n)
    x0 = 0.0
    x1 = ts[1]

    n2 = (n + 1) ÷ 2
    pol3 = pol
    der3 = der

    for kk in 1:n2
        if (ifodd == 1) && (kk == 1)
            ts[kk] = x0
            if itype > 0
                whts[kk] = der
            end
            x0 = x1
            if kk < n2
                x1 = ts[kk + 1]
            end
            pol3 = pol
            der3 = der
            continue
        end

        ifstop = 0
        for _ in 1:10
            hh = x1 - x0
            pol, der = legetayl(pol3, der3, x0, hh, n, k)
            x1 = x1 - pol / der

            if abs(pol) < 1.0e-12
                ifstop += 1
            end
            if ifstop == 3
                break
            end
        end

        ts[kk] = x1
        if itype > 0
            whts[kk] = der
        end

        x0 = x1
        if kk < n2
            x1 = ts[kk + 1]
        end
        pol3 = pol
        der3 = der
    end

    for i in n2:-1:1
        ts[i + half] = ts[i]
    end
    for i in 1:half
        ts[i] = -ts[n - i + 1]
    end
    if itype <= 0
        return ts, whts
    end

    for i in n2:-1:1
        whts[i + half] = whts[i]
    end
    for i in 1:half
        whts[i] = whts[n - i + 1]
    end

    for i in 1:n
        tmp = 1.0 - ts[i] * ts[i]
        whts[i] = 2.0 / tmp / (whts[i] * whts[i])
    end

    return ts, whts
end

function legepols(x::Float64, n::Int)
    pols = zeros(Float64, n + 1)
    if n == 0
        pols[1] = 1.0
        return pols
    elseif n == 1
        pols[1] = 1.0
        pols[2] = x
        return pols
    end

    pols[1] = 1.0
    pols[2] = x
    pkm1 = 1.0
    pk = x
    for k in 1:(n - 1)
        pkp1 = ((2 * k + 1) * x * pk - k * pkm1) / (k + 1)
        pols[k + 2] = pkp1
        pkm1 = pk
        pk = pkp1
    end

    return pols
end

function legeexps(itype::Int, n::Int)
    itype_rts = (itype > 0) ? 1 : 0
    x, whts = legerts(itype_rts, n)

    if itype != 2
        return x, whts, nothing, nothing
    end

    u = zeros(Float64, n, n)
    v = zeros(Float64, n, n)

    for i in 1:n
        pols = legepols(x[i], n - 1)
        for j in 1:n
            u[j, i] = pols[j]
        end
    end

    for i in 1:n
        for j in 1:n
            v[i, j] = u[j, i]
        end
    end

    for i in 1:n
        d = (2 * i - 1) / 2
        for j in 1:n
            u[i, j] = v[j, i] * whts[j] * d
        end
    end

    return x, whts, u, v
end

function legeexev(x::Float64, pexp::AbstractVector{Float64}, n::Int)
    pjm2 = 1.0
    pjm1 = x
    val = pexp[1] * pjm2 + pexp[2] * pjm1

    for j in 2:n
        pj = ((2 * j - 1) * x * pjm1 - (j - 1) * pjm2) / j
        val += pexp[j + 1] * pj
        pjm2 = pjm1
        pjm1 = pj
    end

    return val
end

function legeFDER(x::Float64, pexp::AbstractVector{Float64}, n::Int)
    pjm2 = 1.0
    pjm1 = x
    derjm2 = 0.0
    derjm1 = 1.0

    val = pexp[1] * pjm2 + pexp[2] * pjm1
    der = pexp[2]

    for j in 2:n
        pj = ((2 * j - 1) * x * pjm1 - (j - 1) * pjm2) / j
        val += pexp[j + 1] * pj

        derj = (2 * j - 1) * (pjm1 + x * derjm1) - (j - 1) * derjm2
        derj /= j
        der += pexp[j + 1] * derj

        pjm2 = pjm1
        pjm1 = pj
        derjm2 = derjm1
        derjm1 = derj
    end

    return val, der
end

function legendre(n::Int, x::Float64)
    if n == 0
        return 1.0, 0.0
    elseif n == 1
        return x, 1.0
    end

    pn_minus1 = 1.0
    pn_minus2 = 0.0
    pn = x

    for k in 2:n
        pn_minus2 = pn_minus1
        pn_minus1 = pn
        pn = ((2.0 * k - 1.0) * x * pn_minus1 - (k - 1.0) * pn_minus2) / k
    end

    pn_prime = n * (x * pn - pn_minus1) / (x * x - 1.0)
    return pn, pn_prime
end

function gaussian_quadrature(n::Int)
    nodes, weights = legerts(1, n)
    return copy(nodes), copy(weights)
end

function pseudo_inv(M::Matrix{Float64}, eps::Float64)
    if isempty(M)
        return zeros(Float64, size(M, 2), size(M, 1))
    end

    U, S, V = svd(M)
    max_S = abs(S[1])
    eps_ = max_S * eps

    S_inv = similar(S)
    for i in eachindex(S)
        if S[i] < eps_
            S_inv[i] = 0.0
        else
            S_inv[i] = 1.0 / S[i]
        end
    end

    return V * Diagonal(S_inv) * U'
end

function monomial_nodes_1d(nnodes::Int; a::Float64 = 0.0, b::Float64 = 1.0)
    nodes = zeros(Float64, nnodes)
    for i in 0:(nnodes - 1)
        nodes[i + 1] = i / (nnodes - 1) * (b - a) + a
    end
    return nodes
end

function monomial_basis_1d(order::Int, x::Vector{Float64}; a::Float64 = 0.0, b::Float64 = 1.0)
    n = length(x)
    y = zeros(Float64, order, n)
    for i in 1:n
        for j in 1:order
            y[j, i] = ((x[i] - a) / (b - a))^(order - j)
        end
    end
    return y
end

function cheb_nodes_1d(order::Int, a::Float64, b::Float64)
    nodes = zeros(Float64, order)
    for i in 0:(order - 1)
        nodes[i + 1] = -cos((i + 0.5) * pi / order) * 0.5 + 0.5
        nodes[i + 1] = nodes[i + 1] * (b - a) + a
    end
    return nodes
end

function monomial_interp_1d(order::Int, nnodes::Int, fn_v::Vector{Float64}; a::Float64 = 0.0, b::Float64 = 1.0)
    x = cheb_nodes_1d(nnodes, a, b)
    p = monomial_basis_1d(order, x; a = a, b = b)
    Mp = pseudo_inv(p, eps(Float64))

    dof = length(fn_v) ÷ nnodes
    @assert length(fn_v) == dof * nnodes

    fn_mat = permutedims(reshape(fn_v, (nnodes, dof)))
    coeff = fn_mat * Mp

    return vec(permutedims(coeff))
end
