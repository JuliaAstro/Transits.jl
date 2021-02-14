using Transits
using Limbdark
using LaTeXStrings
using BenchmarkTools
using DataFrames
using Statistics
using StatsPlots, ColorSchemes
using ProgressLogging


N_pts = [100, 316, 1000, 3160, 10000, 31600, 100000]
N_us = [1 2 3 5 8 13 21 34 55 90 144]

function benchmark(N_u, N_pts; r=0.1)
    u = ones(N_u) ./ N_u
    b = [sqrt(((i - N_pts/2) * 2/N_pts * (1 + 2 * r))^2) for i in 1:N_pts]

    ld = PolynomialLimbDark(u)
    output = similar(b)
    bench = @benchmark begin
        @inbounds for i in eachindex($b)
            $output[i] = compute($ld, $b[i], $r)
        end
    end
    return median(bench).time * 1e-9
end
vals = zeros(length(N_pts), length(N_us))
@progress for (j, Nu) in enumerate(N_us)
    @progress for (i, N) in enumerate(N_pts)
        vals[i, j] = benchmark(Nu, N)
    end
end

labels = reduce(hcat, "uₙ: $N" for N in N_us)

plot(
    N_pts,
    vals,
    scale=:log10,
    ls=:solid,
    m=:o,
    palette=palette(:plasma, length(N_us) + 2),
    msw=0,
    lab=labels,
    xlabel="number of points",
    ylabel="timing [s]",
    leg=:topleft
)


function benchmark_limbdark(N_u, N_pts; r=0.1)
    u = ones(N_u) ./ N_u
    b = [sqrt(((i - N_pts/2) * 2/N_pts * (1 + 2 * r))^2) for i in 1:N_pts]

    trans = Limbdark.transit_init(r, b[1], u, true)
    output = similar(b)
    bench = @benchmark begin
        @inbounds for i in eachindex($b)
            $trans.b = $b[i]
            $output[i] = Limbdark.transit_poly!($trans)
        end
    end
    return median(bench).time * 1e-9
end
vals_limbdark = zeros(length(N_pts), length(N_us))
@progress for (j, Nu) in enumerate(N_us)
    @progress for (i, N) in enumerate(N_pts)
        vals_limbdark[i, j] = benchmark_limbdark(Nu, N)
    end
end

plot(
    N_pts,
    vals_limbdark ./ vals,
    xscale=:log10,
    ls=:solid,
    m=:o,
    msw=0,
    palette=palette(:plasma, length(N_us) + 2),
    lab=labels,
    xlabel="number of points",
    ylabel="relative speedup (times faster) - Limbdark.jl"
)

tmed = median(vals ./ vals[:, 1], dims=1)
plot(
    N_us |> vec,
    tmed |> vec,
    scale=:log10,
    ls=:solid,
    m=:o,
    msw=0,
    lab="Transits.jl",
    xlabel="number of limb-darkening coefficients",
    ylabel="relative timing",
    ylims=(0.8, 10),
    leg=:topleft
)
plot!(
    N_us |> vec, 
    median(vals_limbdark ./ vals_limbdark[:, 1], dims=1) |> vec,
    lab="Limbdark.jl",
    c=2,
    ls=:solid,
    m=:o,
    msw=0
)
plot!(N_us |> vec, n -> n^(0.2), ls=:dash, c=:gray, lab=L"N^{0.2}")
plot!(N_us |> vec, n -> tmed[end] * n / N_us[1, end], ls=:dash, c=:gray, lab=L"N^1")

