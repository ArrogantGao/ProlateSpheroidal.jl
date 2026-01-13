using ProlateSpheroidal
using Test

@testset "ProlateSpheroidal.jl" begin
    include(joinpath(@__DIR__, "fixtures", "pswf_reference.jl"))

    @testset "prolc180" begin
        for item in PSWF_REFERENCE[:prolc180]
            @test isapprox(ProlateSpheroidal.prolc180(item.eps), item.value; rtol = 1e-12, atol = 1e-14)
        end
    end

    @testset "prolc180_der3" begin
        for item in PSWF_REFERENCE[:prolc180_der3]
            @test isapprox(ProlateSpheroidal.prolc180_der3(item.eps), item.value; rtol = 1e-12, atol = 1e-14)
        end
    end

    @testset "prolate0_eval" begin
        params_map = Dict{Float64, ProlateSpheroidal.Prolate0Params}()
        for item in PSWF_REFERENCE[:prolate0_eval]
            params = get!(params_map, item.c) do
                ProlateSpheroidal.Prolate0Params(item.c)
            end
            @test isapprox(ProlateSpheroidal.prolate0_eval(params, item.x), item.value; rtol = 1e-9, atol = 1e-11)
        end
    end

    @testset "prolate0_eval_derivative" begin
        params_map = Dict{Float64, ProlateSpheroidal.Prolate0Params}()
        for item in PSWF_REFERENCE[:prolate0_eval_derivative]
            params = get!(params_map, item.c) do
                ProlateSpheroidal.Prolate0Params(item.c)
            end
            @test isapprox(ProlateSpheroidal.prolate0_eval_derivative(params, item.x), item.value;
                           rtol = 1e-8, atol = 1e-10)
        end
    end

    @testset "prolate0_int_eval" begin
        params_map = Dict{Float64, ProlateSpheroidal.Prolate0Params}()
        for item in PSWF_REFERENCE[:prolate0_int_eval]
            params = get!(params_map, item.c) do
                ProlateSpheroidal.Prolate0Params(item.c)
            end
            @test isapprox(ProlateSpheroidal.prolate0_int_eval(params, item.r), item.value;
                           rtol = 1e-8, atol = 1e-10)
        end
    end

    @testset "prolate0_lambda" begin
        for item in PSWF_REFERENCE[:prolate0_lambda]
            params = ProlateSpheroidal.Prolate0Params(item.c)
            @test isapprox(ProlateSpheroidal.prolate0_lambda(params), item.value; rtol = 1e-8, atol = 1e-10)
        end
    end

    @testset "prolate0_eval_from_eps" begin
        params_map = Dict{Float64, ProlateSpheroidal.Prolate0Params}()
        for item in PSWF_REFERENCE[:prolate0_eval_from_eps]
            c = ProlateSpheroidal.prolc180(item.eps)
            params = get!(params_map, c) do
                ProlateSpheroidal.Prolate0Params(c)
            end
            @test isapprox(ProlateSpheroidal.prolate0_eval(params, item.x), item.value; rtol = 1e-9, atol = 1e-11)
        end
    end
end
