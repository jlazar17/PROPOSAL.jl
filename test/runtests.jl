using Test
using PROPOSAL

@testset "PROPOSAL.jl" begin
    # Basic tests will be added once wrapper is functional

    @testset "Module loads" begin
        @test isdefined(PROPOSAL, :PROPOSAL)
    end

    # @testset "Particle definitions" begin
    #     @test MuMinusDef() isa PROPOSAL.ParticleDef
    #     @test TauMinusDef() isa PROPOSAL.ParticleDef
    # end

    # @testset "Basic propagation" begin
    #     # Create a simple propagator and verify it works
    #     # prop = Propagator(MuMinusDef(), config_path)
    #     # state = ParticleState(...)
    #     # track = propagate(prop, state)
    #     # @test length(track) > 0
    # end
end
