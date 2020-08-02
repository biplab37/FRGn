using Test, JLD, FRGn

velocity = load("../example/Bauer/bauer.jld","velocity")
dielectric = load("../example/Bauer/bauer.jld","dielectric")
momentum,cutoff = rand(2)
n, m = size(velocity)

@testset "get_velocity" begin

    @test get_velocity(velocity,rand(),1.0,m,n) == 1.0
    @test get_velocity(velocity,rand(),rand(),1.0,m,n) == (1.0,1.0)

    @test get_velocity(velocity,momentum,cutoff,m,n) == velocity[Int64(round(momentum*n)),m+2- Int64(round(m + 1 - (cutoff*m)))]
end

@testset "get_dielectric" begin

    @test get_dielectric(dielectric,rand(),1.0,m,n) == 1.0
    @test get_dielectric(dielectric,rand(),rand(),1.0,m,n) == (1.0,1.0)

    @test get_dielectric(dielectric,momentum,cutoff,m,n) == dielectric[Int64(round(momentum*n)),m+2- Int64(round(m + 1 - (cutoff*m)))]
end
