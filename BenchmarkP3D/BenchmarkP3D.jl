##############################################################################
#
# BenchmarkP3D.jl
#
# Benchmark the number of interactions a system can model per second
# with and without acceleration.
#
# Copyright HJA Bird 2019
#
##############################################################################

using CVortex
using Dates

let
    # We can print information about the system.
    println("BenchmarkP3D.jl")
	println("\nThis program calculates the velocities and vortex stretching"*
		" terms for different sized groups of vortex particles and times"*
		" times the number of interactions per second. The number of"*
		" measurement points always matches the number of particles.\n")
    println("CVortex found ", number_of_accelerators(), " accelerators "
        *"on this system.")
    println("There are ", number_of_enabled_accelerators(), 
        " enabled accelerators. These are:")
    for i = 1 : number_of_enabled_accelerators()
        println(i, ":\t", accelerator_name(i))
    end

    # We'll only start new runs for a certain amount of time:
    mintime = 30
    println("Mintime is ", mintime, " seconds. The program will finish
        what its working on when it runs out of time.")
    nparticles = sort(vcat(map(i->2^i, 6:30), [1023]))

    # And the regularisation of the particle-particle interaction
    kernel = winckelmans_regularisation()

    r_particles = Vector{Int64}(undef,0)
    r_times = Vector{Float64}(undef, 0)

    starttime = now()
    println("Running...")
    i = 1
    while Float64((now() - starttime).value) / 1000 < mintime
        ninter = nparticles[i] ^2
        particle_pos = rand(nparticles[i], 3)
        particle_vorts = rand(nparticles[i], 3)
        vels, t1,~,~,~ = @timed particle_induced_velocity(particle_pos, particle_vorts, 
            particle_pos, kernel, 0.01)
        dvorts, t2,~,~,~ = @timed particle_induced_dvort(particle_pos, particle_vorts, 
            particle_pos, particle_vorts, kernel, 0.01)
        wallclocktime = t1 + t2
        if wallclocktime > 0
            push!(r_particles, nparticles[i])
            push!(r_times, t1 + t2)
        end
        i += 1
    end
    println("Particles: \tTime(s):\tBandwidth:\n")
    map(i->println(
        r_particles[i],"\t",
        r_times[i],"\t",
        r_particles[i]^2/r_times[i]), 1:length(r_particles))
    return
end
