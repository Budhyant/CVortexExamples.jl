##############################################################################
#
# BenchmarkF3D.jl
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
	println("\nThis program calculates the velocities "*
		" terms for different sized groups of filaments and"*
		" times the number of interactions per second. The number of"*
		" measurement points always matches the number of filaments.\n")
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
    nfilaments = sort(vcat(map(i->2^i, 6:30), [1023]))

    # And the regularisation of the particle-particle interaction
    kernel = winckelmans_regularisation()

    r_filaments = Vector{Int64}(undef,0)
    r_times = Vector{Float64}(undef, 0)

    starttime = now()
    println("Running...")
    i = 1
    while Float64((now() - starttime).value) / 1000 < mintime
        ninter = nfilaments[i] ^2
        fstart = rand(nfilaments[i], 3)
		fend = rand(nfilaments[i], 3)
        fvorts = rand(nfilaments[i])
		mes_points = rand(nfilaments[i], 3)
        vels, t1,~,~,~ = @timed filament_induced_velocity(fstart, fend, fvorts, mes_points)
        wallclocktime = t1
        if wallclocktime > 0
            push!(r_filaments, nfilaments[i])
            push!(r_times, t1)
        end
        i += 1
    end
    println("Filaments: \tTime(s):\tBandwidth:\n")
    map(i->println(
        r_filaments[i],"\t",
        r_times[i],"\t",
        r_filaments[i]^2/r_times[i]), 1:length(r_filaments))
    return
end
