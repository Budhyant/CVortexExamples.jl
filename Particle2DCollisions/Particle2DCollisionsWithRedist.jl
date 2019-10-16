##############################################################################
#
# Particle2DCollisionsWithRedist.jl
#
# Models sets of particles colliding in 2D with particle redistribution.
# Outputs to particles2D_WR_STEP.vtu.
#
# Copyright HJA Bird 2019
#
##############################################################################

import WriteVTK
using CVortex
using Dates

let
    # We can print information about the system.
    println("Running Particle2DCollisionsWithRedist.jl")
    println("CVortex found ", number_of_accelerators(), " accelerators "
        *"on this system.")
    println("There are ", number_of_enabled_accelerators(), 
        " enabled accelerators. These are:")
    for i = 1 : number_of_enabled_accelerators()
        println(i, ":\t", accelerator_name(i))
    end

    n_steps = 10000
    dt = 0.01
    regdist = 0.05
    kernel = winckelmans_regularisation()
    vorticity = 0.001
    period = 1
    offset = 0.5

    xs = collect(1:0.0025:5)
    ys = collect(-0.000 : 0.0025 : 0.0025)
    
    particle_pos = zeros(length(xs) * length(ys)*4, 2)
    particle_vorts = ones(size(particle_pos)[1]) * vorticity
    acc = 1
    for j = 1:length(ys)
        for i = 1:length(xs)
            particle_pos[acc,:] = [xs[i], ys[j]] # Makes save to VTK easier
            particle_vorts[acc] *= -sin(2 * pi * period * xs[i])^4
            acc += 1
        end
    end
    for j = 1:length(ys)
        for i = 1:length(xs)
            particle_pos[acc,:] = [xs[i], ys[j]+offset] # Makes save to VTK easier
            particle_vorts[acc] *= sin(2 * pi * period * xs[i])^4
            acc += 1
        end
    end
    for j = 1:length(ys)
        for i = 1:length(xs)
            particle_pos[acc,:] = [-xs[i], ys[j]] # Makes save to VTK easier
            particle_vorts[acc] *= sin(2 * pi * period * xs[i])^4
            acc += 1
        end
    end
    for j = 1:length(ys)
        for i = 1:length(xs)
            particle_pos[acc,:] = [-xs[i], ys[j]+offset] # Makes save to VTK easier
            particle_vorts[acc] *= -sin(2 * pi * period * xs[i])^4
            acc += 1
        end
    end

    # A method to save the simulation to a file
    function save_particles(step)
		total_particles = size(particle_pos)[1]
		pos3d = hcat(particle_pos, zeros(total_particles, 1))
        cells = Vector{WriteVTK.MeshCell}(undef, total_particles)
        cells = map(
            x->WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_VERTEX, [x]), 
            1:size(pos3d)[1])
        vtkfile = WriteVTK.vtk_grid("particles2D_WR_"*string(step), 
            pos3d', cells)
        WriteVTK.vtk_point_data(vtkfile, particle_vorts', "Vorticity")
        WriteVTK.vtk_save(vtkfile)
    end

    save_particles(0)
    
    print("Started.")
	
    for i = 1 : n_steps
        tstart = now()
		nparticles = length(particle_vorts)
		ninter = nparticles^2
        vels = particle_induced_velocity(particle_pos, particle_vorts, 
            particle_pos, kernel, regdist)
        tend = now()
        print("\rStep:\t", i, "\tInteractions per second: ", 
            1000*round(ninter / Float64((tend - tstart).value)),"\t\t\t")
        particle_pos .+= dt .* vels
		if i % 10 == 0
			particle_pos, particle_vorts, ~ = redistribute_particles_on_grid(
				particle_pos, particle_vorts, lambda3_redistribution(), regdist*0.67)
		end
        if i % 50 == 0
            save_particles(i)
        end
    end
    
end
