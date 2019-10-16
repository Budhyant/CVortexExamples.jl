##############################################################################
#
# ParticleVortexTunnel.jl
#
# A demo of using vortex filament rings to create a tunnel.
# Includes code for particle redistribution and viscocity.
# Euler ODE integration.
# Dumps output to particle_tunnel_output_STEP.vtu
#
# Copyright HJA Bird 2018-2019
#
##############################################################################

import WriteVTK
using CVortex
using Dates

let
    # We can print information about the system.
    println("Running ParticleVortexTunnel.jl")
    println("CVortex found ", number_of_accelerators(), " accelerators "
        *"on this system.")
    println("There are ", number_of_enabled_accelerators(), 
        " enabled accelerators. These are:")
    for i = 1 : number_of_enabled_accelerators()
        println(i, ":\t", accelerator_name(i))
    end

    # Set up the tubes's geometry:
    radius = 0.5
    len = 1
    # The discretisation
    n_rings = 40
    p_per_ring = 120
    # The strength of the vortex particles and the time step parameters
    str = 1 / n_rings
    n_steps = 100
    dt = 0.025
    # And the regularisation of the particle-particle interaction
    kernel = gaussian_regularisation()
	particle_sep = max(len / n_rings, 2 * pi * radius / p_per_ring)
    regdist = particle_sep * 2
	redistribute_every = 10
	kinematic_visc = 0#1/10000
	save_string = "particle_tunnel_output_"

    total_points = n_rings * p_per_ring
    total_particles = n_rings * p_per_ring
    println("Total number of particles is ", total_particles, ".")

    # Create a load of vortex filaments in rings:
    particle_pos = zeros(total_particles, 3)
    particle_vorts = zeros(total_particles, 3)
	particle_vols = particle_sep^3 * ones(total_particles)
    acc = 1
    for ring_idx = 1 : n_rings
        for p_in_ring_idx = 1 : p_per_ring
            xloc = ring_idx * len / n_rings
            yloc = radius * cos(2 * pi * p_in_ring_idx / p_per_ring)
            zloc = radius * sin(2 * pi * p_in_ring_idx / p_per_ring)
            yvort = str * 2 * pi * radius * -sin(2 * pi * p_in_ring_idx / p_per_ring) / p_per_ring
            zvort = str * 2 * pi * radius * cos(2 * pi * p_in_ring_idx / p_per_ring) / p_per_ring
            particle_pos[acc, :] = [xloc, yloc, zloc]
            particle_vorts[acc, :] = [0.0, yvort, zvort] 
            acc += 1
        end
    end

    # A method to save the simulation to a file
    function save_particles(step)
		total_particles = size(particle_pos)[1]
        cells = Vector{WriteVTK.MeshCell}(undef, total_particles)
        cells = map(
            x->WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_VERTEX, [x]), 
            1:size(particle_pos)[1])
        vtkfile = WriteVTK.vtk_grid(save_string*string(step), 
            transpose(particle_pos), cells)
		WriteVTK.vtk_point_data(vtkfile, transpose(particle_vorts), "Vorticity")
        WriteVTK.vtk_save(vtkfile)
    end

    save_particles(0)
    print("Started.")
    for i = 1 : n_steps
		total_particles = size(particle_pos)[1]
		ninter = total_particles^2
        tstart = now()
        vels = particle_induced_velocity(particle_pos, particle_vorts, 
            particle_pos, kernel, regdist)
        dvorts = particle_induced_dvort(particle_pos, particle_vorts, 
            particle_pos, particle_vorts, kernel, regdist)
		if kinematic_visc == 0
			vdvorts = zeros(total_particles, 3)
		else
			vdvorts = particle_visc_induced_dvort(particle_pos, particle_vorts, particle_vols,
				particle_pos, particle_vorts, particle_vols, kernel, regdist, kinematic_visc)
		end
        tend = now()
        print("\rStep:\t", i, "\tInteractions per second: ", 
            1000*round(ninter / Float64((tend - tstart).value)),"\t\t\t")
        particle_pos .+= dt .* vels
        particle_vorts .+= dt * (dvorts + vdvorts)
        if i % 1 == 0
            save_particles(i)
        end
		if i % redistribute_every == 0
			particle_pos, particle_vorts, particle_vols = CVortex.redistribute_particles_on_grid(particle_pos, particle_vorts,
				CVortex.m4p_redistribution(), particle_sep; negligible_vort=0.1)
		end
    end
end
