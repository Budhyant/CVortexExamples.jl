##############################################################################
#
# ParticleRings.jl
#
# A demo of using vortex particle ring interaction based on the work 
# Winckelmans and Leonard, J_Comp_Phys, 1993, Contributions to the 
# Vortex Particle Methods for the Computation of Three-Dimensional 
# Incompressible Unsteady Flows. This should produce something like 
# figure 5.
# Saves to particle_vortex_ring_STEP.vtu.
#
# Copyright HJA Bird 2018-2019
#
##############################################################################

import WriteVTK
using CVortex
using Dates
using DifferentialEquations

let
    # We can print information about the system.
    println("Running particle_vortex_tunnel.jl")
    println("CVortex found ", number_of_accelerators(), " accelerators "
        *"on this system.")
    println("There are ", number_of_enabled_accelerators(), 
        " enabled accelerators. These are:")
    for i = 1 : number_of_enabled_accelerators()
        println(i, ":\t", accelerator_name(i))
    end

    n_steps = 1810
    dt = 0.0001
    regdist = 0.001
    kernel = gaussian_regularisation()
	kinematic_visc = 1e-2

    function make_ring(r_l::Real, radius::Real, vorticity::Real,
        circ_particles::Int, layers::Int; str=ones(10))
        local theta_r = map(i->2 * pi * (i-.5) / circ_particles, 1:circ_particles)
        local particles = zeros(length(theta_r), 3)
        local vorts = zeros(length(theta_r), 3)
        local section_radius = (2 * layers - 1) * r_l
        local cell_area = pi * r_l^2
        local cross_section_area = pi * section_radius^2

        local areas = vcat([cell_area], map(i->pi*r_l^2*((2*i + 1)^2-(2*i-1)^2), 1:layers-1))#
        local vortadj = vorticity / sum(areas .* str[1:length(areas)])
        local cell_vpln = vortadj * cell_area * 2 * pi * radius / circ_particles

        function sub_ring(theta_s, rc, str)
            local p=mapreduce(
                theta->[(radius+rc*cos(theta_s))*cos(theta), (radius+rc*cos(theta_s))*sin(theta), rc*sin(theta_s)]', 
                vcat, theta_r)
            local v=mapreduce(theta->str * [-sin(theta), cos(theta), 0]', 
                vcat, theta_r)
            return p, v
        end

        # Ring centre.
        local particles, vorts = sub_ring(0., 0., cell_vpln * str[1])
        # Outer layers
        for i = 1 : layers-1
            r_c = r_l * (1 + 12*i^2)/(6 * i)
            nr = 8 * i
            theta_s = map(j->2 * pi * j / nr, 1:nr)
            for theta in theta_s
                ps, vs = sub_ring(theta, r_c, cell_vpln * str[i + 1])
                particles=vcat(particles, ps)
                vorts=vcat(vorts, vs)
            end
        end
		volumes = ones(size(particles)[1]) .* cell_area * 2 * pi * radius / circ_particles
        return particles, vorts, volumes
    end
    function translate(points, dx)
        for i in 1 : size(points)[1]
            points[i, :] += dx
        end
        return points
    end
    function rotate(points, vorts, theta_x::Real, centre::Vector{<:Real})
        points = translate(points, -1 .*centre)
        mat = [1 0 0; 0 cos(theta_x) -sin(theta_x); 0 sin(theta_x) cos(theta_x)]
        for i in 1 : size(points)[1]
            points[i, :] = mat * points[i, :]
            vorts[i, :] = mat * vorts[i, :]
        end
        points = translate(points, centre)
        return points, vorts
    end

    particle_pos1, vorticities1, volumes1 = make_ring(0.00002, 0.0267, -3.1e-2, 800, 2)
    particle_pos2 = deepcopy(particle_pos1)
    vorticities2 = deepcopy(vorticities1)
    particle_pos1, vorticities1 = rotate(particle_pos1, vorticities1, deg2rad(20), [0., 0., 0.])
    particle_pos2 = translate(particle_pos2, [0.0267, 0., 0.])

    particle_pos = vcat(particle_pos1, particle_pos2)
    particle_vorts = vcat(vorticities1, vorticities2)
    particle_vols = vcat(volumes1, volumes1)

    total_particles = size(particle_pos)[1]
    println("Total of ", total_particles, " particles.")

    # A method to save the simulation to a file
    function save_particles(step)
        cells = Vector{WriteVTK.MeshCell}(undef, total_particles)
        cells = map(
            x->WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_VERTEX, [x]), 
            1:size(particle_pos)[1])
        vtkfile = WriteVTK.vtk_grid("particle_vortex_ring_"*string(step), 
            particle_pos', cells)
        WriteVTK.vtk_point_data(vtkfile, particle_vorts', "Vorticity")
        WriteVTK.vtk_save(vtkfile)
    end

    save_particles(0)
    
    print("Started.")
    ninter = total_particles^2
    for i = 1 : n_steps
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
        if i % 10 == 0
            save_particles(i)
        end
    end
    
end
