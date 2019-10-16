##############################################################################
#
# Tracer.jl
# Entirely for visual effect. Three sets of vortex particles mix in 3D.
# Because so many particle are used, this allows the shear layers to be
# tracked. Consequentially, a pretty picture of the mixing can be created.
#
# Copyright HJA Bird 2018-2019
#
##############################################################################

import WriteVTK
using CVortex
using Dates

let
    # We can print information about the system.
    println("Running tracer.jl")
    println("CVortex found ", number_of_accelerators(), " accelerators "
        *"on this system.")
    println("There are ", number_of_enabled_accelerators(), 
        " enabled accelerators. These are:")
    for i = 1 : number_of_enabled_accelerators()
        println(i, ":\t", accelerator_name(i))
    end

    n_steps = 3100
    dt = 0.01
    regdist = 0.05
    kernel = winckelmans_regularisation()
    vorticity = 1.
    period = 1
    offset = 0.5

	pstep = 0.0008
    xsr = collect(0:pstep:0.025)
    xsg = collect(0.05:pstep:0.075)
    xsb = collect(0.1:pstep:0.125)
    ys = collect(0:pstep:1)
	nr = length(xsr) * length(ys)
	ng = length(xsg) * length(ys)
	nb = length(xsb) * length(ys)
	ntot = nr + ng + nb
    
    particle_pos = zeros(ntot, 3)
    particle_vorts = ones(ntot) * vorticity / ntot
    acc = 1
    for j = 1:length(ys)
        for i = 1:length(xsr)
            particle_pos[acc,:] = [xsr[i], ys[j], 0] # Makes save to VTK easier
			particle_vorts[acc] *= sin(ys[j] * pi)
            acc += 1
        end
    end
    for j = 1:length(ys)
        for i = 1:length(xsg)
            particle_pos[acc,:] = [xsg[i], ys[j], 0] # Makes save to VTK easier
			particle_vorts[acc] *= sin(ys[j] * pi)
            acc += 1
        end
    end
    for j = 1:length(ys)
        for i = 1:length(xsb)
            particle_pos[acc,:] = [xsb[i], ys[j], 0] # Makes save to VTK easier
			particle_vorts[acc] *= sin(ys[j] * pi)
            acc += 1
        end
    end

    println("Total of ", ntot, " particles.")

    # A method to save the simulation to a file
    function save_particles(step)
        cellsr = map(
            x->WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_VERTEX, [x]), 
            1:nr)
        vtkfile = WriteVTK.vtk_grid("tracer_r_"*string(step), 
            particle_pos[1:nr,:]', cellsr)
        WriteVTK.vtk_save(vtkfile)
		
        cellsg = map(
            x->WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_VERTEX, [x]), 
            1:ng)
        vtkfile = WriteVTK.vtk_grid("tracer_g_"*string(step), 
            particle_pos[nr+1:nr+ng,:]', cellsg)
        WriteVTK.vtk_save(vtkfile)
		
        cellsb = map(
            x->WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_VERTEX, [x]), 
            1:nb)
        vtkfile = WriteVTK.vtk_grid("tracer_b_"*string(step), 
            particle_pos[nr+ng+1:end,:]', cellsb)
        WriteVTK.vtk_save(vtkfile)
    end

    save_particles(0)
    
    print("Started for ", n_steps, " steps.\n")
	ninter = ntot^2
    for i = 1 : n_steps
		tstart = now()
        vels = particle_induced_velocity(particle_pos[:, 1:2], particle_vorts, 
            particle_pos[:, 1:2], kernel, regdist)
        particle_pos[:, 1:2] .+= dt .* vels        
        tend = now()
        print("\rStep:\t", i, "\tInteractions per second: ", 
            1000*round(ninter / Float64((tend - tstart).value)),"\t\t\t")
        if i % 10 == 0
            save_particles(i)
        end
    end
    
end
