if __name__ == '__main__':
    import argparse
    from numpy import pi
    import os

    from cryolike.metadata import ImageDescriptor
    from cryolike.convert_particle_stacks import convert_particle_stacks_from_star_files

    def list_str(values):
        return values.split(',')

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_particles_star_file', type=str, required=True)
    parser.add_argument('--input_particles_stack_files', type=list_str, required=True)
    parser.add_argument('--folder_output', type=str, required=True)

    parser.add_argument('--pixel_size', type=float, required=True)
    parser.add_argument('--box_size', type=int, required=True)
    parser.add_argument('--batch_size', type=int, required=True)

    parser.add_argument('--viewing_distance', type=float, required=True)
    parser.add_argument('--n_inplanes', type=int, required=True)

    parser.add_argument('--use_cuda', required=False, 
                        default=False, action='store_true')

    args = parser.parse_args()

    templates_dir = os.path.join(args.folder_output, "templates")
    parameters_path = os.path.join(templates_dir, "parameters.npz")
    particles_dir = os.path.join(args.folder_output, "particles")
    likelihood_dir = os.path.join(args.folder_output, "likelihood")

    #### Step 1: Write image parameters ####
    os.makedirs(templates_dir, exist_ok=True)
    image_parameters = ImageDescriptor.from_individual_values(
        n_pixels = args.box_size,
        pixel_size = args.pixel_size,
        resolution_factor = 1.0,
        precision = 'single',
        viewing_distance = args.viewing_distance / (4.0 * pi),
        n_inplanes = args.n_inplanes,
        use_protein_residue_model = True,
        atom_shape = 'gaussian'
    )
    image_parameters.save(parameters_path)

    #### Step 2: Convert particles to Fourier ####
    convert_particle_stacks_from_star_files(
        params_input = parameters_path,
        folder_output = particles_dir,
        particle_file_list = args.input_particles_stack_files,
        star_file_list = [args.input_particles_star_file],
        pixel_size = args.pixel_size,
        defocus_angle_is_degree = True,
        phase_shift_is_degree = True,
        skip_exist = False,
        flag_plots = True,
        batch_size = args.batch_size,
        use_cuda = args.use_cuda
    )
