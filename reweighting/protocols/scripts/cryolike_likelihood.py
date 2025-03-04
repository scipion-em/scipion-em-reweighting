if __name__ == '__main__':
    import argparse
    import numpy as np
    import os
    import torch

    from cryolike.stacks.make_templates_from_inputs_api import make_templates_from_inputs
    from cryolike.run_likelihood import run_likelihood

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', type=int, required=True)
    parser.add_argument('--ref', type=str, required=True)
    parser.add_argument('--folder_output', type=str, required=True)
    parser.add_argument('--use_cuda', required=False, 
                        default=False, action='store_true')
    parser.add_argument('--batch_size', type=int, required=True)

    parser.add_argument('--max_displacement_pixels', type=int, required=True)
    parser.add_argument('--n_displacements_x', type=int, required=True)
    parser.add_argument('--n_displacements_y', type=int, required=True)

    args = parser.parse_args()

    templates_dir = os.path.join(args.folder_output, "templates")
    parameters_path = os.path.join(templates_dir, "parameters.npz")
    particles_dir = os.path.join(args.folder_output, "particles")
    likelihood_dir = os.path.join(args.folder_output, "likelihood")

    #### Step 4: Make Fourier reprojection templates (no CTF) ####
    run_likelihood(
        params_input = parameters_path,
        folder_templates = templates_dir,
        folder_particles = particles_dir,
        folder_output = likelihood_dir,
        i_template = args.i,
        n_stacks = 1,
        skip_exist = False,
        n_templates_per_batch = 16,
        n_images_per_batch = args.batch_size,
        search_batch_size = True,
        max_displacement_pixels = args.max_displacement_pixels,
        n_displacements_x = args.n_displacements_x,
        n_displacements_y = args.n_displacements_y,
        return_likelihood_integrated_pose_fourier = True,
        return_likelihood_optimal_pose_physical = False,
        return_likelihood_optimal_pose_fourier = False,
        verbose = True
    )

    ll = torch.load(os.path.join(likelihood_dir,
                                 'template%d/log_likelihood/log_likelihood_integrated_fourier_stack_000000.pt') % args.i)
    np.save(os.path.join(likelihood_dir,
                         'template%d/log_likelihood/log_likelihood_integrated_fourier_stack_000000.npy') % args.i,
                         ll)
