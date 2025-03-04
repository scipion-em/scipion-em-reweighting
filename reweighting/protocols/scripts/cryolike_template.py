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
    args = parser.parse_args()

    templates_dir = os.path.join(args.folder_output, "templates")
    parameters_path = os.path.join(templates_dir, "parameters.npz")
    particles_dir = os.path.join(args.folder_output, "particles")
    likelihood_dir = os.path.join(args.folder_output, "likelihood")


    #### Step 3: Make Fourier reprojection templates (no CTF) for ref volume ####
    make_templates_from_inputs(
        list_of_inputs = [args.ref],
        image_parameters_file = parameters_path,
        folder_output = templates_dir,
        verbose = True,
        use_cuda = args.use_cuda,
        i_start = args.i
    )
