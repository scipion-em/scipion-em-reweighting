if __name__ == '__main__':
    import argparse
    from cryoER.analyze_mcmc import analyze_mcmc
    import numpy as np
    import os

    # Input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_directory', type=str, required=True)
    parser.add_argument('--filename_cluster_counts', type=str, required=True)

    args = parser.parse_args()

    if not args.output_directory.endswith('/'):
        args.output_directory += '/'

    factor_mean_std, rewtprob_mean_std, lp, log_weights_mc_chains = analyze_mcmc(
        output_directory = args.output_directory,
        filename_cluster_counts = args.filename_cluster_counts,
    )

    cluster_counts = np.loadtxt(args.filename_cluster_counts)
    nCluster = len(cluster_counts)
    Nm = cluster_counts.astype(float)
    Nm /= np.sum(Nm)

    factor_mean = factor_mean_std[:,0]
    factor_std = factor_mean_std[:,1]

    np.savetxt(os.path.join(args.output_directory, "mean_weights.txt"), factor_mean)
    np.savetxt(os.path.join(args.output_directory, "std_weights.txt"), factor_std)
