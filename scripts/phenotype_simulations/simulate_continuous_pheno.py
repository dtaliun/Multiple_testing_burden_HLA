import argparse
import gzip
import numpy as np
import statsmodels.api as sm
from scipy.special import expit
from scipy.optimize import brentq
import sys


COVAR_REQUIRED_HEADER = ['IID', 'Sex', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']


argparser = argparse.ArgumentParser(description = 'This script simulates continuous  phenotypes for samples from the input VCF file.')
argparser.add_argument('-i', '--input-vcf', metavar = 'file', dest = 'in_vcf', type = str, required = True, help = 'Input VCF.')
argparser.add_argument('-c', '--covariates', metavar = 'file', dest = 'in_covar_tsv', type = str, required = False, help = f'Input covariate file as PLINK2-compatible TSV. Required columns in order: {COVAR_REQUIRED_HEADER}')
argparser.add_argument('-n', '--n-phenotypes', metavar = 'integer', dest = 'n_pheno', type = int, default = 1, help = 'Number of phenotypes to simulate (default: %(default)s).')
argparser.add_argument('-s', '--seed', metavar = 'integer', dest = 'seed', type = int, default = 1, help = 'Seed for simulations (default: %(default)s).')
argparser.add_argument('-d', '--diagnostic-output', metavar = 'file', dest = 'diag_tsv', type = str, required = True, help = 'Output diagnostic TSV file.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_tsv', type = str, required = True, help = 'Output TSV file with simulated phenotypes.')


if __name__ == '__main__':
    args = argparser.parse_args()

    header = None
    samples = []
    n_samples = 0

    # Get sample names from VCF.
    with gzip.open(args.in_vcf, 'rt') as ifile:
        for line in ifile:
            if line.startswith('#CHROM'):
                header = line.rstrip().split('\t')
                # At least one sample
                assert len(header) > 9
                samples = header[9:]
                break

    assert header is not None

    n_samples = len(samples)
    assert n_samples > 0

    assert args.n_pheno > 0


    # If covariate file is specified, then attempt to model some non-zero covariate effects
    if args.in_covar_tsv is not None:
        covar_sample_index = {}
        with open(args.in_covar_tsv, 'rt') as ifile:
            header = ifile.readline().rstrip().split('\t')
            
            # Header has required column names in required order
            assert len(header) == len(COVAR_REQUIRED_HEADER)
            assert all(x.lower() == y.lower() for x, y in zip(header, COVAR_REQUIRED_HEADER))
            
            n_covars = len(header) - 1
            n_lines = 0
            
            covars = np.zeros((n_samples, n_covars))
            
            for i, line in enumerate(ifile):
                # Number of samples must not exceed that in VCF
                assert i < n_samples

                n_lines = i + 1
                fields = line.rstrip().split('\t')
                
                # No truncated lines
                assert len(fields) == len(COVAR_REQUIRED_HEADER)
                
                # Same sample as in VCF
                assert fields[0] in samples, fields[0]
               
                # No duplicated IIDs
                assert fields[0] not in covar_sample_index, fields[0]
                
                covar_sample_index[fields[0]] = i
 
                for j in range(0, n_covars):
                    covars[i, j] = float(fields[j + 1])

            # Number of samples must match VCF
            assert n_lines == n_samples
            assert len(covar_sample_index) == n_samples
            
            beta = np.zeros(n_covars)

            # Sex effect
            beta[0] = 0.1

            # Standardize PCs
            pcs = covars[:, 1:].copy()
            pc_mean = pcs.mean(axis=0)
            pc_std = pcs.std(axis=0)

            # PCs can't have standard deviation 0
            assert np.all(pc_std > 1e-12)

            covars[:, 1:] = (pcs - pc_mean) / pc_std

            # Generate a sequence of effects for PCs starting with 0.3 and geometrically dicreasing towards 0
            # Then pass them as effects
            for i, pc_effect in enumerate(np.geomspace(0.3, 0.01, num = n_covars - 1)):
                beta[i + 1] = pc_effect
           
           
            # Compute final covariate effects across individuals
            mu = covars @ beta

            # Diagnostics
            sex = covars[:, 0]

    else:
        mu = np.zeros(n_samples)

    diag_file = None
    if args.diag_tsv is not None:
        diag_file = open(args.diag_tsv, 'wt')

        if args.in_covar_tsv is not None:
            diag_file.write(
                "PHENO_ID\tMEAN_Y\tVAR_Y\tSEX_BETA\tSEX_CORR\tR2\tCORR_Y_MU\n"
            )
        else:
            diag_file.write(
                "PHENO_ID\tMEAN_Y\tVAR_Y\n"
            )


    # Simulate Bernoulli trials
    rng = np.random.default_rng(seed = args.seed)

    phenos = []
    for i in range(0, args.n_pheno):
        sigma = 1.0 # environmental SD

        y = rng.normal(
            loc = mu,
            scale = sigma,
            size = n_samples
        )

        phenos.append(y)
        
        if diag_file is not None:
            mean_y = np.mean(y)
            var_y = np.var(y)

            if args.in_covar_tsv is None:
                diag_file.write(
                    f"{i + 1}\t"
                    f"{mean_y:.6f}\t"
                    f"{var_y:.6f}\n"
                )
            else:
                fit = sm.OLS(y, covars).fit()

                mean_y = np.mean(y)
                var_y = np.var(y)

                sex_beta = fit.params[0]
                sex_corr = np.corrcoef(y, sex)[0, 1]

                corr_y_mu = np.corrcoef(y, mu)[0, 1]
                r2 = fit.rsquared
                
                diag_file.write(
                    f"{i+1}\t"
                    f"{mean_y:.6f}\t"
                    f"{var_y:.6f}\t"
                    f"{sex_beta:.6f}\t"
                    f"{sex_corr:.6f}\t"
                    f"{r2:.6f}\t"
                    f"{corr_y_mu:.6f}\n"
                )
               
    
    if diag_file is not None:
        diag_file.close()

    # Write results
    with open(args.out_tsv, 'wt') as ofile:
        ofile.write("IID")
        for i in range(1, len(phenos) + 1):
            ofile.write(f"\tPHENO{i}")
        ofile.write("\n")
    
        if args.in_covar_tsv is not None:
            for sample in samples:
                ofile.write(f"{sample}")
                for pheno in phenos:
                    ofile.write(f"\t{pheno[covar_sample_index[sample]]:.6f}")
                ofile.write("\n")
        else:
            for i, sample in enumerate(samples):
                ofile.write(f"{sample}")
                for pheno in phenos:
                    ofile.write(f"\t{pheno[i]:.6f}")
                ofile.write("\n")
            

    
