import argparse
import gzip
import numpy as np
import statsmodels.api as sm
from scipy.special import expit
from scipy.optimize import brentq
from scipy.stats import bernoulli
import sys


COVAR_REQUIRED_HEADER = ['IID', 'Sex', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']


argparser = argparse.ArgumentParser(description = 'This script simulates binary phenotypes for samples from the input VCF file.')
argparser.add_argument('-i', '--input-vcf', metavar = 'file', dest = 'in_vcf', type = str, required = True, help = 'Input VCF.')
argparser.add_argument('-c', '--covariates', metavar = 'file', dest = 'in_covar_tsv', type = str, required = False, help = f'Input covariate file as PLINK2-compatible TSV. Required columns in order: {COVAR_REQUIRED_HEADER}')
argparser.add_argument('-n', '--n-phenotypes', metavar = 'integer', dest = 'n_pheno', type = int, default = 1, help = 'Number of phenotypes to simulate (default: %(default)s).')
argparser.add_argument('-s', '--seed', metavar = 'integer', dest = 'seed', type = int, default = 1, help = 'Seed for simulations (default: %(default)s).')
argparser.add_argument('-f', '--case-fraction', metavar = 'float', dest = 'case_fraction', type = float, default = 0.5, help = 'Case fraction (default: %(default)s).') 
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

    assert args.case_fraction > 0.0 and args.case_fraction < 1.0

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
            
            covars = sm.add_constant(covars)
            beta = np.zeros(n_covars + 1)

            # Sex effect
            beta[1] = np.log(1.5) # odds-ratio 1.5

            # Standardize PCs
            pcs = covars[:, 2:].copy()
            pc_mean = pcs.mean(axis = 0)
            pc_std = pcs.std(axis = 0)

            # PCs can't have standard deviation 0
            assert np.all(pc_std > 1e-12)

            covars[:, 2:] = (pcs - pc_mean) / pc_std

            # Generate a sequence of odds-ratios for PCs starting with 1.4 and geometrically dicreasing towards 1
            # Then pass them as effects
            for i, pc_or in enumerate(np.geomspace(1.4, 1.01, num = n_covars - 1)):
                beta[i + 2] = np.log(pc_or)
           
            # Estimating intercept effect so that it corresponds to the required case control balance          
            eta = covars[:, 1:] @ beta[1:]
            def f(intercept):
                return np.mean(expit(intercept + eta)) - args.case_fraction
            beta[0] = brentq(f, -20, 20)
            
            # Compute final covariate effects across individuals
            z = covars @ beta

            # Transform to probabilities
            p = expit(z)

            # Diagnostics
            sex = covars[:, 1]

    else:
        p = [args.case_fraction] * n_samples

    print(f"Target prevalence   : {args.case_fraction:.6f}")
    print(f"Expected prevalence : {np.mean(p):.6f}")

    diag_file = None
    if args.diag_tsv is not None:
        diag_file = open(args.diag_tsv, 'wt')

        if args.in_covar_tsv is not None:
            diag_file.write(
                "PHENO_ID\t"
                "N_CASES\t"
                "CASE_FRAC\t"
                "SEX_OR\t"
                "SEX1_CASE_FRAC\t"
                "SEX2_CASE_FRAC\t"
                "CORR_Y_P\t"
                "LOGLIK\t"
                "NULL_CONVERGED\t"
                "NULL_LLF\t"
                "NULL_PSEUDO_R2\t"
                "NULL_SEX_OR\n"
            )
        else:
            diag_file.write(
                "PHENO_ID\tN_CASES\tCASE_FRAC\n"
            )


    # Simulate Bernoulli trials
    rng = np.random.default_rng(seed = args.seed)

    phenos = []
    for i in range(0, args.n_pheno):
        y = bernoulli.rvs(
            p = p,
            size = n_samples,
            random_state = rng
        )

        phenos.append(y)
        
        if diag_file is not None:
            n_cases = int(np.sum(y))
            case_frac = np.mean(y)

            if args.in_covar_tsv is None:
                diag_file.write(
                    f"{i + 1}\t"
                    f"{n_cases}\t"
                    f"{case_frac:.6f}\n"
                )
            else:
                sex1 = (sex == 1)
                sex2 = (sex == 2)

                a = np.sum(sex1 & (y == 1))
                b = np.sum(sex1 & (y == 0))
                c = np.sum(sex2 & (y == 1))
                d = np.sum(sex2 & (y == 0))

                sex_or = (
                    ((c + 0.5) * (b + 0.5)) /
                    ((a + 0.5) * (d + 0.5))
                )

                sex1_case_frac = a / (a + b)
                sex2_case_frac = c / (c + d)

                corr_y_p = np.corrcoef(y, p)[0, 1]

                loglik = np.sum(
                    y * np.log(p) +
                    (1 - y) * np.log(1 - p)
                )

                null_converged = 0
                null_llf = np.nan
                null_pseudo_r2 = np.nan
                null_sex_or = np.nan

                try:
                    null_fit = sm.Logit(y, covars).fit(disp=False)
                    null_converged = int(null_fit.mle_retvals['converged'])
                    null_llf = null_fit.llf
                    null_pseudo_r2 = null_fit.prsquared

                    # Sex coefficient is beta[1]
                    null_sex_or = np.exp(null_fit.params[1])
                except Exception:
                    pass

                diag_file.write(
                    f"{i + 1}\t"
                    f"{n_cases}\t"
                    f"{case_frac:.6f}\t"
                    f"{sex_or:.6f}\t"
                    f"{sex1_case_frac:.6f}\t"
                    f"{sex2_case_frac:.6f}\t"
                    f"{corr_y_p:.6f}\t"
                    f"{loglik:.6f}\t"
                    f"{null_converged}\t"
                    f"{null_llf:.6f}\t"
                    f"{null_pseudo_r2:.6f}\t"
                    f"{null_sex_or:.6f}\n"
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
                    ofile.write(f"\t{pheno[covar_sample_index[sample]] + 1}") # adding one to make PLINK-compatible encoding 0 (no) -> 1 and 1 (yes) -> 2
                ofile.write("\n")
        else:
            for i, sample in enumerate(samples):
                ofile.write(f"{sample}")
                for pheno in phenos:
                    ofile.write(f"\t{pheno[i] + 1}") # adding one to make PLINK-compatible encoding 0 (no) -> 1 and 1 (yes) -> 2
                ofile.write("\n")
            

    
