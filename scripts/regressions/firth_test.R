#!/usr/bin/env Rscript

library(optparse)
library(vcfR)
library(dplyr)
library(logistf)

option_list <- list(
	make_option(c("-v", "--vcf"), type = "character", default = NULL, help = "Path to the input VCF", metavar = "character"),
	make_option(c("-p", "--pheno"), type = "character", default = NULL, help = "Path to the PLINK2-compatible phenotypes file.", metavar = "character"),
	make_option(c("-c", "--covar"), type = "character", default = NULL, help = "Path to the PLINK2-compatible covariates file.", metavar = "character"),
	make_option(c("-o", "--output"), type = "character", default = NULL, help = "Prefix for output files.", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list, description = "This script implements a Firth logistic regression.")
opt <- parse_args(opt_parser)

required_pheno_columns = c("IID")
required_covar_columns = c("IID")

if (is.null(opt$vcf)) {
	print_help(opt_parser)
	stop("Error: Missing required argument '-v'.\n", call. = FALSE)
}

if (is.null(opt$pheno)) {
	print_help(opt_parser)
	stop("Error: Missing required argument '-p'.\n", call. = FALSE)
}

if (is.null(opt$output)) {
	print_help(opt_parser)
	stop("Error: Missing required argument '-o'.\n", call. = FALSE)
}


vcf_filename <- opt$vcf
pheno_filename <- opt$pheno
covar_filename <- opt$covar
output_prefix <- opt$output

pheno <- read.table(pheno_filename, header = TRUE)
if (!all(required_pheno_columns %in% names(pheno))) {
	stop("Phenotype file missing required columns.")
}
pheno_names <- setdiff(colnames(pheno), "IID")

if (anyDuplicated(pheno$IID) != 0) {
        stop("Phenotype file has duplicated IID values.")
}

dat <- pheno

covar_terms <- character(0)

if (!is.null(opt$covar)) {
	covar <- read.table(covar_filename, header = TRUE)
	if (!all(required_covar_columns %in% names(covar))) {
		stop("Covariates file missing required columns.")
	}

        if (!setequal(pheno$IID, covar$IID)) {
                stop("Phenotype and covariate files must contain identical IID sets.")
        }

        if (any(setdiff(names(pheno), "IID") %in% names(covar))) {
                stop("Phenotype and covariate files have at least one column with same name.")
        }

        dat <- merge(pheno, covar, by = "IID", all = FALSE, sort = FALSE)
        if (nrow(dat) != nrow(pheno)) {
                stop("Samples in phenotype and covariate files differ.")
        }

        # Check for PC1, PC2, PC3, ... and standardize
        pc_cols <- grep("^PC[0-9]+$", names(dat), value = TRUE)
        if (length(pc_cols) > 0) {
                dat[pc_cols] <- scale(dat[pc_cols])
        } else {
                message("No PC covariates detected")
        }

	# Check for sex column and factorize
        sex_col <- grep("^sex$", names(dat), ignore.case = TRUE, value = TRUE)
        if (length(sex_col) > 1) {
                stop("More than one sex covariate detected.")
        } else if (length(sex_col) == 0) {
                message("No sex covariate detected.")
        } else {
                dat[[sex_col]] <- as.factor(as.character(dat[[sex_col]]))
        }

	covar_terms <- setdiff(names(covar), "IID")
        if (length(covar_terms) == 0) {
                message("No covariates detected")
        } else {
                message("Covariates detected: ", paste(covar_terms, collapse = ", "))
        }
}

vcf <- read.vcfR(vcf_filename)
gt <- extract.gt(vcf, element = "GT")

if (!setequal(colnames(gt), dat$IID)) {
	stop("Samples in VCF file differ from those in phenotype/covariate files.")
}

# Re-order samples to match VCF
dat <- dat[match(colnames(gt), dat$IID), ]

if (!identical(dat$IID, colnames(gt))) {
        stop("Something went wrong when matching sample order in VCF with phenotype/covariate file.")
}

code_gt <- function(g) {
        if (is.na(g) || g == "./." || g == ".|.") return(NA_integer_)
        alleles <- strsplit(g, "[/|]")[[1]]
        sum(as.integer(alleles), na.rm = TRUE)
}

# Right hand of the regression expression stays same for every phenotype
rhs <- c("G", covar_terms)

for (pheno_name in pheno_names) {
	message("Running phenotype: ", pheno_name)
	message("Regression: ", paste(pheno_name, "~", paste(rhs, collapse = " + ")))

        # Regression formula for the phenotype
        form <- as.formula(paste(pheno_name, "~", paste(rhs, collapse = " + ")))

	# Recode phenotype from 1/2 to 0/1
	dat[[pheno_name]] <- ifelse(dat[[pheno_name]] == 2, 1, 0)

	results <- vector("list", nrow(gt))

	for (i in seq_len(nrow(gt))) {
		allele <- rownames(gt)[i]
		g <- vapply(gt[i, ], code_gt, integer(1))

   		keep <- !is.na(dat[[pheno_name]]) & !is.na(g)

                if (sum(keep) == 0) {
                    stop(paste0("Fatal error: no overlapping non-missing samples for phenotype '", pheno_name, "' and current variant. Check IID alignment / missingness."))
                }

                # drop missing phenotype or genotype
                d_model <- dat[keep, , drop = FALSE]
                d_model$G <- g[keep]

                ac <- sum(d_model$G)
                n_alleles <- 2 * nrow(d_model)
                mac <- min(ac, n_alleles - ac)

                if (n_alleles == 0) {
                        af <- NA_real_
                        maf <- NA_real_
                } else {
                        af <- ac / n_alleles
                        maf <- min(af, 1 - af)
                }
	
		if (mac == 0) {
                        results[[i]] <- data.frame(
                                output = output_prefix,
                                pheno = pheno_name,
                                allele = allele,
                                ac = ac,
                                mac = 0,
                                af = af,
                                maf = maf,
                                n = nrow(d_model),
                                beta = NA_real_,
                                se = NA_real_,
                                p = NA_real_,
                                error = "Monomorphic",
                                stringsAsFactors = FALSE
                        )
                        next
		}
		
		fit <- try(
    			withCallingHandlers(
				logistf(
					form,
         				data = d_model,
					pl = TRUE,
					maxit = 250,
    				),
				warning = function(w) {
            				message("WARNING (", allele, "): ", conditionMessage(w))
					invokeRestart("muffleWarning")
           			}
    			),
			silent = TRUE
  		)

  		if (inherits(fit, "try-error")) {
			err <- conditionMessage(attr(fit, "condition"))
			err <- sub("\\n.*", "", err)   # first line only
                
			results[[i]] <- data.frame(
                                output = output_prefix,
                                pheno = pheno_name,
                                allele  = allele,
                                ac = ac,
                                mac = mac,
                                af = af,
                                maf = maf,
                                n = nrow(d_model),
                                beta = NA_real_,
                                se = NA_real_,
                                p =  NA_real_,
                                error = err,
                                stringsAsFactors = FALSE
                        )
		} else {
			coef_idx <- match("G", names(fit$coefficients))
			se_vec <- sqrt(diag(fit$var))
            		se_val <- if (!is.na(coef_idx) && coef_idx <= length(se_vec)) { se_vec[coef_idx] } else NA_real_			
			results[[i]] <- data.frame(
				output = output_prefix,
				pheno = pheno_name,
				allele  = allele,
				ac = ac,
				mac = mac,
				af = af, 
				maf = maf,
				n = nrow(d_model),
				beta = fit$coefficients[coef_idx],
				se = sqrt(diag(fit$var))[coef_idx],
				p =  fit$prob[coef_idx],
				error = NA_character_,
				stringsAsFactors = FALSE
  			)
		}
	}

	out <- dplyr::bind_rows(results)

    	write.table(
        	out,
        	file = paste0(output_prefix, ".", pheno_name, ".firth.tsv"),
        	sep = "\t",
        	quote = FALSE,
        	row.names = FALSE
    	)
}
