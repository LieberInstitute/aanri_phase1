## This script randomly assigns variable of interest and performs
## limma-voom and mash modeling to get effect sizes. It will
## require that entrez or ensembl ids are provided in the annotation
## file.

library(optparse)

fit_model <- function(counts_file, covs_file, annot_file, target, perm_num, label){
                                        # Load data
    counts <- read.table(counts_file, sep='\t', header=TRUE, row.names=1)
    annot  <- read.table(annot_file, sep='\t', header=TRUE, row.names=1)
    covs   <- read.table(covs_file, sep='\t', header=TRUE, row.names=1)
                                        # Subset features and individuals
    gene_ids   <- intersect(rownames(counts), rownames(annot))
    sample_ids <- intersect(colnames(counts), rownames(covs))
    counts     <- counts[gene_ids, sample_ids]
    annot      <- annot[gene_ids, ]
    covs       <- covs[sample_ids, ]
                                        # Permutate target phenotype
    seed = perm_num + 13131313
    set.seed(seed)
    covs[,target] <- sample(covs[,target])
                                        # Generate edgeR obeject
    x <- edgeR::DGEList(counts=counts, genes=annot, samples=covs)
    x <- edgeR::calcNormFactors(x, method="TMM")
                                        # Extract model
    covs_names <- names(x$samples)[ !names(x$samples) %in%
                                    c("group", "lib.size", "norm.factors")]
    mod <- x$samples[ ,covs_names]
                                        # Build contrasts
    commandstr = paste0("limma::makeContrasts(phenotype=",
                        target,",levels=colnames(mod))")
    contr.matrix <- eval(parse(text=commandstr))
                                        # Fit model
    v            <- limma::voom(x, mod)
    fit0         <- limma::lmFit(v, mod)
    fit <- limma::contrasts.fit(fit0, contrasts=contr.matrix)
    esv <- limma::eBayes(fit)
                                        # Extract data
    top        <- limma::topTable(esv, coef=1, number=Inf, sort.by="P")
    top[,"SE"] <- sqrt(esv$s2.post) * as.vector(esv$stdev.unscaled)
                                        # Save data
    df <- data.frame("feature_ids"=rownames(top), "bhat"=top$logFC, "shat"=top$SE)
    write.table(df, file=paste0(label,"_model_bhat_N_shat_perm",perm_num,".txt"),
                sep="\t", row.names=FALSE, quote=FALSE)
}

#### MAIN ####
                                        # parser
option_list <- list(
    make_option(c("--counts"), type="character", default=NULL,
                help="counts file", metavar="character"),
    make_option(c("--covariates"), type="character", default=NULL,
                help="covariates file", metavar="character"),
    make_option(c("--annotation"), type="character", default=NULL,
                help="feature annotation file"),
    make_option(c("-t", "--target"), type="character", default=NULL,
                help="phenotype of interest in the covariate files"),
    make_option(c("-p", "--perm_num"), type="integer", default=1,
                help="permutation number [default=%default]"),
    make_option(c("-o", "--output"), type="character", default="output",
                help="prefix for output files [default=%default]")
)
opt_parser  <- OptionParser(option_list=option_list)
opt         <- parse_args(opt_parser)
                                        # Check for variables
if(is.null(opt$counts)){
    print_help(opt_parser)
    stop("a counts file (tsv) must be supplied.", call.=FALSE)
} else if(is.null(opt$covariates)){
    print_help(opt_parser)
    stop("a covariates file (tsv) must be supplied.", call.=FALSE)
} else if(is.null(opt$annotation)){
    print_help(opt_parser)
    stop("an annotation file (tsv) must be supplied.", call.=FALSE)
} else if(is.null(opt$target)){
    print_help(opt_parser)
    stop("phenotype of interest (target) must be supplied.", call.=FALSE)
}

outdir = "temp_files/"
dir.create(outdir, showWarnings=FALSE)
fit_model(opt$counts, opt$covariates, opt$annotation, opt$target,
          opt$perm_num, paste0(outdir, opt$output))
