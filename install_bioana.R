# Install Bioana and the packages it uses.
#   Rscript install_bioana.R            # analysis profile: engines, readers, genetics
#   Rscript install_bioana.R full       # + plots and test tools
#   Rscript install_bioana.R suite      # + Simitall
# Run from the bioana folder. The package list is inst/rules/requirements.tsv.
args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args)) args[1] else "analysis"
options(repos = c(CRAN = "https://cloud.r-project.org"))
for (p in c("BiocManager", "remotes", "jsonlite")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
remotes::install_local(".", upgrade = "never", force = TRUE)
bioana::bioana_install_dependencies(profile)
r <- bioana::bioana_check_dependencies(profile)
print(r[, c("package", "source", "profile", "installed", "version")], row.names = FALSE)
