## code to prepare `cps_moisture` dataset goes here

# Central Park soil data at the family level, as prepared for the trac paper:
# running 0create_phyloseq_object.R and 1prep_data_all_levels.R in the
# CentralParkSoil folder of the trac paper's reproducibility code
# (https://github.com/jacobbien/trac-reproducible) creates
# cps_Mois_aggregated.RDS. Its "Family" element holds the counts summed to the
# family level with the matching tree, and the OTU element holds the sample data.
cps_file <- Sys.getenv("CPS_MOIS_AGGREGATED",
                       "../trac-reproducible-main/CentralParkSoil/cps_Mois_aggregated.RDS")
dat <- readRDS(cps_file)
family <- dat$Family
sample_data <- as.data.frame(dat$OTU$sample_data)
stopifnot(identical(as.numeric(sample_data$Moisture), family$y))
stopifnot(identical(colnames(family$x), rownames(family$A)))

x <- as.matrix(family$x)
stopifnot(all(x == round(x)))
storage.mode(x) <- "integer"
rownames(x) <- sample_data$X.SampleID

cps_moisture <- list(
  y = family$y,
  x = x,
  tree = family$tree,
  tax = family$tax,
  A = family$A,
  covariates = data.frame(pH = sample_data$pH, C = sample_data$C,
                          N = sample_data$N, CN = sample_data$CN,
                          CO2_C = sample_data$CO2_C,
                          row.names = sample_data$X.SampleID)
)
usethis::use_data(cps_moisture, overwrite = TRUE, version = 2, compress = "xz")
