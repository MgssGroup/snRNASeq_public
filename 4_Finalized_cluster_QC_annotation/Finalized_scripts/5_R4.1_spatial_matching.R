library(spatialLIBD)

ehub <- ExperimentHub::ExperimentHub()
spe <- fetch_data("spe")

spe_151507 <- spe[, spe$sample_id == "151507"]
counts_151507 <- counts(spe_151507)
coords_151507 <- spatialData(spe_151507)
coords_151507 <- coords_151507[,2:3]
layers_151507 <- spe_151507$layer_guess_reordered_short

save(counts_151507, coords_151507, layers_151507, file = "/home/malosree/scratch/Maynard_data/section_151507.Rda")

spe_151673 <- spe[, spe$sample_id == "151673"]
counts_151673 <- counts(spe_151673)
coords_151673 <- spatialData(spe_151673)
coords_151673 <- coords_151673[,2:3]
layers_151673 <- spe_151673$layer_guess_reordered_short

save(counts_151673, coords_151673, layers_151673, file = "/home/malosree/scratch/Maynard_data/section_151673.Rda")


print(summary(warnings()))
sessionInfo()
