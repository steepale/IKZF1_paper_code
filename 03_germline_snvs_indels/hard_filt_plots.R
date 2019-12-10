library("ggplot2")

# In R, Generate the QD Plot
samples_QD = read.table("./data/hard_filtered_variants/collective_samples_raw_snps_extracted_QD.txt")
qplot(V1, data = samples_QD, geom = "density") + xlab("QD") + ylab("Density")
ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_QD_density_snps.png")

# In R, Generate the FS Plot (x-axis is logged)
samples_FS = read.table("./data/hard_filtered_variants/collective_samples_raw_snps_extracted_FS.txt")
qplot(V1, data = samples_FS, geom = "density", log="x") + xlab("FS") + ylab("Density")
ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_FS_density_snps.png")
# Cutoff right after right peak

# In R, Generate the MQ Plot
samples_MQ = read.table("./data/hard_filtered_variants/collective_samples_raw_snps_extracted_MQ.txt")
qplot(V1, data = samples_MQ, geom = "density") + xlab("MQ") + ylab("Density")
ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_MQ_density_snps.png")

# Zoom in on the MQ peak at 40
qplot(V1, data = samples_MQ, geom = "density") + xlab("MQ") + ylab("Density") + xlim(39.0, 41.0)
ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_MQ_density_40_snps.png")

# Zoom in on the MQ peak at 60
qplot(V1, data = samples_MQ, geom = "density") + xlab("MQ") + ylab("Density") + xlim(59.0, 61.0)
ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_MQ_density_60_snps.png")
# Discard anything that is not 60

# In R, Generate the MQRankSum Plot
samples_MQRankSum = read.table("./data/hard_filtered_variants/collective_samples_raw_snps_extracted_MQRankSum.txt")
qplot(V1, data = samples_MQRankSum, geom = "density") + xlab("MQRankSum") + ylab("Density")
ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_MQRankSum_density_snps.png")
# Remove anything less than -2

# In R, Generate the ReadPosRankSum Plot
samples_ReadPosRankSum = read.table("./data/hard_filtered_variants/collective_samples_raw_snps_extracted_ReadPosRankSum.txt")
qplot(V1, data = samples_ReadPosRankSum, geom = "density") + xlab("ReadPosRankSum") + ylab("Density")
ggsave("./analysis/hard_filtered_variants/R_plots/collective_samples_ReadPosRankSum_density_snps.png")
