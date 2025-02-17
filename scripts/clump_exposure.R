library(pacman)
p_load(tidyverse, data.table, ieugwasr, genetics.binaRies)
#p_load_gh("nicolash2/ggdendroplot")
setwd("/mnt/projects/users/pd_mr/scripts/")
p_load_gh("MRCIEU/TwoSampleMR")

clump_snps <- function(exposure_df, fp) {
    res <- tibble()
    for (id in unique(exposure_df$id.exposure)){
        exp_data <- exposure_df %>%
            dplyr::select(SNP, pval.exposure, id.exposure) %>%
            distinct() %>%
            filter(id.exposure == id) %>%
            rename("rsid" = "SNP",
                   "pval" = "pval.exposure",
                   "id" = "id.exposure")
        res <- res %>%
        bind_rows(tryCatch({
            ieugwasr::ld_clump(
                exp_data,
                bfile = "/mnt/projects/users/parkinson_comorbidities/data/EUR/EUR",
                plink_bin = genetics.binaRies::get_plink_binary()
            )}, error = function(e){}))
    }
    res <- res %>%
        rename("SNP" = "rsid",
               "pval.exposure" = "pval",
               "id.exposure" = "id") %>%
        inner_join(exposure_df,
                   by = c("SNP", "pval.exposure", "id.exposure")) %>%
        write_tsv(file.path(fp, "blood_grn_ld_clumped.txt"))
    # clump_data(df, clump_kb = 10000, clump_r2 = 0.001) # use this for API query
}

out_file_path <- '/mnt/projects/users/pd_mr/analysis/blood'
exposure_data <- read_tsv('/mnt/projects/users/pd_mr/analysis/blood/exposure_df.txt')
clump_expo_snps <- clump_snps(
                        exposure_data,
                        out_file_path
                        )
