library(readxl)
library(BioMonTools)
library(BCGcalc)
library(dplyr)

# Set bioassessment directory and read in data
b_dir   <- "bioassessments_2026/bugs/"
f_dir   <- paste0(b_dir, "data/macroInvert_2021_2022_112125.xlsx")
s_dir   <- paste0(b_dir, "data/macroInvert_2021_2022_sites.xlsx")

df_samps_bugs <- readxl::read_excel(f_dir, guess_max = 10^6)
df_sites_bugs  <- readxl::read_excel(s_dir, guess_max = 10^6)

myDF <- df_samps_bugs

# Calculate Metrics
df_met_val_bugs <- BioMonTools::metric.values(myDF, "bugs")

write.table(df_met_val_bugs, 
            file.path(paste0(b_dir,"results/Metric_Values.tsv")),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

# Import Rules
df_rules <- readxl::read_excel(system.file("extdata/Rules.xlsx", 
                                           package = "BCGcalc"), 
                               sheet="Rules") 

# Calculate Metric Membership
df_met_memb <- BCG.Metric.Membership(df_met_val_bugs, df_rules)

write.table(df_met_memb, 
            file.path(paste0(b_dir,"results/Metric_Membership.tsv")),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

# Calculate Level Memberships
df_lev_memb <- BCG.Level.Membership(df_met_memb, df_rules)

write.table(df_lev_memb,
            file.path(paste0(b_dir,"results/Level_Membership.tsv")),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

# Calculate Level Assignment
df_levels <- BCG.Level.Assignment(df_lev_memb)

write.table(df_levels,
            file.path(paste0(b_dir,"results/Level_Assignment.tsv")),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

# Merge with site data and organize for assessments
df_levels_site <- merge(df_levels,df_sites_bugs,by.x="SampleID",by.y="SAMPLEID")

siteBCG_sum <-  group_by(df_levels_site,staSeq) %>%
                summarize(avgBCG = mean(Continuous_BCG_Level),
                minBCG = min(Continuous_BCG_Level),
                maxBCG = max(Continuous_BCG_Level),
                minBCG1 = min(Primary_BCG_Level),
                maxBCG1 = max(Primary_BCG_Level),
                nBCG   = n())
siteBCG_sum <-  as.data.frame(merge(siteBCG_sum,unique(df_sites_bugs[4:6]),
                                    by="staSeq"))
siteBCG_sum$s_assessment <- ifelse(siteBCG_sum$avgBCG<=4.4,"F",
                                   ifelse(siteBCG_sum$avgBCG>=5,"N","A"))

write.csv(siteBCG_sum, paste0(b_dir,"results/siteBCG_sum.csv"),row.names=FALSE)

