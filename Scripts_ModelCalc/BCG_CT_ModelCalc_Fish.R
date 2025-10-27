# CT BCG Model Calculation
# Fish
#
# Erik.Leppo@tetratech.com
# 2021-07-31
#~~~~~~~~~~~~~~~~~~~~~~~~~

# Packages, Install (if needed)
## Option 1
# if(!require(remotes)){install.packages("remotes")}
# if(!require(BioMonTools)){install_github("leppott/BioMonTools")}
# if(!require(BCGcalc)){install_github("leppott/BCGcalc")}
# if(!require(readxl)){install.packages("readxl")}
# if(!require(dplyr)){install.packages("dplyr")}

# Packages
## Option 2
library(BioMonTools)
library(BCGcalc)
library(readxl)
library(dplyr)

# Global
myCommunity <- "fish"
myDateTime <- format(Sys.time(), "%Y%m%d_%H%M%S")
dn_output <- getwd()
mySiteType <- c("fish01", "fish02", "fish03")

# Data, BCG Rules
df_rules <- read_excel(system.file("./extdata/Rules.xlsx", package="BCGcalc")
                       , sheet = "Rules") 

# METRIC VALUES ----
# Data, Sample Taxa
# Data
fn_SampTaxa <- "Export_R_MetricCalc_Fish.xlsx"
df_SampTaxa <- as.data.frame(read_excel(path = fn_SampTaxa, guess_max = 10^6))

# Munge
names(df_SampTaxa) <- toupper(names(df_SampTaxa))
df_SampTaxa[, "N_TAXA"] <- as.numeric(df_SampTaxa[, "N_TAXA"])
if(myCommunity == "fish") {
  df_SampTaxa[, "N_ANOMALIES"] <- as.numeric(0)
  df_SampTaxa[, "SAMP_WIDTH_M"] <- as.numeric(0)
  df_SampTaxa[, "SAMP_LENGTH_M"] <- as.numeric(1)
}## IF ~ myCommunity ~ END

# Metric Calc
df_Met_Val <- metric.values(fun.DF = df_SampTaxa
                            , fun.Community = myCommunity)

# Save
fn_Met_Val <- file.path(dn_output
                        , paste0("Calc_Met_Val_"
                                 , myCommunity
                                 , "_"
                                 , myDateTime
                                 , ".tsv"))
write.table(df_Met_Val
            , file = fn_Met_Val
            , sep = "\t"
            , row.names = FALSE
            , col.names = TRUE)

# METRIC MEMBERSHIP ----

# Munge
## Match data and R function columns names
df_Met_Val[, "INDEX_REGION"] <- tolower(df_Met_Val[, "INDEX_REGION"])
df_Met_Val[, "SITE_TYPE"] <- df_Met_Val[, "INDEX_REGION"]
## Filter for only relevant index regions
df_Met_Val <- filter(df_Met_Val, SITE_TYPE %in% mySiteType)
## Rules names to upper case
names(df_rules) <- toupper(names(df_rules))
## Filter for only relevant metrics
df_rules <- filter(df_rules, SITE_TYPE %in% mySiteType)

# Calc Metric Membership
df_Met_Mem <- BCG.Metric.Membership(df_Met_Val, df_rules)

# Save
fn_Met_Mem <- file.path(dn_output
                        , paste0("Calc_Met_Mem_"
                                 , myCommunity
                                 , "_"
                                 , myDateTime
                                 , ".tsv"))
write.table(df_Met_Mem
            , file = fn_Met_Mem
            , sep = "\t"
            , row.names = FALSE
            , col.names = TRUE)

# LEVEL MEMBERSHIP ----

# Calc Level Membership
df_Lev_Mem <- BCG.Level.Membership(df_Met_Mem, df_rules)

# Save
fn_Lev_Mem <- file.path(dn_output
                        , paste0("Calc_Lev_Mem_"
                                 , myCommunity
                                 , "_"
                                 , myDateTime
                                 , ".tsv"))
write.table(df_Lev_Mem
            , fn_Lev_Mem
            , sep = "\t"
            , row.names = FALSE
            , col.names = TRUE)

# LEVEL ASSIGNMENT ----

# Calc Level Membership
df_Levels <- BCG.Level.Assignment(df_Lev_Mem)

# Save
fn_Levels <- file.path(dn_output
                       , paste0("Calc_Levels"
                                , myCommunity
                                , "_"
                                , myDateTime
                                , ".tsv"))
write.table(df_Levels
            , fn_Levels
            , sep = "\t"
            , row.names = FALSE
            , col.names = TRUE)