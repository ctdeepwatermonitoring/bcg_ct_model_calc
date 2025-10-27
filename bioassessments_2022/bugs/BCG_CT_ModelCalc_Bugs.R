# CT BCG Model Calculation
# Bugs
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
library(lazyeval)
library(knitr)

# Functions
## Given a dataframe with phylogenetic classes return lowest class in 'rank' col
p_class <- function(data,class_l){
  
  for(i in 1:dim(data)[1]){
    sample = data[i,]
    
    for(rank in class_l){
      if(nchar(sample[rank])==0 | is.na(sample[rank])){
        data[i,"Rank"] = data[i,"Rank"]
      }
      else{
        data[i,"Rank"] = rank
      }
    }
  }
  return(data)
}

# Global
setwd("C:/Users/deepuser/Documents/Projects/BCG/Bioassessments2022/bugs/")
myCommunity <- "bugs"
myDateTime <- format(Sys.time(), "%Y%m%d_%H%M%S")
dn_output <- getwd()
mySiteType <- "bug01"

# Data, BCG Rules
df_rules <- read_excel(system.file("./extdata/Rules.xlsx", package="BCGcalc")
                       , sheet = "Rules") 

# METRIC VALUES ----
# Data, Sample Taxa, SiteInfo
# Data
f_dir <- "/Data/macroInvert_2016_2020.xlsx"
s_name<- "/Data/macroInvert_2016_2020_sites.xlsx"
fn_SampTaxa <- paste0(dn_output,f_dir)
fn_Sites    <- paste0(dn_output,s_name)
df_SampTaxa <- as.data.frame(read_excel(path = fn_SampTaxa, guess_max = 10^6))
df_Sites    <- as.data.frame(read_excel(path = fn_Sites, guess_max = 10^6))

# Update phylogenetic class
rank_l = colnames(df_SampTaxa)[9:21]
df_SampTaxa <- p_class(df_SampTaxa,rank_l)

# Mark Excluded Taxa
df_SampTaxa <- markExcluded(df_SampTaxa
                            , SampID="SampleID"
                            , TaxaID="TAXAID"
                            , TaxaCount = "N_TAXA"
                            , Exclude="EXCLUDE"
                            , TaxaLevels=rank_l
                            , Exceptions=NA)

# Munge
names(df_SampTaxa) <- toupper(names(df_SampTaxa))
df_SampTaxa[, "N_TAXA"] <- as.numeric(df_SampTaxa[, "N_TAXA"])
if(myCommunity == "fish") {
  df_SampTaxa[, "N_ANOMALIES"] <- as.numeric(0)
  df_SampTaxa[, "SAMP_WIDTH_M"] <- as.numeric(0)
  df_SampTaxa[, "SAMP_LENGTH_M"] <- as.numeric(1)
}## IF ~ myCommunity ~ END

# Metric Calc
# 1 missing field for metric that is not calculated UFC Select 'YES'
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

df_Levels_site <- merge(df_Levels,df_Sites,by.x="SampleID",by.y="SAMPLEID")

siteBCG_sum <-  group_by(df_Levels_site,staSeq) %>%
                  summarize(avgBCG = mean(Lev.Prop.Num),
                            minBCG = min(Lev.Prop.Num),
                            maxBCG = max(Lev.Prop.Num),
                            minBCG1 = min(Lev.1.Name),
                            maxBCG1 = max(Lev.1.Name),
                            nBCG   = n())
siteBCG_sum <-  as.data.frame(merge(siteBCG_sum,unique(df_Sites[4:6]),by="staSeq"))
siteBCG_sum$s_assessment <- ifelse(siteBCG_sum$avgBCG<=4.4,"F",
                                   ifelse(siteBCG_sum$avgBCG>=5,"N","A"))

write.csv(siteBCG_sum,"siteBCG_sum.csv",row.names=FALSE)