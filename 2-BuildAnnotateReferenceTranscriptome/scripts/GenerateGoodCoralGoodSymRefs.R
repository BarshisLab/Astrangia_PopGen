#This script reads in different parsed cutoffs (with ReParseBlastbycutoffs.py), which are output from blasting the Trinity-assembled reference against dirty coral,
#clean sym, dirty sym, and clean coral databases. 

#Last updated by Hannah Aichelman (06/2019)


#set working directory
setwd("~/Documents/ODU_MS/MolecularWork/RNASeq_Analysis/MyTrinityAssembly")
getwd()

#read in different parsed outputs and explore the number of contigs in each
dirtycoral95 <- read.delim("TrinityMinusSilva_DC_parsed_95per100bpmatch.txt")
# > dim(dirtycoral95)
# [1] 1117   16
cleansym60 <- read.delim("TrinityMinusSilva_CS_parsed_60per100bpmatch.txt")
# > dim(cleansym60)
# [1] 1779   16
dirtysym95 <- read.delim("TrinityMinusSilva_DS_parsed_95per100bpmatch.txt")
# > dim(dirtysym95)
# [1] 954  16
cleancoral60 <- read.delim("TrinityMinusSilva_CC_parsed_60per100bpmatch.txt")
# > dim(cleancoral60)
# [1] 9383   16

#trying a 50bp match instead of 100bp to see if this helps increase number of contigs:

dirtycoral95_50 <- read.delim("TrinityMinusSilva_DC_parsed_95per50bpmatch.txt")
cleansym95_50 <- read.delim("TrinityMinusSilva_CS_parsed_95per50bpmatch.txt")
dirtysym95_50 <- read.delim("TrinityMinusSilva_DS_parsed_95per50bpmatch.txt")
cleancoral95_50 <- read.delim("TrinityMinusSilva_CC_parsed_95per50bpmatch.txt")

dirtycoral60_50 <- read.delim("TrinityMinusSilva_DC_parsed_60per50bpmatch.txt")
cleansym60_50 <- read.delim("TrinityMinusSilva_CS_parsed_60per50bpmatch.txt")
dirtysym60_50 <- read.delim("TrinityMinusSilva_DS_parsed_60per50bpmatch.txt")
cleancoral60_50 <- read.delim("TrinityMinusSilva_CC_parsed_60per50bpmatch.txt")

#Some comparisons of 50 vs. 100bp matches:
# > sum(!(dirtycoral95_50$Query.Name%in%cleansym60_50$Query.Name))
# [1] 271
# > sum(!(dirtycoral95$Query.Name%in%cleansym60$Query.Name))
# [1] 218


#Write out tables for the clean coral and clean sym references with as many contigs as we can get:
write.table(dirtycoral60[!(dirtycoral60$Query.Name%in%cleansym60$Query.Name), "Query.Name"], file="GoodCoralContig_60vs60_13381.txt", quote = F, row.names = F)

write.table(dirtysym60[!(dirtysym60$Query.Name%in%cleancoral60$Query.Name), "Query.Name"], file="GoodSymbiontContig_60vs60_1698.txt", quote = F, row.names = F)

#REMAINING DIRTY CORAL PARSED
dirtycoral50 <- read.delim("TrinityMinusSilva_DC_parsed_50per100bpmatch.txt")
# > dim(dirtycoral50)
# [1] 15105    16
dirtycoral60 <- read.delim("TrinityMinusSilva_DC_parsed_60per100bpmatch.txt")
# > dim(dirtycoral60)
# [1] 15105    16
dirtycoral70 <- read.delim("TrinityMinusSilva_DC_parsed_70per100bpmatch.txt")
# > dim(dirtycoral70)
# [1] 15105    16
dirtycoral80 <- read.delim("TrinityMinusSilva_DC_parsed_80per100bpmatch.txt")
# > dim(dirtycoral80)
# [1] 14241    16
dirtycoral90 <- read.delim("TrinityMinusSilva_DC_parsed_90per100bpmatch.txt")
# > dim(dirtycoral90)
# [1] 5401   16

#REMAINING CLEAN SYM PARSED
cleansym50 <- read.delim("TrinityMinusSilva_CS_parsed_50per100bpmatch.txt")
# > dim(cleansym50)
# [1] 1779   16
cleansym70 <- read.delim("TrinityMinusSilva_CS_parsed_70per100bpmatch.txt")
# > dim(cleansym70)
# [1] 1779   16
cleansym80 <- read.delim("TrinityMinusSilva_CS_parsed_80per100bpmatch.txt")
# > dim(cleansym80)
# [1] 1748   16
cleansym90 <- read.delim("TrinityMinusSilva_CS_parsed_90per100bpmatch.txt")
# > dim(cleansym90)
# [1] 1533   16
cleansym95 <- read.delim("TrinityMinusSilva_CS_parsed_95per100bpmatch.txt")
# > dim(cleansym95)
# [1] 934  16

#REMAINING DIRTY SYM PARSED
dirtysym50 <- read.delim("TrinityMinusSilva_DS_parsed_50per100bpmatch.txt")
# > dim(dirtysym50)
# [1] 2520   16
dirtysym60 <- read.delim("TrinityMinusSilva_DS_parsed_60per100bpmatch.txt")
# > dim(dirtysym60)
# [1] 2520   16
dirtysym70 <- read.delim("TrinityMinusSilva_DS_parsed_70per100bpmatch.txt")
# > dim(dirtysym70)
# [1] 2520   16
dirtysym80 <- read.delim("TrinityMinusSilva_DS_parsed_80per100bpmatch.txt")
# > dim(dirtysym80)
# [1] 2344   16
dirtysym90 <- read.delim("TrinityMinusSilva_DS_parsed_90per100bpmatch.txt")
# > dim(dirtysym90)
# [1] 1580   16

#REMAINING CLEAN CORAL PARSED
cleancoral50 <- read.delim("TrinityMinusSilva_CC_parsed_50per100bpmatch.txt")
# > dim(cleancoral50)
# [1] 9383   16
cleancoral70 <- read.delim("TrinityMinusSilva_CC_parsed_70per100bpmatch.txt")
# > dim(cleancoral70)
# [1] 9383   16
cleancoral80 <- read.delim("TrinityMinusSilva_CC_parsed_80per100bpmatch.txt")
# > dim(cleancoral80)
# [1] 8689   16
cleancoral90 <- read.delim("TrinityMinusSilva_CC_parsed_90per100bpmatch.txt")
# > dim(cleancoral90)
# [1] 1707   16
cleancoral95 <- read.delim("TrinityMinusSilva_CC_parsed_95per100bpmatch.txt")
# > dim(cleancoral95)
# [1] 96 16


