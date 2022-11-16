rm(list=ls())
setwd("~/Desktop/Cancer_DRP/R/Prepare_Data/")


library(PharmacoGx)
availablePSets()

#GDSC_v2
GDSC_v2 = downloadPSet("GDSC_2020(v2-8.2)")
saveRDS(GDSC_v2, "Processed_data/S34/GDSC_v2_PharmacoGx.rds")


#GDSC
GDSC = downloadPSet("GDSC_2020(v1-8.2)")
saveRDS(GDSC, "Processed_data/S34/GDSC_PharmacoGx.rds")

#gCSI
gCSI = downloadPSet("gCSI_2019")
saveRDS(gCSI, "Processed_data/S34/gCSI_PharmacoGx.rds")

#CTRPv2
CTRPv2 = downloadPSet("CTRPv2_2015")
saveRDS(CTRPv2, "Processed_data/S34/CTRPv2_PharmacoGx.rds")



