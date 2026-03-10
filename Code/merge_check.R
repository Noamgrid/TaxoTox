setwd("c:/Users/Yair/Documents/yair/research/estuary_rehabilitation/yairs_stuff/Students/Noam/TaxoTox/Data")
library(fst)
library(data.table)

k <- read.fst("Known_CAS.fst", as.data.table = TRUE)
t <- read.fst("temp_CAS.fst",  as.data.table = TRUE)

combined <- rbindlist(list(k[, source := "Known_CAS"],
                            t[, source := "temp_CAS"]),
                       use.names = TRUE, fill = TRUE)
combined  <- unique(combined, by = c("PREFERRED_NAME", "CASRN"))

multi_cas <- combined[, .(n = .N, CASRNs  = paste(CASRN,  collapse = " | "),
                                   Sources = paste(source, collapse = " | ")),
                       by = PREFERRED_NAME][n > 1]

cat("=== Names with multiple CASRNs (need decision) ===\n")
print(multi_cas[, .(PREFERRED_NAME, CASRNs, Sources)])
