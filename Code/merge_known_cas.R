setwd("c:/Users/Yair/Documents/yair/research/estuary_rehabilitation/yairs_stuff/Students/Noam/TaxoTox/Data")
library(fst)
library(data.table)

k <- read.fst("Known_CAS.fst", as.data.table = TRUE)
t <- read.fst("temp_CAS.fst",  as.data.table = TRUE)

# ── Combine and remove exact duplicate rows ───────────────────────────────
combined <- rbindlist(list(k, t), use.names = TRUE, fill = TRUE)
combined  <- unique(combined, by = c("PREFERRED_NAME", "CASRN"))

# ── Resolve name → many-CASRN conflicts ──────────────────────────────────
# For conflicts, keep Known_CAS CASRN by default, except for the three
# entries below where a specific correct CASRN has been confirmed.
correct_casrn <- c(
    "2-Aminobenzimidazole" = "934-32-7",
    "Azinphos-methyl oxon" = "961-22-8",
    "Disulfoton sulfoxide"  = "2497-07-6"
)

# Find all names that still have >1 CASRN
multi <- combined[, .N, by = PREFERRED_NAME][N > 1, PREFERRED_NAME]

for (nm in multi) {
    if (nm %in% names(correct_casrn)) {
        keep <- correct_casrn[[nm]]
        cat("Explicit override for", nm, "->", keep, "\n")
    } else {
        keep <- k[PREFERRED_NAME == nm, CASRN][1]
        cat("Keeping Known_CAS   for", nm, "->", keep, "\n")
    }
    combined <- combined[!(PREFERRED_NAME == nm & CASRN != keep)]
}

# ── Final dedup sanity check ──────────────────────────────────────────────
remaining_multi <- combined[, .N, by = PREFERRED_NAME][N > 1]
if (nrow(remaining_multi) > 0) {
    cat("\nWARNING: still has multi-CASRN names:\n")
    print(remaining_multi)
    stop("Aborting — resolve conflicts before saving.")
}

cat("\nFinal row count:", nrow(combined), "\n")
cat("(Known_CAS had", nrow(k), "/ temp_CAS had", nrow(t), ")\n")

# ── Back up old Known_CAS and save merged result ──────────────────────────
file.copy("Known_CAS.fst", "Known_CAS_backup.fst", overwrite = TRUE)
cat("Backup saved to Known_CAS_backup.fst\n")

write.fst(combined, "Known_CAS.fst")
cat("Known_CAS.fst updated.\n")
