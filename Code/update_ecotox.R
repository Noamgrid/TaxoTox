# =============================================================================
# update_ecotox.R
# -----------------------------------------------------------------------------
# Purpose : Download the latest ECOTOX Knowledgebase release from the EPA
#           website and rebuild the local SQLite cache used by TaxoTox.
#
# Usage   : Run once after initial setup, then re-run whenever a new ECOTOX
#           release is available (approximately quarterly).
#           After this script completes, re-run taxotox_install.R to
#           regenerate Data/final_ecotox_data.fst from the updated database.
#
# Note    : Requires an internet connection. The download can be several
#           hundred MB and may take 10–30 minutes depending on bandwidth.
#           The cache location is managed by ECOTOXr:
#             ECOTOXr::get_ecotox_sqlite_file()
# =============================================================================

library(ECOTOXr)

download_ecotox_data()

cat("ECOTOX database update complete.\n")
