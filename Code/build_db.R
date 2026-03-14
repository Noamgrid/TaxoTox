# =============================================================================
# build_db.R
# -----------------------------------------------------------------------------
# Purpose : Build the ECOTOX SQLite database from previously downloaded raw
#           data files (ASCII release).  Use this script only if you already
#           have the raw ECOTOX files on disk and want to (re)build the SQLite
#           cache without re-downloading.
#
#           For a combined download + build, use update_ecotox.R instead.
#
# Input   : ../Data/ecotox_data/  — directory containing ECOTOX ASCII files
#           (results.txt, tests.txt, species.txt, chemicals.txt, etc.)
# Output  : ../Data/ecotox.sqlite — SQLite database for use with ECOTOXr
#
# Note    : The path scheme here places raw data and the SQLite file inside
#           the project Data/ folder, which differs from the default ECOTOXr
#           cache location used by taxotox_install.R.  Adjust the connection
#           string in taxotox_install.R if you use this output file.
# =============================================================================

library(ECOTOXr)

build_ecotox_sqlite(
  source      = "../Data/ecotox_data",
  destination = "../Data/ecotox.sqlite"
)

cat("ECOTOX SQLite database build complete.\n")
