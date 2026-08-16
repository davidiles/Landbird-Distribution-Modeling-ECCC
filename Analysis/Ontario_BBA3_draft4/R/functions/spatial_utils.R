# ============================================================
# spatial_utils.R
#
# Shared spatial helpers for the Ontario_BBA workflow.
#
# Contents
#   get_aea_km_crs   Default kilometre-based projected CRS for spatial modelling.
#
# Sourced by 01 and 02.
# ============================================================

# Return the default kilometre-based projected CRS used for spatial modelling.
# A Lambert conformal conic in km units keeps SPDE mesh edge lengths, buffers,
# and simplification tolerances all expressed directly in kilometres.
get_aea_km_crs <- function() {
  # Single-line PROJ string to avoid newline/whitespace surprises.
  proj <- paste(
    "+proj=lcc",
    "+lat_1=49",
    "+lat_2=77",
    "+lat_0=49",
    "+lon_0=-95",
    "+datum=NAD83",
    "+units=km",
    "+no_defs"
  )
  sf::st_crs(proj)
}
