#!/usr/bin/env Rscript
# =============================================================================
# data_preparation.R
# -----------------------------------------------------------------------------
# Regenerates, from raw public sources, every data file that README3.Rmd reads:
#
#   roads.csv        (from, to, weight)            network edges (weight = time, minutes)
#   coordinates.csv  (ID, X, Y)                    node coordinates (EPSG:2154, meters)
#   roads_aux.csv    (LONGUEUR, length)            per-edge physical length (km)
#   node_commune.csv (com, id_noeud, POPULATION)   one row per commune + population 2021
#   doctor.csv       (CODGEO, NB_D201, id_noeud, POPULATION)  general practitioners (2025)
#   maternity.csv    (CODGEO, NB_D107)             maternity units (2025)
#   commuting.Rds    (object `cmt`: CODGEO, LIBGEO, DCLT, L_DCLT, flux)  home-work flows (2023)
#
# The file `com_simplified_geom.shp` (simplified commune polygons, used only for
# plotting the map) is not regenerated here; it is kept as-is in the repository.
#
# RAW SOURCES (downloaded separately -- raw data is NOT committed to the repo):
#
#   (A) Road network  -- ROUTE500 3.0 - IGN
#       projection LAMB93 / ED211 (EPSG:2154):
#         - TRONCON_ROUTE.shp  road segments (attributes VOCATION, ACCES, SENS, LONGUEUR)
#         - NOEUD_COMMUNE.shp  one representative point per commune (INSEE_COMM)
#       Download: https://data.cquest.org/ign/route500/ROUTE500_3-0__SHP_LAMB93_FXX_2021-11-03.7z
#   (B) Commune population 2021  -- INSEE (open data):
#         donnees_communes.csv   columns COM (5-digit INSEE code) and PTOT (total population)
#       Download : https://www.insee.fr/fr/statistiques/7739582?sommaire=7728826
#   (C) Healthcare facilities -- INSEE "Base de Permanences Et de structures (BPE)", 2025 cube:
#         DS_BPE_2025_data.csv   columns GEO, GEO_OBJECT, FACILITY_TYPE, BPE_MEASURE, VALUE
#           doctor    = FACILITY_TYPE "D265" (médecin généraliste libéral), GEO_OBJECT "COM"
#           maternity = FACILITY_TYPE "D107" (maternité)
#       Download : https://api.insee.fr/melodi/file/DS_BPE/DS_BPE_2025_CSV_FR
#   (D) Home-work commuting flows -- INSEE "Flux domicile-lieu de travail", 2023 census:
#         base-flux-mobilite-domicile-lieu-travail-2023.csv
#           columns CODGEO, LIBGEO (commune codes, 5 digits), DCLT, L_DCLT, flux
#       Download : https://www.insee.fr/fr/statistiques/fichier/8998300/base-flux-mobilite-domicile-lieu-travail-2023_csv.zip
# USAGE:
#   Rscript data_preparation.R [DATA_DIR]
#     DATA_DIR  directory holding the generated files (default: current directory).
#     Raw-file paths can be overridden with the environment variables
#     RAW_TRONCON, RAW_NOEUD, RAW_POP, RAW_BPE, RAW_COMMUTE (defaults below).
#
# Requires: packages sf, data.table, cppRouting (for the network extraction).
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(data.table)
  library(cppRouting)
})

# --------------------------------------------------------------------------- #
# 0. Configuration                                                            #
# --------------------------------------------------------------------------- #
args     <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- if (length(args) >= 1) args[1] else getwd()
setwd(DATA_DIR)
options(scipen = 999)
options(warn = 1)

base_route500 <- "/home/vincent/cppRouting/ROUTE500_3-0__SHP_LAMB93_FXX_2021-11-03/ROUTE500/1_DONNEES_LIVRAISON_2022-01-00175/R500_3-0_SHP_LAMB93_FXX-ED211/RESEAU_ROUTIER"

RAW_TRONCON <- Sys.getenv("RAW_TRONCON", file.path(base_route500, "TRONCON_ROUTE.shp"))
RAW_NOEUD   <- Sys.getenv("RAW_NOEUD",   file.path(base_route500, "NOEUD_COMMUNE.shp"))
RAW_POP     <- Sys.getenv("RAW_POP",     "/home/vincent/cppRouting/donnees_communes.csv")
RAW_BPE     <- Sys.getenv("RAW_BPE",     "/home/vincent/cppRouting/DS_BPE_2025_data.csv")
RAW_COMMUTE <- Sys.getenv("RAW_COMMUTE", "/home/vincent/cppRouting/base-flux-mobilite-domicile-lieu-travail-2023.csv")

# Speed (km/h) per road type, applied to derive edge travel time (minutes).
speed_map <- function(vocation, acces) {
  v <- vocation; a <- acces
  v[!is.na(a)] <- v[!is.na(a)]           # keep
  sp <- rep(NA_real_, length(v))
  sp[v == "Type autoroutier" & a == "A péage"]  <- 110
  sp[v == "Type autoroutier" & a == "Libre"]    <- 90
  sp[v == "Liaison principale"]                 <- 70
  sp[v == "Liaison régionale"]                  <- 60
  sp[v == "Liaison locale"]                     <- 50
  sp[v == "Bretelle"]                           <- 50
  sp
}

# --------------------------------------------------------------------------- #
# 1. Road network: roads.csv, coordinates.csv, roads_aux.csv                 #
# --------------------------------------------------------------------------- #
cat("== 1. Road network (TRONCON_ROUTE.shp -> cpp_extract) ==\n")
troncon <- st_read(RAW_TRONCON, quiet = TRUE)
cat(sprintf("   read %d road segments\n", nrow(troncon)))

# .dbf attributes are latin-1; normalize to UTF-8 so value matching works
troncon$VOCATION <- iconv(as.character(troncon$VOCATION), "latin1", "UTF-8")
troncon$ACCES    <- iconv(as.character(troncon$ACCES),    "latin1", "UTF-8")
troncon$SENS     <- iconv(as.character(troncon$SENS),     "latin1", "UTF-8")

troncon$vitesse <- speed_map(troncon$VOCATION, troncon$ACCES)
stopifnot(!any(is.na(troncon$vitesse)))

tmp <- file.path(tempdir(), "roads_tmp")   # out_prefix; files written as <tmp>_edges.csv etc.

cppRouting::cpp_extract(
  troncon, output = "csv", out_prefix = tmp, coords = TRUE, keep_length = TRUE,
  geom_params = list(
    speed_field     = "vitesse",
    length_field    = "LONGUEUR",
    direction_field = "SENS",
    direction_values = list(forward = "Sens direct",
                             backward = "Sens inverse",
                             twoway = "Double sens"),
    keep_fields = "LONGUEUR"
  )
)

roads <- fread(paste0(tmp, "_edges.csv"))     # from, to, dist
setnames(roads, "dist", "weight")
fwrite(roads, "roads.csv")

coords <- fread(paste0(tmp, "_coords.csv"))    # node_id, X, Y
setnames(coords, "node_id", "ID")
fwrite(coords, "coordinates.csv")

aux <- fread(paste0(tmp, "_aux.csv"))          # LONGUEUR, length
fwrite(aux, "roads_aux.csv")

unlink(paste0(tmp, c("_edges.csv", "_coords.csv", "_aux.csv")))
cat(sprintf("   roads.csv: %d edges | coordinates.csv: %d nodes | roads_aux.csv: %d rows\n",
            nrow(roads), nrow(coords), nrow(aux)))

# --------------------------------------------------------------------------- #
# 2. node_commune.csv  (commune -> representative node + population 2021)    #
# --------------------------------------------------------------------------- #
cat("== 2. node_commune.csv ==\n")
nc <- st_read(RAW_NOEUD, quiet = TRUE)
xy <- st_coordinates(nc)
nc_dt <- data.table(com = as.character(nc$INSEE_COMM), x = xy[, 1], y = xy[, 2])
nc_dt <- nc_dt[!duplicated(com)]                       # one row per commune
cat(sprintf("   %d unique communes\n", nrow(nc_dt)))

pts_nc <- st_as_sf(nc_dt, coords = c("x", "y"), crs = 2154)
pts_co <- st_as_sf(coords, coords = c("X", "Y"), crs = 2154)

raw <- st_nearest_feature(pts_nc, pts_co)              # index of nearest node
id_noeud <- coords$ID[raw]                             # (verified against existing mapping)
d_snap <- sqrt((nc_dt$x - coords$X[raw])^2 + (nc_dt$y - coords$Y[raw])^2)
stopifnot(max(d_snap) < 5e3)                           # sanity: every commune is close to its node
cat(sprintf("   snapped, max distance = %.0f m\n", max(d_snap)))

pop <- fread(RAW_POP)
pop <- pop[, .(COM, PTOT)]
setkey(pop, COM)

node_commune <- merge(nc_dt[, .(com, x, y)], pop, by.x = "com", by.y = "COM", all.x = TRUE)
node_commune[, POPULATION := ifelse(is.na(PTOT), 0L, as.integer(PTOT))]
node_commune <- node_commune[, .(com, id_noeud, POPULATION)]
setcolorder(node_commune, c("com", "id_noeud", "POPULATION"))
fwrite(node_commune, "node_commune.csv")
cat(sprintf("   wrote node_commune.csv (%d rows, total population %d)\n",
            nrow(node_commune), sum(node_commune$POPULATION)))

# --------------------------------------------------------------------------- #
# 3. doctor.csv  (general practitioners, 2025 BPE)                           #
# --------------------------------------------------------------------------- #
cat("== 3. doctor.csv ==\n")
bpe <- fread(RAW_BPE, sep = ";")
doctors <- bpe[GEO_OBJECT == "COM" & FACILITY_TYPE == "D265" & BPE_MEASURE == "FACILITIES",
               .(CODGEO = GEO, NB_D201 = as.integer(OBS_VALUE))]
doctors <- doctors[!duplicated(CODGEO)]
# Keep only mainland communes that have a node in the road network (drop overseas
# and the few mainland communes absent from the network); place GPs at the commune's node.
doctors <- merge(doctors, node_commune[, .(com, id_noeud, POPULATION)],
                 by.x = "CODGEO", by.y = "com", all.x = FALSE)
doctors <- doctors[, .(CODGEO, NB_D201, id_noeud, POPULATION)]
setcolorder(doctors, c("CODGEO", "NB_D201", "id_noeud", "POPULATION"))
fwrite(doctors, "doctor.csv")
cat(sprintf("   wrote doctor.csv (%d communes with GPs, %d GPs total)\n",
            nrow(doctors), sum(doctors$NB_D201)))

# --------------------------------------------------------------------------- #
# 4. maternity.csv  (maternity units, 2025 BPE)                              #
# --------------------------------------------------------------------------- #
cat("== 4. maternity.csv ==\n")
mat <- bpe[GEO_OBJECT == "COM" & FACILITY_TYPE == "D107" & BPE_MEASURE == "FACILITIES",
           .(CODGEO = GEO, NB_D107 = as.integer(OBS_VALUE))]
mat <- mat[!duplicated(CODGEO)]
mat <- mat[CODGEO %in% node_commune$com]   # mainland only (consistent with the network)
fwrite(mat[, .(CODGEO, NB_D107)], "maternity.csv")
cat(sprintf("   wrote maternity.csv (%d communes, %d units)\n", nrow(mat), sum(mat$NB_D107)))

# --------------------------------------------------------------------------- #
# 5. commuting.Rds  (home-work flows, 2023)                                  #
# --------------------------------------------------------------------------- #
cat("== 5. commuting.Rds ==\n")
cmt_raw <- fread(RAW_COMMUTE, sep = ";")
setnames(cmt_raw, "NBFLUX_C23_ACTOCC15P", "flux")   # flux value column
# Keep only flows whose origin is a known French mainland commune (node_commune)
cmt_raw <- cmt_raw[CODGEO %in% node_commune$com]
cmt <- cmt_raw[, .(CODGEO, LIBGEO, DCLT, L_DCLT, flux)]
save(cmt, file = "commuting.Rds", compress = "gzip")
cat(sprintf("   wrote commuting.Rds (%d flows)\n", nrow(cmt)))

cat("\nAll files written to", DATA_DIR, "\n")
