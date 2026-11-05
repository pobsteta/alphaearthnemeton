#' AlphaEarth Foundations satellite embeddings for nemeton
#'
#' Pure-R access to the annual AlphaEarth Foundations Satellite Embedding
#' (Google DeepMind), published as Cloud Optimized GeoTIFFs in the public
#' `gs://alphaearth_foundations` bucket. [aef_fetch()] finds the tiles that
#' cover an area of interest through the dataset's GeoPackage index, reads
#' only the needed window of each COG over HTTP (GDAL `/vsicurl/`),
#' de-quantizes the signed 8-bit values, reprojects and writes one
#' content-addressed COG per year. The result feeds
#' `nemeton::indicateur_n5_distance_reference()` (spec 064).
#'
#' @keywords internal
"_PACKAGE"

# Racine HTTP du jeu (copie Google ; le miroir Source Cooperative porte la
# meme arborescence). Surchargeable : options(alphaearthnemeton.base_url = ).
.AEF_GS_PREFIX <- "gs://alphaearth_foundations/"
.AEF_HTTP_DEFAUT <- "https://storage.googleapis.com/alphaearth_foundations/"
.AEF_INDEX_RELATIF <- "satellite_embedding/v1/annual/aef_index.gpkg"
.AEF_VERSION <- "v1"

.aef_base_url <- function() {
  getOption("alphaearthnemeton.base_url", .AEF_HTTP_DEFAUT)
}

#' Dataset information
#'
#' @return A list: `version`, `years` (published years), `bucket`,
#'   `base_url` (HTTP root actually used), `index` (GeoPackage index URL),
#'   `resolution_m`, `bands`, `license` and the mandatory `attribution` text.
#' @export
#' @examples
#' aef_info()$attribution
aef_info <- function() {
  list(
    version = .AEF_VERSION,
    years = 2017:2025,
    bucket = "gs://alphaearth_foundations/satellite_embedding/v1/annual/",
    base_url = .aef_base_url(),
    index = paste0(.aef_base_url(), .AEF_INDEX_RELATIF),
    resolution_m = 10,
    bands = sprintf("A%02d", 0:63),
    license = "CC-BY-4.0",
    attribution = paste(
      "The AlphaEarth Foundations Satellite Embedding dataset is produced",
      "by Google and Google DeepMind."
    )
  )
}

#' De-quantize raw AlphaEarth values
#'
#' The COGs store each dimension as a signed 8-bit integer; -128 is the
#' no-data value. The analysis-ready value in `[-1, 1]` is
#' `sign(x) * (x / 127.5)^2` (dataset README). Each de-quantized pixel
#' vector has unit Euclidean norm.
#'
#' @param x Numeric vector, matrix or SpatRaster of raw values.
#' @return Same class as `x`, de-quantized, with -128 turned into `NA`.
#' @export
#' @examples
#' aef_dequantize(c(-127, 0, 64, 127, -128))
aef_dequantize <- function(x) {
  if (inherits(x, "SpatRaster")) {
    x <- terra::classify(x, cbind(-128, NA))
    return(sign(x) * (x / 127.5)^2)
  }
  x[!is.na(x) & x == -128] <- NA
  sign(x) * (x / 127.5)^2
}

#' Tiles of the AlphaEarth index covering an area
#'
#' Queries the dataset's GeoPackage index (one row per COG) with a spatial
#' filter, reading only the needed pages over HTTP.
#'
#' @param aoi Area of interest: sf, sfc, SpatVector or bbox, any CRS.
#' @param years Integer vector of years, or `NULL` for every year.
#' @param index Optional sf index already in memory (same columns as the
#'   published one: `path`, `crs`, `year`, `utm_zone`, geometry in
#'   EPSG:4326). `NULL` reads the published index.
#' @return An sf data frame with `path`, `url` (GDAL path to read),
#'   `crs`, `year`, `utm_zone`, one row per COG.
#' @export
aef_index <- function(aoi, years = NULL, index = NULL) {
  zone <- .aef_aoi_4326(aoi)
  if (is.null(index)) {
    src <- paste0("/vsicurl/", aef_info()$index)
    index <- tryCatch(
      sf::st_read(src, wkt_filter = sf::st_as_text(zone), quiet = TRUE),
      error = function(e) {
        cli::cli_abort(c(
          "Cannot read the AlphaEarth index.",
          "i" = "{.url {aef_info()$index}}",
          "x" = conditionMessage(e)
        ))
      }
    )
  } else {
    index <- sf::st_transform(index, 4326)
    index <- index[lengths(sf::st_intersects(index, zone)) > 0, ]
  }
  if (!is.null(years)) index <- index[index$year %in% years, ]
  index$url <- vapply(index$path, .aef_url, "")
  index
}

#' Fetch AlphaEarth embeddings over an area of interest
#'
#' For each year, reads the window of every COG that covers `aoi` (plus
#' `buffer`), de-quantizes it, reprojects it to `crs` on a 10 m grid snapped
#' to multiples of 10 m (nearest neighbour: an embedding is a vector, an
#' interpolation would break its unit norm), mosaics the tiles and writes a
#' float32 COG in `cache_dir`. The cache key hashes the dataset version, the
#' year, the target CRS, the snapped extent and the source tiles: a second
#' call over the same area costs nothing.
#'
#' @param aoi Area of interest: sf, sfc, SpatVector or bbox, any CRS.
#' @param years Integer vector of years. `NULL`: the most recent year
#'   available over `aoi`.
#' @param crs Target CRS (anything [sf::st_crs()] accepts). Default
#'   EPSG:2154 (Lambert-93); use 3035 outside France.
#' @param buffer Numeric, metres added around `aoi`. Default 0.
#' @param cache_dir Directory of the cached COGs. Default
#'   [aef_cache_dir()].
#' @param index Optional in-memory index, see [aef_index()].
#' @param overwrite Logical. Rebuild a cached file.
#'
#' @return One year: a 64-band SpatRaster (`A00` to `A63`). Several years:
#'   a list of SpatRaster named by year, the form
#'   `nemeton::indicateur_n5_distance_reference(embedding =)` accepts.
#'   Each raster's source file is a COG in `cache_dir`.
#'
#' @section Access:
#' The bucket is declared "requester pays" by Google; anonymous HTTP range
#' requests were served on 2026-11-05. If Google starts enforcing it, set
#' `options(alphaearthnemeton.base_url =)` to a mirror (Source Cooperative:
#' `https://data.source.coop/tge-labs/aef/v1/annual/`, same layout under a
#' different root) or configure GDAL's GCS credentials.
#'
#' @export
aef_fetch <- function(aoi, years = NULL, crs = 2154, buffer = 0,
                      cache_dir = aef_cache_dir(), index = NULL,
                      overwrite = FALSE) {
  crs <- sf::st_crs(crs)
  if (is.na(crs)) cli::cli_abort("{.arg crs} is not a valid CRS.")
  tuiles <- aef_index(aoi, years = years, index = index)
  if (!nrow(tuiles)) {
    cli::cli_abort("No AlphaEarth tile covers {.arg aoi} for the years asked.")
  }
  if (is.null(years)) years <- max(tuiles$year)
  manquantes <- setdiff(years, tuiles$year)
  if (length(manquantes)) {
    cli::cli_abort("No AlphaEarth tile over {.arg aoi} for {.val {manquantes}}.")
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  zone <- sf::st_transform(.aef_aoi_sf(aoi), crs)
  bb <- sf::st_bbox(zone)
  pas <- 10
  ext_cible <- c(
    xmin = floor((bb[["xmin"]] - buffer) / pas) * pas,
    xmax = ceiling((bb[["xmax"]] + buffer) / pas) * pas,
    ymin = floor((bb[["ymin"]] - buffer) / pas) * pas,
    ymax = ceiling((bb[["ymax"]] + buffer) / pas) * pas
  )

  res <- lapply(years, function(an) {
    t_an <- tuiles[tuiles$year == an, ]
    cle <- .aef_hash(c(.AEF_VERSION, an, crs$wkt, format(ext_cible, nsmall = 0),
                       sort(t_an$path)))
    fichier <- file.path(cache_dir, sprintf("aef_%s_%s.tif", an, substr(cle, 1, 16)))
    if (file.exists(fichier) && !overwrite) {
      cli::cli_alert_success("AlphaEarth {an}: cache hit ({.file {basename(fichier)}})")
    } else {
      cli::cli_alert_info("AlphaEarth {an}: {nrow(t_an)} tile{?s} to read")
      out <- .aef_assemble(t_an, crs, ext_cible, pas)
      .aef_write_cog(out, fichier)
      cli::cli_alert_success("AlphaEarth {an}: written {.file {basename(fichier)}}")
    }
    r <- terra::rast(fichier)
    names(r) <- sprintf("A%02d", seq_len(terra::nlyr(r)) - 1L)
    attr(r, "embedding_source") <- "aef"
    attr(r, "annee") <- as.integer(an)
    r
  })
  names(res) <- as.character(years)
  if (length(res) == 1L) res[[1L]] else res
}

#' Default cache directory
#'
#' `getOption("alphaearthnemeton.cache_dir")`, else the `NEMETON_CACHE_DIR`
#' environment variable + `/aef`, else the user cache directory.
#' @export
aef_cache_dir <- function() {
  opt <- getOption("alphaearthnemeton.cache_dir")
  if (!is.null(opt)) return(opt)
  env <- Sys.getenv("NEMETON_CACHE_DIR")
  if (nzchar(env)) return(file.path(env, "aef"))
  file.path(tools::R_user_dir("alphaearthnemeton", "cache"))
}

# --------------------------------------------------------------- internes ---

# gs://alphaearth_foundations/x -> /vsicurl/<base_url>x ; chemin local inchange.
.aef_url <- function(path) {
  if (startsWith(path, .AEF_GS_PREFIX)) {
    return(paste0("/vsicurl/", .aef_base_url(),
                  substring(path, nchar(.AEF_GS_PREFIX) + 1L)))
  }
  if (grepl("^https?://", path)) return(paste0("/vsicurl/", path))
  path
}

.aef_aoi_sf <- function(aoi) {
  if (inherits(aoi, "SpatVector")) aoi <- sf::st_as_sf(aoi)
  if (inherits(aoi, "bbox")) aoi <- sf::st_as_sfc(aoi)
  if (inherits(aoi, "sf")) aoi <- sf::st_geometry(aoi)
  if (!inherits(aoi, "sfc")) {
    cli::cli_abort("{.arg aoi} must be an sf, sfc, SpatVector or bbox.")
  }
  if (is.na(sf::st_crs(aoi))) cli::cli_abort("{.arg aoi} has no CRS.")
  sf::st_union(sf::st_make_valid(aoi))
}

.aef_aoi_4326 <- function(aoi) sf::st_transform(.aef_aoi_sf(aoi), 4326)

.aef_hash <- function(x) {
  f <- tempfile()
  on.exit(unlink(f), add = TRUE)
  writeLines(as.character(x), f, useBytes = TRUE)
  unname(tools::md5sum(f))
}

# Lit la fenetre utile de chaque tuile, de-quantifie, reprojette sur la
# grille cible (plus proche voisin) et mosaique (premiere valeur non NA).
.aef_assemble <- function(tuiles, crs, ext_cible, pas) {
  gabarit <- terra::rast(
    xmin = ext_cible[["xmin"]], xmax = ext_cible[["xmax"]],
    ymin = ext_cible[["ymin"]], ymax = ext_cible[["ymax"]],
    resolution = pas, crs = crs$wkt
  )
  poly_cible <- terra::as.polygons(terra::ext(gabarit), crs = terra::crs(gabarit))
  morceaux <- list()
  for (i in seq_len(nrow(tuiles))) {
    src <- terra::rast(tuiles$url[i])
    zone_src <- terra::project(poly_cible, terra::crs(src))
    e <- terra::intersect(terra::ext(zone_src) + 2 * terra::res(src)[1],
                          terra::ext(src))
    if (is.null(e)) next
    brut <- terra::crop(src, e)
    morceaux[[length(morceaux) + 1L]] <-
      terra::project(aef_dequantize(brut), gabarit, method = "near")
  }
  if (!length(morceaux)) {
    cli::cli_abort("The AlphaEarth tiles listed by the index do not overlap {.arg aoi}.")
  }
  out <- Reduce(terra::cover, morceaux)
  names(out) <- sprintf("A%02d", seq_len(terra::nlyr(out)) - 1L)
  out
}

.aef_write_cog <- function(r, fichier) {
  tmp <- paste0(fichier, ".part")
  on.exit(unlink(tmp), add = TRUE)
  ok <- tryCatch({
    # Le pilote COG passe par un raster MEM qui ne connait pas COMPRESS :
    # avertissement GDAL sans consequence, la compression s'applique au COG.
    withCallingHandlers(
      terra::writeRaster(r, tmp, filetype = "COG", datatype = "FLT4S",
                         gdal = c("COMPRESS=DEFLATE"), overwrite = TRUE),
      warning = function(w) {
        if (grepl("driver MEM does not support creation option",
                  conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    # GDAL sans pilote COG : GeoTIFF tuile, meme contenu.
    terra::writeRaster(r, tmp, filetype = "GTiff", datatype = "FLT4S",
                       gdal = c("COMPRESS=DEFLATE", "TILED=YES"), overwrite = TRUE)
  }
  file.rename(tmp, fichier)
  invisible(fichier)
}
