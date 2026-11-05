# Tuiles AlphaEarth simulees : COG locaux int8 (NoData -128) en UTM 32N et un
# index sf au format publie. Aucun reseau.

.vecteurs_unitaires <- function(n, d, seed = 1) {
  set.seed(seed)
  m <- matrix(stats::rnorm(n * d), n, d)
  m / sqrt(rowSums(m^2))
}

# Quantification inverse de la de-quantification du README.
.quantifier <- function(v) {
  q <- round(sign(v) * sqrt(abs(v)) * 127.5)
  q[q > 127] <- 127
  q[q < -127] <- -127
  q
}

# Une tuile UTM 32N de 40 x 40 pixels de 10 m, d bandes, NoData en (1, 1).
.tuile <- function(dossier, nom, xmin, ymin, d = 4, seed = 1) {
  r <- terra::rast(xmin = xmin, xmax = xmin + 400, ymin = ymin, ymax = ymin + 400,
                   resolution = 10, nlyrs = d, crs = "EPSG:32632")
  v <- .quantifier(.vecteurs_unitaires(terra::ncell(r), d, seed))
  v[1, ] <- -128
  terra::values(r) <- v
  names(r) <- sprintf("A%02d", seq_len(d) - 1L)
  f <- file.path(dossier, nom)
  terra::writeRaster(r, f, datatype = "INT1S", NAflag = -128, overwrite = TRUE)
  f
}

.index <- function(fichiers, annees) {
  geoms <- lapply(fichiers, function(f) {
    sf::st_transform(sf::st_as_sfc(sf::st_bbox(terra::rast(f))), 4326)[[1]]
  })
  # `crs` est aussi un argument de st_sf() : on passe par un data.frame.
  df <- data.frame(path = fichiers, crs = "EPSG:32632", year = annees,
                   utm_zone = "32N")
  sf::st_sf(df, geometry = sf::st_sfc(geoms, crs = 4326))
}

# Emprise dans la premiere tuile, en Lambert-93.
.aoi_dans <- function(x0, y0, cote = 200) {
  sf::st_transform(sf::st_as_sfc(sf::st_bbox(
    c(xmin = x0, ymin = y0, xmax = x0 + cote, ymax = y0 + cote),
    crs = sf::st_crs(32632)
  )), 2154)
}

test_that("aef_dequantize applies the README mapping and drops -128", {
  expect_equal(aef_dequantize(c(-127, 0, 64, 127.5)),
               c(-(127 / 127.5)^2, 0, (64 / 127.5)^2, 1))
  expect_true(is.na(aef_dequantize(-128)))
  # un vecteur quantifie puis de-quantifie garde une norme ~ 1
  v <- .vecteurs_unitaires(50, 64)
  dq <- aef_dequantize(.quantifier(v))
  expect_true(all(abs(sqrt(rowSums(dq^2)) - 1) < 0.03))
})

test_that("aef_dequantize works on a SpatRaster", {
  r <- terra::rast(nrows = 1, ncols = 3, vals = c(-128, -64, 127))
  expect_equal(as.vector(terra::values(aef_dequantize(r))),
               c(NA, -(64 / 127.5)^2, (127 / 127.5)^2))
})

test_that(".aef_url maps the bucket to HTTP and leaves local paths alone", {
  expect_equal(
    alphaearthnemeton:::.aef_url("gs://alphaearth_foundations/satellite_embedding/v1/annual/2024/32N/x-0-0.tiff"),
    "/vsicurl/https://storage.googleapis.com/alphaearth_foundations/satellite_embedding/v1/annual/2024/32N/x-0-0.tiff"
  )
  withr::with_options(list(alphaearthnemeton.base_url = "https://mirror.example/aef/"), {
    expect_equal(alphaearthnemeton:::.aef_url("gs://alphaearth_foundations/a/b.tiff"),
                 "/vsicurl/https://mirror.example/aef/a/b.tiff")
  })
  expect_equal(alphaearthnemeton:::.aef_url("/tmp/x.tif"), "/tmp/x.tif")
})

test_that("aef_info carries the licence and the attribution", {
  i <- aef_info()
  expect_equal(i$license, "CC-BY-4.0")
  expect_match(i$attribution, "Google DeepMind")
  expect_length(i$bands, 64)
  expect_match(i$index, "aef_index\\.gpkg$")
})

test_that("aef_index filters an in-memory index by area and year", {
  withr::with_tempdir({
    f1 <- .tuile(".", "a.tif", 400000, 5400000)
    f2 <- .tuile(".", "b.tif", 400000, 5400000, seed = 2)
    idx <- .index(c(f1, f2), c(2023L, 2024L))
    aoi <- .aoi_dans(400100, 5400100)
    expect_equal(nrow(aef_index(aoi, index = idx)), 2L)
    expect_equal(aef_index(aoi, years = 2024, index = idx)$year, 2024L)
    loin <- .aoi_dans(500000, 5400000)
    expect_equal(nrow(aef_index(loin, index = idx)), 0L)
  })
})

test_that("aef_fetch writes a de-quantized raster on a 10 m grid and caches it", {
  withr::with_tempdir({
    f <- .tuile(".", "a.tif", 400000, 5400000)
    idx <- .index(f, 2024L)
    aoi <- .aoi_dans(400100, 5400100)
    r <- suppressMessages(aef_fetch(aoi, index = idx, cache_dir = "cache"))

    expect_s4_class(r, "SpatRaster")
    expect_equal(terra::nlyr(r), 4L)
    expect_equal(names(r), sprintf("A%02d", 0:3))
    expect_equal(terra::res(r), c(10, 10))
    expect_equal(sf::st_crs(terra::crs(r))$epsg, 2154L)
    e <- as.vector(terra::ext(r))
    expect_true(all(e %% 10 == 0))
    expect_equal(attr(r, "embedding_source"), "aef")
    expect_equal(attr(r, "annee"), 2024L)

    v <- terra::values(r)
    v <- v[stats::complete.cases(v), , drop = FALSE]
    expect_gt(nrow(v), 100)
    expect_true(all(abs(v) <= 1))
    expect_true(all(abs(sqrt(rowSums(v^2)) - 1) < 0.05))

    fichiers <- list.files("cache", pattern = "^aef_2024_.*\\.tif$")
    expect_length(fichiers, 1L)
    avant <- file.mtime(file.path("cache", fichiers))
    Sys.sleep(1.1)
    r2 <- suppressMessages(aef_fetch(aoi, index = idx, cache_dir = "cache"))
    expect_equal(file.mtime(file.path("cache", fichiers)), avant)
    expect_equal(terra::values(r2), terra::values(r))
  })
})

test_that("aef_fetch returns a list named by year, latest year by default", {
  withr::with_tempdir({
    f1 <- .tuile(".", "a.tif", 400000, 5400000, seed = 1)
    f2 <- .tuile(".", "b.tif", 400000, 5400000, seed = 2)
    idx <- .index(c(f1, f2), c(2023L, 2024L))
    aoi <- .aoi_dans(400100, 5400100)
    seul <- suppressMessages(aef_fetch(aoi, index = idx, cache_dir = "c"))
    expect_equal(attr(seul, "annee"), 2024L)
    deux <- suppressMessages(aef_fetch(aoi, years = 2023:2024, index = idx, cache_dir = "c"))
    expect_named(deux, c("2023", "2024"))
    expect_false(isTRUE(all.equal(terra::values(deux[["2023"]]),
                                  terra::values(deux[["2024"]]))))
    expect_error(suppressMessages(aef_fetch(aoi, years = 2019, index = idx, cache_dir = "c")),
                 "No AlphaEarth tile")
  })
})

test_that("aef_fetch mosaics two adjacent tiles", {
  withr::with_tempdir({
    f1 <- .tuile(".", "a.tif", 400000, 5400000, seed = 1)
    f2 <- .tuile(".", "b.tif", 400400, 5400000, seed = 2)
    idx <- .index(c(f1, f2), c(2024L, 2024L))
    # a cheval sur la limite x = 400400
    aoi <- .aoi_dans(400300, 5400100)
    r <- suppressMessages(aef_fetch(aoi, index = idx, cache_dir = "c"))
    v <- terra::values(r)
    # pas de trou sur la jointure : toutes les cellules du coeur de l'emprise sont pleines
    coeur <- terra::crop(r, terra::ext(r) - 30)
    expect_false(anyNA(terra::values(coeur)))
  })
})

test_that("aef_fetch refuses an area no tile covers", {
  withr::with_tempdir({
    f <- .tuile(".", "a.tif", 400000, 5400000)
    idx <- .index(f, 2024L)
    expect_error(suppressMessages(aef_fetch(.aoi_dans(600000, 5400000), index = idx,
                                            cache_dir = "c")), "No AlphaEarth tile")
  })
})

test_that("aef_fetch reads the published COGs (network, opt-in)", {
  skip_if_not(identical(Sys.getenv("AEF_TEST_NETWORK"), "true"),
              "set AEF_TEST_NETWORK=true to query the public bucket")
  aoi <- sf::st_transform(sf::st_as_sfc(sf::st_bbox(
    c(xmin = 7.235, ymin = 48.645, xmax = 7.245, ymax = 48.655), crs = 4326)), 2154)
  r <- aef_fetch(aoi, years = 2024, cache_dir = withr::local_tempdir())
  expect_equal(terra::nlyr(r), 64L)
  v <- terra::values(r)
  v <- v[stats::complete.cases(v), , drop = FALSE]
  expect_true(all(abs(sqrt(rowSums(v^2)) - 1) < 0.02))
})
