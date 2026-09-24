library(tidyverse)
library(janitor)
library(sf)
library(stars)
library(starsExtra)

# Util functions ----

#' Find the area of intersection between features of a sfc
#' 
#' @param data sf object of input data
#' @param covar sf object of covariable data
#' @param covar_column name of the column that will be use as covariable
#' @param unit output unit for the area intersection: m2, km2, ha
#' @param column_prefix prefix to use in the intersection columns
#'
#' @return sf object with additional columns for intersection areas
#'
#' @example
#' crop_terrain <- sf::st_read('crops_terrains.shp')
#' maize_and_soy_crops <- sf::st_read('maize_soy_crops.shp')
#' crop_area_for_maize_and_soy <- add_area_covar(crop_terrain, maize_and_soy_crops, 'crop', unit='km2', column_prefix='area')
add_area_covar <- function(data, covar, covar_column, unit = "m2", column_prefix = "area") {
  
  # Ensure both datasets have the same CRS
  if (!sf::st_crs(data) == sf::st_crs(covar)) {
    covar <- sf::st_transform(covar, sf::st_crs(data))
  }
  
  # Define conversion factor based on unit
  conversion_factor <- switch(unit,
                              "m2" = 1,
                              "km2" = 1e-6,
                              "ha" = 1e-4,
                              stop("Unit must be one of: m2, km2, ha"))
  
  # Add row index to data for joining later
  data <- data %>%
    dplyr::mutate(.row_id = dplyr::row_number())
  
  # Perform intersection once for all features
  intersections <- sf::st_intersection(
    data %>% dplyr::select(.row_id),
    covar %>% dplyr::select(dplyr::all_of(covar_column))
  ) %>%
    # Calculate area for each intersection
    dplyr::mutate(
      .area = as.numeric(sf::st_area(geometry)) * conversion_factor
    ) %>%
    # Drop geometry for faster processing
    sf::st_drop_geometry()
  
  # Summarize areas by row_id and covar_column value
  area_summary <- intersections %>%
    dplyr::group_by(.row_id, !!rlang::sym(covar_column)) %>%
    dplyr::summarise(.area = sum(.area), .groups = "drop")
  
  # Pivot to wide format with one column per unique covar value
  area_wide <- area_summary %>%
    tidyr::pivot_wider(
      names_from = !!rlang::sym(covar_column),
      values_from = .area,
      values_fill = 0,
      names_prefix = paste0(column_prefix, "_")
    )
  
  # Join back to original data and remove temporary row_id
  result <- data %>%
    dplyr::left_join(area_wide, by = ".row_id") %>%
    dplyr::select(-.row_id)
  
  # Replace NA with 0 for features with no intersections
  area_columns <- names(result)[grepl(paste0("^", column_prefix, "_"), names(result))]
  result <- result %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(area_columns), ~tidyr::replace_na(., 0))) %>% 
    janitor::clean_names()
  
  return(result)
}

#' Add stats from raster
#' 
#' @param data sf object of input data
#' @param raster_covar stars object covariable data
#' @param funs function or list of functions to apply to pixels intersected by feature
#' @param bands vector of band indices or names to process. NULL processes all bands
#' @param column_prefix prefix to use in the output columns
#' @param na.rm logical, should NA values be removed before applying functions
#'
#' @return sf object with additional columns with the results of applying functions to raster data by feature
#'
#' @examples
#' crop_terrain <- sf::st_read('crops_terrains.shp')
#' gradient <- stars::read_stars('elevation_gradient.tif')
#' 
#' # Single function, single band
#' crop_mean_gradient <- add_stat_raster(crop_terrain, gradient, mean)
#' 
#' # Multiple functions, single band
#' crop_stats <- add_stat_raster(crop_terrain, gradient, 
#'                                funs = list(mean = mean, sd = sd, max = max))
#' 
#' # Multiple functions, multiple bands
#' multi_band_raster <- stars::read_stars('multi_band.tif')
#' crop_multi_stats <- add_stat_raster(crop_terrain, multi_band_raster,
#'                                      funs = list(mean = mean, median = median),
#'                                      bands = c(1, 3))
add_stat_raster <- function(data, 
                            raster_covar, 
                            funs = mean, 
                            bands = NULL,
                            column_prefix = "raster",
                            na.rm = TRUE) {
  
  # Ensure CRS compatibility
  if (!sf::st_crs(data) == sf::st_crs(raster_covar)) {
    raster_covar <- stars::st_warp(raster_covar, crs = sf::st_crs(data))
  }
  
  # Convert single function to named list
  if (!is.list(funs)) {
    fun_name <- deparse(substitute(funs))
    funs <- setNames(list(funs), fun_name)
  }
  
  # Ensure functions are named
  if (is.null(names(funs))) {
    names(funs) <- paste0("fun", seq_along(funs))
  }
  
  # Select bands if specified
  dims <- names(dim(raster_covar))
  band_dim <- dims[!(dims %in% c("x", "y"))]
  
  if (!is.null(bands) && length(band_dim) > 0) {
    raster_covar <- raster_covar %>%
      dplyr::slice(!!rlang::sym(band_dim[1]), bands)
  }
  
  # Add row identifier
  data <- data %>%
    dplyr::mutate(.row_id = dplyr::row_number())
  
  # Use aggregate to extract values for each polygon
  # This is more reliable than st_extract for getting actual pixel values
  result_list <- list()
  
  for (fn_name in names(funs)) {
    fn <- funs[[fn_name]]
    
    # Create wrapper function that handles na.rm
    wrapper_fn <- function(x) {
      if (na.rm) {
        x <- x[!is.na(x)]
      }
      if (length(x) == 0) return(NA_real_)
      fn(x)
    }
    
    # Aggregate raster by polygons
    agg_result <- aggregate(
      raster_covar,
      by = data,
      FUN = wrapper_fn
    )
    
    # Extract values
    col_name <- paste0(column_prefix, "_", fn_name)
    result_list[[col_name]] <- agg_result[[1]]
  }
  
  # Combine all statistics
  result_stats <- tibble::as_tibble(result_list) %>%
    janitor::clean_names()
  
  # Bind to original data
  result <- dplyr::bind_cols(
    data %>% dplyr::select(-.row_id),
    result_stats
  )
  
  return(result)
}

#' Transform and filter a spatial covariate
#'
#' MULTIPOLYGON inputs are split into POLYGON features before filtering so
#' individual components outside the area of interest are discarded.
#'
#' @param data sf object containing the covariate
#' @param aoi sf or sfc object defining the area of interest
#'
#' @return sf object transformed to the AOI CRS and filtered by the AOI
process_spatial_covar <- function(data,
                                  aoi) {
  geometry_types <- unique(as.character(sf::st_geometry_type(data)))

  if ("MULTIPOLYGON" %in% geometry_types) {
    data <- sf::st_cast(
      data,
      "POLYGON",
      group_or_split = FALSE
    )
  }

  data_aoi <- data %>%
    sf::st_transform(sf::st_crs(aoi)) %>%
    sf::st_filter(aoi)

  return(data_aoi)
}

#' Load and process the spatial covariates listed in a catalog
#'
#' Each unique shapefile is read and processed once, even when the catalog uses
#' more than one of its attribute columns.
#'
#' @param catalog validated spatial covariate catalog
#' @param covariates_dir directory containing the catalog shapefiles
#' @param aoi sf or sfc object defining the area of interest
#'
#' @return named list of processed sf objects keyed by shapefile
load_spatial_covariates <- function(catalog,
                                    covariates_dir,
                                    aoi) {
  covariate_shapefiles <- unique(catalog$shapefile)
  spatial_covariates <- stats::setNames(
    vector("list", length(covariate_shapefiles)),
    covariate_shapefiles
  )

  for (shapefile in covariate_shapefiles) {
    covariate <- sf::st_read(
      file.path(covariates_dir, shapefile),
      quiet = TRUE
    )

    spatial_covariates[[shapefile]] <- process_spatial_covar(
      covariate,
      aoi
    )
  }

  return(spatial_covariates)
}


# Load and select geoms ----
INPUT_DIR <- "00-raw_data/02_Entrega_05082026/Shapefiles"
OUTPUT_DIR <- '02-processed_data'
current_date <- now() %>% date()
shp_output <- fs::path_join(c(OUTPUT_DIR, "aoi_shp"))
tiff_output <- fs::path_join(c(OUTPUT_DIR, "aoi_tiff"))

fs::dir_create(shp_output)
fs::dir_create(tiff_output)

nanocuencas_ids <- c("ESTADO_SHT", "CLAVE_SHT")
nanocuecas_data <- sf::read_sf(file.path(INPUT_DIR,'Nanocuencas/SHT_V3_uw.shp'))
aoi_bbox <- nanocuecas_data %>% st_bbox() %>% st_as_sfc() %>% st_buffer(3000) %>% st_bbox() %>% st_as_sfc()

# Read, transform, and filter spatial covariates from the catalog
covariates_dir <- file.path(INPUT_DIR, "Covariables")
covariate_catalog <- readr::read_csv(file.path(INPUT_DIR, "catalogo.csv"))

spatial_covariates <- load_spatial_covariates(
  covariate_catalog,
  covariates_dir,
  aoi_bbox
)

# Add covar data raster

world_dem <- read_stars('00-raw_data/DEM/wc2.1_30s_elev.tif')
dem_crs <- world_dem %>% st_crs()

# Crop to AOI
aoi_bbox_t <- aoi_bbox %>% 
  st_transform(dem_crs)
aoi_dem <- world_dem %>% 
  st_crop(aoi_bbox_t) %>% 
  st_as_stars()

# Transform to projected data and slope calculation
aoi_bbox_grid <- aoi_bbox %>% 
  st_as_stars()
aoi_dem_proj <- aoi_dem %>% 
  st_warp(aoi_bbox_grid) 
aoi_dem_slope <- aoi_dem_proj %>% 
  starsExtra::slope()

# Save raster data
aoi_dem %>% 
  write_stars(fs::path_join(c(tiff_output, paste(current_date, 'dem_aoi.tiff', sep="_"))))

aoi_dem_slope %>% 
  write_stars(fs::path_join(c(tiff_output, paste(current_date, 'dem_slope_aoi.tiff', sep="_"))))

# Create dataset ----
# Add covar data - shapefile
dataset <- nanocuecas_data %>%
  select(nanocuencas_ids)

for (catalog_row in seq_len(nrow(covariate_catalog))) {
  catalog_entry <- covariate_catalog[catalog_row, ]
  dataset <- add_area_covar(
    dataset,
    spatial_covariates[[catalog_entry$shapefile]],
    catalog_entry$variable,
    unit = "ha",
    column_prefix = catalog_entry$prefix
  )
}

# Add covar data - raster 
fn_list <- list(
  "min" = min,
  "max" = max,
  "median" = median,
  "mean" = mean,
  "sd" = sd
)

dataset <- dataset %>% 
  add_stat_raster(aoi_dem_slope, funs = fn_list, column_prefix = "slope", na.rm = TRUE)

# Save dataset shp
dataset %>% st_write(fs::path_join(c(shp_output, paste(current_date, "dataset.gpkg", sep="_"))))

dataset %>% 
  st_drop_geometry() %>% 
  mutate(across(starts_with("slope"), as.numeric)) %>% 
  write_csv(fs::path_join(c(OUTPUT_DIR, paste(current_date, "dataset.csv", sep="_"))))
