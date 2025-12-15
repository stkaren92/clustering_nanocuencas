library(tidyverse)
library(janitor)
library(sf)

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
    raster_covar <- sf::st_transform(raster_covar, sf::st_crs(data))
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
  
  # Get dimension names and handle bands
  dims <- names(dim(raster_covar))
  band_dim <- dims[!(dims %in% c("x", "y"))]
  
  # Select bands if specified
  if (!is.null(bands) && length(band_dim) > 0) {
    raster_covar <- raster_covar %>%
      dplyr::slice(!!rlang::sym(band_dim[1]), bands)
  }
  
  # Get band names/indices
  if (length(band_dim) > 0) {
    band_values <- stars::st_get_dimension_values(raster_covar, band_dim[1])
    if (is.null(band_values)) {
      band_values <- seq_len(dim(raster_covar)[band_dim[1]])
    }
  } else {
    band_values <- 1
    band_dim <- NULL
  }
  
  # Add row identifier
  data <- data %>%
    dplyr::mutate(.row_id = dplyr::row_number())
  
  # Extract raster values for each feature
  extracted_values <- stars::st_extract(raster_covar, data)
  
  # Convert to tibble for easier manipulation
  extracted_df <- extracted_values %>%
    tibble::as_tibble() %>%
    dplyr::mutate(.row_id = dplyr::row_number())
  
  # Get the name of the raster attribute
  raster_attr <- names(extracted_df)[1]
  
  # Process based on whether we have bands or not
  if (!is.null(band_dim) && length(band_values) > 1) {
    # Multiple bands case
    result_stats <- extracted_df %>%
      dplyr::select(.row_id, dplyr::all_of(raster_attr)) %>%
      tidyr::unnest(cols = dplyr::all_of(raster_attr)) %>%
      dplyr::group_by(.row_id, !!rlang::sym(band_dim[1])) %>%
      dplyr::summarise(
        dplyr::across(
          .cols = 1,  # The raster values column
          .fns = funs,
          .names = "{.fn}",
          na.rm = na.rm
        ),
        .groups = "drop"
      ) %>%
      tidyr::pivot_longer(
        cols = -c(.row_id, !!rlang::sym(band_dim[1])),
        names_to = ".stat",
        values_to = ".value"
      ) %>%
      dplyr::mutate(
        .col_name = paste(column_prefix, !!rlang::sym(band_dim[1]), .stat, sep = "_")
      ) %>%
      dplyr::select(.row_id, .col_name, .value) %>%
      tidyr::pivot_wider(
        names_from = .col_name,
        values_from = .value
      )
  } else {
    # Single band case
    result_stats <- extracted_df %>%
      dplyr::select(.row_id, dplyr::all_of(raster_attr)) %>%
      tidyr::unnest(cols = dplyr::all_of(raster_attr)) %>%
      dplyr::group_by(.row_id) %>%
      dplyr::summarise(
        dplyr::across(
          .cols = 1,  # The raster values column
          .fns = funs,
          .names = paste0(column_prefix, "_{.fn}"),
          na.rm = na.rm
        ),
        .groups = "drop"
      )
  }
  
  # Join results back to original data
  result <- data %>%
    dplyr::left_join(result_stats, by = ".row_id") %>%
    dplyr::select(-.row_id)
  
  return(result)
}


# Dataset creation ----
nanocuecas_data <- sf::read_sf('00-raw_data/nanocuencas/SHT_BAconsensuadov2_uw.shp')

geologia_data <- sf::st_read('00-raw_data/geologia/r250k_ccl_dissolveTipo.shp')
suelos_data <- sf::st_read('00-raw_data/suelos/contedafo.shp')
usvVII_data <- sf::read_sf('00-raw_data/uso_suelo_vegetacion/usv250s7cw.shp')

st_crs(suelos_data) <- geologia_data %>% st_crs()

aoi_bbox <- nanocuecas_data %>% st_bbox() %>% st_as_sfc()

geologia_data_poly <- st_cast(geologia_data, "POLYGON", 
                              group_or_split = FALSE) # Convert mutipolygon to polygon

# Set all CRS to my aoi CRS 
aoi_crs <- aoi_bbox %>% st_crs()

geologia_data_poly <- geologia_data_poly %>% st_transform(aoi_crs)
suelos_data <- suelos_data %>% st_transform(aoi_crs)
usvVII_data <- usvVII_data %>% st_transform(aoi_crs)

# Filter by AOI
geologia_data_aoi <- geologia_data_poly %>% 
  st_filter(aoi_bbox)

suelos_data_aoi <- suelos_data %>% 
  st_filter(aoi_bbox)

usvVII_data_aoi <- usvVII_data %>% 
  st_filter(aoi_bbox)

# Select columns of interest
geologia_data_aoi <- geologia_data_aoi %>% 
  select(TIPO)

suelos_data_aoi <- suelos_data_aoi %>% 
  select(CLAVE_IMP)

usvVII_data_aoi <- usvVII_data_aoi %>% 
  select(codigo)

nanocuecas_data <- nanocuecas_data %>%
  select(ESTADO_SHT, CLAVE_SHT)

# 
dataset <- nanocuecas_data %>% 
  add_area_covar(geologia_data_aoi, "TIPO", column_prefix = "geologia")

dataset <- dataset %>% 
  add_area_covar(suelos_data_aoi, "CLAVE_IMP", column_prefix = "suelos")

dataset <- dataset %>% 
  add_area_covar(usvVII_data_aoi, "codigo", column_prefix = "usv")
