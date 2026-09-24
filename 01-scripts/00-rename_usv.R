INPUT_DIR <- "00-raw_data/03_Entrega_21092026/Shapefiles"
VARIABLES_DIR <- file.path(INPUT_DIR, "Variables")

usv_path <- file.path(VARIABLES_DIR, "usv250s7cw.shp")
mapping_path <- file.path(INPUT_DIR, "usv_rename.csv")
output_path <- file.path(VARIABLES_DIR, "usv250s7cw_renamed.shp")
unmatched_label <- "usv_sin_clasificar"

usv_data <- sf::st_read(usv_path, quiet = TRUE)
usv_mapping <- utils::read.csv(
  mapping_path,
  fileEncoding = "MACROMAN",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mapping_descriptions <- stringi::stri_trans_general(
  trimws(as.character(usv_mapping$DESCRIPCIO)),
  "Latin-ASCII"
)
mapping_labels <- trimws(as.character(usv_mapping$ETIQUETA))

usv_descriptions <- stringi::stri_trans_general(
  trimws(as.character(usv_data$DESCRIPCIO)),
  "Latin-ASCII"
)
mapping_matches <- match(usv_descriptions, mapping_descriptions)
unmatched_features <- is.na(mapping_matches)

usv_data$ETIQUETA <- mapping_labels[mapping_matches]
usv_data$ETIQUETA[unmatched_features] <- unmatched_label

unmatched_descriptions <- sort(unique(usv_descriptions[unmatched_features]))

message(
  "Matched ", sum(!unmatched_features), " of ", nrow(usv_data),
  " USV features to an ETIQUETA value."
)
message(
  "Assigned '", unmatched_label, "' to ", sum(unmatched_features),
  " feature(s) across ", length(unmatched_descriptions),
  " unmatched description(s)."
)

if (length(unmatched_descriptions) > 0) {
  message(
    "Unmatched shapefile descriptions:\n- ",
    paste(unmatched_descriptions, collapse = "\n- ")
  )
}

sf::st_write(
  usv_data,
  output_path,
  delete_layer = TRUE,
  quiet = TRUE
)

message("Wrote renamed USV shapefile to: ", output_path)
