#' Holds a list of urls from which to download schema files
#' 
#' @keywords internal, parsing
schema_urls <- function(){
  urls <- list()
  urls$grdi <- "https://raw.githubusercontent.com/cidgoh/pathogen-genomics-package/main/templates/grdi/schema.yaml"
  return(urls)
}

#' Download a schema file from a url
#' 
#' @keywords internal, parsing
download_schema <- function(yaml_url){
  file <- tempfile()
  download.file(url = yaml_url, destfile = file, method = "wget",  extra = "-4")
  return(file)
}

#' Load a schema file
#' 
#' @param file The path to the yaml schema file
#' @returns A schema file -> really just a parsed YAML loaded into R.
#' @keywords internal, parsing
#' @export
load_schema <- function(file){
  schema <- yaml::yaml.load_file(file)
  return(schema)
}

#' Get the slots of a schema
#' 
#' Get the slots (or columns) of a schema
#' 
#' @keywords internal, parsing
#' @param schema A schema file returned from [load_schema]
#' @keywords internal, parsing
slot_names <- function(schema){
  return(names(schema$slots))
}

#' Get the list of slot titles 
#' 
slot_titles <- function(schema, slots){
  unname(sapply(slots, function(x) schema$slots[[x]]$title))
}

#' Get the category of the slots
#' 
#' @inheritParams slot_names
#' @param slots A vector of slots
#' 
#' @return A named list of the slots with their category
#' @keywords internal, parsing
get_category <- function(schema, slots){
  check_slots_in_schema(schema,slots)
  usage <- schema$classes[[2]]$slot_usage
  sapply(FUN=function(x) usage[[x]]$slot_group, X=slots, USE.NAMES = FALSE)
}

#' Get the categories that the schema uses to group its slots
#' 
#' @inheritParams slot_names
#' 
#' @returns a vector of the categories of the schema
#' @keywords internal, parsing
categories <- function(schema){
  unique(get_category(schema, slots = slot_names(schema)))
}

#' Get all the slots that belong to the supplied category
#' 
#' @inheritParams slot_names
#' @param category A single category
#' 
#' @returns a vector of the slots that belong to the category
#' @keywords internal, parsing
get_slots_per_cat <- function(schema, category){
  if ( !all(category %in% categories(schema)) ){
    stop("supplied category not found in schema")
  }
  all_slots <- slot_names(schema)
  cats      <- get_category(schema, slots = all_slots)
  x <- all_slots[cats %in% category]
  return(x)
}

#' Retrieve the ranges of a slot
#' 
#' Returns the 'range' value of a slot.
#' 
#' @inheritParams slot_names
#' @param slot A name of a single slot
#' @returns A vector of the range names
#' @keywords internal, parsing
slot_ranges <- function(schema, slot){
  check_slots_in_schema(schema,slot)
  slots <- schema$slots
  any_ofs <- unlist(slots[[slot]]$any_of)
  ranges <-  unlist(slots[[slot]]$range)
  menus <- unname(c(any_ofs, ranges))
  return(menus)
}

#' Get the permissible values of a range, if they exist
#' 
#' All slots have ranges. Some of these ranges have values
#' Get the permissible values of each range,
#' 
#' @inheritParams slot_ranges
#' @param ranges A list of ranges, returned from [slot_ranges]
#' 
#' @returns A vector of permissible values for each range, NULL if none exist
#' @keywords internal, parsing
get_menu_values <- function(schema, ranges){
  x <- sapply(ranges, function(x) names(schema$enums[[x]]$permissible_values))
  return(x)
}

#' Get all the permissible values of a slot, if they exist.
#' 
#' This gets the permissible values of each range and appends them into 
#' a single vector. This includes the values of the NullValueMenu
#' 
#' @inheritParams slot_ranges
#' 
#' @returns A vector of permissible values for the slot
#' @keywords internal, parsing
get_permissible_values <- function(schema, slot){
  ranges <- slot_ranges(schema, slot)    
  vals <- get_menu_values(schema, ranges)
  return(unname(unlist(vals)))
}

#' Get the importance of a slot 
#' 
#' The standards often define the importance of a slot. 
#' Retrieve this for a single slot.
#' 
#' @inheritParams get_menu_values
#' @return The importance of the field
#' @keywords internal, parsing
get_field_importance <- function(schema, slot){
  check_slots_in_schema(schema,slot)
  slots <- schema$slots
  if (!is.null(slots[[slot]]$required)){
    return("required")
  } else if (!is.null(slots[[slot]]$recommended)){
    return("recommended")
  } else {
    return("normal")
  }
}

#' Get all the enums for each slot
#' 
#' Return a giant list of the slots and their corresponding enum entries
#' 
#' @inheritParams slot_names
#' @return A named list of the slots and their enums
#' @keywords internal, parsing
all_enums_per_col <- function(schema){
  cols <- names(schema$slots)
  # Remove AMR slots, but only if they are present.
  amr_index <- unlist(lapply(FUN = grep, X = amr_regexes(), x = cols))
  if (length(amr_index)>0) cols <- cols[-amr_index]
  vals <- sapply(FUN=get_permissible_values, X=cols, schema = schema)
  menus <- vals[lengths(vals)>0]
  return(menus)
}

#' Get the null value of the schema
#'
#' Returns a valid NULL value. If no query is supplied to search for a specific
#' ontology term, than a list of all the possible NULL values are supplied.
#'
#' @inheritParams slot_names
#' @param x Qery to search for a specific ontology term
#' @return A schema compliant null value
#' @keywords internal
#' @export
get_null_value <- function(schema, x=NULL){
  nulls <- names(schema$enums$NullValueMenu$permissible_values)
  if (is.null(x)){
    return(nulls)
  } else {
    value <- grep(x = nulls, pattern = x, value = TRUE, ignore.case = TRUE )
    if (length(value)==0){
      stop("No value found for query ", x, " in null value menu")
    } else {
      return(value)
    }
  }
}

#' Return the description and comments tags of a GRDI field
#'
#' @inheritParams get_menu_values
#'
#' @export
get_info <- function(schema, slot){
  check_slots_in_schema(schema,slot)
  f <- schema$slots[[slot]]
  cat("Description: ", str_wrap(f$description, width = 80), "\n")
  cat("Comments: ", str_wrap(f$comments, width = 80), "\n")
}

#' Get a list of all the ontology terms used in the standard.
#'
#' This function retrieves a list of all the enums (terms) used in a standard, 
#' returned as a named list. Organisational terms and nulls are removed from 
#' this list, for ease of searching through them
#'
#' @inheritParams slot_names
#' 
#' @return A named vector of all ontology terms from columns with menus, where
#'         the name of each value is the specification field it belongs to.
#'
#' @export
get_all_field_ontology_terms <- function(schema){
  x      <- all_enums_per_col(schema)
  values <- unlist(x, use.names = FALSE)
  nm     <- rep(names(x), times = lengths(x))
  names(values) <- nm
  org_terms <- grepl(x=values, "organizational term")
  nulls     <- values %in% get_null_value(schema)
  filtered  <- values[!(nulls | org_terms)]
  return(filtered)
}

#' Get the regex pattern of the slots
#' 
#' @inheritParams get_category
#' @returns a named vector of the regex that applies to the slot, or NA if none present
#' @keywords internal, parsing
get_pattern <- function(schema, slots){
    x <- sapply(slots, function(x) schema$slots[[x]]$pattern, simplify = T)
    x <- sapply(x    , function(x) if (is.null(x)) return(NA) else return(x))
    return(x)
}

#' Get the type of the slot
#'
#' Get the type of the slot, inferred from the 'range' value when 'Menu' items are removed. 
#' If all the 'Menu' items are removed, it must be a menu type so return that instead.
#' 
#' @inheritParams get_menu_values
#' @returns The type of the slot
#' @keywords internal, parsing
get_slot_type <- function(schema, slot){
  r <- slot_ranges(schema, slot)
  slot_type <- r[!grepl(pattern = 'Menu$', x=r)]
  if (length(slot_type) == 0) slot_type <- 'Menu'
  return(slot_type)
}

#' Determine if the slots are identifier columns.
#' 
#' @inheritParams get_category
#' @returns Logical vector, TRUE if idendifier, FALSE if not.
#' @keywords internal, parsing
is_identifier <- function(schema, slots){
  x <- sapply(slots, 
              function(x){ id <- schema$slots[[x]]$identifier
                           if (!is.null(id)) return(TRUE) else return(FALSE)})
  return(x)
}

#' Check if the slots are multivalued
#' 
#' @inheritParams get_category
#' @returns Logical vector, TRUE if the list is multivalued, FALSE if not.
check_multivalues <- function(schema, slots){
  return(sapply(slots, function(x) if (!is.null(schema$slots[[x]]$multivalued)) TRUE else FALSE))
}

#' Get the date range of a slot
#'
#' This returns the acceptable date range of a date slot. 
#' Slot date ranges appear to be coded in a strange 'todo' section of the slot.
#' Parse this section (will likely error if not present) using '>' and '<' to 
#' determine which is the max value and min value.
#' 
#' @inheritParams slot_ranges
#' @returns A [lubridate::interval] of the acceptable date range.
#' @keywords internal, parsing
get_date_range <- function(schema, slot, dataframe){
  x <- schema$slots[[slot]]$todos
  min_val <- x[grepl(x = x, ">")]
  max_val <- x[grepl(x = x, "<")]
  date_range <- lapply(X=list(min=min_val, max=max_val), function(x) if (length(x)==0) x <- NA else x)
  return(date_range)
}

#' Convert the min or max range tokens to actual dates to check against
#'
#' @param token A string, taken from the "todo" value of a date slot.
#' @param slot The slot name, used for
#' @param dataframe The full dataframe of the data, we need this to extract the dates in the event that
#'                  the min value is acutally another column.
#' @returns description
#' @keywords internal, parsing
convert_range_token_to_date <- function(token, slot, dataframe){
  # Strip >=< characters, this can cause downstream errors in datetime parsing.
  token <- gsub(x=token, pattern = "[><=]{1,2}", "")
  #If brakets, interpret this either as today or get the value from the dataframe
  if        (grepl(x=token, "today")){
    vals <- rep(lubridate::today(), nrow(dataframe))
  } else if (grepl(x=token, "[{}]")){
    col <- stringr::str_extract(string=token, pattern = "\\{([a-zA-Z_]+)\\}", group = 1)
    requested_vals <- dataframe[[col]]
    if (is.null(requested_vals)){
      log_message("warning", slot = slot, "Need slot ", col, " to perform range check, but data missing!")
      vals <- rep(NA, nrow(dataframe))
    } else {
      vals <- suppressWarnings(parse_date(requested_vals))
    }
  } else if (!is.na(suppressWarnings(parse_date(token)))){
    vals <- rep(parse_date(token), nrow(dataframe))
  } else if (is.na(token)){
    log_message("note", slot = slot, "No minimum value found in date range. Setting this to 1900-01-01")
    vals <- rep(lubridate::ymd("1900-01-01"), nrow(dataframe))
  } else {stop("Unhandled token in one of the values in todo, for date range checking: ", token)}
  return(vals) 
}

#' Get the minimum and maximum value of a slot
#' 
#' @inheritParams slot_ranges
#' @returns vector of length 2 (min,max), or NA if there is no range
min_max_value <- function(schema, slot){
  min <- schema$slots[[slot]]$minimum_value
  max <- schema$slots[[slot]]$maximum_value
  if (is.null(min)) min <- NA
  if (is.null(max)) max <- NA
  min_max <- as.numeric(c(min,max))
  if (all(is.na(min_max))) return(NA) else return(min_max)
}

#' Get the identification column, returning the 1st one in there are multiple.
#' 
#' @inheritParams slot_names
get_first_identification_col <- function(schema){
  slots <- slot_names(schema)
  x <- is_identifier(schema, slot_names(schema))
  id_cols <- slots[x]
  if (length(id_cols)>1) message("More then 1 identifier column detected, using first")
  return(id_cols[1])
}
