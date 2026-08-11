suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(shinycssloaders)
  library(cachem)
  library(DT)
})

# ============================================================
# Developing Mouse Visual Cortex Gene Expression Viewer
# ============================================================

app_title <- "Developing Mouse Visual Cortex Gene Expression Viewer"

s3_bucket_name <- "sea-ad-single-cell-profiling"
s3_region <- "us-west-2"
s3_prefix <- paste0(
  "Multiregion_2026/previous_objects/parquet/",
  "TEST_dev_mouse"
)

bucket <- arrow::s3_bucket(
  bucket = s3_bucket_name,
  region = s3_region,
  anonymous = TRUE
)

metadata_source <- bucket$path(
  file.path(s3_prefix, "metadata", "obs_metadata.parquet")
)
gene_map_source <- bucket$path(
  file.path(s3_prefix, "metadata", "gene_map.parquet")
)
cpm_source <- bucket$path(
  file.path(s3_prefix, "metadata", "CPM_scaling.csv")
)
statistics_source <- bucket$path(
  file.path(s3_prefix, "statistics", "mouse_gene_summary_statistics.csv")
)

read_startup_data <- function() {
  tryCatch(
    list(
      metadata = arrow::read_parquet(
        metadata_source,
        as_data_frame = TRUE
      ),
      gene_map = arrow::read_parquet(
        gene_map_source,
        as_data_frame = TRUE
      ),
      cpm = arrow::read_csv_arrow(
        cpm_source,
        as_data_frame = TRUE
      ),
      statistics = arrow::read_csv_arrow(
        statistics_source,
        as_data_frame = TRUE
      )
    ),
    error = function(e) {
      stop(
        paste(
          "Public AWS data connection failed.",
          "The application could not read the mouse visual cortex files.",
          "Please verify that the S3 prefix is public and available.",
          "Original error:",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )
}

startup <- read_startup_data()
metadata <- as.data.frame(startup$metadata, stringsAsFactors = FALSE)
gene_map <- as.data.frame(startup$gene_map, stringsAsFactors = FALSE)
cpm_scaling <- as.data.frame(startup$cpm, stringsAsFactors = FALSE)
gene_statistics <- as.data.frame(
  startup$statistics,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
rm(startup)

# ============================================================
# Validate and align supporting data
# ============================================================

required_metadata <- c("sample_id", "library_label")
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata) > 0) {
  stop(
    "obs_metadata.parquet is missing: ",
    paste(missing_metadata, collapse = ", ")
  )
}

required_gene_map <- c("gene_symbol", "gene_column")
missing_gene_map <- setdiff(required_gene_map, names(gene_map))
if (length(missing_gene_map) > 0) {
  stop(
    "gene_map.parquet is missing: ",
    paste(missing_gene_map, collapse = ", ")
  )
}

required_cpm <- c("library_label", "CPM_scaling_factor")
missing_cpm <- setdiff(required_cpm, names(cpm_scaling))
if (length(missing_cpm) > 0) {
  stop(
    "CPM_scaling.csv is missing: ",
    paste(missing_cpm, collapse = ", ")
  )
}

if (nrow(metadata) != nrow(cpm_scaling)) {
  stop("CPM_scaling.csv and obs_metadata.parquet have different row counts.")
}

if (!identical(
  as.character(metadata$library_label),
  as.character(cpm_scaling$library_label)
)) {
  stop("CPM_scaling.csv library order does not match obs_metadata.parquet.")
}

metadata <- metadata |>
  mutate(
    sample_id = as.character(sample_id),
    library_label = as.character(library_label),
    CPM_scaling_factor = as.numeric(cpm_scaling$CPM_scaling_factor)
  )

if (anyDuplicated(metadata$sample_id)) {
  stop("sample_id must be unique in obs_metadata.parquet.")
}

names(gene_statistics) <- make.unique(
  names(gene_statistics),
  sep = "_"
)

possible_statistics_gene_columns <- c(
  "gene_symbol", "Gene Symbol", "Gene", "gene", "symbol"
)
statistics_gene_column <- possible_statistics_gene_columns[
  possible_statistics_gene_columns %in% names(gene_statistics)
][1]

if (is.na(statistics_gene_column)) {
  stop(
    "mouse_gene_summary_statistics.csv must contain one of: ",
    paste(possible_statistics_gene_columns, collapse = ", ")
  )
}

gene_statistics[[statistics_gene_column]] <- as.character(
  gene_statistics[[statistics_gene_column]]
)

gene_map <- gene_map |>
  transmute(
    gene_symbol = as.character(gene_symbol),
    gene_column = as.character(gene_column)
  ) |>
  filter(
    !is.na(gene_symbol),
    nzchar(gene_symbol),
    !is.na(gene_column),
    nzchar(gene_column)
  ) |>
  distinct(gene_symbol, .keep_all = TRUE) |>
  arrange(gene_symbol)

available_genes <- gene_map$gene_symbol
gene_symbol_to_key <- setNames(
  gene_map$gene_column,
  gene_map$gene_symbol
)

# ============================================================
# Metadata fields and helpers
# ============================================================

first_existing <- function(candidates, data_names = names(metadata)) {
  hit <- candidates[candidates %in% data_names]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

age_field <- first_existing(c(
  "age", "age_label", "donor_age", "developmental_age",
  "age_in_days", "age_days"
))
region_field <- first_existing(c(
  "region_of_interest_label", "Brain Region", "brain_region",
  "region_label", "structure"
))
class_field <- first_existing(c("class", "class_label", "Class"))
subclass_field <- first_existing(c("subclass", "subclass_label", "Subclass"))
cluster_field <- first_existing(c("cluster_alias", "cluster", "cluster_label"))
subcluster_field <- first_existing(c(
  "subcluster_alias", "subcluster", "subcluster_label"
))

if (is.na(age_field)) {
  stop("No age field was found in obs_metadata.parquet.")
}
if (is.na(region_field)) {
  stop("No brain-region field was found in obs_metadata.parquet.")
}
if (is.na(cluster_field)) {
  stop("No cluster field was found in obs_metadata.parquet.")
}

preferred_fields <- unique(c(
  class_field,
  subclass_field,
  cluster_field,
  subcluster_field,
  region_field,
  age_field,
  "sex",
  "donor_sex",
  "library_label",
  "donor_label",
  "donor_id",
  "developmental_compartment",
  "sample_type"
))

preferred_fields <- preferred_fields[
  !is.na(preferred_fields) &
    preferred_fields %in% names(metadata)
]

categorical_fields <- preferred_fields[vapply(
  metadata[preferred_fields],
  function(x) {
    is.character(x) || is.factor(x) ||
      dplyr::n_distinct(x, na.rm = TRUE) <= 100
  },
  logical(1)
)]

plot_fields <- setdiff(
  categorical_fields,
  c("library_label", "donor_label", "donor_id")
)
filter_fields <- categorical_fields
cell_type_fields <- unique(na.omit(c(
  class_field,
  subclass_field,
  cluster_field,
  subcluster_field
)))

parse_age <- function(x) {
  text <- toupper(trimws(as.character(x)))
  out <- suppressWarnings(as.numeric(text))
  
  embryonic <- grepl("^E\\s*[0-9]", text)
  postnatal <- grepl("^P\\s*[0-9]", text)
  
  extract_number <- function(values) {
    suppressWarnings(as.numeric(
      sub(".*?(-?[0-9]+(?:\\.[0-9]+)?).*", "\\1", values)
    ))
  }
  
  out[embryonic] <- extract_number(text[embryonic]) - 100
  out[postnatal] <- extract_number(text[postnatal])
  out
}

capped_age_positions <- function(x, maximum_step = 3) {
  labels <- age_level_order(x)
  raw_positions <- parse_age(labels)
  
  valid <- is.finite(raw_positions)
  labels <- labels[valid]
  raw_positions <- raw_positions[valid]
  
  if (length(raw_positions) == 0) {
    return(setNames(numeric(), character()))
  }
  
  if (length(raw_positions) == 1) {
    return(setNames(0, labels))
  }
  
  raw_steps <- diff(raw_positions)
  capped_steps <- pmin(raw_steps, maximum_step)
  
  positions <- c(0, cumsum(capped_steps))
  setNames(positions, labels)
}

age_level_order <- function(x) {
  values <- unique(as.character(x))
  values <- values[!is.na(values) & nzchar(values)]
  numeric_values <- parse_age(values)
  values[order(numeric_values, values, na.last = TRUE)]
}

field_levels <- function(data, field) {
  values <- unique(as.character(data[[field]]))
  values <- values[!is.na(values) & nzchar(values)]
  if (identical(field, age_field)) {
    age_level_order(values)
  } else {
    sort(values)
  }
}

factor_field <- function(data, field) {
  data[[field]] <- factor(
    as.character(data[[field]]),
    levels = field_levels(data, field)
  )
  data
}

field_colors <- function(data, field) {
  levels <- field_levels(data, field)
  color_column_candidates <- unique(c(
    paste0(field, "_color"),
    paste0(field, "_color_hex_triplet"),
    paste0(gsub(" ", "_", field), "_color"),
    "color_hex_triplet"
  ))
  color_column <- first_existing(
    color_column_candidates,
    names(data)
  )
  
  if (!is.na(color_column)) {
    lookup <- data |>
      transmute(
        level = as.character(.data[[field]]),
        color = as.character(.data[[color_column]])
      ) |>
      filter(
        !is.na(level), nzchar(level),
        !is.na(color), grepl("^#", color)
      ) |>
      distinct(level, .keep_all = TRUE)
    
    colors <- setNames(lookup$color, lookup$level)
    missing <- setdiff(levels, names(colors))
    colors[missing] <- scales::hue_pal()(length(missing))
    return(colors[levels])
  }
  
  setNames(scales::hue_pal()(length(levels)), levels)
}

manual_color_scale <- function(data, field, aesthetic = "color") {
  colors <- field_colors(data, field)
  if (identical(aesthetic, "fill")) {
    scale_fill_manual(values = colors, drop = TRUE, na.value = "#808080")
  } else {
    scale_color_manual(values = colors, drop = TRUE, na.value = "#808080")
  }
}

fields_within_limit <- function(data, fields, maximum) {
  fields[vapply(
    fields,
    function(field) {
      field %in% names(data) &&
        n_distinct(data[[field]], na.rm = TRUE) <= maximum
    },
    logical(1)
  )]
}

make_cache_key <- function(gene_symbol) {
  paste0(
    "gene",
    paste(as.integer(charToRaw(enc2utf8(gene_symbol))), collapse = "")
  )
}

gene_cache <- cachem::cache_mem(
  max_size = 1024 * 1024^2,
  max_age = Inf
)

default_gene <- if ("Reln" %in% available_genes) {
  "Reln"
} else {
  available_genes[[1]]
}
default_x <- if ("donor_age" %in% plot_fields) {
  "donor_age"
} else {
  age_field
}

default_second <- if (!is.na(subclass_field)) {
  subclass_field
} else {
  cluster_field
}
default_facet <- if (!is.na(subclass_field)) subclass_field else cluster_field
default_color <- region_field

# ============================================================
# Header and UI
# ============================================================

app_header <- div(
  class = "mouse-header",
  tags$a(
    class = "mouse-logo-link",
    href = "https://alleninstitute.org/",
    target = "_blank",
    tags$img(
      src = "allen_institute_logo.svg",
      alt = app_title,
      class = "mouse-logo"
    )
  ),
  tags$span(class = "mouse-title", app_title),
  div(
    class = "mouse-header-links",
    tags$a(
      href = "https://alleninstitute.github.io/abc_atlas_access/descriptions/Dev-Mouse-Vis-Cortex-dataset.html",
      target = "_blank",
      icon("brain"),
      tags$span("Dataset")
    ),
    tags$a(
      href = "https://brain-map.org",
      target = "_blank",
      icon("globe"),
      tags$span("Brain Map")
    )
  )
)

ui <- page_sidebar(
  title = app_header,
  theme = bs_theme(
    version = 5,
    primary = "#214E68",
    success = "#4B9B58"
  ),
  
  tags$head(
    tags$title(app_title),
    tags$link(
      rel = "icon",
      type = "image/png",
      href = "allen_institute_logo.svg"
    ),
    tags$style(HTML(
      "
      html, body { width: 100%; overflow-x: hidden; }
      .navbar, header.navbar {
        background: #111827 !important;
        min-height: 64px !important;
        border: 0 !important;
      }
      .navbar-brand { width: 100% !important; padding: 0 !important; }
      .mouse-header {
        position: relative; display: flex; align-items: center;
        justify-content: space-between; width: 100%; height: 64px;
        padding: 0 16px;
      }
      .mouse-header-spacer {
        width: 210px;
      }
      
      .sidebar-logo-container {
        width: 100%;
        text-align: center;
        margin: 0 0 7px 0;
      }
      
      .sidebar-logo {
        display: block;
        width: 100%;
        max-width: 345px;
        max-height: 250px;
        object-fit: contain;
        margin: 0 auto;
      }
      .mouse-logo-link { width: 180px; display: flex; align-items: center; }
      .mouse-logo { height: 52px; max-width: 175px; object-fit: contain; }
      .mouse-title {
        position: absolute; left: 50%; transform: translateX(-50%);
        color: white; font-size: 21px; font-weight: 700; white-space: nowrap;
      }
      .mouse-header-links { display: flex; gap: 6px; width: 210px; justify-content: flex-end; }
      .mouse-header-links a {
        color: white !important; text-decoration: none; font-size: 13px;
        padding: 18px 8px;
      }
      .bslib-sidebar-layout > .sidebar,
      .bslib-sidebar-layout > aside {
        background: #214E68 !important; color: white !important;
        padding: 8px 10px !important;
      }
      .bslib-sidebar-layout > .sidebar label,
      .bslib-sidebar-layout > .sidebar h4,
      .bslib-sidebar-layout > .sidebar .form-check-label,
      .bslib-sidebar-layout > .sidebar .help-block {
        color: white !important;
      }
      .compact-controls h4, .gene-controls h4 {
        font-size: 1rem; margin: 4px 0;
      }
      .compact-controls hr { margin: 5px 0; border-color: rgba(255,255,255,.3); }
      .compact-controls .shiny-input-container { margin-bottom: 3px !important; }
      .compact-controls .shiny-input-container > label,
      .compact-controls .form-check-label,
      .compact-controls input[type='checkbox'] + label {
        font-size: .75rem !important; line-height: 1.08 !important;
      }
      .compact-controls .form-select,
      .compact-controls .form-control,
      .compact-controls .selectize-input {
        min-height: 28px !important; font-size: .78rem !important;
      }
      .repeat-plot-note {
        color: rgba(255, 255, 255, 0.85);
        font-size: 0.67rem;
        line-height: 1.1;
        text-align: center;
        margin: 3px 0 4px 0;
      }
      #filter_values + .selectize-control .selectize-input,
      #filter_values + .selectize-control .selectize-dropdown,
      #filter_values-selectized + .selectize-dropdown {
        background: white !important; color: #111827 !important;
      }
      .gene-status { font-size: .76rem; min-height: 1.2rem; margin: 2px 0; }
      .btn-primary, .btn-outline-primary {
        background: #4B9B58 !important; border-color: white !important;
        color: white !important; font-weight: 600;
      }
      .filter-row { display: grid; grid-template-columns: 1fr 1fr; gap: 6px; }
      .gene-table-panel { padding: 8px; overflow-x: auto; }
      .gene-table-panel table.dataTable { font-size: .82rem; }
      .expression-panel { width: 100%; overflow: hidden; }
      .expression-title {
        min-height: 42px; padding: 9px 13px; background: #111827;
        color: white; font-weight: 600;
      }
      .plot-wrapper {
        width: 100%; height: calc(100vh - 245px); min-height: 380px;
        padding: 6px; overflow: hidden; box-sizing: border-box;
      }
      .plot-wrapper .shiny-plot-output,
      .plot-wrapper .shiny-spinner-output-container,
      .plot-wrapper .load-container {
        width: 100% !important; height: 100% !important;
        min-height: 0 !important; overflow: hidden !important;
      }
      @media (max-width: 1000px) {
        .mouse-title { font-size: 15px; max-width: 55vw; overflow: hidden; text-overflow: ellipsis; }
        .mouse-header-links span { display: none; }
        .mouse-header-links { width: 80px; }
        .mouse-logo-link { width: 110px; }
        .mouse-logo { max-width: 105px; }
        .mouse-header-spacer {
          width: 80px;
        }
      }
      "
    ))
  ),
  
  sidebar = sidebar(
    width = 320,
    div(
      class = "sidebar-logo-container",
      tags$a(
        href = "https://alleninstitute.github.io/abc_atlas_access/descriptions/Dev-Mouse-Vis-Cortex-dataset.html",
        target = "_blank",
        tags$img(
          src = "Developing_Mouse_Visual_Cortex_logo.png",
          alt = app_title,
          class = "sidebar-logo"
        )
      )
    ),
    helpText(
      "Choose a gene, retrieve its data, select plot options, and generate the plot."
    ),
    div(
      class = "gene-controls",
      h4("Gene"),
      selectizeInput(
        "gene",
        "Gene symbol",
        choices = NULL,
        selected = NULL,
        options = list(
          placeholder = "Type a gene symbol",
          maxOptions = 50,
          create = FALSE
        )
      ),
      actionButton(
        "get_gene_data",
        "Get gene data",
        class = "btn-outline-primary",
        width = "100%"
      ),
      div(class = "gene-status", textOutput("gene_status"))
    ),
    div(
      class = "compact-controls",
      hr(),
      h4("Filter and scale"),
      div(
        class = "filter-row",
        selectInput(
          "filter_field",
          "Filter metadata",
          choices = c("No filter" = "none", filter_fields),
          selected = "none"
        ),
        selectizeInput(
          "filter_values",
          "Include values",
          choices = NULL,
          selected = NULL,
          multiple = TRUE,
          options = list(
            placeholder = "Select values",
            maxOptions = 1000,
            closeAfterSelect = TRUE
          )
        )
      ),
      checkboxInput(
        "omit_zero_values",
        "Discard observations with zero counts",
        FALSE
      ),
      checkboxInput(
        "log_scale",
        "Plot ln(CPM + 1)",
        TRUE
      ),
      hr(),
      h4("Plot"),
      selectInput(
        "plot_type",
        "Plot type",
        choices = c(
          "Trajectory across development" = "trajectory",
          "Heatmap of mean expression" = "heatmap",
          "Dot plot of mean expression" = "dot",
          "Violin plot with observations" = "violin"
        ),
        selected = "trajectory"
      ),
      conditionalPanel(
        condition = "input.plot_type == 'trajectory'",
        selectInput(
          "progression_variable",
          "Developmental progression axis",
          choices = age_field,
          selected = age_field
        ),
        selectInput(
          "facet_variable",
          "Facet by (maximum 30 values)",
          choices = c(cell_type_fields, region_field),
          selected = default_facet
        ),
        selectInput(
          "color_variable",
          "Color by (maximum 15 values)",
          choices = region_field,
          selected = default_color
        ),
        selectInput(
          "smoother",
          "Trend line",
          choices = c(
            "LOESS smoother" = "loess",
            "Linear fit" = "lm",
            "None" = "none"
          ),
          selected = "loess"
        ),
        checkboxInput(
          "show_points",
          "Show individual observations",
          TRUE
        )
      ),
      conditionalPanel(
        condition = "input.plot_type != 'trajectory'",
        selectInput(
          "x_variable",
          "Horizontal axis",
          choices = plot_fields,
          selected = default_x
        ),
        selectInput(
          "second_dimension",
          "Second dimension",
          choices = setdiff(plot_fields, default_x),
          selected = default_second
        )
      ),
      actionButton(
        "make_plot",
        "Generate plot",
        class = "btn-primary",
        width = "100%"
      ),
      div(
        class = "repeat-plot-note",
        "If nothing happens, press ^ again."
      ),
      actionButton(
        "reset_defaults",
        "Reset plot options",
        class = "btn-outline-light btn-sm mt-1",
        width = "100%"
      )
    )
  ),
  
  navset_card_tab(
    id = "main_tabs",
    nav_panel(
      "Select gene of interest",
      div(
        class = "gene-table-panel",
        DT::DTOutput("gene_statistics_table")
      )
    ),
    nav_panel(
      "Plot gene of interest",
      div(
        class = "expression-panel",
        div(
          class = "expression-title",
          textOutput("plot_title")
        ),
        div(
          class = "plot-wrapper",
          shinycssloaders::withSpinner(
            plotOutput(
              "expression_plot",
              width = "100%",
              height = "100%"
            ),
            type = 8,
            color = "#4B9B58",
            size = 1
          )
        )
      )
    )
  )
)

# ============================================================
# Server
# ============================================================

server <- function(input, output, session) {
  gene_data <- reactiveVal(NULL)
  loaded_gene <- reactiveVal(NULL)
  requested_gene <- reactiveVal(default_gene)
  gene_status_message <- reactiveVal("No gene retrieved yet.")
  
  output$gene_statistics_table <- DT::renderDT({
    table_data <- gene_statistics
    
    numeric_columns <- names(table_data)[vapply(
      table_data,
      is.numeric,
      logical(1)
    )]
    
    first_numeric_index <- if (length(numeric_columns) > 0) {
      match(numeric_columns[[1]], names(table_data)) - 1L
    } else {
      0L
    }
    
    table <- DT::datatable(
      table_data,
      rownames = FALSE,
      filter = "top",
      selection = list(
        mode = "single",
        target = "row"
      ),
      options = list(
        pageLength = 10,
        lengthMenu = c(10, 20, 50, 100),
        scrollX = TRUE,
        autoWidth = TRUE,
        stateSave = FALSE,
        searchHighlight = TRUE,
        order = list(
          list(first_numeric_index, "desc")
        )
      ),
      class = "compact stripe hover"
    )
    
    if (length(numeric_columns) > 0) {
      table <- DT::formatRound(
        table,
        columns = numeric_columns,
        digits = 3
      )
    }
    
    table
  }, server = TRUE)
  
  observeEvent(input$gene_statistics_table_rows_selected, {
    selected_row <- input$gene_statistics_table_rows_selected
    req(length(selected_row) == 1)
    
    selected_gene <- gene_statistics[[statistics_gene_column]][selected_row]
    req(length(selected_gene) == 1, !is.na(selected_gene), nzchar(selected_gene))
    
    if (selected_gene %in% available_genes) {
      requested_gene(selected_gene)
      updateSelectizeInput(
        session,
        "gene",
        choices = available_genes,
        selected = selected_gene,
        server = TRUE
      )
      gene_status_message(
        paste("Selected", selected_gene, "from the statistics table.")
      )
    } else {
      showNotification(
        paste("This gene is not available:", selected_gene),
        type = "warning"
      )
    }
  }, ignoreInit = TRUE)
  
  session$onFlushed(function() {
    updateSelectizeInput(
      session,
      "gene",
      choices = available_genes,
      selected = default_gene,
      server = TRUE
    )
  }, once = TRUE)
  
  observeEvent(input$gene, {
    if (
      length(input$gene) == 1 &&
      nzchar(input$gene) &&
      input$gene %in% available_genes
    ) {
      requested_gene(input$gene)
    }
  }, ignoreInit = TRUE)
  
  read_gene <- function(gene_symbol) {
    cache_key <- make_cache_key(gene_symbol)
    cached <- gene_cache$get(cache_key, missing = NULL)
    if (!is.null(cached)) return(cached)
    
    gene_key <- unname(gene_symbol_to_key[[gene_symbol]])
    if (is.null(gene_key) || is.na(gene_key)) {
      stop("No gene partition key was found for: ", gene_symbol)
    }
    
    gene_source <- bucket$path(
      file.path(
        s3_prefix,
        "counts_by_gene",
        paste0("gene=", gene_key)
      )
    )
    
    # open_dataset() deliberately supports one or many Parquet files in the
    # selected gene partition.
    gene_counts <- arrow::open_dataset(
      gene_source,
      format = "parquet",
      partitioning = NULL,
      unify_schemas = FALSE
    ) |>
      select(sample_id, summed_counts) |>
      collect() |>
      transmute(
        sample_id = as.character(sample_id),
        summed_counts = as.numeric(summed_counts)
      )
    
    if (nrow(gene_counts) == 0) {
      stop("No count data were found for gene: ", gene_symbol)
    }
    
    joined <- metadata |>
      inner_join(gene_counts, by = "sample_id") |>
      filter(
        is.finite(summed_counts),
        is.finite(CPM_scaling_factor),
        CPM_scaling_factor >= 0
      ) |>
      mutate(
        CPM = summed_counts * CPM_scaling_factor
      )
    
    if (nrow(joined) == 0) {
      stop("Gene data were found, but no rows matched valid metadata and CPM scaling.")
    }
    
    gene_cache$set(cache_key, joined)
    joined
  }
  
  observeEvent(input$get_gene_data, {
    gene_to_load <- requested_gene()
    req(gene_to_load %in% available_genes)
    gene_status_message(paste("Retrieving", gene_to_load, "..."))
    
    tryCatch(
      {
        loaded <- withProgress(
          message = paste("Retrieving", gene_to_load),
          value = 0.25,
          {
            result <- read_gene(gene_to_load)
            incProgress(0.65)
            result
          }
        )
        gene_data(loaded)
        loaded_gene(gene_to_load)
        updateSelectizeInput(
          session,
          "gene",
          choices = available_genes,
          selected = gene_to_load,
          server = TRUE
        )
        gene_status_message(
          paste0(
            "Loaded ", gene_to_load, " (",
            comma(nrow(loaded)), " observations)."
          )
        )
      },
      error = function(e) {
        gene_status_message(
          paste("Unable to load", gene_to_load, ":", conditionMessage(e))
        )
        showNotification(
          paste("Unable to load gene data:", conditionMessage(e)),
          type = "error",
          duration = NULL
        )
      }
    )
  })
  
  output$gene_status <- renderText(gene_status_message())
  
  observeEvent(input$filter_field, {
    if (is.null(input$filter_field) || input$filter_field == "none") {
      updateSelectizeInput(
        session,
        "filter_values",
        choices = character(),
        selected = character(),
        server = TRUE
      )
      return()
    }
    
    choices <- field_levels(metadata, input$filter_field)
    updateSelectizeInput(
      session,
      "filter_values",
      choices = choices,
      selected = NULL,
      server = TRUE
    )
  }, ignoreInit = FALSE)
  
  control_metadata <- reactive({
    if (
      is.null(input$filter_field) ||
      input$filter_field == "none" ||
      length(input$filter_values) == 0
    ) {
      return(metadata)
    }
    
    out <- metadata |>
      filter(
        as.character(.data[[input$filter_field]]) %in%
          as.character(input$filter_values)
      )
    if (nrow(out) == 0) metadata else out
  })
  
  observeEvent(control_metadata(), {
    data <- control_metadata()
    facet_choices <- fields_within_limit(
      data,
      unique(c(cell_type_fields, region_field)),
      30
    )
    color_choices <- fields_within_limit(
      data,
      plot_fields,
      15
    )
    
    if (length(facet_choices) == 0) facet_choices <- default_facet
    if (length(color_choices) == 0) color_choices <- default_color
    
    current_facet <- isolate(input$facet_variable)
    if (!current_facet %in% facet_choices) current_facet <- facet_choices[[1]]
    
    current_color <- isolate(input$color_variable)
    if (!current_color %in% color_choices) current_color <- color_choices[[1]]
    
    updateSelectInput(
      session,
      "facet_variable",
      choices = facet_choices,
      selected = current_facet
    )
    updateSelectInput(
      session,
      "color_variable",
      choices = color_choices,
      selected = current_color
    )
  }, ignoreInit = FALSE)
  
  observeEvent(
    list(input$plot_type, input$x_variable, control_metadata()),
    {
      req(input$plot_type)
      data <- control_metadata()
      
      dimension_choices <- if (input$plot_type %in% c("heatmap", "dot")) {
        plot_fields
      } else {
        fields_within_limit(data, plot_fields, 30)
      }
      if (length(dimension_choices) < 2) dimension_choices <- plot_fields
      
      current_x <- isolate(input$x_variable)
      if (!current_x %in% dimension_choices) current_x <- dimension_choices[[1]]
      
      updateSelectInput(
        session,
        "x_variable",
        choices = dimension_choices,
        selected = current_x
      )
      
      second_choices <- setdiff(dimension_choices, current_x)
      current_second <- isolate(input$second_dimension)
      if (!current_second %in% second_choices) {
        current_second <- if (region_field %in% second_choices) {
          region_field
        } else {
          second_choices[[1]]
        }
      }
      
      updateSelectInput(
        session,
        "second_dimension",
        choices = second_choices,
        selected = current_second
      )
    },
    ignoreInit = FALSE
  )
  
  observeEvent(input$reset_defaults, {
    updateSelectInput(session, "filter_field", selected = "none")
    updateSelectizeInput(
      session,
      "filter_values",
      choices = character(),
      selected = character(),
      server = TRUE
    )
    updateCheckboxInput(session, "omit_zero_values", value = FALSE)
    updateCheckboxInput(session, "log_scale", value = TRUE)
    updateSelectInput(session, "plot_type", selected = "trajectory")
    updateSelectInput(session, "progression_variable", selected = age_field)
    updateSelectInput(session, "facet_variable", selected = default_facet)
    updateSelectInput(session, "color_variable", selected = default_color)
    updateSelectInput(session, "smoother", selected = "loess")
    updateCheckboxInput(session, "show_points", value = TRUE)
    updateSelectInput(session, "x_variable", selected = default_x)
    updateSelectInput(session, "second_dimension", selected = default_second)
  })
  
  observeEvent(input$make_plot, {
    bslib::nav_select(
      id = "main_tabs",
      selected = "Plot gene of interest",
      session = session
    )
  }, ignoreInit = TRUE)
  
  plot_settings <- eventReactive(input$make_plot, {
    req(gene_data(), loaded_gene(), input$plot_type)
    list(
      plot_type = input$plot_type,
      filter_field = input$filter_field,
      filter_values = input$filter_values,
      omit_zero_values = isTRUE(input$omit_zero_values),
      log_scale = isTRUE(input$log_scale),
      progression_variable = input$progression_variable,
      facet_variable = input$facet_variable,
      color_variable = input$color_variable,
      smoother = input$smoother,
      show_points = isTRUE(input$show_points),
      x_variable = input$x_variable,
      second_dimension = input$second_dimension
    )
  }, ignoreInit = TRUE)
  
  filtered_data <- reactive({
    settings <- plot_settings()
    data <- gene_data()
    req(data)
    
    if (
      !is.null(settings$filter_field) &&
      settings$filter_field != "none" &&
      length(settings$filter_values) > 0
    ) {
      data <- data |>
        filter(
          as.character(.data[[settings$filter_field]]) %in%
            as.character(settings$filter_values)
        )
    }
    
    if (settings$omit_zero_values) {
      data <- data |> filter(summed_counts > 0)
    }
    
    validate(need(
      nrow(data) >= 2,
      "Fewer than two observations remain after filtering."
    ))
    
    data |>
      mutate(
        plotted_expression = if (settings$log_scale) {
          log1p(CPM)
        } else {
          CPM
        }
      )
  })
  
  summarized_data <- reactive({
    settings <- plot_settings()
    data <- filtered_data()
    req(settings$x_variable, settings$second_dimension)
    validate(need(
      settings$x_variable != settings$second_dimension,
      "Choose two different plot dimensions."
    ))
    
    data <- factor_field(data, settings$x_variable)
    data <- factor_field(data, settings$second_dimension)
    
    out <- data |>
      filter(
        !is.na(.data[[settings$x_variable]]),
        !is.na(.data[[settings$second_dimension]]),
        is.finite(plotted_expression)
      ) |>
      group_by(across(all_of(c(
        settings$x_variable,
        settings$second_dimension
      )))) |>
      summarise(
        expression = mean(plotted_expression, na.rm = TRUE),
        .groups = "drop"
      )
    
    validate(need(nrow(out) > 0, "No valid groups remain for this plot."))
    out
  })
  
  output$expression_plot <- renderPlot({
    settings <- plot_settings()
    data <- filtered_data()
    
    y_label <- if (settings$log_scale) "ln(CPM + 1)" else "Counts per million"
    
    if (settings$plot_type == "trajectory") {
      data <- factor_field(data, settings$facet_variable)
      data <- factor_field(data, settings$color_variable)
      age_position_lookup <- capped_age_positions(
        data[[settings$progression_variable]],
        maximum_step = 3
      )
      
      age_labels <- names(age_position_lookup)
      age_breaks <- unname(age_position_lookup)
      
      data$progression_value <- unname(
        age_position_lookup[
          as.character(data[[settings$progression_variable]])
        ]
      )
      
      plot_data <- data |>
        filter(
          is.finite(progression_value),
          is.finite(plotted_expression),
          !is.na(.data[[settings$facet_variable]]),
          !is.na(.data[[settings$color_variable]])
        )
      
      validate(
        need(nrow(plot_data) >= 2, "Too few observations remain."),
        need(
          n_distinct(plot_data$progression_value) >= 2,
          "At least two ages are required."
        )
      )
      
      plot <- ggplot(
        plot_data,
        aes(
          x = progression_value,
          y = plotted_expression,
          color = .data[[settings$color_variable]]
        )
      )
      
      if (settings$show_points) {
        plot <- plot + geom_point(alpha = 0.42, size = 1.0)
      }
      if (settings$smoother == "loess") {
        plot <- plot + geom_smooth(
          method = "loess",
          formula = y ~ x,
          se = FALSE,
          linewidth = 0.9,
          na.rm = TRUE
        )
      } else if (settings$smoother == "lm") {
        plot <- plot + geom_smooth(
          method = "lm",
          formula = y ~ x,
          se = FALSE,
          linewidth = 0.9,
          na.rm = TRUE
        )
      }
      
      plot +
        manual_color_scale(plot_data, settings$color_variable, "color") +
        scale_x_continuous(breaks = age_breaks, labels = age_labels) +
        facet_wrap(
          vars(.data[[settings$facet_variable]]),
          scales = "fixed",
          drop = TRUE
        ) +
        labs(
          x = "Developmental age",
          y = y_label,
          color = settings$color_variable
        ) +
        theme_minimal(base_size = 11) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold", size = 9),
          legend.position = "bottom"
        )
      
    } else if (settings$plot_type == "heatmap") {
      plot_data <- summarized_data()
      ggplot(
        plot_data,
        aes(
          x = .data[[settings$x_variable]],
          y = .data[[settings$second_dimension]],
          fill = expression
        )
      ) +
        geom_tile(color = "white", linewidth = 0.2) +
        scale_fill_viridis_c(option = "C", name = y_label) +
        labs(x = settings$x_variable, y = settings$second_dimension) +
        theme_minimal(base_size = 11) +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 55, hjust = 1)
        )
      
    } else if (settings$plot_type == "dot") {
      plot_data <- summarized_data()
      ggplot(
        plot_data,
        aes(
          x = .data[[settings$x_variable]],
          y = .data[[settings$second_dimension]],
          color = expression
        )
      ) +
        geom_point(size = 3.2, alpha = 0.9) +
        scale_color_viridis_c(option = "C", name = y_label) +
        labs(x = settings$x_variable, y = settings$second_dimension) +
        theme_minimal(base_size = 11) +
        theme(axis.text.x = element_text(angle = 55, hjust = 1))
      
    } else {
      data <- factor_field(data, settings$x_variable)
      data <- factor_field(data, settings$second_dimension)
      plot_data <- data |>
        filter(
          !is.na(.data[[settings$x_variable]]),
          !is.na(.data[[settings$second_dimension]]),
          is.finite(plotted_expression)
        ) |>
        group_by(across(all_of(c(
          settings$x_variable,
          settings$second_dimension
        )))) |>
        filter(n() >= 2) |>
        ungroup()
      
      validate(need(
        nrow(plot_data) >= 2,
        "No violin group has at least two observations."
      ))
      
      ggplot(
        plot_data,
        aes(
          x = .data[[settings$x_variable]],
          y = plotted_expression,
          fill = .data[[settings$x_variable]]
        )
      ) +
        geom_violin(
          scale = "width",
          trim = TRUE,
          color = "#214E68",
          alpha = 0.75,
          na.rm = TRUE
        ) +
        geom_jitter(
          width = 0.10,
          color = "black",
          alpha = 0.20,
          size = 0.35
        ) +
        stat_summary(
          fun.data = mean_se,
          geom = "errorbar",
          width = 0.18,
          color = "black"
        ) +
        stat_summary(
          fun = mean,
          geom = "point",
          shape = 21,
          fill = "white",
          color = "black",
          size = 2
        ) +
        manual_color_scale(plot_data, settings$x_variable, "fill") +
        facet_wrap(
          vars(.data[[settings$second_dimension]]),
          scales = "fixed",
          drop = TRUE
        ) +
        labs(x = settings$x_variable, y = y_label) +
        theme_minimal(base_size = 11) +
        theme(
          axis.text.x = element_text(angle = 55, hjust = 1, size = 7),
          legend.position = "none",
          strip.text = element_text(face = "bold", size = 9)
        )
    }
  }, res = 96, execOnResize = TRUE)
  
  output$plot_title <- renderText({
    if (is.null(loaded_gene())) {
      return("Retrieve a gene, then generate a plot")
    }
    settings <- plot_settings()
    zero_note <- if (settings$omit_zero_values) {
      "; zero-count observations removed"
    } else {
      ""
    }
    
    if (settings$plot_type == "trajectory") {
      paste0(
        loaded_gene(),
        " across developmental age by ",
        settings$facet_variable,
        " (colored by ",
        settings$color_variable,
        ")",
        zero_note
      )
    } else {
      paste0(
        loaded_gene(),
        " expression by ",
        settings$x_variable,
        " and ",
        settings$second_dimension,
        zero_note
      )
    }
  })
}

shinyApp(ui = ui, server = server)
