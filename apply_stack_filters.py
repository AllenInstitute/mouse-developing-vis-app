#!/usr/bin/env python3
from pathlib import Path
import sys

source = Path(sys.argv[1] if len(sys.argv) > 1 else "app.R")
target = Path(sys.argv[2] if len(sys.argv) > 2 else "app_stacked_filters.R")
s = source.read_text(encoding="utf-8")

def replace_once(old, new, label):
    global s
    count = s.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected one matching block, found {count}")
    s = s.replace(old, new, 1)

replace_once(
'''filter_fields <- setdiff(
  categorical_fields,
  omit_filter_fields
)
cell_type_fields <- unique(na.omit(c(
''',
'''numeric_filter_fields <- names(metadata)[vapply(
  metadata,
  is.numeric,
  logical(1)
)]
numeric_filter_fields <- setdiff(
  numeric_filter_fields,
  c(
    "CPM_scaling_factor"
  )
)
filter_fields <- setdiff(
  unique(c(
    categorical_fields,
    numeric_filter_fields
  )),
  omit_filter_fields
)

filter_input_id <- function(field) {
  field_index <- match(field, filter_fields)
  if (is.na(field_index)) {
    stop("Unknown filter field: ", field)
  }
  paste0("stack_filter_", field_index)
}

apply_filter_specification <- function(data, filters) {
  if (length(filters) == 0) {
    return(data)
  }
  
  for (specification in filters) {
    field <- specification$field
    if (!field %in% names(data)) {
      next
    }
    
    if (identical(specification$type, "numeric")) {
      bounds <- as.numeric(specification$value)
      if (
        length(bounds) == 2 &&
        all(is.finite(bounds))
      ) {
        numeric_values <- suppressWarnings(
          as.numeric(data[[field]])
        )
        data <- data[
          is.finite(numeric_values) &
            numeric_values >= min(bounds) &
            numeric_values <= max(bounds),
          ,
          drop = FALSE
        ]
      }
    } else {
      included_values <- as.character(specification$value)
      if (length(included_values) > 0) {
        data <- data[
          as.character(data[[field]]) %in% included_values,
          ,
          drop = FALSE
        ]
      }
    }
  }
  
  data
}

cell_type_fields <- unique(na.omit(c(
''',
"global filter helpers"
)

replace_once(
'''      #filter_values + .selectize-control .selectize-input,
      #filter_values + .selectize-control .selectize-dropdown,
      #filter_values-selectized + .selectize-dropdown {
        background: white !important; color: #111827 !important;
      }
''',
'''      #filter_fields + .selectize-control .selectize-input,
      #filter_fields + .selectize-control .selectize-dropdown,
      #filter_fields-selectized + .selectize-dropdown,
      .stacked-filter-control .selectize-input,
      .stacked-filter-control .selectize-dropdown {
        background: white !important; color: #111827 !important;
      }
      .stacked-filter-panel {
        border: 1px solid rgba(255,255,255,.28);
        border-radius: 5px;
        padding: 5px 6px 2px;
        margin: 3px 0 5px;
        background: rgba(255,255,255,.06);
      }
      .stacked-filter-control {
        margin-bottom: 5px;
      }
      .stacked-filter-control:last-child {
        margin-bottom: 1px;
      }
      .stacked-filter-control .irs--shiny .irs-bar,
      .stacked-filter-control .irs--shiny .irs-single,
      .stacked-filter-control .irs--shiny .irs-from,
      .stacked-filter-control .irs--shiny .irs-to {
        background: #4B9B58;
        border-color: #4B9B58;
      }
''',
"stack filter CSS"
)

replace_once(
'''      div(
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
''',
'''      selectizeInput(
        "filter_fields",
        "Filter metadata",
        choices = filter_fields,
        selected = character(),
        multiple = TRUE,
        options = list(
          placeholder = "Add one or more filters",
          maxOptions = 1000,
          closeAfterSelect = TRUE,
          plugins = list("remove_button")
        )
      ),
      div(
        class = "stacked-filter-panel",
        uiOutput("stacked_filter_controls")
      ),
''',
"stack filter UI"
)

start = s.index('''  observeEvent(input$filter_field, {\n''')
end = s.index('''  observeEvent(control_metadata(), {\n''', start)
old = s[start:end]
new = '''  output$stacked_filter_controls <- renderUI({
    selected_fields <- input$filter_fields
    if (is.null(selected_fields) || length(selected_fields) == 0) {
      return(
        div(
          class = "small text-white-50",
          "No metadata filters selected."
        )
      )
    }
    
    selected_fields <- selected_fields[
      selected_fields %in% filter_fields
    ]
    
    tagList(lapply(selected_fields, function(field) {
      input_id <- filter_input_id(field)
      current_value <- isolate(input[[input_id]])
      
      if (field %in% numeric_filter_fields) {
        numeric_values <- suppressWarnings(
          as.numeric(metadata[[field]])
        )
        numeric_values <- numeric_values[is.finite(numeric_values)]
        
        if (length(numeric_values) == 0) {
          return(NULL)
        }
        
        minimum <- min(numeric_values)
        maximum <- max(numeric_values)
        selected_range <- if (
          length(current_value) == 2 &&
          all(is.finite(as.numeric(current_value)))
        ) {
          pmax(
            minimum,
            pmin(maximum, as.numeric(current_value))
          )
        } else {
          c(minimum, maximum)
        }
        
        step_size <- if (
          all(abs(numeric_values - round(numeric_values)) < 1e-9)
        ) {
          1
        } else {
          max((maximum - minimum) / 100, .Machine$double.eps)
        }
        
        div(
          class = "stacked-filter-control",
          sliderInput(
            input_id,
            field,
            min = minimum,
            max = maximum,
            value = selected_range,
            step = step_size,
            separator = ""
          )
        )
      } else {
        choices <- field_levels(metadata, field)
        selected_values <- as.character(current_value)
        selected_values <- selected_values[
          selected_values %in% choices
        ]
        
        div(
          class = "stacked-filter-control",
          selectizeInput(
            input_id,
            field,
            choices = choices,
            selected = selected_values,
            multiple = TRUE,
            options = list(
              placeholder = "Include values",
              maxOptions = 1000,
              closeAfterSelect = TRUE,
              plugins = list("remove_button")
            )
          )
        )
      }
    }))
  })
  
  active_filter_specification <- reactive({
    selected_fields <- input$filter_fields
    if (is.null(selected_fields) || length(selected_fields) == 0) {
      return(list())
    }
    
    selected_fields <- selected_fields[
      selected_fields %in% filter_fields
    ]
    
    filters <- lapply(selected_fields, function(field) {
      value <- input[[filter_input_id(field)]]
      
      if (field %in% numeric_filter_fields) {
        if (
          is.null(value) ||
          length(value) != 2 ||
          !all(is.finite(as.numeric(value)))
        ) {
          return(NULL)
        }
        list(
          field = field,
          type = "numeric",
          value = as.numeric(value)
        )
      } else {
        value <- as.character(value)
        value <- value[!is.na(value) & nzchar(value)]
        if (length(value) == 0) {
          return(NULL)
        }
        list(
          field = field,
          type = "categorical",
          value = value
        )
      }
    })
    
    Filter(Negate(is.null), filters)
  })
  
  control_metadata <- reactive({
    filtered <- apply_filter_specification(
      metadata,
      active_filter_specification()
    )
    
    if (nrow(filtered) == 0) {
      metadata[0, , drop = FALSE]
    } else {
      filtered
    }
  })
  
'''
s = s[:start] + new + s[end:]

replace_once(
'''  observeEvent(input$reset_defaults, {
    updateSelectInput(session, "filter_field", selected = "none")
    updateSelectizeInput(
      session,
      "filter_values",
      choices = character(),
      selected = character(),
      server = TRUE
    )
''',
'''  observeEvent(input$reset_defaults, {
    updateSelectizeInput(
      session,
      "filter_fields",
      selected = character(),
      server = TRUE
    )
''',
"reset filters"
)

replace_once(
'''      plot_type = input$plot_type,
      filter_field = input$filter_field,
      filter_values = input$filter_values,
      omit_zero_values = isTRUE(input$omit_zero_values),
''',
'''      plot_type = input$plot_type,
      filters = active_filter_specification(),
      omit_zero_values = isTRUE(input$omit_zero_values),
''',
"plot settings filters"
)

replace_once(
'''    if (
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
''',
'''    data <- apply_filter_specification(
      data,
      settings$filters
    )
''',
"apply stacked filters"
)

target.write_text(s, encoding="utf-8")
print(f"Wrote {target} from {source}")
