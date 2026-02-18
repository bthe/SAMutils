#' Convert SAM Input Object to Table-Friendly Data Frame
#'
#' @param x Object read by `stockassessment::read.ices()`.
#'
#' @return A data frame suitable for display in dashboard tables.
#' @keywords internal
sam_dashboard_as_table <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    out <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
    rn <- rownames(out)
    if (!is.null(rn) && any(nzchar(rn))) {
      out <- cbind(year = rn, out)
      rownames(out) <- NULL
    }
    return(out)
  }

  if (is.list(x)) {
    item_names <- names(x)
    if (is.null(item_names)) {
      item_names <- paste0("item_", seq_along(x))
    }
    preview <- vapply(x, function(el) {
      if (length(el) == 0) {
        return("")
      }
      paste(utils::head(as.vector(el), 5), collapse = ", ")
    }, FUN.VALUE = character(1))
    n <- vapply(x, length, FUN.VALUE = integer(1))
    return(data.frame(field = item_names, length = n, preview = preview, stringsAsFactors = FALSE))
  }

  if (is.vector(x)) {
    return(data.frame(value = x, stringsAsFactors = FALSE))
  }

  data.frame(value = as.character(x), stringsAsFactors = FALSE)
}

#' Read NSHER SAM Input Files
#'
#' @param model_dir Directory containing NSHER `.dat` input files.
#'
#' @return A named list of raw ICES inputs.
#' @keywords internal
sam_dashboard_read_inputs <- function(model_dir) {
  dat_files <- sort(list.files(model_dir, pattern = "\\.dat$", full.names = TRUE))
  if (length(dat_files) == 0) {
    stop("No .dat files found in `model_dir`.")
  }

  out <- lapply(dat_files, stockassessment::read.ices)
  names(out) <- tools::file_path_sans_ext(basename(dat_files))
  out
}

#' Prepare NSHER Model Bundle for Dashboard
#'
#' @param model_dir Directory containing NSHER inputs (defaults to `testmore/nsher`).
#' @param retro_year Number of retrospective peels.
#'
#' @return A list with raw inputs, table inputs, fit bundle, and standard outputs.
#' @export
sam_dashboard_prepare <- function(model_dir = "testmore/nsher", retro_year = 5) {
  raw_inputs <- sam_dashboard_read_inputs(model_dir)

  needed <- c("cn", "cw", "dw", "lf", "lw", "mo", "nm", "pf", "pm", "sw", "survey")
  missing <- setdiff(needed, names(raw_inputs))
  if (length(missing) > 0) {
    stop("Missing required NSHER input files: ", paste(missing, collapse = ", "))
  }

  dat <- stockassessment::setup.sam.data(
    surveys = raw_inputs$survey,
    residual.fleets = raw_inputs$cn,
    prop.mature = raw_inputs$mo,
    stock.mean.weight = raw_inputs$sw,
    catch.mean.weight = raw_inputs$cw,
    dis.mean.weight = raw_inputs$dw,
    land.mean.weight = raw_inputs$lw,
    prop.f = raw_inputs$pf,
    prop.m = raw_inputs$pm,
    natural.mortality = raw_inputs$nm,
    land.frac = raw_inputs$lf
  )

  conf <- stockassessment::defcon(dat)
  conf$fbarRange <- c(2, 6)
  conf$corFlag <- 1
  conf$keyLogFpar <- matrix(
    c(
      -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, 0, 1, 2, 3, 4, 5, 6, -1,
      -1, 7, -1, -1, -1, -1, -1, -1, -1,
      8, -1, -1, -1, -1, -1, -1, -1, -1
    ),
    nrow = 4,
    byrow = TRUE
  )

  res <- {
    out <- NULL
    utils::capture.output({
      out <- full_sam_fit(dat, conf)
    })
    out
  }

  if (!is.null(retro_year) && !inherits(res$fit, "try-error")) {
    res$retro <- {
      out <- NULL
      utils::capture.output({
        out <- try(stockassessment::retro(res$fit, year = retro_year, ncores = 1))
      })
      out
    }
  }

  model_dat <- as.data.frame.table(raw_inputs$sw, stringsAsFactors = FALSE, responseName = "stock_weight") |>
    dplyr::rename(year = Var1, age = Var2) |>
    dplyr::mutate(
      year = suppressWarnings(as.numeric(as.character(year))),
      age = suppressWarnings(as.numeric(as.character(age)))
    ) |>
    tibble::as_tibble()

  std_out <- rby.sam(res$fit) |>
    dplyr::filter(variable %in% c("fbar", "ssb", "rec", "catch"))

  list(
    raw_inputs = raw_inputs,
    input_tables = lapply(raw_inputs, sam_dashboard_as_table),
    dat = dat,
    conf = conf,
    par = res$par,
    model_dat = model_dat,
    res = res,
    standard_output = std_out
  )
}

#' Build Input Tables from a Fitted SAM Object
#'
#' @param fit A fitted SAM object (`sam_fit$fit` from `full_sam_fit`).
#'
#' @return A named list of data frames for dashboard input tabs.
#' @keywords internal
sam_dashboard_input_tables_from_fit <- function(fit) {
  dat <- fit$data
  keep <- c(
    "stockMeanWeight", "catchMeanWeight", "natMor", "propMat",
    "propF", "propM", "landFrac", "disMeanWeight", "landMeanWeight"
  )
  keep <- keep[keep %in% names(dat)]
  out <- lapply(dat[keep], sam_dashboard_as_table)
  names(out) <- keep
  out
}

#' Convert `full_sam_fit` Output to a Dashboard Bundle
#'
#' @param sam_fit Output object from [full_sam_fit()].
#' @param input_tables Optional named list of input tables for the Input Data tab.
#' @param model_dat Optional long data frame with `year`, `age`, `stock_weight`.
#'
#' @return A dashboard bundle consumed by [dashboard_app()].
#' @keywords internal
sam_dashboard_bundle <- function(sam_fit, input_tables = NULL, model_dat = NULL) {
  if (is.null(sam_fit$fit)) {
    stop("`sam_fit` must be output from `full_sam_fit()` and include `$fit`.")
  }

  if (is.null(model_dat)) {
    model_dat <- as.data.frame.table(sam_fit$fit$data$stockMeanWeight, stringsAsFactors = FALSE, responseName = "stock_weight") |>
      dplyr::rename(year = Var1, age = Var2) |>
      dplyr::mutate(
        year = suppressWarnings(as.numeric(as.character(year))),
        age = suppressWarnings(as.numeric(as.character(age)))
      ) |>
      tibble::as_tibble()
  }

  if (is.null(input_tables)) {
    input_tables <- sam_dashboard_input_tables_from_fit(sam_fit$fit)
  }

  std_out <- rby.sam(sam_fit$fit) |>
    dplyr::filter(variable %in% c("fbar", "ssb", "rec", "catch"))

  list(
    input_tables = input_tables,
    model_dat = model_dat,
    res = sam_fit,
    standard_output = std_out
  )
}

#' Launch Dashboard from `full_sam_fit` Output
#'
#' @param sam_fit Output object from [full_sam_fit()].
#' @param input_tables Optional named list of input tables for the Input Data tab.
#' @param model_dat Optional long data frame with `year`, `age`, `stock_weight`.
#' @param title Dashboard title.
#'
#' @return A Shiny app object.
#' @export
dashboard_app <- function(
  sam_fit,
  input_tables = NULL,
  model_dat = NULL,
  title = "SAMutils Dashboard"
) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required. Install with install.packages('shiny').")
  }

  bundle <- sam_dashboard_bundle(
    sam_fit = sam_fit,
    input_tables = input_tables,
    model_dat = model_dat
  )
  input_ids <- make.names(names(bundle$input_tables), unique = TRUE)

  input_tabs <- lapply(seq_along(input_ids), function(i) {
    shiny::tabPanel(
      title = names(bundle$input_tables)[i],
      shiny::tableOutput(outputId = paste0("table_", input_ids[i]))
    )
  })

  input_data_panel <- if (length(input_tabs) == 0) {
    shiny::tabPanel(
      title = "Input Data",
      shiny::h4("No input tables supplied")
    )
  } else {
    shiny::tabPanel(
      title = "Input Data",
      do.call(shiny::tabsetPanel, input_tabs)
    )
  }

  ui <- shiny::navbarPage(
    title = title,
    input_data_panel,
    shiny::tabPanel(
      title = "Model Fit",
      shiny::fluidRow(
        shiny::column(6, shiny::h4("Observation residuals"), shiny::plotOutput("obs_resid", height = "600px")),
        shiny::column(6, shiny::h4("Process residuals"), shiny::plotOutput("proc_resid", height = "600px"))
      )
    ),
    shiny::tabPanel(
      title = "Standard Output",
      shiny::plotOutput("standard_plot", height = "700px"),
      shiny::tableOutput("standard_table")
    ),
    shiny::tabPanel(
      title = "Combined Fit",
      shiny::plotOutput("combined_fit_plot", height = "700px")
    ),
    shiny::tabPanel(
      title = "Leave-Out",
      shiny::plotOutput("leaveout_plot", height = "700px")
    ),
    shiny::tabPanel(
      title = "Likelihood",
      shiny::plotOutput("lik_plot", height = "900px")
    ),
    shiny::tabPanel(
      title = "Selectivity",
      shiny::plotOutput("selectivity_plot", height = "900px")
    ),
    shiny::tabPanel(
      title = "Stock-Recruit",
      shiny::plotOutput("sr_plot", height = "700px")
    ),
    shiny::tabPanel(
      title = "Retro",
      shiny::plotOutput("retro_plot", height = "700px")
    )
  )

  server <- function(input, output, session) {
    for (i in seq_along(input_ids)) {
      local({
        j <- i
        out_id <- paste0("table_", input_ids[j])
        tab_name <- names(bundle$input_tables)[j]
        output[[out_id]] <- shiny::renderTable({
          bundle$input_tables[[tab_name]]
        }, rownames = FALSE)
      })
    }

    output$obs_resid <- shiny::renderPlot({
      model_resid_plot(bundle$res)
    })

    output$proc_resid <- shiny::renderPlot({
      model_pres_resid_plot(bundle$res)
    })

    output$standard_plot <- shiny::renderPlot({
      bundle$standard_output |>
        ggplot2::ggplot(ggplot2::aes(year, median)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "steelblue") +
        ggplot2::geom_line(linewidth = 0.5) +
        ggplot2::facet_wrap(~variable, scales = "free_y") +
        ggplot2::labs(x = "Year", y = "Estimate")
    })

    output$standard_table <- shiny::renderTable({
      bundle$standard_output |>
        dplyr::select(year, variable, median, lower, upper) |>
        tidyr::pivot_wider(
          names_from = variable,
          values_from = c(median, lower, upper),
          names_glue = "{variable}_{.value}"
        ) |>
        dplyr::arrange(year)
    }, rownames = FALSE)

    output$combined_fit_plot <- shiny::renderPlot({
      model_combfit(bundle$res, bundle$model_dat)
    })

    output$leaveout_plot <- shiny::renderPlot({
      if (inherits(bundle$res$lo, "try-error")) {
        return(
          ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0, y = 0, label = "Leave-out failed to compute") +
            ggplot2::theme_void()
        )
      }
      model_lo_plot(bundle$res)
    })

    output$lik_plot <- shiny::renderPlot({
      model_lik_plot(bundle$res)
    })

    output$selectivity_plot <- shiny::renderPlot({
      model_selectivity_plot(bundle$res)
    })

    output$sr_plot <- shiny::renderPlot({
      SR_plot(bundle$res)
    })

    output$retro_plot <- shiny::renderPlot({
      if (inherits(bundle$res$retro, "try-error")) {
        return(
          ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0, y = 0, label = "Retro failed to compute") +
            ggplot2::theme_void()
        )
      }
      model_retro_plot(bundle$res)
    })
  }

  shiny::shinyApp(ui = ui, server = server)
}

#' Build SAM Dashboard App
#'
#' @param model_dir Directory containing NSHER inputs.
#' @param retro_year Number of retrospective peels.
#'
#' @return A Shiny app object.
#' @export
sam_dashboard_app <- function(model_dir = "testmore/nsher", retro_year = 5) {
  bundle <- sam_dashboard_prepare(model_dir = model_dir, retro_year = retro_year)
  dashboard_app(
    sam_fit = bundle$res,
    input_tables = bundle$input_tables,
    model_dat = bundle$model_dat
  )
}

#' Run SAM Dashboard
#'
#' @param model_dir Directory containing NSHER inputs.
#' @param retro_year Number of retrospective peels.
#' @param launch.browser Passed to `shiny::runApp()`.
#'
#' @return Invisibly returns the Shiny app run result.
#' @export
run_sam_dashboard <- function(
  model_dir = "testmore/nsher",
  retro_year = 5,
  launch.browser = interactive()
) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required. Install with install.packages('shiny').")
  }

  shiny::runApp(
    app = sam_dashboard_app(model_dir = model_dir, retro_year = retro_year),
    launch.browser = launch.browser
  )
}
