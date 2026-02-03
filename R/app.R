#' csbCLUST: Automated cluster detection from brightfield microscopy images
#'
#' The csbCLUST R package provides an interactive Shiny application
#' and a robust image-processing backend for automated detection,
#' quantification, and export of clusters from brightfield microscopy images.
#'
#' @keywords internal
"_PACKAGE"

#' Run the csbCLUST application
#'
#' @description
#' Launches the csbCLUST Shiny application.
#'
#' @export
runCsbCLUST <- function() {
  shiny::shinyApp(
    ui = ui,
    server = server
  )
}


# Imports ----------------------------------------------------------------
#' @import shiny
#' @import shinyjs
#' @import DT
#' @import fontawesome
#' @import colourpicker
#' @import shinyFiles
#' @import tibble
#' @import EBImage
#' @import dplyr

#' @importFrom utils write.csv
#' @importFrom tools file_path_sans_ext
#' @importFrom stats sd median

NULL

csbCLUST_section4_detect_clusters <- function(metadata, um_per_pixel = 0.65, progress_cb = NULL) {
  ## -------- GLOBAL SIZE CONSTRAINTS --------
  min_diameter_um <- 30
  max_diameter_um <- 100

  ## -------- CENTRAL REGION DETECTION PARAMETERS --------
  central_radius_frac        <- 0.65
  central_bg_stat            <- "mean"
  central_bg_offset          <- 0.05
  central_seed_brush_radius  <- 60

  ## -------- BORDER AND CORNER EXCLUSION PARAMETERS --------
  border_edge_frac           <- 0.98
  border_dark_mad_multiplier <- 4

  corner_edge_frac           <- 0.90
  corner_rim_thickness_frac  <- 0.10
  corner_min_dark_fraction   <- 0.05
  corner_max_core_white_frac <- 0.40
  corner_min_intensity_sd    <- 0.30

  ## -------- OUTPUT CONTAINERS (LIGHTWEIGHT ONLY) --------
  cluster_summary  <- vector("list", nrow(metadata))
  cluster_geometry <- vector("list", nrow(metadata))

  total_steps <- nrow(metadata) * 5
  current_step <- 0

  for (i in seq_len(nrow(metadata))) {

    message(
      "---- Processing image ", i, " / ", nrow(metadata),
      " : ", metadata$file_name[i]
    )

    ## -------- LOAD IMAGE --------
    img <- suppressWarnings(
      suppressWarnings(
        readImage(metadata$file_path[i])
      )
    )

    if (!is.null(progress_cb)) {
      current_step <<- current_step + 1
      progress_cb(current_step)
    }

    d  <- dim(img)
    w  <- d[1]
    h  <- d[2]
    nc <- if (length(d) == 3) d[3] else 1

    ## -------- GRAYSCALE CONVERSION AND THRESHOLDING --------
    img_gray <- if (nc == 3) channel(img, "gray") else img
    img_gray <- normalize(img_gray)

    thr    <- otsu(img_gray)
    img_bw <- (img_gray > thr)

    ## -------- WELL DETECTION --------
    ds <- 6

    small <- resize(
      img_bw,
      w = round(w / ds),
      h = round(h / ds)
    )

    small <- fillHull(small)
    small <- closing(small, suppressWarnings(makeBrush(25, "disc")))

    lab_small <- bwlabel(small)

    if (max(lab_small) == 0) {
      warning("No well detected in ", metadata$file_name[i])
      well_mask <- NULL
    } else {
      id <- which.max(computeFeatures.shape(lab_small)[,"s.area"])

      well_mask <- resize(
        lab_small == id,
        w = w,
        h = h
      ) > 0.5
    }

    ## -------- WELL MASK REFINEMENT --------
    if (!is.null(well_mask)) {
      well_mask <- opening(well_mask, suppressWarnings(makeBrush(6,  "disc")))
      well_mask <- closing(well_mask, suppressWarnings(makeBrush(10, "disc")))
    }

    ## -------- BACKGROUND NORMALIZATION (OUTSIDE WELL) --------
    img_bw_clean <- img_bw
    if (!is.null(well_mask)) {
      img_bw_clean[!well_mask] <- 1
    }

    ## -------- DARK PIXEL MASK AND DENSITY MAP --------
    dark_mask <- !img_bw_clean

    density <- filter2(
      dark_mask,
      suppressWarnings(makeBrush(15, "disc")),
      boundary = "replicate"
    )

    ## -------- CLUSTER PASS 1: STRICT --------
    seed_mask_strict <- density > 0.6 * max(density)
    seed_mask_strict <- opening(seed_mask_strict, suppressWarnings(makeBrush(5, "disc")))

    seed_lab_strict <- bwlabel(seed_mask_strict)
    message("  Strict pass seeds: ", max(seed_lab_strict))

    if (!is.null(progress_cb)) {
      current_step <<- current_step + 1
      progress_cb(current_step)
    }

    lab_strict <- propagate(density, seed_lab_strict, dark_mask)

    shape_strict  <- computeFeatures.shape(lab_strict)
    moment_strict <- computeFeatures.moment(lab_strict)

    labels_strict <- as.integer(rownames(shape_strict))
    diam_um_s <- 2 * sqrt(shape_strict[, "s.area"] / pi) * um_per_pixel

    labels_strict <- labels_strict[
      diam_um_s > min_diameter_um &
        diam_um_s < max_diameter_um
    ]

    ## -------- CLUSTER PASS 2: LOOSE --------
    seed_mask_loose <- density > 0.45 * max(density)
    seed_mask_loose <- opening(seed_mask_loose, suppressWarnings(makeBrush(5, "disc")))

    seed_lab_loose <- bwlabel(seed_mask_loose)
    message("  Loose pass seeds: ", max(seed_lab_loose))

    if (!is.null(progress_cb)) {
      current_step <<- current_step + 1
      progress_cb(current_step)
    }

    lab_loose <- propagate(density, seed_lab_loose, dark_mask)

    shape_loose  <- computeFeatures.shape(lab_loose)
    moment_loose <- computeFeatures.moment(lab_loose)

    labels_loose <- as.integer(rownames(shape_loose))
    diam_um_l <- 2 * sqrt(shape_loose[, "s.area"] / pi) * um_per_pixel

    labels_loose <- labels_loose[
      diam_um_l > min_diameter_um &
        diam_um_l < max_diameter_um
    ]

    ## -------- CLUSTER PASS 3: CENTRAL REGION --------
    cx_well <- w / 2
    cy_well <- h / 2
    well_radius <- min(w, h) / 2

    xx <- matrix(rep(seq_len(w), h), nrow = w)
    yy <- matrix(rep(seq_len(h), each = w), nrow = w)

    central_mask <- (xx - cx_well)^2 + (yy - cy_well)^2 <=
      (central_radius_frac * well_radius)^2

    bg_val <- if (central_bg_stat == "median")
      median(dark_mask[central_mask]) else mean(dark_mask[central_mask])

    seed_mask_central <- density > (bg_val + central_bg_offset)
    seed_mask_central <- seed_mask_central & central_mask
    seed_mask_central <- opening(
      seed_mask_central,
      suppressWarnings(makeBrush(central_seed_brush_radius, "disc"))
    )

    seed_lab_central <- bwlabel(seed_mask_central)
    message("  Central pass seeds: ", max(seed_lab_central))

    if (!is.null(progress_cb)) {
      current_step <<- current_step + 1
      progress_cb(current_step)
    }

    labels_central <- integer(0)

    if (max(seed_lab_central) > 0) {

      lab_central <- propagate(
        density,
        seed_lab_central,
        dark_mask & central_mask
      )

      shape_central  <- computeFeatures.shape(lab_central)
      moment_central <- computeFeatures.moment(lab_central)

      labels_central <- as.integer(rownames(shape_central))
      diam_um_c <- 2 * sqrt(shape_central[, "s.area"] / pi) * um_per_pixel

      labels_central <- labels_central[
        diam_um_c > min_diameter_um &
          diam_um_c < max_diameter_um
      ]
    }

    ## -------- MERGE CLUSTERS FROM ALL PASSES --------
    clusters <- list()

    add_clusters <- function(labels, shape, moment) {
      for (id in labels) {
        clusters[[length(clusters) + 1]] <<- list(
          cx = moment[as.character(id), "m.cx"],
          cy = moment[as.character(id), "m.cy"],
          r  = sqrt(shape[as.character(id), "s.area"] / pi)
        )
      }
    }

    add_clusters(labels_strict,  shape_strict,  moment_strict)
    add_clusters(labels_loose,   shape_loose,   moment_loose)
    add_clusters(labels_central, shape_central, moment_central)

    ## -------- EXCLUDE CORNER ARTIFACT CLUSTERS --------
    if (length(clusters) > 0) {

      keep_idx <- sapply(clusters, function(cl) {

        dist_edge <- sqrt((cl$cx - cx_well)^2 + (cl$cy - cy_well)^2)
        if (dist_edge < corner_edge_frac * well_radius)
          return(TRUE)

        r <- ceiling(cl$r)

        full_mask <- (xx - cl$cx)^2 + (yy - cl$cy)^2 <= r^2
        core_r <- ceiling(r * (1 - corner_rim_thickness_frac))
        core_mask <- (xx - cl$cx)^2 + (yy - cl$cy)^2 <= core_r^2
        rim_mask  <- full_mask & !core_mask

        rim_dark_frac   <- mean(img_bw_clean[rim_mask] == 0)
        core_white_frac <- mean(img_bw_clean[core_mask] == 1)
        intensity_sd    <- sd(img_bw_clean[full_mask])

        !(rim_dark_frac > corner_min_dark_fraction &
            core_white_frac > corner_max_core_white_frac &
            intensity_sd > corner_min_intensity_sd)
      })

      message("  Removed corner fake clusters: ", sum(!keep_idx))
      clusters <- clusters[keep_idx]
    }

    ## -------- RESOLVE OVERLAPPING CLUSTERS --------
    if (length(clusters) > 1) {

      keep_idx <- rep(TRUE, length(clusters))

      for (i1 in seq_along(clusters)) {
        if (!keep_idx[i1]) next
        for (i2 in seq_along(clusters)) {
          if (i2 <= i1 || !keep_idx[i2]) next

          dx <- clusters[[i1]]$cx - clusters[[i2]]$cx
          dy <- clusters[[i1]]$cy - clusters[[i2]]$cy
          dist <- sqrt(dx^2 + dy^2)

          if (dist < (clusters[[i1]]$r + clusters[[i2]]$r)) {
            if (clusters[[i1]]$r >= clusters[[i2]]$r) {
              keep_idx[i1] <- FALSE
            } else {
              keep_idx[i2] <- FALSE
            }
          }
        }
      }

      message("  Removed nested / overlapping clusters: ",
              sum(!keep_idx))

      clusters <- clusters[keep_idx]
    }

    message("  Final clusters to draw: ", length(clusters))

    if (!is.null(progress_cb)) {
      current_step <<- current_step + 1
      progress_cb(current_step)
    }

    ## -------- STORE CLUSTER METRICS AND GEOMETRY --------
    if (length(clusters) > 0) {

      areas_px <- sapply(clusters, function(cl) pi * cl$r^2)
      areas_um <- areas_px * (um_per_pixel^2)

      cluster_summary[[i]] <- tibble(
        file_name          = metadata$file_name[i],
        n_clusters         = length(clusters),
        mean_cluster_area  = mean(areas_um),
        total_cluster_area = sum(areas_um)
      )

      cluster_geometry[[i]] <- tibble(
        cx = sapply(clusters, `[[`, "cx"),
        cy = sapply(clusters, `[[`, "cy"),
        r  = sapply(clusters, `[[`, "r")
      )

    } else {

      cluster_summary[[i]] <- tibble(
        file_name          = metadata$file_name[i],
        n_clusters         = 0,
        mean_cluster_area  = NA_real_,
        total_cluster_area = 0
      )

      cluster_geometry[[i]] <- tibble()
    }

    rm(img, img_gray, img_bw, img_bw_clean, dark_mask, density,
       lab_strict, lab_loose, lab_central)
    gc(verbose = FALSE)
  }

  list(
    cluster_summary  = cluster_summary,
    cluster_geometry = cluster_geometry,
    metadata         = metadata
  )
}

css_tabs <- "
body { background-color: #ffffff; }
.nav-tabs { border-bottom: 1px solid #000000; }
.nav-tabs > li > a { color: #000000 !important; }
.nav-tabs > li.active > a {
  border: 1px solid #000000;
  border-bottom-color: transparent;
}
.btn-primary {
  background-color: #000000 !important;
  border-color: #000000 !important;
}
/* Disable ALL tabs by default */
.nav-tabs > li > a {
  pointer-events: none;
  cursor: default;
  color: #bdbdbd !important;  /* light grey text for disabled tabs */
}

/* Active tab: restore original text color + clickability */
.nav-tabs > li.active > a {
  pointer-events: auto;
  cursor: pointer;
  color: #000000 !important;  /* original active text colour */
}
"

ui <- fluidPage(
  useShinyjs(),

  tags$head(
    tags$style(HTML(css_tabs)),
    tags$style("#about {border-color:white; font-size:12px}"),
    tags$style(HTML("
  #input_dir,
  #export_dir {
    font-size: 13px !important;
    padding: 6px 12px !important;
    line-height: 1.2 !important;
  }

  #input_dir::before {
    font-family: 'Font Awesome 6 Free';
    font-weight: 900;
    content: '\\f07c'; /* folder-open */
    margin-right: 6px;
  }

  #export_dir::before {
    font-family: 'Font Awesome 6 Free';
    font-weight: 900;
    content: '\\f07b'; /* folder-tree */
    margin-right: 6px;
  }
"))
  ),

  headerPanel(
    fluidRow(
      column(
        width = 3,
        tags$div(
          style = "display:flex; align-items:center; gap:8px;",
          tags$img(
            src = base64enc::dataURI(
              file = system.file("icons", "Logo.png", package = "csbCLUST"),
              mime = "image/png"
            ),
            height = "50px"
          ),
          tags$span(
            "csbCLUST",
            style = "font-weight:600; font-size:30px;"
          )
        )
      ),
      column(
        width = 2,
        offset = 7,
        align = "right",
        actionButton(
          inputId = "about",
          label = "About csbCLUST"
        )
      )
    )
  ),

  tabsetPanel(
    id = "main_tabs",

    ## ---- TAB 1: Processing ----
    tabPanel(
      value = "Processing_tab",
      title = tagList(icon("cogs"), "Processing"),
      div(
        style = "height:70vh; display:flex; flex-direction:column; align-items:center; justify-content:center;",
        shinyDirButton(
          "input_dir",
          "Select input folder",
          "Select input folder"
        ),
        br(),
        actionButton("Processing_proceed", "Proceed", class = "btn-primary", disabled = TRUE)
      )
    ),


    ## ---- TAB 3: Clusters ----
    tabPanel(
      value = "clusters_tab",
      title = tagList(icon("circle-nodes"), "Clusters"),
      br(),
      fluidRow(
        column(
          3,
          div(style="border:1px solid #000; padding:10px;",
              selectInput(
                "cluster_image",
                tagList(icon("image"), "Select image"),
                choices = NULL
              ),
              checkboxInput(
                "add_scalebar",
                tagList(icon("ruler-horizontal"), "Add scale bar"),
                TRUE
              ),
              sliderInput("scalebar_offset", "Scale bar offset",
                          min = 0, max = 0.1, value = 0.05),
              sliderInput("scalebar_width", "Scale bar width (um)",
                          min = 50, max = 400, value = 200),
              actionButton("clusters_proceed", "Proceed", class="btn-primary")
          )
        ),
        column(
          6,
          div(
            style = "
      width:120%;
      height:100vh;
    ",
            imageOutput(
              "cluster_plot",
              width  = "100%",
              height = "100%"
            )
          )
        ),
        column(3)
      )
    ),

    ## ---- TAB 4: Table ----
    tabPanel(
      value = "table_tab",
      title = tagList(icon("table"), "Table"),
      br(),
      fluidRow(
        column(12, align="right",
               actionButton("table_proceed", "Proceed", class="btn-primary")),
        br(), br(),
        column(12, DTOutput("cluster_table"))
      )
    ),

    ## ---- TAB 5: Exports ----
    tabPanel(
      value = "exports_tab",
      title = tagList(icon("file-export"), "Exports"),
      div(
        style = "height:70vh; display:flex; flex-direction:column; align-items:center; justify-content:center;",
        shinyDirButton(
          "export_dir",
          "Select export folder",
          "Select export folder"
        ),
        br(),
        actionButton(
          "export_run",
          tagList("Export images and stats"),
          class = "btn-primary",
          disabled = TRUE
        )
      )
    )
  )
)

server <- function(input, output, session) {

  observeEvent(input$about, {
    showModal(
      modalDialog(
        title = "csbCLUST",
        "This is an open-source project that provides an automated pipeline for
      robust cluster detection and quantification from brightfield microscopy images.",
        br(),
        br(),
        strong("Project website:"),
        "https://github.com/BonilhaCaio",
        footer = NULL,
        size = "m",
        easyClose = TRUE
      )
    )
  })

  volumes <- c(Home = normalizePath("~"), Root = "/")

  shinyDirChoose(
    input,
    "input_dir",
    roots = volumes,
    session = session
  )

  shinyDirChoose(
    input,
    "export_dir",
    roots = volumes,
    session = session
  )

  rv <- reactiveValues(
    input_dir = NULL,
    export_dir = NULL,
    metadata = NULL,
    results = NULL
  )

  ## ---- Folder selection ----
  observeEvent(input$input_dir, {

    path <- parseDirPath(volumes, input$input_dir)

    if (length(path) == 0) {
      disable("Processing_proceed")
      return()
    }

    if (!dir.exists(path)) {
      disable("Processing_proceed")
      return()
    }

    rv$input_dir <- path
    enable("Processing_proceed")
  })

  observeEvent(input$Processing_proceed, {

    ## ---- DISABLE INPUTS WHILE RUNNING ----
    disable("input_dir")
    disable("Processing_proceed")

    ## ---- Build metadata FIRST (guaranteed) ----
    files <- list.files(rv$input_dir, pattern = "\\.tif(f)?$", full.names = TRUE)

    rv$metadata <- tibble(
      file_path = files,
      file_name = basename(files)
    )

    ## ---- Progress handling ----
    n <- nrow(rv$metadata)

    step_counter <- 0
    total_steps  <- n * 5

    withProgress(message = "Progress: ", value = 0, {

      rv$results <- csbCLUST_section4_detect_clusters(
        metadata = rv$metadata,
        um_per_pixel = 0.65,
        progress_cb = function(step) {

          step_counter <<- step_counter + 1

          pct <- round((step_counter / total_steps) * 100)

          incProgress(
            amount = 1 / total_steps,
            detail = paste0(pct, "%")
          )
        }
      )

    })

    ## ---- RE-ENABLE AFTER COMPLETION ----
    enable("input_dir")

    updateTabsetPanel(session, "main_tabs", "clusters_tab")

  })

  ## ---- Clusters tab ----
  observe({
    updateSelectInput(session,"cluster_image",
                      choices = rv$results$metadata$file_name)
  })

  output$cluster_plot <- renderImage({
    req(input$cluster_image, rv$results)

    idx <- which(rv$results$metadata$file_name == input$cluster_image)
    img <- readImage(rv$results$metadata$file_path[idx])

    geom <- rv$results$cluster_geometry[[idx]]

    d  <- dim(img)
    w0 <- d[1]
    h0 <- d[2]
    nc <- if (length(d) == 3) d[3] else 1

    ## ---- TEMPORARY DISPLAY DOWNSCALING ----
    max_dim <- 1200
    scale <- min(1, max_dim / max(w0, h0))

    if (scale < 1) {
      img <- resize(
        img,
        w = round(w0 * scale),
        h = round(h0 * scale)
      )
    }

    d <- dim(img)
    w <- d[1]
    h <- d[2]

    ## ---- DRAW CLUSTERS ----
    if (nrow(geom) > 0) {

      xx <- matrix(rep(seq_len(w), h), nrow = w)
      yy <- matrix(rep(seq_len(h), each = w), nrow = w)

      for (k in seq_len(nrow(geom))) {

        cx <- geom$cx[k] * scale
        cy <- geom$cy[k] * scale
        r  <- geom$r[k]  * scale

        ring <- abs((xx - cx)^2 + (yy - cy)^2 - r^2) <= (2 * r)

        if (nc == 1) {
          img[ring] <- 1
        } else {
          img[ring, 1] <- 1
          img[ring, 2] <- 0
          img[ring, 3] <- 0
        }
      }
    }

    ## ---- DRAW SCALE BAR ----
    if (isTRUE(input$add_scalebar)) {

      um_per_pixel <- 0.65
      px_per_um <- 1 / um_per_pixel

      scale_um <- input$scalebar_width
      scale_px <- round(scale_um * px_per_um * scale)

      bar_height <- round(40 * scale)
      margin_px  <- round(input$scalebar_offset * min(w, h))

      x_start <- w - margin_px - scale_px
      x_end   <- x_start + scale_px
      y_start <- h - margin_px
      y_end   <- y_start + bar_height

      x_start <- max(1, x_start)
      y_start <- max(1, y_start)
      x_end   <- min(w, x_end)
      y_end   <- min(h, y_end)

      if (nc == 1) {
        img[x_start:x_end, y_start:y_end] <- 1
      } else {
        img[x_start:x_end, y_start:y_end, 1:3] <- 1
      }
    }

    tmp <- tempfile(fileext = ".png")
    writeImage(img, tmp)

    list(
      src = tmp,
      contentType = "image/png",
      width = "100%",
      deleteFile = TRUE
    )
  }, deleteFile = TRUE)

  observeEvent(input$clusters_proceed,{
    updateTabsetPanel(session, "main_tabs", "table_tab")
  })

  ## ---- Table ----
  output$cluster_table <- renderDT({

    tbl <- dplyr::bind_rows(rv$results$cluster_summary)

    tbl <- tbl %>%
      dplyr::mutate(
        dplyr::across(where(is.numeric), ~ as.integer(round(.)))
      )

    datatable(
      tbl,
      rownames = FALSE,
      options = list(
        dom = "t",
        paging = FALSE,
        ordering = FALSE,
        searching = FALSE,
        info = FALSE
      )
    )
  })

  observeEvent(input$table_proceed,{
    updateTabsetPanel(session, "main_tabs", "exports_tab")
  })

  ## ---- Export ----

  observeEvent(input$export_dir, {

    path <- parseDirPath(volumes, input$export_dir)

    if (length(path) == 0) return()
    if (!dir.exists(path)) return()
    if (identical(path, rv$input_dir)) return()

    rv$export_dir <- path
    enable("export_run")
  })

  observeEvent(input$export_run, {

    req(rv$export_dir, rv$results)

    ## ---- DISABLE ALL ACTION BUTTONS DURING EXPORT ----
    shinyjs::disable(selector = "button")

    ## ---- SHOW WAITING MODAL ----
    showModal(
      modalDialog(
        title = "Export in progress",
        div(
          style = "font-size:15px;",
          icon("hourglass-half"),
          " Please wait while images and statistics are being exported.",
          br(), br(),
          "This may take a few moments depending on the number and size of images."
        ),
        footer = NULL,
        easyClose = FALSE
      )
    )

    ## ---- EXPORT SUMMARY TABLE ----
    summary_tbl <- dplyr::bind_rows(rv$results$cluster_summary)

    write.csv(
      summary_tbl,
      file = file.path(rv$export_dir, "cluster_summary.csv"),
      row.names = FALSE
    )

    ## ---- EXPORT IMAGES + PER-IMAGE CLUSTER GEOMETRY ----
    for (i in seq_len(nrow(rv$results$metadata))) {

      img <- readImage(rv$results$metadata$file_path[i])
      geom <- rv$results$cluster_geometry[[i]]

      d  <- dim(img)
      w  <- d[1]
      h  <- d[2]
      nc <- if (length(d) == 3) d[3] else 1

      if (nrow(geom) > 0) {

        xx <- matrix(rep(seq_len(w), h), nrow = w)
        yy <- matrix(rep(seq_len(h), each = w), nrow = w)

        for (k in seq_len(nrow(geom))) {

          cx <- geom$cx[k]
          cy <- geom$cy[k]
          r  <- geom$r[k]

          ring <- abs((xx - cx)^2 + (yy - cy)^2 - r^2) <= (2 * r)

          if (nc == 1) {
            img[ring] <- 1
          } else {
            img[ring, 1] <- 1
            img[ring, 2] <- 0
            img[ring, 3] <- 0
          }
        }
      }

      if (isTRUE(input$add_scalebar)) {

        um_per_pixel <- 0.65
        px_per_um <- 1 / um_per_pixel

        scale_um <- input$scalebar_width
        scale_px <- round(scale_um * px_per_um)

        bar_height <- 40
        margin_px  <- round(input$scalebar_offset * min(w, h))

        x_start <- w - margin_px - scale_px
        x_end   <- x_start + scale_px
        y_start <- h - margin_px
        y_end   <- y_start + bar_height

        x_start <- max(1, x_start)
        y_start <- max(1, y_start)
        x_end   <- min(w, x_end)
        y_end   <- min(h, y_end)

        if (nc == 1) {
          img[x_start:x_end, y_start:y_end] <- 1
        } else {
          img[x_start:x_end, y_start:y_end, 1:3] <- 1
        }
      }

      img_name <- paste0(
        tools::file_path_sans_ext(rv$results$metadata$file_name[i]),
        "_clusters.png"
      )

      writeImage(
        img,
        file.path(rv$export_dir, img_name)
      )

      if (nrow(geom) > 0) {

        geom_out <- geom |>
          dplyr::rename(
            `x coord` = cx,
            `y coord` = cy,
            r = r
          )

        csv_name <- paste0(
          tools::file_path_sans_ext(rv$results$metadata$file_name[i]),
          "_clusters.csv"
        )

        write.csv(
          geom_out,
          file = file.path(rv$export_dir, csv_name),
          row.names = FALSE
        )
      }
    }

    ## ---- CLOSE WAIT MODAL ----
    removeModal()

    ## ---- SHOW SUCCESS MODAL ----
    showModal(
      modalDialog(
        title = "Export completed",
        div(
          style = "font-size:15px;",
          icon("circle-check"),
          " All images and statistics were successfully exported.",
          br(), br(),
          "You may now safely close the application."
        ),
        easyClose = TRUE,
        footer = modalButton("OK")
      )
    )

    ## ---- RE-ENABLE BUTTONS ----
    shinyjs::enable(selector = "button")

  })

}


