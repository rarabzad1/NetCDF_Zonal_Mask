# app.R - Shiny wrapper for nc_spatial_mask() - output filename forced to NULL and improved variable selection
required_packages <- c(
  "shiny", "leaflet", "sf", "DT", "shinyWidgets", "zip", "shinyjs",
  "ncdf4", "geosphere", "dplyr", "sp", "lwgeom", "rmapshaper", "raster"
)

install_log <- ""
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install_log <- paste0(install_log, Sys.time(), " — Installing missing package: ", pkg, "\n")
    tryCatch({
      install.packages(pkg, repos = "https://cloud.r-project.org")
      install_log <- paste0(install_log, Sys.time(), " — ✅ Successfully installed: ", pkg, "\n")
    }, error = function(e) {
      install_log <- paste0(install_log, Sys.time(), " — ❌ Failed to install ", pkg, ": ", e$message, "\n")
    })
  }
}
if (install_log == "") install_log <- paste0(Sys.time(), " — ✅ All required packages are already installed.\n")

library(shiny)
library(leaflet)
library(sf)
library(DT)
library(shinyWidgets)
library(zip)
library(shinyjs)
library(ncdf4)
library(geosphere)
library(dplyr)
library(sp)
library(lwgeom)
library(rmapshaper)
library(raster)

options(shiny.maxRequestSize = 200 * 1024^2)  # 200 MB
initial_log <- install_log

# ----- Top-level helper: same heuristics as in your function to find spatio-temporal vars -----
find_spt_vars <- function(nc) {
  if (is.null(nc) || !"var" %in% names(nc)) stop("`nc` must be an already-opened ncdf4 object.")
  norm <- function(x) tolower(gsub("[^a-z0-9]", "", as.character(x)))
  time_words <- c("time","t","tim","date","day","hour","minute","second","record","step","tau")
  lat_words  <- c("lat","latitude","y","northing")
  lon_words  <- c("lon","longitude","x","easting")
  time_regex <- paste0("\\b(", paste(time_words, collapse="|"), ")\\b")
  lat_regex  <- paste0("\\b(", paste(lat_words, collapse="|"), ")\\b")
  lon_regex  <- paste0("\\b(", paste(lon_words, collapse="|"), ")\\b")
  units_indicate_time <- function(un) {
    if (is.null(un)) return(FALSE)
    un <- tolower(un)
    grepl("since", un) || grepl("\\b(sec|second|seconds|min|minute|minutes|hour|hr|day|doy|julian|year|month)\\b", un)
  }
  units_indicate_lat <- function(un) {
    if (is.null(un)) return(FALSE)
    un <- tolower(un)
    grepl("degree_north|degrees_north|deg_north|degree_n|degrees_n", un) ||
      (grepl("degree|deg", un) && grepl("north|n", un))
  }
  units_indicate_lon <- function(un) {
    if (is.null(un)) return(FALSE)
    un <- tolower(un)
    grepl("degree_east|degrees_east|deg_east|degree_e", un) ||
      (grepl("degree|deg", un) && grepl("east|e", un))
  }
  
  vars <- nc$var
  dims_info <- nc$dim
  found <- character(0)
  
  for (vname in names(vars)) {
    v <- vars[[vname]]
    dim_objs <- v$dim
    if (is.null(dim_objs) || length(dim_objs) == 0) next
    dim_names <- sapply(dim_objs, function(d) d$name)
    
    is_time <- is_lat <- is_lon <- logical(length(dim_names))
    
    # robustly fetch coordinates attribute (if present)
    coords_attr_raw <- try(ncatt_get(nc, vname, "coordinates"), silent = TRUE)
    coord_tokens <- NULL
    if (!inherits(coords_attr_raw, "try-error") && !is.null(coords_attr_raw)) {
      if (is.list(coords_attr_raw) && !is.null(coords_attr_raw$value) && is.character(coords_attr_raw$value)) {
        coord_tokens <- unlist(strsplit(coords_attr_raw$value, "\\s+"))
        coord_tokens <- tolower(coord_tokens[nzchar(coord_tokens)])
      } else if (is.character(coords_attr_raw) && length(coords_attr_raw) == 1) {
        coord_tokens <- unlist(strsplit(coords_attr_raw, "\\s+"))
        coord_tokens <- tolower(coord_tokens[nzchar(coord_tokens)])
      }
    }
    
    for (i in seq_along(dim_names)) {
      dn <- dim_names[i]
      if (is.na(dn) || nchar(dn) == 0) next
      dn_norm <- norm(dn)
      if (grepl(time_regex, dn_norm, perl = TRUE)) is_time[i] <- TRUE
      if (grepl(lat_regex, dn_norm, perl = TRUE))  is_lat[i]  <- TRUE
      if (grepl(lon_regex, dn_norm, perl = TRUE))  is_lon[i]  <- TRUE
      
      dimobj <- if (!is.null(dims_info[[dn]])) dims_info[[dn]] else NULL
      if (!is.null(dimobj) && !is.null(dimobj$units)) {
        if (units_indicate_time(dimobj$units)) is_time[i] <- TRUE
        if (units_indicate_lat(dimobj$units))  is_lat[i]  <- TRUE
        if (units_indicate_lon(dimobj$units))  is_lon[i]  <- TRUE
      }
      
      coordvar <- NULL
      if (!is.null(vars[[dn]])) coordvar <- vars[[dn]]
      if (!is.null(coordvar) && !is.null(coordvar$units)) {
        if (units_indicate_time(coordvar$units)) is_time[i] <- TRUE
        if (units_indicate_lat(coordvar$units))  is_lat[i]  <- TRUE
        if (units_indicate_lon(coordvar$units))  is_lon[i]  <- TRUE
      }
      
      if (!is.null(coord_tokens) && dn %in% coord_tokens) {
        if (any(grepl("lat", coord_tokens))) is_lat[i] <- TRUE
        if (any(grepl("lon", coord_tokens))) is_lon[i] <- TRUE
      }
      
      vals_try <- try(ncvar_get(nc, dn), silent = TRUE)
      if (!inherits(vals_try, "try-error") && is.numeric(vals_try)) {
        rng <- range(as.numeric(vals_try), na.rm = TRUE)
        if (is.finite(rng[1]) && is.finite(rng[2])) {
          if (rng[1] >= -90 && rng[2] <= 90) is_lat[i] <- TRUE
          if ((rng[1] >= -180 && rng[2] <= 180) || (rng[1] >= 0 && rng[2] <= 360)) is_lon[i] <- TRUE
        }
      }
    } # end dims loop
    
    n_time <- sum(is_time, na.rm = TRUE)
    n_lat  <- sum(is_lat, na.rm = TRUE)
    n_lon  <- sum(is_lon, na.rm = TRUE)
    
    if ((n_time >= 1) && ((n_lat + n_lon) >= 2)) found <- c(found, vname)
  } # end vars loop
  
  unique(found)
}

ui <- fluidPage(
  useShinyjs(),
  tags$div(
    style = "display: flex; align-items: center; gap: 15px; margin-bottom: 20px;",
    tags$img(src = "logo.png", width = "150px", style = "border-radius: 20px;"),
    tags$h2("NetCDF Zonal Mask", style = "margin: 0;"),
    tags$p(
      "The NetCDF mask Shiny App provides a user-friendly interface to spatially mask NetCDF datasets using a shapefile polygon. The app allows users to selectively process spatio-temporal variables in a NetCDF file, apply a polygonal mask from a shapefile, and generate a masked or aggregated NetCDF output. It also provides an interactive visualization of the resulting data, including time-averaged maps for each variable.",
      style = "font-size: 1.2em; color: #555;"
    )
  ),
  sidebarLayout(
    sidebarPanel(
      fileInput("ncfile", "Upload NetCDF file (.nc)", accept = ".nc"),
      fileInput("maskzip", "Upload HRU shapefile (zip containing .shp/.shx/.dbf/.prj)", accept = ".zip"),
      uiOutput("var_ui"),   # dynamic variable multiselect (spatio-temporal only)
      actionButton("generate", "Generate Masked netCDF", icon = icon("play")),
      checkboxInput("show_plot_tab", "Show Plot tab", value = TRUE),
      conditionalPanel(
        "output.has_nc == true",
        downloadButton("download_nc", "Download masked NetCDF", class = "btn-primary")
      ),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Log", verbatimTextOutput("log", placeholder = TRUE)),
        # Plot tab is only rendered on the client when the checkbox is checked
        conditionalPanel(
          condition = "input.show_plot_tab == true",
          tabPanel("Plot", value = "plot", uiOutput("plot_tabs_ui"))
        )
      ),
      width = 9
    )
  )
)

server <- function(input, output, session) {
  source("https://raw.githubusercontent.com/rarabzad/NetCDF_Zonal_Mask/refs/heads/main/nc_spatial_mask.R")
  logs <- reactiveVal(initial_log)
  append_log <- function(msg) {
    logs(paste0(logs(), Sys.time(), " — ", msg, "\n"))
  }
  output$log <- renderText({ logs() })
  output_nc_path <- reactiveVal(NULL)
  output$has_nc <- reactive({ !is.null(output_nc_path()) })
  outputOptions(output, "has_nc", suspendWhenHidden = FALSE)
  output$status <- renderText({
    if (is.null(input$ncfile)) return("Waiting for NetCDF upload...")
    paste0("NetCDF: ", input$ncfile$name)
  })

  # ---------- Helper functions for plotting ----------
  # find coordinate variable names in the nc object (lon & lat)
  find_coord_vars <- function(nc) {
    vars <- names(nc$var)
    norm <- function(x) tolower(gsub("[^a-z0-9]", "", as.character(x)))
    lon_words <- c("lon","longitude","x","easting")
    lat_words <- c("lat","latitude","y","northing")
    lon_name <- NULL
    lat_name <- NULL
    for (v in vars) {
      n <- norm(v)
      if (is.null(lon_name) && any(sapply(lon_words, function(w) grepl(w, n)))) lon_name <- v
      if (is.null(lat_name) && any(sapply(lat_words, function(w) grepl(w, n)))) lat_name <- v
    }
    list(lon = lon_name, lat = lat_name)
  }
  
  # detect time-dimension index for a variable (returns integer index or NULL)
  detect_time_dim_index <- function(nc, varname) {
    time_words <- c("time","t","tim","date","day","hour","minute","second","record","step","tau")
    time_regex <- paste0("\\b(", paste(time_words, collapse="|"), ")\\b")
    dimnames_var <- sapply(nc$var[[varname]]$dim, function(d) d$name)
    # check dim names first
    hit <- which(grepl(time_regex, tolower(dimnames_var), perl = TRUE))
    if (length(hit) >= 1) return(hit[1])
    # check dim units for "since" or time words
    for (i in seq_along(dimnames_var)) {
      dn <- dimnames_var[i]
      dimobj <- nc$dim[[dn]]
      if (!is.null(dimobj) && !is.null(dimobj$units)) {
        un <- tolower(dimobj$units)
        if (grepl("since", un) || grepl("\\b(sec|second|minute|min|hour|day|julian|year|month)\\b", un)) return(i)
      }
    }
    # fallback heuristics: if var has 3+ dims assume the last is time
    var_dimlen <- length(dimnames_var)
    if (var_dimlen >= 3) return(var_dimlen)
    return(NULL)
  }
  
  # ---------- Build Plot UI and plot outputs when a masked file is available ----------
  observeEvent(output_nc_path(), {
    req(output_nc_path())
    ncpath <- output_nc_path()
    # try to open; if fails, make UI empty
    nc <- tryCatch(nc_open(ncpath), error = function(e) { append_log(paste0("Error opening masked nc for plotting: ", e$message)); NULL })
    if (is.null(nc)) {
      output$plot_tabs_ui <- renderUI({ NULL })
      return()
    }
    
    # Try to identify spatio-temporal variables in the output using same heuristic
    vars_to_plot <- tryCatch(find_spt_vars(nc), error = function(e) { names(nc$var) })
    
    # as a fallback exclude coordinate vars (lon/lat) from plotting
    coords <- find_coord_vars(nc)
    vars_to_plot <- setdiff(vars_to_plot, c(coords$lon, coords$lat))
    if (length(vars_to_plot) == 0) {
      # fallback: use all non-dimension-only variables
      vars_to_plot <- setdiff(names(nc$var), c(coords$lon, coords$lat))
    }
    
    # Clip to reasonable number to avoid creating too many tabs (optional)
    # (you can remove this if you want all variables)
    # vars_to_plot <- vars_to_plot[1:min(length(vars_to_plot), 20)]
    
    # Create UI: a tab for each variable
    output$plot_tabs_ui <- renderUI({
      if (is.null(output_nc_path()) || !isTruthy(input$show_plot_tab)) return(NULL)
      tabs <- lapply(vars_to_plot, function(v) {
        plotid <- paste0("plot_", v)
        tabPanel(title = v,
                 # allow some vertical room for map/plot
                 plotOutput(plotid, height = "450px")
        )
      })
      do.call(tabsetPanel, c(tabs, id = "plot_var_tabs"))
    })
    
    # For each variable create a renderPlot that computes time-average and plots it
    for (v in vars_to_plot) {
      local({
        vv <- v
        pid <- paste0("plot_", vv)
        output[[pid]] <- renderPlot({
          req(output_nc_path())
          # open the netcdf inside the render (keeps lifetime correct)
          nc2 <- nc_open(output_nc_path())
          on.exit(nc_close(nc2), add = TRUE)
          
          # get variable array
          varok <- tryCatch(ncvar_get(nc2, vv), error = function(e) NULL)
          if (is.null(varok)) {
            plot.new(); title(main = paste0("Unable to read variable: ", vv)); return()
          }
          arr <- varok
          darr <- dim(arr)
          if (is.null(darr)) {
            # scalar - draw message
            plot.new(); title(main = paste0("Variable ", vv, " is scalar."))
            return()
          }
          
          # find time dim index
          t_idx <- detect_time_dim_index(nc2, vv)
          # Compute time average: collapse the time dimension if found, else if array has 3 dims guess the 3rd
          if (!is.null(t_idx) && length(darr) >= t_idx) {
            spatial_idx <- setdiff(seq_along(darr), t_idx)
            # apply mean over the t_idx
            avg <- apply(arr, spatial_idx, mean, na.rm = TRUE)
            # avg may be vector (if spatial_idx length 1) or matrix
            # ensure matrix orientation: we want matrix with dims [nrow, ncol]
            if (is.null(dim(avg))) {
              # vector -> reshape
              avg <- matrix(avg, nrow = darr[spatial_idx[1]], ncol = ifelse(length(spatial_idx) > 1, darr[spatial_idx[2]], 1))
            }
          } else if (length(darr) == 3) {
            # assume the third is time
            avg <- apply(arr, c(1,2), mean, na.rm = TRUE)
          } else if (length(darr) == 2) {
            avg <- arr
          } else {
            # Unexpected shape - show message
            plot.new(); title(main = paste0("Cannot compute time-average for variable: ", vv)); return()
          }
          
          # find lon/lat
          coords2 <- find_coord_vars(nc2)
          lon <- NULL; lat <- NULL
          if (!is.null(coords2$lon)) {
            lon_val <- tryCatch(ncvar_get(nc2, coords2$lon), error = function(e) NULL)
            if (!is.null(lon_val)) lon <- as.numeric(lon_val)
          }
          if (!is.null(coords2$lat)) {
            lat_val <- tryCatch(ncvar_get(nc2, coords2$lat), error = function(e) NULL)
            if (!is.null(lat_val)) lat <- as.numeric(lat_val)
          }
          
          # If lon/lat are still NULL, create simple axes
          if (is.null(lon)) lon <- seq_len(dim(avg)[2])
          if (is.null(lat)) lat <- seq_len(dim(avg)[1])
          
          # Try to align orientation: many netCDFs have dimension order [y,x], so avg is [nrow, ncol]
          # use image() — users can tweak orientation if needed
          # flip lat if needed so map looks upright
          if (length(lat) == nrow(avg) && length(lon) == ncol(avg)) {
            # OK dims match
            image(x = lon, y = lat, z = t(avg[nrow(avg):1, , drop = FALSE]), 
                  xlab = "Longitude", ylab = "Latitude", main = paste0("Time-mean of ", vv))
          } else {
            # fallback: just plot matrix with axis indices
            image(t(avg[nrow(avg):1, , drop = FALSE]), axes = FALSE, main = paste0("Time-mean of ", vv))
            axis(1, at = seq(0,1,length.out = 5), labels = round(seq(min(lon), max(lon), length.out = 5),2))
            axis(2, at = seq(0,1,length.out = 5), labels = round(seq(min(lat), max(lat), length.out = 5),2))
          }
        }, bg = "white")
      })
    }
    
    # close the file we opened earlier in this observe
    nc_close(nc)
  })

  output$download_nc <- downloadHandler(
    filename = function() {
      req(output_nc_path())                 # require an available path
      basename(output_nc_path())            # present original filename to user
    },
    content = function(file) {
      req(output_nc_path())
      # copy the file to the destination path 'file' (required by downloadHandler)
      file.copy(output_nc_path(), file, overwrite = TRUE)
    },
    contentType = "application/x-netcdf"    # appropriate mime type
  )
  # Inspect NetCDF and populate only spatio-temporal variables
  observeEvent(input$ncfile, {
    req(input$ncfile)
    append_log("Opening NetCDF to examine variables...")
    nc_path <- input$ncfile$datapath
    nc <- tryCatch(nc_open(nc_path), error = function(e) { append_log(paste("Error opening nc:", e$message)); NULL })
    if (is.null(nc)) return()
    
    spt_vars <- character(0)
    try({
      spt_vars <- find_spt_vars(nc)
    }, silent = TRUE)
    nc_close(nc)
    
    if (length(spt_vars) == 0) {
      append_log("No spatio-temporal variables detected. Showing all variables as fallback.")
      nc2 <- nc_open(nc_path)
      spt_vars <- names(nc2$var)
      nc_close(nc2)
    } else {
      append_log(paste0("Spatio-temporal variables detected: ", paste(spt_vars, collapse = ", ")))
    }
    
    output$var_ui <- renderUI({
      selectizeInput("varnames", "Select variable(s) to mask (spatio-temporal only)",
                     choices = spt_vars, multiple = TRUE,
                     options = list(placeholder = "Select one or more variables (optional)"))
    })
  })
  
  # Extract shapefile zip and keep path to .shp
  shp_path_r <- reactiveVal(NULL)
  observeEvent(input$maskzip, {
    req(input$maskzip)
    append_log("Extracting shapefile ZIP...")
    tmpdir <- tempfile("shp_unzip_")
    dir.create(tmpdir)
    unzip(input$maskzip$datapath, exdir = tmpdir)
    shp_files <- list.files(tmpdir, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)
    if (length(shp_files) == 0) {
      append_log("No .shp file found in the uploaded ZIP.")
      shp_path_r(NULL)
    } else {
      shp_path_r(shp_files[1])
      append_log(paste0("Found shapefile: ", basename(shp_files[1])))
    }
  })
  
  # Run nc_spatial_mask when user clicks generate; always pass ncFileOut = NULL
  observeEvent(input$generate, {
    req(input$ncfile)
    req(shp_path_r())   # extracted .shp path
    append_log("Starting nc_spatial_mask...")
    
    # vars to pass
    vars_arg <- if (!is.null(input$varnames) && length(input$varnames) > 0) input$varnames else NULL
    
    # create work dir
    work_dir <- tempfile("nc_run_")
    dir.create(work_dir)
    append_log(paste0("Created work directory: ", work_dir))
    
    # copy uploaded netcdf into work_dir using original filename
    orig_nc_name <- input$ncfile$name
    nc_copy_path <- file.path(work_dir, orig_nc_name)
    if (!file.copy(from = input$ncfile$datapath, to = nc_copy_path, overwrite = TRUE)) {
      append_log("❌ Failed to copy uploaded NetCDF into work_dir.")
      return()
    }
    append_log(paste0("Copied uploaded NetCDF to: ", nc_copy_path))
    
    # copy shapefile components into work_dir
    shp_fullpath <- shp_path_r()
    shp_basename <- sub("\\.shp$", "", basename(shp_fullpath), ignore.case = TRUE)
    shp_dir <- dirname(shp_fullpath)
    shp_components <- list.files(shp_dir, pattern = paste0("^", shp_basename, "\\.(shp|shx|dbf|prj|cpg)$"),
                                 full.names = TRUE, ignore.case = TRUE)
    if (length(shp_components) == 0) {
      # fallback: copy every file from the extracted shapefile folder (safe)
      shp_components <- list.files(shp_dir, full.names = TRUE)
    }
    for (f in shp_components) {
      file.copy(f, file.path(work_dir, basename(f)), overwrite = TRUE)
    }
    shp_work_shp <- file.path(work_dir, paste0(shp_basename, ".shp"))
    append_log(paste0("Copied shapefile components to work dir. shp used: ", shp_work_shp))
    
    # prepare naming info
    nc_basename <- sub("\\.\\w+$", "", basename(orig_nc_name))
    inpdir <- dirname(nc_copy_path)
    
    append_log(paste("Calling nc_spatial_mask with ncFile:", basename(nc_copy_path),
                     " maskFile:", basename(shp_work_shp),
                     " ncFileOut: <NULL>",
                     " var:", ifelse(is.null(vars_arg), "<NULL>", paste(vars_arg, collapse = ","))))
    
    # list files before
    before_files <- list.files(inpdir, pattern = "\\.nc$", full.names = TRUE)
    if (length(before_files) == 0) {
      append_log(paste0("Files in work dir before run: <none> (", inpdir, ")"))
    } else {
      append_log(paste0("Files in work dir before run: ", paste(basename(before_files), collapse = "; ")))
    }
    
    # run function and capture stdout messages/warnings
    captured_output <- character(0)
    returned_value <- NULL
    run_ok <- FALSE
    
    run_res <- tryCatch({
      out_capture <- capture.output({
        returned_value <<- withCallingHandlers({
          nc_spatial_mask(ncFile = nc_copy_path,
                          maskFile = shp_work_shp,
                          ncFileOut = NULL,    # ALWAYS NA as requested
                          var = vars_arg)
        }, message = function(m) {
          msg <- paste0("[message] ", m$message)
          append_log(msg); captured_output <<- c(captured_output, msg); invokeRestart("muffleMessage")
        }, warning = function(w) {
          msg <- paste0("[warning] ", w$message)
          append_log(msg); captured_output <<- c(captured_output, msg); invokeRestart("muffleWarning")
        })
      }, type = "output")
      if (length(out_capture) > 0) {
        append_log(paste0("[stdout] ", paste(out_capture, collapse = "\n")))
        captured_output <<- c(captured_output, out_capture)
      }
      run_ok <<- TRUE
      list(success = TRUE)
    }, error = function(e) {
      append_log(paste0("❌ nc_spatial_mask raised an error: ", e$message))
      list(success = FALSE, error = e$message)
    })
    
    # list files after
    after_files <- list.files(inpdir, pattern = "\\.nc$", full.names = TRUE)
    if (length(after_files) == 0) {
      append_log(paste0("Files in work dir after run: <none> (", inpdir, ")"))
    } else {
      append_log(paste0("Files in work dir after run: ", paste(basename(after_files), collapse = "; ")))
    }
    
    # Decide final output:
    final_out <- NULL
    
    # 1) If function returned a path (character), prefer that
    if (!is.null(returned_value) && is.character(returned_value) && length(returned_value) == 1) {
      # if relative path returned, make absolute relative to work_dir (some functions may return basename only)
      candidate_path <- if (file.exists(returned_value)) returned_value else file.path(inpdir, returned_value)
      if (file.exists(candidate_path)) {
        final_out <- normalizePath(candidate_path, winslash = "/")
        append_log(paste0("Function returned output path: ", final_out))
      } else {
        append_log(paste0("Function returned a path but file not found: ", returned_value))
      }
    }
    
    # 2) If no returned path, look for predictable filenames: *_aggregated.nc then *_masked.nc
    if (is.null(final_out) && length(after_files) > 0) {
      candidate_patterns <- c(
        paste0("^", nc_basename, "_aggregated\\.nc$"),
        paste0("^", nc_basename, "_masked\\.nc$")
      )
      bn_after <- basename(after_files)
      for (pat in candidate_patterns) {
        hits <- after_files[grepl(pat, bn_after, ignore.case = TRUE)]
        if (length(hits) > 0) {
          final_out <- hits[which.max(file.info(hits)$mtime)]
          append_log(paste0("Found output by pattern '", pat, "': ", basename(final_out)))
          break
        }
      }
      # 3) fallback: most recent .nc in work dir
      if (is.null(final_out)) {
        final_out <- after_files[which.max(file.info(after_files)$mtime)]
        append_log(paste0("No pattern match. Falling back to most recent .nc: ", basename(final_out)))
      }
    }
    
    # set result for download
    if (!is.null(final_out) && file.exists(final_out)) {
      append_log(paste0("✅ Masked/aggregated netCDF located: ", final_out))
      # copy to tempdir if you want the file to persist after work_dir cleanup (optional)
      safe_dst <- file.path(tempdir(), basename(final_out))
      file.copy(final_out, safe_dst, overwrite = TRUE)
      output_nc_path(safe_dst)
      append_log(paste0("Output copied to: ", safe_dst))
    } else {
      append_log("⚠️ Could not locate output netCDF file after running nc_spatial_mask.")
      if (length(captured_output) > 0) {
        append_log("Captured output from nc_spatial_mask():")
        for (line in captured_output) append_log(paste0("  ", line))
      }
      output_nc_path(NULL)
    }
  })
  
  # download handler for the created netcdf
  output$download_nc <- downloadHandler(
    filename = function() {
      req(output_nc_path())
      basename(output_nc_path())
    },
    content = function(file) {
      req(output_nc_path())
      file.copy(output_nc_path(), file, overwrite = TRUE)
    },
    contentType = "application/x-netcdf"
  )
}

shinyApp(ui = ui, server = server)

