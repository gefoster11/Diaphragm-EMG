
# Shiny App: ECG & EMG Processing with RMS Peak Detection, Interactive Plots, Peak Toggling,
# Reset, Tooltips, Suggested Peak Height, Clustering Logic, Adjustable Cluster Window
library(shiny)
library(tidyverse)
library(signal)
library(RcppRoll)
library(plotly)
library(pracma)

# Package Dependency
packages = c("remotes", "shiny", "shinyBS", "tidyverse", "thematic", "shinythemes", "signal", "RcppRoll", "pracma")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# Install plotly from GitHub if needed (keeping your original logic)
if (!require("plotly", character.only = TRUE)) {
  remotes::install_github("ropensci/plotly")
  library("plotly", character.only = TRUE)
}

options(shiny.maxRequestSize = 50 * 1024^2, shiny.launch.browser = .rs.invokeShinyWindowExternal, scipen = 999)

ui <- fluidPage(
  titlePanel("Diaphragm EMG Processing Tool"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload Tab-Delimited File", accept = ".txt"),
      
      # ---- Time window controls (auto-filled to min/max after upload)
      numericInput("start_time", "Start time (s)", value = NA_real_, step = 0.1),
      numericInput("end_time",   "End time (s)",   value = NA_real_, step = 0.1),
      
      numericInput("threshold",    "R-wave Detection Threshold", value = NA,   step = 0.05),
      numericInput("clean_window", "EMG Cleaning Window (s)",    value = 0.2,  step = 0.05),
      numericInput("rms_window",   "RMS Window (ms)",            value = 100,  step = 10),
      
      sliderInput("low_cut",  "Low Cutoff Frequency (Hz)",  min = 0.05, max = 1, value = 0.1, step = 0.01),
      sliderInput("high_cut", "High Cutoff Frequency (Hz)", min = 0.2,  max = 2, value = 0.5, step = 0.01),
      
      # ---- Channel selector for visualization
      selectizeInput(
        "viz_channels", "Channels to visualize",
        choices = NULL, selected = NULL, multiple = TRUE,
        options = list(plugins = list("remove_button"))
      ),
      
      # ---- NEW: Toggle to analyze selected channels only
      checkboxInput("analyze_selected", "Analyze selected channels only", value = TRUE),
      
      tags$hr(),
      tags$div(HTML("
2025; Created by Glen Foster
"))
    ),
    mainPanel(
      tabsetPanel(id = "mainTabs", type = "tabs",
                  
                  # --- Tab 1: Overview (MWI + EMG)
                  tabPanel("ECG & EMG",
                           fluidRow(column(12, plotlyOutput("mwiPlot", height = "20vh"))),
                           fluidRow(column(12, plotlyOutput("emgPlot", height = "60vh"))),
                           downloadButton("downloadCombined", "Download Cleaned EMG + RMS")
                  ),
                  
                  # --- Tab 2: RMS & Clusters
                  tabPanel("RMS Analysis",
                           fluidRow(
                             column(3, numericInput("peak_height",   "Minimum RMS Peak Height",       value = 0.1,  step = 0.05)),
                             column(3, numericInput("peak_distance", "Minimum RMS Peak Distance (ms)", value = 1000, step = 100)),
                             column(3, numericInput("cluster_window","Cluster Window (ms)",            value = 1000, step = 50))
                           ),
                           fluidRow(
                             column(3, textOutput("suggestedPeakHeight")),
                             column(3, actionButton("resetPeaks", "Reset Toggled Peaks"))
                           ),
                           fluidRow(tags$hr()),
                           fluidRow(column(12, plotlyOutput("rmsPlot", height = "60vh"))),
                           fluidRow(tags$hr()),
                           fluidRow(column(12, tableOutput("rmsTable"))),
                           downloadButton("downloadRMSPeaks", "Download RMS Peak Summary")
                  )
      )
    )
  )
)

server <- function(input, output, session) {
  
  session$onSessionEnded(function() { stopApp() })
  

# ---- Helper: robust reader that auto-detects delimiter & normalizes ----
read_emg_file <- function(path) {
  # Peek first line to detect delimiter
  line1 <- tryCatch(readr::read_lines(path, n_max = 1), error = function(e) "")
  delim <- if (grepl("\t", line1)) "\t" else if (grepl(",", line1)) "," else if (grepl(";", line1)) ";" else "ws"

  # Read with the detected delimiter
  d <- tryCatch(
    {
      if (delim == "ws") {
        readr::read_table(path, col_names = FALSE, progress = FALSE, trim_ws = TRUE)
      } else {
        readr::read_delim(path, delim = delim, col_names = FALSE, show_col_types = FALSE, trim_ws = TRUE)
      }
    },
    error = function(e) tibble::tibble()  # empty tibble on hard failure
  )

  # Fallback: if only 1 column was read, try generic whitespace once more
  if (ncol(d) < 2L) {
    d <- tryCatch(readr::read_table(path, col_names = FALSE, progress = FALSE, trim_ws = TRUE),
                  error = function(e) tibble::tibble())
  }

  # If still <2 columns, return empty (caller will validate and message the user)
  if (ncol(d) < 2L) return(d)

  # Name columns: time, ECG, EMG1..EMGk
  nm <- seq_len(ncol(d))
  colnames(d)[nm[1]] <- "time"
  colnames(d)[nm[2]] <- "ECG"
  if (ncol(d) >= 3L) {
    for (i in nm[3:length(nm)]) colnames(d)[i] <- paste0("EMG", i - 2)
  }

  # Coerce all columns to numeric (non-numeric -> NA)
  for (nm in names(d)) {
    d[[nm]] <- suppressWarnings(as.numeric(d[[nm]]))
  }

  # Drop rows where time is NA; order & de-duplicate by time
  if (!("time" %in% names(d))) return(tibble::tibble())  # safety guard
  d <- d[!is.na(d$time), , drop = FALSE]
  if (nrow(d) > 0) d <- d[order(d$time), , drop = FALSE]
  if (nrow(d) > 0) d <- d[!duplicated(d$time), , drop = FALSE]

  d}

  
  # --- Toggle store for peak markers
  toggled_keys <- reactiveVal(character(0))
  observeEvent(input$resetPeaks, { toggled_keys(character(0)) })
  
  # --- Raw data (read once per upload)
  
  raw_data <- reactive({
    req(input$file)
    
    d <- read_emg_file(input$file$datapath)
    
    # Validate we have at least time + ECG
    validate(
      shiny::need(ncol(d) >= 2, "Unable to parse file: expected at least 2 columns (time, ECG)."),
      shiny::need("time" %in% names(d), "Unable to parse file: missing 'time' column after import."),
      shiny::need("ECG"  %in% names(d), "Unable to parse file: missing 'ECG' column after import."),
      shiny::need(nrow(d) > 0, "No usable rows found (check delimiter and header).")
    )
    
    d
  })
  
  
  
  # --- Initialize start/end time and channel selector after file upload
  
  
  observeEvent(raw_data(), {
    d <- raw_data()
    req(nrow(d) > 0, "time" %in% names(d))
    
    tmin <- min(d$time, na.rm = TRUE)
    tmax <- max(d$time, na.rm = TRUE)
    if (!is.finite(tmin) || !is.finite(tmax) || tmin >= tmax) {
      showNotification("Cannot establish a valid time range from the file.", type = "error")
      return()
    }
    
    updateNumericInput(session, "start_time", value = tmin)
    updateNumericInput(session, "end_time",   value = tmax)
    
    ch <- names(d)[grepl("^EMG", names(d))]
    updateSelectizeInput(session, "viz_channels", choices = ch, selected = ch, server = TRUE)
  }, ignoreInit = FALSE)
  
  
  
  # --- Full processing pipeline (respects time window + analyzer toggle)
  
  processed_data <- reactive({
    d <- raw_data()
    req(d)
    
    # EMG columns present (require at least 1 overall)
    emg_cols_all <- names(d)[grepl("^EMG", names(d))]
    req(length(emg_cols_all) >= 1)
    
    
    
    # --- Time window (robust defaults) ---
    st <- suppressWarnings(as.numeric(input$start_time))
    en <- suppressWarnings(as.numeric(input$end_time))
    
    if (!is.finite(st) || !is.finite(en)) {
      st <- min(d$time, na.rm = TRUE)
      en <- max(d$time, na.rm = TRUE)
    }
    
    validate(shiny::need(is.finite(st) && is.finite(en) && st < en,
                         "Start time must be less than End time."))
    
    # Base subsetting to avoid NSE/masking
    idx <- which(d$time >= st & d$time <= en)
    validate(shiny::need(length(idx) > 1, "No data in the selected time window (or window too small)."))
    d <- d[idx, , drop = FALSE]
    
    # --- Sampling rate guards ---
    dt <- diff(d$time)
    validate(shiny::need(length(dt) > 0 && any(is.finite(dt)),
                         "Selected time window is too small to estimate sampling rate."))
    fs <- 1 / mean(dt, na.rm = TRUE)
    validate(shiny::need(is.finite(fs) && fs > 0, "Invalid sampling rate in selected window."))
    
  
    
    # Channels to ANALYZE (toggle)
    viz <- input$viz_channels
    channels_to_analyze <- if (isTRUE(input$analyze_selected)) {
      intersect(viz, emg_cols_all)   # <-- fixed
    } else {
      emg_cols_all
    }
    req(length(channels_to_analyze) >= 1)
    
    
    # ECG MWI → R-peaks → Cleaning → RMS → Filter (unchanged from your logic)
    bp <- signal::butter(2, c(5, 15)/(fs/2), type = "pass")
    ecg_filt <- signal::filtfilt(bp, d$ECG)
    ecg_diff <- c(0, diff(ecg_filt))
    ecg_sq   <- ecg_diff^2
    win_size <- round(0.15 * fs)
    ecg_mwi  <- as.numeric(stats::filter(ecg_sq, rep(1/win_size, win_size), sides = 1))
    ecg_mwi[is.na(ecg_mwi)] <- 0
    
    observeEvent(ecg_mwi, {
      if (is.na(input$threshold)) {
        updateNumericInput(session, "threshold", value = max(ecg_mwi) * 0.5)
      }
    })
    
    candidate_peaks <- which(ecg_mwi > input$threshold)
    refractory <- 0.25 * fs
    r_peaks <- c()
    if (length(candidate_peaks) > 0) {
      regions <- split(candidate_peaks, cumsum(c(1, diff(candidate_peaks) > 1)))
      peak_indices <- sapply(regions, function(r) r[which.max(ecg_mwi[r])])
      r_peaks <- peak_indices[c(1, which(diff(peak_indices) > refractory) + 1)]
    }
    
    half_window <- floor((input$clean_window * fs)/2)
    cleaned_emg <- d %>% dplyr::select(dplyr::all_of(channels_to_analyze))
    for (ch in names(cleaned_emg)) {
      emg <- cleaned_emg[[ch]]
      for (rp in r_peaks) {
        start_idx <- max(1, rp - half_window)
        end_idx   <- min(nrow(d), rp + half_window)
        mid <- floor((start_idx + end_idx)/2)
        if (start_idx > half_window) {
          target_len  <- mid - start_idx + 1
          source_vals <- emg[(start_idx - half_window):(start_idx - 1)]
          emg[start_idx:mid] <- source_vals[seq_len(min(target_len, length(source_vals)))]
        }
        if (end_idx + half_window <= nrow(d)) {
          target_len  <- end_idx - mid
          source_vals <- emg[(end_idx + 1):(end_idx + half_window)]
          emg[(mid + 1):end_idx] <- source_vals[seq_len(min(target_len, length(source_vals)))]
        }
      }
      cleaned_emg[[ch]] <- emg
    }
    
    rms_window_samples <- round((input$rms_window/1000) * fs)
    rms_results <- cleaned_emg %>%
      dplyr::mutate(dplyr::across(dplyr::everything(),
                                  ~ RcppRoll::roll_mean(.^2, rms_window_samples, fill = NA) |> sqrt()))
    bp_resp <- signal::butter(2, c(input$low_cut, input$high_cut)/(fs/2), type = "pass")
    rms_filtered <- rms_results %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ {
        vals <- na.omit(.)
        if (length(vals) > 3) {
          filtered_vals <- signal::filtfilt(bp_resp, vals)
          out <- rep(NA, length(.))
          out[!is.na(.)] <- filtered_vals
          out
        } else rep(NA, length(.))
      }))
    
    list(
      data         = d,
      ecg_mwi      = ecg_mwi,
      r_peaks      = r_peaks,
      cleaned_emg  = cleaned_emg,
      rms          = rms_results,
      rms_filtered = rms_filtered,
      fs           = fs
    )
  })
  
  # --- Peak detection on filtered RMS (channels from processed_data)
  rms_peaks_df <- reactive({
    dp <- processed_data()
    req(dp$rms, dp$rms_filtered, dp$data$time, dp$fs)
    req(input$peak_height, input$peak_distance, input$low_cut, input$high_cut)
    
    channels <- names(dp$rms)[grepl("^EMG", names(dp$rms))]
    req(length(channels) >= 1)
    
    min_distance_samples <- max(1L, as.integer(round((input$peak_distance / 1000) * dp$fs)))
    
    out_list <- lapply(channels, function(ch) {
      filt <- as.numeric(dp$rms_filtered[[ch]])
      raw  <- as.numeric(dp$rms[[ch]])
      time <- dp$data$time
      
      if (length(filt) < 3 || length(raw) != length(filt) || length(time) != length(filt)) {
        return(tibble(Channel = ch, sample_idx = integer(0), time = numeric(0), RMS = numeric(0)))
      }
      
      good <- which(is.finite(filt))
      if (length(good) < 3) {
        return(tibble(Channel = ch, sample_idx = integer(0), time = numeric(0), RMS = numeric(0)))
      }
      
      filt_good <- filt[good]
      pk <- tryCatch(
        pracma::findpeaks(filt_good, minpeakheight = input$peak_height, minpeakdistance = min_distance_samples),
        error = function(e) NULL
      )
      if (is.null(pk) || nrow(pk) == 0) {
        tibble(Channel = ch, sample_idx = integer(0), time = numeric(0), RMS = numeric(0))
      } else {
        idx_rel <- as.integer(pk[, 2])
        idx     <- good[idx_rel]
        tibble(Channel = ch, sample_idx = idx, time = time[idx], RMS = raw[idx])
      }
    })
    
    dplyr::bind_rows(out_list)
  })
  
  # --- Suggested peak height
  output$suggestedPeakHeight <- renderText({
    dp <- processed_data()
    all_filtered_vals <- unlist(dp$rms_filtered)
    guessed_height <- suppressWarnings(quantile(all_filtered_vals, 0.60, na.rm = TRUE))
    paste("Suggested peak height:", round(guessed_height, 3))
  })
  
  # --- ECG MWI plot
  output$mwiPlot <- renderPlotly({
    dp <- processed_data()
    diagnostic_data <- tibble(time = dp$data$time, ecg_mwi = dp$ecg_mwi)
    p <- ggplot(diagnostic_data, aes(x = time, y = ecg_mwi)) +
      geom_line(color = "orange") +
      geom_hline(yintercept = input$threshold, color = "red", linetype = "dashed") +
      geom_point(data = diagnostic_data[dp$r_peaks, ], aes(x = time, y = ecg_mwi), color = "blue", size = 2) +
      labs(title = "ECG MWI with User Threshold", x = "Time (s)", y = "Integrated Signal") +
      theme_minimal()
    ggplotly(p)
  })
  
  # --- EMG plot (Raw vs Cleaned) — visualization respects selected channels
  output$emgPlot <- renderPlotly({
    dp <- processed_data()
    
    ch_all <- names(dp$cleaned_emg)
    viz <- input$viz_channels
    if (is.null(viz) || length(viz) == 0) {
      return(
        plotly_empty(type = "scatter", mode = "lines") %>%
          layout(title = "Select at least one channel to visualize", xaxis = list(title = "Time (s)"),
                 yaxis = list(title = "EMG"))
      )
    }
    viz <- intersect(viz, ch_all)
    if (length(viz) == 0) {
      return(
        plotly_empty(type = "scatter", mode = "lines") %>%
          layout(title = "Selected channels not found in analyzed set", xaxis = list(title = "Time (s)"),
                 yaxis = list(title = "EMG"))
      )
    }
    
    raw_emg <- dp$data %>% dplyr::select(time, dplyr::all_of(viz)) %>% dplyr::mutate(type = "Raw")
    cleaned_emg <- dp$cleaned_emg %>% dplyr::select(dplyr::all_of(viz)) %>%
      dplyr::mutate(time = dp$data$time, type = "Cleaned")
    
    emg_long <- dplyr::bind_rows(raw_emg, cleaned_emg) %>%
      tidyr::pivot_longer(cols = viz, names_to = "Channel", values_to = "Value")
    
    y_min <- suppressWarnings(min(emg_long$Value, na.rm = TRUE))
    y_max <- suppressWarnings(max(emg_long$Value, na.rm = TRUE))
    
    plots <- lapply(viz, function(ch) {
      raw_df   <- emg_long %>% dplyr::filter(Channel == ch, type == "Raw")
      clean_df <- emg_long %>% dplyr::filter(Channel == ch, type == "Cleaned")
      plot_ly() %>%
        add_trace(data = raw_df,   x = ~time, y = ~Value, type = 'scatter', mode = 'lines',
                  line = list(color = 'blue'),   name = 'Raw') %>%
        add_trace(data = clean_df, x = ~time, y = ~Value, type = 'scatter', mode = 'lines',
                  line = list(color = 'orange'), name = 'Cleaned') %>%
        layout(title = ch, xaxis = list(title = "Time (s)"),
               yaxis = list(title = "EMG", range = c(y_min, y_max)))
    })
    
    fig <- subplot(plots, nrows = length(plots), shareX = TRUE, titleX = TRUE) %>%
      layout(title = "Raw vs Cleaned EMG", showlegend = FALSE)
    
    if (length(plots) > 1) {
      for (i in 2:length(plots)) {
        fig$x$layout[[paste0("yaxis", i)]]$matches <- "y"
      }
    }
    fig
  })
  
  # --- RMS plot — visualization respects selected channels; analysis respects toggle
  output$rmsPlot <- renderPlotly({
    dp <- processed_data()
    
    req(dp, dp$rms, dp$rms_filtered, dp$data$time)
    
    pk <- rms_peaks_df()
    keys <- toggled_keys()
    
    cat("[render] toggled_keys:", if (length(keys)) paste(keys, collapse = ", ") else "<empty>", "\n")
    
    ch_all <- names(dp$rms)
    emg_mask <- grepl("^EMG", ch_all)
    if (!any(emg_mask)) {
      return(
        plotly_empty(type = "scatter", mode = "lines") %>%
          layout(
            title = "No EMG channels found (need at least one)",
            xaxis = list(title = "Time (s)"),
            yaxis = list(title = "RMS")
          )
      )
    }
    
    viz <- input$viz_channels
    if (is.null(viz) || length(viz) == 0) {
      return(
        plotly_empty(type = "scatter", mode = "lines") %>%
          layout(
            title = "Select at least one channel to visualize",
            xaxis = list(title = "Time (s)"),
            yaxis = list(title = "RMS")
          )
      )
    }
    channels <- intersect(viz, ch_all[emg_mask])
    if (length(channels) == 0) {
      return(
        plotly_empty(type = "scatter", mode = "lines") %>%
          layout(
            title = "Selected channels not found in analyzed set",
            xaxis = list(title = "Time (s)"),
            yaxis = list(title = "RMS")
          )
      )
    }
    
    
    
    plots <- lapply(channels, function(ch) {
      pk_ch <- pk %>% dplyr::filter(Channel == ch)
      if (nrow(pk_ch) > 0) {
        pk_ch <- pk_ch %>% dplyr::mutate(key = paste(Channel, sample_idx, sep = "||"))
        cat("[render] ch:", ch, "matched:", sum(pk_ch$key %in% keys), "/", nrow(pk_ch), "\n")
        # ... split to red/grey and add traces as above ...
      }
      
      df_raw  <- tibble::tibble(time = dp$data$time, value = as.numeric(dp$rms[[ch]]))
      df_filt <- tibble::tibble(time = dp$data$time, value = as.numeric(dp$rms_filtered[[ch]]))
      
      p <- plot_ly(source = "rms") %>%
        add_lines(data = df_filt, x = ~time, y = ~value, name = paste(ch, "Filtered RMS"),
                  line = list(color = "grey"), legendgroup = ch, showlegend = FALSE, hoverinfo = "none") %>%
        add_lines(data = df_raw,  x = ~time, y = ~value, name = paste(ch, "Raw RMS"),
                  line = list(color = "purple"), legendgroup = ch, showlegend = FALSE, hoverinfo = "none")
      
      if (nrow(pk) > 0) {
        pk_ch <- pk %>% dplyr::filter(Channel == ch)
        if (nrow(pk_ch) > 0 && all(c("time", "RMS") %in% names(pk_ch))) {
          pk_ch <- pk_ch %>%
            dplyr::mutate(
              time = as.numeric(time),
              RMS  = as.numeric(RMS),
              key = paste(Channel, sample_idx, sep = "||"),
              color = ifelse(key %in% keys, "grey", "red")
            )
          good_pts <- which(is.finite(pk_ch$time) & is.finite(pk_ch$RMS))
          if (length(good_pts) > 0) {
            pk_good <- pk_ch[good_pts, , drop = FALSE]
            pk_red  <- pk_good[pk_good$color == "red",  , drop = FALSE]
            pk_grey <- pk_good[pk_good$color == "grey", , drop = FALSE]
            
            
            # Red markers (included)
            if (nrow(pk_red) > 0) {
              p <- p %>% add_markers(
                data = pk_red, x = ~time, y = ~RMS,
                name = paste(ch, "Peaks (included)"),
                marker = list(size = 12, color = "red"),
                key = ~key,             # <-- ensure character keys
                customdata = ~key,
                text = ~paste("Channel:", Channel, "\nTime:", round(time, 3), "\nRMS:", round(RMS, 3)),
                hoverinfo = "text", legendgroup = ch, showlegend = FALSE
              )
            }
            
            # Grey markers (excluded)
            if (nrow(pk_grey) > 0) {
              p <- p %>% add_markers(
                data = pk_grey, x = ~time, y = ~RMS,
                name = paste(ch, "Peaks (excluded)"),
                marker = list(size = 6, color = "grey"),
                key = ~key,             # <-- ensure character keys
                customdata = ~key,
                text = ~paste("Channel:", Channel, "\nTime:", round(time, 3), "\nRMS:", round(RMS, 3)),
                hoverinfo = "text", legendgroup = ch, showlegend = FALSE
              )
            }
            
          }
        }
      }
      
      p %>% layout(title = list(text = ch),
                   xaxis = list(title = "Time (s)"),
                   yaxis = list(title = "RMS"))
    })
    
    
    fig <- subplot(plots, nrows = length(plots), shareX = TRUE, titleY = TRUE) %>%
      layout(title = "RMS vs Filtered RMS with Detected Peaks", margin = list(t = 35))
    fig$x$source <- "rms"
    fig <- fig %>% plotly::event_register("plotly_click")
    fig
    
  })
  
  # --- Plotly click to toggle peak inclusion
  
  observeEvent(plotly::event_data("plotly_click", source = "rms"), {
    d <- plotly::event_data("plotly_click", source = "rms")
    if (is.null(d)) return()
    
    # Prefer key, fall back to customdata — they are identical now
    k <- if (!is.null(d$key) && length(d$key) && nzchar(d$key)) d$key
    else if (!is.null(d$customdata) && length(d$customdata) && nzchar(d$customdata)) d$customdata
    else return()
    
    cur <- toggled_keys()
    if (k %in% cur) toggled_keys(setdiff(cur, k)) else toggled_keys(c(cur, k))

  }, ignoreInit = TRUE)
  
  
  
  # --- Excluded (grey) and included (red) peaks
  
  excluded_peaks <- reactive({
    pk   <- rms_peaks_df()
    keys <- toggled_keys()
    if (nrow(pk) == 0 || length(keys) == 0) return(pk[0, ])
    
    pk_key <- dplyr::mutate(pk, key = paste(Channel, sample_idx, sep = "||"))
    dplyr::semi_join(pk_key, tibble::tibble(key = keys), by = "key") %>%
      dplyr::select(Channel, sample_idx, time, RMS)
  })
  
  included_peaks <- reactive({
    pk <- rms_peaks_df()
    if (nrow(pk) == 0) return(pk)
    
    keys <- toggled_keys()
    if (length(keys) == 0) return(pk)
    
    pk %>%
      dplyr::mutate(key = paste(Channel, sample_idx, sep = "||")) %>%
      dplyr::anti_join(tibble::tibble(key = keys), by = "key") %>%
      dplyr::select(-key)
  })
  
  # --- Clustered summary (across ANALYZED channels only)
  
  clustered_rms <- reactive({
    df <- included_peaks()
    if (nrow(df) == 0) {
      return(tibble::tibble(
        cluster_id   = integer(0),
        cluster_time = numeric(0),
        RMSmax       = numeric(0)
      ))
    }
    
    # Sort by time
    df <- dplyr::arrange(df, time)
    
    # Cluster window (ms -> s)
    cw_ms <- if (is.null(input$cluster_window)) 200 else input$cluster_window
    cw <- cw_ms / 1000
    
    # Cluster IDs across all channels by inter-event gap
    cluster_id <- cumsum(c(1L, diff(df$time) > cw))
    df <- dplyr::mutate(df, cluster_id = cluster_id)
    
    # Per cluster & channel: max RMS
    per_ch <- df %>%
      dplyr::group_by(cluster_id, Channel) %>%
      dplyr::summarize(RMS = max(RMS, na.rm = TRUE), .groups = "drop")
    
    # One row per cluster with earliest time
    base <- df %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarize(cluster_time = min(time), .groups = "drop")
    
    # Wide format: one column per (analyzed) EMG channel
    
    # Wide format: one column per (analyzed) EMG channel
    wide <- base %>%
      dplyr::left_join(per_ch, by = "cluster_id") %>%
      tidyr::pivot_wider(names_from = Channel, values_from = RMS)
    
    # Compute RMSmax across all EMG columns present (base-safe, no '.')
    emg_cols <- grep("^EMG\\d+$", names(wide), value = TRUE)
    
    if (length(emg_cols) == 0) {
      wide$RMSmax <- NA_real_
    } else {
      # Sort EMG columns by numeric suffix: EMG1, EMG2, ... EMG10, ...
      emg_order <- emg_cols[order(as.integer(sub("^EMG", "", emg_cols)))]
      # pmax across the (sorted) EMG columns
      wide$RMSmax <- do.call(pmax, c(as.data.frame(wide[emg_order]), na.rm = TRUE))
      # Reorder columns: cluster_id, cluster_time, EMG1..EMGk (sorted), RMSmax
      wide <- wide %>%
        dplyr::relocate(cluster_id, cluster_time) %>%
        dplyr::select(cluster_id, cluster_time, dplyr::all_of(emg_order), RMSmax)
    }
    
  })
  
  
  output$rmsTable <- renderTable({
    clustered_rms()
  }, digits = 3, striped = TRUE, bordered = TRUE, hover = TRUE)
  
  # --- Downloads (respect time window + analyzer toggle)
  output$downloadRMSPeaks <- downloadHandler(
    filename = function() { paste0("RMS_clusters_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv") },
    content  = function(file) { readr::write_csv(clustered_rms(), file) }
  )
  
  output$downloadCombined <- downloadHandler(
    filename = function() { paste0("cleaned_emg_plus_rms_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv") },
    content  = function(file) {
      dp <- processed_data()
      req(dp$cleaned_emg, dp$rms, dp$data$time)
      
      time    <- dp$data$time
      cleaned <- dp$cleaned_emg %>% dplyr::mutate(time = time) %>% dplyr::relocate(time)
      rms_raw <- dp$rms         %>% dplyr::mutate(time = time) %>% dplyr::relocate(time)
      rms_cols <- setdiff(names(rms_raw), "time")
      names(rms_raw)[names(rms_raw) != "time"] <- paste0(rms_cols, "_RMS")
      
      combined <- cleaned %>% dplyr::left_join(rms_raw, by = "time")
      readr::write_csv(combined, file)
    }
  )
}

shinyApp(ui, server)
