
# Shiny App: ECG & EMG Processing with RMS Peak Detection, Interactive Plots, Peak Toggling, Reset, Tooltips, Suggested Peak Height, Clustering Logic, Adjustable Cluster Window
library(shiny)
library(tidyverse)
library(signal)
library(RcppRoll)
library(plotly)
library(pracma)

# Package Dependency
packages = c("remotes",
             "shiny",
             "shinyBS",
             "tidyverse",
             "thematic",
             "shinythemes",
             "signal",
             "RcppRoll",
             "pracma"
)

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Package Dependency - Install plotly from github
if (!require("plotly", character.only = TRUE)) {
  remotes::install_github("ropensci/plotly")
  library("plotly", character.only = TRUE)
}

options(shiny.maxRequestSize = 50 * 1024^2,
shiny.launch.browser = .rs.invokeShinyWindowExternal,
scipen = 999)



ui <- fluidPage(
  titlePanel("Diaphragm EMG Processing Tool"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload Tab-Delimited File", accept = ".txt"),
      numericInput("threshold", "R-wave Detection Threshold", value = NA, step = 0.05),
      numericInput("clean_window", "EMG Cleaning Window (s)", value = 0.2, step = 0.05),
      numericInput("rms_window", "RMS Window (ms)", value = 100, step = 10),
      # numericInput("peak_height", "Minimum RMS Peak Height", value = 0.1, step = 0.05),
      # textOutput("suggestedPeakHeight"),
      # numericInput("peak_distance", "Minimum RMS Peak Distance (ms)", value = 1000, step = 100),
      # numericInput("cluster_window", "Cluster Window (ms)", value = 1000, step = 50),
      sliderInput("low_cut", "Low Cutoff Frequency (Hz)", min = 0.05, max = 1, value = 0.1, step = 0.01),
      sliderInput("high_cut", "High Cutoff Frequency (Hz)", min = 0.2, max = 2, value = 0.5, step = 0.01),
      # downloadButton("downloadCombined", "Download Cleaned EMG + RMS"),
      # downloadButton("downloadRMSPeaks", "Download RMS Peak Summary"),
      # ---- Horizontal line ----
      tags$hr(),
      
      # ---- Reference ----
      tags$div(
        HTML("<p>2025; Created by Glen Foster</p><br>")),
    ),
    mainPanel(
      tabsetPanel(id = "mainTabs", type = "tabs",
                  
                  # --- Tab 1: Overview (MWI fixed height + EMG fills rest) ---
                  
                  tabPanel("ECG & EMG",
                           fluidRow(
                             column(
                               width = 12,
                               # ECG/MWI takes 30% of viewport height
                               plotlyOutput("mwiPlot", height = "20vh")
                             )
                           ),
                           fluidRow(
                             column(
                               width = 12,
                               # EMG takes the remaining 70%
                               plotlyOutput("emgPlot", height = "60vh")
                             )
                           ),
                           downloadButton("downloadCombined", "Download Cleaned EMG + RMS"),
                  ),
                  
                  # --- Tab 2: RMS & Clusters ---
                  
                  
                  tabPanel("RMS Analysis",
                           fluidRow(
                             column(width = 3,
                                    numericInput("peak_height", "Minimum RMS Peak Height", value = 0.1, step = 0.05)
                             ),
                             column(width = 3,
                                    numericInput("peak_distance", "Minimum RMS Peak Distance (ms)", value = 1000, step = 100)
                             ),
                             column(width = 3,
                                    numericInput("cluster_window", "Cluster Window (ms)", value = 1000, step = 50)
                             ),
                           ),
                           fluidRow(
                             column(width = 3,
                                    textOutput("suggestedPeakHeight")
                             ),
                             column(width = 3,
                                    actionButton("resetPeaks", "Reset Toggled Peaks")
                             ),
                           ),
                           # ---- Horizontal line ----
                           fluidRow(
                             tags$hr(),  
                           ),
                           fluidRow(
                             column(width = 12, plotlyOutput("rmsPlot", height = "60vh"))
                           ),
                           # ---- Horizontal line ----
                           fluidRow(
                             tags$hr(),  
                           ),
                           fluidRow(
                             column(width = 12, tableOutput("rmsTable"))
                           ),
                           downloadButton("downloadRMSPeaks", "Download RMS Peak Summary")
                  )
      )
    )
  )
)



server <- function(input, output, session) {
  
  session$onSessionEnded(function() {
    stopApp()
  })
  
  # ReactiveVal that stores toggled keys (grey points)
  toggled_keys <- reactiveVal(character(0))
  
  # If you want to keep the old table too, you can maintain both:
  # toggled_peaks <- reactiveVal(tibble(Channel = character(), sample_idx = integer(), time = numeric()))
   
   # Optional: reset button wiring (you already have input$resetPeaks)
   observeEvent(input$resetPeaks, {
     toggled_keys(character(0))
   })
   
  processed_data <- reactive({
    req(input$file)
    data <- read_delim(input$file$datapath, delim = "	", col_names = FALSE)
    colnames(data)[1] <- "time"
    colnames(data)[2] <- "ECG"
    for (i in 3:ncol(data)) colnames(data)[i] <- paste0("EMG", i - 2)
    
    fs <- 1 / mean(diff(data$time))
    
    bp <- butter(2, c(5, 15)/(fs/2), type = "pass")
    ecg_filt <- filtfilt(bp, data$ECG)
    ecg_diff <- c(0, diff(ecg_filt))
    ecg_sq <- ecg_diff^2
    win_size <- round(0.15 * fs)
    ecg_mwi <- as.numeric(stats::filter(ecg_sq, rep(1/win_size, win_size), sides = 1))
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
    cleaned_emg <- data %>% select(starts_with("EMG"))
    for (ch in names(cleaned_emg)) {
      emg <- cleaned_emg[[ch]]
      for (rp in r_peaks) {
        start_idx <- max(1, rp - half_window)
        end_idx <- min(nrow(data), rp + half_window)
        mid <- floor((start_idx + end_idx)/2)
        if (start_idx > half_window) {
          target_len <- mid - start_idx + 1
          source_vals <- emg[(start_idx - half_window):(start_idx - 1)]
          emg[start_idx:mid] <- source_vals[seq_len(min(target_len, length(source_vals)))]
        }
        if (end_idx + half_window <= nrow(data)) {
          target_len <- end_idx - mid
          source_vals <- emg[(end_idx + 1):(end_idx + half_window)]
          emg[(mid + 1):end_idx] <- source_vals[seq_len(min(target_len, length(source_vals)))]
        }
      }
      cleaned_emg[[ch]] <- emg
    }
    
    rms_window_samples <- round((input$rms_window/1000) * fs)
    rms_results <- cleaned_emg %>% mutate(across(everything(), ~ RcppRoll::roll_mean(.^2, rms_window_samples, fill = NA) |> sqrt()))
    
    bp_resp <- signal::butter(2, c(input$low_cut, input$high_cut)/(fs/2), type = "pass")
    rms_filtered <- rms_results %>% mutate(across(everything(), ~ {
      vals <- na.omit(.)
      if (length(vals) > 3) {
        filtered_vals <- signal::filtfilt(bp_resp, vals)
        out <- rep(NA, length(.))
        out[!is.na(.)] <- filtered_vals
        out
      } else rep(NA, length(.))
    }))
    
    list(data = data, ecg_mwi = ecg_mwi, r_peaks = r_peaks, cleaned_emg = cleaned_emg, rms = rms_results, rms_filtered = rms_filtered, fs = fs)
  })
  
  
  
  
  rms_peaks_df <- reactive({
    dp <- processed_data()
    req(dp$rms, dp$rms_filtered, dp$data$time, dp$fs)
    req(input$peak_height, input$peak_distance, input$low_cut, input$high_cut)
    
    # Limit to first 5 EMG channels (or fewer if present)
    ch_all   <- names(dp$rms)
    channels <- ch_all[grepl("^EMG", ch_all)]
    channels <- channels[seq_len(min(5, length(channels)))]
    
    # Peak detection parameters (in samples)
    min_distance_samples <- max(1L, as.integer(round((input$peak_distance / 1000) * dp$fs)))
    
    # We do NOT re-filter here; we use dp$rms_filtered as provided.
    out_list <- lapply(channels, function(ch) {
      filt <- as.numeric(dp$rms_filtered[[ch]])
      raw  <- as.numeric(dp$rms[[ch]])
      time <- dp$data$time
      
      # Guard clauses
      if (length(filt) < 3 || length(raw) != length(filt) || length(time) != length(filt)) {
        return(tibble(Channel = ch, sample_idx = integer(0), time = numeric(0), RMS = numeric(0)))
      }
      
      # Work only on non-NA samples of the filtered RMS
      good <- which(is.finite(filt))
      if (length(good) < 3) {
        return(tibble(Channel = ch, sample_idx = integer(0), time = numeric(0), RMS = numeric(0)))
      }
      
      filt_good <- filt[good]
      
      # Find peaks on filtered RMS (indices are relative to 'filt_good')
      pk <- tryCatch(
        pracma::findpeaks(filt_good,
                          minpeakheight   = input$peak_height,
                          minpeakdistance = min_distance_samples),
        error = function(e) NULL
      )
      
      if (is.null(pk) || nrow(pk) == 0) {
        tibble(Channel = ch, sample_idx = integer(0), time = numeric(0), RMS = numeric(0))
      } else {
        # pk[,2] are indices into filt_good; map them back to global sample indices
        idx_rel <- as.integer(pk[, 2])
        idx     <- good[idx_rel]
        
        tibble(
          Channel    = ch,
          sample_idx = idx,          # global sample index
          time       = time[idx],
          RMS        = raw[idx]      # raw RMS at the filtered-peak locations
        )
      }
    })
    
    dplyr::bind_rows(out_list)
  })
  
  
    
  
  
  output$suggestedPeakHeight <- renderText({
    dp <- processed_data()
    all_filtered_vals <- unlist(dp$rms_filtered)
    guessed_height <- quantile(all_filtered_vals, 0.60, na.rm = TRUE)
    paste("Suggested peak height:", round(guessed_height, 3))
  })
  
  output$mwiPlot <- renderPlotly({
    dp <- processed_data()
    diagnostic_data <- tibble(time = dp$data$time, ecg_mwi = dp$ecg_mwi)
    p <- ggplot(diagnostic_data, aes(x = time, y = ecg_mwi)) +
      geom_line(color = "orange") +
      geom_hline(yintercept = input$threshold, color = "red", linetype = "dashed") +
      geom_point(data = diagnostic_data[dp$r_peaks, ], aes(x = time, y = ecg_mwi), color = "blue", size = 2) +
      labs(title = "ECG MWI with User Threshold", x = "Time (s)", y = "Integrated Signal") + theme_minimal()
    ggplotly(p)
  })
  
  
  
  output$emgPlot <- renderPlotly({
    dp <- processed_data()
    
    # Prepare data
    raw_emg <- dp$data %>%
      select(time, starts_with("EMG")) %>%
      mutate(type = "Raw")
    
    cleaned_emg <- dp$cleaned_emg %>%
      mutate(time = dp$data$time, type = "Cleaned")
    
    emg_long <- bind_rows(raw_emg, cleaned_emg) %>%
      pivot_longer(cols = starts_with("EMG"), names_to = "Channel", values_to = "Value")
    
    # Compute global y-axis range
    y_min <- min(emg_long$Value, na.rm = TRUE)
    y_max <- max(emg_long$Value, na.rm = TRUE)
    
    # Create individual plots for each channel
    channels <- unique(emg_long$Channel)
    plots <- lapply(channels, function(ch) {
      raw_df <- emg_long %>% filter(Channel == ch, type == "Raw")
      clean_df <- emg_long %>% filter(Channel == ch, type == "Cleaned")
      
      plot_ly() %>%
        add_trace(data = raw_df, x = ~time, y = ~Value,
                  type = 'scatter', mode = 'lines',
                  line = list(color = 'blue'),
                  name = 'Raw') %>%
        add_trace(data = clean_df, x = ~time, y = ~Value,
                  type = 'scatter', mode = 'lines',
                  line = list(color = 'orange'),
                  name = 'Cleaned') %>%
        layout(title = ch,
               xaxis = list(title = "Time (s)"),
               yaxis = list(title = "EMG", range = c(y_min, y_max)))
    })
    
    # Combine plots and enforce axis matching
    fig <- subplot(plots, nrows = length(channels), shareX = TRUE, titleX = TRUE) %>%
      layout(title = "Raw vs Cleaned EMG", showlegend = FALSE)
    
    # Force all y-axes to match the first one
    for (i in 2:length(channels)) {
      fig$x$layout[[paste0("yaxis", i)]]$matches <- "y"
    }
    
    fig
  })
  
  
  
  output$rmsPlot <- renderPlotly({
    # --- Upstream data and reactive deps ---
    dp <- processed_data()
    req(dp, dp$rms, dp$rms_filtered, dp$data$time)
    
    # Peaks: find on filtered RMS, index raw RMS
    pk <- rms_peaks_df()  # or rms_peaks_df_debounced() if you added debounce
    
    # Toggle keys (grey points) — plot depends on this so it re-renders on click
    keys <- if (exists("toggled_keys")) {
      toggled_keys()
    } else if (exists("toggled_peaks") && nrow(toggled_peaks()) > 0) {
      with(toggled_peaks(), paste(Channel, sample_idx, sep = "|"))
    } else {
      character(0)
    }
    
    # --- Define channels safely (limit to 5 EMG columns) ---
    ch_all   <- names(dp$rms)
    emg_mask <- grepl("^EMG", ch_all)
    if (!any(emg_mask)) {
      return(
        plotly_empty(type = "scatter", mode = "lines") %>%
          layout(
            title = "No EMG channels found (expected names like EMG1, EMG2, ...)",
            xaxis = list(title = "Time (s)"),
            yaxis = list(title = "RMS")
          )
      )
    }
    channels <- ch_all[emg_mask]
    channels <- channels[seq_len(min(5, length(channels)))]
    
    # --- Build per-channel subplots (each will now auto-scale its own Y) ---
    plots <- lapply(channels, function(ch) {
      df_raw  <- tibble::tibble(time = dp$data$time, value = as.numeric(dp$rms[[ch]]))
      df_filt <- tibble::tibble(time = dp$data$time, value = as.numeric(dp$rms_filtered[[ch]]))
      
      # Base lines (always added)
      p <- plot_ly(source = "rms") %>%
        add_lines(
          data = df_filt, x = ~time, y = ~value,
          name = paste(ch, "Filtered RMS"),
          line = list(color = "grey"),
          legendgroup = ch, showlegend = FALSE
        ) %>%
        add_lines(
          data = df_raw, x = ~time, y = ~value,
          name = paste(ch, "Raw RMS"),
          line = list(color = "purple"),
          legendgroup = ch, showlegend = FALSE
        )
      
      # Add peaks only if we have them for this channel
      if (nrow(pk) > 0) {
        pk_ch <- pk %>% dplyr::filter(Channel == ch)
        
        if (nrow(pk_ch) > 0 && all(c("time", "RMS") %in% names(pk_ch))) {
          pk_ch <- pk_ch %>%
            dplyr::mutate(
              time = as.numeric(time),
              RMS  = as.numeric(RMS),
              key  = paste(Channel, sample_idx, sep = "|"),
              color = ifelse(key %in% keys, "grey", "red")
            )
          
          # keep only finite points
          good_pts <- which(is.finite(pk_ch$time) & is.finite(pk_ch$RMS))
          if (length(good_pts) > 0) {
            pk_good <- pk_ch[good_pts, , drop = FALSE]
            pk_red  <- pk_good[pk_good$color == "red",  , drop = FALSE]
            pk_grey <- pk_good[pk_good$color == "grey", , drop = FALSE]
            
            # Red markers (included)
            if (nrow(pk_red) > 0) {
              p <- p %>%
                add_markers(
                  data = pk_red,
                  x = ~time, y = ~RMS,
                  name = paste(ch, "Peaks (included)"),
                  marker = list(size = 12, color = "red"),
                  key = ~key,
                  customdata = ~sample_idx,
                  text = ~paste("Channel:", Channel,
                                "<br>Time:", round(time, 3),
                                "<br>RMS:", round(RMS, 3)),
                  hoverinfo = "text",
                  legendgroup = ch, showlegend = FALSE
                )
            }
            
            # Grey markers (excluded)
            if (nrow(pk_grey) > 0) {
              p <- p %>%
                add_markers(
                  data = pk_grey,
                  x = ~time, y = ~RMS,
                  name = paste(ch, "Peaks (excluded)"),
                  marker = list(size = 6, color = "grey"),
                  key = ~key,
                  customdata = ~sample_idx,
                  text = ~paste("Channel:", Channel,
                                "<br>Time:", round(time, 3),
                                "<br>RMS:", round(RMS, 3)),
                  hoverinfo = "text",
                  legendgroup = ch, showlegend = FALSE
                )
            }
          }
        }
      }
      
      # y-axis (per-plot auto range — no shared/global range)
      p %>% layout(
        title = list(text = ch),
        xaxis = list(title = "Time (s)"),
        yaxis = list(title = "RMS")  # no fixed 'range' applied
      )
    })
    
    subplot(plots, nrows = length(plots), shareX = TRUE, titleY = TRUE) %>%
      layout(title = "RMS vs Filtered RMS with Detected Peaks",
             margin = list(t = 35))
  })
  
  
  
  
  
  observeEvent(plotly::event_data("plotly_click", source = "rms"), {
    d <- plotly::event_data("plotly_click", source = "rms")
    if (is.null(d)) return()
    
    k <- d$key  # "Channel|sample_idx"
    if (is.null(k) || !nzchar(k)) return()
    
    cur <- toggled_keys()
    if (k %in% cur) {
      toggled_keys(setdiff(cur, k))  # Untoggle → back to red
    } else {
      toggled_keys(c(cur, k))        # Toggle → grey
    }
  })
  
  
  
  excluded_peaks <- reactive({
    pk <- rms_peaks_df()
    keys <- toggled_keys()
    if (nrow(pk) == 0 || length(keys) == 0) return(pk[0, ])
    pk_key <- transform(pk, key = paste(Channel, sample_idx, sep = "|"))
    dplyr::semi_join(pk_key, tibble::tibble(key = keys), by = "key") %>%
      dplyr::select(Channel, sample_idx, time, RMS)
  })
  
  
  # Included (red) peaks after toggling
  included_peaks <- reactive({
    pk <- rms_peaks_df()
    if (nrow(pk) == 0) return(pk)
    
    # Preferred: toggled_keys()
    if (exists("toggled_keys")) {
      keys <- toggled_keys()
    } else if (exists("toggled_peaks") && nrow(toggled_peaks()) > 0) {
      # Fallback: derive keys from toggled_peaks table
      keys <- with(toggled_peaks(), paste(Channel, sample_idx, sep = "|"))
    } else {
      keys <- character(0)
    }
    
    if (length(keys) == 0) return(pk)
    
    pk %>%
      dplyr::mutate(key = paste(Channel, sample_idx, sep = "|")) %>%
      dplyr::anti_join(tibble::tibble(key = keys), by = "key") %>%
      dplyr::select(-key)
  })
  
  
  clustered_rms <- reactive({
    df <- included_peaks()
    if (nrow(df) == 0) {
      # Create an empty table with expected columns
      return(tibble::tibble(
        cluster_id   = integer(0),
        cluster_time = numeric(0),
        EMG1 = numeric(0), EMG2 = numeric(0), EMG3 = numeric(0),
        EMG4 = numeric(0), EMG5 = numeric(0),
        RMSmax = numeric(0)
      ))
    }
    
    # Sort by time
    df <- dplyr::arrange(df, time)
    
    # Cluster window in seconds
    cw <- (input$cluster_window %||% 200) / 1000  # default 200 ms if NULL
    
    # Compute cluster IDs across all channels by time gap
    cluster_id <- cumsum(c(1L, diff(df$time) > cw))
    df <- dplyr::mutate(df, cluster_id = cluster_id)
    
    # For each cluster & channel, keep the **max RMS** within that window
    per_ch <- df %>%
      dplyr::group_by(cluster_id, Channel) %>%
      dplyr::summarize(RMS = max(RMS, na.rm = TRUE), .groups = "drop")
    
    # One row per cluster with cluster_time = min(time) in window
    base <- df %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarize(cluster_time = min(time), .groups = "drop")
    
    # Wide format: columns EMG1..EMG5 (NA if channel absent in that cluster)
    wide <- base %>%
      dplyr::left_join(per_ch, by = "cluster_id") %>%
      tidyr::pivot_wider(names_from = Channel, values_from = RMS)
    
    # Ensure all EMG1..EMG5 columns exist (fill missing ones with NA)
    target_ch <- paste0("EMG", 1:5)
    for (ch in target_ch) {
      if (!ch %in% names(wide)) wide[[ch]] <- NA_real_
    }
    wide <- wide %>%
      dplyr::select(cluster_id, cluster_time, dplyr::all_of(target_ch))
    
    # RMSmax across EMG1..EMG5 (row-wise maximum, ignoring NA)
    rms_cols <- target_ch[target_ch %in% names(wide)]
    wide <- wide %>%
      dplyr::mutate(RMSmax = do.call(pmax, c(dplyr::select(., dplyr::all_of(rms_cols)), na.rm = TRUE)))
    
    wide %>% dplyr::arrange(cluster_time)
  })
  
  
  # Show the clustered summary
  output$rmsTable <- renderTable({
    clustered_rms()
  }, digits = 3, striped = TRUE, bordered = TRUE, hover = TRUE)
  
  # Optional: download as CSV
  output$downloadRMSPeaks <- downloadHandler(
    filename = function() {
      paste0("RMS_clusters_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      readr::write_csv(clustered_rms(), file)
    }
  )
  
  
  # ---- Download: Cleaned EMG + Raw RMS (CSV only) ----
  output$downloadCombined <- downloadHandler(
    filename = function() {
      paste0("cleaned_emg_plus_rms_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      dp <- processed_data()
      req(dp$cleaned_emg, dp$rms, dp$data$time)
      
      # Base time vector
      time <- dp$data$time
      
      # Cleaned EMG + time
      cleaned <- dp$cleaned_emg %>%
        dplyr::mutate(time = time) %>%
        dplyr::relocate(time)
      
      # Raw RMS + time, rename to *_RMS
      rms_raw <- dp$rms %>%
        dplyr::mutate(time = time) %>%
        dplyr::relocate(time)
      rms_cols <- setdiff(names(rms_raw), "time")
      names(rms_raw)[names(rms_raw) != "time"] <- paste0(rms_cols, "_RMS")
      
      # Join by time (same sampling grid)
      combined <- cleaned %>%
        dplyr::left_join(rms_raw, by = "time")
      
      # Write CSV
      readr::write_csv(combined, file)
    }
  )
  
  
}

shinyApp(ui, server)