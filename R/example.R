#' Example of using \pkg{destim} package
#'
#' @description This is just an example on how to compute the location and joint
#' location probabilities #' using simulated data. All the files used in this
#' example are supposed to be produced using the simulation software.
#' The "simulation.xml" file is an exception and it is an input file for the
#' simulation software. The files used in this example are provided with
#' the \pkg{destim} package.
#'
#' @details This is just an example on how to compute the location and joint
#' location probabilities #' using simulated data. All the files used in this
#' example are supposed to be produced using the simulation software.
#' The "simulation.xml" file is an exception and it is an input file for the
#' simulation software. The files used in this example are provided with
#' the \pkg{destim} package.
#' @references \url{https://github.com/MobilePhoneESSnetBigData}
#'
#' @examples
#'
#'library(data.table)
#'library(tidyr)
#'library(stringr)
#'library(Matrix)
#'library(xml2)

#'# set file names
#'path_root = 'extdata'
#'fileGridName       <- system.file(path_root, 'grid.csv', package = 'destim')
#'fileEventsInfoName <- system.file(path_root, 'AntennaInfo_MNO_MNO1.csv', package = 'destim')
#'signalFileName <- system.file(path_root, 'SignalMeasure_MNO1.csv', package = 'destim')
#'simulationFileName <- system.file(path_root, 'simulation.xml', package = 'destim')
#'# read simulation params
#'simulation.xml  <- as_list(read_xml(simulationFileName))
#'simulation.xml <- simulation.xml$simulation
#'start_time <- as.numeric(simulation.xml$start_time)
#'end_time <- as.numeric(simulation.xml$end_time)
#'time_increment <- as.numeric(simulation.xml$time_increment)
#'times <-
#'  seq(from = start_time,
#'      to = (end_time - time_increment),
#'      by = time_increment)
#'sigMin <- as.numeric(simulation.xml$conn_threshold)
#'
#'# read grid params
#'gridParam <-
#'  fread(
#'    fileGridName,
#'    sep = ',',
#'    header = TRUE,
#'    stringsAsFactors = FALSE
#'  )
#'ncol_grid  <- gridParam[['No Tiles Y']]
#'nrow_grid  <- gridParam[['No Tiles X']]
#'tile_sizeX <- gridParam[['X Tile Dim']]
#'tile_sizeY <- gridParam[['Y Tile Dim']]
#'ntiles     <- ncol_grid * nrow_grid
#'
#'# tile-rasterCell equivalence
#'tileEquiv.dt <- data.table(tileEquivalence(ncol_grid, nrow_grid))
#'
#'# read Received Signal Strength file and compute emission probabilities
#'RSS <-
#'  fread(
#'    signalFileName,
#'    sep = ",",
#'    header = TRUE,
#'    stringsAsFactors = FALSE
#'  )
#'setnames(RSS, c('antennaID', 0:(ntiles - 1)))
#'RSS <- melt(
#'  RSS,
#'  id.vars = 'antennaID',
#'  variable.name = 'tile',
#'  variable.factor = FALSE,
#'  value.name = 'RSS'
#')
#'RSS[, RSS := ifelse(RSS < sigMin, NA, RSS)]
#'
#'# compute event location (emission probabilities)
#'RSS <-
#'  RSS[, eventLoc := 10 ** RSS / sum(10 ** RSS, na.rm = TRUE), by = 'tile']
#'RSS <- RSS[is.na(eventLoc), eventLoc := 0]
#'RSS[, tile := as.numeric(tile)]
#'RSS <- RSS[tileEquiv.dt, on = 'tile'][, tile := NULL]
#'RSS <-
#'  dcast(RSS, rasterCell ~ antennaID, value.var = 'eventLoc')[, rasterCell := NULL]
#'emissionProbs <- Matrix(data = as.matrix(RSS))
#'dimnames(emissionProbs)[[1]] <-
#'  as.character(1:dim(emissionProbs)[1])
#'
#'# read and process network event data
#'allEvents.dt <-
#'  fread(
#'    fileEventsInfoName,
#'    sep = ',',
#'    stringsAsFactors = FALSE,
#'    colClasses = c(
#'      'integer',
#'      'character',
#'      'character',
#'      'character',
#'      'numeric',
#'      'numeric',
#'      'character'
#'    )
#'  )
#'allEvents.dt <- allEvents.dt[!duplicated(allEvents.dt)]
#'setnames(allEvents.dt ,
#'         c('time', 'antennaID', 'eventCode', 'device', 'x', 'y', 'tile'))
#'allEvents.dt[, obsVar := do.call(paste, c(.SD, sep = "-")),
#'             .SDcols = c('antennaID', 'eventCode')]
#'events.dt <- allEvents.dt[eventCode %in% c('0', '2', '3')]
#'events.dt_noDup <-
#'  copy(events.dt)[, list(eventCode = as.character(min(as.numeric(eventCode)))),
#'                  by = c("time", "device")]
#'events.dt <-
#'  merge(events.dt_noDup,
#'        events.dt,
#'        by = names(events.dt_noDup),
#'        all.x = TRUE)
#'events.dt <-
#'  events.dt[!duplicated(events.dt, by = c("time", "device", "eventCode"))][,
#'.(time, device, eventCode, antennaID, obsVar)][order(time)]
#'
#'# Set maximum velocity (from an external source)
#'vMax_ms <- 16
#'# Set time padding params
#'distMax <- vMax_ms * time_increment
#'pad_coef <-
#'  as.integer(ceiling(distMax / max(tile_sizeX, tile_sizeY)))
#'pad_coef <- pad_coef + 1
#'
#'# Initial state distribution (PRIOR)
#'
#'# Prepare prior_network distribution (uniform prior)
#'prior_network <- rep(1 / ntiles, ntiles)
#'
#'# Initialize HMM
#'model <- HMMrectangle(nrow_grid, ncol_grid)
#'emissions(model) <- emissionProbs
#'
#'model <- initparams(model)  # initialize transition prob
#'model <-
#'  minparams(model)   # parameter reduction according to restrictions
#'istates(model) <- prior_network
#'
#'# comute posterior location probabilities
#'deviceIDs <- sort(unique(events.dt$device))
#'
#'# for each device
#'for (i in seq(along = deviceIDs)) {
#'  devID <- deviceIDs[i]
#'  cat(paste0('    device ', devID, '...\n'))
#'  cat(' Selecting network events...')
#'  events_device.dt <- events.dt[device == devID, .(device, time, antennaID)][
#'                                                          order(device, time)]
#'
#'
#'  antennas_deviceID  <- unlist(events_device.dt[, c("antennaID")])
#'  if (!all(is.na(antennas_deviceID))) {
#'    # Fit and compute HMM model
#'    observedValues_pad <-
#'      rep(NA, pad_coef * length(antennas_deviceID))
#'    observedValues_pad[seq(1, length(observedValues_pad), by = pad_coef)] <-
#'      antennas_deviceID
#'    colEvents <- sapply(observedValues_pad,
#'                        function(x)
#'                          ifelse(!is.na(x), which(x == colnames(emissionProbs)), NA))
#'
#'    # Fit HMM - ML estimation of transition probabilities
#'    fitTry <-
#'      try(model_devID <- fit(model, colEvents, init = TRUE))
#'    if (inherits(fitTry, "try-error")) {
#'      stop("Fit model fails")
#'    }
#'    ssTry <- try(A <- sstates(model_devID, colEvents))
#'    if (inherits(ssTry, "try-error")) {
#'      stop("[compute_HMM] Smooth States fails")
#'    }
#'    B <- scpstates(model_devID, colEvents)
#'    # Transform output of the HMM model to sparse matrix file format
#'    transform_output <- transform_postLoc(
#'      postLocP = A,
#'      postLocJointP = B,
#'      observedValues = antennas_deviceID,
#'      times = times,
#'      t_increment = time_increment,
#'      ntiles = ntiles,
#'      pad_coef = pad_coef,
#'      tileEquiv.dt = tileEquiv.dt,
#'      devID = devID,
#'      sparse_postLocP = TRUE,
#'      sparse_postLocJointP = TRUE
#'    )
#'    rm(A, B)
#'    gc()
#'    fwrite(
#'      transform_output$postLocProb[, .(tile , time, postLocProb)],
#'      paste0('postLocProb_', devID, '.csv'),
#'      col.names = FALSE,
#'      row.names = FALSE,
#'      sep = ','
#'    )
#'    transform_output$postLocJointProb[, time_to := time_from + time_increment]
#'    fwrite(
#'      transform_output$postLocJointProb[, .(time_from, time_to, tile_from,
#'                                            tile_to, postLocProb)],
#'      paste0('postLocJointProb_', devID, '.csv'),
#'      col.names = FALSE,
#'      row.names = FALSE,
#'      sep = ','
#'    )
#'}
#'}
#'
example <- function() {}
