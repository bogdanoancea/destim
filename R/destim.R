#' destim: A package for mobile devices position estimation.
#' @title A package for mobile devices position estimation using HMM.
#'
#' @description  This package contains functions to compute the posterior
#' location probability for each device over a grid of tiles covering the
#' geographical area under consideration. It uses Hidden Markov Models.
#' The theory behind the method is described in detail in
#' \href{https://webgate.ec.europa.eu/fpfis/mwikis/essnetbigdata/images/f/fb/WPI_Deliverable_I3_A_proposed_production_framework_with_mobile_network_data_2020_05_31_draft.pdf}{WPI
#' Deliverable 3} and in the paper \emph{An end-to-end statistical process
#' with mobile network data for Official Statistics}. For an example on how to
#' use this package please read \link{example1} and \link{example2}.
#' @docType package
#' @name destim
#' @useDynLib destim
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
NULL
