# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

createBCT <- function(TL, S) {
    .Call('_destim_createBCT', PACKAGE = 'destim', TL, S)
}

updateCTaddtran <- function(CT, newt, stillt, sizeCT, tCT) {
    .Call('_destim_updateCTaddtran', PACKAGE = 'destim', CT, newt, stillt, sizeCT, tCT)
}

createrectangleCT <- function(TL, x, y) {
    .Call('_destim_createrectangleCT', PACKAGE = 'destim', TL, x, y)
}

is_sortedTL <- function(TL) {
    .Call('_destim_is_sortedTL', PACKAGE = 'destim', TL)
}

orderTL <- function(TL) {
    .Call('_destim_orderTL', PACKAGE = 'destim', TL)
}

findTorder <- function(TL, T) {
    .Call('_destim_findTorder', PACKAGE = 'destim', TL, T)
}

frbind <- function(mat1, mat2) {
    .Call('_destim_frbind', PACKAGE = 'destim', mat1, mat2)
}

createEQ <- function(tran, ncol) {
    .Call('_destim_createEQ', PACKAGE = 'destim', tran, ncol)
}

faddconstraint <- function(tran, TL, CT, stillt, S) {
    .Call('_destim_faddconstraint', PACKAGE = 'destim', tran, TL, CT, stillt, S)
}

faddconstraintm <- function(CT, CTS) {
    .Call('_destim_faddconstraintm', PACKAGE = 'destim', CT, CTS)
}

createTM <- function(TL, transitions, states) {
    .Call('_destim_createTM', PACKAGE = 'destim', TL, transitions, states)
}

createrectangleTL <- function(x, y) {
    .Call('_destim_createrectangleTL', PACKAGE = 'destim', x, y)
}

createsteady <- function(TM) {
    .Call('_destim_createsteady', PACKAGE = 'destim', TM)
}

createsteadypattern <- function(TM) {
    .Call('_destim_createsteadypattern', PACKAGE = 'destim', TM)
}

createsteadyfrompattern <- function(TM, ptrn) {
    .Call('_destim_createsteadyfrompattern', PACKAGE = 'destim', TM, ptrn)
}

createtransmatrix <- function(CT) {
    .Call('_destim_createtransmatrix', PACKAGE = 'destim', CT)
}

createlconmatrix <- function(CT) {
    .Call('_destim_createlconmatrix', PACKAGE = 'destim', CT)
}

floglik <- function(TM, values, rparams, transmatrix, init, EM, obs) {
    .Call('_destim_floglik', PACKAGE = 'destim', TM, values, rparams, transmatrix, init, EM, obs)
}

fforward <- function(TM, init, EM, obs) {
    .Call('_destim_fforward', PACKAGE = 'destim', TM, init, EM, obs)
}

fbackward <- function(TM, EM, obs, sfactors) {
    .Call('_destim_fbackward', PACKAGE = 'destim', TM, EM, obs, sfactors)
}

fscpstates <- function(TM, alpha, beta, sfactors, EM, obs) {
    .Call('_destim_fscpstates', PACKAGE = 'destim', TM, alpha, beta, sfactors, EM, obs)
}

cppfunique <- function(mat, tol) {
    .Call('_destim_cppfunique', PACKAGE = 'destim', mat, tol)
}

cppfuniqueind <- function(mat, tol) {
    .Call('_destim_cppfuniqueind', PACKAGE = 'destim', mat, tol)
}

inittransitions <- function(CT) {
    .Call('_destim_inittransitions', PACKAGE = 'destim', CT)
}

