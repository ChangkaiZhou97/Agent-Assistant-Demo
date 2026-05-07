# =============================================================================
# LIVEtools_extensions.R
# -----------------------------------------------------------------------------
# Add-on functions for the LIVEtools R package
# (https://github.com/johnmurraylab/LIVE_tools)
#
# Provides three new functions that operate on any CD-like data.frame
# (the standard object returned by readEmbryoTable()[[1]] or totalRePosition()):
#
#   1. extractSubset()        - pull a subset of cells (by names / lineages /
#                               time range) out of an embryo data.frame.
#   2. rotateSubset()          - rigid-body rotate (and optionally translate)
#                               a subset in 3D, either by explicit Euler
#                               angles, by AP/DV indicators, or by an
#                               arbitrary 3x3 rotation matrix.
#   3. alignSubsetToRelatives() - align a subset (translation +/- rotation) to
#                               the SAME cells' parents at an earlier time
#                               point or children at a later time point.
#                               Uses the standard C. elegans naming
#                               convention (Sulston) to find relatives.
#
# Usage:
#   library(LIVEtools)
#   source("LIVEtools_extensions.R")
#
# Author: generated as add-on to LIVEtools v1.0
# =============================================================================


# -----------------------------------------------------------------------------
# Internal helpers (re-implemented because LIVEtools does NOT export them)
# These mirror the logic in LIVEtools/R/tree_plots.R
# -----------------------------------------------------------------------------

.getParent <- function(x) {
  if (is.na(x) || x == "P" || x == "P0") return(NA_character_)
  last <- substr(x, nchar(x), nchar(x))
  if (last %in% c("a", "p", "d", "v", "l", "r")) {
    return(substr(x, 1, nchar(x) - 1))
  }
  if (x %in% c("AB", "P1"))   return("P0")
  if (x %in% c("EMS", "P2"))  return("P1")
  if (x %in% c("E", "MS"))    return("EMS")
  if (x %in% c("C", "P3"))    return("P2")
  if (x %in% c("D", "P4"))    return("P3")
  if (x %in% c("Z2", "Z3"))   return("P4")
  return("P")
}

# Children of a given cell (the inverse of .getParent). Each cell has up to 2.
.getChildren <- function(x) {
  if (is.na(x)) return(character(0))
  # Special founder cells
  special <- list(
    P0  = c("AB",  "P1"),
    P1  = c("EMS", "P2"),
    P2  = c("C",   "P3"),
    P3  = c("D",   "P4"),
    P4  = c("Z2",  "Z3"),
    EMS = c("MS",  "E")
  )
  if (x %in% names(special)) return(special[[x]])
  # AB / MS / E / C / D follow the lettered naming
  # Default: append a/p (anterior/posterior). Note: in C. elegans, axes used
  # for naming are context-dependent (a/p, d/v, l/r). We use a/p as the
  # canonical default because the StarryNite / AceTree convention is to
  # preserve whatever last letter the parent was named with by appending
  # one of the six axis letters. We list a/p first; user code that needs
  # d/v or l/r daughters can supply explicit cell names instead.
  return(c(paste0(x, "a"), paste0(x, "p")))
}


# -----------------------------------------------------------------------------
# 1. extractSubset()
# -----------------------------------------------------------------------------

#' extractSubset
#'
#' Pull a subset of cells out of an embryo CD-like data.frame.
#' Thin, well-documented wrapper around \code{LIVEtools::grepCells()} plus
#' an optional time-range filter and a "drop empty" safety check.
#'
#' @param CDFrame  A CD-like data.frame (must contain at least the columns
#'                 \code{cell}, \code{time}, \code{x}, \code{y}, \code{z}).
#' @param cells    Character vector of cell names. \code{"x"} acts as a
#'                 single-letter wildcard (e.g. \code{"MSxa"}). May be
#'                 \code{NULL} if \code{lineages} is given.
#' @param lineages Character vector of lineage / mother cell names. The
#'                 mother cell and ALL its descendants are returned. May be
#'                 \code{NULL} if \code{cells} is given.
#' @param times    Either \code{"ALL"} (default) for no time filter, a
#'                 numeric vector of specific time points, or a length-2
#'                 numeric vector \code{c(t_min, t_max)} when
#'                 \code{time_is_range = TRUE}.
#' @param time_is_range Logical. If \code{TRUE}, \code{times} is interpreted
#'                 as a closed range \code{[t_min, t_max]}.
#' @param verbose  Print a summary of what was extracted.
#'
#' @return A data.frame with the same columns as \code{CDFrame}, containing
#'         only the matching rows. Row names are reset.
#'
#' @examples
#' \dontrun{
#'   sub <- extractSubset(embryo, lineages = "MS", times = c(80, 140),
#'                        time_is_range = TRUE)
#' }
extractSubset <- function(CDFrame,
                          cells          = NULL,
                          lineages       = NULL,
                          times          = "ALL",
                          time_is_range  = FALSE,
                          verbose        = TRUE) {
  if (!is.data.frame(CDFrame)) {
    stop("`CDFrame` must be a data.frame (CD-like).")
  }
  required_cols <- c("cell", "time", "x", "y", "z")
  missing_cols  <- setdiff(required_cols, names(CDFrame))
  if (length(missing_cols) > 0) {
    stop("`CDFrame` is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  if (is.null(cells) && is.null(lineages)) {
    stop("Either `cells` or `lineages` (or both) must be supplied.")
  }

  # ---- Time filter ---------------------------------------------------------
  if (identical(times, "ALL")) {
    df <- CDFrame
  } else if (time_is_range) {
    if (length(times) != 2) {
      stop("When time_is_range = TRUE, `times` must have length 2: c(t_min, t_max).")
    }
    df <- CDFrame[CDFrame$time >= times[1] & CDFrame$time <= times[2], ]
  } else {
    df <- CDFrame[CDFrame$time %in% times, ]
  }

  # ---- Cell / lineage filter (delegated to LIVEtools) ----------------------
  out <- LIVEtools::grepCells(
    CDData   = df,
    cells    = cells,
    lineages = lineages,
    times    = "ALL"   # we already filtered by time above
  )
  rownames(out) <- NULL

  if (verbose) {
    message(sprintf(
      "extractSubset: %d rows kept (%d unique cells, time range [%s, %s]).",
      nrow(out),
      length(unique(out$cell)),
      ifelse(nrow(out) > 0, min(out$time), NA),
      ifelse(nrow(out) > 0, max(out$time), NA)
    ))
  }
  out
}


# -----------------------------------------------------------------------------
# 2. rotateSubset()
# -----------------------------------------------------------------------------

# Build a 3x3 rotation matrix from extrinsic Euler angles in radians
# (rotate around X, then Y, then Z) suitable for post-multiplication of
# row-vector coordinates: rotated_coords = coords %*% R.
# (This is the same convention LIVEtools uses internally.)
.eulerToMatrix <- function(rx, ry, rz) {
  Rx <- matrix(c(1, 0,        0,
                 0, cos(rx), -sin(rx),
                 0, sin(rx),  cos(rx)),
               nrow = 3, byrow = TRUE)
  Ry <- matrix(c( cos(ry), 0, sin(ry),
                  0,       1, 0,
                 -sin(ry), 0, cos(ry)),
               nrow = 3, byrow = TRUE)
  Rz <- matrix(c(cos(rz), -sin(rz), 0,
                 sin(rz),  cos(rz), 0,
                 0,        0,       1),
               nrow = 3, byrow = TRUE)
  M <- Rz %*% Ry %*% Rx          # standard column-vector convention
  t(M)                           # transpose -> row-vector convention
}

#' rotateSubset
#'
#' Rigid-body rotate (and optionally translate) a subset of an embryo in 3D.
#' The rotation is applied ONLY to the rows in \code{subsetCD}; nothing
#' outside the subset is touched. Three modes are supported:
#'
#' \enumerate{
#'   \item \strong{angles}:       supply Euler angles (radians or degrees).
#'   \item \strong{matrix}:       supply a 3x3 rotation matrix directly.
#'   \item \strong{indicators}:   supply AP / DV / LR indicator lineages
#'                                (same idea as \code{LIVEtools::RePosition})
#'                                and the function computes a rotation that
#'                                aligns the subset's principal axes to the
#'                                global x / y / z axes.
#' }
#'
#' @param subsetCD  A CD-like data.frame (typically from \code{extractSubset}).
#' @param mode      One of \code{"angles"}, \code{"matrix"}, \code{"indicators"}.
#' @param angles    Numeric vector of length 3, c(rx, ry, rz). Used when
#'                  \code{mode = "angles"}.
#' @param degrees   Logical. If \code{TRUE} (default) \code{angles} are in
#'                  degrees; if \code{FALSE}, in radians.
#' @param rotMatrix 3x3 numeric matrix. Used when \code{mode = "matrix"}.
#'                  The rotation is applied as
#'                  \code{coords \%*\% rotMatrix} (LIVEtools convention,
#'                  i.e. row vectors are post-multiplied).
#' @param indicatorP,indicatorD,indicatorV,indicatorL,indicatorR
#'                  Lineage names used to define AP and DV/LR axes.
#'                  Same semantics as \code{LIVEtools::RePosition}.
#' @param ref_time  Time point used to compute the rotation matrix when
#'                  \code{mode = "indicators"}. Defaults to the median time
#'                  in \code{subsetCD}.
#' @param recenter  Logical. If \code{TRUE} (default), translate the subset
#'                  centroid to the origin BEFORE rotating, then leave it
#'                  there. Set to \code{FALSE} to rotate around the global
#'                  origin in place.
#' @param translate Optional length-3 numeric vector. After rotation, add
#'                  this vector to (x, y, z). Useful for placing the subset
#'                  somewhere specific.
#'
#' @return A data.frame with the same columns and row order as
#'         \code{subsetCD}, with rotated coordinates. The function also
#'         attaches the rotation matrix that was applied as
#'         \code{attr(result, "rotation")} and the centering offset as
#'         \code{attr(result, "center")} for full reproducibility.
#'
#' @examples
#' \dontrun{
#'   sub <- extractSubset(embryo, lineages = "MS", times = 100)
#'   sub_rot <- rotateSubset(sub, mode = "angles", angles = c(0, 0, 90))
#' }
rotateSubset <- function(subsetCD,
                         mode       = c("angles", "matrix", "indicators"),
                         angles     = c(0, 0, 0),
                         degrees    = TRUE,
                         rotMatrix  = NULL,
                         indicatorP = "C",
                         indicatorD = "Cxa",
                         indicatorV = "MSxxp",
                         indicatorL = NULL,
                         indicatorR = NULL,
                         ref_time   = NULL,
                         recenter   = TRUE,
                         translate  = NULL) {
  mode <- match.arg(mode)
  if (!is.data.frame(subsetCD) || nrow(subsetCD) == 0) {
    stop("`subsetCD` must be a non-empty data.frame.")
  }

  # ---- 1. Build the rotation matrix ---------------------------------------
  R <- switch(mode,
    angles = {
      if (length(angles) != 3) stop("`angles` must have length 3.")
      if (degrees) angles <- angles * pi / 180
      .eulerToMatrix(angles[1], angles[2], angles[3])
    },
    matrix = {
      if (is.null(rotMatrix) || !all(dim(rotMatrix) == c(3, 3))) {
        stop("`rotMatrix` must be a 3x3 matrix when mode = 'matrix'.")
      }
      rotMatrix
    },
    indicators = {
      # Reproduce the LIVEtools::rotationVec logic in-line (it's not
      # exported). This builds a 3x3 rotation matrix M such that
      #   coords_rotated = coords_centered %*% M
      # has AP aligned to +x, DV to +z, LR to +y.
      if (is.null(ref_time)) ref_time <- median(subsetCD$time)
      ref_slice <- subsetCD[abs(subsetCD$time - ref_time) < 1e-6, ]
      if (nrow(ref_slice) < 4) {
        stop("Need at least 4 cells at ref_time = ", ref_time,
             " in the subset to fit AP/DV indicators.")
      }
      # Center the reference slice before fitting principal axes
      ref_centered <- ref_slice
      ref_centered$x <- ref_centered$x - mean(ref_centered$x)
      ref_centered$y <- ref_centered$y - mean(ref_centered$y)
      ref_centered$z <- ref_centered$z - mean(ref_centered$z)

      .orthProj  <- function(v, u) v - sum(v * u) * u
      .crossProd <- function(v, u) c(v[2]*u[3] - v[3]*u[2],
                                     v[3]*u[1] - v[1]*u[3],
                                     v[1]*u[2] - v[2]*u[1])

      # AP axis = first PC, sign-flipped so positive is toward indicatorP
      pc <- prcomp(ref_centered[, c("x", "y", "z")])
      AP <- pc$rotation[, 1]
      Pcells <- LIVEtools::grepCells(ref_centered, lineages = indicatorP)
      if (nrow(Pcells) > 0) {
        Pvec <- c(mean(Pcells$x), mean(Pcells$y), mean(Pcells$z))
        if (sum(Pvec * AP) < 0) AP <- -AP
      }
      # DV axis: weighted contributions from D/V/L/R indicators
      DV <- c(0, 0, 0)
      if (!is.null(indicatorD)) {
        Dc <- LIVEtools::grepCells(ref_centered, lineages = indicatorD)
        if (nrow(Dc) > 0) {
          DV <- DV + .orthProj(c(mean(Dc$x), mean(Dc$y), mean(Dc$z)), AP)
        }
      }
      if (!is.null(indicatorV)) {
        Vc <- LIVEtools::grepCells(ref_centered, lineages = indicatorV)
        if (nrow(Vc) > 0) {
          DV <- DV - .orthProj(c(mean(Vc$x), mean(Vc$y), mean(Vc$z)), AP)
        }
      }
      if (!is.null(indicatorR)) {
        Rc <- LIVEtools::grepCells(ref_centered, lineages = indicatorR)
        if (nrow(Rc) > 0) {
          DV <- DV - .crossProd(c(mean(Rc$x), mean(Rc$y), mean(Rc$z)), AP)
        }
      }
      if (!is.null(indicatorL)) {
        Lc <- LIVEtools::grepCells(ref_centered, lineages = indicatorL)
        if (nrow(Lc) > 0) {
          DV <- DV + .crossProd(c(mean(Lc$x), mean(Lc$y), mean(Lc$z)), AP)
        }
      }
      if (sqrt(sum(DV^2)) < 1e-10) {
        stop("Indicator-based DV axis is degenerate. ",
             "Check that the indicator lineages exist in the subset.")
      }
      DV <- DV / sqrt(sum(DV^2))
      RL <- .crossProd(-AP, DV)
      cbind(AP, RL, DV)
    }
  )

  # ---- 2. Optionally translate to origin ----------------------------------
  coords <- as.matrix(subsetCD[, c("x", "y", "z")])
  center <- if (recenter) colMeans(coords) else c(0, 0, 0)
  coords <- sweep(coords, 2, center)

  # ---- 3. Apply rotation --------------------------------------------------
  rotated <- coords %*% R

  # ---- 4. Optional translation back ---------------------------------------
  if (!is.null(translate)) {
    if (length(translate) != 3) stop("`translate` must have length 3.")
    rotated <- sweep(rotated, 2, translate, FUN = "+")
  } else if (!recenter) {
    # leave coords where they are - already computed without centering
  }

  out <- subsetCD
  out$x <- rotated[, 1]
  out$y <- rotated[, 2]
  out$z <- rotated[, 3]

  attr(out, "rotation")  <- R
  attr(out, "center")    <- center
  attr(out, "translate") <- if (is.null(translate)) c(0, 0, 0) else translate
  out
}


# -----------------------------------------------------------------------------
# 3. alignSubsetToRelatives()
# -----------------------------------------------------------------------------

# Procrustes-style rigid alignment of point set A onto point set B
# (no scaling). Both A and B are n x 3 matrices, rows already paired.
# Returns list(R = 3x3, t = length-3 translation, fit = aligned A).
.procrustesRT <- function(A, B) {
  centroidA <- colMeans(A)
  centroidB <- colMeans(B)
  Ac <- sweep(A, 2, centroidA)
  Bc <- sweep(B, 2, centroidB)
  H  <- t(Ac) %*% Bc
  svd_H <- svd(H)
  d <- sign(det(svd_H$v %*% t(svd_H$u)))
  R <- svd_H$v %*% diag(c(1, 1, d)) %*% t(svd_H$u)
  t <- as.numeric(centroidB - R %*% centroidA)
  list(R = R, t = t,
       fit = sweep(A %*% t(R), 2, t, FUN = "+"))
}

#' alignSubsetToRelatives
#'
#' Rigid-body align a subset of cells to the spatial positions of their
#' parents (at an earlier time point) or their children (at a later time
#' point). Useful for tracking how a group of nuclei moves between time
#' frames, or for stitching a subset back into a previous/future embryo
#' state for visualization.
#'
#' Two alignment modes:
#'
#' \describe{
#'   \item{\code{"translate"}}{Shift only - move the subset centroid to
#'         match the relatives' centroid. No rotation.}
#'   \item{\code{"rigid"}}{Procrustes-style rigid alignment - find the
#'         best rotation + translation (no scaling) that maps each subset
#'         cell onto its corresponding relative. Requires at least 3
#'         matching pairs.}
#' }
#'
#' @param subsetCD     A CD-like data.frame (output of \code{extractSubset}).
#'                     Must contain a single time point OR you must specify
#'                     \code{subset_time}.
#' @param fullCD       The full embryo data.frame from which to look up
#'                     parents/children (typically the same data.frame the
#'                     subset was extracted from).
#' @param relative     One of \code{"parent"} (use parents at an earlier
#'                     time) or \code{"child"} (use children at a later
#'                     time).
#' @param subset_time  Time point of the subset. Required if subset spans
#'                     multiple times; defaults to \code{unique(subsetCD$time)}
#'                     when there is only one.
#' @param target_time  Optional. If supplied (a single number), all relatives
#'                     are looked up at exactly this time. If \code{NULL}
#'                     (default), each relative is looked up at its OWN
#'                     last-observed time before \code{subset_time} (for
#'                     \code{relative = "parent"}) or first-observed time
#'                     after \code{subset_time} (for \code{relative = "child"}).
#'                     This per-cell default is much more robust because in
#'                     C. elegans different lineages divide at very different
#'                     times, so there is rarely a single time point where
#'                     every relative is alive simultaneously.
#' @param mode         \code{"translate"} or \code{"rigid"}.
#' @param require_all  Logical. If \code{TRUE}, the function errors out
#'                     unless EVERY cell in the subset has a matching
#'                     relative in \code{fullCD}.
#'                     If \code{FALSE} (default), unmatched cells are
#'                     dropped from the alignment fit but still moved by
#'                     the same transform.
#'
#' @return A list with elements:
#'   \item{aligned}{The aligned subset data.frame (rotated + translated).}
#'   \item{R}{The 3x3 rotation matrix (identity when
#'           \code{mode = "translate"}).}
#'   \item{t}{The 3-vector translation.}
#'   \item{pairs}{A data.frame of matched (subset_cell, relative_cell,
#'           relative_time) used to fit the transform.}
#'   \item{rmsd}{Root-mean-square distance between aligned subset cells
#'           and their relatives (in input units).}
#'
#' @examples
#' \dontrun{
#'   sub  <- extractSubset(embryo, lineages = "MS", times = 100)
#'   res  <- alignSubsetToRelatives(sub, embryo, relative = "parent",
#'                                  mode = "rigid")
#'   res$rmsd
#' }
alignSubsetToRelatives <- function(subsetCD,
                                   fullCD,
                                   relative    = c("parent", "child"),
                                   subset_time = NULL,
                                   target_time = NULL,
                                   mode        = c("rigid", "translate"),
                                   require_all = FALSE) {
  relative <- match.arg(relative)
  mode     <- match.arg(mode)

  if (!is.data.frame(subsetCD) || nrow(subsetCD) == 0) {
    stop("`subsetCD` must be a non-empty data.frame.")
  }
  if (is.null(subset_time)) {
    ut <- unique(subsetCD$time)
    if (length(ut) != 1) {
      stop("Subset spans multiple times; please specify `subset_time`.")
    }
    subset_time <- ut
  }
  sub <- subsetCD[abs(subsetCD$time - subset_time) < 1e-6, ]
  if (nrow(sub) == 0) {
    stop("No subset rows at subset_time = ", subset_time)
  }
  # Aggregate to one row per cell (handle duplicates conservatively)
  sub_agg <- aggregate(cbind(x, y, z) ~ cell, data = sub, FUN = mean)

  # ---- 1. Per-cell lookup of one relative position -----------------------
  # Helper: for a given relative cell name, return the row of fullCD at the
  # relevant time (latest before subset_time for parent; earliest after for
  # child). Returns NULL if no match.
  lookup_one <- function(rel_name, side) {
    rows <- fullCD[fullCD$cell == rel_name, , drop = FALSE]
    if (nrow(rows) == 0) return(NULL)
    if (!is.null(target_time)) {
      hit <- rows[abs(rows$time - target_time) < 1e-6, , drop = FALSE]
      if (nrow(hit) == 0) return(NULL)
      return(list(time = target_time,
                  x = mean(hit$x), y = mean(hit$y), z = mean(hit$z)))
    }
    if (side == "before") {
      cand <- rows[rows$time < subset_time, , drop = FALSE]
      if (nrow(cand) == 0) return(NULL)
      tt <- max(cand$time)
    } else {  # "after"
      cand <- rows[rows$time > subset_time, , drop = FALSE]
      if (nrow(cand) == 0) return(NULL)
      tt <- min(cand$time)
    }
    hit <- rows[abs(rows$time - tt) < 1e-6, , drop = FALSE]
    list(time = tt,
         x = mean(hit$x), y = mean(hit$y), z = mean(hit$z))
  }

  # ---- 2. Build matched coordinate pairs ---------------------------------
  rows <- list()
  if (relative == "parent") {
    parents_per_cell <- vapply(sub_agg$cell, .getParent, character(1))
    for (i in seq_len(nrow(sub_agg))) {
      pname <- parents_per_cell[i]
      if (is.na(pname)) next
      hit <- lookup_one(pname, "before")
      if (is.null(hit)) next
      rows[[length(rows) + 1L]] <- data.frame(
        cell           = sub_agg$cell[i],
        relative       = pname,
        relative_time  = hit$time,
        x.sub = sub_agg$x[i], y.sub = sub_agg$y[i], z.sub = sub_agg$z[i],
        x.rel = hit$x,        y.rel = hit$y,        z.rel = hit$z,
        stringsAsFactors = FALSE
      )
    }
  } else {  # "child"
    for (i in seq_len(nrow(sub_agg))) {
      kids <- .getChildren(sub_agg$cell[i])
      if (length(kids) == 0) next
      kid_hits <- lapply(kids, lookup_one, side = "after")
      kid_hits <- Filter(Negate(is.null), kid_hits)
      if (length(kid_hits) == 0) next
      hit_x <- mean(vapply(kid_hits, function(h) h$x, numeric(1)))
      hit_y <- mean(vapply(kid_hits, function(h) h$y, numeric(1)))
      hit_z <- mean(vapply(kid_hits, function(h) h$z, numeric(1)))
      kid_t <- vapply(kid_hits, function(h) h$time, numeric(1))
      rows[[length(rows) + 1L]] <- data.frame(
        cell           = sub_agg$cell[i],
        relative       = paste(kids[seq_along(kid_hits)], collapse = "+"),
        relative_time  = mean(kid_t),
        x.sub = sub_agg$x[i], y.sub = sub_agg$y[i], z.sub = sub_agg$z[i],
        x.rel = hit_x,        y.rel = hit_y,        z.rel = hit_z,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) {
    stop("No subset cell had any ", relative,
         " in `fullCD` (subset_time = ", subset_time, ").")
  }
  matched <- do.call(rbind, rows)

  # Informative message about the time spread of relatives we matched to
  rt <- matched$relative_time
  message(sprintf(
    "alignSubsetToRelatives: %d/%d subset cells matched to a %s. ",
    nrow(matched), nrow(sub_agg), relative),
    sprintf("relative_time range: [%g, %g] (median %g)",
            min(rt), max(rt), stats::median(rt)))

  if (require_all && nrow(matched) < nrow(sub_agg)) {
    missing_cells <- setdiff(sub_agg$cell, matched$cell)
    stop(sprintf("require_all = TRUE: %d subset cells had no %s match: %s",
                 length(missing_cells), relative,
                 paste(head(missing_cells, 10), collapse = ", ")))
  }
  if (nrow(matched) < ifelse(mode == "rigid", 3, 1)) {
    stop("Need at least ", ifelse(mode == "rigid", 3, 1),
         " matched pairs; got ", nrow(matched), ".")
  }

  # ---- 3. Fit the rigid transform -----------------------------------------
  A <- as.matrix(matched[, c("x.sub", "y.sub", "z.sub")])
  B <- as.matrix(matched[, c("x.rel", "y.rel", "z.rel")])

  if (mode == "translate") {
    R <- diag(3)
    t <- as.numeric(colMeans(B) - colMeans(A))
  } else {
    pr <- .procrustesRT(A, B)
    R <- pr$R
    t <- pr$t
  }

  # ---- 5. Apply transform to ALL rows of the original subset --------------
  all_coords <- as.matrix(subsetCD[, c("x", "y", "z")])
  new_coords <- sweep(all_coords %*% t(R), 2, t, FUN = "+")
  aligned <- subsetCD
  aligned$x <- new_coords[, 1]
  aligned$y <- new_coords[, 2]
  aligned$z <- new_coords[, 3]

  # ---- 6. Compute fit RMSD on the matched pairs ---------------------------
  fit_A <- sweep(A %*% t(R), 2, t, FUN = "+")
  rmsd  <- sqrt(mean(rowSums((fit_A - B)^2)))

  list(
    aligned = aligned,
    R       = R,
    t       = t,
    pairs   = matched[, c("cell", "relative", "relative_time")],
    rmsd    = rmsd
  )
}


# -----------------------------------------------------------------------------
# Convenience: print summary on source()
# -----------------------------------------------------------------------------
message("LIVEtools_extensions.R loaded. New functions: ",
        "extractSubset(), rotateSubset(), alignSubsetToRelatives().")
