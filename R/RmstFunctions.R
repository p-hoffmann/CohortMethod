# Copyright 2025 Observational Health Data Sciences and Informatics
#
# This file is part of CohortMethod
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Calculate Restricted Mean Survival Time (RMST)
#'
#' @description
#' Calculate RMST for treatment groups with appropriate variance estimation for matched,
#' stratified, or weighted designs.
#'
#' @details
#' Automatically detects study design and applies the appropriate variance estimation:
#' matched (pseudovalue-based paired differences), stratified (stratified variance),
#' or IPTW (weighted Kaplan-Meier).
#'
#' @param population            Population from [createStudyPopulation()]. Required columns:
#'                              rowId, treatment, survivalTime, outcomeCount. Optional: stratumId,
#'                              matchId, iptw.
#' @param tau                   Time horizon(s) for RMST. If NULL, derived using tauMethod.
#' @param tauMethod             Auto tau selection: "max", "median", or "quantile". Default "max".
#' @param tauQuantile           Quantile for tauMethod = "quantile". Default 0.75.
#' @param confLevel             Confidence level. Default 0.95.
#' @param returnCurves          Return survival curves for plotting. Default FALSE.
#'
#' @return
#' An object of class `RmstResult` containing:
#' \describe{
#'   \item{tau}{The time horizon(s) used}
#'   \item{rmst}{A data frame with treatment-specific RMST estimates, SE, and CI}
#'   \item{contrast}{A data frame with RMST differences between treatment groups}
#'   \item{settings}{A list with calculation settings and metadata}
#'   \item{survivalCurves}{Optional survival curves (if returnCurves = TRUE)}
#' }
#'
#' @examples
#' \dontrun{
#' rmst <- computeRmst(population, tau = 365)
#' rmst <- computeRmst(population, tau = c(90, 180, 365))
#' }
#'
#' @export
computeRmst <- function(population,
                        tau = NULL,
                        tauMethod = "max",
                        tauQuantile = 0.75,
                        confLevel = 0.95,
                        returnCurves = FALSE) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(population, add = errorMessages)
  checkmate::assertNumeric(tau,
    lower = 0, finite = TRUE, any.missing = FALSE,
    min.len = 1, null.ok = TRUE, add = errorMessages
  )
  checkmate::assertChoice(tauMethod, c("max", "median", "quantile"), add = errorMessages)
  checkmate::assertNumber(tauQuantile, lower = 0, upper = 1, add = errorMessages)
  checkmate::assertNumber(confLevel, lower = 0, upper = 1, add = errorMessages)
  checkmate::assertLogical(returnCurves, len = 1, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)

  alpha <- 1 - confLevel

  requiredFields <- c("rowId", "treatment", "survivalTime", "outcomeCount")
  if (!all(requiredFields %in% names(population))) {
    stop(sprintf(
      "Population must contain fields: %s",
      paste(requiredFields, collapse = ", ")
    ))
  }

  if (any(is.na(population$survivalTime)) || any(is.na(population$outcomeCount)) ||
    any(is.na(population$treatment))) {
    stop("Missing values detected in required fields")
  }

  hasMatching <- "matchId" %in% names(population)
  hasStratification <- "stratumId" %in% names(population)
  hasIPTW <- "iptw" %in% names(population)

  nAdjustments <- sum(hasMatching, hasIPTW)
  if (nAdjustments > 1) {
    methods <- c()
    if (hasMatching) methods <- c(methods, "matching (matchId)")
    if (hasIPTW) methods <- c(methods, "IPTW (iptw)")
    stop(sprintf(
      "Multiple adjustment methods detected: %s. Population should contain only one adjustment method.",
      paste(methods, collapse = " and ")
    ))
  }

  adjustmentMethod <- "none"
  if (hasIPTW) {
    adjustmentMethod <- "iptw"
    population$iptw <- .handleExtremeWeights(population$iptw)
  } else if (hasMatching) {
    adjustmentMethod <- "matching"
  } else if (hasStratification) {
    adjustmentMethod <- "stratification"
  }

  warnings <- character()

  nByTreatment <- table(population$treatment)
  if (any(nByTreatment < 20)) {
    stop("Insufficient subjects: at least 20 subjects required per treatment group")
  }

  eventsByTreatment <- aggregate(population$outcomeCount,
    by = list(treatment = population$treatment),
    FUN = sum
  )
  if (any(eventsByTreatment$x < 5)) {
    msg <- "Fewer than 5 events in at least one treatment group. RMST estimates may be unreliable."
    warning(msg)
    warnings <- c(warnings, msg)
  }

  maxFollowUp <- max(population$survivalTime)
  medianFollowUp <- median(population$survivalTime)

  if (maxFollowUp < 30) {
    msg <- sprintf("Very short follow-up detected (max %.1f days). RMST may not be meaningful.", maxFollowUp)
    warning(msg)
    warnings <- c(warnings, msg)
  } else if (medianFollowUp < 30) {
    msg <- sprintf("Short median follow-up (%.1f days). Consider longer observation period.", medianFollowUp)
    warning(msg)
    warnings <- c(warnings, msg)
  }

  if (is.null(tau)) {
    tau <- switch(tauMethod,
      "max" = maxFollowUp,
      "median" = median(population$survivalTime),
      "quantile" = quantile(population$survivalTime, tauQuantile, names = FALSE),
      stop("Invalid tauMethod: ", tauMethod)
    )
    ParallelLogger::logInfo(sprintf("Derived tau = %.2f using method '%s'", tau, tauMethod))
  }

  tauResult <- .validateTau(tau, maxFollowUp)
  tau <- tauResult$tau
  tauTruncated <- tauResult$tauTruncated

  if (any(tauTruncated)) {
    ParallelLogger::logWarn(sprintf(
      "Tau truncated for %d value(s) exceeding maximum follow-up (%.2f days)",
      sum(tauTruncated), maxFollowUp
    ))
  }

  ParallelLogger::logInfo(sprintf(
    "Computing RMST for %d tau value(s) with %s adjustment (N=%d, events=%d)",
    length(tau), adjustmentMethod, nrow(population), sum(population$outcomeCount)
  ))

  rmstEstimates <- data.frame()
  contrastEstimates <- data.frame()
  survivalCurves <- if (returnCurves) list() else NULL

  kmFits <- list()

  for (trt in c(0, 1)) {
    popSubset <- population[population$treatment == trt, ]

    if (hasIPTW) {
      kmFits[[as.character(trt)]] <- survival::survfit(
        survival::Surv(survivalTime, outcomeCount) ~ 1,
        data = popSubset,
        weights = popSubset$iptw
      )
    } else {
      kmFits[[as.character(trt)]] <- survival::survfit(
        survival::Surv(survivalTime, outcomeCount) ~ 1,
        data = popSubset
      )
    }
  }

  for (t in tau) {
    for (trt in c(0, 1)) {
      kmFit <- kmFits[[as.character(trt)]]

      rmst <- tryCatch(
        {
          .integrateKmCurve(kmFit, t)
        },
        error = function(e) {
          warning(sprintf(
            "Failed to calculate RMST for treatment=%d, tau=%.2f: %s",
            trt, t, e$message
          ))
          return(NA)
        }
      )

      se <- tryCatch(
        {
          .computeRmstVariance(kmFit, t)
        },
        error = function(e) {
          warning(sprintf(
            "Failed to calculate RMST variance for treatment=%d, tau=%.2f: %s",
            trt, t, e$message
          ))
          return(NA)
        }
      )

      z <- qnorm(1 - alpha / 2)
      lower <- rmst - z * se
      upper <- rmst + z * se

      rmstEstimates <- rbind(rmstEstimates, data.frame(
        tau = t,
        treatment = trt,
        rmst = rmst,
        se = se,
        lower95ci = lower,
        upper95ci = upper
      ))

      if (returnCurves) {
        trtLabel <- if (trt == 0) "comparator" else "target"
        survivalCurves[[paste0(trtLabel, "_tau", t)]] <- list(
          time = kmFit$time,
          surv = kmFit$surv,
          lower = kmFit$lower,
          upper = kmFit$upper
        )
      }
    }

    rmst0 <- rmstEstimates$rmst[rmstEstimates$tau == t & rmstEstimates$treatment == 0]
    rmst1 <- rmstEstimates$rmst[rmstEstimates$tau == t & rmstEstimates$treatment == 1]
    se0 <- rmstEstimates$se[rmstEstimates$tau == t & rmstEstimates$treatment == 0]
    se1 <- rmstEstimates$se[rmstEstimates$tau == t & rmstEstimates$treatment == 1]

    rmstDiff <- rmst1 - rmst0

    if (adjustmentMethod == "matching") {
      pairedResult <- .computeMatchedRmstContrast(population, t, kmFits)

      if (pairedResult$method == "failed" || is.na(pairedResult$diff)) {
        msg <- sprintf(
          "Matched RMST pseudovalue computation failed for tau=%.2f. Falling back to independent-group variance. %s",
          t,
          if (!is.null(pairedResult$error)) paste("Error:", pairedResult$error) else ""
        )
        warning(msg)
        warnings <- c(warnings, msg)
        ParallelLogger::logWarn(msg)

        rmstDiff <- rmst1 - rmst0
        seDiff <- sqrt(se0^2 + se1^2)

        if (any(population$outcomeCount == 0 & population$survivalTime < t)) {
          censorMsg <- sprintf(
            "Censoring detected before tau=%.2f. Independent variance may be conservative. Consider investigating censoring patterns.",
            t
          )
          warning(censorMsg)
          warnings <- c(warnings, censorMsg)
          ParallelLogger::logWarn(censorMsg)
        }
      } else if (pairedResult$nPairs < 10) {
        rmstDiff <- pairedResult$diff
        seDiff <- pairedResult$se

        msg <- sprintf(
          "Matched RMST based on only %d pairs for tau=%.2f. Estimates may be unstable.",
          pairedResult$nPairs, t
        )
        warning(msg)
        warnings <- c(warnings, msg)
      } else {
        rmstDiff <- pairedResult$diff
        seDiff <- pairedResult$se

        ParallelLogger::logInfo(sprintf(
          "Matched RMST computed using pseudovalues (%d pairs) for tau=%.2f",
          pairedResult$nPairs, t
        ))
      }
    } else {
      seDiff <- sqrt(se0^2 + se1^2)
    }

    z <- qnorm(1 - alpha / 2)
    lowerDiff <- rmstDiff - z * seDiff
    upperDiff <- rmstDiff + z * seDiff
    pValue <- 2 * (1 - pnorm(abs(rmstDiff / seDiff)))

    contrastEstimates <- rbind(contrastEstimates, data.frame(
      tau = t,
      rmstDiff = rmstDiff,
      se = seDiff,
      lower95ci = lowerDiff,
      upper95ci = upperDiff,
      pValue = pValue
    ))
  }

  varianceMethod <- if (adjustmentMethod == "matching") {
    "paired"
  } else if (adjustmentMethod == "stratification") {
    "stratified"
  } else if (adjustmentMethod == "iptw") {
    "weighted"
  } else {
    "greenwood"
  }

  design <- switch(adjustmentMethod,
    "matching" = "matched",
    "stratification" = "stratified",
    "iptw" = "weighted",
    "none" = "none"
  )

  result <- list(
    tau = tau,
    rmst = rmstEstimates,
    contrast = contrastEstimates,
    design = design,
    settings = list(
      alpha = alpha,
      nSubjects = nrow(population),
      nEvents = sum(population$outcomeCount),
      maxFollowUp = maxFollowUp,
      medianFollowUp = medianFollowUp,
      adjustmentMethod = adjustmentMethod,
      nStrata = if (hasStratification) length(unique(population$stratumId)) else NA,
      nPairs = if (hasMatching) length(unique(population$matchId)) else NA,
      tauTruncated = tauTruncated,
      varianceMethod = varianceMethod,
      warnings = if (length(warnings) > 0) warnings else NULL
    ),
    survivalCurves = survivalCurves
  )

  class(result) <- "RmstResult"
  return(result)
}

#' Print method for RmstResult
#'
#' @param x An object of class `RmstResult`.
#' @param ... Additional arguments (not used).
#'
#' @export
print.RmstResult <- function(x, ...) {
  cat("Restricted Mean Survival Time (RMST) Results\n")
  cat("=============================================\n\n")

  cat(sprintf(
    "Sample size: %d subjects, %d events\n",
    x$settings$nSubjects, x$settings$nEvents
  ))
  cat(sprintf("Adjustment method: %s\n", x$settings$adjustmentMethod))
  cat(sprintf("Variance method: %s\n\n", x$settings$varianceMethod))

  for (t in unique(x$rmst$tau)) {
    cat(sprintf("Time horizon (tau) = %.1f:\n", t))

    rmst_tau <- x$rmst[x$rmst$tau == t, ]
    for (i in seq_len(nrow(rmst_tau))) {
      trt <- rmst_tau$treatment[i]
      trtLabel <- if (trt == 0) "Comparator" else "Target    "
      cat(sprintf(
        "  %s: RMST = %6.2f (95%% CI: %6.2f to %6.2f)\n",
        trtLabel, rmst_tau$rmst[i],
        rmst_tau$lower95ci[i], rmst_tau$upper95ci[i]
      ))
    }

    contrast_tau <- x$contrast[x$contrast$tau == t, ]
    cat(sprintf(
      "  Difference: %6.2f (95%% CI: %6.2f to %6.2f), p-value = %.4f\n\n",
      contrast_tau$rmstDiff,
      contrast_tau$lower95ci, contrast_tau$upper95ci,
      contrast_tau$pValue
    ))
  }

  if (x$settings$tauTruncated) {
    cat("Note: One or more tau values were truncated to maximum follow-up time\n")
  }

  invisible(x)
}

.validateTau <- function(tau, maxFollowUp) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertNumeric(tau,
    lower = 0, finite = TRUE, any.missing = FALSE,
    min.len = 1, add = errorMessages
  )
  checkmate::reportAssertions(collection = errorMessages)

  tauTruncated <- any(tau > maxFollowUp)

  if (tauTruncated) {
    originalTau <- tau
    tau <- pmin(tau, maxFollowUp)
    warning(sprintf(
      "tau exceeds maximum follow-up time (%.2f). Truncating to maximum: %s",
      maxFollowUp, paste(sprintf("%.2f", tau), collapse = ", ")
    ))
  }

  if (any(duplicated(tau))) {
    tau <- unique(tau)
    warning("Duplicate tau values after truncation. Using unique values only.")
  }

  return(list(
    tau = tau,
    tauTruncated = tauTruncated
  ))
}

.handleExtremeWeights <- function(weights, minWeight = 0.1, maxWeight = 10) {
  if (is.null(weights)) {
    return(NULL)
  }

  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertNumeric(weights,
    lower = 0, finite = TRUE, any.missing = FALSE,
    add = errorMessages
  )
  checkmate::reportAssertions(collection = errorMessages)

  originalWeights <- weights
  needsTruncation <- any(weights < minWeight | weights > maxWeight)

  if (needsTruncation) {
    nLow <- sum(weights < minWeight)
    nHigh <- sum(weights > maxWeight)

    weights <- pmax(pmin(weights, maxWeight), minWeight)

    msg <- sprintf(
      paste(
        "Extreme IPTW weights detected:",
        "%d weights < %.2f, %d weights > %.2f.",
        "Truncating to [%.2f, %.2f]"
      ),
      nLow, minWeight, nHigh, maxWeight, minWeight, maxWeight
    )
    warning(msg)
    ParallelLogger::logWarn(msg)
  }

  return(weights)
}

# RMST variance using influence function (Uno et al. 2014)
.computeRmstVariance <- function(kmFit, tau) {
  idx <- kmFit$time <= tau
  wk.time <- sort(c(kmFit$time[idx], tau))
  wk.surv <- kmFit$surv[idx]
  wk.n.risk <- kmFit$n.risk[idx]
  wk.n.event <- kmFit$n.event[idx]

  if (length(wk.time) <= 1) {
    return(0)
  }

  time.diff <- diff(c(0, wk.time))
  areas <- time.diff * c(1, wk.surv)

  # Greenwood variance
  wk.var <- ifelse((wk.n.risk - wk.n.event) == 0, 0,
    wk.n.event / (wk.n.risk * (wk.n.risk - wk.n.event))
  )
  wk.var <- c(wk.var, 0)

  cum_rev_areas <- cumsum(rev(areas[-1]))
  rmst.var <- sum(cum_rev_areas^2 * rev(wk.var)[-1])

  return(sqrt(rmst.var))
}

.computeRmstPseudovalues <- function(kmFit, population, tau) {
  rmst_all <- .integrateKmCurve(kmFit, tau)

  pseudovalues <- computeRmstPseudovaluesInternal(
    subjectTimes = population$survivalTime,
    subjectEvents = as.integer(population$outcomeCount),
    kmTimes = kmFit$time,
    kmSurv = kmFit$surv,
    kmNRisk = as.integer(kmFit$n.risk),
    kmNEvent = as.integer(kmFit$n.event),
    tau = tau,
    rmstAll = rmst_all
  )

  return(data.frame(
    rowId = population$rowId,
    pseudovalue = pseudovalues
  ))
}

.integrateKmCurve <- function(kmFit, tau) {
  times <- c(0, kmFit$time[kmFit$time <= tau])
  surv <- c(1, kmFit$surv[kmFit$time <= tau])

  if (max(times) < tau) {
    times <- c(times, tau)
    surv <- c(surv, tail(surv, 1))
  }

  uniqueIdx <- !duplicated(times)
  times <- times[uniqueIdx]
  surv <- surv[uniqueIdx]

  time.diff <- diff(times)
  areas <- time.diff * surv[-length(surv)]

  return(sum(areas))
}

.computeMatchedRmstContrast <- function(population, tau, kmFits) {
  if (!"matchId" %in% names(population)) {
    stop("matchId column not found in population for matched RMST")
  }

  matchIds <- unique(population$matchId[!is.na(population$matchId)])

  if (length(matchIds) == 0) {
    stop("No valid match pairs found in population")
  }

  censoringBeforeTau <- any(population$outcomeCount == 0 & population$survivalTime < tau)

  tryCatch(
    {
      pop_target <- population[population$treatment == 1, ]
      pv_target <- .computeRmstPseudovalues(kmFits[["1"]], pop_target, tau)

      pop_comparator <- population[population$treatment == 0, ]
      pv_comparator <- .computeRmstPseudovalues(kmFits[["0"]], pop_comparator, tau)

      population <- merge(population, pv_target, by = "rowId", all.x = TRUE, suffixes = c("", "_target"))
      population <- merge(population, pv_comparator, by = "rowId", all.x = TRUE, suffixes = c("", "_comparator"))

      population$pseudovalue <- ifelse(population$treatment == 1,
        population$pseudovalue,
        population$pseudovalue_comparator
      )
      population$pseudovalue_comparator <- NULL

      pairStats <- population %>%
        filter(!is.na(.data$matchId), !is.na(.data$pseudovalue)) %>%
        group_by(.data$matchId) %>%
        summarise(
          nSubjects = n(),
          hasBothArms = all(c(0, 1) %in% .data$treatment),
          diff = if (n() == 2 && all(c(0, 1) %in% .data$treatment)) {
            .data$pseudovalue[.data$treatment == 1] - .data$pseudovalue[.data$treatment == 0]
          } else {
            NA_real_
          },
          .groups = "drop"
        )

      skippedWrongSize <- sum(pairStats$nSubjects != 2)
      skippedMissingArm <- sum(pairStats$nSubjects == 2 & !pairStats$hasBothArms)
      totalSkipped <- skippedWrongSize + skippedMissingArm

      if (totalSkipped > 0) {
        skipReasons <- c()
        if (skippedWrongSize > 0) {
          skipReasons <- c(skipReasons, sprintf("%d wrong size", skippedWrongSize))
        }
        if (skippedMissingArm > 0) {
          skipReasons <- c(skipReasons, sprintf("%d missing treatment arm", skippedMissingArm))
        }
        warning(sprintf(
          "Skipped %d of %d matched pairs: %s",
          totalSkipped, length(matchIds), paste(skipReasons, collapse = ", ")
        ))
      }

      pairDifferences <- pairStats$diff[!is.na(pairStats$diff)]

      if (length(pairDifferences) == 0) {
        stop("No valid matched pairs found after pseudovalue computation")
      }

      if (length(pairDifferences) < 10) {
        msg <- sprintf(
          "Only %d valid matched pairs available. Consider using independent-group RMST for more stable estimates.",
          length(pairDifferences)
        )
        warning(msg)
        ParallelLogger::logWarn(msg)
      }

      meanDiff <- mean(pairDifferences)
      seDiff <- sd(pairDifferences) / sqrt(length(pairDifferences))

      return(list(
        diff = meanDiff,
        se = seDiff,
        nPairs = length(pairDifferences),
        method = "pseudovalue",
        censoringBeforeTau = censoringBeforeTau
      ))
    },
    error = function(e) {
      msg <- sprintf("Matched RMST pseudovalue computation failed: %s", e$message)
      warning(msg)
      ParallelLogger::logWarn(msg)

      return(list(
        diff = NA,
        se = NA,
        nPairs = 0,
        method = "failed",
        error = e$message
      ))
    }
  )
}
