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

#' Create diagnostic plots for proportional hazards assumption
#'
#' @description
#' Generate diagnostic plots showing Schoenfeld residuals over time to visually assess
#' proportional hazards assumption violations.
#'
#' @details
#' This function creates a plot of Schoenfeld residuals versus time for a specified
#' covariate from a proportional hazards test. A non-constant pattern in the residuals
#' indicates violation of the proportional hazards assumption for that covariate.
#'
#' The plot includes:
#' - Scatter plot of Schoenfeld residuals over time
#' - Reference line at y=0
#' - Optional smooth line (LOWESS) to show trend
#' - Optional confidence bands
#' - Test statistics and p-value in subtitle
#'
#' @param phTestResult          An object of class `PhTestResult` from [testProportionalHazards()].
#' @param covariateId           The covariate ID to plot. If NULL, plots the treatment effect.
#'                              Default is NULL.
#' @param showSmooth            Whether to add a LOWESS smooth line. Default is TRUE.
#' @param showConfidenceBand    Whether to add ±2 SE confidence band. Default is FALSE.
#' @param fileName              Optional file name to save the plot. If NULL, plot is returned
#'                              but not saved.
#' @param width                 Width of saved plot in inches. Default is 8.
#' @param height                Height of saved plot in inches. Default is 6.
#'
#' @return
#' A ggplot2 plot object. If fileName is provided, the plot is also saved to file.
#'
#' @examples
#' \dontrun{
#' # After testing proportional hazards:
#' phTest <- testProportionalHazards(
#'   population = population,
#'   outcomeModel = outcomeModel
#' )
#'
#' # Create diagnostic plot
#' plot <- plotSchoenfeld(
#'   phTestResult = phTest,
#'   showSmooth = TRUE
#' )
#'
#' print(plot)
#'
#' # Save to file
#' plotSchoenfeld(
#'   phTestResult = phTest,
#'   fileName = "ph_diagnostics.png",
#'   width = 10,
#'   height = 6
#' )
#' }
#'
#' @references
#' Grambsch, P. M., & Therneau, T. M. (1994). Proportional hazards tests and diagnostics
#' based on weighted residuals. Biometrika, 81(3), 515-526.
#'
#' @export
plotSchoenfeld <- function(phTestResult,
                           covariateId = NULL,
                           showSmooth = TRUE,
                           showConfidenceBand = FALSE,
                           fileName = NULL,
                           width = 8,
                           height = 6) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertClass(phTestResult, "PhTestResult", add = errorMessages)
  checkmate::assertNumeric(covariateId, null.ok = TRUE, len = 1, add = errorMessages)
  checkmate::assertLogical(showSmooth, len = 1, add = errorMessages)
  checkmate::assertLogical(showConfidenceBand, len = 1, add = errorMessages)
  checkmate::assertCharacter(fileName, null.ok = TRUE, len = 1, add = errorMessages)
  checkmate::assertNumber(width, lower = 1, add = errorMessages)
  checkmate::assertNumber(height, lower = 1, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)

  if (is.null(phTestResult$schoenfeld) || is.null(phTestResult$schoenfeld$residuals)) {
    stop("No Schoenfeld residuals available in PhTestResult object")
  }

  if (is.null(covariateId)) {
    covIdx <- 1
    covName <- "Treatment"
  } else {
    covRow <- phTestResult$covariates[phTestResult$covariates$covariateId == covariateId, ]
    if (nrow(covRow) == 0) {
      stop(sprintf("Covariate ID %d not found in PhTestResult", covariateId))
    }
    covName <- covRow$covariateName
    covIdx <- which(phTestResult$covariates$covariateId == covariateId)
  }

  time <- phTestResult$schoenfeld$time
  residuals <- phTestResult$schoenfeld$residuals

  if (is.matrix(residuals)) {
    if (covIdx > ncol(residuals)) {
      stop(sprintf(
        "Covariate index %d exceeds number of columns in residuals (%d)",
        covIdx, ncol(residuals)
      ))
    }
    resid <- residuals[, covIdx]
  } else {
    resid <- residuals
  }

  if (length(resid) < 10) {
    warning("Fewer than 10 residuals available. Plot may not be informative.")
  }

  plotData <- data.frame(
    time = time,
    residual = resid
  )

  covStats <- phTestResult$covariates[covIdx, ]

  p <- ggplot2::ggplot(plotData, ggplot2::aes(x = time, y = residual)) +
    ggplot2::geom_point(alpha = 0.5, size = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::labs(
      title = sprintf("Schoenfeld Residuals: %s", covName),
      subtitle = sprintf(
        "Chi-square = %.2f, p-value = %.4f%s",
        covStats$testStatistic,
        covStats$pValue,
        if (covStats$violated) " *PH VIOLATION*" else ""
      ),
      x = "Time",
      y = "Schoenfeld Residual"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10)
    )

  if (showSmooth) {
    p <- p + ggplot2::geom_smooth(
      method = "loess",
      se = showConfidenceBand,
      color = "blue",
      fill = "lightblue",
      alpha = 0.3
    )
  }

  if (showConfidenceBand && !showSmooth) {
    se <- sd(resid) / sqrt(length(resid))
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = 0 - 2 * se, ymax = 0 + 2 * se),
        alpha = 0.2,
        fill = "gray"
      )
  }

  if (!is.null(fileName)) {
    ggplot2::ggsave(
      filename = fileName,
      plot = p,
      width = width,
      height = height,
      dpi = 300
    )
    message(sprintf("Plot saved to: %s", fileName))
  }

  return(p)
}
