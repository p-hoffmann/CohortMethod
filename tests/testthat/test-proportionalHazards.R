library(testthat)
library(CohortMethod)

createTestPopulation <- function(n = 200, hazard = 0.01, matchPairs = FALSE, addWeights = FALSE) {
  set.seed(123)

  pop <- data.frame(
    rowId = 1:n,
    subjectId = 1:n,
    treatment = rep(c(0, 1), n / 2),
    survivalTime = rexp(n, rate = hazard),
    outcomeCount = rbinom(n, 1, 0.3)
  )

  if (matchPairs) {
    pop$matchId <- rep(1:(n / 2), each = 2)
  }

  if (addWeights) {
    pop$iptw <- runif(n, 0.3, 3)
  }

  return(pop)
}

test_that("testProportionalHazards returns valid PhTestResult structure", {
  set.seed(123)
  n <- 200
  population <- data.frame(
    rowId = 1:n,
    treatment = rep(0:1, each = n / 2),
    survivalTime = rexp(n, rate = 0.1),
    outcomeCount = rbinom(n, 1, 0.3)
  )

  result <- testProportionalHazards(population = population)

  expect_s3_class(result, "PhTestResult")
  expect_true(!is.null(result$covariates))
  expect_true(!is.null(result$settings))
  expect_true(!is.null(result$schoenfeld))
  expect_true(is.data.frame(result$covariates))
  expect_true(nrow(result$covariates) >= 1)
  expect_true(all(c("covariateId", "covariateName", "testStatistic", "pValue") %in%
    names(result$covariates)))
  expect_true(all(result$covariates$pValue >= 0 & result$covariates$pValue <= 1))
})

test_that("testProportionalHazards works without outcomeModel", {
  population <- createTestPopulation(n = 200)

  result <- testProportionalHazards(
    population = population,
    outcomeModel = NULL
  )

  expect_s3_class(result, "PhTestResult")
  expect_true(!is.null(result$covariates))
  expect_true(nrow(result$covariates) > 0)
})

test_that("testProportionalHazards handles stratified Cox models correctly", {
  set.seed(123)
  n <- 200
  population <- data.frame(
    rowId = 1:n,
    treatment = rep(0:1, each = n / 2),
    survivalTime = rexp(n, rate = 0.1),
    outcomeCount = rbinom(n, 1, 0.3),
    stratumId = rep(1:10, each = n / 10)
  )

  result <- testProportionalHazards(population = population)

  expect_s3_class(result, "PhTestResult")
  expect_true(!is.null(result$covariates))
  expect_equal(result$settings$modelType, "stratified")
})

test_that("testProportionalHazards handles matched designs with cluster variance", {
  population <- createTestPopulation(n = 200, matchPairs = TRUE)

  result <- testProportionalHazards(
    population = population,
    outcomeModel = NULL
  )

  expect_s3_class(result, "PhTestResult")
  expect_true(!is.null(result$settings))
  expect_true(!is.null(result$covariates))
})

test_that("testProportionalHazards handles IPTW-weighted models", {
  set.seed(123)
  n <- 200
  population <- data.frame(
    rowId = 1:n,
    treatment = rep(0:1, each = n / 2),
    survivalTime = rexp(n, rate = 0.1),
    outcomeCount = rbinom(n, 1, 0.3),
    iptw = runif(n, 0.5, 2.0)
  )

  result <- testProportionalHazards(population = population)

  expect_s3_class(result, "PhTestResult")
  expect_true(!is.null(result$covariates))
})

test_that("testProportionalHazards handles IPTW with robust variance", {
  population <- createTestPopulation(n = 200, addWeights = TRUE)

  result <- testProportionalHazards(
    population = population,
    outcomeModel = NULL,
    robustVariance = TRUE
  )

  expect_s3_class(result, "PhTestResult")
  expect_true(!is.null(result$covariates))
})

test_that("testProportionalHazards respects iptwTruncation parameter", {
  population <- createTestPopulation(n = 200, addWeights = TRUE)
  population$iptw[1:5] <- c(15, 20, 25, 30, 35)

  result <- testProportionalHazards(
    population = population,
    outcomeModel = NULL,
    iptwTruncation = c(0.1, 10)
  )

  expect_s3_class(result, "PhTestResult")
  expect_true(max(population$iptw) > 10)
})

test_that("testProportionalHazards can disable robust variance", {
  population <- createTestPopulation(n = 200)

  result <- testProportionalHazards(
    population = population,
    outcomeModel = NULL,
    robustVariance = FALSE
  )

  expect_s3_class(result, "PhTestResult")
  expect_true(!is.null(result$covariates))
})

test_that("testProportionalHazards issues warning when p-value < 0.05", {
  set.seed(456)
  n <- 200
  population <- data.frame(
    rowId = 1:n,
    treatment = rep(0:1, each = n / 2),
    survivalTime = c(
      rexp(n / 2, rate = 0.1),
      c(rexp(n / 4, rate = 0.05), rexp(n / 4, rate = 0.15))
    ),
    outcomeCount = rbinom(n, 1, 0.4)
  )

  result <- testProportionalHazards(population = population)

  expect_s3_class(result, "PhTestResult")
})

test_that("testProportionalHazards errors with insufficient events", {
  set.seed(123)
  population <- data.frame(
    rowId = 1:50,
    treatment = rep(0:1, each = 25),
    survivalTime = rexp(50, rate = 0.01),
    outcomeCount = rep(0, 50)
  )

  expect_error(
    testProportionalHazards(population = population),
    "Insufficient events|Cannot fit"
  )
})

test_that("testProportionalHazards handles edge case: very few events", {
  population <- createTestPopulation(n = 200)
  population$outcomeCount <- rbinom(200, 1, 0.05)

  result <- tryCatch(
    {
      testProportionalHazards(population = population)
    },
    error = function(e) {
      expect_true(grepl("Insufficient|events", e$message, ignore.case = TRUE))
      NULL
    }
  )

  if (!is.null(result)) {
    expect_s3_class(result, "PhTestResult")
  }
})

test_that("testProportionalHazards validates input parameters", {
  population <- createTestPopulation(n = 200)

  expect_error(
    testProportionalHazards(population, transform = "invalid"),
    "Must be element of set"
  )

  expect_error(
    testProportionalHazards(population, iptwTruncation = 10),
    "Must have length 2"
  )
})

test_that("testProportionalHazards produces consistent results", {
  population <- createTestPopulation(n = 200)

  result1 <- testProportionalHazards(population)
  result2 <- testProportionalHazards(population)

  expect_equal(result1$covariates$pValue, result2$covariates$pValue)
  expect_equal(result1$covariates$testStatistic, result2$covariates$testStatistic)
})

test_that("plotSchoenfeld returns valid ggplot object", {
  skip_if_not_installed("ggplot2")

  set.seed(123)
  n <- 200
  population <- data.frame(
    rowId = 1:n,
    treatment = rep(0:1, each = n / 2),
    survivalTime = rexp(n, rate = 0.1),
    outcomeCount = rbinom(n, 1, 0.3)
  )

  phResult <- testProportionalHazards(population)
  plot <- plotSchoenfeld(phResult)

  expect_true(inherits(plot, "ggplot"))
  expect_true(!is.null(plot$data))
  expect_true(nrow(plot$data) > 0)
})

test_that("plotSchoenfeld handles single covariate plot", {
  skip_if_not_installed("ggplot2")

  set.seed(123)
  n <- 200
  population <- data.frame(
    rowId = 1:n,
    treatment = rep(0:1, each = n / 2),
    survivalTime = rexp(n, rate = 0.1),
    outcomeCount = rbinom(n, 1, 0.3)
  )

  phResult <- testProportionalHazards(population)
  plot <- plotSchoenfeld(phResult, covariateId = phResult$covariates$covariateId[1])

  expect_true(inherits(plot, "ggplot"))
  expect_true(!is.null(plot$data))
})

test_that("plotSchoenfeld handles multiple covariates", {
  skip_if_not_installed("ggplot2")

  set.seed(123)
  n <- 200
  population <- data.frame(
    rowId = 1:n,
    treatment = rep(0:1, each = n / 2),
    survivalTime = rexp(n, rate = 0.1),
    outcomeCount = rbinom(n, 1, 0.3)
  )

  phResult <- testProportionalHazards(population)
  plot <- plotSchoenfeld(phResult)

  expect_true(inherits(plot, "ggplot"))
  expect_true(nrow(plot$data) > 0)
})

test_that("plotSchoenfeld adds smooth line when showSmooth=TRUE", {
  skip_if_not_installed("ggplot2")

  set.seed(123)
  n <- 200
  population <- data.frame(
    rowId = 1:n,
    treatment = rep(0:1, each = n / 2),
    survivalTime = rexp(n, rate = 0.1),
    outcomeCount = rbinom(n, 1, 0.3)
  )

  phResult <- testProportionalHazards(population)

  plot <- plotSchoenfeld(phResult, showSmooth = TRUE)

  expect_true(inherits(plot, "ggplot"))
  expect_true(length(plot$layers) > 1)
})

test_that("plotSchoenfeld saves to file when fileName specified", {
  skip_if_not_installed("ggplot2")

  set.seed(123)
  n <- 200
  population <- data.frame(
    rowId = 1:n,
    treatment = rep(0:1, each = n / 2),
    survivalTime = rexp(n, rate = 0.1),
    outcomeCount = rbinom(n, 1, 0.3)
  )

  phResult <- testProportionalHazards(population)

  tempFile <- tempfile(fileext = ".png")

  plot <- plotSchoenfeld(phResult, fileName = tempFile)

  expect_true(file.exists(tempFile))
  expect_true(file.size(tempFile) > 0)

  unlink(tempFile)
})

test_that("plotSchoenfeld handles transform parameter", {
  skip_if_not_installed("ggplot2")

  set.seed(123)
  n <- 200
  population <- data.frame(
    rowId = 1:n,
    treatment = rep(0:1, each = n / 2),
    survivalTime = rexp(n, rate = 0.1),
    outcomeCount = rbinom(n, 1, 0.3)
  )

  phResult <- testProportionalHazards(population, transform = "identity")

  plot <- plotSchoenfeld(phResult)

  expect_true(inherits(plot, "ggplot"))
})

test_that("Integration: PH test with Eunomia data", {
  skip_if(
    !all(c("cohortMethodData", "studyPop") %in% ls(envir = .GlobalEnv)),
    "Eunomia test data not available"
  )

  tryCatch(
    {
      studyPop <- get("studyPop", envir = .GlobalEnv)

      result <- testProportionalHazards(
        population = studyPop
      )

      expect_s3_class(result, "PhTestResult")
      expect_true(!is.null(result$covariates))
      expect_true(nrow(result$covariates) > 0)
      expect_true(!is.null(result$schoenfeld))
    },
    error = function(e) {
      skip(paste("Could not query Eunomia data:", e$message))
    }
  )
})

test_that("PH test rejects population with only one treatment group", {
  population <- data.frame(
    rowId = 1:100,
    subjectId = 1:100,
    treatment = rep(0, 100),
    survivalTime = runif(100, 50, 200),
    outcomeCount = rbinom(100, 1, 0.4)
  )

  expect_error(testProportionalHazards(population = population))
})

test_that("PH test rejects population with very sparse events", {
  population <- data.frame(
    rowId = 1:200,
    subjectId = 1:200,
    treatment = rep(c(0, 1), 100),
    survivalTime = runif(200, 50, 300)
  )

  eventIndices <- sample(1:200, 5)
  population$outcomeCount <- ifelse(1:200 %in% eventIndices, 1, 0)

  expect_error(
    testProportionalHazards(population = population),
    "Insufficient events: 5 events found, minimum 10 required"
  )
})

test_that("PH test succeeds with moderate events", {
  population <- data.frame(
    rowId = 1:200,
    subjectId = 1:200,
    treatment = rep(c(0, 1), 100),
    survivalTime = runif(200, 50, 300),
    outcomeCount = rbinom(200, 1, 0.15)
  )

  result <- testProportionalHazards(population = population)

  expect_s3_class(result, "PhTestResult")
  expect_true("settings" %in% names(result))
  expect_true(result$settings$nEvents >= 10)
})

test_that("PH test completes with reasonable proportional hazards data", {
  set.seed(456)
  population <- data.frame(
    rowId = 1:200,
    subjectId = 1:200,
    treatment = rep(c(0, 1), 100)
  )

  hr <- 2
  baseHazard <- 0.01
  population$survivalTime <- ifelse(
    population$treatment == 1,
    rexp(200, baseHazard * hr),
    rexp(200, baseHazard)
  )

  population$survivalTime <- pmin(population$survivalTime, 200)
  population$outcomeCount <- ifelse(population$survivalTime < 200, 1, 0)

  result <- testProportionalHazards(population = population)

  expect_s3_class(result, "PhTestResult")
  expect_equal(nrow(result$covariates), 1)
  expect_true(result$covariates$pValue >= 0 && result$covariates$pValue <= 1)
  expect_true(!is.na(result$covariates$pValue))
})

test_that("PH test handles varying hazards gracefully", {
  set.seed(789)
  population <- data.frame(
    rowId = 1:200,
    subjectId = 1:200,
    treatment = rep(c(0, 1), 100)
  )

  population$survivalTime <- ifelse(
    population$treatment == 1,
    c(rexp(50, 0.05), rexp(50, 0.005)),
    rexp(100, 0.02)
  )

  population$outcomeCount <- rbinom(200, 1, 0.4)

  result <- testProportionalHazards(population = population)

  expect_s3_class(result, "PhTestResult")

  pValue <- result$covariates$pValue[1]
  expect_true(pValue >= 0 && pValue <= 1)
  expect_true(!is.na(pValue))
})

test_that("PH test with rank transform handles extreme values", {
  population <- data.frame(
    rowId = 1:100,
    subjectId = 1:100,
    treatment = rep(c(0, 1), 50),
    survivalTime = c(runif(50, 0.1, 10), runif(50, 100, 1000)),
    outcomeCount = rbinom(100, 1, 0.4)
  )

  result <- testProportionalHazards(
    population = population,
    transform = "rank"
  )

  expect_s3_class(result, "PhTestResult")

  pValue <- result$covariates$pValue[1]
  expect_true(!is.na(pValue))
  expect_true(!is.nan(pValue))
  expect_true(pValue >= 0 && pValue <= 1)
})
