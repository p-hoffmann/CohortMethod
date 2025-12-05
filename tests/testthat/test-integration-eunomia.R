library(testthat)
library(CohortMethod)

skip_on_cran()

test_that("Eunomia integration test: PH testing with real data", {
  skip_if_not_installed("Eunomia")
  skip_if_not_installed("DatabaseConnector")

  library(Eunomia)
  library(DatabaseConnector)

  connectionDetails <- tryCatch(
    {
      getEunomiaConnectionDetails()
    },
    error = function(e) {
      skip("Eunomia database not available")
    }
  )

  connection <- connect(connectionDetails)
  on.exit(disconnect(connection), add = TRUE)

  sql <- "
  SELECT
    cohort_definition_id AS treatment,
    subject_id,
    cohort_start_date,
    cohort_end_date
  FROM main.cohort
  WHERE cohort_definition_id IN (1, 2)
  LIMIT 200
  "

  cohortData <- querySql(connection, sql, snakeCaseToCamelCase = TRUE)

  expect_gt(nrow(cohortData), 0)
  expect_true(all(c("treatment", "subjectId", "cohortStartDate", "cohortEndDate") %in% colnames(cohortData)))

  outcomeSql <- "
  SELECT
    person_id AS subject_id,
    condition_start_date AS outcome_date
  FROM main.condition_occurrence
  WHERE condition_concept_id = 192671
  "

  outcomeData <- querySql(connection, outcomeSql, snakeCaseToCamelCase = TRUE)

  population <- cohortData %>%
    dplyr::left_join(outcomeData, by = "subjectId") %>%
    dplyr::mutate(
      outcomeCount = ifelse(is.na(outcomeDate), 0, 1),
      survivalTime = ifelse(
        is.na(outcomeDate),
        as.numeric(difftime(cohortEndDate, cohortStartDate, units = "days")),
        as.numeric(difftime(outcomeDate, cohortStartDate, units = "days"))
      ),
      treatment = ifelse(treatment == 1, 1, 0)
    ) %>%
    dplyr::filter(survivalTime > 0) %>%
    dplyr::select(subjectId, treatment, outcomeCount, survivalTime) %>%
    as.data.frame()

  population$rowId <- seq_len(nrow(population))

  phResult <- testProportionalHazards(population)

  expect_s3_class(phResult, "PhTestResult")
  expect_true(!is.null(phResult$global))
  expect_true(!is.null(phResult$covariates))
  expect_true("pValue" %in% names(phResult$global))
  expect_true(phResult$global$pValue >= 0 && phResult$global$pValue <= 1)
})

test_that("Eunomia integration test: RMST computation with real data", {
  skip_if_not_installed("Eunomia")
  skip_if_not_installed("DatabaseConnector")

  library(Eunomia)
  library(DatabaseConnector)

  connectionDetails <- tryCatch(
    {
      getEunomiaConnectionDetails()
    },
    error = function(e) {
      skip("Eunomia database not available")
    }
  )

  connection <- connect(connectionDetails)
  on.exit(disconnect(connection), add = TRUE)

  sql <- "
  SELECT
    cohort_definition_id AS treatment,
    subject_id,
    cohort_start_date,
    cohort_end_date
  FROM main.cohort
  WHERE cohort_definition_id IN (1, 2)
  LIMIT 200
  "

  cohortData <- querySql(connection, sql, snakeCaseToCamelCase = TRUE)

  outcomeSql <- "
  SELECT
    person_id AS subject_id,
    condition_start_date AS outcome_date
  FROM main.condition_occurrence
  WHERE condition_concept_id = 192671
  "

  outcomeData <- querySql(connection, outcomeSql, snakeCaseToCamelCase = TRUE)

  population <- cohortData %>%
    dplyr::left_join(outcomeData, by = "subjectId") %>%
    dplyr::mutate(
      outcomeCount = ifelse(is.na(outcomeDate), 0, 1),
      survivalTime = ifelse(
        is.na(outcomeDate),
        as.numeric(difftime(cohortEndDate, cohortStartDate, units = "days")),
        as.numeric(difftime(outcomeDate, cohortStartDate, units = "days"))
      ),
      treatment = ifelse(treatment == 1, 1, 0)
    ) %>%
    dplyr::filter(survivalTime > 0) %>%
    dplyr::select(subjectId, treatment, outcomeCount, survivalTime) %>%
    as.data.frame()

  population$rowId <- seq_len(nrow(population))

  rmstResult <- computeRmst(
    population = population,
    tau = NULL,
    tauMethod = "median"
  )

  expect_s3_class(rmstResult, "RmstResult")
  expect_true(!is.null(rmstResult$rmst))
  expect_true(!is.null(rmstResult$contrast))
  expect_equal(nrow(rmstResult$rmst), 2)
  expect_true(all(c("treatment", "rmst", "se", "lower", "upper") %in% colnames(rmstResult$rmst)))
  expect_true(rmstResult$settings$adjustmentMethod == "none")

  maxFollow <- max(population$survivalTime)
  rmstMulti <- computeRmst(
    population = population,
    tau = c(30, 90, 180)
  )

  expect_equal(length(rmstMulti$tau), 3)
  expect_equal(nrow(rmstMulti$rmst), 6)
})

test_that("Eunomia integration test: End-to-end workflow", {
  skip_if_not_installed("Eunomia")
  skip_if_not_installed("DatabaseConnector")

  library(Eunomia)
  library(DatabaseConnector)

  connectionDetails <- tryCatch(
    {
      getEunomiaConnectionDetails()
    },
    error = function(e) {
      skip("Eunomia database not available")
    }
  )

  connection <- connect(connectionDetails)
  on.exit(disconnect(connection), add = TRUE)

  sql <- "
  SELECT
    cohort_definition_id AS treatment,
    subject_id,
    cohort_start_date,
    cohort_end_date
  FROM main.cohort
  WHERE cohort_definition_id IN (1, 2)
  LIMIT 100
  "

  cohortData <- querySql(connection, sql, snakeCaseToCamelCase = TRUE)

  outcomeSql <- "
  SELECT
    person_id AS subject_id,
    condition_start_date AS outcome_date
  FROM main.condition_occurrence
  WHERE condition_concept_id = 192671
  "

  outcomeData <- querySql(connection, outcomeSql, snakeCaseToCamelCase = TRUE)

  population <- cohortData %>%
    dplyr::left_join(outcomeData, by = "subjectId") %>%
    dplyr::mutate(
      outcomeCount = ifelse(is.na(outcomeDate), 0, 1),
      survivalTime = ifelse(
        is.na(outcomeDate),
        as.numeric(difftime(cohortEndDate, cohortStartDate, units = "days")),
        as.numeric(difftime(outcomeDate, cohortStartDate, units = "days"))
      ),
      treatment = ifelse(treatment == 1, 1, 0)
    ) %>%
    dplyr::filter(survivalTime > 0) %>%
    dplyr::select(subjectId, treatment, outcomeCount, survivalTime) %>%
    as.data.frame()

  population$rowId <- seq_len(nrow(population))

  phResult <- testProportionalHazards(population)

  if (phResult$covariates$violated[1]) {
    result <- computeRmst(population = population, tau = NULL, tauMethod = "median")
    expect_s3_class(result, "RmstResult")
  } else {
    result <- computeRmst(population = population, tau = 365)
    expect_s3_class(result, "RmstResult")
  }

  expect_true(!is.null(result$settings))
  expect_true(!is.null(result$settings$nSubjects))
  expect_true(!is.null(result$settings$adjustmentMethod))
})
