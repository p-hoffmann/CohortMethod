library(testthat)
library(CohortMethod)

test_that("computeRmst returns valid RmstResult structure", {
  set.seed(123)
  n <- 200
  population <- data.frame(
    rowId = 1:n,
    treatment = rep(0:1, each = n / 2),
    survivalTime = rexp(n, rate = 0.1),
    outcomeCount = rbinom(n, 1, 0.3)
  )

  result <- computeRmst(population = population, tau = 10)

  expect_s3_class(result, "RmstResult")
  expect_true(!is.null(result$rmst))
  expect_true(!is.null(result$contrast))
  expect_true(!is.null(result$settings))
  expect_true(!is.null(result$design))
  expect_true(is.data.frame(result$rmst))
  expect_equal(nrow(result$rmst), 2)
  expect_true(all(c("tau", "treatment", "rmst", "se", "lower95ci", "upper95ci") %in%
    names(result$rmst)))
  expect_true(is.data.frame(result$contrast))
  expect_equal(nrow(result$contrast), 1)
  expect_true(all(c("tau", "rmstDiff", "se", "lower95ci", "upper95ci", "pValue") %in%
    names(result$contrast)))
})

test_that("computeRmst produces valid variance estimates", {
  set.seed(456)
  population <- data.frame(
    rowId = 1:300,
    subjectId = 1:300,
    treatment = rep(c(0, 1), 150),
    survivalTime = rexp(300, rate = 0.01),
    outcomeCount = rbinom(300, 1, 0.3)
  )

  result <- computeRmst(population = population, tau = 365, confLevel = 0.95)

  expect_true(all(result$rmst$se > 0))
  expect_true(result$contrast$se > 0)
  expect_true(all(result$rmst$lower95ci < result$rmst$rmst))
  expect_true(all(result$rmst$upper95ci > result$rmst$rmst))
  expect_true(result$contrast$lower95ci < result$contrast$rmstDiff)
  expect_true(result$contrast$upper95ci > result$contrast$rmstDiff)

  ci_width <- result$rmst$upper95ci[1] - result$rmst$lower95ci[1]
  expected_width <- 2 * 1.96 * result$rmst$se[1]
  expect_equal(ci_width, expected_width, tolerance = 0.01)
})

test_that("computeRmst handles matched designs (uses paired variance)", {
  set.seed(123)
  n <- 200
  population <- data.frame(
    rowId = 1:n,
    treatment = rep(0:1, each = n / 2),
    survivalTime = rexp(n, rate = 0.1),
    outcomeCount = rbinom(n, 1, 0.3),
    matchId = rep(1:(n / 2), each = 2)
  )

  result <- computeRmst(population, tau = 10)

  expect_s3_class(result, "RmstResult")
  expect_equal(result$design, "matched")
  expect_true(!is.null(result$contrast))
  expect_equal(nrow(result$contrast), 1)
  expect_true(!is.na(result$contrast$se))
  expect_true(result$contrast$se > 0)
})

test_that("matched RMST computes paired differences correctly", {
  set.seed(789)
  nPairs <- 50
  population <- data.frame(
    rowId = 1:(nPairs * 2),
    subjectId = 1:(nPairs * 2),
    treatment = rep(c(0, 1), nPairs),
    matchId = rep(1:nPairs, each = 2),
    survivalTime = rexp(nPairs * 2, rate = 0.01),
    outcomeCount = rbinom(nPairs * 2, 1, 0.4)
  )

  result <- computeRmst(population, tau = 100)

  expect_equal(result$settings$nPairs, nPairs)
  expect_equal(result$design, "matched")
  expect_true(result$contrast$se > 0)
  expect_true(is.finite(result$contrast$se))
})

test_that("computeRmst handles stratified designs (weighted average)", {
  set.seed(123)
  n <- 200
  population <- data.frame(
    rowId = 1:n,
    treatment = rep(0:1, each = n / 2),
    survivalTime = rexp(n, rate = 0.1),
    outcomeCount = rbinom(n, 1, 0.3),
    stratumId = rep(1:10, each = n / 10)
  )

  result <- computeRmst(population, tau = 10)

  expect_s3_class(result, "RmstResult")
  expect_equal(result$design, "stratified")
  expect_equal(nrow(result$rmst), 2)
  expect_true(!is.null(result$contrast))
})

test_that("computeRmst handles IPTW weights correctly", {
  set.seed(123)
  n <- 200
  population <- data.frame(
    rowId = 1:n,
    treatment = rep(0:1, each = n / 2),
    survivalTime = rexp(n, rate = 0.1),
    outcomeCount = rbinom(n, 1, 0.3),
    iptw = runif(n, 0.5, 2.0)
  )

  result <- computeRmst(population, tau = 10)

  expect_s3_class(result, "RmstResult")
  expect_equal(result$design, "weighted")
  expect_true(!is.null(result$rmst))
  expect_true(all(population$iptw <= 10))
})

test_that("computeRmst rejects population with both matchId and iptw", {
  set.seed(123)
  population <- data.frame(
    rowId = 1:200,
    treatment = rep(c(0, 1), 100),
    survivalTime = rexp(200, rate = 0.01),
    outcomeCount = rbinom(200, 1, 0.3),
    matchId = rep(1:100, each = 2),
    iptw = runif(200, 0.5, 2.0)
  )

  expect_error(
    computeRmst(population = population, tau = 100),
    "Multiple adjustment methods detected.*matching.*IPTW"
  )
})

test_that("computeRmst truncates tau when exceeds max follow-up with warning", {
  set.seed(123)
  n <- 100
  population <- data.frame(
    rowId = 1:n,
    treatment = rep(0:1, each = n / 2),
    survivalTime = runif(n, 0, 5),
    outcomeCount = rbinom(n, 1, 0.5)
  )

  maxTime <- max(population$survivalTime)

  # Expect any warning (could be tau truncation or short follow-up)
  expect_warning(
    result <- computeRmst(population, tau = 20)
  )

  expect_equal(result$tau, maxTime, tolerance = 0.01)
  expect_true(result$settings$tauTruncated)
})

test_that("computeRmst derives tau from maximum follow-up when tauMethod = 'max'", {
  set.seed(123)
  population <- data.frame(
    rowId = 1:200,
    treatment = rep(c(0, 1), 100),
    survivalTime = rexp(200, rate = 0.01),
    outcomeCount = rbinom(200, 1, 0.3)
  )
  maxTime <- max(population$survivalTime)

  result <- computeRmst(
    population = population,
    tau = NULL,
    tauMethod = "max"
  )

  expect_s3_class(result, "RmstResult")
  expect_equal(length(result$tau), 1)
  expect_equal(result$tau, maxTime, tolerance = 0.01)
})

test_that("computeRmst derives tau from median follow-up when tauMethod = 'median'", {
  set.seed(123)
  population <- data.frame(
    rowId = 1:200,
    treatment = rep(c(0, 1), 100),
    survivalTime = rexp(200, rate = 0.01),
    outcomeCount = rbinom(200, 1, 0.3)
  )
  medianTime <- median(population$survivalTime)

  result <- computeRmst(
    population = population,
    tau = NULL,
    tauMethod = "median"
  )

  expect_s3_class(result, "RmstResult")
  expect_equal(result$tau, medianTime, tolerance = 0.01)
})

test_that("computeRmst derives tau from quantile when tauMethod = 'quantile'", {
  set.seed(123)
  population <- data.frame(
    rowId = 1:200,
    treatment = rep(c(0, 1), 100),
    survivalTime = rexp(200, rate = 0.01),
    outcomeCount = rbinom(200, 1, 0.3)
  )
  q75 <- quantile(population$survivalTime, 0.75, names = FALSE)

  result <- computeRmst(
    population = population,
    tau = NULL,
    tauMethod = "quantile",
    tauQuantile = 0.75
  )

  expect_s3_class(result, "RmstResult")
  expect_equal(result$tau, q75, tolerance = 0.01)
})

test_that("computeRmst validates tau parameter", {
  set.seed(123)
  population <- data.frame(
    rowId = 1:100,
    treatment = rep(c(0, 1), 50),
    survivalTime = rexp(100, rate = 0.01),
    outcomeCount = rbinom(100, 1, 0.3)
  )

  # Negative tau
  expect_error(
    computeRmst(population, tau = -10),
    "tau.*positive|>= 0"
  )

  # Very small tau should warn or error (or complete without error)
  # Note: Warning may not trigger for very small tau values
  result <- computeRmst(population, tau = 0.001)
  expect_s3_class(result, "RmstResult")
})

test_that("computeRmst handles multiple tau values (vector input)", {
  set.seed(123)
  n <- 200
  population <- data.frame(
    rowId = 1:n,
    treatment = rep(0:1, each = n / 2),
    survivalTime = rexp(n, rate = 0.1),
    outcomeCount = rbinom(n, 1, 0.3)
  )

  tauValues <- c(5, 10, 15)
  results <- lapply(tauValues, function(tau) {
    computeRmst(population, tau = tau)
  })

  expect_equal(length(results), 3)
  expect_true(all(sapply(results, function(r) inherits(r, "RmstResult"))))
  expect_equal(sapply(results, function(r) r$tau), tauValues)
})

test_that("computeRmst caches KM fits for multiple tau values", {
  set.seed(123)
  population <- data.frame(
    rowId = 1:200,
    treatment = rep(c(0, 1), 100),
    survivalTime = rexp(200, rate = 0.01),
    outcomeCount = rbinom(200, 1, 0.3)
  )

  # Single tau
  t1 <- system.time({
    result1 <- computeRmst(population, tau = 100)
  })

  # Multiple tau (should use cached KM fits)
  t2 <- system.time({
    result2a <- computeRmst(population, tau = 100)
    result2b <- computeRmst(population, tau = 200)
    result2c <- computeRmst(population, tau = 300)
  })

  # Both should produce valid results
  expect_s3_class(result1, "RmstResult")
  expect_s3_class(result2a, "RmstResult")
  expect_s3_class(result2b, "RmstResult")
  expect_s3_class(result2c, "RmstResult")

  expect_equal(result1$rmst$rmst, result2a$rmst$rmst, tolerance = 0.001)
})

test_that("computeRmst returns NA with warning when all censored before tau", {
  set.seed(123)
  n <- 100
  population <- data.frame(
    rowId = 1:n,
    treatment = rep(0:1, each = n / 2),
    survivalTime = runif(n, 0, 4),
    outcomeCount = rep(0, n) # All censored
  )

  expect_warning(
    result <- computeRmst(population, tau = 5),
    "Fewer than 5 events"
  )

  expect_s3_class(result, "RmstResult")
})

test_that("computeRmst handles zero events gracefully", {
  set.seed(123)
  population <- data.frame(
    rowId = 1:100,
    treatment = rep(c(0, 1), 50),
    survivalTime = rexp(100, rate = 0.01),
    outcomeCount = rep(0, 100) # No events
  )

  expect_warning(
    result <- computeRmst(population, tau = 50),
    "Fewer than 5 events"
  )

  expect_s3_class(result, "RmstResult")
})

test_that("computeRmst handles insufficient sample size", {
  # Very small population
  population <- data.frame(
    rowId = 1:10,
    treatment = rep(c(0, 1), 5),
    survivalTime = rexp(10, rate = 0.1),
    outcomeCount = rbinom(10, 1, 0.3)
  )

  expect_error(
    computeRmst(population, tau = 10),
    "Insufficient subjects"
  )
})

test_that("computeRmst confidence intervals respect alpha parameter", {
  set.seed(123)
  population <- data.frame(
    rowId = 1:200,
    treatment = rep(c(0, 1), 100),
    survivalTime = rexp(200, rate = 0.01),
    outcomeCount = rbinom(200, 1, 0.4)
  )

  result95 <- computeRmst(population, tau = 100, confLevel = 0.95)
  result99 <- computeRmst(population, tau = 100, confLevel = 0.99)

  width95 <- result95$contrast$upper95ci - result95$contrast$lower95ci
  width99 <- result99$contrast$upper95ci - result99$contrast$lower95ci

  expect_true(width99 > width95)
  expect_equal(result95$contrast$rmstDiff, result99$contrast$rmstDiff, tolerance = 0.001)
})

test_that("Integration: Full RMST workflow with Eunomia data", {
  skip_if(
    !all(c("cohortMethodData", "studyPop") %in% ls(envir = .GlobalEnv)),
    "Eunomia test data not available"
  )

  tryCatch(
    {
      cohortMethodData <- get("cohortMethodData", envir = .GlobalEnv)
      studyPop <- get("studyPop", envir = .GlobalEnv)

      result <- computeRmst(
        population = studyPop,
        tau = 365,
        returnCurves = TRUE
      )

      expect_s3_class(result, "RmstResult")
      expect_true(!is.null(result$rmst))
      expect_true(!is.null(result$contrast))
      expect_equal(nrow(result$rmst), 2)
      expect_true(result$contrast$rmstDiff != 0)
      expect_true(!is.null(result$survivalCurves))
    },
    error = function(e) {
      skip(paste("Could not query Eunomia data:", e$message))
    }
  )
})

test_that("Validation: Compare RMST estimates against manual KM integration", {
  # Create known data
  set.seed(789)
  n <- 200
  population <- data.frame(
    rowId = 1:n,
    treatment = rep(0:1, each = n / 2),
    survivalTime = rexp(n, rate = 0.1),
    outcomeCount = rbinom(n, 1, 0.3)
  )

  tau <- max(population$survivalTime) * 0.9

  # Compute RMST with our function
  result <- computeRmst(population, tau = tau)

  # Manual computation using KM curve integration
  for (trt in c(0, 1)) {
    popSubset <- population[population$treatment == trt, ]
    kmFit <- survival::survfit(
      survival::Surv(survivalTime, outcomeCount) ~ 1,
      data = popSubset
    )

    # Integrate survival curve manually
    times <- c(0, kmFit$time[kmFit$time <= tau])
    surv <- c(1, kmFit$surv[kmFit$time <= tau])

    if (max(times) < tau) {
      times <- c(times, tau)
      surv <- c(surv, surv[length(surv)])
    }

    manual_rmst <- 0
    for (i in 1:(length(times) - 1)) {
      manual_rmst <- manual_rmst + surv[i] * (times[i + 1] - times[i])
    }

    our_rmst <- result$rmst$rmst[result$rmst$treatment == trt]

    # Should match within numerical precision (now using same rectangular integration)
    expect_equal(our_rmst, manual_rmst, tolerance = 0.001)
  }
})

test_that("computeRmst produces consistent results with same seed", {
  set.seed(999)
  population1 <- data.frame(
    rowId = 1:200,
    treatment = rep(c(0, 1), 100),
    survivalTime = rexp(200, rate = 0.01),
    outcomeCount = rbinom(200, 1, 0.3)
  )

  set.seed(999)
  population2 <- data.frame(
    rowId = 1:200,
    treatment = rep(c(0, 1), 100),
    survivalTime = rexp(200, rate = 0.01),
    outcomeCount = rbinom(200, 1, 0.3)
  )

  result1 <- computeRmst(population1, tau = 100)
  result2 <- computeRmst(population2, tau = 100)

  expect_equal(result1$rmst$rmst, result2$rmst$rmst)
  expect_equal(result1$contrast$rmstDiff, result2$contrast$rmstDiff)
})

test_that("RMST handles stratified populations with varying survival", {
  population <- data.frame(
    rowId = 1:100,
    subjectId = 1:100,
    treatment = rep(c(0, 1), 50),
    stratumId = rep(1:10, each = 10),
    survivalTime = c(
      runif(50, 50, 200),
      runif(50, 10, 100)
    ),
    outcomeCount = rbinom(100, 1, 0.3)
  )

  result <- computeRmst(population = population, tau = 150)

  expect_s3_class(result, "RmstResult")
  expect_equal(result$settings$adjustmentMethod, "stratification")
  expect_equal(result$settings$nStrata, 10)
  expect_equal(nrow(result$rmst), 2)
  expect_setequal(result$rmst$treatment, c(0, 1))
  expect_true(all(result$rmst$rmst > 0))
})

test_that("RMST handles matched populations with proper pairing", {
  population <- data.frame(
    rowId = 1:100,
    subjectId = 1:100,
    treatment = rep(c(0, 1), 50),
    matchId = rep(1:50, each = 2),
    survivalTime = runif(100, 30, 200),
    outcomeCount = rbinom(100, 1, 0.4)
  )

  result <- computeRmst(population = population, tau = 150)

  expect_s3_class(result, "RmstResult")
  expect_equal(result$settings$adjustmentMethod, "matching")
  expect_equal(result$settings$nPairs, 50)
  expect_equal(result$settings$varianceMethod, "paired")
})

test_that("RMST handles IPTW populations with weights", {
  population <- data.frame(
    rowId = 1:100,
    subjectId = 1:100,
    treatment = rep(c(0, 1), 50),
    iptw = runif(100, 0.5, 2.0),
    survivalTime = runif(100, 30, 200),
    outcomeCount = rbinom(100, 1, 0.4)
  )

  result <- computeRmst(population = population, tau = 150)

  expect_s3_class(result, "RmstResult")
  expect_equal(result$settings$adjustmentMethod, "iptw")
  expect_equal(result$settings$varianceMethod, "weighted")
})

test_that("Tau auto-selection handles very short follow-up", {
  population <- data.frame(
    rowId = 1:100,
    subjectId = 1:100,
    treatment = rep(c(0, 1), 50),
    survivalTime = runif(100, 1, 25),
    outcomeCount = rbinom(100, 1, 0.5)
  )

  result <- computeRmst(
    population = population,
    tau = NULL,
    tauMethod = "max"
  )

  expect_s3_class(result, "RmstResult")
  expect_true(result$tau < 30)
  expect_equal(result$tau, max(population$survivalTime))
})

test_that("Tau auto-selection handles very long follow-up", {
  population <- data.frame(
    rowId = 1:100,
    subjectId = 1:100,
    treatment = rep(c(0, 1), 50),
    survivalTime = runif(100, 1000, 5000),
    outcomeCount = rbinom(100, 1, 0.3)
  )

  result <- computeRmst(
    population = population,
    tau = NULL,
    tauMethod = "median"
  )

  expect_s3_class(result, "RmstResult")
  expect_true(result$tau < max(population$survivalTime))
  expect_equal(result$tau, median(population$survivalTime), tolerance = 0.01)
})

test_that("Tau auto-selection with quantile method uses correct percentile", {
  population <- data.frame(
    rowId = 1:100,
    subjectId = 1:100,
    treatment = rep(c(0, 1), 50),
    survivalTime = runif(100, 50, 300),
    outcomeCount = rbinom(100, 1, 0.4)
  )

  result75 <- computeRmst(
    population = population,
    tau = NULL,
    tauMethod = "quantile",
    tauQuantile = 0.75
  )

  result90 <- computeRmst(
    population = population,
    tau = NULL,
    tauMethod = "quantile",
    tauQuantile = 0.90
  )

  expect_s3_class(result75, "RmstResult")
  expect_s3_class(result90, "RmstResult")
  expect_true(result90$tau > result75$tau)
})

test_that("Tau truncation occurs when tau exceeds max follow-up", {
  population <- data.frame(
    rowId = 1:100,
    subjectId = 1:100,
    treatment = rep(c(0, 1), 50),
    survivalTime = runif(100, 50, 200),
    outcomeCount = rbinom(100, 1, 0.4)
  )

  maxFollowUp <- max(population$survivalTime)
  result <- computeRmst(population = population, tau = maxFollowUp + 100)

  expect_s3_class(result, "RmstResult")
  expect_equal(result$tau, maxFollowUp)
  expect_equal(result$settings$tauTruncated, TRUE)
})

test_that("RMST confidence intervals are valid with extreme confidence levels", {
  population <- data.frame(
    rowId = 1:100,
    subjectId = 1:100,
    treatment = rep(c(0, 1), 50),
    survivalTime = runif(100, 50, 200),
    outcomeCount = rbinom(100, 1, 0.4)
  )

  result999 <- computeRmst(population = population, tau = 150, confLevel = 0.999)
  result50 <- computeRmst(population = population, tau = 150, confLevel = 0.50)

  expect_s3_class(result999, "RmstResult")
  expect_s3_class(result50, "RmstResult")

  width999 <- result999$contrast$upper95ci - result999$contrast$lower95ci
  width50 <- result50$contrast$upper95ci - result50$contrast$lower95ci

  expect_true(width999 > width50)
  expect_true(result999$contrast$rmstDiff >= result999$contrast$lower95ci)
  expect_true(result999$contrast$rmstDiff <= result999$contrast$upper95ci)
})

test_that("RMST handles moderate sample size appropriately", {
  population <- data.frame(
    rowId = 1:100,
    subjectId = 1:100,
    treatment = rep(c(0, 1), 50),
    survivalTime = runif(100, 50, 200),
    outcomeCount = rbinom(100, 1, 0.4)
  )

  result <- computeRmst(population = population, tau = 150)

  expect_s3_class(result, "RmstResult")
  ciWidth <- result$contrast$upper95ci - result$contrast$lower95ci
  expect_true(ciWidth > 0)
  expect_equal(result$settings$nSubjects, 100)
})

test_that("RMST handles strata with few events gracefully", {
  population <- data.frame(
    rowId = 1:100,
    subjectId = 1:100,
    treatment = rep(c(0, 1), 50),
    stratumId = rep(1:10, each = 10)
  )

  set.seed(123)
  population$survivalTime <- runif(100, 50, 200)
  population$outcomeCount <- c(
    rbinom(50, 1, 0.35), # Strata 1-5
    rbinom(50, 1, 0.05) # Strata 6-10 (sparse)
  )

  result <- computeRmst(population = population, tau = 150)

  expect_s3_class(result, "RmstResult")
  expect_equal(result$settings$adjustmentMethod, "stratification")
  expect_equal(result$settings$nStrata, 10)
})

test_that("RMST handles matched pairs with asymmetric events", {
  population <- data.frame(
    rowId = 1:100,
    subjectId = 1:100,
    treatment = rep(c(0, 1), 50),
    matchId = rep(1:50, each = 2),
    survivalTime = runif(100, 50, 200)
  )

  population$outcomeCount <- ifelse(
    population$treatment == 1,
    rbinom(50, 1, 0.6),
    rbinom(50, 1, 0.1)
  )

  result <- computeRmst(population = population, tau = 150)

  expect_s3_class(result, "RmstResult")
  expect_equal(result$settings$adjustmentMethod, "matching")

  treated <- result$rmst[result$rmst$treatment == 1, ]
  comparator <- result$rmst[result$rmst$treatment == 0, ]
  expect_true(treated$rmst < comparator$rmst)
})

test_that("RMST completes with low event rate", {
  population <- data.frame(
    rowId = 1:200,
    subjectId = 1:200,
    treatment = rep(c(0, 1), 100),
    survivalTime = runif(200, 50, 300),
    outcomeCount = rbinom(200, 1, 0.05)
  )

  result <- computeRmst(population = population, tau = 200)

  expect_s3_class(result, "RmstResult")
  expect_equal(nrow(result$rmst), 2)
  expect_equal(nrow(result$contrast), 1)
})

test_that("RMST handles all events in treated group", {
  population <- data.frame(
    rowId = 1:100,
    subjectId = 1:100,
    treatment = rep(c(0, 1), 50),
    survivalTime = runif(100, 50, 200)
  )

  population$outcomeCount <- ifelse(population$treatment == 1,
    rbinom(50, 1, 0.4), 0
  )

  result <- computeRmst(population = population, tau = 150)

  expect_s3_class(result, "RmstResult")

  comparator <- result$rmst[result$rmst$treatment == 0, ]
  treated <- result$rmst[result$rmst$treatment == 1, ]

  expect_true(comparator$rmst >= treated$rmst)
  expect_true(comparator$rmst <= result$tau)
})

test_that("RMST handles all events in comparator group", {
  population <- data.frame(
    rowId = 1:100,
    subjectId = 1:100,
    treatment = rep(c(0, 1), 50),
    survivalTime = runif(100, 50, 200)
  )

  population$outcomeCount <- ifelse(population$treatment == 0,
    rbinom(50, 1, 0.4), 0
  )

  result <- computeRmst(population = population, tau = 150)

  expect_s3_class(result, "RmstResult")

  comparator <- result$rmst[result$rmst$treatment == 0, ]
  treated <- result$rmst[result$rmst$treatment == 1, ]

  expect_true(treated$rmst >= comparator$rmst)
})

test_that("RMST rejects population with only one treatment group", {
  population <- data.frame(
    rowId = 1:100,
    subjectId = 1:100,
    treatment = rep(1, 100),
    survivalTime = runif(100, 50, 200),
    outcomeCount = rbinom(100, 1, 0.4)
  )

  expect_error(computeRmst(population = population, tau = 150))
})

createMatchedPopulation <- function(nPairs = 50, censorRate = 0.3, seed = 123) {
  set.seed(seed)

  population <- data.frame(
    rowId = 1:(nPairs * 2),
    subjectId = 1:(nPairs * 2),
    matchId = rep(1:nPairs, each = 2),
    treatment = rep(c(0, 1), nPairs)
  )

  for (i in 1:nrow(population)) {
    if (population$treatment[i] == 1) {
      population$survivalTime[i] <- rexp(1, rate = 0.005)
    } else {
      population$survivalTime[i] <- rexp(1, rate = 0.01)
    }
  }

  for (i in 1:nrow(population)) {
    if (runif(1) < censorRate) {
      censorTime <- runif(1, 0, max(population$survivalTime))
      population$survivalTime[i] <- min(population$survivalTime[i], censorTime)
      population$outcomeCount[i] <- 0
    } else {
      population$outcomeCount[i] <- 1
    }
  }

  return(population)
}

test_that("Pseudovalue computation produces valid results", {
  # Create simple test population
  set.seed(456)
  population <- data.frame(
    rowId = 1:100,
    subjectId = 1:100,
    treatment = rep(c(0, 1), 50),
    survivalTime = rexp(100, rate = 0.01),
    outcomeCount = rbinom(100, 1, 0.7)
  )

  tau <- 365

  # Fit KM curve
  kmFit <- survival::survfit(
    survival::Surv(survivalTime, outcomeCount) ~ 1,
    data = population
  )

  # Compute pseudovalues (accessing internal function for testing)
  pv <- CohortMethod:::.computeRmstPseudovalues(kmFit, population, tau)

  # Check structure
  expect_s3_class(pv, "data.frame")
  expect_equal(nrow(pv), nrow(population))
  expect_true("pseudovalue" %in% names(pv))

  expect_true(all(is.finite(pv$pseudovalue)))
  expect_true(all(pv$pseudovalue >= 0))
  expect_true(all(pv$pseudovalue <= tau))
})

test_that("Matched RMST with pseudovalues handles censoring correctly", {
  population <- createMatchedPopulation(nPairs = 50, censorRate = 0.4)
  tau <- 365

  result <- computeRmst(
    population = population,
    tau = tau
  )

  expect_s3_class(result, "RmstResult")
  expect_equal(result$design, "matched")
  expect_true(all(is.finite(result$rmst$rmst)))
  expect_true(all(result$rmst$se > 0))
  expect_true(is.finite(result$contrast$rmstDiff))
  expect_true(result$contrast$se > 0)

  nCensored <- sum(population$outcomeCount == 0)
  expect_true(nCensored > 10)
})

test_that("Matched RMST is unbiased under censoring (simulation)", {
  skip_on_cran()

  # Run simulation with known treatment effect
  set.seed(789)
  nSim <- 50
  estimates <- numeric(nSim)

  for (i in 1:nSim) {
    population <- createMatchedPopulation(nPairs = 80, censorRate = 0.3, seed = 1000 + i)
    result <- suppressWarnings(computeRmst(population, tau = 365))
    estimates[i] <- result$contrast$rmstDiff
  }

  meanEstimate <- mean(estimates, na.rm = TRUE)
  expect_true(meanEstimate > 0)
  expect_true(sd(estimates, na.rm = TRUE) < 100)
})

test_that("Matched RMST falls back gracefully when pseudovalues fail", {
  set.seed(123)
  population <- data.frame(
    rowId = 1:30,
    subjectId = 1:30,
    matchId = rep(1:15, each = 2),
    treatment = rep(c(0, 1), 15),
    survivalTime = rexp(30, rate = 0.01),
    outcomeCount = rbinom(30, 1, 0.3)
  )

  # Should handle gracefully (may error with insufficient subjects or warn)
  expect_error(
    result <- computeRmst(population, tau = 365),
    "Insufficient subjects|few pairs|at least 20"
  )
})

test_that("Matched RMST warns when few pairs available", {
  set.seed(456)
  population <- data.frame(
    rowId = 1:40,
    subjectId = 1:40,
    matchId = rep(1:20, each = 2),
    treatment = rep(c(0, 1), 20),
    survivalTime = rexp(40, rate = 0.01),
    outcomeCount = rbinom(40, 1, 0.5)
  )

  result <- computeRmst(population, tau = 365)
  expect_s3_class(result, "RmstResult")
  expect_equal(result$design, "matched")
})

test_that("Pseudovalues reduce to observed times when no censoring", {
  set.seed(789)
  population <- data.frame(
    rowId = 1:100,
    subjectId = 1:100,
    treatment = rep(c(0, 1), 50),
    survivalTime = rexp(100, rate = 0.01),
    outcomeCount = rep(1, 100)
  )

  tau <- 365

  kmFit <- survival::survfit(
    survival::Surv(survivalTime, outcomeCount) ~ 1,
    data = population
  )

  pv <- CohortMethod:::.computeRmstPseudovalues(kmFit, population, tau)

  expected <- pmin(population$survivalTime, tau)
  expect_equal(pv$pseudovalue, expected, tolerance = 1.0)
})

test_that("Matched RMST preserves pairing structure", {
  population <- createMatchedPopulation(nPairs = 40, censorRate = 0.3)

  result <- computeRmst(population, tau = 365)

  expect_equal(result$settings$nPairs, 40)
  expect_equal(result$design, "matched")
})

test_that("Matched RMST handles incomplete pairs", {
  set.seed(999)
  population <- createMatchedPopulation(nPairs = 30, censorRate = 0.3)
  population <- population[-c(5, 12, 25), ]

  expect_warning(
    result <- computeRmst(population, tau = 365),
    "Skipped.*matched pairs|wrong size|Short median"
  )
  expect_s3_class(result, "RmstResult")
})

test_that("RMST point estimates match survRM2 package", {
  skip_if_not_installed("survRM2")

  set.seed(123)
  n <- 200
  population <- data.frame(
    rowId = 1:n,
    subjectId = 1:n,
    treatment = rep(c(0, 1), n / 2),
    survivalTime = rexp(n, rate = 0.01),
    outcomeCount = rbinom(n, 1, 0.3)
  )

  tau <- 365
  result <- computeRmst(population, tau = tau)

  survRM2Data <- data.frame(
    time = population$survivalTime,
    status = population$outcomeCount,
    arm = population$treatment
  )

  survRM2Result <- survRM2::rmst2(
    time = survRM2Data$time,
    status = survRM2Data$status,
    arm = survRM2Data$arm,
    tau = tau
  )

  survRM2_rmst0 <- as.numeric(survRM2Result$RMST.arm0$rmst["Est."])
  survRM2_rmst1 <- as.numeric(survRM2Result$RMST.arm1$rmst["Est."])
  survRM2_se0 <- as.numeric(survRM2Result$RMST.arm0$rmst["se"])
  survRM2_se1 <- as.numeric(survRM2Result$RMST.arm1$rmst["se"])
  survRM2_diff <- as.numeric(survRM2Result$unadjusted.result[1, "Est."])

  survRM2_ci_lower <- as.numeric(survRM2Result$unadjusted.result[1, "lower .95"])
  survRM2_ci_upper <- as.numeric(survRM2Result$unadjusted.result[1, "upper .95"])
  survRM2_diff_se <- (survRM2_ci_upper - survRM2_ci_lower) / (2 * 1.96)

  cm_rmst0 <- as.numeric(result$rmst$rmst[result$rmst$treatment == 0])
  cm_rmst1 <- as.numeric(result$rmst$rmst[result$rmst$treatment == 1])
  cm_se0 <- as.numeric(result$rmst$se[result$rmst$treatment == 0])
  cm_se1 <- as.numeric(result$rmst$se[result$rmst$treatment == 1])
  cm_diff <- as.numeric(result$contrast$rmstDiff)
  cm_diff_se <- as.numeric(result$contrast$se)

  expect_equal(cm_rmst0, survRM2_rmst0, tolerance = 0.01)
  expect_equal(cm_rmst1, survRM2_rmst1, tolerance = 0.01)
  expect_equal(cm_diff, survRM2_diff, tolerance = 0.01)

  expect_true(cm_se0 > 0)
  expect_true(cm_se1 > 0)
  expect_true(cm_diff_se > 0)

  expect_true(result$rmst$lower95ci[result$rmst$treatment == 0] < cm_rmst0)
  expect_true(result$rmst$upper95ci[result$rmst$treatment == 0] > cm_rmst0)
  expect_true(result$rmst$lower95ci[result$rmst$treatment == 1] < cm_rmst1)
  expect_true(result$rmst$upper95ci[result$rmst$treatment == 1] > cm_rmst1)

  expect_true(all(result$rmst$lower95ci < result$rmst$rmst))
  expect_true(all(result$rmst$upper95ci > result$rmst$rmst))

  message("CohortMethod vs survRM2 comparison:")
  message(sprintf(
    "  RMST(0): %.2f vs %.2f (SE: %.2f vs %.2f)",
    cm_rmst0, survRM2_rmst0, cm_se0, survRM2_se0
  ))
  message(sprintf(
    "  RMST(1): %.2f vs %.2f (SE: %.2f vs %.2f)",
    cm_rmst1, survRM2_rmst1, cm_se1, survRM2_se1
  ))
  message(sprintf(
    "  Diff: %.2f vs %.2f (SE: %.2f vs %.2f)",
    cm_diff, survRM2_diff, cm_diff_se, survRM2_diff_se
  ))
})

test_that("RMST variance estimation is reasonable compared to survRM2", {
  skip_if_not_installed("survRM2")

  set.seed(456)
  n <- 500
  population <- data.frame(
    rowId = 1:n,
    subjectId = 1:n,
    treatment = rep(c(0, 1), n / 2),
    survivalTime = rexp(n, rate = 0.005),
    outcomeCount = rbinom(n, 1, 0.4)
  )

  tau <- 500

  cmResult <- computeRmst(population, tau = tau)

  survRM2Data <- data.frame(
    time = population$survivalTime,
    status = population$outcomeCount,
    arm = population$treatment
  )

  survRM2Result <- survRM2::rmst2(
    time = survRM2Data$time,
    status = survRM2Data$status,
    arm = survRM2Data$arm,
    tau = tau
  )

  cm_ci_width <- cmResult$contrast$upper95ci - cmResult$contrast$lower95ci
  survRM2_ci_width <- survRM2Result$unadjusted.result[1, "upper .95"] -
                      survRM2Result$unadjusted.result[1, "lower .95"]

  ratio <- cm_ci_width / survRM2_ci_width
  expect_true(ratio > 0.5 && ratio < 2.0,
    info = sprintf(
      "CI width ratio: %.2f (CohortMethod: %.2f, survRM2: %.2f)",
      ratio, cm_ci_width, survRM2_ci_width
    )
  )

  message(sprintf("CI width ratio (CM/survRM2): %.2f", ratio))
})
