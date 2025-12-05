/**
 * @file RcppWrapper.cpp
 *
 * This file is part of CohortMethod
 *
 * Copyright 2025 Observational Health Data Sciences and Informatics
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * @author Observational Health Data Sciences and Informatics
 * @author Martijn Schuemie
 * @author Marc Suchard
 */

#include <Rcpp.h>
#include <map>
#include "Match.h"
#include "Auc.h"
#include "AdjustedKm.h"
#include "RmstPseudovalues.h"

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame matchPsInternal(std::vector<double> propensityScores, std::vector<int> treatment, unsigned int maxRatio, double caliper) {

	using namespace ohdsi::cohortMethod;

	try {
	  Match match(propensityScores, treatment, maxRatio, caliper);
		std::vector<int64_t> stratumIds = match.match();
		return DataFrame::create(_["propensityScore"] = propensityScores, _["treatment"] = treatment,_["stratumId"] = stratumIds);
	} catch (std::exception &e) {
		forward_exception_to_r(e);
	} catch (...) {
		::Rf_error("c++ exception (unknown reason)");
	}
	return DataFrame::create();
}

// [[Rcpp::export]]
std::vector<double> aucWithCi(std::vector<double> propensityScores, std::vector<int> treatment) {

  using namespace ohdsi::cohortMethod;

	try {
		std::vector<double> auc = Auc::aucWithCi(propensityScores, treatment);
		return auc;
	} catch (std::exception &e) {
		forward_exception_to_r(e);
	} catch (...) {
		::Rf_error("c++ exception (unknown reason)");
	}
  std::vector<double> auc(3,0);
	return auc;
}

// [[Rcpp::export]]
double aucWithoutCi(std::vector<double> propensityScores, std::vector<int> treatment) {

  using namespace ohdsi::cohortMethod;

  try {
		double auc = Auc::auc(propensityScores, treatment);
		return auc;
	} catch (std::exception &e) {
		forward_exception_to_r(e);
	} catch (...) {
		::Rf_error("c++ exception (unknown reason)");
	}
	return 0.0;
}

//' Compute a weight-adjusted Kaplan-Meier curve
//'
//' @param weight      Vector of observation weights
//' @param time        Vector of event times
//' @param y           Vector outcomes (0 indicates censoring, 1 indicates event-of-interest)
//'
//' @export
//'
// [[Rcpp::export]]
DataFrame adjustedKm(const std::vector<double> &weight, const std::vector<int> &time, const std::vector<int> &y) {

  using namespace ohdsi::cohortMethod;

  try {
    Surv surv = AdjustedKm::surv(weight, time, y);

    return DataFrame::create(_["time"] = surv.time, _["s"] = surv.s, _["var"] = surv.var);
  } catch (std::exception &e) {
    forward_exception_to_r(e);
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return DataFrame::create();
}

//' Compute RMST pseudovalues using Infinitesimal Jackknife
//'
//' @param subjectTimes    Subject survival times
//' @param subjectEvents   Event indicators (0=censored, 1=event)
//' @param kmTimes         KM curve event times
//' @param kmSurv          KM survival probabilities
//' @param kmNRisk         Number at risk at each event time
//' @param kmNEvent        Number of events at each event time
//' @param tau             Restriction time
//' @param rmstAll         Overall RMST
//'
//' @export
//'
// [[Rcpp::export]]
std::vector<double> computeRmstPseudovaluesInternal(
    const std::vector<double> &subjectTimes,
    const std::vector<int> &subjectEvents,
    const std::vector<double> &kmTimes,
    const std::vector<double> &kmSurv,
    const std::vector<int> &kmNRisk,
    const std::vector<int> &kmNEvent,
    double tau,
    double rmstAll) {

  using namespace ohdsi::cohortMethod;

  try {
    std::vector<double> pseudovalues = RmstPseudovalues::computePseudovalues(
      subjectTimes,
      subjectEvents,
      kmTimes,
      kmSurv,
      kmNRisk,
      kmNEvent,
      tau,
      rmstAll
    );
    return pseudovalues;
  } catch (std::exception &e) {
    forward_exception_to_r(e);
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return std::vector<double>();
}
