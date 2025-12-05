/**
 * @file RmstPseudovalues.cpp
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
 */

#include <cmath>
#include <algorithm>
#include "RmstPseudovalues.h"

namespace ohdsi {
namespace cohortMethod {

std::vector<double> RmstPseudovalues::computeAreaRemaining(
    const std::vector<double>& kmTimes,
    const std::vector<double>& kmSurv,
    double tau) {

  int nTimes = kmTimes.size();
  std::vector<double> areaRemaining(nTimes, 0.0);

  if (nTimes == 0) {
    return areaRemaining;
  }

  int lastIdx = nTimes - 1;
  if (kmTimes[lastIdx] < tau) {
    areaRemaining[lastIdx] = kmSurv[lastIdx] * (tau - kmTimes[lastIdx]);
  } else {
    areaRemaining[lastIdx] = 0.0;
  }

  for (int j = nTimes - 2; j >= 0; j--) {
    double nextTime = std::min(kmTimes[j + 1], tau);
    double segmentArea = kmSurv[j] * (nextTime - kmTimes[j]);
    areaRemaining[j] = segmentArea + areaRemaining[j + 1];
  }

  return areaRemaining;
}

int RmstPseudovalues::findEventTimeIndex(
    double subjectTime,
    const std::vector<double>& kmTimes,
    double tolerance) {

  if (kmTimes.empty()) {
    return -1;
  }

  auto it = std::lower_bound(kmTimes.begin(), kmTimes.end(), subjectTime - tolerance);

  if (it != kmTimes.end() && std::abs(*it - subjectTime) < tolerance) {
    return static_cast<int>(std::distance(kmTimes.begin(), it));
  }

  if (it != kmTimes.begin()) {
    --it;
    if (std::abs(*it - subjectTime) < tolerance) {
      return static_cast<int>(std::distance(kmTimes.begin(), it));
    }
  }

  return -1;
}

std::vector<double> RmstPseudovalues::computePseudovalues(
    const std::vector<double>& subjectTimes,
    const std::vector<int>& subjectEvents,
    const std::vector<double>& kmTimes,
    const std::vector<double>& kmSurv,
    const std::vector<int>& kmNRisk,
    const std::vector<int>& kmNEvent,
    double tau,
    double rmstAll) {

  int n = subjectTimes.size();
  int nEventTimes = kmTimes.size();
  std::vector<double> pseudovalues(n);

  if (nEventTimes == 0) {
    for (int i = 0; i < n; i++) {
      pseudovalues[i] = std::min(subjectTimes[i], tau);
    }
    return pseudovalues;
  }

  std::vector<double> areaRemaining = computeAreaRemaining(kmTimes, kmSurv, tau);

  std::vector<double> hazardIncrement(nEventTimes);
  for (int j = 0; j < nEventTimes; j++) {
    hazardIncrement[j] = static_cast<double>(kmNEvent[j]) / static_cast<double>(kmNRisk[j]);
  }

  int maxEventIdx = nEventTimes;
  for (int j = 0; j < nEventTimes; j++) {
    if (kmTimes[j] > tau) {
      maxEventIdx = j;
      break;
    }
  }

  // Infinitesimal Jackknife pseudovalues
  for (int i = 0; i < n; i++) {
    double T_i = subjectTimes[i];
    int delta_i = subjectEvents[i];
    double IF_i = 0.0;

    int eventTimeIdx = (delta_i == 1) ? findEventTimeIndex(T_i, kmTimes) : -1;

    auto upperIt = std::upper_bound(kmTimes.begin(), kmTimes.begin() + maxEventIdx, T_i);
    int jMax = static_cast<int>(std::distance(kmTimes.begin(), upperIt));

    for (int j = 0; j < jMax; j++) {
      double M_ij = (eventTimeIdx == j) ? (1.0 - hazardIncrement[j]) : (-hazardIncrement[j]);
      IF_i += areaRemaining[j] * M_ij / static_cast<double>(kmNRisk[j]);
    }

    pseudovalues[i] = rmstAll + IF_i;
  }

  return pseudovalues;
}

} // namespace cohortMethod
} // namespace ohdsi
