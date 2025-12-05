/**
 * @file RmstPseudovalues.h
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

#ifndef __RmstPseudovalues_h__
#define __RmstPseudovalues_h__

#include <vector>

namespace ohdsi {
namespace cohortMethod {

class RmstPseudovalues {
public:
  static std::vector<double> computePseudovalues(
    const std::vector<double>& subjectTimes,
    const std::vector<int>& subjectEvents,
    const std::vector<double>& kmTimes,
    const std::vector<double>& kmSurv,
    const std::vector<int>& kmNRisk,
    const std::vector<int>& kmNEvent,
    double tau,
    double rmstAll
  );

private:
  static std::vector<double> computeAreaRemaining(
    const std::vector<double>& kmTimes,
    const std::vector<double>& kmSurv,
    double tau
  );

  static int findEventTimeIndex(
    double subjectTime,
    const std::vector<double>& kmTimes,
    double tolerance = 1e-10
  );
};

} // namespace cohortMethod
} // namespace ohdsi

#endif // __RmstPseudovalues_h__
