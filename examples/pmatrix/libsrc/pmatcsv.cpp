// --------------------------------------------
// Copyright KAPSARC. Open source MIT License.
// --------------------------------------------
// The MIT License (MIT)
//
// Copyright (c) 2015 King Abdullah Petroleum Studies and Research Center
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software
// and associated documentation files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom
// the Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or
// substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
// BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
// --------------------------------------------
// This is a very simple demonstration of reading CSV files that are in the format
// which SMPC expects.
// --------------------------------------------

#include "pmatcsv.h"
#include <easylogging++.h>

// --------------------------------------------

namespace PMatDemo {

using KBase::KMatrix;
using KBase::KException;

// --------------------------------------------

FittingParameters pccCSV(const string fs) {


  csv::ifstream inStream(fs.c_str());
  inStream.set_delimiter(',', "$$");
  inStream.enable_trim_quote_on_str(true, '\"');

  assert(inStream.is_open());
  string dummy = "";
  unsigned int numAct = 0;
  unsigned int numScen = 0;
  unsigned int numCase = 0;


  inStream.read_line();
  inStream >> dummy >> numAct;
  assert(KBase::Model::minNumActor <= numAct);

  inStream.read_line();
  inStream >> dummy >> numScen;
  assert(1 < numScen);

  inStream.read_line();
  inStream >> dummy >> numCase;
  assert(1 <= numCase);

  LOG(DEBUG) << "Actors" << numAct << ", Scenarios" << numScen << ", Cases" << numCase;

  // skip headers
  inStream.read_line();

  vector<string> aNames = {};
  auto caseWeights = KMatrix(numCase, numAct, 1.0);
  vector<bool> maxVect = {};
  auto outcomes = KMatrix(numAct, numScen);
  auto probWeight = KMatrix(numScen, numCase, 1.0);

  for (unsigned int i = 0; i < numAct; i++) {
    inStream.read_line();
    string ani = "";
    inStream >> ani;
    aNames.push_back(ani);
    for (unsigned int j = 0; j < numCase; j++) {
      double cw = 0.0;
      inStream >> cw;
      assert(0.0 <= cw);
      assert(cw <= 100.0);
      caseWeights(j, i) = cw / 100.0;
    } // loop over j, cases
    inStream >> dummy;
    if ("up" == dummy) {
      maxVect.push_back(true);
    }
    else if ("down" == dummy) {
      maxVect.push_back(false);
    }
    else {
      throw (KException("Unrecognized group-optimization direction"));
    }

    for (unsigned int j = 0; j < numScen; j++) {
      double sv = 0.0;
      inStream >> sv;
      outcomes(i, j) = sv;
    } // loop over j, scenarios
  } // loop over i, actors

  LOG(DEBUG) << "Actor names (min/max)";
  for (unsigned int i = 0; i < numAct; i++) {
    LOG(DEBUG) << (maxVect[i] ? "max" : "min") << "   " << aNames[i];
  }

  LOG(DEBUG) << "Case Weights:";
  caseWeights.mPrintf("%7.3f ");

  LOG(DEBUG) << "Outcomes:";
  outcomes.mPrintf(" %+.4e  ");


  vector<double> threshVal = {};
  vector<bool> overThresh = {};
  for (unsigned int j = 0; j < numCase; j++) {
    inStream.read_line();
    double tv = 0.0;
    string dir;
    inStream >> dummy; // skip "prob-n"
    inStream >> tv; // read a threshold value
    assert(0.0 <= tv);
    assert(tv <= 1.0);
    threshVal.push_back(tv);

    inStream >> dir; // read dir
    if ("higher" == dir) {
      overThresh.push_back(true);
    }
    else if ("lower" == dir) {
      overThresh.push_back(false);
    }
    else {
      throw (KException("Unrecognized threshold direction"));
    }

    for (unsigned int k = 0; k < numCase - 1; k++) {
      inStream >> dummy; // skip blanks, if any
    }

    for (unsigned int k = 0; k < numScen; k++) {
      double pw = 0.0;
      inStream >> pw;
      assert(0.0 <= pw);
      assert(pw <= 100.0);
      probWeight(k, j) = pw / 100.0;
    }
  } // loop over j, cases

  LOG(DEBUG) << "Prob threshholds:";
  for (unsigned int k = 0; k < numCase; k++) {
    if (overThresh[k]) {
      LOG(DEBUG) << KBase::getFormattedString("Over  %.3f", threshVal[k]);
    }
    else {
      LOG(DEBUG) << KBase::getFormattedString("Under %.3f", threshVal[k]);
    }
  }

  LOG(DEBUG) << "ProbWeights:";
  probWeight.mPrintf(" %5.3f ");

  auto rslt = FittingParameters(aNames, maxVect,
                                outcomes, caseWeights,
                                probWeight, threshVal, overThresh);

  return rslt;
}



}
// end of namespace

// --------------------------------------------
// Copyright KAPSARC. Open source MIT License.
// --------------------------------------------
