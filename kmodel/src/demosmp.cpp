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
//
// Demonstrate a very basic Spatial Model of Politics.
// 
// --------------------------------------------


#include "kutils.h"
#include "kmodel.h"
#include "demo.h"
#include "demosmp.h"
#include "csv_parser.hpp"


using KBase::PRNG;
using KBase::KMatrix;
using KBase::Actor;
using KBase::Model;
using KBase::Position;
using KBase::State;
using KBase::VotingRule;


namespace DemoSMP {
  using std::cout;
  using std::endl;
  using std::flush;
  using std::function;
  using std::get;
  using std::string;

  using KBase::ReportingLevel;

  BargainSMP::BargainSMP(const SMPActor* ai, const SMPActor* ar, const VctrPstn & pi, const VctrPstn & pr) {
    assert(nullptr != ai);
    assert(nullptr != ar);
    actInit = ai;
    actRcvr = ar;
    posInit = pi;
    posRcvr = pr;
  }

  BargainSMP::~BargainSMP() {
    actInit = nullptr;
    actRcvr = nullptr;
    posInit = KMatrix(0, 0);
    posRcvr = KMatrix(0, 0);
  }

  SMPActor::SMPActor(string n, string d) : Actor(n, d){
    vr = VotingRule::Proportional; // just a default 
  }

  SMPActor::~SMPActor() {
  }

  double SMPActor::vote(unsigned int, unsigned int, const State*) const {
    assert(false);  // TDO: finish this
    return 0;
  }


  double SMPActor::vote(const Position * ap1, const Position * ap2, const SMPState* ast) const {
    double u1 = posUtil(ap1, ast);
    double u2 = posUtil(ap2, ast);
    double v = Model::vote(vr, sCap, u1, u2);
    return v;
  }


  double SMPActor::posUtil(const Position * ap1, const SMPState* as) const {
    assert(nullptr != as);
    int ai = as->model->actrNdx(this);
    double ri = as->nra(ai, 0);
    assert(0 <= ai);
    const VctrPstn* p0 = ((const VctrPstn*)(as->pstns[ai]));
    assert(nullptr != p0);
    auto p1 = ((const VctrPstn*)ap1);
    assert(nullptr != p1);
    double u1 = SMPModel::bvUtil((*p0) - (*p1), vSal, ri);
    return u1;
  }


  void SMPActor::randomize(PRNG* rng, unsigned int numD) {
    sCap = rng->uniform(10.0, 200.0);

    // assign an overall salience, and then by-component saliences
    double s = rng->uniform(0.75, 0.99);
    vSal = KMatrix::uniform(rng, numD, 1, 0.1, 1.0);
    vSal = (s * vSal) / sum(vSal);
    assert(fabs(s - sum(vSal)) < 1E-4);

    // I could randomly assign different voting rules
    // to different actors, but that would just be too cute.
    vr = VotingRule::Proportional;

    return;
  }



  void SMPActor::interpBrgnSnPm(unsigned int n, unsigned int m,
    double tik, double sik, double prbI,
    double tjk, double sjk, double prbJ,
    double & bik, double & bjk) {
    assert((1 == n) || (2 == n));
    assert((1 == m) || (2 == m));

    double wsi = pow(sik, n);
    double wpi = pow(prbI, m);
    double wik = wsi * wpi;

    double wsj = pow(sjk, n);
    double wpj = pow(prbJ, m);
    double wjk = wsj * wpj;

    // imagine that either neither actor cares, or neither actor can coerce the other,
    // so that wik = 0 = wjk. We need to avoid 0/0 error, and have bi=ti and bj=tj.
    const double minW = 1e-6;
    bik = ((wik + minW)*tik + wjk*tjk) / (wik + minW + wjk);
    bjk = (wik*tik + (minW + wjk)*tjk) / (wik + minW + wjk);
    return;
  }


  void SMPActor::interpBrgnS2PMax(double tik, double sik, double prbI,
    double tjk, double sjk, double prbJ,
    double & bik, double & bjk) {
    double di = (prbJ > prbI) ? (prbJ - prbI) : 0;  // max(0, prbJ - prbI);
    double dj = (prbI > prbJ) ? (prbI - prbJ) : 0;  // max(0, prbI - prbJ);
    double sik2 = sik * sik;
    double sjk2 = sjk * sjk;

    const double minW = 1e-6;
    double dik = (di * sjk2) / ((di * sjk2) + minW + ((1 - di) * sik2));
    double djk = (dj * sik2) / ((dj * sik2) + minW + ((1 - dj) * sjk2));


    bik = tik + dik * (tjk - tik);
    bjk = tjk + djk * (tik - tjk);
    return;
  }


  BargainSMP* SMPActor::interpolateBrgn(const SMPActor* ai, const SMPActor* aj,
    const VctrPstn* posI, const VctrPstn * posJ,
    double prbI, double prbJ, InterVecBrgn ivb)  {
    assert((1 == posI->numC()) && (1 == posJ->numC()));
    unsigned int numD = posI->numR();
    assert(numD == posJ->numR());
    auto brgnI = VctrPstn(numD, 1);
    auto brgnJ = VctrPstn(numD, 1);

    for (unsigned int k = 0; k < numD; k++) {
      double tik = (*posI)(k, 0);
      double sik = ai->vSal(k, 0);

      double tjk = (*posJ)(k, 0);
      double sjk = aj->vSal(k, 0);
      double & bik = tik;
      double & bjk = tjk;
      switch (ivb) {
      case InterVecBrgn::S1P1:
        interpBrgnSnPm(1, 1, tik, sik, prbI, tjk, sjk, prbJ, bik, bjk);
        break;
      case InterVecBrgn::S2P2:
        interpBrgnSnPm(2, 2, tik, sik, prbI, tjk, sjk, prbJ, bik, bjk);
        break;
      case InterVecBrgn::S2PMax:
        interpBrgnS2PMax(tik, sik, prbI, tjk, sjk, prbJ, bik, bjk);
        break;
      default:
        throw KBase::KException("interpolateBrgn: unrecognized InterVecBrgn value");
        break;
      }
      brgnI(k, 0) = bik;
      brgnJ(k, 0) = bjk;
    }

    auto brgn = new BargainSMP(ai, aj, brgnI, brgnJ);
    return brgn;
  }


  KMatrix SMPState::bigRfromProb(const KMatrix& p, BigRRange rr)   {
    double pMin = 1.0;
    double pMax = 0.0;
    for (double pi : p) {
      assert(0.0 <= pi);
      assert(pi <= 1.0);
      pMin = (pi < pMin) ? pi : pMin;
      pMax = (pi > pMax) ? pi : pMax;
    }

    const double pTol = 1E-8;
    assert(fabs(1 - KBase::sum(p)) < pTol);

    function<double(unsigned int, unsigned int)> rfn = nullptr;
    switch (rr) {
    case BigRRange::Min:
      rfn = [pMin, pMax, p](unsigned int i, unsigned int j) {
        return (p(i, j) - pMin) / (pMax - pMin);
      };
      break;
    case BigRRange::Mid:
      rfn = [pMin, pMax, p](unsigned int i, unsigned int j) {
        return (3 * p(i, j) - (pMax + 2 * pMin)) / (2 * (pMax - pMin));
      };
      break;
    case BigRRange::Max:
      rfn = [pMin, pMax, p](unsigned int i, unsigned int j) {
        return (2 * p(i, j) - (pMax + pMin)) / (pMax - pMin);
      };
      break;
    }
    auto rMat = KMatrix::map(rfn, p.numR(), p.numC());
    return rMat;
  }


  SMPState::SMPState(Model * m) : State(m) {
    nra = KMatrix();
  }


  SMPState::~SMPState() {
    nra = KMatrix();
  }



  void SMPState::setDiff() {
    auto dfn = [this](unsigned int i, unsigned int j) {
      auto ai = ((const SMPActor*)(model->actrs[i]));
      KMatrix si = ai->vSal;
      auto pi = ((const VctrPstn*)(pstns[i]));
      auto pj = ((const VctrPstn*)(pstns[j]));
      const double d = SMPModel::bvDiff((*pi) - (*pj), si);
      return d;
    };

    const unsigned int na = model->numAct;
    diff = KMatrix::map(dfn, na, na);
    return;
  }


  double SMPState::estNRA(unsigned int h, unsigned int i, SMPState::BigRAdjust ra) const {
    double rh = nra(h, 0);
    double ri = nra(i, 0);
    double rhi = 0.0;
    switch (ra) {
    case SMPState::BigRAdjust::Full:
      rhi = ri;
      break;
    case SMPState::BigRAdjust::Half:
      rhi = (rh + ri) / 2;
      break;
    case SMPState::BigRAdjust::None:
      rhi = rh;
      break;
    }
    return rhi;
  }

  KMatrix SMPState::actrCaps() const {
    auto wFn = [this](unsigned int i, unsigned int j) {
      auto aj = ((SMPActor*)(model->actrs[j]));
      return aj->sCap;
    };

    auto w = KMatrix::map(wFn, 1, model->numAct);
    return w;
  }

  void SMPState::setAUtil(ReportingLevel rl) {
    // you can change these parameters
    auto vr = VotingRule::Proportional;
    auto ra = SMPState::BigRAdjust::Half;
    auto rr = BigRRange::Mid; // use [-0.5, +1.0] scale
    auto vpm = Model::VPModel::Linear;

    const unsigned int na = model->numAct;
    auto w_j = actrCaps();
    setDiff();
    nra = KMatrix(na, 1); // zero-filled, i.e. risk neutral
    auto uFn1 = [this](unsigned int i, unsigned int j) {
      return  SMPModel::bsUtil(diff(i, j), nra(i, 0));
    };

    auto rnUtil_ij = KMatrix::map(uFn1, na, na);

    if (ReportingLevel::Silent < rl) {
      cout << "Raw actor-pos value matrix (risk neutral)" << endl;
      rnUtil_ij.printf(" %+.3f ");
      cout << endl << flush;
    }

    auto pv_ij = Model::vProb(vr, vpm, w_j, rnUtil_ij);
    auto p_i = Model::probCE(pv_ij);
    nra = bigRfromProb(p_i, rr);


    if (ReportingLevel::Silent < rl) {
      cout << "Inferred risk attitudes: " << endl;
      nra.printf(" %+.3f ");
      cout << endl << flush;
    }

    auto raUtil_ij = KMatrix::map(uFn1, na, na);

    if (ReportingLevel::Silent < rl) {
      cout << "Risk-aware actor-pos utility matrix (objective):" << endl;
      raUtil_ij.printf(" %+.4f ");
      cout << endl;
      cout << "RMS change in value vs utility: " << norm(rnUtil_ij - raUtil_ij) / na << endl << flush;
    }

    const double duTol = 1E-6;
    assert(duTol < norm(rnUtil_ij - raUtil_ij)); // I've never seen it below 0.07


    if (ReportingLevel::Silent < rl) {
      switch (ra) {
      case SMPState::BigRAdjust::Full:
        cout << "Using Full adjustment of ra, r^h_i = ri" << endl;
        break;
      case SMPState::BigRAdjust::Half:
        cout << "Using Half adjustment of ra, r^h_i = (rh + ri)/2" << endl;
        break;
      case SMPState::BigRAdjust::None:
        cout << "Using None adjustment of ra, r^h_i = rh " << endl;
        break;
      default:
        cout << "Unrecognized SMPState::BigRAdjust" << endl;
        assert(false);
        break;
      }
    }

    aUtil = vector<KMatrix>();
    for (unsigned int h = 0; h < na; h++) {
      auto u_h_ij = KMatrix(na, na);
      for (unsigned int i = 0; i < na; i++) {
        double rhi = estNRA(h, i, ra);
        for (unsigned int j = 0; j < na; j++) {
          double dij = diff(i, j);
          u_h_ij(i, j) = SMPModel::bsUtil(dij, rhi);
        }
      }
      aUtil.push_back(u_h_ij);


      if (ReportingLevel::Silent < rl) {
        cout << "Estimate by " << h << " of risk-aware utility matrix:" << endl;
        u_h_ij.printf(" %+.4f ");
        cout << endl;

        cout << "RMS change in util^h vs utility: " << norm(u_h_ij - raUtil_ij) / na << endl;
        cout << endl;
      }

      assert(duTol < norm(u_h_ij - raUtil_ij)); // I've never seen it below 0.03
    }
    return;
  }


  void SMPState::showBargains(const vector < vector < BargainSMP* > > & brgns) const {
    for (unsigned int i = 0; i < brgns.size(); i++) {
      printf("Bargains involving actor %i: ", i);
      for (unsigned int j = 0; j < brgns[i].size(); j++) {
        BargainSMP* bij = brgns[i][j];
        if (nullptr != bij) {
          int a1 = model->actrNdx(bij->actInit);
          int a2 = model->actrNdx(bij->actRcvr);
          printf(" [%i:%i] ", a1, a2);
        }
        else {
          printf(" SQ ");
        }
      }
      cout << endl << flush;
    }
    return;
  }


  void SMPState::addPstn(Position* ap) {
    auto sp = (VctrPstn*)ap;
    auto sm = (SMPModel*)model;
    assert(1 == sp->numC());
    assert(sm->numDim == sp->numR());

    State::addPstn(ap);
    return;
  }


  // set the diff matrix, do probCE for risk neutral,
  // estimate Ri, and set all the aUtil[h] matrices 
  SMPState* SMPState::stepBCN() {
    setAUtil(ReportingLevel::Low);
    auto s2 = doBCN();
    s2->step = [s2]() {return s2->stepBCN(); };
    return s2;
  }


  SMPState* SMPState::doBCN() const {
    auto brgns = vector< vector < BargainSMP* > >();
    const unsigned int na = model->numAct;
    brgns.resize(na);
    for (unsigned int i = 0; i < na; i++) {
      brgns[i] = vector<BargainSMP*>();
      brgns[i].push_back(nullptr); // null bargain is SQ
    }

    auto ivb = SMPActor::InterVecBrgn::S2P2;
    // For each actor, identify good targets, and propose bargains to them.
    // (This loop would be an excellent place for high-level parallelism)
    for (unsigned int i = 0; i < na; i++) {
      auto chlgI = bestChallenge(i);
      int bestJ = get<0>(chlgI);
      double piJ = get<1>(chlgI);
      double bestEU = get<2>(chlgI);
      if (0 < bestEU) {
        assert(0 <= bestJ);

        printf("Actor %i has most advantageous target %i worth %.3f\n", i, bestJ, bestEU);

        auto ai = ((const SMPActor*)(model->actrs[i]));
        auto aj = ((const SMPActor*)(model->actrs[bestJ]));
        auto posI = ((const VctrPstn*)pstns[i]);
        auto posJ = ((const VctrPstn*)pstns[bestJ]);
        BargainSMP* brgnIJ = SMPActor::interpolateBrgn(ai, aj, posI, posJ, piJ, 1 - piJ, ivb);
        auto nai = model->actrNdx(brgnIJ->actInit);
        auto naj = model->actrNdx(brgnIJ->actRcvr);

        brgns[i].push_back(brgnIJ); // initiator's copy, delete only it later
        brgns[bestJ].push_back(brgnIJ); // receiver's copy, just null it out later

        printf(" %2i proposes %2i adopt: ", nai, nai);
        KBase::trans(brgnIJ->posInit).printf(" %.3f ");
        printf(" %2i proposes %2i adopt: ", nai, naj);
        KBase::trans(brgnIJ->posRcvr).printf(" %.3f ");
      }
      else {
        printf("Actor %i has no advantageous targets \n", i);
      }
    }


    cout << endl << "Bargains to be resolved" << endl << flush;
    showBargains(brgns);

    auto w = actrCaps();
    cout << "w:" << endl;
    w.printf(" %6.2f ");

    // of course, you  can change these two parameters
    auto vr = VotingRule::Proportional;
    auto vpm = Model::VPModel::Linear;

    auto ndxMaxProb = [](const KMatrix & cv) {

      const double pTol = 1E-8;
      assert(fabs(KBase::sum(cv) - 1.0) < pTol);

      assert(0 < cv.numR());
      assert(1 == cv.numC());

      auto ndxIJ = ndxMaxAbs(cv);
      unsigned int iMax = get<0>(ndxIJ);
      return iMax;
    };

    // what is the utility to actor nai of the state resulting after
    // the nbj-th bargain of the k-th actor is implemented?
    auto brgnUtil = [this, brgns](unsigned int nk, unsigned int nai, unsigned int nbj) {
      const unsigned int na = model->numAct;
      BargainSMP * b = brgns[nk][nbj];
      double uAvrg = 0.0;

      if (nullptr == b) { // SQ bargain
        uAvrg = 0.0;
        for (unsigned int n = 0; n < na; n++) {
          // nai's estimate of the utility to nai of position n, i.e. the true value
          uAvrg = uAvrg + aUtil[nai](nai, n);
        }
      }

      if (nullptr != b) { // all positions unchanged, except Init and Rcvr
        uAvrg = 0.0;
        auto ndxInit = model->actrNdx(b->actInit);
        assert((0 <= ndxInit) && (ndxInit < na)); // must find it
        double uPosInit = ((SMPActor*)(model->actrs[nai]))->posUtil(&(b->posInit), this);
        uAvrg = uAvrg + uPosInit;

        auto ndxRcvr = model->actrNdx(b->actRcvr);
        assert((0 <= ndxRcvr) && (ndxRcvr < na)); // must find it
        double uPosRcvr = ((SMPActor*)(model->actrs[nai]))->posUtil(&(b->posRcvr), this);
        uAvrg = uAvrg + uPosRcvr;

        for (unsigned int n = 0; n < na; n++) {
          if ((ndxInit != n) && (ndxRcvr != n)) {
            // again, nai's estimate of the utility to nai of position n, i.e. the true value
            uAvrg = uAvrg + aUtil[nai](nai, n);
          }
        }
      }

      uAvrg = uAvrg / na;

      assert(0.0 < uAvrg); // none negative, at least own is positive
      assert(uAvrg <= 1.0); // can not all be over 1.0
      return uAvrg;
    };
    // end of λ-fn


    // TODO: finish this
    // For each actor, assess what bargains result from CDMP, and put it into s2

    // The key is to build the usual matrix of U_ai (Brgn_m) for all bargains in brgns[k],
    // making sure to divide the sum of the utilities of positions by 1/N
    // so 0 <= Util(state after Brgn_m) <= 1, then do the standard scalarPCE for bargains involving k.


    SMPState* s2 = new SMPState(model);
    // (This loop would be a good place for high-level parallelism)
    for (unsigned int k = 0; k < na; k++) {
      unsigned int nb = brgns[k].size();
      auto buk = [brgnUtil, k](unsigned int nai, unsigned int nbj) { return brgnUtil(k, nai, nbj); };
      auto u_im = KMatrix::map(buk, na, nb);

      cout << "u_im: " << endl;
      u_im.printf(" %.5f ");

      cout << "Doing probCE for the " << nb << " bargains of actor " << k << " ... " << flush;
      auto p = Model::scalarPCE(na, nb, w, u_im, vr, vpm, ReportingLevel::Medium);
      assert(nb == p.numR());
      assert(1 == p.numC());
      cout << "done" << endl << flush;
      unsigned int mMax = ndxMaxProb(p); // indexing actors by i, bargains by m
      cout << "Chosen bargain: " << mMax << endl;



      // TODO: create a fresh position for k, from the selected bargain mMax.
      VctrPstn * pk = nullptr;
      auto bkm = brgns[k][mMax];
      if (nullptr == bkm) {
        auto oldPK = (VctrPstn *)pstns[k];
        pk = new VctrPstn(*oldPK);
      }
      else {
        const unsigned int ndxInit = model->actrNdx(bkm->actInit);
        const unsigned int ndxRcvr = model->actrNdx(bkm->actRcvr);
        if (ndxInit == k) {
          pk = new VctrPstn(bkm->posInit);
        }
        else if (ndxRcvr == k) {
          pk = new VctrPstn(bkm->posRcvr);
        }
        else {
          cout << "SMPState::doBCN: unrecognized actor in bargain" << endl;
          assert(false);
        }
      }
      assert(nullptr != pk);

      assert(k == s2->pstns.size());
      s2->pstns.push_back(pk);

      cout << endl << flush;
    }


    // Some bargains are nullptr, and there are two copies of every non-nullptr randomly
    // arranged. If we delete them as we find them, then the second occurance will be corrupted,
    // so the code crashes when it tries to access the memory to see if it matches something
    // already deleted. Hence, we scan them all, building a list of unique bargains,
    // then delete those in order. 
    auto uBrgns = vector<BargainSMP*>();

    for (unsigned int i = 0; i < brgns.size(); i++) {
      auto ai = ((const SMPActor*)(model->actrs[i]));
      for (unsigned int j = 0; j < brgns[i].size(); j++) {
        BargainSMP* bij = brgns[i][j];
        if (nullptr != bij) {
          if (ai == bij->actInit) {
            uBrgns.push_back(bij); // this is the initiator's copy, so save it for deletion
          }
        }
        brgns[i][j] = nullptr; // either way, null it out.
      }
    }

    for (auto b : uBrgns) {
      int aI = model->actrNdx(b->actInit);
      int aR = model->actrNdx(b->actRcvr);
      //printf("Delete bargain [%2i:%2i] \n", aI, aR);
      delete b;
    }

    return s2;
  }


  // h's estimate of the victory probability and expected change in utility for k from i challenging j,
  // compared to status quo.
  // Note that the  aUtil vector of KMatrix must be set before starting this.
  // TODO: offer a choice the different ways of estimating value-of-a-state: even sum or expected value.
  // TODO: we may need to separate euConflict from this at some point
  tuple<double, double> SMPState::probEduChlg(unsigned int h, unsigned int k, unsigned int i, unsigned int j) const {

    // you could make other choices for these two sub-models
    auto vr = VotingRule::Proportional;
    auto tpc = KBase::ThirdPartyCommit::Semi;

    double uii = aUtil[h](i, i);
    double uij = aUtil[h](i, j);
    double uji = aUtil[h](j, i);
    double ujj = aUtil[h](j, j);

    // h's estimate of utility to k of status-quo positions of i and j
    double euSQ = aUtil[h](k, i) + aUtil[h](k, j);
    assert((0.0 <= euSQ) && (euSQ <= 2.0));

    // h's estimate of utility to k of i defeating j, so j adopts i's position
    double uhkij = aUtil[h](k, i) + aUtil[h](k, i);
    assert((0.0 <= uhkij) && (uhkij <= 2.0));

    // h's estimate of utility to k of j defeating i, so i adopts j's position
    double uhkji = aUtil[h](k, j) + aUtil[h](k, j);
    assert((0.0 <= uhkji) && (uhkji <= 2.0));

    auto ai = ((const SMPActor*)(model->actrs[i]));
    double si = KBase::sum(ai->vSal);
    double ci = ai->sCap;
    auto aj = ((const SMPActor*)(model->actrs[j]));
    double sj = KBase::sum(aj->vSal);
    assert((0 < sj) && (sj <= 1));
    double cj = aj->sCap;
    double minCltn = 1E-10;

    // h's estimate of i's unilateral influence contribution to (i:j), hence positive
    double contrib_i_ij = Model::vote(vr, si*ci, uii, uij);
    assert(0 <= contrib_i_ij);
    double chij = minCltn + contrib_i_ij; // strength of complete coalition supporting i over j
    assert (0.0 < chij);

    // h's estimate of j's unilateral influence contribution to (i:j), hence negative
    double contrib_j_ij = Model::vote(vr, sj*cj, uji, ujj);
    assert(contrib_j_ij <= 0);
    double chji = minCltn - contrib_j_ij; // strength of complete coalition supporting j over i
    assert (0.0 < chji);

    const unsigned int na = model->numAct;

    // we assess the overall coalition strengths by adding up the contribution of
    // individual actors (including i and j, above). We assess the contribution of third 
    // parties by looking at little coalitions in the hypothetical (in:j) or (i:nj) contests.
    for (unsigned int n = 0; n < na; n++) {
      if ((n != i) && (n != j)) { // already got their influence-contributions
        auto an = ((const SMPActor*)(model->actrs[n]));
        double cn = an->sCap;
        double sn = KBase::sum(an->vSal);
        double uni = aUtil[h](n, i);
        double unj = aUtil[h](n, j);
        double unn = aUtil[h](n, n);

        double pin = Actor::vProbLittle(vr, sn*cn, uni, unj, chij, chji);
        assert(0.0 <= pin);
        assert(pin <= 1.0);
        double pjn = 1.0 - pin;

        double vnij = Actor::thirdPartyVoteSU(sn*cn, vr, tpc, pin, pjn, uni, unj, unn);

        chij = (vnij > 0) ? (chij + vnij) : chij;
        assert(0 < chij);
        chji = (vnij < 0) ? (chji - vnij) : chji;
        assert(0 < chji);
      }
    }

    const double phij = chij / (chij + chji);
    const double phji = chji / (chij + chji);

    const double euCh = (1 - sj)*uhkij + sj*(phij*uhkij + phji*uhkji);
    const double euChlg = euCh - euSQ;
    // printf ("SMPState::probEduChlg(%2i, %2i, %2i, %i2) = %+6.4f - %+6.4f = %+6.4f\n", h, k, i, j, euCh, euSQ, euChlg);
    auto rslt = tuple<double, double>(phij, euChlg);
    return rslt;
  }



  tuple<int, double, double> SMPState::bestChallenge(unsigned int i) const {
    int bestJ = -1;
    double piJ = 0;
    double bestEU = 0;

    // for SMP, positive ej are typically in the 0.5 to 0.01 range, so I take 1/1000 of the minimum, 
    const double minSig = 1e-5;

    for (unsigned int j = 0; j < model->numAct; j++) {
      if (j != i) {
        auto pej = probEduChlg(i, i, i, j);
        double pj = get<0>(pej); // i's estimate of the victory-Prob for i challengeing j
        double ej = get<1>(pej); // i's estimate of the change in utility to i of i challengeing j, compared to SQ
        if ((minSig < ej) && (bestEU < ej)) {
          bestJ = j;
          piJ = pj;
          bestEU = ej;
        }
      }
    }
    auto rslt = tuple<int, double, double>(bestJ, piJ, bestEU);
    return rslt;
  }



  KMatrix SMPState::pDist(int persp) const {
    auto na = model->numAct;
    auto w = actrCaps();
    auto vr = VotingRule::Proportional;
    auto rl = ReportingLevel::Silent;
    auto vpm = Model::VPModel::Linear;
    auto uij = KMatrix(na, na);

    if ((0 <= persp) && (persp < na)) {
      uij = aUtil[persp];
    }
    else if (-1 == persp) {
      for (unsigned int i = 0; i < na; i++) {
        for (unsigned int j = 0; j < na; j++) {
          auto ui = aUtil[i];
          uij(i, j) = aUtil[i](i, j);
        }
      }
    }
    else {
      cout << "SMPState::pDist: unrecognized perspective, " << persp << endl << flush;
      assert(false);
    }
    auto pd = Model::scalarPCE(na, na, w, uij, vr, vpm, rl);
    return pd;
  }


  // -------------------------------------------------


  SMPModel::SMPModel(PRNG * r) : Model(r) {
    numDim = 0;
    dimName = vector<string>();
  }

  SMPModel::~SMPModel() {
    //
  }

  void SMPModel::addDim(string dn) {
    dimName.push_back(dn);
    numDim = dimName.size();
    return;
  }

double SMPModel::stateDist(const SMPState* s1 , const SMPState* s2 ) {
  unsigned int n = s1->pstns.size();
  assert (n == s2->pstns.size());
  double dSum = 0;
  for (unsigned int i=0; i<n; i++) {
    auto vp1i = ((const VctrPstn*) (s1->pstns[i]));
    auto vp2i = ((const VctrPstn*) (s2->pstns[i]));
    dSum = dSum + KBase::norm((*vp1i) - (*vp2i));
  }
  return dSum;
}


  // 0 <= d <= 1 is the difference in normalized position
  // -1 <= R <= +1 is normalized risk-aversion
  double SMPModel::bsUtil(double d, double R){
    double u = 0;
    assert(0 <= d);
    if (d <= 1) {
      u = (1 - d)*(1 + d*R);  //  (0 <= u) && (u <= 1)
    }
    else {
      // linearly interpolate with last util-slope at d=1,
      // This is "unphysical", but we have to deal with it anyway
      // because VHCSearch will vary components out of physical limits
      // for both scalar and vector cases.
      double uSlope = -(R + 1); // obviously, uSlope <= 0
      u = uSlope*(d - 1); // because u=0 at d=1, regardless of R, this u<0
    }
    return u;
  }

  double SMPModel::bvDiff(const  KMatrix & d, const  KMatrix & s) {
    assert(KBase::sameShape(d, s));
    double dsSqr = 0;
    double ssSqr = 0;
    for (unsigned int i = 0; i < d.numR(); i++) {
      for (unsigned int j = 0; j < d.numC(); j++) {
        const double dij = d(i, j);
        const double sij = s(i, j);
        assert(0 <= sij);
        const double ds = dij * sij;
        const double s = sij;
        dsSqr = dsSqr + (ds*ds);
        ssSqr = ssSqr + (s*s);
      }
    }
    assert(0 < ssSqr);
    double sd = sqrt(dsSqr / ssSqr);
    return sd;
  };

  double SMPModel::bvUtil(const  KMatrix & d, const  KMatrix & s, double R) {
    const double sd = bvDiff(d, s);
    const double u = bsUtil(sd, R);
    return u;
  };


  void SMPModel::showVPHistory() const {
    const string cs = " , ";

    // show positions over time
    for (unsigned int i = 0; i < numAct; i++) {
      for (unsigned int k = 0; k < numDim; k++) {
        cout << actrs[i]->name << cs;
        printf("Dim-%02i %s", k, cs.c_str());
        //cout << "Dim-"<< k << cs;
        for (unsigned int t = 0; t < history.size(); t++) {
          State* st = history[t];
          Position* pit = st->pstns[i];
          VctrPstn* vpit = (VctrPstn*)pit;
          assert(1 == vpit->numC());
          assert(numDim == vpit->numR());
          printf("%7.3f %s", 100 * (*vpit)(k, 0), cs.c_str());
        }
        cout << endl;
      }
    }
    cout << endl;
    // show probabilities over time.
    // Note that we have to set the aUtil matrices for the last one.
    auto prbHist = vector<KMatrix>();
    for (unsigned int t = 0; t < history.size(); t++) {
      auto sst = (SMPState*)history[t];
      if (t == history.size() - 1) {
        sst->setAUtil(ReportingLevel::Silent);
      }
      auto pdt = sst->pDist(-1);
      prbHist.push_back(pdt);
    }

    for (unsigned int i = 0; i < numAct; i++) {
      cout << actrs[i]->name << cs;
      printf("prob %s", cs.c_str());
      for (unsigned int t = 0; t < history.size(); t++) {
        printf("%.4f %s", prbHist[t](i, 0), cs.c_str());
      }
      cout << endl << flush;
    }
    return;
  }


  SMPModel * readCSV(string fName, PRNG * rng) {
    using KBase::KException;
    auto newChars = [](unsigned int len) {
      auto rslt = new char[len];
      for (unsigned int i = 0; i < len; i++) {
        rslt[i] = ((char)0);
      }
      return rslt;
    };

    const unsigned int minNumActor = 3;
    const unsigned int maxNumActor = 100; // It's just a demo 
    char * errBuff = newChars(100); // as sprintf requires

    csv_parser csv(fName);

    // Get values according to row number and column number.
    // Remember it starts from (1,1) and not (0,0)
    string scenName = csv.get_value(1, 1);
    cout << "Scenario name: |" << scenName << "|" << endl;
    cout << flush;
    string numActorString = csv.get_value(1, 3);
    unsigned int numActor = atoi(numActorString.c_str());
    string numDimString = csv.get_value(1, 4);
    int numDim = atoi(numDimString.c_str());
    printf("Number of actors: %i \n", numActor);
    printf("Number of dimensions: %i \n", numDim);
    cout << endl << flush;

    if (numDim < 1) { // lower limit
      throw(KBase::KException("SMPModel::readCSV: Invalid number of dimensions"));
    }
    assert(0 < numDim);
    if ((numActor < minNumActor) || (maxNumActor < numActor)) { // avoid impossibly low or ridiculously large
      throw(KBase::KException("SMPModel::readCSV: Invalid number of actors"));
    }
    assert(minNumActor <= numActor);
    assert(numActor <= maxNumActor);

    // Read actor data
    auto actorNames = vector<string>();
    auto actorDescs = vector<string>();
    auto cap = KMatrix(numActor, 1);
    auto nra = KMatrix(numActor, 1);
    for (unsigned int i = 0; i < numActor; i++) {
      // get short names
      string nis = csv.get_value(3 + i, 1);
      assert(0 < nis.length());
      actorNames.push_back(nis);
      printf("Actor %3i name: %s \n", i, actorNames[i].c_str());

      // get long descriptions
      string descsi = csv.get_value(3 + i, 2);
      actorDescs.push_back(descsi);
      printf("Actor %3i desc: %s \n", i, actorDescs[i].c_str());

      // get capability/power, often on 0-100 scale
      string psi = csv.get_value(3 + i, 3);
      double pi = atof(psi.c_str());
      printf("Actor %3i power: %5.1f \n", i, pi);
      assert(0 <= pi); // zero weight is pointless, but not incorrect
      assert(pi < 1E8); // no real upper limit, so this is just a sanity-check
      cap(i, 0) = pi;

   
      cout << endl << flush;

    } // loop over actors, i


    // get issue names
    auto dNames = vector<string>();
    for (unsigned int j = 0; j < numDim; j++) {
      string insi = csv.get_value(2, 4 + 2 * j);
      dNames.push_back(insi);
      printf("Dimension %2i: %s \n", j, dNames[j].c_str());
    }
    cout << endl;

    // get position/salience data
    auto pos = KMatrix(numActor, numDim);
    auto sal = KMatrix(numActor, numDim);
    for (unsigned int i = 0; i < numActor; i++) {
      double salI = 0.0;
      for (unsigned int j = 0; j < numDim; j++) {
        string posSIJ = csv.get_value(3 + i, 4 + 2 * j);
        double posIJ = atof(posSIJ.c_str());
        printf("pos[%3i , %3i] =  %5.3f \n", i, j, posIJ);
        cout << flush;
        if ((posIJ < 0.0) || (+100.0 < posIJ)) { // lower and upper limit
          errBuff = newChars(100);
          sprintf(errBuff, "SMPModel::readCSV: Out-of-bounds position for actor %i on dimension %i:  %f",
            i, j, posIJ);
          throw(KException(errBuff));
        }
        assert(0.0 <= posIJ);
        assert(posIJ <= 100.0);
        pos(i, j) = posIJ;

        string salSIJ = csv.get_value(3 + i, 5 + 2 * j);
        double salIJ = atof(salSIJ.c_str());
        //printf("sal[%3i , %3i] = %5.3f \n", i, j, salIJ);
        //cout << flush;
        if ((salIJ < 0.0) || (+100.0 < salIJ)) { // lower and upper limit
          errBuff = newChars(100);
          sprintf(errBuff, "SMPModel::readCSV: Out-of-bounds salience for actor %i on dimension %i:  %f",
            i, j, salIJ);
          throw(KException(errBuff));
        }
        assert(0.0 <= salIJ);
        salI = salI + salIJ;
        //printf("sal[%3i] = %5.3f \n", i, salI);
        //cout << flush;
        if (+100.0 < salI) { // upper limit: no more than 100% of attention to all issues
          errBuff = newChars(100);
          sprintf(errBuff,
            "SMPModel::readCSV: Out-of-bounds total salience for actor %i:  %f",
            i, salI);
          throw(KException(errBuff));
        }
        assert(salI <= 100.0);
        sal(i, j) = (double)salIJ;
        //cout << endl << flush;
      }
    }

    cout << "Position matrix:" << endl;
    pos.printf("%5.1f  ");
    cout << endl << endl << flush;
    cout << "Salience matrix:" << endl;
    sal.printf("%5.1f  ");
    cout << endl << flush;

    // get them into the proper internal scale:
    pos = pos / 100.0;
    sal = sal / 100.0;

    // now that it is read and verified, use the data
    auto sm0 = SMPModel::initModel(actorNames, actorDescs, dNames, cap, pos, sal, rng);
    return sm0;
  }

  SMPModel * SMPModel::initModel(vector<string> aName, vector<string> aDesc, vector<string> dName,
    KMatrix cap, KMatrix pos, KMatrix sal, PRNG * rng) {
    SMPModel * sm0 = new SMPModel(rng);
    SMPState * st0 = new SMPState(sm0);
    st0->step = [st0]() {return st0->stepBCN(); };
    sm0->addState(st0);


    const unsigned int na = aName.size();
    const unsigned int nd = dName.size();


    for (auto dn : dName) {
      sm0->addDim(dn);
    }

    for (unsigned int i = 0; i < na; i++) {
      auto ai = new SMPActor(aName[i], aDesc[i]);
      ai->sCap = cap(i, 0);
      ai->vSal = KMatrix(nd, 1);
      auto vpi = new VctrPstn(nd, 1);
      for (unsigned int j = 0; j < nd; j++) {
        ai->vSal(j, 0) = sal(i, j);
        (*vpi)(j, 0) = pos(i, j);
      }
      sm0->addActor(ai);
      st0->addPstn(vpi);
    }

    return sm0;
  }

  // -------------------------------------------------
  void demoActorUtils(uint64_t s, PRNG* rng) {

    printf("Using PRNG seed: %020lu \n", s);
    rng->setSeed(s);

    cout << "Demonstrate simple voting by spatial actors (scalar capability)" << endl;
    cout << "and by economic actors (vector capability)" << endl;

    unsigned int sDim = 3;
    auto sp1 = new VctrPstn(KMatrix::uniform(rng, sDim, 1, 0.0, 1.0));
    cout << "Random spatial position, sp1:" << endl;
    sp1->printf(" %.3f ");
    cout << endl << flush;
    auto sp2 = new VctrPstn(KMatrix::uniform(rng, sDim, 1, 0.0, 1.0));
    cout << "Random spatial position, sp2:" << endl;
    sp2->printf(" %.3f ");
    cout << endl << flush;
    auto alice = new SMPActor("Alice", "first cryptographer");
    alice->randomize(rng, sDim);
    auto md0 = new SMPModel(rng);
    auto st0 = new SMPState(md0);
    st0->nra = KMatrix(1, 1); // 
    auto iPos = new VctrPstn((2 * (*sp1) + (*sp2)) / 3);
    md0->addActor(alice);
    st0->addPstn(iPos);
    cout << "Alice's position is (2*sp1 + sp2)/3:" << endl;
    iPos->printf(" %.3f ");
    cout << endl << flush;
    printf("Alice's scalar capability: %.3f \n", alice->sCap);
    cout << "Alice's voting rule (overall): " << vrName(alice->vr) << endl;
    printf("Alice's risk attitude: %.3f \n", st0->nra(0, 0));
    printf("Alice's total salience %.4f \n", sum(alice->vSal));
    printf("Alice's vector salience: \n");
    alice->vSal.printf(" %.3f ");


    double va = alice->vote(sp1, sp2, st0);
    printf("A's vote on [sp1:sp2] is %+.3f \n", va);
    printf("Her vote should always be positive \n");
    cout << flush;
    assert(va > 0); // by construction
    cout << endl << flush;

    delete sp1;
    delete sp2;
    delete alice;
    delete st0;
    delete md0;


    return;
  }


  void demoEUSpatial(unsigned int numA, unsigned int sDim, uint64_t s, PRNG* rng) {
    assert(0 < sDim);
    assert(2 < numA);

    printf("Using PRNG seed: %020lu \n", s);
    rng->setSeed(s);

    cout << "EU State for SMP actors with scalar capabilities" << endl;
    printf("Number of actors; %i \n", numA);
    printf("Number of SMP dimensions %i \n", sDim);

    // note that because all actors use the same scale for capability, utility, etc,
    // their 'votes' are on the same scale and influence can be added up meaningfully

    const unsigned int maxIter = 5000;
      double qf = 100.0;
    auto md0 = new SMPModel(rng);
    md0->stop = [maxIter] (unsigned int iter, const State * s) { return (maxIter <= iter); };
    md0->stop = [maxIter, qf] (unsigned int iter, const State * s) { 
      bool tooLong = (maxIter <= iter);
      bool quiet = false;
      if (1 < iter) {
        auto sf = [](unsigned int i1, unsigned int i2, double d12) {
        printf("sDist [%2i,%2i] = %.2E   ", i1, i2, d12);
        return;
        };
        
        auto s0 = ((const SMPState*) (s->model->history[0]));
        auto s1 = ((const SMPState*) (s->model->history[1]));
        auto d01 = SMPModel::stateDist(s0, s1);
        sf(0, 1, d01);
        
        auto sx = ((const SMPState*) (s->model->history[iter-0]));
        auto sy = ((const SMPState*) (s->model->history[iter-1]));
        auto dxy = SMPModel::stateDist(sx, sy);
        sf( iter-0, iter-1, dxy);
        
        quiet = (dxy < d01/qf); 
        if (quiet) 
          printf("Quiet \n");
        else
          printf("Not Quiet \n");
        
        cout << endl<<flush;
      }
      return tooLong || quiet;
      
    };
    md0->numDim = sDim;

    auto st0 = new SMPState(md0);
    md0->addState(st0); // now state 0 of the history

    st0->step = [st0]() {return st0->stepBCN(); };

    for (unsigned int i = 0; i < numA; i++){
      //string ni = "SActor-";
      //ni.append(std::to_string(i));
      unsigned int nbSize = 15;
      char * nameBuff = new char[nbSize];
      for (unsigned int j = 0; j < nbSize; j++) { nameBuff[j] = (char)0; }
      sprintf(nameBuff, "SActor-%02i", i);
      auto ni = string(nameBuff);
      delete nameBuff;
      nameBuff = nullptr;
      string di = "Random spatial actor";


      auto ai = new SMPActor(ni, di);
      ai->randomize(rng, sDim);
      auto iPos = new VctrPstn(KMatrix::uniform(rng, sDim, 1, 0.0, 1.0)); // SMP is always on [0,1] scale
      md0->addActor(ai);
      st0->addPstn(iPos);
    }
    st0->nra = KMatrix::uniform(rng, numA, 1, -0.5, 1.0);

    for (unsigned int i = 0; i < numA; i++){
      auto ai = ((SMPActor*)(md0->actrs[i]));
      auto ri = st0->nra(i, 0);
      printf("%2i: %s , %s \n", i, ai->name.c_str(), ai->desc.c_str());
      cout << "voting rule: " << vrName(ai->vr) << endl;
      cout << "Pos vector: ";
      VctrPstn * pi = ((VctrPstn*)(st0->pstns[i]));
      trans(*pi).printf(" %+7.4f ");
      cout << "Sal vector: ";
      trans(ai->vSal).printf(" %+7.4f ");
      printf("Capability: %.3f \n", ai->sCap);
      printf("Risk attitude: %+.4f \n", ri);
      cout << endl;
    }

    // with SMP actors, we can always read their ideal position.
    // with strategic voting, they might want to advocate positions
    // separate from their ideal, but this simple demo skips that.
    auto uFn1 = [st0](unsigned int i, unsigned int j) {
      auto ai = ((SMPActor*)(st0->model->actrs[i]));
      auto pj = ((VctrPstn*)(st0->pstns[j])); // aj->iPos;
      double uij = ai->posUtil(pj, st0);
      return uij;
    };

    auto u = KMatrix::map(uFn1, numA, numA);
    cout << "Raw actor-pos util matrix" << endl;
    u.printf(" %.4f ");
    cout << endl << flush;

    auto w = st0->actrCaps(); //  KMatrix::map(wFn, 1, numA);

    // arbitrary but illustrates that we can do an election with arbitrary
    // voting rules - not necessarily the same as the actors would do.
    auto vr = VotingRule::Binary;
    cout << "Using voting rule " << vrName(vr) << endl;

    auto vpm = Model::VPModel::Linear;

    KMatrix p = Model::scalarPCE(numA, numA, w, u, vr, vpm, ReportingLevel::Medium);

    cout << "Expected utility to actors: " << endl;
    (u*p).printf(" %.3f ");
    cout << endl << flush;

    cout << "Net support for positions: " << endl;
    (w*u).printf(" %.3f ");
    cout << endl << flush;

    auto aCorr = [](const KMatrix & x, const KMatrix & y) {
      using KBase::lCorr;
      using KBase::mean;
      return  lCorr(x - mean(x), y - mean(y));
    };

    // for nearly flat distributions, and nearly flat net support,
    // one can sometimes see negative affine-correlations because of
    // random variations in 3rd or 4th decimal places.
    printf("L-corr of prob and net support: %+.4f \n", KBase::lCorr((w*u), trans(p)));
    printf("A-corr of prob and net support: %+.4f \n", aCorr((w*u), trans(p)));

    // no longer need external reference to the state
    st0 = nullptr;

    cout << "Starting model run" << endl << flush;
    md0->run();

    cout << "Completed model run" << endl << endl;

    cout << "History of actor positions over time" << endl;
    md0->showVPHistory();

    cout << endl;
    cout << "Delete model (actors, states, positions, etc.)" << endl << flush;
    delete md0;
    md0 = nullptr;

    return;
  }

  void readEUSpatial(uint64_t seed, string inputCSV, PRNG* rng){
    auto md0 = readCSV(inputCSV, rng);

    const unsigned int maxIter = 5;
    md0->stop = [maxIter](unsigned int iter, const State * s) { return (maxIter <= iter); };

    cout << "Starting model run" << endl << flush;
    md0->run();

    cout << "Completed model run" << endl << endl;

    cout << "History of actor positions over time" << endl;
    md0->showVPHistory();

    delete md0;
    return;
  }

} // namespace


int main(int ac, char **av) {
  using std::cout;
  using std::endl;
  using std::string;

  auto sTime = KBase::displayProgramStart();
  uint64_t dSeed = 0x89B9E56731783372; // 0x9B3BDF88E9E466D3; // arbitrary 
  uint64_t seed = dSeed;
  bool run = true;
  bool euSmpP = false;
  bool csvP = false;
  string inputCSV = "";

  cout << "smpApp version " << DemoSMP::appVersion << endl << endl;

  auto showHelp = [dSeed]() {
    printf("\n");
    printf("Usage: specify one or more of these options\n");
    printf("--help            print this message\n");
    printf("--euSMP           exp. util. of spatial model of politics\n");
    printf("--csv <f>         read a scenario from CSV\n");
    printf("--seed <n>        set a 64bit seed\n");
    printf("                  0 means truly random\n");
    printf("                  default: %020lu \n", dSeed);
  };

  // tmp args
  //seed = 0xC4CB0A3E298562EA;
  //euSmpP = true;
  //csvP = true;
  //inputCSV = "doc/dummyData_3Dim.csv";
  seed =  0x10299382392D026D;//  01164624166166266477;

  if (ac > 1) {
    for (int i = 1; i < ac; i++) {
      if (strcmp(av[i], "--seed") == 0) {
        i++;
        seed = std::stoull(av[i]);
      }
      else if (strcmp(av[i], "--csv") == 0) {
        csvP = true;
        i++;
        inputCSV = av[i];
      }
      else if (strcmp(av[i], "--euSMP") == 0) {
        euSmpP = true;
      }
      else if (strcmp(av[i], "--help") == 0) {
        run = false;
      }
      else {
        run = false;
        printf("Unrecognized argument %s\n", av[i]);
      }
    }
  }

  if (!run){
    showHelp();
    return 0;
  }

  PRNG * rng = new PRNG();
  seed = rng->setSeed(seed); // 0 == get a random number
  printf("Using PRNG seed:  %020lu \n", seed);
  printf("Same seed in hex:   0x%016lX \n", seed);

  // note that we reset the seed every time, so that in case something 
  // goes wrong, we need not scroll back too far to find the
  // seed required to reproduce the bug.  
  if (euSmpP) {
    cout << "-----------------------------------" << endl;
    DemoSMP::demoEUSpatial(7, 3, seed, rng);
  }
  if (csvP) {
    cout << "-----------------------------------" << endl;
    DemoSMP::readEUSpatial(seed, inputCSV, rng);
  }
  cout << "-----------------------------------" << endl;


  delete rng;
  KBase::displayProgramEnd(sTime);
  return 0;
}

// --------------------------------------------
// Copyright KAPSARC. Open source MIT License.
// --------------------------------------------
