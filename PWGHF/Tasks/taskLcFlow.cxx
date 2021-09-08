// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file taskLcFlow.cxx
/// \task to calculate flow of Lambda_c
///
/// \author Katarina Krizkova Gajdosova <katarina.gajdosova@cern.ch>, Czech Technical University in Prague

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include <TComplex.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_prong3;
using namespace o2::analysis::hf_cuts_lc_topkpi;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, false, {"Fill MC histograms."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// Λc± → p± K∓ π± flow analysis task
struct TaskLcFlow {
  HistogramRegistry registry{
    "registry",
    {
      {"hPtHadrons", "charged hadrons; #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"hPtHadronsSelected", "charged hadrons after cuts; #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"hEtaHadrons", "charged hadrons; #eta; entries", {HistType::kTH1F, {{100, -4., 4.}}}},
      {"hEtaHadronsSelected", "charged hadrons after cuts; #eta; entries", {HistType::kTH1F, {{100, -4., 4.}}}},
      {"hPhiHadronsSelected", "charged hadrons after cuts; #varphi; entries", {HistType::kTH1F, {{60, 0., 2*TMath::Pi()}}}},
      {"hPhiEtaHadronsSelected", "charged hadrons after cuts; #varphi; #eta", {HistType::kTH2F, {{30, 0., 2*TMath::Pi()}, {40, -2., 2.}}}},
      {"hDCAzHadronsSelected", "charged hadrons after cuts; DCA_{z} (cm); entries", {HistType::kTH1F, {{60, -3., 3.}}}},
      {"hPhiEtaCand", "candidates; #varphi; #eta", {HistType::kTH2F, {{30, 0., 2*TMath::Pi()}, {40, -2., 2.}}}},
      {"hPtCand", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"hNtracks", "charged hadrons; Number of tracks; entries", {HistType::kTH1F, {{100, 0., 10000}}}},
      {"hcorrRef", "charged hadron <<2>>; Number of tracks; <<2>>", {HistType::kTProfile, {{100, 0., 10000}}}},
      {"hcorrRefGap", "charged hadron <<2>> |#Delta #eta| > 0.0; Number of tracks; <<2>> |#Delta #eta| > 0.0", {HistType::kTProfile, {{100, 0., 10000}}}}
      }};

  Configurable<int> d_selectionFlagLc{"d_selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_lc_topkpi::pTBins_v}, "pT bin limits"};

  Filter filterSelectCandidates = (aod::hf_selcandidate_lc::isSelLcpKpi >= d_selectionFlagLc || aod::hf_selcandidate_lc::isSelLcpiKp >= d_selectionFlagLc);

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)bins;
    registry.add("hmass", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 1.6, 3.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hmass2", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 1.6, 3.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hmassEtaPos", "3-prong candidates, |#eta|>0;inv. mass (p K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 1.6, 3.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hmassEtaNeg", "3-prong candidates, |#eta|<0;inv. mass (p K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 1.6, 3.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hcorrDiff", "charged hadrons <<2'>>; p_{T} (GeV/#it{c})", {HistType::kTProfile, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hcorrDiffGapNeg", "charged hadrons <<2'>> |#Delta#eta|>0, #eta(p-vector)<0; p_{T} (GeV/#it{c})", {HistType::kTProfile, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hcorrDiffGapPos", "charged hadrons <<2'>> |#Delta#eta|>0, #eta(p-vector)>0; p_{T} (GeV/#it{c})", {HistType::kTProfile, {{vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hcorrLcDiffGapNeg", "Lc <<2'>> |#Delta#eta|>0, #eta(p-vector)<0;inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T} (GeV/#it{c})", {HistType::kTProfile2D, {{30, 1.6, 3.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hcorrLcDiffGapPos", "Lc <<2'>> |#Delta#eta|>0, #eta(p-vector)>0;inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T} (GeV/#it{c})", {HistType::kTProfile2D, {{30, 1.6, 3.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  /// aod::BigTracks is not soa::Filtered, should be added when filters are added
  void process(aod::Collision const& collisions, soa::Join<aod::FullTracks, aod::TracksExtended> const& tracks, soa::Filtered<soa::Join<aod::HfCandProng3, aod::HFSelLcCandidate>> const& candidates)
  {

    auto vbins = (std::vector<double>)bins;
    const int sizePt = 10;
    const int sizeMinv = 30;
    double minvarray[31] = {0};
    for(int i = 0; i < 31; i++) {
      minvarray[i] = 1.6 + 0.05*i;
    }

    int Ntracks = 0;
    int NtracksGapNegative = 0;
    int NtracksGapPositive = 0;
    double Qcos[3] = {0};
    double Qsin[3] = {0};
    TComplex Qvector[3] = {0};
    double QcosGapNegative[3] = {0};
    double QcosGapPositive[3] = {0};
    double QsinGapNegative[3] = {0};
    double QsinGapPositive[3] = {0};
    TComplex QvectorGapNegative[3] = {0};
    TComplex QvectorGapPositive[3] = {0};

    int NtracksDiff[sizePt] = {0};
    int NtracksDiffGapNegative[sizePt] = {0};
    int NtracksDiffGapPositive[sizePt] = {0};
    double pcos[3][sizePt] = {0};
    double psin[3][sizePt] = {0};
    double pcosGapNegative[3][sizePt] = {0};
    double pcosGapPositive[3][sizePt] = {0};
    double psinGapNegative[3][sizePt] = {0};
    double psinGapPositive[3][sizePt] = {0};
    TComplex pvector[3] = {0};
    TComplex pvectorGapNegative[3] = {0};
    TComplex pvectorGapPositive[3] = {0};

    for(auto& track : tracks) {
      registry.fill(HIST("hPtHadrons"), track.pt());
      registry.fill(HIST("hEtaHadrons"), track.eta());

      if(track.pt() < 0.2) continue;
      if(track.pt() > 5.0) continue;
      if(TMath::Abs(track.eta()) > 2.0) continue;

      registry.fill(HIST("hPtHadronsSelected"), track.pt());
      registry.fill(HIST("hEtaHadronsSelected"), track.eta());
      registry.fill(HIST("hPhiHadronsSelected"), track.phi());
      registry.fill(HIST("hPhiEtaHadronsSelected"), track.phi(), track.eta());
      registry.fill(HIST("hDCAzHadronsSelected"), track.dcaZ());

      Ntracks++;

      int indexPt = getPtIndex(track.pt(), vbins);
      if(indexPt == -1) continue;

      for(int iharm=0; iharm<3; iharm++){
        Qcos[iharm] += TMath::Cos(iharm*track.phi());
        Qsin[iharm] += TMath::Sin(iharm*track.phi());
        pcos[iharm][indexPt] += TMath::Cos(iharm*track.phi());
        psin[iharm][indexPt] += TMath::Sin(iharm*track.phi());
      }

      NtracksDiff[indexPt]++;

      if(track.eta() < 0.) {
        for(int iharm=0; iharm<3;iharm++){
          QcosGapNegative[iharm] += TMath::Cos(iharm*track.phi());
          QsinGapNegative[iharm] += TMath::Sin(iharm*track.phi());
          pcosGapNegative[iharm][indexPt] += TMath::Cos(iharm*track.phi());
          psinGapNegative[iharm][indexPt] += TMath::Sin(iharm*track.phi());
        }
        NtracksGapNegative++;
        NtracksDiffGapNegative[indexPt]++;
      }

      if(track.eta() > 0.) {
        for(int iharm=0; iharm<3; iharm++){
          QcosGapPositive[iharm] += TMath::Cos(iharm*track.phi());
          QsinGapPositive[iharm] += TMath::Sin(iharm*track.phi());
          pcosGapPositive[iharm][indexPt] += TMath::Cos(iharm*track.phi());
          psinGapPositive[iharm][indexPt] += TMath::Sin(iharm*track.phi());
        }
        NtracksGapPositive++;
        NtracksDiffGapPositive[indexPt]++;
      }
    } // charged hadrons

    registry.fill(HIST("hNtracks"), Ntracks);

    for(int iharm=0; iharm<3; iharm++) {
      Qvector[iharm] = TComplex(Qcos[iharm], Qsin[iharm]);
      QvectorGapNegative[iharm] = TComplex(QcosGapNegative[iharm], QsinGapNegative[iharm]);
      QvectorGapPositive[iharm] = TComplex(QcosGapPositive[iharm], QsinGapPositive[iharm]);
    }
 
    double numerator = getCorrelation(Qvector, 2, 2);
    double denominator = getCorrelation(Qvector, 0, 0);
    double corr = numerator/denominator;

    if(Ntracks > 1 && denominator != 0) {
      registry.fill(HIST("hcorrRef"), Ntracks, corr, denominator);
    }

    double numeratorGap = getCorrelationGap(QvectorGapNegative, QvectorGapPositive, 2, 2);
    double denominatorGap = getCorrelationGap(QvectorGapNegative, QvectorGapPositive, 0, 0);
    double corrGap = numeratorGap/denominatorGap;

    if(NtracksGapNegative > 0 && NtracksGapPositive > 0 && denominatorGap != 0) {
      registry.fill(HIST("hcorrRefGap"), Ntracks, corrGap, denominatorGap);
    }

    for(int ipt = 0; ipt < sizePt; ipt ++) {
      for(int iharm=0; iharm<3; iharm++){
        pvector[iharm] = TComplex(pcos[iharm][ipt], psin[iharm][ipt]);
        pvectorGapNegative[iharm] = TComplex(pcosGapNegative[iharm][ipt], psinGapNegative[iharm][ipt]);
        pvectorGapPositive[iharm] = TComplex(pcosGapPositive[iharm][ipt], psinGapPositive[iharm][ipt]);
      }

      //  no eta gap
      double numeratorDiff = getCorrelationDiff(pvector, Qvector, 2, 2);
      double denominatorDiff = getCorrelationDiff(pvector, Qvector, 0, 0);
      double corrDiff = numeratorDiff/denominatorDiff;

      if(NtracksDiff[ipt] > 0 && Ntracks > 0 && denominatorDiff != 0) {
        registry.fill(HIST("hcorrDiff"), (vbins[ipt]+vbins[ipt+1])/2.0, corrDiff, denominatorDiff);
      }

      // case 1 (hadrons from negative eta)
      double numeratorDiffGappNeg = getCorrelationGap(pvectorGapNegative, QvectorGapPositive, 2, 2);
      double denominatorDiffGappNeg = getCorrelationGap(pvectorGapNegative, QvectorGapPositive, 0, 0);
      double corrDiffGappNeg = numeratorDiffGappNeg/denominatorDiffGappNeg;

      if(NtracksDiffGapNegative[ipt] > 0 && NtracksGapPositive > 0 && denominatorDiffGappNeg != 0) {
        registry.fill(HIST("hcorrDiffGapNeg"), (vbins[ipt]+vbins[ipt+1])/2.0, corrDiffGappNeg, denominatorDiffGappNeg);
      }

      // case 2 (hadrons from positive eta)
      double numeratorDiffGappPos = getCorrelationGap(QvectorGapNegative, pvectorGapPositive, 2, 2);
      double denominatorDiffGappPos = getCorrelationGap(QvectorGapNegative, pvectorGapPositive, 0, 0);
      double corrDiffGappPos = numeratorDiffGappPos/denominatorDiffGappPos;

      if(NtracksDiffGapPositive[ipt] > 0 && NtracksGapNegative > 0 && denominatorDiffGappPos != 0) {
        registry.fill(HIST("hcorrDiffGapPos"), (vbins[ipt]+vbins[ipt+1])/2.0, corrDiffGappPos, denominatorDiffGappPos);
      }
    } // pT

    //  Lambda candidates
    int NtracksLcNegative[sizePt][sizeMinv] = {0};
    int NtracksLcPositive[sizePt][sizeMinv] = {0};
    double pcosLcNegative[3][sizePt][sizeMinv] = {0};
    double pcosLcPositive[3][sizePt][sizeMinv] = {0};
    double psinLcNegative[3][sizePt][sizeMinv] = {0};
    double psinLcPositive[3][sizePt][sizeMinv] = {0};
    TComplex pvectorLc[3] = {0};
    TComplex pvectorLcNegative[3] = {0};
    TComplex pvectorLcPositive[3] = {0};

    for (auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << DecayType::LcToPKPi)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YLc(candidate)) > cutYCandMax) {
        continue;
      }
      if (candidate.isSelLcpKpi() >= d_selectionFlagLc) {
        registry.fill(HIST("hmass"), InvMassLcpKpi(candidate), candidate.pt());
      }
      if (candidate.isSelLcpiKp() >= d_selectionFlagLc) {
        registry.fill(HIST("hmass2"), InvMassLcpiKp(candidate), candidate.pt());
      }
      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPhiEtaCand"), candidate.phi(), candidate.eta());
     
      int indexPt = getPtIndex(candidate.pt(), vbins);
      int indexMinv = getMinvIndex(InvMassLcpKpi(candidate), minvarray);

      if(indexPt == -1) continue;
      if(indexMinv == -1) continue;

      if(candidate.eta() < 0.) {
        for(int iharm=0; iharm<3; iharm++){
          pcosLcNegative[iharm][indexPt][indexMinv] += TMath::Cos(iharm*candidate.phi());
          psinLcNegative[iharm][indexPt][indexMinv] += TMath::Sin(iharm*candidate.phi());
        }
        NtracksLcNegative[indexPt][indexMinv]++;
        registry.fill(HIST("hmassEtaNeg"), InvMassLcpKpi(candidate), candidate.pt());
      }

      if(candidate.eta() > 0.) {
        for(int iharm=0; iharm<3; iharm++){
          pcosLcPositive[iharm][indexPt][indexMinv] += TMath::Cos(iharm*candidate.phi());
          psinLcPositive[iharm][indexPt][indexMinv] += TMath::Sin(iharm*candidate.phi());
        }
        NtracksLcPositive[indexPt][indexMinv]++;
        registry.fill(HIST("hmassEtaPos"), InvMassLcpKpi(candidate), candidate.pt());
      }
    } // Lc candidates

    for(int ipt = 0; ipt < sizePt; ipt ++) {
      for(int iminv = 0; iminv < sizeMinv; iminv ++) {
        double binPt = vbins[ipt] + (vbins[ipt+1]-vbins[ipt])/2.0;
        double binMinv = minvarray[iminv] + (minvarray[iminv+1]-minvarray[iminv])/2.0;

        for(int iharm=0; iharm<3; iharm++){
          pvectorLcNegative[iharm] = TComplex(pcosLcNegative[iharm][ipt][iminv], psinLcNegative[iharm][ipt][iminv]);
          pvectorLcPositive[iharm] = TComplex(pcosLcPositive[iharm][ipt][iminv], psinLcPositive[iharm][ipt][iminv]);
        }

        // case 1 (Lambda from negative eta)
        double numeratorLcNeg = getCorrelationGap(pvectorLcNegative, QvectorGapPositive, 2, 2);
        double denominatorLcNeg = getCorrelationGap(pvectorLcNegative, QvectorGapPositive, 0, 0);
        double corrLcNeg = numeratorLcNeg/denominatorLcNeg;

        if(NtracksLcNegative[ipt][iminv] > 0 && NtracksGapPositive > 0 && denominatorLcNeg != 0) {
          registry.fill(HIST("hcorrLcDiffGapNeg"), binMinv, binPt, corrLcNeg, denominatorLcNeg);
        }

        // case 2 (Lambda from positive eta)
        double numeratorLcPos = getCorrelationGap(QvectorGapNegative, pvectorLcPositive, 2, 2);
        double denominatorLcPos = getCorrelationGap(QvectorGapNegative, pvectorLcPositive, 0, 0);
        double corrLcPos = numeratorLcPos/denominatorLcPos;

        if(NtracksLcPositive[ipt][iminv] > 0 && NtracksGapNegative > 0 && denominatorLcPos != 0) {
          registry.fill(HIST("hcorrLcDiffGapPos"), binMinv, binPt, corrLcPos, denominatorLcPos);
        }
      } // minv
    } // pT

  } // process()

  int getPtIndex(double candidatePt, std::vector<double> vbins) {
    int indexPt = -1;
    for(int ipt=0; ipt<10; ipt++) {
      if(candidatePt > vbins[ipt] && candidatePt <= vbins[ipt+1]) {
        indexPt = ipt;
      }
    }
    return indexPt;
  }

  int getMinvIndex(double candidateMinv, double minvarray[]) {
    int indexMinv = -1;
    for(int iminv = 0; iminv < 30; iminv++) {
      if(candidateMinv > minvarray[iminv] && candidateMinv <= minvarray[iminv+1]) {
        indexMinv = iminv;
      }
    }
    return indexMinv;
  }

  double getCorrelation(TComplex vector[], int n1, int n2)
  {
    TComplex correlation = TComplex(0,0);
    correlation = vector[n1]*TComplex::Conjugate(vector[n2]) - vector[0];
    return correlation.Re();
  }
  
  double getCorrelationDiff(TComplex pvector[], TComplex Qvector[], int n1, int n2)
  {
    TComplex correlation = TComplex(0,0);
    correlation = pvector[n1]*TComplex::Conjugate(Qvector[n2]) - pvector[0];
    return correlation.Re();
  }
  
  double getCorrelationGap(TComplex vector1[], TComplex vector2[], int n1, int n2)
  {
    TComplex correlation = TComplex(0,0);
    correlation = vector1[n1]*TComplex::Conjugate(vector2[n2]);
    return correlation.Re();
  }
  
}; // end of struct of Lc flow task

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<TaskLcFlow>(cfgc, TaskName{"hf-task-lc-flow"})};
  return workflow;
}
