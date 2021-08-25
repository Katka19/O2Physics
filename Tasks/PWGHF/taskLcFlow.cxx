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
#include "AnalysisDataModel/HFSecondaryVertex.h"
#include "AnalysisDataModel/HFCandidateSelectionTables.h"
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
      {"hPhiHadronsSelected", "charged hadrons after cuts; #varphi; entries", {HistType::kTH1F, {{180, 0., 2*TMath::Pi()}}}},
      {"hPtCand", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c}); entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"Ntracks", "charged hadrons; Number of tracks; entries", {HistType::kTH1F, {{100, 0., 100}}}},
      {"hcorrelations", "charged hadron <<2>>; Number of tracks; <<2>>", {HistType::kTProfile, {{100, 0., 100}}}},
      {"hcorrelationsGap", "charged hadron <<2>> |#Delta #eta| > 0.4; Number of tracks; <<2>> |#Delta #eta| > 0.4", {HistType::kTProfile, {{100, 0., 100}}}}
      }};

  Configurable<int> d_selectionFlagLc{"d_selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_lc_topkpi::pTBins_v}, "pT bin limits"};

  Filter filterSelectCandidates = (aod::hf_selcandidate_lc::isSelLcpKpi >= d_selectionFlagLc || aod::hf_selcandidate_lc::isSelLcpiKp >= d_selectionFlagLc);

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)bins;
    registry.add("hmass", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 1.6, 3.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  /// aod::BigTracks is not soa::Filtered, should be added when filters are added
  void process(aod::Collision const& collisions, aod::BigTracks const& tracks, soa::Filtered<soa::Join<aod::HfCandProng3, aod::HFSelLcCandidate>> const& candidates)
  {

    int Ntracks = 0;
    int NtracksGapNegative = 0;
    int NtracksGapPositive = 0;
    double Qcos0 = 0;
    double Qsin0 = 0;
    double Qcos2 = 0;
    double Qsin2 = 0;
    TComplex Qvector2;
    TComplex Qvector0;
    double Qcos0GapNegative = 0;
    double Qcos0GapPositive = 0;
    double Qsin0GapNegative = 0;
    double Qsin0GapPositive = 0;
    double Qcos2GapNegative = 0;
    double Qcos2GapPositive = 0;
    double Qsin2GapNegative = 0;
    double Qsin2GapPositive = 0;
    TComplex Qvector2GapNegative;
    TComplex Qvector2GapPositive;
    TComplex Qvector0GapNegative;
    TComplex Qvector0GapPositive;

    for(auto& track : tracks) {
      registry.fill(HIST("hPtHadrons"), track.pt());
      registry.fill(HIST("hEtaHadrons"), track.eta());

      if(track.pt() < 0.2) continue;
      if(track.pt() > 3.0) continue;
      if(TMath::Abs(track.eta()) > 0.8) continue;

      registry.fill(HIST("hPtHadronsSelected"), track.pt());
      registry.fill(HIST("hEtaHadronsSelected"), track.eta());
      registry.fill(HIST("hPhiHadronsSelected"), track.phi());

      Ntracks++;

      Qcos0 += TMath::Cos(0*track.phi());
      Qsin0 += TMath::Sin(0*track.phi());
      Qcos2 += TMath::Cos(2*track.phi());
      Qsin2 += TMath::Sin(2*track.phi());

      if(track.eta() < -0.2) {
        Qcos0GapNegative += TMath::Cos(0*track.phi());
        Qsin0GapNegative += TMath::Sin(0*track.phi());
        Qcos2GapNegative += TMath::Cos(2*track.phi());
        Qsin2GapNegative += TMath::Sin(2*track.phi());
        NtracksGapNegative++;
      }

      if(track.eta() > 0.2) {
        Qcos0GapPositive += TMath::Cos(0*track.phi());
        Qsin0GapPositive += TMath::Sin(0*track.phi());
        Qcos2GapPositive += TMath::Cos(2*track.phi());
        Qsin2GapPositive += TMath::Sin(2*track.phi());
        NtracksGapPositive++;
      }
    } // charged hadrons

    registry.fill(HIST("Ntracks"), Ntracks);

    Qvector0 = TComplex(Qcos0, Qsin0);
    Qvector2 = TComplex(Qcos2, Qsin2);
    Qvector0GapNegative = TComplex(Qcos0GapNegative, Qsin0GapNegative);
    Qvector0GapPositive = TComplex(Qcos0GapPositive, Qsin0GapPositive);
    Qvector2GapNegative = TComplex(Qcos2GapNegative, Qsin2GapNegative);
    Qvector2GapPositive = TComplex(Qcos2GapPositive, Qsin2GapPositive);

    TComplex numerator = Qvector2*TComplex::Conjugate(Qvector2) - Qvector0;
    TComplex denominator = Qvector0*TComplex::Conjugate(Qvector0) - Qvector0;
    double num = numerator.Re();
    double den = denominator.Re();
    double corr = num/den;

    if(Ntracks > 1) {
      registry.fill(HIST("hcorrelations"), Ntracks, corr, den);
    }

    TComplex numeratorGap = Qvector2GapNegative*TComplex::Conjugate(Qvector2GapPositive);
    TComplex denominatorGap = Qvector0GapNegative*TComplex::Conjugate(Qvector0GapPositive);
    double numGap = numeratorGap.Re();
    double denGap = denominatorGap.Re();
    double corrGap = numGap/denGap;

    if(NtracksGapNegative > 0 && NtracksGapPositive > 0) {
      registry.fill(HIST("hcorrelationsGap"), Ntracks, corrGap, denGap);
    }

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
        registry.fill(HIST("hmass"), InvMassLcpiKp(candidate), candidate.pt());
      }
      registry.fill(HIST("hPtCand"), candidate.pt());
    } // Lc candidates
  } // process()
}; // end of struct of Lc flow task

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<TaskLcFlow>(cfgc, TaskName{"hf-task-lc-flow"})};
  return workflow;
}
