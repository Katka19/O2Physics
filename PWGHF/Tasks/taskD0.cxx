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

/// \file taskD0.cxx
/// \brief D0 analysis task
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::analysis::hf_cuts_d0_topik;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, true, {"Fill MC histograms."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// D0 analysis task
struct TaskD0 {
  HistogramRegistry registry{
    "registry",
    {{"hptcand", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hptprong0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hptprong1", "2-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}}}};

  Configurable<int> d_selectionFlagD0{"d_selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> d_selectionFlagD0bar{"d_selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> bins{"pTBins", std::vector<double>{hf_cuts_d0_topik::pTBins_v}, "pT bin limits"};

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)bins;
    registry.add("hmass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hdeclength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hdeclengthxy", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0d0", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCTS", "2-prong candidates;cos #it{#theta}* (D^{0});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCt", "2-prong candidates;proper lifetime (D^{0}) * #it{c} (cm);entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hselectionstatus", "2-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr", "2-prong candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "2-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenXYErr", "2-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  Filter filterSelectCandidates = (aod::hf_selcandidate_d0::isSelD0 >= d_selectionFlagD0 || aod::hf_selcandidate_d0::isSelD0bar >= d_selectionFlagD0bar);

  void process(soa::Filtered<soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate>> const& candidates)
  {
    for (auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << DecayType::D0ToPiK)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YD0(candidate)) > cutYCandMax) {
        continue;
      }

      if (candidate.isSelD0() >= d_selectionFlagD0) {
        registry.fill(HIST("hmass"), InvMassD0(candidate), candidate.pt());
      }
      if (candidate.isSelD0bar() >= d_selectionFlagD0bar) {
        registry.fill(HIST("hmass"), InvMassD0bar(candidate), candidate.pt());
      }

      registry.fill(HIST("hptcand"), candidate.pt());
      registry.fill(HIST("hptprong0"), candidate.ptProng0());
      registry.fill(HIST("hptprong1"), candidate.ptProng1());
      registry.fill(HIST("hdeclength"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("hdeclengthxy"), candidate.decayLengthXY(), candidate.pt());
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("hd0d0"), candidate.impactParameterProduct(), candidate.pt());
      registry.fill(HIST("hCTS"), CosThetaStarD0(candidate), candidate.pt());
      registry.fill(HIST("hCt"), CtD0(candidate), candidate.pt());
      registry.fill(HIST("hCPA"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("hEta"), candidate.eta(), candidate.pt());
      registry.fill(HIST("hselectionstatus"), candidate.isSelD0() + (candidate.isSelD0bar() * 2), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), candidate.pt());
    }
  }
};

/// Fills MC histograms.
struct TaskD0MC {
  HistogramRegistry registry{
    "registry",
    {{"hPtRecSig", "2-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtRecSigPrompt", "2-prong candidates (matched, prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtRecSigNonPrompt", "2-prong candidates (matched, non-prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtRecBg", "2-prong candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtGen", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtGenPrompt", "MC particles (matched, prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtGenNonPrompt", "MC particles (matched, non-prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtGenSig", "2-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hCPARecSig", "2-prong candidates (matched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hCPARecBg", "2-prong candidates (unmatched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hEtaRecSig", "2-prong candidates (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -5., 5.}}}},
     {"hEtaRecBg", "2-prong candidates (unmatched);#it{#eta};entries", {HistType::kTH1F, {{100, -5., 5.}}}},
     {"hEtaGen", "MC particles (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -5., 5.}}}},
     {"hptProng0Sig", "prong0 pt (matched); #it{y}", {HistType::kTH2F, {{300, 0., 30.}, {10, -5., 5.}}}},
     {"hptProng1Sig", "prong1 pt (matched); #it{y}", {HistType::kTH2F, {{300, 0., 30.}, {10, -5., 5.}}}},
     {"hdeclengthSig", "2-prong candidates (matched);decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hdeclengthxySig", "2-prong candidates (matched);decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hnormaliseddeclengthSig", "2-prong candidates (matched);normalised decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hnormaliseddeclengthxySig", "2-prong candidates (matched);normalised decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hd0Prong0Sig", "2-prong candidates (matched);prong 0 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0Prong1Sig", "2-prong candidates (matched);prong 1 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0d0Sig", "2-prong candidates (matched);product of DCAxy to prim. vertex (cm^{2}); #it{y}", {HistType::kTH2F, {{500, -1., 1.}, {10, -5., 5.}}}},
     {"hCTSSig", "2-prong candidates (matched);cos #it{#theta}* (D^{0}); #it{y}", {HistType::kTH2F, {{110, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCtSig", "2-prong candidates (matched);proper lifetime (D^{0}) * #it{c} (cm); #it{y}", {HistType::kTH2F, {{120, -20., 100.}, {10, -5., 5.}}}},
     {"hCPASig", "2-prong candidates (matched);cosine of pointing angle; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCPAxySig", "2-prong candidates (matched);cosine of pointing angle xy; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hptProng0Bkg", "prong0 pt (matched); #it{y}", {HistType::kTH2F, {{300, 0., 30.}, {10, -5., 5.}}}},
     {"hptProng1Bkg", "prong1 pt (matched); #it{y}", {HistType::kTH2F, {{300, 0., 30.}, {10, -5., 5.}}}},
     {"hdeclengthBkg", "2-prong candidates (checked);decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hdeclengthxyBkg", "2-prong candidates (checked);decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hnormaliseddeclengthBkg", "2-prong candidates (checked);normalised decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hnormaliseddeclengthxyBkg", "2-prong candidates (checked);normalised decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hd0Prong0Bkg", "2-prong candidates (checked);prong 0 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0Prong1Bkg", "2-prong candidates (checked);prong 1 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0d0Bkg", "2-prong candidates (checked);product of DCAxy to prim. vertex (cm^{2}); #it{y}", {HistType::kTH2F, {{500, -1., 1.}, {10, -5., 5.}}}},
     {"hCTSBkg", "2-prong candidates (checked);cos #it{#theta}* (D^{0}); #it{y}", {HistType::kTH2F, {{110, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCtBkg", "2-prong candidates (checked);proper lifetime (D^{0}) * #it{c} (cm); #it{y}", {HistType::kTH2F, {{120, -20., 100.}, {10, -5., 5.}}}},
     {"hCPABkg", "2-prong candidates (checked);cosine of pointing angle; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCPAxyBkg", "2-prong candidates (checked);cosine of pointing angle xy; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hPtvsYRecSig_RecoPID", "2-prong candidates (RecoPID - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{100, 0., 10.}, {100, -5., 5.}}}},
     {"hPtvsYRecSigPrompt_RecoPID", "2-prong candidates (RecoPID - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{100, 0., 10.}, {100, -5., 5.}}}},
     {"hPtvsYRecSigNonPrompt_RecoPID", "2-prong candidates (RecoPID - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{100, 0., 10.}, {100, -5., 5.}}}},
     {"hPtvsYRecSig_RecoCand", "2-prong candidates (RecoCand - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{100, 0., 10.}, {100, -5., 5.}}}},
     {"hPtvsYRecSigPrompt_RecoCand", "2-prong candidates (RecoCand - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{100, 0., 10.}, {100, -5., 5.}}}},
     {"hPtvsYRecSigNonPrompt_RecoCand", "2-prong candidates (RecoCand - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{100, 0., 10.}, {100, -5., 5.}}}},
     {"hPtvsYRecSig_RecoTopol", "2-prong candidates (RecoTopol - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{100, 0., 10.}, {100, -5., 5.}}}},
     {"hPtvsYRecSigPrompt_RecoTopol", "2-prong candidates (RecoTopol - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{100, 0., 10.}, {100, -5., 5.}}}},
     {"hPtvsYRecSigNonPrompt_RecoTopol", "2-prong candidates (RecoTopol - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{100, 0., 10.}, {100, -5., 5.}}}},
     {"hPtvsYRecSig_RecoHFFlag", "2-prong candidates (RecoHFFlag - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{100, 0., 10.}, {100, -5., 5.}}}},
     {"hPtvsYRecSigPrompt_RecoHFFlag", "2-prong candidates (RecoHFFlag - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{100, 0., 10.}, {100, -5., 5.}}}},
     {"hPtvsYRecSigNonPrompt_RecoHFFlag", "2-prong candidates (RecoHFFlag - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{100, 0., 10.}, {100, -5., 5.}}}},
     {"hPtvsYGen", "2-prong candidates (matched);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{100, 0., 10.}, {100, -5., 5.}}}},
     {"hPtvsYGenPrompt", "2-prong candidates (matched, prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{100, 0., 10.}, {100, -5., 5.}}}},
     {"hPtvsYGenNonPrompt", "2-prong candidates (matched, non-prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{100, 0., 10.}, {100, -5., 5.}}}},
     {"hMassSigD0", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassBkgD0", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassReflBkgD0", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassSigBkgD0", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassSigD0bar", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassBkgD0bar", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassReflBkgD0bar", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassSigBkgD0bar", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}}}};

  Configurable<int> d_selectionFlagD0{"d_selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> d_selectionFlagD0bar{"d_selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<int> d_selectionHFFlag{"d_selectionHFFlag", 1, "Selection Flag for HF flagged candidates"};
  Configurable<int> d_selectionTopol{"d_selectionTopol", 1, "Selection Flag for topologically selected candidates"};
  Configurable<int> d_selectionCand{"d_selectionCand", 1, "Selection Flag for conj. topol. selected candidates"};
  Configurable<int> d_selectionPID{"d_selectionPID", 1, "Selection Flag for reco PID candidates"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};

  Filter filterSelectCandidates = (aod::hf_selcandidate_d0::isRecoHFFlag >= d_selectionHFFlag);
  //Filter filterSelectCandidates = (aod::hf_selcandidate_d0::isSelD0 >= d_selectionFlagD0 || aod::hf_selcandidate_d0::isSelD0bar >= d_selectionFlagD0bar);

  void process(soa::Filtered<soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate, aod::HfCandProng2MCRec>> const& candidates,
               soa::Join<aod::McParticles, aod::HfCandProng2MCGen> const& particlesMC, aod::BigTracksMC const& tracks)
  {
    // MC rec.
    //Printf("MC Candidates: %d", candidates.size());
    for (auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << DecayType::D0ToPiK)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YD0(candidate)) > cutYCandMax) {
        continue;
      }
      if (std::abs(candidate.flagMCMatchRec()) == 1 << DecayType::D0ToPiK) {
        // Get the corresponding MC particle.
        auto indexMother = RecoDecay::getMother(particlesMC, candidate.index0_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandProng2MCGen>>(), pdg::Code::kD0, true);
        auto particleMother = particlesMC.iteratorAt(indexMother);
        registry.fill(HIST("hPtGenSig"), particleMother.pt()); // gen. level pT
        auto ptRec = candidate.pt();
        auto yRec = YD0(candidate);
        registry.fill(HIST("hPtRecSig"), ptRec); // rec. level pT
        if (candidate.isRecoHFFlag() >= d_selectionHFFlag) {
          registry.fill(HIST("hPtvsYRecSig_RecoHFFlag"), ptRec, yRec);
        }
        if (candidate.isRecoTopol() >= d_selectionTopol) {
          registry.fill(HIST("hPtvsYRecSig_RecoTopol"), ptRec, yRec);
        }
        if (candidate.isRecoCand() >= d_selectionCand) {
          registry.fill(HIST("hPtvsYRecSig_RecoCand"), ptRec, yRec);
        }
        if (candidate.isRecoPID() >= d_selectionPID) {
          registry.fill(HIST("hPtvsYRecSig_RecoPID"), ptRec, yRec);
        }

        if (candidate.originMCRec() == OriginType::Prompt) {
          registry.fill(HIST("hPtRecSigPrompt"), ptRec); // rec. level pT, prompt
          if (candidate.isRecoHFFlag() >= d_selectionHFFlag) {
            registry.fill(HIST("hPtvsYRecSigPrompt_RecoHFFlag"), ptRec, yRec);
          }
          if (candidate.isRecoTopol() >= d_selectionTopol) {
            registry.fill(HIST("hPtvsYRecSigPrompt_RecoTopol"), ptRec, yRec);
          }
          if (candidate.isRecoCand() >= d_selectionCand) {
            registry.fill(HIST("hPtvsYRecSigPrompt_RecoCand"), ptRec, yRec);
          }
          if (candidate.isRecoPID() >= d_selectionPID) {
            registry.fill(HIST("hPtvsYRecSigPrompt_RecoPID"), ptRec, yRec);
          }
        } else {
          registry.fill(HIST("hPtRecSigNonPrompt"), ptRec); // rec. level pT, non-prompt
          if (candidate.isRecoHFFlag() >= d_selectionHFFlag) {
            registry.fill(HIST("hPtvsYRecSigNonPrompt_RecoHFFlag"), ptRec, yRec);
          }
          if (candidate.isRecoTopol() >= d_selectionTopol) {
            registry.fill(HIST("hPtvsYRecSigNonPrompt_RecoTopol"), ptRec, yRec);
          }
          if (candidate.isRecoCand() >= d_selectionCand) {
            registry.fill(HIST("hPtvsYRecSigNonPrompt_RecoCand"), ptRec, yRec);
          }
          if (candidate.isRecoPID() >= d_selectionPID) {
            registry.fill(HIST("hPtvsYRecSigNonPrompt_RecoPID"), ptRec, yRec);
          }
        }
        registry.fill(HIST("hCPARecSig"), candidate.cpa());
        registry.fill(HIST("hEtaRecSig"), candidate.eta());
      } else {
        registry.fill(HIST("hPtRecBg"), candidate.pt());
        registry.fill(HIST("hCPARecBg"), candidate.cpa());
        registry.fill(HIST("hEtaRecBg"), candidate.eta());
      }
      auto massD0 = InvMassD0(candidate);
      auto massD0bar = InvMassD0bar(candidate);
      auto ptCandidate = candidate.pt();
      auto ptProng0 = candidate.ptProng0();
      auto ptProng1 = candidate.ptProng1();
      auto rapidityCandidate = YD0(candidate);
      auto declengthCandidate = candidate.decayLength();
      auto declengthxyCandidate = candidate.decayLengthXY();
      auto normaliseddeclengthCandidate = candidate.decayLengthNormalised();
      auto normaliseddeclengthxyCandidate = candidate.decayLengthXYNormalised();
      auto d0Prong0 = candidate.impactParameter0();
      auto d0Prong1 = candidate.impactParameter1();
      auto d0d0Candidate = candidate.impactParameterProduct();
      auto ctsCandidate = CosThetaStarD0(candidate);
      auto ctCandidate = CtD0(candidate);
      auto cpaCandidate = candidate.cpa();
      auto cpaxyCandidate = candidate.cpaXY();
      if (candidate.isSelD0() >= d_selectionFlagD0) {
        registry.fill(HIST("hMassSigBkgD0"), massD0, ptCandidate, rapidityCandidate);
        if (candidate.flagMCMatchRec() == (1 << DecayType::D0ToPiK)) {
          registry.fill(HIST("hptProng0Sig"), ptProng0, rapidityCandidate);
          registry.fill(HIST("hptProng1Sig"), ptProng1, rapidityCandidate);
          registry.fill(HIST("hdeclengthSig"), declengthCandidate, rapidityCandidate);
          registry.fill(HIST("hdeclengthxySig"), declengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("hnormaliseddeclengthSig"), normaliseddeclengthCandidate, rapidityCandidate);
          registry.fill(HIST("hnormaliseddeclengthxySig"), normaliseddeclengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("hd0Prong0Sig"), d0Prong0, rapidityCandidate);
          registry.fill(HIST("hd0Prong1Sig"), d0Prong1, rapidityCandidate);
          registry.fill(HIST("hd0d0Sig"), d0d0Candidate, rapidityCandidate);
          registry.fill(HIST("hCTSSig"), ctsCandidate, rapidityCandidate);
          registry.fill(HIST("hCtSig"), ctCandidate, rapidityCandidate);
          registry.fill(HIST("hCPASig"), cpaCandidate, rapidityCandidate);
          registry.fill(HIST("hCPAxySig"), cpaxyCandidate, rapidityCandidate);
          registry.fill(HIST("hMassSigD0"), massD0, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hptProng0Bkg"), ptProng0, rapidityCandidate);
          registry.fill(HIST("hptProng1Bkg"), ptProng1, rapidityCandidate);
          registry.fill(HIST("hdeclengthBkg"), declengthCandidate, rapidityCandidate);
          registry.fill(HIST("hdeclengthxyBkg"), declengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("hnormaliseddeclengthBkg"), normaliseddeclengthCandidate, rapidityCandidate);
          registry.fill(HIST("hnormaliseddeclengthxyBkg"), normaliseddeclengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("hd0Prong0Bkg"), d0Prong0, rapidityCandidate);
          registry.fill(HIST("hd0Prong1Bkg"), d0Prong1, rapidityCandidate);
          registry.fill(HIST("hd0d0Bkg"), d0d0Candidate, rapidityCandidate);
          registry.fill(HIST("hCTSBkg"), ctsCandidate, rapidityCandidate);
          registry.fill(HIST("hCtBkg"), ctCandidate, rapidityCandidate);
          registry.fill(HIST("hCPABkg"), cpaCandidate, rapidityCandidate);
          registry.fill(HIST("hCPAxyBkg"), cpaxyCandidate, rapidityCandidate);
          registry.fill(HIST("hMassBkgD0"), massD0, ptCandidate, rapidityCandidate);
          if (candidate.flagMCMatchRec() == -(1 << DecayType::D0ToPiK)) {
            registry.fill(HIST("hMassReflBkgD0"), massD0, ptCandidate, rapidityCandidate);
          }
        }
      }
      if (candidate.isSelD0bar() >= d_selectionFlagD0) {
        registry.fill(HIST("hMassSigBkgD0bar"), massD0bar, ptCandidate, rapidityCandidate);
        if (candidate.flagMCMatchRec() == -(1 << DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassSigD0bar"), massD0bar, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgD0bar"), massD0bar, ptCandidate, rapidityCandidate);
          if (candidate.flagMCMatchRec() == (1 << DecayType::D0ToPiK)) {
            registry.fill(HIST("hMassReflBkgD0bar"), massD0bar, ptCandidate, rapidityCandidate);
          }
        }
      }
    }
    // MC gen.
    //Printf("MC Particles: %d", particlesMC.size());
    for (auto& particle : particlesMC) {
      if (std::abs(particle.flagMCMatchGen()) == 1 << DecayType::D0ToPiK) {
        if (cutYCandMax >= 0. && std::abs(RecoDecay::Y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()))) > cutYCandMax) {
          continue;
        }
        auto ptGen = particle.pt();
        auto yGen = RecoDecay::Y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()));
        registry.fill(HIST("hPtGen"), ptGen);
        registry.fill(HIST("hPtvsYGen"), ptGen, yGen);
        if (particle.originMCGen() == OriginType::Prompt) {
          registry.fill(HIST("hPtGenPrompt"), ptGen);
          registry.fill(HIST("hPtvsYGenPrompt"), ptGen, yGen);
        } else {
          registry.fill(HIST("hPtGenNonPrompt"), ptGen);
          registry.fill(HIST("hPtvsYGenNonPrompt"), ptGen, yGen);
        }
        registry.fill(HIST("hEtaGen"), particle.eta());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<TaskD0>(cfgc, TaskName{"hf-task-d0"})};
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<TaskD0MC>(cfgc, TaskName{"hf-task-d0-mc"}));
  }
  return workflow;
}
