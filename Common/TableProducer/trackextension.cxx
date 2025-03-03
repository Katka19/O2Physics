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

//
// Task performing basic track selection.
//

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

//****************************************************************************************
/**
 * Produce the more complicated derived track quantities needed for track selection.
 * FIXME: we shall run this only if all other selections are passed to avoid
 * FIXME: computing overhead and errors in calculations
 */
//****************************************************************************************
struct TrackExtensionTask {

  Produces<aod::TracksExtended> extendedTrackQuantities;

  void process(aod::FullTracks const& tracks, aod::Collisions const&)
  {
    for (auto& track : tracks) {

      std::array<float, 2> dca{1e10f, 1e10f};
      if (track.has_collision()) {
        if ((track.trackType() == o2::aod::track::TrackTypeEnum::Track) ||
            (track.trackType() == o2::aod::track::TrackTypeEnum::Run2Track && track.itsChi2NCl() != 0.f && track.tpcChi2NCl() != 0.f && std::abs(track.x()) < 10.f)) {
          float magField = 5.0; // in kG (FIXME: get this from CCDB)
          auto trackPar = getTrackPar(track);
          auto const& collision = track.collision();
          trackPar.propagateParamToDCA({collision.posX(), collision.posY(), collision.posZ()}, magField, &dca);
        }
      }
      extendedTrackQuantities(dca[0], dca[1]);

      // TODO: add realtive pt resolution sigma(pt)/pt \approx pt * sigma(1/pt)
      // TODO: add geometrical length / fiducial volume
    }
  }
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<TrackExtensionTask>(cfgc, TaskName{"track-extension"})};
  return workflow;
}
