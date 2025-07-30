#pragma once

#include <functional>
#include <cmath>
#include <string>
#include "duneanaobj/StandardRecord/StandardRecord.h"
#include <vector>
#include <utility>
#include "TVector3.h"

bool kIsNC(const caf::StandardRecord* sr, int idx);

bool kIsCC(const caf::StandardRecord* sr,int idx);

bool kIsNumu(const caf::StandardRecord* sr,int idx);

bool kIsAntiNumu(const caf::StandardRecord* sr,int idx);

bool kIsVtxCont(const caf::StandardRecord* sr,int idx);

bool kIsVtxFV(const caf::SRVector3D& vertex);

bool kIsArgon(const caf::StandardRecord* sr,int idx);

bool AngleCrossMx2(const caf::StandardRecord* sr,int idx);

bool EnergyCrossMx2(const caf::StandardRecord* sr,int idx);

bool IsBasicCCNuAr(const caf::StandardRecord* sr, int idx,const char* mode = "RHC" );

bool ContLargeProtons(const caf::StandardRecord* sr,int idx);

bool QELikeSignal(const caf::StandardRecord* sr,int idx, const char* mode = "RHC");

// Declare the track_minerva_match function
bool track_minerva_match(int nixn, int npart, caf::StandardRecord* sr);


std::pair<int, int> getNumRecoSecPrim(const caf::StandardRecord* sr, int idx);
std::pair<bool, bool> has_valid_muon(const caf::StandardRecord* sr, int idx);
bool RecoIsExitingMu(const caf::StandardRecord* sr, int idx);
bool RecoIsBasicCCNuAr(const caf::StandardRecord* sr, int idx);
bool RecoHasOnePrimary(const caf::StandardRecord* sr, int idx);
bool RecoIsQELikeSignal(const caf::StandardRecord* sr, int idx, int prim_limit = 1, int sec_limit = 3);