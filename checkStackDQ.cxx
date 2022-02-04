// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// Executable to check functioning of stack
// Analyses kinematics and track references of a kinematics file

// Usage: root checkStackDQ.cxx++
//          .x checkStackDQ.cxx++

#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCUtils.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/Stack.h"
#include "SimulationDataFormat/TrackReference.h"
#include "Steer/MCKinematicsReader.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include "FairLogger.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include "ITSMFTSimulation/Hit.h"
#include <unordered_map>
#include <string>

////////////////////////////
// Breit-Wigner Function ///
////////////////////////////

//                2              Γ^2 * M^2
//B(m;M,Γ) = N * ---  -------------------------------       ----> Relativistic Breit-Wigner Distribution
//                π    (m^2 - M^2)^2 + m^4(Γ^2 / M^2)

//Γ = Gamma
//π = Pi
Double_t mybwMass(Double_t* x, Double_t* par)
{
  Double_t arg1 = 14.0/22.0; // 2 over pi
  Double_t arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma=par[1]  M=par[2]
  Double_t arg3 = ((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
  Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[1]*par[1])/(par[2]*par[2]));
  return par[0]*arg1*arg2/(arg3 + arg4);
}

int checkStackDQ()
{
  // Define Rapidity, Lepton PDG and Sim Rapidity Type
  Int_t leptonPDG;
  Double_t rapMin, rapMax;
  Bool_t midysim = kFALSE; Bool_t forwardysim = kFALSE;

  // User Guide
  std::cout << "============================================================================================" << '\n';
  std::cout << "          checkStackDQ.cxx is Executable macro to check functioning of stack"                 << '\n';
  std::cout << "          This macro has been developed for the Charmonium Analysis of the PWG-DQ group"      << '\n';
  std::cout << "          It provides analyses kinematics and track references of a kinematics file"          << '\n';
  std::cout << "============================================================================================" << '\n';

  // Enter PDG code daughter for selection
  std::cout << "============================================================================================" << '\n';
  std::cout << "                Selection for daughters (dilepton pairs). You should enter a PDG code."       << '\n';
  std::cout << "                      For electron 11 or For Muon 13 and then press enter."                   << '\n';
  std::cout << "============================================================================================" << '\n';
  std::cin >> leptonPDG;
  assert(leptonPDG == 13 || leptonPDG == 11);

  // Settings for Selections
  if(leptonPDG == 11) { rapMin = -1.5; rapMax = 1.5;  midysim     = kTRUE; }
  if(leptonPDG == 13) { rapMin = -4.3; rapMax = -2.3; forwardysim = kTRUE; }

  /// define histos prompt jpsi / psi2s before cuts
  // TODO : Check the histogram ranges and optimize them.

  // J/Psi
  TH1F *histPtPairJpsiPrompt = new TH1F("PtPairJpsiPromptBeforeCuts",";p_{T}(GeV/c);Entries",100,0.,20.);
  TH1F *histMassPairJpsiPrompt = new TH1F("MassJpsiPromptBeforeCuts",";m_{ll}(GeV/c^{2});Entries",5./0.02,0.,5.);
  TH1F *histYJpsiPrompt = new TH1F("YJpsiPromptBeforeCuts",";y;Entries",100,rapMin,rapMax);
  TH1F *histPtJpsiPrompt = new TH1F("PtJpsiPromptBeforeCuts",";p_{T}(GeV/c);Entries",100,0.,15.);
  TH1F *histEtaJpsiPrompt = new TH1F("EtaJpsiPromptBeforeCuts",";|#eta|;Entries",100,rapMin,rapMax);
  TH1F *histThetaJpsiPrompt = new TH1F("ThetaJpsiPromptBeforeCuts",";|#theta|;Entries",100,rapMin,rapMax);
  TH1F *histPhiJpsiPrompt = new TH1F("PhiJpsiPromptBeforeCuts",";|#phi|;Entries",100,rapMin,rapMax);
  // Psi(2s)
  TH1F *histPtPairPsi2sPrompt = new TH1F("PtPairPsi2sPromptBeforeCuts",";p_{T}(GeV/c);Entries",100,0.,20.);
  TH1F *histMassPairPsi2sPrompt = new TH1F("MassPsi2sPromptBeforeCuts",";m_{ll}(GeV/c^{2});Entries",5./0.02,0.,5.);
  TH1F *histYPsi2sPrompt = new TH1F("YPsi2sPromptBeforeCuts",";y;Entries",100,rapMin,rapMax);
  TH1F *histPtPsi2sPrompt = new TH1F("PtPsi2sPromptBeforeCuts",";p_{T}(GeV/c);Entries",100,0.,15.);
  TH1F *histEtaPsi2sPrompt = new TH1F("EtaPsi2sPromptBeforeCuts",";|#eta|;Entries",100,-20,20);
  TH1F *histThetaPsi2sPrompt = new TH1F("ThetaPsi2sPromptBeforeCuts",";|#theta|;Entries",100,-20,20);
  TH1F *histPhiPsi2sPrompt = new TH1F("PhiPsi2sPromptBeforeCuts",";|#phi|;Entries",100,-20,20);

  /// define histos prompt jpsi / psi2s after cuts

  // J/Psi
  TH1F *histPtPairJpsiPromptAfterCuts = new TH1F("PtPairJpsiPromptAfterCuts",";p_{T}(GeV/c);Entries",100,0.,20.);
  TH1F *histMassPairJpsiPromptAfterCuts = new TH1F("MassJpsiPromptAfterCuts",";m_{ll}(GeV/c^{2});Entries",5./0.02,0.,5.);
  TH1F *histYJpsiPromptAfterCuts = new TH1F("YJpsiPromptAfterCuts",";y;Entries",100,rapMin,rapMax);
  TH1F *histPtJpsiPromptAfterCuts = new TH1F("PtJpsiPromptAfterCuts",";p_{T}(GeV/c);Entries",100,0.,15.);
  TH1F *histEtaJpsiPromptAfterCuts = new TH1F("EtaJpsiPromptAfterCuts",";|#eta|;Entries",100,rapMin,rapMax);
  TH1F *histThetaJpsiPromptAfterCuts = new TH1F("ThetaJpsiPromptAfterCuts",";|#theta|;Entries",100,rapMin,rapMax);
  TH1F *histPhiJpsiPromptAfterCuts = new TH1F("PhiJpsiPromptAfterCuts",";|#phi|;Entries",100,rapMin,rapMax);
  // Psi(2s)
  TH1F *histPtPairPsi2sPromptAfterCuts = new TH1F("PtPairPsi2sPromptAfterCuts",";p_{T}(GeV/c);Entries",100,0.,20.);
  TH1F *histMassPairPsi2sPromptAfterCuts = new TH1F("MassPsi2sPromptAfterCuts",";m_{ll}(GeV/c^{2});Entries",5./0.02,0.,5.);
  TH1F *histYPsi2sPromptAfterCuts = new TH1F("YPsi2sPromptAfterCuts",";y;Entries",100,rapMin,rapMax);
  TH1F *histPtPsi2sPromptAfterCuts = new TH1F("PtPsi2sPromptAfterCuts",";p_{T}(GeV/c);Entries",100,0.,15.);
  TH1F *histEtaPsi2sPromptAfterCuts = new TH1F("EtaPsi2sPromptAfterCuts",";|#eta|;Entries",100,rapMin,rapMax);
  TH1F *histThetaPsi2sPromptAfterCuts = new TH1F("ThetaPsi2sPromptAfterCuts",";|#theta|;Entries",100,rapMin,rapMax);
  TH1F *histPhiPsi2sPromptAfterCuts = new TH1F("PhiPsi2sPromptAfterCuts",";|#phi|;Entries",100,rapMin,rapMax);

  /// define histos non-prompt jpsi / psi2s before cuts

  // J/Psi
  TH1F *histPtPairJpsiNonPrompt = new TH1F("PtPairNonPromptJpsiBeforeCuts",";p_{T}(GeV/c);Entries",100,0.,20.);
  TH1F *histMassPairJpsiNonPrompt = new TH1F("MassNonPromptJpsiBeforeCuts",";m_{ll}(GeV/c^{2});Entries",5./0.02,0.,5.);
  TH1F *histYJpsiNonPrompt = new TH1F("YJpsiNonPromptBeforeCuts",";y;Entries",100,rapMin,rapMax);
  TH1F *histPtJpsiNonPrompt = new TH1F("PtJpsiNonPromptBeforeCuts",";p_{T}(GeV/c);Entries",100,0.,15.);
  TH1F *histEtaJpsiNonPrompt = new TH1F("EtaJpsiNonPromptBeforeCuts",";|#eta|;Entries",100,-20,20);
  TH1F *histThetaJpsiNonPrompt = new TH1F("ThetaJpsiNonPromptBeforeCuts",";|#theta|;Entries",100,-20,20);
  TH1F *histPhiJpsiNonPrompt = new TH1F("PhiJpsiNonPromptBeforeCuts",";|#phi|;Entries",100,20,-20);
  // Psi(2s)
  TH1F *histPtPairPsi2sNonPrompt = new TH1F("PtPairPsi2sNonPromptBeforeCuts",";p_{T}(GeV/c);Entries",100,0.,20.);
  TH1F *histMassPairPsi2sNonPrompt = new TH1F("MassPsi2sNonPromptBeforeCuts",";m_{ll}(GeV/c^{2});Entries",5./0.02,0.,5.);
  TH1F *histYPsi2sNonPrompt = new TH1F("YPsi2sNonPromptBeforeCuts",";y;Entries",100,rapMin,rapMax);
  TH1F *histPtPsi2sNonPrompt = new TH1F("PtPsi2sNonPromptBeforeCuts",";p_{T}(GeV/c);Entries",100,0.,15.);
  TH1F *histEtaPsi2sNonPrompt = new TH1F("EtaPsi2sNonPromptBeforeCuts",";|#eta|;Entries",100,-9,-1);
  TH1F *histThetaPsi2sNonPrompt = new TH1F("ThetaPsi2sNonPromptBeforeCuts",";|#theta|;Entries",100,2.5,4.5);
  TH1F *histPhiPsi2sNonPrompt = new TH1F("PhiPsi2sNonPromptBeforeCuts",";|#phi|;Entries",100,-1,7);

  /// define histos non-prompt jpsi / psi2s after cuts

  // J/Psi
  TH1F *histPtPairJpsiNonPromptAfterCuts = new TH1F("PtPairNonPromptJpsiAfterCuts",";p_{T}(GeV/c);Entries",100,0.,20.);
  TH1F *histMassPairJpsiNonPromptAfterCuts = new TH1F("MassNonPromptJpsiAfterCuts",";m_{ll}(GeV/c^{2});Entries",5./0.02,0.,5.);
  TH1F *histYJpsiNonPromptAfterCuts = new TH1F("YJpsiNonPromptJpsiAfterCuts",";y;Entries",100,rapMin,rapMax);
  TH1F *histPtJpsiNonPromptAfterCuts = new TH1F("PtJpsiNonPromptJpsiAfterCuts",";p_{T}(GeV/c);Entries",100,0.,15.);
  TH1F *histEtaJpsiNonPromptAfterCuts = new TH1F("EtaJpsiNonPromptAfterCuts",";|#eta|;Entries",100,rapMin,rapMax);
  TH1F *histThetaJpsiNonPromptAfterCuts = new TH1F("ThetaJpsiNonPromptAfterCuts",";|#theta|;Entries",100,rapMin,rapMax);
  TH1F *histPhiJpsiNonPromptAfterCuts = new TH1F("PhiJpsiNonPromptAfterCuts",";|#phi|;Entries",100,rapMin,rapMax);
  // Psi(2s)
  TH1F *histPtPairPsi2sNonPromptAfterCuts = new TH1F("PtPairPsi2sNonPromptAfterCuts",";p_{T}(GeV/c);Entries",100,0.,20.);
  TH1F *histMassPairPsi2sNonPromptAfterCuts = new TH1F("MassPsi2sNonPromptAfterCuts",";m_{ll}(GeV/c^{2});Entries",5./0.02,0.,5.);
  TH1F *histYPsi2sNonPromptAfterCuts = new TH1F("YPsi2sNonPromptAfterCuts",";y;Entries",100,rapMin,rapMax);
  TH1F *histPtPsi2sNonPromptAfterCuts = new TH1F("PtPsi2sNonPromptAfterCuts",";p_{T}(GeV/c);Entries",100,0.,15.);
  TH1F *histEtaPsi2sNonPromptAfterCuts = new TH1F("EtaPsi2sNonPromptAfterCuts",";|#eta|;Entries",100,rapMin,rapMax);
  TH1F *histThetaPsi2sNonPromptAfterCuts = new TH1F("ThetaPsi2sNonPromptAfterCuts",";|#theta|;Entries",100,rapMin,rapMax);
  TH1F *histPhiPsi2sNonPromptAfterCuts = new TH1F("PhiPsi2sNonPromptAfterCuts",";|#phi|;Entries",100,rapMin,rapMax);

  /// define histos for Origin before cuts
  TH1D *histIsPrompt = new TH1D("jpsiOrigBeforeCuts","",4,0,4);
  histIsPrompt->GetXaxis()->SetBinLabel(1,"IsJpsiPrompt");
  histIsPrompt->GetXaxis()->SetBinLabel(2,"IsJpsiFromB");
  histIsPrompt->GetXaxis()->SetBinLabel(3,"IsPsi2sPrompt");
  histIsPrompt->GetXaxis()->SetBinLabel(4,"IsPsi2sFromB");

  /// define histos for Origin after cuts
  TH1D *histIsPromptAfterCuts = new TH1D("jpsiOrigAfterCuts","",4,0,4);
  histIsPromptAfterCuts->GetXaxis()->SetBinLabel(1,"IsJpsiPrompt");
  histIsPromptAfterCuts->GetXaxis()->SetBinLabel(2,"IsJpsiFromB");
  histIsPromptAfterCuts->GetXaxis()->SetBinLabel(3,"IsPsi2sPrompt");
  histIsPromptAfterCuts->GetXaxis()->SetBinLabel(4,"IsPsi2sFromB");

  /// define canvas for drawing functions in same histogram

  // Prompt Canvas
  TCanvas *canvasPromptMassCompareAfterCuts = new TCanvas("canvasPromptMassCompareAfterCuts","Compare InvMass J/psi and Psi(2S)",800,600);
  TCanvas *canvasPromptPtPairCompareAfterCuts = new TCanvas("canvasPromptPtPairCompareAfterCuts","Compare Pt Pair J/psi and Psi(2S)",800,600);
  TCanvas *canvasPromptPtCompareAfterCuts = new TCanvas("canvasPromptPtCompareAfterCuts","Compare Pt J/psi and Psi(2S)",800,600);
  TCanvas *canvasPromptyCompareAfterCuts = new TCanvas("canvasPromptyCompareAfterCuts","Compare Y J/psi and Psi(2S)",800,600);

  // Non-Prompt Canvas
  TCanvas *canvasNonPromptMassCompareAfterCuts = new TCanvas("canvasNonPromptMassCompareAfterCuts","Compare InvMass J/psi and Psi(2S)",800,600);
  TCanvas *canvasNonPromptPtPairCompareAfterCuts = new TCanvas("canvasNonPromptPtPairCompareAfterCuts","Compare Pt Pair J/psi and Psi(2S)",800,600);
  TCanvas *canvasNonPromptPtCompareAfterCuts = new TCanvas("canvasNonPromptPtCompareAfterCuts","Compare Pt J/psi and Psi(2S)",800,600);
  TCanvas *canvasNonPromptyCompareAfterCuts = new TCanvas("canvasNonPromptyCompareAfterCuts","Compare Y J/psi and Psi(2S)",800,600);

  /// define histos for Efficiency

  // Prompt
  TH1F *histPtEfficiencyJpsiPrompt = new TH1F("Prompt J/#it{#psi} Acceptance for p_{T}",";p_{T}(GeV/c);Acceptance",100,0.,15.);
  TH1F *histPtPairEfficiencyJpsiPrompt = new TH1F("Prompt J/#it{#psi} Acceptance for p_{T} Pairs",";p_{T}(GeV/c);Acceptance",100,0.,20.);
  TH1F *histPtEfficiencyPsi2sPrompt = new TH1F("Prompt #it{#psi}(2S) Acceptance for p_{T}",";p_{T}(GeV/c);Acceptance",100,0.,15.);
  TH1F *histPtPairEfficiencyPsi2sPrompt = new TH1F("Prompt #it{#psi}(2S) Acceptance for p_{T} Pairs",";p_{T}(GeV/c);Acceptance",100,0.,20.);

  // Non-Prompt
  TH1F *histPtEfficiencyJpsiNonPrompt = new TH1F("NonPrompt J/#it{#psi} Acceptance for p_{T}",";p_{T}(GeV/c);Acceptance",100,0.,15.);
  TH1F *histPtPairEfficiencyJpsiNonPrompt = new TH1F("NonPrompt J/#it{#psi} Acceptance for p_{T} Pairs",";p_{T}(GeV/c);Acceptance",100,0.,20.);
  TH1F *histPtEfficiencyPsi2sNonPrompt = new TH1F("NonPrompt #it{#psi}(2S) Acceptance for p_{T}",";p_{T}(GeV/c);Acceptance",100,0.,15.);
  TH1F *histPtPairEfficiencyPsi2sNonPrompt = new TH1F("NonPrompt #it{#psi}(2S) Acceptance for p_{T} Pairs",";p_{T}(GeV/c);Acceptance",100,0.,20.);

  ////////////////////////////////////////////////////
  // Error bars for Efficiency and pT Distributions //
  ////////////////////////////////////////////////////

  // For Efficiencies
  histPtEfficiencyJpsiPrompt->Sumw2();
  histPtPairEfficiencyJpsiPrompt->Sumw2();
  histPtEfficiencyPsi2sPrompt->Sumw2();
  histPtPairEfficiencyPsi2sPrompt->Sumw2();

  histPtEfficiencyJpsiNonPrompt->Sumw2();
  histPtPairEfficiencyJpsiNonPrompt->Sumw2();
  histPtEfficiencyPsi2sNonPrompt->Sumw2();
  histPtPairEfficiencyPsi2sNonPrompt->Sumw2();

  // For pT Distributions
  histPtJpsiPrompt->Sumw2();
  histPtPairJpsiPrompt->Sumw2();
  histPtPsi2sPrompt->Sumw2();
  histPtPairPsi2sPrompt->Sumw2();

  histPtJpsiPromptAfterCuts->Sumw2();
  histPtPairJpsiPromptAfterCuts->Sumw2();
  histPtPsi2sPromptAfterCuts->Sumw2();
  histPtPairPsi2sPromptAfterCuts->Sumw2();

  histPtJpsiNonPrompt->Sumw2();
  histPtPairJpsiNonPrompt->Sumw2();
  histPtPsi2sNonPrompt->Sumw2();
  histPtPairPsi2sNonPrompt->Sumw2();

  histPtJpsiNonPromptAfterCuts->Sumw2();
  histPtPairJpsiNonPromptAfterCuts->Sumw2();
  histPtPsi2sNonPromptAfterCuts->Sumw2();
  histPtPairPsi2sNonPromptAfterCuts->Sumw2();

  ////////////////////////////////////////////////
  // Declare some bin parameters for BW Fitting //
  ////////////////////////////////////////////////

  // For J/Psi Mass
  int   divisionMass  = histMassPairJpsiPromptAfterCuts->GetNbinsX();
  float massMIN       = histMassPairJpsiPromptAfterCuts->GetBinLowEdge(1);
  float massMAX       = histMassPairJpsiPromptAfterCuts->GetBinLowEdge(divisionMass+1);
  float BIN_SIZE_MASS = histMassPairJpsiPromptAfterCuts->GetBinWidth(1);

  // For Psi2s Mass
  int   divisionMass_psi2s  = histMassPairPsi2sPromptAfterCuts->GetNbinsX();
  float massMIN_psi2s       = histMassPairPsi2sPromptAfterCuts->GetBinLowEdge(1);
  float massMAX_psi2s       = histMassPairPsi2sPromptAfterCuts->GetBinLowEdge(divisionMass_psi2s+1);
  float BIN_SIZE_MASS_psi2s = histMassPairPsi2sPromptAfterCuts->GetBinWidth(1);

  /// TPave Settings for Statistics and Fitting Box
  gStyle->SetOptStat("KSiouRMen");
  gStyle->SetOptFit(11112);

  /// Variables for Reading Kinematics
  Double_t invMassPair, ptPair;
  Double_t m1, m2;
  Double_t px, py;
  Double_t ptMoth, yMoth, etaMoth, pMoth, phiMoth, thetaMoth, massMoth;

  /// Variables For Applying cuts at Numerator Level
  Bool_t ptCut;
  Bool_t etaCut;
  Bool_t yCut;

  /// Variables For Applying cuts at Generator Level
  Bool_t generatorLevelYCut;

  /// Variables for Checking Prompt / Non-Prompt
  Bool_t PromptCheck    = kFALSE;
  Bool_t NonPromptCheck = kFALSE;

  /// J/Psi and Psi(2S) Pdg Codes
  Int_t mothPDG[] = {443, 100443};

  // Prefix for signal or background
  const char* nameprefix = "sgn_1";

  FairLogger::GetLogger()->SetLogScreenLevel("DEBUG");
  TFile f(o2::base::NameConf::getMCKinematicsFileName(nameprefix).c_str());

  LOG(debug) << "Checking input file :" << f.GetPath();

  std::vector<o2::MCTrack>* mctracks = nullptr;
  auto tr = (TTree*)f.Get("o2sim");
  if(!tr) printf("Warning!!! File not Found. \n");
  assert(tr);

  auto mcbr = tr->GetBranch("MCTrack");
  assert(mcbr);
  mcbr->SetAddress(&mctracks);

  std::vector<o2::TrackReference>* trackrefs = nullptr;
  auto refbr = tr->GetBranch("TrackRefs");
  assert(refbr);
  refbr->SetAddress(&trackrefs);

  o2::steer::MCKinematicsReader mcreader(nameprefix, o2::steer::MCKinematicsReader::Mode::kMCKine);

  // when present we also read some hits for ITS to test consistency of trackID assignments
  TFile hitf(o2::base::DetectorNameConf::getHitsFileName(o2::detectors::DetID::ITS, nameprefix).c_str());
  auto hittr = (TTree*)hitf.Get("o2sim");
  auto hitbr = hittr ? hittr->GetBranch("ITSHit") : nullptr;
  std::vector<o2::itsmft::Hit>* hits = nullptr;
  if (hitbr) {
    hitbr->SetAddress(&hits);
  }

  for (int eventID = 0; eventID < mcbr->GetEntries(); ++eventID) {
    mcbr->GetEntry(eventID);
    refbr->GetEntry(eventID);
    LOG(debug) << "-- Entry --" << eventID;
    LOG(debug) << "Have " << mctracks->size() << " tracks";

    std::unordered_map<int, bool> trackidsinITS_fromhits;
    if (hitbr) {
      hitbr->GetEntry(eventID);
      //LOG(debug) << "Have " << hits->size() << " hits";

      // check that trackIDs from the hits are within range
      int maxid = 0;
      for (auto& h : *hits) {
        maxid = std::max(maxid, h.GetTrackID());
        trackidsinITS_fromhits[h.GetTrackID()] = true;
        assert(maxid < mctracks->size());
      }
    }

    int ti = 0;
    int tiprimary = 0;

    int primaries = 0;
    int physicalprimaries = 0;
    int secondaries = 0;

    // record tracks that left a hit in TPC and ITS
    // (we know that these tracks should then have a TrackRef)
    std::vector<int> trackidsinTPC;
    std::vector<int> trackidsinITS;


    //                                                            Beauty Hadrons PDGs                                                       //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 511: B0 and anti-B0 || 521: B- and B+ || 531: B_s0 and anti-B_s0 || 541: B_c+ and B_c- || 5112: anti-Sigma_b+ and Sigma_b-           //
    // 5122: anti-Lambda_b0 and Lambda_b0 || 5132: Xi_b- and anti-Xi_b+ || 5232: anti-Xi_b0  and Xi_b0 || 5332: anti-Omega_b+ and Omega_b-  //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Int_t bpdgs[] = {511, 521, 531, 541, 5112, 5122, 5132, 5232, 5332};
    Int_t sizePdg = sizeof(bpdgs)/sizeof(Int_t);

    // Main Loop For MC Tracks
    for (auto& t : *mctracks) {
      if (t.isSecondary()) {
        // check that mother indices are monotonic
        // for primaries, this may be different (for instance with Pythia8)
        assert(ti > t.getMotherTrackId());
      }
      if (t.leftTrace(o2::detectors::DetID::TPC)) {
        trackidsinTPC.emplace_back(ti);
      }
      if (t.leftTrace(o2::detectors::DetID::ITS)) {
        trackidsinITS.emplace_back(ti);
      }
      bool physicalPrim = o2::mcutils::MCTrackNavigator::isPhysicalPrimary(t, *mctracks);
      //LOG(debug) << " track " << ti << "\t" << t.getMotherTrackId() << " hits " << t.hasHits() << " isPhysicalPrimary " << physicalPrim;
      if (t.isPrimary()) {
        primaries++;
      } else {
        secondaries++;
      }
      if (physicalPrim) {
        physicalprimaries++;
      }
      ti++;
      if(t.isPrimary()) tiprimary++; else continue;
      // check that mother indices are reasonable
      Int_t ifirst = t.getFirstDaughterTrackId();
      Int_t ilast = t.getLastDaughterTrackId();
      Bool_t hasBeautyMoth = kFALSE;
      if( (TMath::Abs(t.GetPdgCode()) == mothPDG[0]) || (TMath::Abs(t.GetPdgCode()) == mothPDG[1] ) ){
        o2::MCTrack *lepton=0x0; o2::MCTrack *antilepton=0x0;
        Int_t idMoth = t.getMotherTrackId();
        Bool_t isPrompt = idMoth <= 0 ? kTRUE : kFALSE;
        if(!isPrompt){ //  check beauty mother
          auto tdM = mcreader.getTrack(eventID, idMoth);
          for(int i=0; i<sizePdg; i++){ if (TMath::Abs(tdM->GetPdgCode()) == bpdgs[i] ) hasBeautyMoth = kTRUE; }
        }

        /// Get Mother Tracks. Use getter.
        ptMoth = t.GetPt();
        pMoth = t.GetP();
        phiMoth = t.GetPhi();
        thetaMoth = t.GetTheta();
        yMoth = t.GetRapidity();
        etaMoth = t.GetEta();
        //massMoth = t.GetMass();

        /// Check Daugthers is Primary
        ifirst  = (mcreader.getTrack(eventID, ifirst))->isPrimary() ? ifirst : -1;
        ilast  = (mcreader.getTrack(eventID, ilast))->isPrimary() ? ilast : -1;

        //LOG(DEBUG) << " mother track - pdg " << t.GetPdgCode() << " first daughter "<<  ifirst <<" last daughter " << ilast << " position " << ti;

        for(int idaugh=ifirst; idaugh<ilast+1; idaugh++ ){
          auto td = mcreader.getTrack(eventID, idaugh);
          if(td->GetPdgCode() == -1*leptonPDG) lepton = (o2::MCTrack*)td;
          if(td->GetPdgCode() == leptonPDG) antilepton = (o2::MCTrack*)td;
        }

        if((!lepton) || (!antilepton)) continue;
        // evaluate inv mass, pt, y of pairs before cuts
        m1 = TDatabasePDG::Instance()->GetParticle(leptonPDG)->Mass();
        m2 = TDatabasePDG::Instance()->GetParticle(leptonPDG)->Mass();
        invMassPair =  m1*m1+m2*m2 + 2.0*(TMath::Sqrt(m1*m1+lepton->GetP()*lepton->GetP())*TMath::Sqrt(m2*m2+antilepton->GetP()*antilepton->GetP()) - lepton->Px()*antilepton->Px() - lepton->Py()*antilepton->Py() - lepton->Pz()*antilepton->Pz());
        invMassPair = TMath::Sqrt(invMassPair);
        ////
        px = lepton->Px()+antilepton->Px();
        py = lepton->Py()+antilepton->Py();
        ptPair = TMath::Sqrt(px*px + py*py);

        /////////////////////////////////
        // Fill Histograms Before Cuts //
        /////////////////////////////////

        /// fiducial cut on generator level for the definition of the acceptance before the application of the single track cuts.
        if(midysim     == kTRUE){generatorLevelYCut = TMath::Abs(yMoth) < 0.9;}
        if(forwardysim == kTRUE){generatorLevelYCut = yMoth < -2.5 && yMoth > -4.0;}

        // Applying to cuts at generator level
        if(generatorLevelYCut == kFALSE) continue;

        //// fill prompt
        if(isPrompt){
          if(TMath::Abs(t.GetPdgCode()) == mothPDG[0]){
            // Before Cuts
            histPtPairJpsiPrompt->Fill(ptPair);
            histMassPairJpsiPrompt->Fill(invMassPair);
            histPtJpsiPrompt->Fill(ptMoth);
            histYJpsiPrompt->Fill(yMoth);
            histEtaJpsiPrompt->Fill(etaMoth);
            histThetaJpsiPrompt->Fill(thetaMoth);
            histPhiJpsiPrompt->Fill(phiMoth);
            histIsPrompt->Fill(0.5);
          }else{
            // Before Cuts
            histPtPairPsi2sPrompt->Fill(ptPair);
            histMassPairPsi2sPrompt->Fill(invMassPair);
            histPtPsi2sPrompt->Fill(ptMoth);
            histYPsi2sPrompt->Fill(yMoth);
            histEtaPsi2sPrompt->Fill(etaMoth);
            histThetaPsi2sPrompt->Fill(thetaMoth);
            histPhiPsi2sPrompt->Fill(phiMoth);
            histIsPrompt->Fill(2.5);
          }
        }else if(!isPrompt && hasBeautyMoth){
          //// fill non prompt
          if(TMath::Abs(t.GetPdgCode()) == mothPDG[0]){
            histPtPairJpsiNonPrompt->Fill(ptMoth);
            histMassPairJpsiNonPrompt->Fill(invMassPair);
            histPtJpsiNonPrompt->Fill(ptMoth);
            histYJpsiNonPrompt->Fill(yMoth);
            histEtaJpsiNonPrompt->Fill(etaMoth);
            histThetaJpsiNonPrompt->Fill(thetaMoth);
            histPhiJpsiNonPrompt->Fill(phiMoth);
            histIsPrompt->Fill(1.5);
            //}
          }else{
            histPtPairPsi2sNonPrompt->Fill(ptPair);
            histMassPairPsi2sNonPrompt->Fill(invMassPair);
            histPtPsi2sNonPrompt->Fill(ptMoth);
            histYPsi2sNonPrompt->Fill(yMoth);
            histEtaPsi2sNonPrompt->Fill(etaMoth);
            histThetaPsi2sNonPrompt->Fill(thetaMoth);
            histPhiPsi2sNonPrompt->Fill(phiMoth);
            histIsPrompt->Fill(3.5);
          }
        }
        //////////////////////////////////////
        // Cut Definition and Applying Part //
        //////////////////////////////////////

        //                        DEFAULT CUTS                   //
        ///////////////////////////////////////////////////////////
        //  For Mid Rapidity,     pT > 1   &&      |eta| <  0.9  //
        //  For Forward Rapidity, pT > 1   &&  -4 < eta  < -2.5  //
        ///////////////////////////////////////////////////////////

        // Cut Definition Part. You can add specific kinematic cuts here.
        ptCut = lepton->GetPt() > 1.0 && antilepton->GetPt() > 1.0;
        if(midysim == kTRUE){
          etaCut    =   TMath::Abs(lepton->GetEta()) < 0.9 && TMath::Abs(antilepton->GetEta()) < 0.9;
          yCut      =   TMath::Abs(lepton->GetRapidity()) < 0.9 && TMath::Abs(antilepton->GetRapidity()) < 0.9;
        }
        if(forwardysim == kTRUE){
          etaCut    =   lepton->GetEta() < -2.5 && antilepton->GetEta() < -2.5 && lepton->GetEta() > -4.0 && antilepton->GetEta() > -4.0 ;
          yCut      =   lepton->GetRapidity() < -2.5 && antilepton->GetRapidity() < -2.5 && lepton->GetRapidity() > -4.0 && antilepton->GetRapidity() > -4.0;
        }

        // Applying The Cuts
        if(ptCut  == kFALSE) continue;
        if(etaCut == kFALSE) continue;
        if(yCut   == kFALSE) continue;

        //// re-evaluate inv mass, pt, y of pairs after cuts
        m1 = TDatabasePDG::Instance()->GetParticle(leptonPDG)->Mass();
        m2 = TDatabasePDG::Instance()->GetParticle(leptonPDG)->Mass();
        invMassPair =  m1*m1+m2*m2 + 2.0*(TMath::Sqrt(m1*m1+lepton->GetP()*lepton->GetP())*TMath::Sqrt(m2*m2+antilepton->GetP()*antilepton->GetP()) - lepton->Px()*antilepton->Px() - lepton->Py()*antilepton->Py() - lepton->Pz()*antilepton->Pz());
        invMassPair = TMath::Sqrt(invMassPair);
        ////
        px = lepton->Px()+antilepton->Px();
        py = lepton->Py()+antilepton->Py();
        ptPair = TMath::Sqrt(px*px + py*py);

        ////////////////////////////////
        // Fill Histograms After Cuts //
        ////////////////////////////////

        //// fill prompt
        if(isPrompt){
          if(TMath::Abs(t.GetPdgCode()) == mothPDG[0]){
            histPtPairJpsiPromptAfterCuts->Fill(ptPair);
            histMassPairJpsiPromptAfterCuts->Fill(invMassPair);
            histPtJpsiPromptAfterCuts->Fill(ptMoth);
            histYJpsiPromptAfterCuts->Fill(yMoth);
            histEtaJpsiPromptAfterCuts->Fill(etaMoth);
            histThetaJpsiPromptAfterCuts->Fill(thetaMoth);
            histPhiJpsiPromptAfterCuts->Fill(phiMoth);
            histPtEfficiencyJpsiPrompt->Fill(ptMoth);
            histPtPairEfficiencyJpsiPrompt->Fill(ptPair);
            histIsPromptAfterCuts->Fill(0.5);
          }else{
            histPtPairPsi2sPromptAfterCuts->Fill(ptPair);
            histMassPairPsi2sPromptAfterCuts->Fill(invMassPair);
            histPtPsi2sPromptAfterCuts->Fill(ptMoth);
            histYPsi2sPromptAfterCuts->Fill(yMoth);
            histEtaPsi2sPromptAfterCuts->Fill(etaMoth);
            histThetaPsi2sPromptAfterCuts->Fill(thetaMoth);
            histPhiPsi2sPromptAfterCuts->Fill(phiMoth);
            histPtEfficiencyPsi2sPrompt->Fill(ptMoth);
            histPtPairEfficiencyPsi2sPrompt->Fill(ptPair);
            histIsPromptAfterCuts->Fill(2.5);
          }
        }else if(!isPrompt && hasBeautyMoth){
          //// fill non prompts
          if(TMath::Abs(t.GetPdgCode()) == mothPDG[0]){
            histPtPairJpsiNonPromptAfterCuts->Fill(ptPair);
            histMassPairJpsiNonPromptAfterCuts->Fill(invMassPair);
            histPtJpsiNonPromptAfterCuts->Fill(ptMoth);
            histYJpsiNonPromptAfterCuts->Fill(yMoth);
            histEtaJpsiNonPromptAfterCuts->Fill(etaMoth);
            histThetaJpsiNonPromptAfterCuts->Fill(thetaMoth);
            histPhiJpsiNonPromptAfterCuts->Fill(phiMoth);
            histPtEfficiencyJpsiNonPrompt->Fill(ptMoth);
            histPtPairEfficiencyJpsiNonPrompt->Fill(ptPair);
            histIsPromptAfterCuts->Fill(1.5);
            //}
          }else{
            histPtPairPsi2sNonPromptAfterCuts->Fill(ptPair);
            histMassPairPsi2sNonPromptAfterCuts->Fill(invMassPair);
            histPtPsi2sNonPromptAfterCuts->Fill(ptMoth);
            histYPsi2sNonPromptAfterCuts->Fill(yMoth);
            histEtaPsi2sNonPromptAfterCuts->Fill(etaMoth);
            histThetaPsi2sNonPromptAfterCuts->Fill(thetaMoth);
            histPhiPsi2sNonPromptAfterCuts->Fill(phiMoth);
            histPtEfficiencyPsi2sNonPrompt->Fill(ptMoth);
            histPtPairEfficiencyPsi2sNonPrompt->Fill(ptPair);
            histIsPromptAfterCuts->Fill(3.5);
          }
        }
      }
      //Delete some pointer variables for don't get to memory leak
      //o2::MCTrack lepton = n
    }

    if (hitbr) {
      assert(trackidsinITS.size() == trackidsinITS_fromhits.size());
      for (auto id : trackidsinITS) {
        assert(trackidsinITS_fromhits[id] == true);
      }
    }
    // Count TPC Tracks and Debugging
    LOG(debug) << "Have " << trackidsinTPC.size() << " tracks with hits in TPC";
    LOG(debug) << "Have " << trackidsinITS.size() << " tracks with hits in ITS";
    //LOG(debug) << "Have " << trackrefs->size() << " track refs";
    //LOG(info) << "Have " << primaries << " primaries and " << physicalprimaries << " physical primaries";

    // check correct working of MCKinematicsReader
    bool havereferences = trackrefs->size();
    if (havereferences) {
      for (auto& trackID : trackidsinTPC) {
        auto trackrefs = mcreader.getTrackRefs(eventID, trackID);
        assert(trackrefs.size() > 0);
        //LOG(debug) << " Track " << trackID << " has " << trackrefs.size() << " TrackRefs";
        for (auto& ref : trackrefs) {
          assert(ref.getTrackID() == trackID);
        }
      }
    }
  }

  ////////////////////////////////////////////
  // Histogram Filling Checking for Drawing //
  ////////////////////////////////////////////

  if(histPtJpsiPrompt->GetEntries()    > 0)   {  PromptCheck    = kTRUE; }
  if(histPtJpsiNonPrompt->GetEntries() > 0)   {  NonPromptCheck = kTRUE; }

  /// Set Parameters to Breit-Wigner Fit Functions for Mass, pT and Y

  TF1 *bwFuncMass_jpsi = new TF1("mybwMass_jpsi",mybwMass,massMIN, massMAX,3);
  bwFuncMass_jpsi->SetParameter(0,1.0);     bwFuncMass_jpsi->SetParName(0,"const");
  bwFuncMass_jpsi->SetParameter(2,5.0);     bwFuncMass_jpsi->SetParName(1,"sigma");
  bwFuncMass_jpsi->SetParameter(1,95.0);    bwFuncMass_jpsi->SetParName(2,"mean");

  TF1 *bwFuncMass_psi2s = new TF1("mybwMass_psi2s",mybwMass,massMIN_psi2s, massMAX_psi2s,3);
  bwFuncMass_psi2s->SetParameter(0,10.0);    bwFuncMass_psi2s->SetParName(0,"const");
  bwFuncMass_psi2s->SetParameter(2,5.0);     bwFuncMass_psi2s->SetParName(1,"sigma");
  bwFuncMass_psi2s->SetParameter(1,95.0);    bwFuncMass_psi2s->SetParName(2,"mean");


  ////////////////////////////////////////
  // Efficiencies Part Calculation Part //
  ////////////////////////////////////////

  // Efficiency = N_Accepted (After Cuts) / N_Recorded (Before Cuts)

  if(PromptCheck == kTRUE){
    // Efficiency for Prompt
    histPtEfficiencyJpsiPrompt->Divide(histPtJpsiPrompt);
    histPtPairEfficiencyJpsiPrompt->Divide(histPtPairJpsiPrompt);
    histPtEfficiencyPsi2sPrompt->Divide(histPtPsi2sPrompt);
    histPtPairEfficiencyPsi2sPrompt->Divide(histPtPairPsi2sPrompt);

    // Set Range Efficiencies for Prompt
    histPtEfficiencyJpsiPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
    histPtPairEfficiencyJpsiPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
    histPtEfficiencyPsi2sPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
    histPtPairEfficiencyPsi2sPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
  }

  if(NonPromptCheck == kTRUE){
    // Efficiency for Non-Prompt
    histPtEfficiencyJpsiNonPrompt->Divide(histPtJpsiNonPrompt);
    histPtPairEfficiencyJpsiNonPrompt->Divide(histPtPairJpsiNonPrompt);
    histPtEfficiencyPsi2sNonPrompt->Divide(histPtPsi2sNonPrompt);
    histPtPairEfficiencyPsi2sNonPrompt->Divide(histPtPairPsi2sNonPrompt);

    // Set Range Efficiencies for Non-Prompt
    histPtEfficiencyJpsiNonPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
    histPtPairEfficiencyJpsiNonPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
    histPtEfficiencyPsi2sNonPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
    histPtPairEfficiencyPsi2sNonPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
  }

  /// Writing Histograms to Disk (Output will histKine.root file) --> For Before Cuts
  TFile foutput1("histKine.root","RECREATE");

  /// Write Histograms for Prompt J/Psi
  if(PromptCheck == kTRUE){
    histMassPairJpsiPrompt->Write();
    histPtPairJpsiPrompt->Write();
    histPtJpsiPrompt->Write();
    histYJpsiPrompt->Write();
    histEtaJpsiPrompt->Write();
    histThetaJpsiPrompt->Write();
    histPhiJpsiPrompt->Write();
    histMassPairPsi2sPrompt->Write();
    histPtPairPsi2sPrompt->Write();
    histPtPsi2sPrompt->Write();
    histYPsi2sPrompt->Write();
    histEtaPsi2sPrompt->Write();
    histThetaPsi2sPrompt->Write();
    histPhiPsi2sPrompt->Write();
  }
  /// Write Histograms for Non - Prompt J/Psi
  if(NonPromptCheck == kTRUE){
    histMassPairJpsiNonPrompt->Write();
    histPtPairJpsiNonPrompt->Write();
    histPtJpsiNonPrompt->Write();
    histYJpsiNonPrompt->Write();
    histEtaJpsiNonPrompt->Write();
    histThetaJpsiNonPrompt->Write();
    histPhiJpsiNonPrompt->Write();
    histMassPairPsi2sNonPrompt->Write();
    histPtPairPsi2sNonPrompt->Write();
    histPtPsi2sNonPrompt->Write();
    histYPsi2sNonPrompt->Write();
    histEtaPsi2sNonPrompt->Write();
    histThetaPsi2sNonPrompt->Write();
    histPhiPsi2sNonPrompt->Write();
  }
  histIsPrompt->Write();

  /// Writing Histograms to Disk (Output will histCutKine.root file) --> For After Cuts
  TFile foutput2("histCutKine.root","RECREATE");

  /// Write Histograms for Prompt J/Psi and Psi(2S)
  if(PromptCheck == kTRUE){
    histMassPairJpsiPromptAfterCuts->Write();
    histPtPairJpsiPromptAfterCuts->Write();
    histPtJpsiPromptAfterCuts->Write();
    histYJpsiPromptAfterCuts->Write();
    histEtaJpsiPromptAfterCuts->Write();
    histThetaJpsiPromptAfterCuts->Write();
    histPhiJpsiPromptAfterCuts->Write();
    histMassPairPsi2sPromptAfterCuts->Write();
    histPtPairPsi2sPromptAfterCuts->Write();
    histPtPsi2sPromptAfterCuts->Write();
    histYPsi2sPromptAfterCuts->Write();
    histEtaPsi2sPromptAfterCuts->Write();
    histThetaPsi2sPromptAfterCuts->Write();
    histPhiPsi2sPromptAfterCuts->Write();
  }
  /// Write Histograms for Non - Prompt J/Psi and Psi(2S)
  if(NonPromptCheck == kTRUE){
    histMassPairJpsiNonPromptAfterCuts->Write();
    histPtPairJpsiNonPromptAfterCuts->Write();
    histPtJpsiNonPromptAfterCuts->Write();
    histYJpsiNonPromptAfterCuts->Write();
    histEtaJpsiNonPromptAfterCuts->Write();
    histThetaJpsiNonPromptAfterCuts->Write();
    histPhiJpsiNonPromptAfterCuts->Write();
    histMassPairPsi2sNonPromptAfterCuts->Write();
    histPtPairPsi2sNonPromptAfterCuts->Write();
    histPtPsi2sNonPromptAfterCuts->Write();
    histYPsi2sNonPromptAfterCuts->Write();
    histEtaPsi2sNonPromptAfterCuts->Write();
    histThetaPsi2sNonPromptAfterCuts->Write();
    histPhiPsi2sNonPromptAfterCuts->Write();
  }
  histIsPromptAfterCuts->Write();

  /// Writing Histograms to Disk (Output will histCompareAfterCuts.root file) --> For Comparing J/Psi and Psi(2S)
  TFile foutput3("histCompareAfterCuts.root","RECREATE");

  // Setting the Line Colours for comparing J/psi and Psi(2S)
  histMassPairPsi2sPromptAfterCuts->SetLineColor(kGreen+3);
  histPtPairPsi2sPromptAfterCuts->SetLineColor(kGreen+3);
  histPtPsi2sPromptAfterCuts->SetLineColor(kGreen+3);
  histYPsi2sPromptAfterCuts->SetLineColor(kGreen+3);

  histMassPairPsi2sNonPromptAfterCuts->SetLineColor(kGreen+3);
  histPtPairPsi2sNonPromptAfterCuts->SetLineColor(kGreen+3);
  histPtPsi2sNonPromptAfterCuts->SetLineColor(kGreen+3);
  histYPsi2sNonPromptAfterCuts->SetLineColor(kGreen+3);

  // Draw Histograms on Canvas

  // Prompt J/Psi and Psi(2S) Compare Plots
  canvasPromptMassCompareAfterCuts->Divide(1,1);
  canvasPromptMassCompareAfterCuts->cd(1);
  histMassPairJpsiPromptAfterCuts->Draw();
  histMassPairPsi2sPromptAfterCuts->Draw("same");

  canvasPromptPtPairCompareAfterCuts->Divide(1,1);
  canvasPromptPtPairCompareAfterCuts->cd(1);
  histPtPairJpsiPromptAfterCuts->Draw();
  histPtPairPsi2sPromptAfterCuts->Draw("same");

  canvasPromptPtCompareAfterCuts->Divide(1,1);
  canvasPromptPtCompareAfterCuts->cd(1);
  histPtJpsiPromptAfterCuts->Draw();
  histPtPsi2sPromptAfterCuts->Draw("same");

  canvasPromptyCompareAfterCuts->Divide(1,1);
  canvasPromptyCompareAfterCuts->cd(1);
  histYJpsiPromptAfterCuts->Draw();
  histYPsi2sPromptAfterCuts->Draw("same");

  // Non-Prompt J/Psi and Psi(2S) Compare Plots
  canvasNonPromptMassCompareAfterCuts->Divide(1,1);
  canvasNonPromptMassCompareAfterCuts->cd(1);
  histMassPairJpsiNonPromptAfterCuts->Draw();
  histMassPairPsi2sNonPromptAfterCuts->Draw("same");

  canvasNonPromptPtPairCompareAfterCuts->Divide(1,1);
  canvasNonPromptPtPairCompareAfterCuts->cd(1);
  histPtPairJpsiNonPromptAfterCuts->Draw();
  histPtPairPsi2sNonPromptAfterCuts->Draw("same");

  canvasNonPromptPtCompareAfterCuts->Divide(1,1);
  canvasNonPromptPtCompareAfterCuts->cd(1);
  histPtJpsiNonPromptAfterCuts->Draw();
  histPtPsi2sNonPromptAfterCuts->Draw("same");

  canvasNonPromptyCompareAfterCuts->Divide(1,1);
  canvasNonPromptyCompareAfterCuts->cd(1);
  histYJpsiNonPromptAfterCuts->Draw();
  histYPsi2sNonPromptAfterCuts->Draw("same");


  /// Write Histograms for Prompt J/Psi and Psi(2S)

  // Prompt
  canvasPromptMassCompareAfterCuts->Write();
  canvasPromptPtCompareAfterCuts->Write();
  canvasPromptPtPairCompareAfterCuts->Write();
  canvasPromptyCompareAfterCuts->Write();

  // Non-Prompt
  canvasNonPromptMassCompareAfterCuts->Write();
  canvasNonPromptPtCompareAfterCuts->Write();
  canvasNonPromptPtPairCompareAfterCuts->Write();
  canvasNonPromptyCompareAfterCuts->Write();

  //////////////////////////////////////////////
  // Applying Breit Wigner Fit For After Cuts //
  //////////////////////////////////////////////

  if(PromptCheck == kTRUE){
    ///Fit Inv Mass Histograms for Prompt J/Psi and Psi(2S)
    histMassPairJpsiPromptAfterCuts->Fit("mybwMass_jpsi","QR");       TF1 *fitMassPairJpsiPromptAfterCuts   = histMassPairJpsiPrompt->GetFunction("mybwMass_jpsi");
    histMassPairPsi2sPromptAfterCuts->Fit("mybwMass_psi2s","QR");     TF1 *fitMassPairPsi2sPromptAfterCuts  = histMassPairPsi2sPromptAfterCuts->GetFunction("mybwMass_psi2s");
  }
  if(NonPromptCheck == kTRUE){
    /// Fit Inv Mass Histograms for Non - Prompt J/Psi and Non-Prompt Psi(2S)
    histMassPairJpsiNonPromptAfterCuts->Fit("mybwMass_jpsi","QR");    TF1 *fitMassPairJpsiNonPromptAfterCuts   = histMassPairJpsiNonPromptAfterCuts->GetFunction("mybwMass_jpsi");
    histMassPairPsi2sNonPromptAfterCuts->Fit("mybwMass_psi2s","QR");  TF1 *fitMassPairPsi2sNonPromptAfterCuts  = histMassPairPsi2sNonPromptAfterCuts->GetFunction("mybwMass_psi2s");
  }

  /// Writing Histograms to Disk (Output will histFittingAfterCuts.root file) --> For Fitted (Breit - Wigner) Histos After Cuts
  TFile foutput4("histFittingAfterCuts.root","RECREATE");

  /// Write Histograms for Prompt J/Psi and Psi(2S)
  if(PromptCheck == kTRUE){
    //histEtaPromptAfterCuts->Write("PE0");
    histMassPairJpsiPromptAfterCuts->Write("PE0");
    histPtPairJpsiPromptAfterCuts->Write("PE0");
    histPtJpsiPromptAfterCuts->Write("PE0");
    histYJpsiPromptAfterCuts->Write("PE0");
    histMassPairPsi2sPromptAfterCuts->Write("PE0");
    histPtPairPsi2sPromptAfterCuts->Write("PE0");
    histPtPsi2sPromptAfterCuts->Write("PE0");
    histYPsi2sPromptAfterCuts->Write("PE0");
  }
  /// Write Histograms for Non - Prompt J/Psi and Psi(2S)
  if(NonPromptCheck == kTRUE){
    //histEtaNonPromptAfterCuts->Write("PE0");
    histMassPairJpsiNonPromptAfterCuts->Write("PE0");
    histPtPairJpsiNonPromptAfterCuts->Write("PE0");
    histPtJpsiNonPromptAfterCuts->Write("PE0");
    histYJpsiNonPromptAfterCuts->Write("PE0");
    histMassPairPsi2sNonPromptAfterCuts->Write("PE0");
    histPtPairPsi2sNonPromptAfterCuts->Write("PE0");
    histPtPsi2sNonPromptAfterCuts->Write("PE0");
    histYPsi2sNonPromptAfterCuts->Write("PE0");
  }
  /// Writing Histograms to Disk (Output will histAcceptance.root file) --> For Acceptance Factor
  TFile foutput5("histAcceptance.root","RECREATE");

  // Prompt J/Psi and Psi(2S)
  if(PromptCheck == kTRUE){
  histPtEfficiencyJpsiPrompt->Write();
  histPtPairEfficiencyJpsiPrompt->Write();
  histPtEfficiencyPsi2sPrompt->Write();
  histPtPairEfficiencyPsi2sPrompt->Write();
  }
  // Non-Prompt J/psi and Psi(2S)
  if(NonPromptCheck == kTRUE){
  histPtEfficiencyJpsiNonPrompt->Write();
  histPtPairEfficiencyJpsiNonPrompt->Write();
  histPtEfficiencyPsi2sNonPrompt->Write();
  histPtPairEfficiencyPsi2sNonPrompt->Write();
  }

  LOG(info) << "STACK TEST SUCCESSFULL\n";
  return 0;
}
