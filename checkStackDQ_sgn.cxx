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
#include "TH2F.h"
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

int checkStackDQ_sgn()
{
	// Define Rapidity, Lepton PDG and Sim Rapidity Type
	Int_t kSelectionPdg;
	Double_t kYRangeMin, kYRangeMax;
	Bool_t kMidYSimulation = kFALSE; Bool_t kForwardYSimulation = kFALSE;

	// User Guide
	std::cout << "============================================================================================" << '\n';
	std::cout << "          checkStackDQ.cxx is Executable macro to check functioning of stack"                 << '\n';
	std::cout << "          This macro has been developed for the Charmonium Analysis of the PWG-DQ group"      << '\n';
	std::cout << "          It provides analyses kinematics and track references of a kinematics file"          << '\n';
	std::cout << "                             Usage of this macro for signals"                                 << '\n';
	std::cout << "============================================================================================" << '\n';

	// Enter PDG code daughter for selection
	std::cout << "============================================================================================" << '\n';
	std::cout << "                Selection for daughters (dilepton pairs). You should enter a PDG code."       << '\n';
	std::cout << "                      For electron 11 or For Muon 13 and then press enter."                   << '\n';
	std::cout << "============================================================================================" << '\n';
	std::cin >> kSelectionPdg;
	//if(kSelectionPdg != 13 || kSelectionPdg != 11){printf("Your LeptonPDG selection is wrong!");}
	assert(kSelectionPdg == 13 || kSelectionPdg == 11);

	// set eta, rapidity, phi and theta ranges
	Int_t kEtaRangeMin, kEtaRangeMax;      	kEtaRangeMin = -5; kEtaRangeMax = 5;
	Int_t kPhiRangeMin, kPhiRangeMax;      	kPhiRangeMin = -6; kPhiRangeMax = 6;
	Int_t kThetaRangeMin, kThetaRangeMax;  	kThetaRangeMin = -5; kThetaRangeMax = 6;
	Int_t kVertexRangeMin, kVertexRangeMax;	kVertexRangeMin = -15; kVertexRangeMax = 15;

	// set pT and pT Pair Range
	Int_t kPtRangeMin, kPtRangeMax;					kPtRangeMin = 0;	kPtRangeMax = 15;
	Int_t kPtPairRangeMin, kPtPairRangeMax;	kPtPairRangeMin = 0; kPtPairRangeMax = 20;

	// set histogram bins
	// TODO: HİSTOGRAMLAR İÇİN EN UYGUN BİNLERİ BUL.
	Int_t kBinRange_100  = 100;
	Int_t kBinRange_250  = 250;
	Int_t kBinRange_500  = 500;
	Int_t kBinRange_750  = 750;
	Int_t kBinRange_1000 = 1000;

	// Settings for Selections
	if(kSelectionPdg == 11) { kYRangeMin = -1.5; kYRangeMax = 1.5;  kMidYSimulation     = kTRUE; }
	if(kSelectionPdg == 13) { kYRangeMin = -4.3; kYRangeMax = -2.3; kForwardYSimulation = kTRUE; }

	/// TPave Settings for Statistics and Fitting Box
	gStyle->SetOptStat("KSiouRMen");
	gStyle->SetOptFit(11112);

	/// Variables for Reading Kinematics
	Float_t kInvMassPair, kPtPair, kVertexPair;
	Float_t kMassPair1, kMassPair2;
	Float_t kPx, kPy;
	Float_t kVx, kVy, kVz;
	Float_t kPt, kY, kEta, kP, kPhi, kTheta;
	Float_t kPrimaryVx, kPrimaryVy, kPrimaryVz;
	Float_t kMass;

	/// Variables for Reading Kinematics For Beauties Vertexes
	Float_t kCandidateBzeroVx, kCandidateBzeroVy;
	Float_t kCandidateBplusVx, kCandidateBplusVy;
	Float_t kCandidateBszeroVx, kCandidateBszeroVy;
	Float_t kCandidateBcplusVx, kCandidateBcplusVy;
	Float_t kCandidateSigmaBVx, kCandidateSigmaBVy;
	Float_t kCandidateLambdaBzeroVx, kCandidateLambdaBzeroVy;
	Float_t kCandidateXiBminusVx, kCandidateXiBminusVy;
	Float_t kCandidateXiBzeroVx, kCandidateXiBzeroVy;
	Float_t kCandidateOmegaBminusVx, kCandidateOmegaBminusVy;

	/// Variables for Reading Kinematics For Quarkonium
	Float_t kCandidateJpsiVx, kCandidateJpsiVy;
	Float_t kCandidatePsi2SVx, kCandidatePsi2SVy;


	/// Variables for Reaing Kinematics After Cuts
	Float_t kInvMassPairBranchCut, kPtPairBranchCut;
	Float_t kPtBranchCut, kYBranchCut, kEtaBranchCut, kPhiBranchCut, kThetaBranchCut;

	/// Variables For Applying cuts at Numerator Level
	Bool_t kPtCut;
	Bool_t kEtaCut;
	Bool_t kYCut;
	Bool_t kInvMassCutJpsi, kInvMassCutPsi2s;

	// Set pT Cut Range
	Float_t kPtCutRangeMin = 1.0;
	Float_t kPtCutRangeMax = 100;

	// For Forward Rapidity Set eta and Y Cut Ranges
	Float_t kMuonQualityEtaCutMin = -4.0;
	Float_t kMuonQualityEtaCutMax = -2.5;

	Float_t kMuonQualityYCutMin = -4.0;
	Float_t kMuonQualityYCutMax = -2.5;

	// For Mid Rapidity Set eta and Y Cut Ranges
	Float_t kTrackBarrelYCutMin = -0.9;
	Float_t kTrackBarrelYCutMax = 0.9;

	Float_t kTrackBarrelEtaCutMin = -0.9;
	Float_t kTrackBarrelEtaCutMax = 0.9;

	/// Variables For Applying cuts at Generator Level
	Bool_t kGeneratorLevelYCut;
	Bool_t kGeneratorLevelEtaCut;

	// For Forward Rapidity Set Generator Level Y Cut
	Float_t kGeneratorLevelYCutRangeMinMuonQuality = -4.0;
	Float_t kGeneratorLevelYCutRangeMaxMuonQuality = -2.5;

	// For Mid Rapidity Set Generator Level Y Cut
	Float_t kGeneratorLevelYCutRangeMinTrackBarrel = 0.9;
	Float_t kGeneratorLevelYCutRangeMaxTrackBarrel = 0.9;

	/// Variables for Checking Prompt / Non-Prompt
	Bool_t kPromptCheck    = kFALSE;
	Bool_t kNonPromptCheck = kFALSE;

	// TODO : Check the histogram ranges and optimize them.
	// TODO: Px Py Pz ve vertex için histogramlar ekle.

	////////////////////////////////////////////////////////////
	/// Histogram Definitions for Basic Kinematic Properties ///
	////////////////////////////////////////////////////////////

	/// define histos prompt jpsi / psi2s before cuts

	// J/Psi
	TH1F *histPtPairJpsiPrompt = new TH1F("PtPairJpsiPromptBeforeCuts","p_{T} Pair Distribution;p_{T}(GeV/c);Entries",kBinRange_100,kPtPairRangeMin,kPtPairRangeMax);
	TH1F *histMassPairJpsiPrompt = new TH1F("MassJpsiPromptBeforeCuts","Invariant Mass;m_{ll}(GeV/c^{2});Entries",5./0.02,0.,5.);
	TH1F *histYJpsiPrompt = new TH1F("YJpsiPromptBeforeCuts","Y Distribution;y;Entries",kBinRange_100,kYRangeMin,kYRangeMax);
	TH1F *histPtJpsiPrompt = new TH1F("PtJpsiPromptBeforeCuts","p_{T} Distribution;p_{T}(GeV/c);Entries",kBinRange_100,kPtRangeMin,kPtRangeMax);
	TH1F *histEtaJpsiPrompt = new TH1F("EtaJpsiPromptBeforeCuts","#eta Distribution;#eta;Entries",kBinRange_100,kEtaRangeMin,kEtaRangeMax);
	TH1F *histThetaJpsiPrompt = new TH1F("ThetaJpsiPromptBeforeCuts","#theta Distribution;#theta (rad.);Entries",kBinRange_100,kThetaRangeMin,kThetaRangeMax);
	TH1F *histPhiJpsiPrompt = new TH1F("PhiJpsiPromptBeforeCuts","#varphi Distribution;#varphi (rad.);Entries",kBinRange_100,kPhiRangeMin,kPhiRangeMax);
	// TODO : Vertex için histolar oluştur böyle.
	//TH1F *histVxPairJpsiPrompt = new TH1F("VxPairJpsiPromptBeforeCuts",";p_{T}(GeV/c);Entries",kBinRange_100,kPtRangeMin,kPtPairRangeMax);
	//TH1F *histVyPairJpsiPrompt = new TH1F("VyPairJpsiPromptBeforeCuts",";p_{T}(GeV/c);Entries",kBinRange_100,kPtRangeMin,kPtPairRangeMax);
	//TH1F *histVzPairJpsiPrompt = new TH1F("VzPairJpsiPromptBeforeCuts",";p_{T}(GeV/c);Entries",kBinRange_100,kPtRangeMin,kPtPairRangeMax);



	// Psi(2s)
	TH1F *histPtPairPsi2sPrompt = new TH1F("PtPairPsi2sPromptBeforeCuts","p_{T} Pair Distribution;p_{T}(GeV/c);Entries",kBinRange_100,kPtPairRangeMin,kPtPairRangeMax);
	TH1F *histMassPairPsi2sPrompt = new TH1F("MassPsi2sPromptBeforeCuts","Invariant Mass;m_{ll}(GeV/c^{2});Entries",5./0.02,0.,5.);
	TH1F *histYPsi2sPrompt = new TH1F("YPsi2sPromptBeforeCuts","Y Distribution;y;Entries",kBinRange_100,kYRangeMin,kYRangeMax);
	TH1F *histPtPsi2sPrompt = new TH1F("PtPsi2sPromptBeforeCuts","p_{T} Distribution;p_{T}(GeV/c);Entries",kBinRange_100,kPtRangeMin,kPtRangeMax);
	TH1F *histEtaPsi2sPrompt = new TH1F("EtaPsi2sPromptBeforeCuts","#eta Distribution;#eta;Entries",kBinRange_100,kEtaRangeMin,kEtaRangeMax);
	TH1F *histThetaPsi2sPrompt = new TH1F("ThetaPsi2sPromptBeforeCuts","#theta Distribution;#theta (rad.);Entries",kBinRange_100,kThetaRangeMin,kThetaRangeMax);
	TH1F *histPhiPsi2sPrompt = new TH1F("PhiPsi2sPromptBeforeCuts","#varphi Distribution;#varphi (rad.);Entries",kBinRange_100,kPhiRangeMin,kPhiRangeMax);

	/// define histos prompt jpsi / psi2s after cuts

	// J/Psi
	TH1F *histPtPairJpsiPromptAfterCuts = new TH1F("PtPairJpsiPromptAfterCuts","p_{T} Pair Distribution;p_{T}(GeV/c);Entries",kBinRange_100,kPtPairRangeMin,kPtPairRangeMax);
	TH1F *histMassPairJpsiPromptAfterCuts = new TH1F("MassJpsiPromptAfterCuts","Invariant Mass;m_{ll}(GeV/c^{2});Entries",5./0.02,0.,5.);
	TH1F *histYJpsiPromptAfterCuts = new TH1F("YJpsiPromptAfterCuts","Y Distribution;y;Entries",kBinRange_100,kYRangeMin,kYRangeMax);
	TH1F *histPtJpsiPromptAfterCuts = new TH1F("PtJpsiPromptAfterCuts","p_{T} Distribution;p_{T}(GeV/c);Entries",kBinRange_100,kPtRangeMin,kPtRangeMax);
	TH1F *histEtaJpsiPromptAfterCuts = new TH1F("EtaJpsiPromptAfterCuts","#eta Distribution;#eta;Entries",kBinRange_100,kEtaRangeMin,kEtaRangeMax);
	TH1F *histThetaJpsiPromptAfterCuts = new TH1F("ThetaJpsiPromptAfterCuts","#theta Distribution;#theta (rad.);Entries",kBinRange_100,kThetaRangeMin,kThetaRangeMax);
	TH1F *histPhiJpsiPromptAfterCuts = new TH1F("PhiJpsiPromptAfterCuts","#varphi Distribution;#varphi (rad.);Entries",kBinRange_100,kPhiRangeMin,kPhiRangeMax);
	// Psi(2s)
	TH1F *histPtPairPsi2sPromptAfterCuts = new TH1F("PtPairPsi2sPromptAfterCuts","p_{T} Pair Distribution;p_{T}(GeV/c);Entries",kBinRange_100,kPtPairRangeMin,kPtPairRangeMax);
	TH1F *histMassPairPsi2sPromptAfterCuts = new TH1F("MassPsi2sPromptAfterCuts","Invariant Mass;m_{ll}(GeV/c^{2});Entries",5./0.02,0.,5.);
	TH1F *histYPsi2sPromptAfterCuts = new TH1F("YPsi2sPromptAfterCuts","Y Distribution;y;Entries",kBinRange_100,kYRangeMin,kYRangeMax);
	TH1F *histPtPsi2sPromptAfterCuts = new TH1F("PtPsi2sPromptAfterCuts","p_{T} Distribution;p_{T}(GeV/c);Entries",kBinRange_100,kPtRangeMin,kPtRangeMax);
	TH1F *histEtaPsi2sPromptAfterCuts = new TH1F("EtaPsi2sPromptAfterCuts","#eta Distribution;#eta;Entries",kBinRange_100,kEtaRangeMin,kEtaRangeMax);
	TH1F *histThetaPsi2sPromptAfterCuts = new TH1F("ThetaPsi2sPromptAfterCuts","#theta Distribution;#theta (rad.);Entries",kBinRange_100,kThetaRangeMin,kThetaRangeMax);
	TH1F *histPhiPsi2sPromptAfterCuts = new TH1F("PhiPsi2sPromptAfterCuts","#varphi Distribution;#varphi (rad.);Entries",kBinRange_100,kPhiRangeMin,kPhiRangeMax);

	/// define histos non-prompt jpsi / psi2s before cuts

	// J/Psi
	TH1F *histPtPairJpsiNonPrompt = new TH1F("PtPairNonPromptJpsiBeforeCuts","p_{T} Pair Distribution;p_{T}(GeV/c);Entries",kBinRange_100,kPtPairRangeMin,kPtPairRangeMax);
	TH1F *histMassPairJpsiNonPrompt = new TH1F("MassNonPromptJpsiBeforeCuts","Invariant Mass;m_{ll}(GeV/c^{2});Entries",5./0.02,0.,5.);
	TH1F *histYJpsiNonPrompt = new TH1F("YJpsiNonPromptBeforeCuts","Y Distribution;y;Entries",kBinRange_100,kYRangeMin,kYRangeMax);
	TH1F *histPtJpsiNonPrompt = new TH1F("PtJpsiNonPromptBeforeCuts","p_{T} Distribution;p_{T}(GeV/c);Entries",kBinRange_100,kPtRangeMin,kPtRangeMax);
	TH1F *histEtaJpsiNonPrompt = new TH1F("EtaJpsiNonPromptBeforeCuts","#eta Distribution;#eta (rad.);Entries",kBinRange_100,kEtaRangeMin,kEtaRangeMax);
	TH1F *histThetaJpsiNonPrompt = new TH1F("ThetaJpsiNonPromptBeforeCuts","#theta Distribution;#theta (rad.);Entries",kBinRange_100,kThetaRangeMin,kThetaRangeMax);
	TH1F *histPhiJpsiNonPrompt = new TH1F("PhiJpsiNonPromptBeforeCuts","#varphi Distribution;#varphi (rad.);Entries",kBinRange_100,kPhiRangeMin,kPhiRangeMax);
	// BURAYA VERTEX EKLEDİM CHECK ET.
	TH1F *histVxPairJpsiNonPrompt = new TH1F("VxPairNonPromptJpsiBeforeCuts",";p_{T}(GeV/c);Entries",1000,kVertexRangeMin,kVertexRangeMax);
	TH1F *histVyPairJpsiNonPrompt = new TH1F("VyPairNonPromptJpsiBeforeCuts",";p_{T}(GeV/c);Entries",1000,kVertexRangeMin,kVertexRangeMax);
	TH1F *histVzPairJpsiNonPrompt = new TH1F("VzPairNonPromptJpsiBeforeCuts",";p_{T}(GeV/c);Entries",1000,kVertexRangeMin,kVertexRangeMax);
	TH1F *histVertexPairJpsiNonPrompt = new TH1F("VertexPairNonPromptJpsiBeforeCuts",";p_{T}(GeV/c);Entries",1000,kVertexRangeMin,kVertexRangeMax);
	// Psi(2s)
	TH1F *histPtPairPsi2sNonPrompt = new TH1F("PtPairPsi2sNonPromptBeforeCuts","p_{T} Pair Distribution;p_{T}(GeV/c);Entries",kBinRange_100,kPtPairRangeMin,kPtPairRangeMax);
	TH1F *histMassPairPsi2sNonPrompt = new TH1F("MassPsi2sNonPromptBeforeCuts","Invariant Mass;m_{ll}(GeV/c^{2});Entries",5./0.02,0.,5.);
	TH1F *histYPsi2sNonPrompt = new TH1F("YPsi2sNonPromptBeforeCuts","Y Distribution;y;Entries",kBinRange_100,kYRangeMin,kYRangeMax);
	TH1F *histPtPsi2sNonPrompt = new TH1F("PtPsi2sNonPromptBeforeCuts","p_{T} Distribution;p_{T}(GeV/c);Entries",kBinRange_100,kPtRangeMin,kPtRangeMax);
	TH1F *histEtaPsi2sNonPrompt = new TH1F("EtaPsi2sNonPromptBeforeCuts","#eta Distribution;#eta;Entries",kBinRange_100,kEtaRangeMin,kEtaRangeMax);
	TH1F *histThetaPsi2sNonPrompt = new TH1F("ThetaPsi2sNonPromptBeforeCuts","#theta Distribution;#theta (rad.);Entries",kBinRange_100,kThetaRangeMin,kThetaRangeMax);
	TH1F *histPhiPsi2sNonPrompt = new TH1F("PhiPsi2sNonPromptBeforeCuts","#varphi Distribution;#varphi (rad.);Entries",kBinRange_100,kPhiRangeMin,kPhiRangeMax);

	/// define histos non-prompt jpsi / psi2s after cuts

	// J/Psi
	TH1F *histPtPairJpsiNonPromptAfterCuts = new TH1F("PtPairNonPromptJpsiAfterCuts","p_{T} Pair Distribution;p_{T}(GeV/c);Entries",kBinRange_100,kPtPairRangeMin,kPtPairRangeMax);
	TH1F *histMassPairJpsiNonPromptAfterCuts = new TH1F("MassNonPromptJpsiAfterCuts","Invariant Mass;m_{ll}(GeV/c^{2});Entries",5./0.02,0.,5.);
	TH1F *histYJpsiNonPromptAfterCuts = new TH1F("YJpsiNonPromptJpsiAfterCuts","Y Distribution;y;Entries",kBinRange_100,kYRangeMin,kYRangeMax);
	TH1F *histPtJpsiNonPromptAfterCuts = new TH1F("PtJpsiNonPromptJpsiAfterCuts","p_{T} Distribution;p_{T}(GeV/c);Entries",kBinRange_100,kPtRangeMin,kPtRangeMax);
	TH1F *histEtaJpsiNonPromptAfterCuts = new TH1F("EtaJpsiNonPromptAfterCuts","#eta Distribution;#eta;Entries",kBinRange_100,kEtaRangeMin,kEtaRangeMax);
	TH1F *histThetaJpsiNonPromptAfterCuts = new TH1F("ThetaJpsiNonPromptAfterCuts","#theta Distribution;#theta (rad.);Entries",kBinRange_100,kThetaRangeMin,kThetaRangeMax);
	TH1F *histPhiJpsiNonPromptAfterCuts = new TH1F("PhiJpsiNonPromptAfterCuts","#varphi Distribution;#varphi (rad.);Entries",kBinRange_100,kPhiRangeMin,kPhiRangeMax);
	// Psi(2s)
	TH1F *histPtPairPsi2sNonPromptAfterCuts = new TH1F("PtPairPsi2sNonPromptAfterCuts","p_{T} Pair Distribution;p_{T}(GeV/c);Entries",kBinRange_100,kPtPairRangeMin,kPtPairRangeMax);
	TH1F *histMassPairPsi2sNonPromptAfterCuts = new TH1F("MassPsi2sNonPromptAfterCuts","Invariant Mass;m_{ll}(GeV/c^{2});Entries",5./0.02,0.,5.);
	TH1F *histYPsi2sNonPromptAfterCuts = new TH1F("YPsi2sNonPromptAfterCuts","Y Distribution;y;Entries",kBinRange_100,kYRangeMin,kYRangeMax);
	TH1F *histPtPsi2sNonPromptAfterCuts = new TH1F("PtPsi2sNonPromptAfterCuts","p_{T} Distribution;p_{T}(GeV/c);Entries",kBinRange_100,kPtRangeMin,kPtRangeMax);
	TH1F *histEtaPsi2sNonPromptAfterCuts = new TH1F("EtaPsi2sNonPromptAfterCuts","#eta Distribution;#eta;Entries",kBinRange_100,kEtaRangeMin,kEtaRangeMax);
	TH1F *histThetaPsi2sNonPromptAfterCuts = new TH1F("ThetaPsi2sNonPromptAfterCuts","#theta Distribution;#theta (rad.);Entries",kBinRange_100,kThetaRangeMin,kThetaRangeMax);
	TH1F *histPhiPsi2sNonPromptAfterCuts = new TH1F("PhiPsi2sNonPromptAfterCuts","#varphi Distribution;#varphi (rad.);Entries",kBinRange_100,kPhiRangeMin,kPhiRangeMax);

	/// define histos for Origin before cuts
	TH1D *histOrigin = new TH1D("histOriginBeforeCuts","",5,0,5);
	histOrigin->GetXaxis()->SetBinLabel(1,"IsJpsiPrompt");
	histOrigin->GetXaxis()->SetBinLabel(2,"IsJpsiFromB");
	histOrigin->GetXaxis()->SetBinLabel(3,"IsPsi2sPrompt");
	histOrigin->GetXaxis()->SetBinLabel(4,"IsPsi2sFromB");
	histOrigin->GetXaxis()->SetBinLabel(5,"IsOtherParticle");
	histOrigin->SetTitle("Charmonia Origin Before Cuts");
	histOrigin->GetXaxis()->SetTitle("Quarkonium Candidates");
	histOrigin->GetYaxis()->SetTitle("Entries");

	/// define histos for Origin after cuts
	TH1D *histOriginAfterCuts = new TH1D("histOriginAfterCuts","",5,0,5);
	histOriginAfterCuts->GetXaxis()->SetBinLabel(1,"IsJpsiPrompt");
	histOriginAfterCuts->GetXaxis()->SetBinLabel(2,"IsJpsiFromB");
	histOriginAfterCuts->GetXaxis()->SetBinLabel(3,"IsPsi2sPrompt");
	histOriginAfterCuts->GetXaxis()->SetBinLabel(4,"IsPsi2sFromB");
	histOriginAfterCuts->GetXaxis()->SetBinLabel(5,"IsOtherParticle");
	histOriginAfterCuts->SetTitle("Charmonia Origin After Cuts");
	histOriginAfterCuts->GetXaxis()->SetTitle("Quarkonium Candidates");
	histOriginAfterCuts->GetYaxis()->SetTitle("Entries");

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

	/////////////////////////////////////////////////////////
	///  Histogram Definitions for Extract to Acceptance  ///
	/////////////////////////////////////////////////////////

	// Prompt
	TH1F *histPtEfficiencyJpsiPrompt = new TH1F("Prompt J/#it{#psi} Acceptance for p_{T}",";p_{T}(GeV/c);Acceptance",kBinRange_100,kPtRangeMin,kPtRangeMax);
	TH1F *histPtPairEfficiencyJpsiPrompt = new TH1F("Prompt J/#it{#psi} Acceptance for p_{T} Pairs",";p_{T}(GeV/c);Acceptance",kBinRange_100,kPtPairRangeMin,kPtPairRangeMax);
	TH1F *histPtEfficiencyPsi2sPrompt = new TH1F("Prompt #it{#psi}(2S) Acceptance for p_{T}",";p_{T}(GeV/c);Acceptance",kBinRange_100,kPtRangeMin,kPtRangeMax);
	TH1F *histPtPairEfficiencyPsi2sPrompt = new TH1F("Prompt #it{#psi}(2S) Acceptance for p_{T} Pairs",";p_{T}(GeV/c);Acceptance",kBinRange_100,kPtPairRangeMin,kPtPairRangeMax);
	TH1F *histYEfficiencyJpsiPrompt = new TH1F("Prompt J/#it{#psi} Acceptance for Y",";y;Acceptance",kBinRange_100,kYRangeMin,kYRangeMax);
	TH1F *histYEfficiencyPsi2sPrompt = new TH1F("Prompt #it{#psi}(2S) Acceptance for Y",";y;Acceptance",kBinRange_100,kYRangeMin,kYRangeMax);

	// Non-Prompt
	TH1F *histPtEfficiencyJpsiNonPrompt = new TH1F("NonPrompt J/#it{#psi} Acceptance for p_{T}",";p_{T}(GeV/c);Acceptance",kBinRange_100,kPtRangeMin,kPtRangeMax);
	TH1F *histPtPairEfficiencyJpsiNonPrompt = new TH1F("NonPrompt J/#it{#psi} Acceptance for p_{T} Pairs",";p_{T}(GeV/c);Acceptance",kBinRange_100,kPtPairRangeMin,kPtPairRangeMax);
	TH1F *histPtEfficiencyPsi2sNonPrompt = new TH1F("NonPrompt #it{#psi}(2S) Acceptance for p_{T}",";p_{T}(GeV/c);Acceptance",kBinRange_100,kPtRangeMin,kPtRangeMax);
	TH1F *histPtPairEfficiencyPsi2sNonPrompt = new TH1F("NonPrompt #it{#psi}(2S) Acceptance for p_{T} Pairs",";p_{T}(GeV/c);Acceptance",kBinRange_100,kPtPairRangeMin,kPtPairRangeMax);
	TH1F *histYEfficiencyJpsiNonPrompt = new TH1F("NonPrompt J/#it{#psi} Acceptance for Y",";y;Acceptance",kBinRange_100,kYRangeMin,kYRangeMax);
	TH1F *histYEfficiencyPsi2sNonPrompt = new TH1F("NonPrompt #it{#psi}(2S) Acceptance for Y",";y;Acceptance",kBinRange_100,kYRangeMin,kYRangeMax);

	///////////////////////////////////////////////////////////
	/// Histogram Definitions for Check Beauties Statistics ///
	///////////////////////////////////////////////////////////

	// From J/psi before cuts
	TH1D *histBeautyMesonsFromJpsi = new TH1D("histBeautyMesonsFromJpsi","h_{B};Entries",4,0,4);
	histBeautyMesonsFromJpsi->GetXaxis()->SetBinLabel(1,"B_{0}");
	histBeautyMesonsFromJpsi->GetXaxis()->SetBinLabel(2,"B^{+}");
	histBeautyMesonsFromJpsi->GetXaxis()->SetBinLabel(3,"B_{s}^{0}");
	histBeautyMesonsFromJpsi->GetXaxis()->SetBinLabel(4,"B_{c}^{+}");
	histBeautyMesonsFromJpsi->SetTitle("Beauty Mesons From J/#it{#psi} Before Cuts");
	histBeautyMesonsFromJpsi->GetXaxis()->SetTitle("h_{B} Candidates");
	histBeautyMesonsFromJpsi->GetYaxis()->SetTitle("Entries");

	TH1D *histBeautyBaryonsFromJpsi = new TH1D("histBeautyBaryonsFromJpsi","h_{B};Entries",5,0,5);
	histBeautyBaryonsFromJpsi->GetXaxis()->SetBinLabel(1,"#Sigma_{b}^{-}");
	histBeautyBaryonsFromJpsi->GetXaxis()->SetBinLabel(2,"#Lambda_{B}^{-}");
	histBeautyBaryonsFromJpsi->GetXaxis()->SetBinLabel(3,"#Xi_{B}^{-}");
	histBeautyBaryonsFromJpsi->GetXaxis()->SetBinLabel(4,"#Xi_{B}^{0}");
	histBeautyBaryonsFromJpsi->GetXaxis()->SetBinLabel(5,"#Omega_{B}^{-}");
	histBeautyBaryonsFromJpsi->SetTitle("Beauty Baryons From J/#it{#psi} Before Cuts");
	histBeautyBaryonsFromJpsi->GetXaxis()->SetTitle("h_{B} Candidates");
	histBeautyBaryonsFromJpsi->GetYaxis()->SetTitle("Entries");

	// From psi(2S) before cuts
	TH1D *histBeautyMesonsFromPsi2s = new TH1D("histBeautyMesonsFromPsi2s","h_{B};Entries",4,0,4);
	histBeautyMesonsFromPsi2s->GetXaxis()->SetBinLabel(1,"B_{0}");
	histBeautyMesonsFromPsi2s->GetXaxis()->SetBinLabel(2,"B^{+}");
	histBeautyMesonsFromPsi2s->GetXaxis()->SetBinLabel(3,"B_{s}^{0}");
	histBeautyMesonsFromPsi2s->GetXaxis()->SetBinLabel(4,"B_{c}^{+}");
	histBeautyMesonsFromPsi2s->SetTitle("Beauty Mesons From #psi(2S) Before Cuts");
	histBeautyMesonsFromPsi2s->GetXaxis()->SetTitle("h_{B} Candidates");
	histBeautyMesonsFromPsi2s->GetYaxis()->SetTitle("Entries");

	TH1D *histBeautyBaryonsFromPsi2s = new TH1D("histBeautyBaryonsFromPsi2s","h_{B};Entries",5,0,5);
	histBeautyBaryonsFromPsi2s->GetXaxis()->SetBinLabel(1,"#Sigma_{b}^{-}");
	histBeautyBaryonsFromPsi2s->GetXaxis()->SetBinLabel(2,"#Lambda_{B}^{-}");
	histBeautyBaryonsFromPsi2s->GetXaxis()->SetBinLabel(3,"#Xi_{B}^{-}");
	histBeautyBaryonsFromPsi2s->GetXaxis()->SetBinLabel(4,"#Xi_{B}^{0}");
	histBeautyBaryonsFromPsi2s->GetXaxis()->SetBinLabel(5,"#Omega_{B}^{-}");
	histBeautyBaryonsFromPsi2s->SetTitle("Beauty Baryons From #psi(2S) Before Cuts");
	histBeautyBaryonsFromPsi2s->GetXaxis()->SetTitle("h_{B} Candidates");
	histBeautyBaryonsFromPsi2s->GetYaxis()->SetTitle("Entries");

	// From J/psi after cuts
	TH1D *histBeautyMesonsFromJpsiAfterCuts = new TH1D("histBeautyMesonsFromJpsiAfterCuts","h_{B};Entries",4,0,4);
	histBeautyMesonsFromJpsiAfterCuts->GetXaxis()->SetBinLabel(1,"B_{0}");
	histBeautyMesonsFromJpsiAfterCuts->GetXaxis()->SetBinLabel(2,"B^{+}");
	histBeautyMesonsFromJpsiAfterCuts->GetXaxis()->SetBinLabel(3,"B_{s}^{0}");
	histBeautyMesonsFromJpsiAfterCuts->GetXaxis()->SetBinLabel(4,"B_{c}^{+}");
	histBeautyMesonsFromJpsiAfterCuts->SetTitle("Beauty Mesons From J/#it{#psi} After Cuts");
	histBeautyMesonsFromJpsiAfterCuts->GetXaxis()->SetTitle("h_{B} Candidates");
	histBeautyMesonsFromJpsiAfterCuts->GetYaxis()->SetTitle("Entries");

	TH1D *histBeautyBaryonsFromJpsiAfterCuts = new TH1D("histBeautyBaryonsFromJpsiAfterCuts","h_{B};Entries",5,0,5);
	histBeautyBaryonsFromJpsiAfterCuts->GetXaxis()->SetBinLabel(1,"#Sigma_{b}^{-}");
	histBeautyBaryonsFromJpsiAfterCuts->GetXaxis()->SetBinLabel(2,"#Lambda_{B}^{-}");
	histBeautyBaryonsFromJpsiAfterCuts->GetXaxis()->SetBinLabel(3,"#Xi_{B}^{-}");
	histBeautyBaryonsFromJpsiAfterCuts->GetXaxis()->SetBinLabel(4,"#Xi_{B}^{0}");
	histBeautyBaryonsFromJpsiAfterCuts->GetXaxis()->SetBinLabel(5,"#Omega_{B}^{-}");
	histBeautyBaryonsFromJpsiAfterCuts->SetTitle("Beauty Baryons From J/#it{#psi} After Cuts");
	histBeautyBaryonsFromJpsiAfterCuts->GetXaxis()->SetTitle("h_{B} Candidates");
	histBeautyBaryonsFromJpsiAfterCuts->GetYaxis()->SetTitle("Entries");

	// From psi(2S) before cuts
	TH1D *histBeautyMesonsFromPsi2sAfterCuts = new TH1D("histBeautyMesonsFromPsi2sAfterCuts","h_{B};Entries",4,0,4);
	histBeautyMesonsFromPsi2sAfterCuts->GetXaxis()->SetBinLabel(1,"B_{0}");
	histBeautyMesonsFromPsi2sAfterCuts->GetXaxis()->SetBinLabel(2,"B^{+}");
	histBeautyMesonsFromPsi2sAfterCuts->GetXaxis()->SetBinLabel(3,"B_{s}^{0}");
	histBeautyMesonsFromPsi2sAfterCuts->GetXaxis()->SetBinLabel(4,"B_{c}^{+}");
	histBeautyMesonsFromPsi2sAfterCuts->SetTitle("Beauty Mesons From #psi(2S) Before Cuts");
	histBeautyMesonsFromPsi2sAfterCuts->GetXaxis()->SetTitle("h_{B} Candidates");
	histBeautyMesonsFromPsi2sAfterCuts->GetYaxis()->SetTitle("Entries");

	TH1D *histBeautyBaryonsFromPsi2sAfterCuts = new TH1D("histBeautyBaryonsFromPsi2sAfterCuts","h_{B};Entries",5,0,5);
	histBeautyBaryonsFromPsi2sAfterCuts->GetXaxis()->SetBinLabel(1,"#Sigma_{b}^{-}");
	histBeautyBaryonsFromPsi2sAfterCuts->GetXaxis()->SetBinLabel(2,"#Lambda_{B}^{-}");
	histBeautyBaryonsFromPsi2sAfterCuts->GetXaxis()->SetBinLabel(3,"#Xi_{B}^{-}");
	histBeautyBaryonsFromPsi2sAfterCuts->GetXaxis()->SetBinLabel(4,"#Xi_{B}^{0}");
	histBeautyBaryonsFromPsi2sAfterCuts->GetXaxis()->SetBinLabel(5,"#Omega_{B}^{-}");
	histBeautyBaryonsFromPsi2sAfterCuts->SetTitle("Beauty Baryons From #psi(2S) Before Cuts");
	histBeautyBaryonsFromPsi2sAfterCuts->GetXaxis()->SetTitle("h_{B} Candidates");
	histBeautyBaryonsFromPsi2sAfterCuts->GetYaxis()->SetTitle("Entries");

	//////////////////////////////////////////////////////////////////
	/// Histogram Definetions for Check Daughter Pairs Statistics ///
	////////////////////////////////////////////////////////////////

	// From NonPrompt + Prompt J/psi before cuts
	TH1D *histJpsiPairStatistics = new TH1D("histJpsiPairStatistics","",4,0,4);
	histJpsiPairStatistics->GetXaxis()->SetBinLabel(1,"e^{-}");
	histJpsiPairStatistics->GetXaxis()->SetBinLabel(2,"e^{+}");
	histJpsiPairStatistics->GetXaxis()->SetBinLabel(3,"#mu^{-}");
	histJpsiPairStatistics->GetXaxis()->SetBinLabel(4,"#mu^{+}");
	histJpsiPairStatistics->SetTitle("J/#it{#psi} Pair Statistics Before Cuts");
	histJpsiPairStatistics->GetXaxis()->SetTitle("Pairs");
	histJpsiPairStatistics->GetYaxis()->SetTitle("Entries");


	// From NonPrompt + Prompt psi(2S) before cuts
	TH1D *histPsi2sPairStatistics = new TH1D("histPsi2sPairStatistics","",4,0,4);
	histPsi2sPairStatistics->GetXaxis()->SetBinLabel(1,"e^{-}");
	histPsi2sPairStatistics->GetXaxis()->SetBinLabel(2,"e^{+}");
	histPsi2sPairStatistics->GetXaxis()->SetBinLabel(3,"#mu^{-}");
	histPsi2sPairStatistics->GetXaxis()->SetBinLabel(4,"#mu^{+}");
	histPsi2sPairStatistics->SetTitle("#psi(2S) Pair Statistics Before Cuts");
	histPsi2sPairStatistics->GetXaxis()->SetTitle("Pairs");
	histPsi2sPairStatistics->GetYaxis()->SetTitle("Entries");

	// From NonPrompt + Prompt J/psi after cuts
	TH1D *histJpsiPairStatisticsAfterCuts = new TH1D("histJpsiPairStatisticsAfterCuts","",4,0,4);
	histJpsiPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(1,"e^{-}");
	histJpsiPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(2,"e^{+}");
	histJpsiPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(3,"#mu^{-}");
	histJpsiPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(4,"#mu^{+}");
	histJpsiPairStatisticsAfterCuts->SetTitle("J/#it{#psi} Pair Statistics After Cuts");
	histJpsiPairStatisticsAfterCuts->GetXaxis()->SetTitle("Pairs");
	histJpsiPairStatisticsAfterCuts->GetYaxis()->SetTitle("Entries");

	// From NonPrompt + Prompt psi(2S) after cuts
	TH1D *histPsi2sPairStatisticsAfterCuts = new TH1D("histPsi2sPairStatisticsAfterCuts","",4,0,4);
	histPsi2sPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(1,"e^{-}");
	histPsi2sPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(2,"e^{+}");
	histPsi2sPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(3,"#mu^{-}");
	histPsi2sPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(4,"#mu^{+}");
	histPsi2sPairStatisticsAfterCuts->SetTitle("#psi(2S) Pair Statistics After Cuts");
	histPsi2sPairStatisticsAfterCuts->GetXaxis()->SetTitle("Pairs");
	histPsi2sPairStatisticsAfterCuts->GetYaxis()->SetTitle("Entries");

	// From Prompt J/psi before cuts
	TH1D *histPromptJpsiPairStatistics = new TH1D("histPromptJpsiPairStatistics","",4,0,4);
	histPromptJpsiPairStatistics->GetXaxis()->SetBinLabel(1,"e^{-}");
	histPromptJpsiPairStatistics->GetXaxis()->SetBinLabel(2,"e^{+}");
	histPromptJpsiPairStatistics->GetXaxis()->SetBinLabel(3,"#mu^{-}");
	histPromptJpsiPairStatistics->GetXaxis()->SetBinLabel(4,"#mu^{+}");
	histPromptJpsiPairStatistics->SetTitle("Prompt J/#it{#psi} Pair Statistics Before Cuts");
	histPromptJpsiPairStatistics->GetXaxis()->SetTitle("Pairs");
	histPromptJpsiPairStatistics->GetYaxis()->SetTitle("Entries");

	// From Prompt psi(2S) before cuts
	TH1D *histPromptPsi2sPairStatistics = new TH1D("histPromptPsi2sPairStatistics","",4,0,4);
	histPromptPsi2sPairStatistics->GetXaxis()->SetBinLabel(1,"e^{-}");
	histPromptPsi2sPairStatistics->GetXaxis()->SetBinLabel(2,"e^{+}");
	histPromptPsi2sPairStatistics->GetXaxis()->SetBinLabel(3,"#mu^{-}");
	histPromptPsi2sPairStatistics->GetXaxis()->SetBinLabel(4,"#mu^{+}");
	histPromptPsi2sPairStatistics->SetTitle("Prompt #psi(2S) Pair Statistics Before Cuts");
	histPromptPsi2sPairStatistics->GetXaxis()->SetTitle("Pairs");
	histPromptPsi2sPairStatistics->GetYaxis()->SetTitle("Entries");

	// From Prompt J/psi after cuts
	TH1D *histPromptJpsiPairStatisticsAfterCuts = new TH1D("histPromptJpsiPairStatisticsAfterCuts","",4,0,4);
	histPromptJpsiPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(1,"e^{-}");
	histPromptJpsiPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(2,"e^{+}");
	histPromptJpsiPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(3,"#mu^{-}");
	histPromptJpsiPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(4,"#mu^{+}");
	histPromptJpsiPairStatisticsAfterCuts->SetTitle("Prompt J/#it{#psi} Pair Statistics After Cuts");
	histPromptJpsiPairStatisticsAfterCuts->GetXaxis()->SetTitle("Pairs");
	histPromptJpsiPairStatisticsAfterCuts->GetYaxis()->SetTitle("Entries");

	// From Prompt psi(2S) after cuts
	TH1D *histPromptPsi2sPairStatisticsAfterCuts = new TH1D("histPromptPsi2sPairStatisticsAfterCuts","",4,0,4);
	histPromptPsi2sPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(1,"e^{-}");
	histPromptPsi2sPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(2,"e^{+}");
	histPromptPsi2sPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(3,"#mu^{-}");
	histPromptPsi2sPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(4,"#mu^{+}");
	histPromptPsi2sPairStatisticsAfterCuts->SetTitle("Prompt #psi(2S) Pair Statistics After Cuts");
	histPromptPsi2sPairStatisticsAfterCuts->GetXaxis()->SetTitle("Pairs");
	histPromptPsi2sPairStatisticsAfterCuts->GetYaxis()->SetTitle("Entries");

	// From NonPrompt J/psi before cuts
	TH1D *histNonPromptJpsiPairStatistics = new TH1D("histNonPromptJpsiPairStatistics","",4,0,4);
	histNonPromptJpsiPairStatistics->GetXaxis()->SetBinLabel(1,"e^{-}");
	histNonPromptJpsiPairStatistics->GetXaxis()->SetBinLabel(2,"e^{+}");
	histNonPromptJpsiPairStatistics->GetXaxis()->SetBinLabel(3,"#mu^{-}");
	histNonPromptJpsiPairStatistics->GetXaxis()->SetBinLabel(4,"#mu^{+}");
	histNonPromptJpsiPairStatistics->SetTitle("Non Prompt J/#it{#psi} Pair Statistics Before Cuts");
	histNonPromptJpsiPairStatistics->GetXaxis()->SetTitle("Pairs");
	histNonPromptJpsiPairStatistics->GetYaxis()->SetTitle("Entries");

	// From NonPrompt psi(2S) before cuts
	TH1D *histNonPromptPsi2sPairStatistics = new TH1D("histNonPromptPsi2sPairStatistics","",4,0,4);
	histNonPromptPsi2sPairStatistics->GetXaxis()->SetBinLabel(1,"e^{-}");
	histNonPromptPsi2sPairStatistics->GetXaxis()->SetBinLabel(2,"e^{+}");
	histNonPromptPsi2sPairStatistics->GetXaxis()->SetBinLabel(3,"#mu^{-}");
	histNonPromptPsi2sPairStatistics->GetXaxis()->SetBinLabel(4,"#mu^{+}");
	histNonPromptPsi2sPairStatistics->SetTitle("Non Prompt #psi(2S) Pair Statistics Before Cuts");
	histNonPromptPsi2sPairStatistics->GetXaxis()->SetTitle("Pairs");
	histNonPromptPsi2sPairStatistics->GetYaxis()->SetTitle("Entries");

	// From NonPrompt J/psi after cuts
	TH1D *histNonPromptJpsiPairStatisticsAfterCuts = new TH1D("histNonPromptJpsiPairStatisticsAfterCuts","",4,0,4);
	histNonPromptJpsiPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(1,"e^{-}");
	histNonPromptJpsiPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(2,"e^{+}");
	histNonPromptJpsiPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(3,"#mu^{-}");
	histNonPromptJpsiPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(4,"#mu^{+}");
	histNonPromptJpsiPairStatisticsAfterCuts->SetTitle("Non Prompt J/#it{#psi} Pair Statistics After Cuts");
	histNonPromptJpsiPairStatisticsAfterCuts->GetXaxis()->SetTitle("Pairs");
	histNonPromptJpsiPairStatisticsAfterCuts->GetYaxis()->SetTitle("Entries");

	// From NonPrompt psi(2S) after cuts
	TH1D *histNonPromptPsi2sPairStatisticsAfterCuts = new TH1D("histNonPromptPsi2sPairStatisticsAfterCuts","",4,0,4);
	histNonPromptPsi2sPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(1,"e^{-}");
	histNonPromptPsi2sPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(2,"e^{+}");
	histNonPromptPsi2sPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(3,"#mu^{-}");
	histNonPromptPsi2sPairStatisticsAfterCuts->GetXaxis()->SetBinLabel(4,"#mu^{+}");
	histNonPromptPsi2sPairStatisticsAfterCuts->SetTitle("Non Prompt #psi(2s) Pair Statistics After Cuts");
	histNonPromptPsi2sPairStatisticsAfterCuts->GetXaxis()->SetTitle("Pairs");
	histNonPromptPsi2sPairStatisticsAfterCuts->GetYaxis()->SetTitle("Entries");

	/////////////////////////////////////////////////////////
	/// Histogram Definitions for Check Versus Kinematics ///
	/////////////////////////////////////////////////////////

	// Prompt J/psi and psi(2s) before cuts
	TH2F *histPromptJpsiPhiVersusEta = new TH2F("histPromptJpsiPhiVersusEta","#varphi vs #eta;#varphi;#eta",100,-5,5,40,-5,5);
	TH2F *histPromptJpsiEtaVersusPt = new TH2F("histPromptJpsiEtaVersusPt","#eta vs p_{T};#eta;p_{T}",100,-5,5,100,0,20);
	TH2F *histPromptJpsiEtaVersusPtPair = new TH2F("histPromptJpsiEtaVersusPtPair","#eta vs p_{T} Pairs;#eta;p_{T} Pairs",100,-5,5,100,0,20);
	TH2F *histPromptJpsiMassVersusPt  = new TH2F("histPromptJpsiMassVersusPt ","InvMass vs p_{T};InvMass;p_{T}",100,0,5,100,0,20);
	TH2F *histPromptJpsiMassVersusPtPair  = new TH2F("histPromptJpsiMassVersusPtPair","InvMass vs p_{T} Pairs;InvMass;p_{T} Pairs",100,0,5,100,0,20);
	TH2F *histPromptJpsiMassVersusY  = new TH2F("histPromptJpsiMassVersusY","InvMass vs Y;InvMass;Y",100,0,5,100,-5,5);
	TH2F *histPromptJpsiPtVersusY  = new TH2F("histPromptJpsiPtVersusY","p_{T} vs Y;p_{T};Y",100,0,20,100,-5,5);
	TH2F *histPromptJpsiPtPairVersusY  = new TH2F("histPromptJpsiPtPairVersusY","p_{T} Pairs vs Y;p_{T} Pairs;Y",100,0,20,100,-5,5);
	TH2F *histPromptJpsiPtPairVersusPt  = new TH2F("histPromptJpsiPtPairVersusPt","p_{T} Pairs vs p_{T};p_{T} Pairs;p_{T}",100,0,20,100,0,20);

	TH2F *histPromptPsi2sPhiVersusEta = new TH2F("histPromptPsi2sPhiVersusEta","#varphi vs #eta;#varphi;#eta",100,-5,5,40,-5,5);
	TH2F *histPromptPsi2sEtaVersusPt = new TH2F("histPromptPsi2sEtaVersusPt","#eta vs p_{T};#eta;p_{T}",100,-5,5,100,0,20);
	TH2F *histPromptPsi2sEtaVersusPtPair = new TH2F("histPromptPsi2sEtaVersusPtPair","#eta vs p_{T} Pairs;#eta;p_{T} Pairs",100,-5,5,100,0,20);
	TH2F *histPromptPsi2sMassVersusPt  = new TH2F("histPromptPsi2sMassVersusPt ","InvMass vs p_{T};InvMass;p_{T}",100,0,5,100,0,20);
	TH2F *histPromptPsi2sMassVersusPtPair  = new TH2F("histPromptPsi2sMassVersusPtPair","InvMass vs p_{T} Pairs;InvMass;p_{T} Pairs",100,0,5,100,0,20);
	TH2F *histPromptPsi2sMassVersusY  = new TH2F("histPromptPsi2sMassVersusY","InvMass vs Y;InvMass;Y",100,0,5,100,-5,5);
	TH2F *histPromptPsi2sPtVersusY  = new TH2F("histPromptPsi2sPtVersusY","p_{T} vs Y;p_{T};Y",100,0,20,100,-5,5);
	TH2F *histPromptPsi2sPtPairVersusY  = new TH2F("histPromptPsi2sPtPairVersusY","p_{T} Pairs vs Y;p_{T} Pairs;Y",100,0,20,100,-5,5);
	TH2F *histPromptPsi2sPtPairVersusPt  = new TH2F("histPromptPsi2sPtPairVersusPt","p_{T} Pairs vs p_{T};p_{T} Pairs;p_{T}",100,0,20,100,0,20);

	// Prompt J/psi and psi(2s) after cuts
	TH2F *histPromptJpsiPhiVersusEtaAfterCuts = new TH2F("histPromptJpsiPhiVersusEtaAfterCuts","#varphi vs #eta;#varphi;#eta",100,-5,5,100,-5,5);
	TH2F *histPromptJpsiEtaVersusPtAfterCuts = new TH2F("histPromptJpsiEtaVersusPtAfterCuts","#eta vs p_{T};#eta;p_{T}",100,-5,5,100,0,20);
	TH2F *histPromptJpsiEtaVersusPtPairAfterCuts = new TH2F("histPromptJpsiEtaVersusPtPairAfterCuts","#eta vs p_{T} Pairs;#eta;p_{T} Pairs",100,-5,5,100,0,20);
	TH2F *histPromptJpsiMassVersusPtAfterCuts  = new TH2F("histPromptJpsiMassVersusPtAfterCuts","InvMass vs p_{T};InvMass;p_{T}",100,0,5,100,0,20);
	TH2F *histPromptJpsiMassVersusPtPairAfterCuts  = new TH2F("histPromptJpsiMassVersusPtPairAfterCuts","InvMass vs p_{T} Pairs;InvMass;p_{T} Pairs",100,0,5,100,0,20);
	TH2F *histPromptJpsiMassVersusYAfterCuts  = new TH2F("histPromptJpsiMassVersusYAfterCuts","InvMass vs Y;InvMass;Y",100,0,5,100,-5,5);
	TH2F *histPromptJpsiPtVersusYAfterCuts  = new TH2F("histPromptJpsiPtVersusYAfterCuts","p_{T} vs Y;p_{T};Y",100,0,20,100,-5,5);
	TH2F *histPromptJpsiPtPairVersusYAfterCuts  = new TH2F("histPromptJpsiPtPairVersusYAfterCuts","p_{T} Pairs vs Y;p_{T} Pairs;Y",100,0,20,100,-5,5);
	TH2F *histPromptJpsiPtPairVersusPtAfterCuts  = new TH2F("histPromptJpsiPtPairVersusPtAfterCuts","p_{T} Pairs vs p_{T};p_{T} Pairs;p_{T}",100,0,20,100,0,20);

	TH2F *histPromptPsi2sPhiVersusEtaAfterCuts = new TH2F("histPromptPsi2sPhiVersusEtaAfterCuts","#varphi vs #eta;#varphi;#eta",100,-5,5,100,-5,5);
	TH2F *histPromptPsi2sEtaVersusPtAfterCuts = new TH2F("histPromptPsi2sEtaVersusPtAfterCuts","#eta vs p_{T};#eta;p_{T}",100,-5,5,100,0,20);
	TH2F *histPromptPsi2sEtaVersusPtPairAfterCuts = new TH2F("histPromptPsi2sEtaVersusPtPairAfterCuts","#eta vs p_{T} Pairs;#eta;p_{T} Pairs",100,0,5,100,0,20);
	TH2F *histPromptPsi2sMassVersusPtAfterCuts  = new TH2F("histPromptPsi2sMassVersusPtAfterCuts","InvMass vs p_{T};InvMass;p_{T}",100,0,5,100,0,20);
	TH2F *histPromptPsi2sMassVersusPtPairAfterCuts  = new TH2F("histPromptPsi2sMassVersusPtPairAfterCuts","InvMass vs p_{T} Pairs;InvMass;p_{T} Pairs",100,0,5,100,0,20);
	TH2F *histPromptPsi2sMassVersusYAfterCuts  = new TH2F("histPromptPsi2sMassVersusYAfterCuts","InvMass vs Y;InvMass;Y",100,0,5,100,-5,5);
	TH2F *histPromptPsi2sPtVersusYAfterCuts  = new TH2F("histPromptPsi2sPtVersusYAfterCuts ","p_{T} vs Y;p_{T};Y",100,0,20,100,-5,5);
	TH2F *histPromptPsi2sPtPairVersusYAfterCuts  = new TH2F("histPromptPsi2sPtPairVersusYAfterCuts","p_{T} Pairs vs Y;p_{T} Pairs;Y",100,0,20,100,-5,5);
	TH2F *histPromptPsi2sPtPairVersusPtAfterCuts  = new TH2F("histPromptPsi2sPtPairVersusPtAfterCuts","p_{T} Pairs vs p_{T};p_{T} Pairs;p_{T}",100,0,20,100,0,20);

	// Non Prompt J/psi and psi(2s) before cuts
	TH2F *histNonPromptJpsiPhiVersusEta = new TH2F("histNonPromptJpsiPhiVersusEta","#varphi vs #eta;#varphi;#eta",100,-5,5,100,-5,5);
	TH2F *histNonPromptJpsiEtaVersusPt = new TH2F("histNonPromptJpsiEtaVersusPt","#eta vs p_{T};#eta;p_{T}",100,-5,5,100,0,20);
	TH2F *histNonPromptJpsiEtaVersusPtPair = new TH2F("histNonPromptJpsiEtaVersusPtPair","#eta v. p_{T} Pairs;#eta;p_{T} Pairs",100,-5,5,100,0,20);
	TH2F *histNonPromptJpsiMassVersusPt  = new TH2F("histNonPromptJpsiMassVersusPt","InvMass vs p_{T};InvMass;p_{T}",100,0,5,100,0,20);
	TH2F *histNonPromptJpsiMassVersusPtPair  = new TH2F("histNonPromptJpsiMassVersusPtPair","InvMass vs p_{T} Pairs;InvMass;p_{T} Pairs",100,0,5,100,0,20);
	TH2F *histNonPromptJpsiMassVersusY  = new TH2F("histNonPromptJpsiMassVersusY","InvMass vs Y;InvMass;Y",100,0,5,100,-5,5);
	TH2F *histNonPromptJpsiPtVersusY  = new TH2F("histNonPromptJpsiPtVersusY","p_{T} vs Y;p_{T};Y",100,0,20,100,-5,5);
	TH2F *histNonPromptJpsiPtPairVersusY  = new TH2F("histNonPromptJpsiPtPairVersusY","p_{T} Pairs vs Y;p_{T} Pairs;Y",100,0,20,100,-5,5);
	TH2F *histNonPromptJpsiPtPairVersusPt  = new TH2F("histNonPromptJpsiPtPairVersusPt","p_{T} Pairs vs p_{T};p_{T} Pairs;p_{T}",100,0,20,100,0,20);

	TH2F *histNonPromptPsi2sPhiVersusEta = new TH2F("histNonPromptPsi2sPhiVersusEta","#varphi vs #eta;#varphi;#eta",100,-5,5,100,-5,5);
	TH2F *histNonPromptPsi2sEtaVersusPt = new TH2F("histNonPromptPsi2sEtaVersusPt","#eta vs p_{T};#eta;p_{T}",100,-5,5,100,0,20);
	TH2F *histNonPromptPsi2sEtaVersusPtPair = new TH2F("histNonPromptPsi2sEtaVersusPtPair","#eta v. p_{T} Pairs;#eta;p_{T} Pairs",100,-5,5,100,0,20);
	TH2F *histNonPromptPsi2sMassVersusPt  = new TH2F("histNonPromptPsi2sMassVersusPt","InvMass vs p_{T};InvMass;p_{T}",100,-5,5,100,0,20);
	TH2F *histNonPromptPsi2sMassVersusPtPair  = new TH2F("histNonPromptPsi2sMassVersusPtPair","InvMass vs p_{T} Pairs;InvMass;p_{T} Pairs",100,-5,5,100,0,20);
	TH2F *histNonPromptPsi2sMassVersusY  = new TH2F("histNonPromptPsi2sMassVersusY","InvMass vs Y;InvMass;Y",100,-5,5,100,-5,5);
	TH2F *histNonPromptPsi2sPtVersusY  = new TH2F("histNonPromptPsi2sPtVersusY","p_{T} vs Y;p_{T};Y",100,0,20,100,-5,5);
	TH2F *histNonPromptPsi2sPtPairVersusY  = new TH2F("histnONPromptPsi2sPtPairVersusY","p_{T} Pairs vs Y;p_{T} Pairs;Y",100,0,20,100,-5,5);
	TH2F *histNonPromptPsi2sPtPairVersusPt  = new TH2F("histNonPromptPsi2sPtPairVersusPt","p_{T} Pairs vs p_{T};p_{T} Pairs;p_{T}",100,0,20,100,0,20);

	// Non Prompt J/psi and psi(2s) after cuts
	TH2F *histNonPromptJpsiPhiVersusEtaAfterCuts = new TH2F("histNonPromptJpsiPhiVersusEtaAfterCuts","#varphi vs #eta;#varphi;#eta",100,-5,5,100,-5,5);
	TH2F *histNonPromptJpsiEtaVersusPtAfterCuts = new TH2F("histNonPromptJpsiEtaVersusPtAfterCuts","#eta vs p_{T};#eta;p_{T}",100,-5,5,100,0,20);
	TH2F *histNonPromptJpsiEtaVersusPtPairAfterCuts = new TH2F("histNonPromptJpsiEtaVersusPtPairAfterCuts","#eta v. p_{T} Pairs;#eta;p_{T} Pairs",100,-5,5,100,0,20);
	TH2F *histNonPromptJpsiMassVersusPtAfterCuts  = new TH2F("histNonPromptJpsiMassVersusPtAfterCuts","InvMass vs p_{T};InvMass;p_{T}",100,-5,5,100,0,20);
	TH2F *histNonPromptJpsiMassVersusPtPairAfterCuts  = new TH2F("histNonPromptJpsiMassVersusPtPairAfterCuts","InvMass vs p_{T} Pairs;InvMass;p_{T} Pairs",100,0,5,100,0,20);
	TH2F *histNonPromptJpsiMassVersusYAfterCuts  = new TH2F("histNonPromptJpsiMassVersusYAfterCuts","InvMass vs Y;InvMass;Y",100,-5,5,100,-5,5);
	TH2F *histNonPromptJpsiPtVersusYAfterCuts  = new TH2F("histNonPromptJpsiPtVersusYAfterCuts","p_{T} vs Y;p_{T};Y",100,0,20,100,-5,5);
	TH2F *histNonPromptJpsiPtPairVersusYAfterCuts  = new TH2F("histNonPromptJpsiPtPairVersusYAfterCuts","p_{T} Pairs vs Y;p_{T} Pairs;Y",100,0,20,100,-5,5);
	TH2F *histNonPromptJpsiPtPairVersusPtAfterCuts  = new TH2F("histNonPromptJpsiPtPairVersusPtAfterCuts","p_{T} Pairs vs p_{T};p_{T} Pairs;p_{T}",100,0,20,100,0,20);

	TH2F *histNonPromptPsi2sPhiVersusEtaAfterCuts = new TH2F("histNonPromptPsi2sPhiVersusEtaAfterCuts","#varphi vs #eta;#varphi;#eta",100,-5,5,100,-5,5);
	TH2F *histNonPromptPsi2sEtaVersusPtAfterCuts = new TH2F("histNonPromptPsi2sEtaVersusPtAfterCuts","#eta vs p_{T};#eta;p_{T}",100,-5,5,100,0,20);
	TH2F *histNonPromptPsi2sEtaVersusPtPairAfterCuts = new TH2F("histNonPromptPsi2sEtaVersusPtPairAfterCuts","#eta v. p_{T} Pairs;#eta;p_{T} Pairs",100,-5,5,100,0,20);
	TH2F *histNonPromptPsi2sMassVersusPtAfterCuts  = new TH2F("histNonPromptPsi2sMassVersusPtAfterCuts","InvMass vs p_{T};InvMass;p_{T}",100,0,5,100,0,20);
	TH2F *histNonPromptPsi2sMassVersusPtPairAfterCuts  = new TH2F("histNonPromptPsi2sMassVersusPtPairAfterCuts","InvMass vs p_{T} Pairs;InvMass;p_{T} Pairs",100,-0,5,100,0,20);
	TH2F *histNonPromptPsi2sMassVersusYAfterCuts  = new TH2F("histNonPromptPsi2sMassVersusYAfterCuts","InvMass vs Y;InvMass;Y",100,-5,5,100,-5,5);
	TH2F *histNonPromptPsi2sPtVersusYAfterCuts  = new TH2F("histNonPromptPsi2sPtVersusYAfterCuts","p_{T} vs Y;p_{T};Y",100,0,20,100,-5,5);
	TH2F *histNonPromptPsi2sPtPairVersusYAfterCuts  = new TH2F("histNonPromptPsi2sPtPairVersusYAfterCuts","p_{T} Pairs vs Y;p_{T} Pairs;Y",100,0,20,100,-5,5);
	TH2F *histNonPromptPsi2sPtPairVersusPtAfterCuts  = new TH2F("histNonPromptPsi2sPtPairVersusPtAfterCuts","p_{T} Pairs vs p_{T};p_{T} Pairs;p_{T}",100,0,20,100,0,20);

	/// Special Corelations for Detector Analysis

	/////////////////////////////////////////////////
	/// Histograms For Checking hits on Detectors ///
	/////////////////////////////////////////////////

	TH1F *histTPCHits = new TH1F("histTPCHits","TPC Hits;Hits;Entries",kBinRange_250,0,1000);
	TH1F *histITSHits = new TH1F("histITSHits","ITS Hits;Hits;Entries",kBinRange_250,0,1000);
	TH2F *histTPCVersusITS = new TH2F("histTPCVersusITS","TPC Hits vs ITS Hits;TPC Hits;ITS Hits",250,0,1000,250,0,1000);

	//////////////////////////////////////////////////
	/// Histograms For h_{B} Analysis pseudoproper ///
	//////////////////////////////////////////////////

  //pseudoproper decay length

	TH1F *histPseudoProperDecayLengthB0 = new TH1F("histPseudoProperDecayLengthB0",";Pseudo-proper decay length;Entries",kBinRange_100,-0.2,0.2);
	TH1F *histPseudoProperDecayLengthBplus = new TH1F("histPseudoProperDecayLengthBplus",";Pseudo-proper decay length;Entries",kBinRange_100,-2,2);
	TH1F *histPseudoProperDecayLengthBszero = new TH1F("histPseudoProperDecayLengthBszero",";Pseudo-proper decay length;Entries",kBinRange_100,-0.2,0.2);
	TH1F *histPseudoProperDecayLengthBcplus = new TH1F("histPseudoProperDecayLengthBcplus",";Pseudo-proper decay length;Entries",kBinRange_100,-0.2,0.2);
	TH1F *histPseudoProperDecayLengthSigmaB = new TH1F("histPseudoProperDecayLengthSigmaB ",";Pseudo-proper decay length;Entries",kBinRange_100,-0.2,0.2);
	TH1F *histPseudoProperDecayLengthLambdaBzero = new TH1F("histPseudoProperDecayLengthLambdaBzero",";Pseudo-proper decay length;Entries",kBinRange_100,-0.2,0.2);
	TH1F *histPseudoProperDecayLengthXiBminus = new TH1F("histPseudoProperDecayLengthXiBminus",";Pseudo-proper decay length;Entries",kBinRange_100,-0.2,0.2);
	TH1F *histPseudoProperDecayLengthXiBzero = new TH1F("histPseudoProperDecayLengthXiBzero",";Pseudo-proper decay length;Entries",kBinRange_100,-0.2,0.2);
	TH1F *histPseudoProperDecayLengthOmegaBminus = new TH1F("histPseudoProperDecayLengthOmegaBminus",";Pseudo-proper decay length;Entries",kBinRange_100,-0.2,0.2);

	//pseudoproper decaytime
	TH1F *histPseudoProperDecayTimeB0 = new TH1F("histPseudoProperDecayTimeB0",";Pseudo-proper decay time;Entries",kBinRange_100,-0.2,0.2);
	TH1F *histPseudoProperDecayTimeBplus = new TH1F("histPseudoProperDecayTimeBplus",";Pseudo-proper decay time;Entries",kBinRange_100,-0.2,0.2);


	///////////////////////////////////////////////////////
	/// Histograms For Quarkonium Analysis pseudoproper ///
	///////////////////////////////////////////////////////

	//pseudoproper decay length

	TH1F *histPseudoProperDecayLengthJpsi = new TH1F("histPseudoProperDecayLengthJpsi",";p_{T}(GeV/c);Entries",kBinRange_100,-0.2,0.2);

	//pseudoproper decaytime

	TH1F *histPseudoProperDecayTimeJpsi = new TH1F("histPseudoProperDecayTimeJpsi",";p_{T}(GeV/c);Entries",kBinRange_100,-0.2,0.2);


	////////////////////////////////////////////////////
	// Error bars for Efficiency and pT Distributions //
	////////////////////////////////////////////////////

	// For Efficiencies
	histPtEfficiencyJpsiPrompt->Sumw2();
	histPtPairEfficiencyJpsiPrompt->Sumw2();
	histPtEfficiencyPsi2sPrompt->Sumw2();
	histPtPairEfficiencyPsi2sPrompt->Sumw2();
	histYEfficiencyJpsiPrompt->Sumw2();
	histYEfficiencyPsi2sPrompt->Sumw2();

	histPtEfficiencyJpsiNonPrompt->Sumw2();
	histPtPairEfficiencyJpsiNonPrompt->Sumw2();
	histPtEfficiencyPsi2sNonPrompt->Sumw2();
	histPtPairEfficiencyPsi2sNonPrompt->Sumw2();
	histYEfficiencyJpsiNonPrompt->Sumw2();
	histYEfficiencyPsi2sNonPrompt->Sumw2();

	// For pT and pT Pair Distributions
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

	// For Eta Distributions
	histEtaJpsiPrompt->Sumw2();
	histEtaPsi2sPrompt->Sumw2();
	histEtaJpsiPromptAfterCuts->Sumw2();
	histEtaPsi2sPromptAfterCuts->Sumw2();

	histEtaJpsiNonPrompt->Sumw2();
	histEtaPsi2sNonPrompt->Sumw2();
	histEtaJpsiNonPromptAfterCuts->Sumw2();
	histEtaPsi2sNonPromptAfterCuts->Sumw2();

	// For Phi Distributions
	histPhiJpsiPrompt->Sumw2();
	histPhiPsi2sPrompt->Sumw2();
	histPhiJpsiPromptAfterCuts->Sumw2();
	histPhiPsi2sPromptAfterCuts->Sumw2();

	histPhiJpsiNonPrompt->Sumw2();
	histPhiPsi2sNonPrompt->Sumw2();
	histPhiJpsiNonPromptAfterCuts->Sumw2();
	histPhiPsi2sNonPromptAfterCuts->Sumw2();

	// For Theta Distrubutions
	histThetaJpsiPrompt->Sumw2();
	histThetaPsi2sPrompt->Sumw2();
	histThetaJpsiPromptAfterCuts->Sumw2();
	histThetaPsi2sPromptAfterCuts->Sumw2();

	histThetaJpsiNonPrompt->Sumw2();
	histThetaPsi2sNonPrompt->Sumw2();
	histThetaJpsiNonPromptAfterCuts->Sumw2();
	histThetaPsi2sNonPromptAfterCuts->Sumw2();

	///////////////////////////////////////////////////////////////////
	/// Histogram Definitions for Skimmed Data and Machine Learning ///
	///////////////////////////////////////////////////////////////////

	/// Define Branches Before Cuts

	TTree *PromptJpsiTree = new TTree("PromptJpsiTree", "TTree with a structure before cuts");
	PromptJpsiTree->Branch("kInvMassPair", &kInvMassPair, "kInvMassPair/F");
	PromptJpsiTree->Branch("kPt", &kPt, "kPt/F");
	PromptJpsiTree->Branch("kPtPair", &kPtPair, "kPtPair/F");
	PromptJpsiTree->Branch("kY", &kY, "kY/F");
	PromptJpsiTree->Branch("kEta", &kEta, "kEta/F");
	PromptJpsiTree->Branch("kTheta", &kTheta, "kTheta/F");
	PromptJpsiTree->Branch("kPhi", &kPhi, "kPhi/F");

	TTree *PromptPsi2STree = new TTree("PromptPsi2STree", "TTree with a structure before cuts");
	PromptPsi2STree->Branch("kInvMassPair", &kInvMassPair, "kInvMassPair/F");
	PromptPsi2STree->Branch("kPt", &kPt, "kPt/F");
	PromptPsi2STree->Branch("kPtPair", &kPtPair, "kPtPair/F");
	PromptPsi2STree->Branch("kY", &kY, "kY/F");
	PromptPsi2STree->Branch("kEta", &kEta, "kEta/F");
	PromptPsi2STree->Branch("kTheta", &kTheta, "kTheta/F");
	PromptPsi2STree->Branch("kPhi", &kPhi, "kPhi/F");

	TTree *NonPromptJpsiTree = new TTree("NonPromptJpsiTree", "TTree with a structure before cuts");
	NonPromptJpsiTree->Branch("kInvMassPair", &kInvMassPair, "kInvMassPair/F");
	NonPromptJpsiTree->Branch("kPt", &kPt, "kPt/F");
	NonPromptJpsiTree->Branch("kPtPair", &kPtPair, "kPtPair/F");
	NonPromptJpsiTree->Branch("kY", &kY, "kY/F");
	NonPromptJpsiTree->Branch("kEta", &kEta, "kEta/F");
	NonPromptJpsiTree->Branch("kTheta", &kTheta, "kTheta/F");
	NonPromptJpsiTree->Branch("kPhi", &kPhi, "kPhi/F");

	TTree *NonPromptPsi2STree = new TTree("NonPromptPsi2STree", "TTree with a structure before cuts");
	NonPromptPsi2STree->Branch("kInvMassPair", &kInvMassPair, "kInvMassPair/F");
	NonPromptPsi2STree->Branch("kPt", &kPt, "kPt/F");
	NonPromptPsi2STree->Branch("kPtPair", &kPtPair, "kPtPair/F");
	NonPromptPsi2STree->Branch("kY", &kY, "kY/F");
	NonPromptPsi2STree->Branch("kEta", &kEta, "kEta/F");
	NonPromptPsi2STree->Branch("kTheta", &kTheta, "kTheta/F");
	NonPromptPsi2STree->Branch("kPhi", &kPhi, "kPhi/F");

	// Define Branches After Cuts
	TTree *PromptJpsiTreeAfterCuts = new TTree("PromptJpsiTreeAfterCuts", "TTree with a structure before cuts");
	PromptJpsiTreeAfterCuts->Branch("kInvMassPair", &kInvMassPairBranchCut, "kInvMassPair/F");
	PromptJpsiTreeAfterCuts->Branch("kPt", &kPtBranchCut, "kPt/F");
	PromptJpsiTreeAfterCuts->Branch("kPtPair", &kPtPairBranchCut, "kPtPair/F");
	PromptJpsiTreeAfterCuts->Branch("kY", &kYBranchCut, "kY/F");
	PromptJpsiTreeAfterCuts->Branch("kEta", &kEtaBranchCut, "kEta/F");
	PromptJpsiTreeAfterCuts->Branch("kTheta", &kThetaBranchCut, "kTheta/F");
	PromptJpsiTreeAfterCuts->Branch("kPhi", &kPhiBranchCut, "kPhi/F");

	TTree *PromptPsi2STreeAfterCuts = new TTree("PromptPsi2STreeAfterCuts", "TTree with a structure before cuts");
	PromptPsi2STreeAfterCuts->Branch("kInvMassPair", &kInvMassPairBranchCut, "kInvMassPair/F");
	PromptPsi2STreeAfterCuts->Branch("kPt", &kPtBranchCut, "kPt/F");
	PromptPsi2STreeAfterCuts->Branch("kPtPair", &kPtPairBranchCut, "kPtPair/F");
	PromptPsi2STreeAfterCuts->Branch("kY", &kYBranchCut, "kY/F");
	PromptPsi2STreeAfterCuts->Branch("kEta", &kEtaBranchCut, "kEta/F");
	PromptPsi2STreeAfterCuts->Branch("kTheta", &kThetaBranchCut, "kTheta/F");
	PromptPsi2STreeAfterCuts->Branch("kPhi", &kPhiBranchCut, "kPhi/F");

	TTree *NonPromptJpsiTreeAfterCuts = new TTree("NonPromptJpsiTreeAfterCuts", "TTree with a structure before cuts");
	NonPromptJpsiTreeAfterCuts->Branch("kInvMassPair", &kInvMassPairBranchCut, "kInvMassPair/F");
	NonPromptJpsiTreeAfterCuts->Branch("kPt", &kPtBranchCut, "kPt/F");
	NonPromptJpsiTreeAfterCuts->Branch("kPtPair", &kPtPairBranchCut, "kPtPair/F");
	NonPromptJpsiTreeAfterCuts->Branch("kY", &kYBranchCut, "kY/F");
	NonPromptJpsiTreeAfterCuts->Branch("kEta", &kEtaBranchCut, "kEta/F");
	NonPromptJpsiTreeAfterCuts->Branch("kTheta", &kThetaBranchCut, "kTheta/F");
	NonPromptJpsiTreeAfterCuts->Branch("kPhi", &kPhiBranchCut, "kPhi/F");

	TTree *NonPromptPsi2STreeAfterCuts = new TTree("NonPromptPsi2STreeAfterCuts", "TTree with a structure before cuts");
	NonPromptPsi2STreeAfterCuts->Branch("kInvMassPair", &kInvMassPairBranchCut, "kInvMassPair/F");
	NonPromptPsi2STreeAfterCuts->Branch("kPt", &kPtBranchCut, "kPt/F");
	NonPromptPsi2STreeAfterCuts->Branch("kPtPair", &kPtPairBranchCut, "kPtPair/F");
	NonPromptPsi2STreeAfterCuts->Branch("kY", &kYBranchCut, "kY/F");
	NonPromptPsi2STreeAfterCuts->Branch("kEta", &kEtaBranchCut, "kEta/F");
	NonPromptPsi2STreeAfterCuts->Branch("kTheta", &kThetaBranchCut, "kTheta/F");
	NonPromptPsi2STreeAfterCuts->Branch("kPhi", &kPhiBranchCut, "kPhi/F");

	/// J/Psi and Psi(2S) Pdg Codes
	Int_t kMotherPdgCode[] = {443, 100443};
	Double_t kCandidateJpsiMass = TDatabasePDG::Instance()->GetParticle(kMotherPdgCode[0])->Mass();
	Double_t kCandidatePsi2SMass = TDatabasePDG::Instance()->GetParticle(kMotherPdgCode[1])->Mass();

	/// Dilepton Pairs Pdg Codes
	Int_t kPairsPdgCode[] = {11,-11,13,-13};

	Int_t kPairElectronPdg  = kPairsPdgCode[0];
	Int_t kPairPositronPdg  = kPairsPdgCode[1];
	Int_t kPairMuonPdg      = kPairsPdgCode[2];
	Int_t kPairAntiMuonPdg  = kPairsPdgCode[3];

	Bool_t isPairElectron;
	Bool_t isPairPositron;
	Bool_t isPairMuon;
	Bool_t isPairAntiMuon;

	// Prefix for signal or background
	//	const char* nameprefix = "sgn";
	const char* nameprefix = "sgn_1";

	FairLogger::GetLogger()->SetLogScreenLevel("DEBUG");
	TFile f(o2::base::NameConf::getMCKinematicsFileName(nameprefix).c_str());

	LOG(debug) << "Checking input file :" << f.GetPath();

	std::vector<o2::MCTrack>* mctracks = nullptr;
	auto getTTree_o2sim = (TTree*)f.Get("o2sim");
	if(!getTTree_o2sim) printf("Warning!!! File not Found. \n");
	assert(getTTree_o2sim);

	auto getBranchMCTrack = getTTree_o2sim->GetBranch("MCTrack");
	assert(getBranchMCTrack);
	getBranchMCTrack->SetAddress(&mctracks);

	std::vector<o2::TrackReference>* trackrefs = nullptr;
	auto getBranchTrackRefs = getTTree_o2sim->GetBranch("TrackRefs");
	assert(getBranchTrackRefs);
	getBranchTrackRefs->SetAddress(&trackrefs);

	o2::steer::MCKinematicsReader mcreader(nameprefix, o2::steer::MCKinematicsReader::Mode::kMCKine);

	// when present we also read some hits for ITS to test consistency of trackID assignments
	TFile hitf(o2::base::DetectorNameConf::getHitsFileName(o2::detectors::DetID::ITS, nameprefix).c_str());
	auto hittr = (TTree*)hitf.Get("o2sim");
	auto hitbr = hittr ? hittr->GetBranch("ITSHit") : nullptr;
	std::vector<o2::itsmft::Hit>* hits = nullptr;
	if (hitbr) {
		hitbr->SetAddress(&hits);
	}

	for (int eventID = 0; eventID < getBranchMCTrack->GetEntries(); ++eventID) {
		getBranchMCTrack->GetEntry(eventID);
		getBranchTrackRefs->GetEntry(eventID);
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

		Int_t kBeautyPdgCode[] = {511, 521, 531, 541, 5112, 5122, 5132, 5232, 5332};
		Int_t kSizeBeautyPdgCode = sizeof(kBeautyPdgCode)/sizeof(Int_t);

		/// Variables for Count histBeauties (PDG CODES)
		Int_t FromBeauty;
		Int_t kCandidateBzeroPdg 				= kBeautyPdgCode[0];
		Int_t kCandidateBplusPdg 				= kBeautyPdgCode[1];
		Int_t kCandidateBszeroPdg 			= kBeautyPdgCode[2];
		Int_t kCandidateBcplusPdg 			= kBeautyPdgCode[3];
		Int_t kCandidateSigmaBPdg 			= kBeautyPdgCode[4];
		Int_t kCandidateLambdaBzeroPdg 	= kBeautyPdgCode[5];
		Int_t kCandidateXiBminusPdg 		= kBeautyPdgCode[6];
		Int_t kCandidateXiBzeroPdg 			= kBeautyPdgCode[7];
		Int_t kCandidateOmegaBminusPdg 	= kBeautyPdgCode[8];

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
			Int_t kFirstDaughterTrackID = t.getFirstDaughterTrackId();
			Int_t kLastDaughterTrackID = t.getLastDaughterTrackId();
			Bool_t hasBeautyMother = kFALSE;
			Bool_t isJpsi =  TMath::Abs(t.GetPdgCode()) == kMotherPdgCode[0];
			Bool_t isPsi2S =  TMath::Abs(t.GetPdgCode()) == kMotherPdgCode[1];
			if( isJpsi || isPsi2S ) {
				// Daughters
				o2::MCTrack *kPair=0x0; o2::MCTrack *kAntiPair=0x0;
				// Charmonium
				o2::MCTrack *kCandidateJpsi=0x0; o2::MCTrack *kCandidatePsi2S=0x0;
				// Beauties
				o2::MCTrack *kCandidateBzero=0x0;
				o2::MCTrack *kCandidateBplus=0x0;
				o2::MCTrack *kCandidateBszero=0x0;
				o2::MCTrack *kCandidateBcplus=0x0;
				o2::MCTrack *kCandidateSigmaB=0x0;
				o2::MCTrack *kCandidateLambdaBzero=0x0;
				o2::MCTrack *kCandidateXiBminus=0x0;
				o2::MCTrack *kCandidateXiBzero=0x0;
				o2::MCTrack *kCandidateOmegaBminus=0x0;

				////////////////////////////////////////////////////////////////////////
				// CONSTRUCTOR DENEME
				//o2::MCTrackT(Int_t pdgCode = 1, Int_t motherID, Int_t secondMotherID, Int_t firstDaughterID, Int_t lastDaughterID,
        //   Double_t px, Double_t py, Double_t pz, Double_t x, Double_t y, Double_t z, Double_t t,
        //   Int_t nPoints);
				//o2::MCTrackT myObj = default;

				//o2::MCTrackT(pdgCode, motherID, secondMotherID, firstDaughterID, lastDaughterID,
				 //px,  py, pz,  x,  y,  z,  t,
				//	nPoints);
				///  Standard constructor
				//o2::MCTrackT(pdgCode, motherID, secondMotherID, firstDaughterID, lastDaughterID,
				// px,  py, pz,  x,  y,  z,  t,
				  //nPoints);
				//o2::MCTrackT();
				//o2::MCTrackT();
				////////////////////////////////////////////////////////////////////////

				Int_t kMotherTrackID = t.getMotherTrackId();
				Int_t kSecondMotherTrackID = t.getSecondMotherTrackId();
				Bool_t isPrompt = kMotherTrackID <= 0 ? kTRUE : kFALSE;
				auto tdPromptQuarkonium = mcreader.getTrack(eventID, kMotherTrackID);
				auto td1 =	mcreader.getTrack(eventID, kFirstDaughterTrackID);
				auto td2 =	mcreader.getTrack(eventID, kLastDaughterTrackID);
				auto td3 = 	mcreader.getTrack(eventID, kSecondMotherTrackID);
				std::cout << "tdPromptQuarkonium için PDG: " << TMath::Abs(tdPromptQuarkonium->GetPdgCode()) << '\n';
				std::cout << "td1 için PDG: " << TMath::Abs(td1->GetPdgCode()) << '\n';
				std::cout << "td2 için PDG: " << TMath::Abs(td2->GetPdgCode()) << '\n';
				std::cout << "td3 için PDG: " << TMath::Abs(td3->GetPdgCode()) << '\n';


				if(isPrompt) {
					 //auto tdPromptQuarkonium = mcreader.getTrack(eventID, kMotherTrackID);
					 //std::cout << "tdPromptQuarkonium için PDG: " << TMath::Abs(tdPromptQuarkonium->GetPdgCode()) << '\n';
					 //if( TMath::Abs(t.GetPdgCode()) == kMotherPdgCode[0] ) kCandidateJpsi = (o2::MCTrack*)tdPromptQuarkonium;
					 //if( TMath::Abs(t.GetPdgCode()) == kMotherPdgCode[1] ) kCandidatePsi2S = (o2::MCTrack*)tdPromptQuarkonium;
				 }
				if(!isPrompt) { //  check beauty mother
					auto tdM = mcreader.getTrack(eventID, kMotherTrackID);
					// DENEME
					//auto tdM2 = mcreader.getTrack(eventID, kSecondMotherTrackID);
					//if (TMath::Abs(tdM->GetPdgCode()) == kMotherPdgCode[0] ) kCandidateJpsi  = (o2::MCTrack*)tdM2; std::cout << "zaa" << '\n';
					//if (TMath::Abs(tdM->GetPdgCode()) == kMotherPdgCode[1] ) kCandidatePsi2S = (o2::MCTrack*)tdM2;
					//if( TMath::Abs(t.GetPdgCode()) == kMotherPdgCode[0] ) kCandidateJpsi = (o2::MCTrack*)tdM;; // NonPrompt Jpsi
					//if( TMath::Abs(t.GetPdgCode()) == kMotherPdgCode[1] ) kCandidatePsi2S = (o2::MCTrack*)tdM; //NonPrompt Psi2s
					for(int i=0; i<kSizeBeautyPdgCode; i++) {
						if (TMath::Abs(tdM->GetPdgCode()) == kBeautyPdgCode[i] ){
							 hasBeautyMother = kTRUE;
							 //std::cout << "tdM için PDG: " << TMath::Abs(tdM->GetPdgCode()) << '\n';
						if (TMath::Abs(tdM->GetPdgCode()) == kBeautyPdgCode[0] ) kCandidateBzero = (o2::MCTrack*)tdM;
						if (TMath::Abs(tdM->GetPdgCode()) == kBeautyPdgCode[1] ) kCandidateBplus = (o2::MCTrack*)tdM;
						if (TMath::Abs(tdM->GetPdgCode()) == kBeautyPdgCode[2] ) kCandidateBszero = (o2::MCTrack*)tdM;
						if (TMath::Abs(tdM->GetPdgCode()) == kBeautyPdgCode[3] ) kCandidateBcplus = (o2::MCTrack*)tdM;
						if (TMath::Abs(tdM->GetPdgCode()) == kBeautyPdgCode[4] ) kCandidateSigmaB = (o2::MCTrack*)tdM;
						if (TMath::Abs(tdM->GetPdgCode()) == kBeautyPdgCode[5] ) kCandidateLambdaBzero = (o2::MCTrack*)tdM;
						if (TMath::Abs(tdM->GetPdgCode()) == kBeautyPdgCode[6] ) kCandidateXiBminus = (o2::MCTrack*)tdM;
						if (TMath::Abs(tdM->GetPdgCode()) == kBeautyPdgCode[7] ) kCandidateXiBzero = (o2::MCTrack*)tdM;
						if (TMath::Abs(tdM->GetPdgCode()) == kBeautyPdgCode[8] ) kCandidateOmegaBminus = (o2::MCTrack*)tdM;

					  }   // if (TMath::Abs(tdM->GetPdgCode()) == kBeautyPdgCode[i] )
					} // for(int i=0; i<kSizeBeautyPdgCode; i++)
				}  // if(!isPrompt)

				/// Get Mother Tracks with Getter Method
				kPt = t.GetPt();
				kP = t.GetP();
				kPhi = t.GetPhi();
				kTheta = t.GetTheta();
				kY = t.GetRapidity();
				kEta = t.GetEta();
				kPrimaryVx = t.Vx();
				kPrimaryVy = t.Vy();
				kPrimaryVz = t.Vz();
				kMass = t.GetMass();
				std::cout << "Track Mass: " << kMass << '\n';
				/// Checking Daughters from Primary Vertex
				kFirstDaughterTrackID  = (mcreader.getTrack(eventID, kFirstDaughterTrackID))->isPrimary() ? kFirstDaughterTrackID : -1;
				kLastDaughterTrackID   = (mcreader.getTrack(eventID, kLastDaughterTrackID))->isPrimary() ? kLastDaughterTrackID : -1;

				/// Checking Mothers From Primary Vertex (SONRA EKLENDİ CHECK ET)
				//kMotherTrackID = (mcreader.getTrack(eventID, kMotherTrackID))->isPrimary() ? kMotherTrackID : -1;
				//if(kMotherTrackID->GetPdgCode() == kMotherPdgCode[0]) kCandidateJpsi = (o2::MCTrack*)kMotherTrackID;
				//LOG(DEBUG) << " mother track - pdg " << t.GetPdgCode() << " first daughter "<<  kFirstDaughterTrackID <<" last daughter " << kLastDaughterTrackID << " position " << ti;

				for(int idaugh=kFirstDaughterTrackID; idaugh<kLastDaughterTrackID+1; idaugh++ ) {
					// kDaughterTrack track daughter demek. t ise motherları temsil eder.
					auto kDaughterTrack = mcreader.getTrack(eventID, idaugh);
					if(kDaughterTrack->GetPdgCode() == kSelectionPdg) 		kPair = (o2::MCTrack*)kDaughterTrack;
					if(kDaughterTrack->GetPdgCode() == -1*kSelectionPdg) kAntiPair = (o2::MCTrack*)kDaughterTrack;
				}

				if((!kPair) || (!kAntiPair)) continue;

				///////////////////////
				/// Cut Definitions ///
				///////////////////////

				// TODO: Bunu cut definition için optimize et.
				/*
				Bool_t TrackBarrelCut = TMath::Abs(kPair->GetEta()) < 0.9 && TMath::Abs(kAntiPair->GetEta()) < 0.9 &&
						                    TMath::Abs(kPair->GetRapidity()) < 0.9 && TMath::Abs(kAntiPair->GetRapidity()) < 0.9 &&
																kPair->GetPt() > 1.0 && kAntiPair->GetPt() > 1.0;

				Bool_t ForwardYCut 		=	kPair->GetEta() < -2.5 && kAntiPair->GetEta() < -2.5 && kPair->GetEta() > -4.0 && kAntiPair->GetEta() > -4.0 &&
																kPair->GetRapidity() < -2.5 && kAntiPair->GetRapidity() < -2.5 && kPair->GetRapidity() > -4.0 && kAntiPair->GetRapidity() > -4.0 &&
																kPair->GetPt() > 1.0 && kAntiPair->GetPt() > 1.0;
																*/

				//if(kMidYSimulation) TrackBarrelCut == kTRUE; ForwardYCut = kFALSE;
				//if(kForwardYSimulation) TrackBarrelCut == kFALSE; ForwardYCut = kTRUE;

				//TrackBarrelCut ? -->

				// if(!TrackBarrelCut) continue;
				// if(!kForwardYSimulation) continue;

				// Check the Pairs
				isPairElectron = kPair->GetPdgCode()     == kPairElectronPdg;
				isPairPositron = kAntiPair->GetPdgCode() == kPairPositronPdg;
				isPairMuon     = kPair->GetPdgCode() 		 == kPairMuonPdg;
				isPairAntiMuon = kAntiPair->GetPdgCode() == kPairAntiMuonPdg;

				/////

				//Double_t vertexDistance = vertexB0->DistanceToVertex(primaryVertex);
				//Double_t normDecayLength = kCandidateBzero->NormalizedDecayLength();
				//Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1 / ((pdgMassMotherBzero * pdgMassMotherBzero / (momentumMotherBzero * momentumMotherBzero)) + 1)))
				//Double_t pseudoProperDecayLength = ((vertexB0->GetX() - primaryVertex->GetX()) * kCandidateBzero->Px() / TMath::Abs(kCandidateBzero->Pt())) + ((vertexB0->GetY() - primaryVertex->GetY()) * kCandidateBzero->Py() / TMath::Abs(kCandidateBzero->Pt()));

				//////

				////////////////////////////////////////////
				/// Beauty Vertex and Decay Calculations ///
				////////////////////////////////////////////

				if(kCandidateBzero){
					// Get Secondary Vertex
					//std::cout << "B zero Mass: " << kCandidateBzero->GetMass() << '\n';
					if(kCandidateBzero->isSecondary()) kCandidateBzeroVx = kCandidateBzero->Vx(); kCandidateBzeroVy = kCandidateBzero->Vy();
					Double_t pdgMassMotherBzero = TDatabasePDG::Instance()->GetParticle(kBeautyPdgCode[0])->Mass();
					Double_t momentumMotherBzero = kCandidateBzero->GetP();
					Double_t ptMotherBzero = kCandidateBzero->GetPt();
					Double_t pseudoProperDecayLengthB0 = ((kCandidateBzeroVx - kPrimaryVx) * kCandidateBzero->Px() / TMath::Abs(kCandidateBzero->GetPt())) + ((kCandidateBzeroVy - kPrimaryVy) * kCandidateBzero->Py() / TMath::Abs(kCandidateBzero->GetPt()));
					Double_t pseudoProperDecayTimeB0 = pseudoProperDecayLengthB0 * pdgMassMotherBzero / ptMotherBzero;

					histPseudoProperDecayLengthB0->Fill(pseudoProperDecayLengthB0);
					histPseudoProperDecayTimeB0->Fill(pseudoProperDecayTimeB0);
				}
				// TODO: bug var fixle.
				if(kCandidateBplus){
					if(kCandidateBplus->isSecondary()) { kCandidateBplusVx = kCandidateBplus->Vx(); kCandidateBplusVy = kCandidateBplus->Vy(); }
					Double_t pdgMassMotherBplus = TDatabasePDG::Instance()->GetParticle(kBeautyPdgCode[1])->Mass();
					Double_t momentumMotherBplus = kCandidateBplus->GetP();
					Double_t ptMotherBplus = kCandidateBplus->GetPt();
					Double_t pseudoProperDecayLengthBplus = ((kCandidateBplusVx - kPrimaryVx) * kCandidateBplus->Px() / TMath::Abs(kCandidateBplus->GetPt())) + ((kCandidateBplusVy - kPrimaryVy) * kCandidateBplus->Py() / TMath::Abs(kCandidateBplus->GetPt()));
					Double_t pseudoProperDecayTimeBplus = pseudoProperDecayLengthBplus * pdgMassMotherBplus / ptMotherBplus;

					histPseudoProperDecayLengthBplus->Fill(pseudoProperDecayLengthBplus);
					histPseudoProperDecayTimeBplus->Fill(pseudoProperDecayTimeBplus);
				}
				if(kCandidateBszero){
					if(kCandidateBszero->isSecondary()) kCandidateBszeroVx = kCandidateBszero->Vx(); kCandidateBszeroVy = kCandidateBszero->Vy();
					Double_t pdgMassMotherBszero = TDatabasePDG::Instance()->GetParticle(kBeautyPdgCode[2])->Mass();
					Double_t momentumMotherBszero = kCandidateBszero->GetP();
					Double_t ptMotherBszero = kCandidateBszero->GetPt();
					Double_t pseudoProperDecayLengthBszero = ((kCandidateBszeroVx - kPrimaryVx) * kCandidateBszero->Px() / TMath::Abs(kCandidateBszero->GetPt())) + ((kCandidateBszeroVy - kPrimaryVy) * kCandidateBszero->Py() / TMath::Abs(kCandidateBszero->GetPt()));
					histPseudoProperDecayLengthBszero->Fill(pseudoProperDecayLengthBszero);
				}
				if(kCandidateBcplus){
					if(kCandidateBcplus->isSecondary()) kCandidateBcplusVx = kCandidateBcplus->Vx(); kCandidateBcplusVy = kCandidateBcplus->Vy();
					Double_t pdgMassMotherBcplus = TDatabasePDG::Instance()->GetParticle(kBeautyPdgCode[3])->Mass();
					Double_t momentumMotherBcplus = kCandidateBcplus->GetP();
					Double_t ptMotherBcplus = kCandidateBcplus->GetPt();
					Double_t pseudoProperDecayLengthBcplus = ((kCandidateBcplusVx - kPrimaryVx) * kCandidateBcplus->Px() / TMath::Abs(kCandidateBcplus->GetPt())) + ((kCandidateBcplusVy - kPrimaryVy) * kCandidateBcplus->Py() / TMath::Abs(kCandidateBcplus->GetPt()));
					histPseudoProperDecayLengthBcplus->Fill(pseudoProperDecayLengthBcplus);
				}
				if(kCandidateSigmaB){
					if(kCandidateSigmaB->isSecondary()) kCandidateSigmaBVx = kCandidateSigmaB->Vx(); kCandidateSigmaBVy = kCandidateSigmaB->Vy();
					Double_t pdgMassMotherSigmaB = TDatabasePDG::Instance()->GetParticle(kBeautyPdgCode[4])->Mass();
					Double_t momentumMotherSigmaB = kCandidateSigmaB->GetP();
					Double_t ptMotherSigmaB = kCandidateSigmaB->GetPt();
					Double_t pseudoProperDecayLengthSigmaB = ((kCandidateSigmaBVx - kPrimaryVx) * kCandidateSigmaB->Px() / TMath::Abs(kCandidateSigmaB->GetPt())) + ((kCandidateSigmaBVy - kPrimaryVy) * kCandidateSigmaB->Py() / TMath::Abs(kCandidateSigmaB->GetPt()));
					histPseudoProperDecayLengthSigmaB->Fill(pseudoProperDecayLengthSigmaB);
				}
				if(kCandidateLambdaBzero){
					if(kCandidateLambdaBzero->isSecondary()) kCandidateLambdaBzeroVx = kCandidateLambdaBzero->Vx(); kCandidateLambdaBzeroVy = kCandidateLambdaBzero->Vy();
					Double_t pdgMassMotherLambdaBzero = TDatabasePDG::Instance()->GetParticle(kBeautyPdgCode[5])->Mass();
					Double_t momentumMotherLambdaBzero = kCandidateLambdaBzero->GetP();
					Double_t ptMotherLambdaBzero = kCandidateLambdaBzero->GetPt();
					Double_t pseudoProperDecayLengthLambdaBzero = ((kCandidateLambdaBzeroVx - kPrimaryVx) * kCandidateLambdaBzero->Px() / TMath::Abs(kCandidateLambdaBzero->GetPt())) + ((kCandidateLambdaBzeroVy - kPrimaryVy) * kCandidateLambdaBzero->Py() / TMath::Abs(kCandidateLambdaBzero->GetPt()));
					histPseudoProperDecayLengthLambdaBzero->Fill(pseudoProperDecayLengthLambdaBzero);
				}
				if(kCandidateXiBminus){
					if(kCandidateXiBminus->isSecondary()) kCandidateXiBminusVx = kCandidateXiBminus->Vx(); kCandidateXiBminusVy = kCandidateXiBminus->Vy();
					Double_t pdgMassMotherXiBminus = TDatabasePDG::Instance()->GetParticle(kBeautyPdgCode[6])->Mass();
					Double_t momentumMotherXiBminus = kCandidateXiBminus->GetP();
					Double_t ptMotherXiBminus = kCandidateXiBminus->GetPt();
					Double_t pseudoProperDecayLengthXiBminus = ((kCandidateXiBminusVx - kPrimaryVx) * kCandidateXiBminus->Px() / TMath::Abs(kCandidateXiBminus->GetPt())) + ((kCandidateXiBminusVy - kPrimaryVy) * kCandidateXiBminus->Py() / TMath::Abs(kCandidateXiBminus->GetPt()));
					histPseudoProperDecayLengthXiBminus->Fill(pseudoProperDecayLengthXiBminus);
				}

				if(kCandidateXiBzero){
					//std::cout << "XiBzero Mass: " << kCandidateJpsi->GetMass() << '\n';
					if(kCandidateXiBzero->isSecondary()) kCandidateXiBzeroVx = kCandidateXiBzero->Vx(); kCandidateXiBzeroVy = kCandidateXiBzero->Vy();
					Double_t pdgMassMotherXiBzero = TDatabasePDG::Instance()->GetParticle(kBeautyPdgCode[7])->Mass();
					Double_t momentumMotherXiBzero = kCandidateXiBzero->GetP();
					Double_t ptMotherXiBzero = kCandidateXiBzero->GetPt();
					Double_t pseudoProperDecayLengthXiBzero = ((kCandidateXiBzeroVx - kPrimaryVx) * kCandidateXiBzero->Px() / TMath::Abs(kCandidateXiBzero->GetPt())) + ((kCandidateXiBzeroVy - kPrimaryVy) * kCandidateXiBzero->Py() / TMath::Abs(kCandidateXiBzero->GetPt()));
					histPseudoProperDecayLengthXiBzero->Fill(pseudoProperDecayLengthXiBzero);
				}
				if(kCandidateOmegaBminus){
					//std::cout << "OmegaBminus Mass: " << kCandidateJpsi->GetMass() << '\n';
					if(kCandidateOmegaBminus->isSecondary()) kCandidateOmegaBminusVx = kCandidateOmegaBminus->Vx(); kCandidateOmegaBminusVy = kCandidateOmegaBminus->Vy();
					Double_t pdgMassMotherOmegaBminus = TDatabasePDG::Instance()->GetParticle(kBeautyPdgCode[8])->Mass();
					Double_t momentumMotherOmegaBminus = kCandidateOmegaBminus->GetP();
					Double_t ptMotherOmegaBminus = kCandidateOmegaBminus->GetPt();
					Double_t pseudoProperDecayLengthOmegaBminus = ((kCandidateOmegaBminusVx - kPrimaryVx) * kCandidateOmegaBminus->Px() / TMath::Abs(kCandidateOmegaBminus->GetPt())) + ((kCandidateOmegaBminusVy - kPrimaryVy) * kCandidateOmegaBminus->Py() / TMath::Abs(kCandidateOmegaBminus->GetPt()));
					histPseudoProperDecayLengthOmegaBminus->Fill(pseudoProperDecayLengthOmegaBminus);
				}

				//////////////////////////////////////////////////////////
				/// Quarkonium Vertex and Decay Calculations 					 ///
				/////////////////////////////////////////////////////////
				/*
				if(kCandidateJpsi){
					std::cout << "J/psi Mass: " << kCandidateJpsi->GetMass() << '\n';
					// Get Secondary Vertex
					if(kCandidateJpsi->isSecondary()) kCandidateJpsiVx = kCandidateJpsi->Vx(); kCandidateJpsiVy = kCandidateJpsi->Vy();
					Double_t pdgMassMotherJpsi = TDatabasePDG::Instance()->GetParticle(443)->Mass();
					Double_t momentumMotherJpsi = kCandidateJpsi->GetP();
					Double_t ptMotherJpsi = kCandidateJpsi->GetPt();
					Double_t pseudoProperDecayLengthJpsi = ((kCandidateJpsiVx - kPrimaryVx) * kCandidateJpsi->Px() / TMath::Abs(kCandidateJpsi->GetPt())) + ((kCandidateJpsiVy - kPrimaryVy) * kCandidateJpsi->Py() / TMath::Abs(kCandidateJpsi->GetPt()));
					Double_t pseudoProperDecayTimeJpsi = pseudoProperDecayLengthJpsi * pdgMassMotherJpsi / ptMotherJpsi;

					histPseudoProperDecayLengthJpsi->Fill(pseudoProperDecayLengthJpsi);
					histPseudoProperDecayTimeJpsi->Fill(pseudoProperDecayTimeJpsi);
				}
				*/


				///////////////////////////////////
				/// Pair Statistics Before Cuts ///
				///////////////////////////////////

				// J/psi
				if(isPairElectron && isJpsi) histJpsiPairStatistics->Fill(0.5);
				if(isPairPositron && isJpsi) histJpsiPairStatistics->Fill(1.5);
				if(isPairMuon && isJpsi) histJpsiPairStatistics->Fill(2.5);
				if(isPairAntiMuon && isJpsi) histJpsiPairStatistics->Fill(3.5);

				// Psi(2S)
				if(isPairElectron && isPsi2S) histPsi2sPairStatistics->Fill(0.5);
				if(isPairPositron && isPsi2S) histPsi2sPairStatistics->Fill(1.5);
				if(isPairMuon && isPsi2S) histPsi2sPairStatistics->Fill(2.5);
				if(isPairAntiMuon && isPsi2S) histPsi2sPairStatistics->Fill(3.5);

				// Prompt J/psi and psi(2S)
				if(isPairElectron && isJpsi && isPrompt) histPromptJpsiPairStatistics->Fill(0.5);
				if(isPairPositron && isJpsi && isPrompt) histPromptJpsiPairStatistics->Fill(1.5);
				if(isPairMuon && isJpsi && isPrompt) histPromptJpsiPairStatistics->Fill(2.5);
				if(isPairAntiMuon && isJpsi && isPrompt) histPromptJpsiPairStatistics->Fill(3.5);

				if(isPairElectron && isPsi2S && isPrompt) histPromptPsi2sPairStatistics->Fill(0.5);
				if(isPairPositron && isPsi2S && isPrompt) histPromptPsi2sPairStatistics->Fill(1.5);
				if(isPairMuon && isPsi2S && isPrompt) histPromptPsi2sPairStatistics->Fill(2.5);
				if(isPairAntiMuon && isPsi2S && isPrompt) histPromptPsi2sPairStatistics->Fill(3.5);

				// NonPrompt J/psi and psi(2S)
				if(isPairElectron && isJpsi && !isPrompt) histNonPromptJpsiPairStatistics->Fill(0.5);
				if(isPairPositron && isJpsi && !isPrompt) histNonPromptJpsiPairStatistics->Fill(1.5);
				if(isPairMuon && isJpsi && !isPrompt) histNonPromptJpsiPairStatistics->Fill(2.5);
				if(isPairAntiMuon && isJpsi && !isPrompt) histNonPromptJpsiPairStatistics->Fill(3.5);

				if(isPairElectron && isPsi2S && !isPrompt) histNonPromptPsi2sPairStatistics->Fill(0.5);
				if(isPairPositron && isPsi2S && !isPrompt) histNonPromptPsi2sPairStatistics->Fill(1.5);
				if(isPairMuon && isPsi2S && !isPrompt) histNonPromptPsi2sPairStatistics->Fill(2.5);
				if(isPairAntiMuon && isPsi2S && !isPrompt) histNonPromptPsi2sPairStatistics->Fill(3.5);


				// evaluate inv mass, pt, y of pairs before cuts
				// TODO: ANTİPAİR VE PAİRİN PDG KODLARI MC TRACKLERDEN ALIYOM. AMA BAZEN 1 İ FALSE DÖNÜP BREAK SEGMENT ... VEREBİLİR. HATA HANDLER YAP.
				kMassPair1 = TDatabasePDG::Instance()->GetParticle(kAntiPair->GetPdgCode())->Mass();
				kMassPair2 = TDatabasePDG::Instance()->GetParticle(kPair->GetPdgCode())->Mass();
				kInvMassPair =  kMassPair1*kMassPair1+kMassPair2*kMassPair2 + 2.0*(TMath::Sqrt(kMassPair1*kMassPair1+kPair->GetP()*kPair->GetP())*TMath::Sqrt(kMassPair2*kMassPair2+kAntiPair->GetP()*kAntiPair->GetP()) - kPair->Px()*kAntiPair->Px() - kPair->Py()*kAntiPair->Py() - kPair->Pz()*kAntiPair->Pz());
				kInvMassPair = TMath::Sqrt(kInvMassPair);
				////
				kPx = kPair->Px()+kAntiPair->Px();
				kPy = kPair->Py()+kAntiPair->Py();
				kPtPair = TMath::Sqrt(kPx*kPx + kPy*kPy);
				//kPtPair = kPtPair.GetPt(); TODO: BU GETTERİ KULLANMANIN YOLUNU BUL.
				////
				// TODO: Vertexin birimini kontro let.
				kVx = (kPair->Vx()+kAntiPair->Vx());
				kVy = (kPair->Vy()+kAntiPair->Vy());
				kVz = (kPair->Vz()+kAntiPair->Vz());
				kVertexPair = TMath::Sqrt(kVx*kVx+kVy*kPy);


				/// JPSİ DENEME PSEUDO VERTEX
				// TODO: J/psi için bu doğru mu? Check et.
				if(isJpsi){
				// Get Secondary Vertex
				//if(kCandidateJpsi->isSecondary()) kCandidateJpsiVx = kCandidateJpsi->Vx(); kCandidateJpsiVy = kCandidateJpsi->Vy();
				Double_t pdgMassMotherJpsi = TDatabasePDG::Instance()->GetParticle(443)->Mass();
				//Double_t momentumMotherJpsi = kCandidateJpsi->GetP();
				//Double_t ptMotherJpsi = kCandidateJpsi->GetPt();
				Double_t pseudoProperDecayLengthJpsi = ((kVx - kPrimaryVx) * kPx / kPtPair) + ((kVy - kPrimaryVy) * kPy / kPtPair);
				Double_t pseudoProperDecayTimeJpsi = pseudoProperDecayLengthJpsi * kInvMassPair / kPtPair;

				histPseudoProperDecayLengthJpsi->Fill(pseudoProperDecayLengthJpsi);
				histPseudoProperDecayTimeJpsi->Fill(pseudoProperDecayTimeJpsi);
				}

				/////////////////////////////////
				// Fill Histograms Before Cuts //
				/////////////////////////////////

				/// fiducial cut on generator level for the definition of the acceptance before the application of the single track cuts.
				if(kMidYSimulation     == kTRUE) kGeneratorLevelYCut = TMath::Abs(kY) < 0.9;
				if(kForwardYSimulation == kTRUE) kGeneratorLevelYCut = kY < -2.5 && kY > -4.0;

				//if(kMidYSimulation     == kTRUE) kGeneratorLevelEtaCut = TMath::Abs(kEta) < 0.9;
				//if(kForwardYSimulation == kTRUE) kGeneratorLevelEtaCut = kEta < -2.5 && kEta > -4.0;

				// Applying to cuts at generator level
				if(kGeneratorLevelYCut == kFALSE) continue;
				//if(kGeneratorLevelEtaCut == kFALSE) continue;

				//// fill prompt
				if(isPrompt && isJpsi) {
					histPtPairJpsiPrompt->Fill(kPtPair);
					histMassPairJpsiPrompt->Fill(kInvMassPair);
					histPtJpsiPrompt->Fill(kPt);
					histYJpsiPrompt->Fill(kY);
					histEtaJpsiPrompt->Fill(kEta);
					histThetaJpsiPrompt->Fill(kTheta);
					histPhiJpsiPrompt->Fill(kPhi);
					histOrigin->Fill(0.5);
					PromptJpsiTree->Fill();
					// Analysis Fill
					histPromptJpsiPhiVersusEta->Fill(kPhi,kEta);
					histPromptJpsiEtaVersusPt->Fill(kEta,kPt);
					histPromptJpsiEtaVersusPtPair->Fill(kEta,kPtPair);
					histPromptJpsiMassVersusPt->Fill(kInvMassPair,kPt);
					histPromptJpsiMassVersusPtPair->Fill(kInvMassPair,kPtPair);
					histPromptJpsiMassVersusY->Fill(kInvMassPair,kY);
					histPromptJpsiPtVersusY->Fill(kPt,kY);
					histPromptJpsiPtPairVersusY->Fill(kPtPair,kY);
					histPromptJpsiPtPairVersusPt->Fill(kPtPair,kPt);
				}else if(isPrompt && isPsi2S) {
					histPtPairPsi2sPrompt->Fill(kPtPair);
					histMassPairPsi2sPrompt->Fill(kInvMassPair);
					histPtPsi2sPrompt->Fill(kPt);
					histYPsi2sPrompt->Fill(kY);
					histEtaPsi2sPrompt->Fill(kEta);
					histThetaPsi2sPrompt->Fill(kTheta);
					histPhiPsi2sPrompt->Fill(kPhi);
					histOrigin->Fill(2.5);
					PromptPsi2STree->Fill();
					// Analysis Fill
					histPromptPsi2sPhiVersusEta->Fill(kPhi,kEta);
					histPromptPsi2sEtaVersusPt->Fill(kEta,kPt);
					histPromptPsi2sEtaVersusPtPair->Fill(kEta,kPtPair);
					histPromptPsi2sMassVersusPt->Fill(kInvMassPair,kPt);
					histPromptPsi2sMassVersusPtPair->Fill(kInvMassPair,kPtPair);
					histPromptPsi2sMassVersusY->Fill(kInvMassPair,kY);
					histPromptPsi2sPtVersusY->Fill(kPt,kY);
					histPromptPsi2sPtPairVersusY->Fill(kPtPair,kY);
					histPromptPsi2sPtPairVersusPt->Fill(kPtPair,kPt);
				}
				//// fill non prompt
				else if(!isPrompt && isJpsi) {
					histPtPairJpsiNonPrompt->Fill(kPtPair);
					histMassPairJpsiNonPrompt->Fill(kInvMassPair);
					histPtJpsiNonPrompt->Fill(kPt);
					histYJpsiNonPrompt->Fill(kY);
					histEtaJpsiNonPrompt->Fill(kEta);
					histThetaJpsiNonPrompt->Fill(kTheta);
					histPhiJpsiNonPrompt->Fill(kPhi);
					histVxPairJpsiNonPrompt->Fill(kVx);
					histVyPairJpsiNonPrompt->Fill(kVy);
					histVzPairJpsiNonPrompt->Fill(kVz);
					histVertexPairJpsiNonPrompt->Fill(kVertexPair);
					histOrigin->Fill(1.5);
					NonPromptJpsiTree->Fill();
					// Analysis Fill
					histNonPromptJpsiPhiVersusEta->Fill(kPhi,kEta);
					histNonPromptJpsiEtaVersusPt->Fill(kEta,kPt);
					histNonPromptJpsiEtaVersusPtPair->Fill(kEta,kPtPair);
					histNonPromptJpsiMassVersusPt->Fill(kInvMassPair,kPt);
					histNonPromptJpsiMassVersusPtPair->Fill(kInvMassPair,kPtPair);
					histNonPromptJpsiMassVersusY->Fill(kInvMassPair,kY);
					histNonPromptJpsiPtVersusY->Fill(kPt,kY);
					histNonPromptJpsiPtPairVersusY->Fill(kPtPair,kY);
					histNonPromptJpsiPtPairVersusPt->Fill(kPtPair,kPt);
				}else if(!isPrompt && isPsi2S) {
					histPtPairPsi2sNonPrompt->Fill(kPtPair);
					histMassPairPsi2sNonPrompt->Fill(kInvMassPair);
					histPtPsi2sNonPrompt->Fill(kPt);
					histYPsi2sNonPrompt->Fill(kY);
					histEtaPsi2sNonPrompt->Fill(kEta);
					histThetaPsi2sNonPrompt->Fill(kTheta);
					histPhiPsi2sNonPrompt->Fill(kPhi);
					histOrigin->Fill(3.5);
					NonPromptPsi2STree->Fill();
					// Analysis Fill
					histNonPromptPsi2sPhiVersusEta->Fill(kPhi,kEta);
					histNonPromptPsi2sEtaVersusPt->Fill(kEta,kPt);
					histNonPromptPsi2sEtaVersusPtPair->Fill(kEta,kPtPair);
					histNonPromptPsi2sMassVersusPt->Fill(kInvMassPair,kPt);
					histNonPromptPsi2sMassVersusPtPair->Fill(kInvMassPair,kPtPair);
					histNonPromptPsi2sMassVersusY->Fill(kInvMassPair,kY);
					histNonPromptPsi2sPtVersusY->Fill(kPt,kY);
					histNonPromptPsi2sPtPairVersusY->Fill(kPtPair,kY);
					histNonPromptPsi2sPtPairVersusPt->Fill(kPtPair,kPt);
				}

				//////////////////////////////////
				/// Count Beauty Filling Part ///
				///				Before Cuts					///
				/////////////////////////////////

				// TODO : Make a Loop over tdM Length
				for(int i=0; i<kSizeBeautyPdgCode; i++) {
					auto tdM = mcreader.getTrack(eventID, kMotherTrackID);
					//std::cout << "My TdM Values:" << '\n' << tdM->GetPdgCode();
					//std::cout << "-----------------------" << '\n';
					if (TMath::Abs(tdM->GetPdgCode()) == kBeautyPdgCode[i]) {
						FromBeauty = TMath::Abs(tdM->GetPdgCode());
						if(!isPrompt && isJpsi) {
							if(FromBeauty == kCandidateBzeroPdg) {
								histBeautyMesonsFromJpsi->Fill(0.5);
							}
							if(FromBeauty == kCandidateBplusPdg) {
								histBeautyMesonsFromJpsi->Fill(1.5);
							}
							if(FromBeauty == kCandidateBszeroPdg) {
								histBeautyMesonsFromJpsi->Fill(2.5);
							}
							if(FromBeauty == kCandidateBcplusPdg) {
								histBeautyMesonsFromJpsi->Fill(3.5);
							}
							if(FromBeauty == kCandidateSigmaBPdg) {
								histBeautyBaryonsFromJpsi->Fill(0.5);
							}
							if(FromBeauty == kCandidateLambdaBzeroPdg) {
								histBeautyBaryonsFromJpsi->Fill(1.5);
							}
							if(FromBeauty == kCandidateXiBminusPdg) {
								histBeautyBaryonsFromJpsi->Fill(2.5);
							}
							if(FromBeauty == kCandidateXiBzeroPdg) {
								histBeautyBaryonsFromJpsi->Fill(3.5);
							}
							if(FromBeauty == kCandidateOmegaBminusPdg) {
								histBeautyBaryonsFromPsi2s->Fill(4.5);
							}
						}else if(!isPrompt && isPsi2S) {
							if(FromBeauty == kCandidateBzeroPdg) {
								histBeautyMesonsFromPsi2s->Fill(0.5);
							}
							if(FromBeauty == kCandidateBplusPdg) {
								histBeautyMesonsFromPsi2s->Fill(1.5);
							}
							if(FromBeauty == kCandidateBszeroPdg) {
								histBeautyMesonsFromPsi2s->Fill(2.5);
							}
							if(FromBeauty == kCandidateBcplusPdg) {
								histBeautyMesonsFromPsi2s->Fill(3.5);
							}
							if(FromBeauty == kCandidateSigmaBPdg) {
								histBeautyBaryonsFromPsi2s->Fill(0.5);
							}
							if(FromBeauty == kCandidateLambdaBzeroPdg) {
								histBeautyBaryonsFromPsi2s->Fill(1.5);
							}
							if(FromBeauty == kCandidateXiBminusPdg) {
								histBeautyBaryonsFromPsi2s->Fill(2.5);
							}
							if(FromBeauty == kCandidateXiBzeroPdg) {
								histBeautyBaryonsFromPsi2s->Fill(3.5);
							}
							if(FromBeauty == kCandidateOmegaBminusPdg) {
								histBeautyBaryonsFromPsi2s->Fill(4.5);
							}
						} // if(!isPrompt && isJpsi)
					} // if (TMath::Abs(tdM->GetPdgCode()) == kBeautyPdgCode[i])
				} // for(int i=0; i<kSizeBeautyPdgCode; i++)

				//////////////////////////////////////
				// Cut Definition and Applying Part //
				//////////////////////////////////////

				//                       			 DEFAULT CUTS                   			 //
				//////////////////////////////////////////////////////////////////////
				//--------------------------------------------------------------------
				//                 					TRACK BARREL CUTS												//
				//--------------------------------------------------------------------
				//  For Mid Rapidity,     pT > 1   &&      |eta| <  0.9  						//
				//------------------------------------------------------------------//
				//                  				MUON QUALITY CUTYS											//
				//------------------------------------------------------------------//
				//  For Forward Rapidity, pT > 1   &&  -4 < eta  < -2.5  						//
				//------------------------------------------------------------------//
				//													MASS SELECTIONS													//
				//------------------------------------------------------------------//
				//  For Jpsi/Psi(2S), 3.3 > invMass > 2.8  &&  3.9 > invMass > 3.4  //
				//------------------------------------------------------------------//

				// Cut Definition Part. You can add specific kinematic cuts here.
				kPtCut = kPair->GetPt() > 1.0 && kAntiPair->GetPt() > 1.0;
				if(kMidYSimulation == kTRUE) {
					kEtaCut    =   TMath::Abs(kPair->GetEta()) < 0.9 && TMath::Abs(kAntiPair->GetEta()) < 0.9;
					kYCut      =   TMath::Abs(kPair->GetRapidity()) < 0.9 && TMath::Abs(kAntiPair->GetRapidity()) < 0.9;
				}
				if(kForwardYSimulation == kTRUE) {
					kEtaCut    =   kPair->GetEta() < -2.5 && kAntiPair->GetEta() < -2.5 && kPair->GetEta() > -4.0 && kAntiPair->GetEta() > -4.0;
					kYCut      =   kPair->GetRapidity() < -2.5 && kAntiPair->GetRapidity() < -2.5 && kPair->GetRapidity() > -4.0 && kAntiPair->GetRapidity() > -4.0;
				}

				// Applying The Cuts
				if(kPtCut  == kFALSE) continue;
				if(kEtaCut == kFALSE) continue;
				if(kYCut   == kFALSE) continue;

				// J/psi
				if(isPairElectron && isJpsi) histJpsiPairStatisticsAfterCuts->Fill(0.5);
				if(isPairPositron && isJpsi) histJpsiPairStatisticsAfterCuts->Fill(1.5);
				if(isPairMuon && isJpsi) histJpsiPairStatisticsAfterCuts->Fill(2.5);
				if(isPairAntiMuon && isJpsi) histJpsiPairStatisticsAfterCuts->Fill(3.5);

				// Psi(2S)
				if(isPairElectron && isPsi2S) histPsi2sPairStatisticsAfterCuts->Fill(0.5);
				if(isPairPositron && isPsi2S) histPsi2sPairStatisticsAfterCuts->Fill(1.5);
				if(isPairMuon && isPsi2S) histPsi2sPairStatisticsAfterCuts->Fill(2.5);
				if(isPairAntiMuon && isPsi2S) histPsi2sPairStatisticsAfterCuts->Fill(3.5);

				// Prompt J/psi and psi(2S)
				if(isPairElectron && isJpsi && isPrompt) histPromptJpsiPairStatisticsAfterCuts->Fill(0.5);
				if(isPairPositron && isJpsi && isPrompt) histPromptJpsiPairStatisticsAfterCuts->Fill(1.5);
				if(isPairMuon && isJpsi && isPrompt) histPromptJpsiPairStatisticsAfterCuts->Fill(2.5);
				if(isPairAntiMuon && isJpsi && isPrompt) histPromptJpsiPairStatisticsAfterCuts->Fill(3.5);

				if(isPairElectron && isPsi2S && isPrompt) histPromptPsi2sPairStatisticsAfterCuts->Fill(0.5);
				if(isPairPositron && isPsi2S && isPrompt) histPromptPsi2sPairStatisticsAfterCuts->Fill(1.5);
				if(isPairMuon && isPsi2S && isPrompt) histPromptPsi2sPairStatisticsAfterCuts->Fill(2.5);
				if(isPairAntiMuon && isPsi2S && isPrompt) histPromptPsi2sPairStatisticsAfterCuts->Fill(3.5);

				// NonPrompt J/psi and psi(2S)
				if(isPairElectron && isJpsi && !isPrompt) histNonPromptJpsiPairStatisticsAfterCuts->Fill(0.5);
				if(isPairPositron && isJpsi && !isPrompt) histNonPromptJpsiPairStatisticsAfterCuts->Fill(1.5);
				if(isPairMuon && isJpsi && !isPrompt) histNonPromptJpsiPairStatisticsAfterCuts->Fill(2.5);
				if(isPairAntiMuon && isJpsi && !isPrompt) histNonPromptJpsiPairStatisticsAfterCuts->Fill(3.5);

				if(isPairElectron && isPsi2S && !isPrompt) histNonPromptPsi2sPairStatisticsAfterCuts->Fill(0.5);
				if(isPairPositron && isPsi2S && !isPrompt) histNonPromptPsi2sPairStatisticsAfterCuts->Fill(1.5);
				if(isPairMuon && isPsi2S && !isPrompt) histNonPromptPsi2sPairStatisticsAfterCuts->Fill(2.5);
				if(isPairAntiMuon && isPsi2S && !isPrompt) histNonPromptPsi2sPairStatisticsAfterCuts->Fill(3.5);

				//// re-evaluate inv mass, pt, y of pairs after cuts
				kMassPair1 = TDatabasePDG::Instance()->GetParticle(kSelectionPdg)->Mass();
				kMassPair2 = TDatabasePDG::Instance()->GetParticle(kSelectionPdg)->Mass();
				kInvMassPair =  kMassPair1*kMassPair1+kMassPair2*kMassPair2 + 2.0*(TMath::Sqrt(kMassPair1*kMassPair1+kPair->GetP()*kPair->GetP())*TMath::Sqrt(kMassPair2*kMassPair2+kAntiPair->GetP()*kAntiPair->GetP()) - kPair->Px()*kAntiPair->Px() - kPair->Py()*kAntiPair->Py() - kPair->Pz()*kAntiPair->Pz());
				kInvMassPair = TMath::Sqrt(kInvMassPair);
				////
				kPx = kPair->Px()+kAntiPair->Px();
				kPy = kPair->Py()+kAntiPair->Py();
				kPtPair = TMath::Sqrt(kPx*kPx + kPy*kPy);

				// invMass cut
				/*
				if(isJpsi){
			  	kInvMassCutJpsi  = kInvMassPair > 2.8 && kInvMassPair < 3.3;
					histMassPairJpsiPromptAfterCuts->GetXaxis()->SetRangeUser(2.8,3.3);
					histMassPairJpsiNonPromptAfterCuts->GetXaxis()->SetRangeUser(2.8,3.3);
				}

				if(isPsi2S){
					kInvMassCutPsi2s = kInvMassPair > 3.4 && kInvMassPair < 3.9;
					histMassPairPsi2sPromptAfterCuts->GetXaxis()->SetRangeUser(3.4,3.9);
					histMassPairPsi2sNonPromptAfterCuts->GetXaxis()->SetRangeUser(3.4,3.9);
				}

				if(kInvMassCutJpsi   == kFALSE) continue;
				if(kInvMassCutPsi2s  == kFALSE) continue;
				*/

				// some adressing for re-fillable branches after Cuts
				kInvMassPairBranchCut = kInvMassPair;
				kPtPairBranchCut = kPtPair;
				kPtBranchCut = kPt;
				kYBranchCut = kY;
				kEtaBranchCut = kEta;
				kThetaBranchCut = kTheta;
				kPhiBranchCut = kPhi;

				////////////////////////////////
				// Fill Histograms After Cuts //
				////////////////////////////////

				//// fill prompt
				if(isPrompt && isJpsi) {
					histPtPairJpsiPromptAfterCuts->Fill(kPtPair);
					histMassPairJpsiPromptAfterCuts->Fill(kInvMassPair);
					histPtJpsiPromptAfterCuts->Fill(kPt);
					histYJpsiPromptAfterCuts->Fill(kY);
					histEtaJpsiPromptAfterCuts->Fill(kEta);
					histThetaJpsiPromptAfterCuts->Fill(kTheta);
					histPhiJpsiPromptAfterCuts->Fill(kPhi);
					histPtEfficiencyJpsiPrompt->Fill(kPt);
					histPtPairEfficiencyJpsiPrompt->Fill(kPtPair);
					histYEfficiencyJpsiPrompt->Fill(kY);
					histOriginAfterCuts->Fill(0.5);
					PromptJpsiTreeAfterCuts->Fill();
					// Analysis Fill
					histPromptJpsiPhiVersusEtaAfterCuts->Fill(kPhi,kEta);
					histPromptJpsiEtaVersusPtAfterCuts->Fill(kEta,kPt);
					histPromptJpsiEtaVersusPtPairAfterCuts->Fill(kEta,kPtPair);
					histPromptJpsiMassVersusPtAfterCuts->Fill(kInvMassPair,kPt);
					histPromptJpsiMassVersusPtPairAfterCuts->Fill(kInvMassPair,kPtPair);
					histPromptJpsiMassVersusYAfterCuts->Fill(kInvMassPair,kY);
					histPromptJpsiPtVersusYAfterCuts->Fill(kPt,kY);
					histPromptJpsiPtPairVersusYAfterCuts->Fill(kPtPair,kY);
					histPromptJpsiPtPairVersusPtAfterCuts->Fill(kPtPair,kPt);
				}else if(isPrompt && isPsi2S) {
					histPtPairPsi2sPromptAfterCuts->Fill(kPtPair);
					histMassPairPsi2sPromptAfterCuts->Fill(kInvMassPair);
					histPtPsi2sPromptAfterCuts->Fill(kPt);
					histYPsi2sPromptAfterCuts->Fill(kY);
					histEtaPsi2sPromptAfterCuts->Fill(kEta);
					histThetaPsi2sPromptAfterCuts->Fill(kTheta);
					histPhiPsi2sPromptAfterCuts->Fill(kPhi);
					histPtEfficiencyPsi2sPrompt->Fill(kPt);
					histPtPairEfficiencyPsi2sPrompt->Fill(kPtPair);
					histYEfficiencyPsi2sPrompt->Fill(kY);
					histOriginAfterCuts->Fill(2.5);
					PromptPsi2STreeAfterCuts->Fill();
					// Analysis Fill
					histPromptPsi2sPhiVersusEtaAfterCuts->Fill(kPhi,kEta);
					histPromptPsi2sEtaVersusPtAfterCuts->Fill(kEta,kPt);
					histPromptPsi2sEtaVersusPtPairAfterCuts->Fill(kEta,kPtPair);
					histPromptPsi2sMassVersusPtAfterCuts->Fill(kInvMassPair,kPt);
					histPromptPsi2sMassVersusPtPairAfterCuts->Fill(kInvMassPair,kPtPair);
					histPromptPsi2sMassVersusYAfterCuts->Fill(kInvMassPair,kY);
					histPromptPsi2sPtVersusYAfterCuts->Fill(kPt,kY);
					histPromptPsi2sPtPairVersusYAfterCuts->Fill(kPtPair,kY);
					histPromptPsi2sPtPairVersusPtAfterCuts->Fill(kPtPair,kPt);
				}
				//// fill non prompt
				else if(!isPrompt && isJpsi) {
					histPtPairJpsiNonPromptAfterCuts->Fill(kPtPair);
					histMassPairJpsiNonPromptAfterCuts->Fill(kInvMassPair);
					histPtJpsiNonPromptAfterCuts->Fill(kPt);
					histYJpsiNonPromptAfterCuts->Fill(kY);
					histEtaJpsiNonPromptAfterCuts->Fill(kEta);
					histThetaJpsiNonPromptAfterCuts->Fill(kTheta);
					histPhiJpsiNonPromptAfterCuts->Fill(kPhi);
					histPtEfficiencyJpsiNonPrompt->Fill(kPt);
					histPtPairEfficiencyJpsiNonPrompt->Fill(kPtPair);
					histYEfficiencyJpsiNonPrompt->Fill(kY);
					histOriginAfterCuts->Fill(1.5);
					NonPromptJpsiTreeAfterCuts->Fill();
					// Analysis Fill
					histNonPromptJpsiPhiVersusEtaAfterCuts->Fill(kPhi,kEta);
					histNonPromptJpsiEtaVersusPtAfterCuts->Fill(kEta,kPt);
					histNonPromptJpsiEtaVersusPtPairAfterCuts->Fill(kEta,kPtPair);
					histNonPromptJpsiMassVersusPtAfterCuts->Fill(kInvMassPair,kPt);
					histNonPromptJpsiMassVersusPtPairAfterCuts->Fill(kInvMassPair,kPtPair);
					histNonPromptJpsiMassVersusYAfterCuts->Fill(kInvMassPair,kY);
					histNonPromptJpsiPtVersusYAfterCuts->Fill(kPt,kY);
					histNonPromptJpsiPtPairVersusYAfterCuts->Fill(kPtPair,kY);
					histNonPromptJpsiPtPairVersusPtAfterCuts->Fill(kPtPair,kPt);
				}else if(!isPrompt && isPsi2S) {
					histPtPairPsi2sNonPromptAfterCuts->Fill(kPtPair);
					histMassPairPsi2sNonPromptAfterCuts->Fill(kInvMassPair);
					histPtPsi2sNonPromptAfterCuts->Fill(kPt);
					histYPsi2sNonPromptAfterCuts->Fill(kY);
					histEtaPsi2sNonPromptAfterCuts->Fill(kEta);
					histThetaPsi2sNonPromptAfterCuts->Fill(kTheta);
					histPhiPsi2sNonPromptAfterCuts->Fill(kPhi);
					histPtEfficiencyPsi2sNonPrompt->Fill(kPt);
					histPtPairEfficiencyPsi2sNonPrompt->Fill(kPtPair);
					histYEfficiencyPsi2sPrompt->Fill(kY);
					histOriginAfterCuts->Fill(3.5);
					NonPromptPsi2STreeAfterCuts->Fill();
					// Analysis Fill
					histNonPromptPsi2sPhiVersusEtaAfterCuts->Fill(kPhi,kEta);
					histNonPromptPsi2sEtaVersusPtAfterCuts->Fill(kEta,kPt);
					histNonPromptPsi2sEtaVersusPtPairAfterCuts->Fill(kEta,kPtPair);
					histNonPromptPsi2sMassVersusPtAfterCuts->Fill(kInvMassPair,kPt);
					histNonPromptPsi2sMassVersusPtPairAfterCuts->Fill(kInvMassPair,kPtPair);
					histNonPromptPsi2sMassVersusYAfterCuts->Fill(kInvMassPair,kY);
					histNonPromptPsi2sPtVersusYAfterCuts->Fill(kPt,kY);
					histNonPromptPsi2sPtPairVersusYAfterCuts->Fill(kPtPair,kY);
					histNonPromptPsi2sPtPairVersusPtAfterCuts->Fill(kPtPair,kPt);
				}

				//////////////////////////////////
				/// Count Beauty Filling Part ///
				///				  After Cuts			  ///
				/////////////////////////////////

				// TODO : Make a Loop over tdM Length
				for(int i=0; i<kSizeBeautyPdgCode; i++) {
					auto tdM = mcreader.getTrack(eventID, kMotherTrackID);
					//std::cout << "My TdM Values:" << '\n' << tdM->GetPdgCode();
					//std::cout << "-----------------------" << '\n';
					if (TMath::Abs(tdM->GetPdgCode()) == kBeautyPdgCode[i]) {
						FromBeauty = TMath::Abs(tdM->GetPdgCode());
						if(!isPrompt && isJpsi) {
							if(FromBeauty == kCandidateBzeroPdg) {
								histBeautyMesonsFromJpsiAfterCuts->Fill(0.5);
							}
							if(FromBeauty == kCandidateBplusPdg) {
								histBeautyMesonsFromJpsiAfterCuts->Fill(1.5);
							}
							if(FromBeauty == kCandidateBszeroPdg) {
								histBeautyMesonsFromJpsiAfterCuts->Fill(2.5);
							}
							if(FromBeauty == kCandidateBcplusPdg) {
								histBeautyMesonsFromJpsiAfterCuts->Fill(3.5);
							}
							if(FromBeauty == kCandidateSigmaBPdg) {
								histBeautyBaryonsFromJpsiAfterCuts->Fill(0.5);
							}
							if(FromBeauty == kCandidateLambdaBzeroPdg) {
								histBeautyBaryonsFromJpsiAfterCuts->Fill(1.5);
							}
							if(FromBeauty == kCandidateXiBminusPdg) {
								histBeautyBaryonsFromJpsiAfterCuts->Fill(2.5);
							}
							if(FromBeauty == kCandidateXiBzeroPdg) {
								histBeautyBaryonsFromJpsiAfterCuts->Fill(3.5);
							}
							if(FromBeauty == kCandidateOmegaBminusPdg) {
								histBeautyBaryonsFromPsi2sAfterCuts->Fill(4.5);
							}
						}else if(!isPrompt && isPsi2S) {
							if(FromBeauty == kCandidateBzeroPdg) {
								histBeautyMesonsFromPsi2sAfterCuts->Fill(0.5);
							}
							if(FromBeauty == kCandidateBplusPdg) {
								histBeautyMesonsFromPsi2sAfterCuts->Fill(1.5);
							}
							if(FromBeauty == kCandidateBszeroPdg) {
								histBeautyMesonsFromPsi2sAfterCuts->Fill(2.5);
							}
							if(FromBeauty == kCandidateBcplusPdg) {
								histBeautyMesonsFromPsi2sAfterCuts->Fill(3.5);
							}
							if(FromBeauty == kCandidateSigmaBPdg) {
								histBeautyBaryonsFromPsi2sAfterCuts->Fill(0.5);
							}
							if(FromBeauty == kCandidateLambdaBzeroPdg) {
								histBeautyBaryonsFromPsi2sAfterCuts->Fill(1.5);
							}
							if(FromBeauty == kCandidateXiBminusPdg) {
								histBeautyBaryonsFromPsi2sAfterCuts->Fill(2.5);
							}
							if(FromBeauty == kCandidateXiBzeroPdg) {
								histBeautyBaryonsFromPsi2sAfterCuts->Fill(3.5);
							}
							if(FromBeauty == kCandidateOmegaBminusPdg) {
								histBeautyBaryonsFromPsi2sAfterCuts->Fill(4.5);
							}
						} // else if(!isPrompt && isPsi2S)
					} //if (TMath::Abs(tdM->GetPdgCode()) == kBeautyPdgCode[i])
				}  // for(int i=0; i<kSizeBeautyPdgCode; i++)
			} // if( isJpsi || isPsi2S )
			// TODO : ADD DESTRUCTOR
			//~o2::MCTrack<float>::kPair() = default;
			//~o2::MCTrack() = default;
			//~o2::MCTrackT();
			//~o2::MCTrack();
			//o2::MCTrack *kPair
			//~MCTrackT() = default;
			//~mctracks() = default;
			//~*mctracks() = default;
			//~kPtPair() = default;
			//~kAntiPair() = default;
			//~(o2::MCTrack*)kDaughterTrack() = default;
			//std::cout << "DENEME\n" << kPtPair << '\n';

		} // for (auto& t : *mctracks)

		if (hitbr) {
			assert(trackidsinITS.size() == trackidsinITS_fromhits.size());
			for (auto id : trackidsinITS) {
				assert(trackidsinITS_fromhits[id] == true);
			}
		}
		// Count TPC Tracks and Debugging
		LOG(debug) << "Have " << trackidsinTPC.size() << " tracks with hits in TPC";
		histTPCHits->Fill(trackidsinTPC.size());
		histITSHits->Fill(trackidsinITS.size());
		histTPCVersusITS->Fill(trackidsinTPC.size(),trackidsinITS.size());
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
	} // for (int eventID = 0; eventID < getBranchMCTrack->GetEntries(); ++eventID)

	// CHECK FOR MEMORY lEAK
	/*
	std::cout << trackrefs << '\n';
	*/

	////////////////////////////////////////////
	// Histogram Filling Checking for Drawing //
	////////////////////////////////////////////

	if(histPtJpsiPrompt->GetEntries()    > 0)   {  kPromptCheck    = kTRUE; }
	if(histPtJpsiNonPrompt->GetEntries() > 0)   {  kNonPromptCheck = kTRUE; }

	////////////////////////////////////////////////
	// Declare some bin parameters for BW Fitting //
	////////////////////////////////////////////////

	// For J/Psi Mass
	int divisionMass  = histMassPairJpsiPromptAfterCuts->GetNbinsX();
	float massMIN       = histMassPairJpsiPromptAfterCuts->GetBinLowEdge(1);
	float massMAX       = histMassPairJpsiPromptAfterCuts->GetBinLowEdge(divisionMass+1);
	float BIN_SIZE_MASS = histMassPairJpsiPromptAfterCuts->GetBinWidth(1);

	// For Psi2s Mass
	int divisionMass_psi2s  = histMassPairPsi2sPromptAfterCuts->GetNbinsX();
	float massMIN_psi2s       = histMassPairPsi2sPromptAfterCuts->GetBinLowEdge(1);
	float massMAX_psi2s       = histMassPairPsi2sPromptAfterCuts->GetBinLowEdge(divisionMass_psi2s+1);
	float BIN_SIZE_MASS_psi2s = histMassPairPsi2sPromptAfterCuts->GetBinWidth(1);

	/// Set Parameters to Breit-Wigner Fit Functions for Mass, pT and Y

	TF1 *bwFuncMass_jpsi = new TF1("mybwMass_jpsi",mybwMass,massMIN, massMAX,3);
	bwFuncMass_jpsi->SetParameter(0,1.0);     bwFuncMass_jpsi->SetParName(0,"const");
	bwFuncMass_jpsi->SetParameter(0,5.0);     bwFuncMass_jpsi->SetParName(1,"sigma");
	bwFuncMass_jpsi->SetParameter(1,95.0);     bwFuncMass_jpsi->SetParName(2,"mean");

	TF1 *bwFuncMass_psi2s = new TF1("mybwMass_psi2s",mybwMass,massMIN_psi2s, massMAX_psi2s,3);
	bwFuncMass_psi2s->SetParameter(0,1.0);     bwFuncMass_psi2s->SetParName(0,"const");
	bwFuncMass_psi2s->SetParameter(2,5.0);     bwFuncMass_psi2s->SetParName(1,"sigma");
	bwFuncMass_psi2s->SetParameter(1,95.0);    bwFuncMass_psi2s->SetParName(2,"mean");

	////////////////////////////////////////
	// Efficiencies Part Calculation Part //
	////////////////////////////////////////

	// Efficiency = N_Accepted (After Cuts) / N_Recorded (Before Cuts)

	if(kPromptCheck == kTRUE) {
		// Efficiency for Prompt
		histPtEfficiencyJpsiPrompt->Divide(histPtJpsiPrompt);
		histPtPairEfficiencyJpsiPrompt->Divide(histPtPairJpsiPrompt);
		histPtEfficiencyPsi2sPrompt->Divide(histPtPsi2sPrompt);
		histPtPairEfficiencyPsi2sPrompt->Divide(histPtPairPsi2sPrompt);
		histYEfficiencyJpsiPrompt->Divide(histYJpsiPrompt);
		histYEfficiencyPsi2sPrompt->Divide(histYPsi2sPrompt);

		// Set Range Efficiencies for Prompt
		histPtEfficiencyJpsiPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
		histPtPairEfficiencyJpsiPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
		histPtEfficiencyPsi2sPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
		histPtPairEfficiencyPsi2sPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
		histYEfficiencyJpsiPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
		histYEfficiencyPsi2sPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
	}

	if(kNonPromptCheck == kTRUE) {
		// Efficiency for Non-Prompt
		histPtEfficiencyJpsiNonPrompt->Divide(histPtJpsiNonPrompt);
		histPtPairEfficiencyJpsiNonPrompt->Divide(histPtPairJpsiNonPrompt);
		histPtEfficiencyPsi2sNonPrompt->Divide(histPtPsi2sNonPrompt);
		histPtPairEfficiencyPsi2sNonPrompt->Divide(histPtPairPsi2sNonPrompt);
		histYEfficiencyJpsiNonPrompt->Divide(histYJpsiNonPrompt);
		histYEfficiencyPsi2sNonPrompt->Divide(histYPsi2sNonPrompt);

		// Set Range Efficiencies for Non-Prompt
		histPtEfficiencyJpsiNonPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
		histPtPairEfficiencyJpsiNonPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
		histPtEfficiencyPsi2sNonPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
		histPtPairEfficiencyPsi2sNonPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
		histYEfficiencyJpsiNonPrompt->GetYaxis()->SetRangeUser(0.2,1.0);
		histYEfficiencyPsi2sNonPrompt->GetYaxis()->SetRangeUser(0.2,1.0);;
	}

	/// Writing Histograms to Disk (Output will histKine.root file) --> For Before Cuts
	TFile foutput1("histSgnKine.root","RECREATE");
	PromptJpsiTree->Write();
	PromptPsi2STree->Write();
	NonPromptJpsiTree->Write();
	NonPromptPsi2STree->Write();

	/// Write Histograms for Prompt J/Psi
	if(kPromptCheck == kTRUE) {
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
	if(kNonPromptCheck == kTRUE) {
		histMassPairJpsiNonPrompt->Write();
		histPtPairJpsiNonPrompt->Write();
		histPtJpsiNonPrompt->Write();
		histYJpsiNonPrompt->Write();
		histEtaJpsiNonPrompt->Write();
		histThetaJpsiNonPrompt->Write();
		histPhiJpsiNonPrompt->Write();
		histVxPairJpsiNonPrompt->Write();								// new
		histVyPairJpsiNonPrompt->Write();        				// new
		histVzPairJpsiNonPrompt->Write();								// new
		histVertexPairJpsiNonPrompt->Write(); 					// new
		histMassPairPsi2sNonPrompt->Write();
		histPtPairPsi2sNonPrompt->Write();
		histPtPsi2sNonPrompt->Write();
		histYPsi2sNonPrompt->Write();
		histEtaPsi2sNonPrompt->Write();
		histThetaPsi2sNonPrompt->Write();
		histPhiPsi2sNonPrompt->Write();
	}

	/// Writing Histograms to Disk (Output will histCutKine.root file) --> For After Cuts
	TFile foutput2("histSgnCutKine.root","RECREATE");
	PromptJpsiTreeAfterCuts->Write();
	PromptPsi2STreeAfterCuts->Write();
	NonPromptJpsiTreeAfterCuts->Write();
	NonPromptPsi2STreeAfterCuts->Write();

	/// Write Histograms for Prompt J/Psi and Psi(2S)
	if(kPromptCheck == kTRUE) {
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
	if(kNonPromptCheck == kTRUE) {
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

	/// Writing Histograms to Disk (Output will histCompareAfterCuts.root file) --> For Comparing J/Psi and Psi(2S)
	TFile foutput3("histSgnCompareAfterCuts.root","RECREATE");

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

	if(kPromptCheck == kTRUE) {
		///Fit Inv Mass Histograms for Prompt J/Psi and Psi(2S)
		histMassPairJpsiPromptAfterCuts->Fit("mybwMass_jpsi","QR");       TF1 *fitMassPairJpsiPromptAfterCuts   = histMassPairJpsiPrompt->GetFunction("mybwMass_jpsi");
		histMassPairPsi2sPromptAfterCuts->Fit("mybwMass_psi2s","QR");     TF1 *fitMassPairPsi2sPromptAfterCuts  = histMassPairPsi2sPromptAfterCuts->GetFunction("mybwMass_psi2s");
	}
	if(kNonPromptCheck == kTRUE) {
		/// Fit Inv Mass Histograms for Non - Prompt J/Psi and Non-Prompt Psi(2S)
		histMassPairJpsiNonPromptAfterCuts->Fit("mybwMass_jpsi","QR");    TF1 *fitMassPairJpsiNonPromptAfterCuts   = histMassPairJpsiNonPromptAfterCuts->GetFunction("mybwMass_jpsi");
		histMassPairPsi2sNonPromptAfterCuts->Fit("mybwMass_psi2s","QR");  TF1 *fitMassPairPsi2sNonPromptAfterCuts  = histMassPairPsi2sNonPromptAfterCuts->GetFunction("mybwMass_psi2s");
	}


	/// Writing Histograms to Disk (Output will histFittingAfterCuts.root file) --> For Fitted (Breit - Wigner) Histos After Cuts
	TFile foutput4("histSgnAfterCuts.root","RECREATE");

	/// Write Histograms for Prompt J/Psi and Psi(2S)
	if(kPromptCheck == kTRUE) {
		//histEtaPromptAfterCuts->Write("PE0");
		histMassPairJpsiPromptAfterCuts->Write("PE0");
		histMassPairPsi2sPromptAfterCuts->Write("PE0");

	}
	/// Write Histograms for Non - Prompt J/Psi and Psi(2S)
	if(kNonPromptCheck == kTRUE) {
		//histEtaNonPromptAfterCuts->Write("PE0");
		histMassPairJpsiNonPromptAfterCuts->Write("PE0");
		histMassPairPsi2sNonPromptAfterCuts->Write("PE0");
	}
	/// Writing Histograms to Disk (Output will histAcceptance.root file) --> For Acceptance Factor
	TFile foutput5("histSgnAcceptance.root","RECREATE");

	// Prompt J/Psi and Psi(2S)
	if(kPromptCheck == kTRUE) {
		histPtEfficiencyJpsiPrompt->Write();
		histPtPairEfficiencyJpsiPrompt->Write();
		histYEfficiencyJpsiPrompt->Write();
		histPtEfficiencyPsi2sPrompt->Write();
		histPtPairEfficiencyPsi2sPrompt->Write();
		histYEfficiencyPsi2sPrompt->Write();
	}
	// Non-Prompt J/psi and Psi(2S)
	if(kNonPromptCheck == kTRUE) {
		histPtEfficiencyJpsiNonPrompt->Write();
		histPtPairEfficiencyJpsiNonPrompt->Write();
		histYEfficiencyJpsiNonPrompt->Write();
		histPtEfficiencyPsi2sNonPrompt->Write();
		histPtPairEfficiencyPsi2sNonPrompt->Write();
		histYEfficiencyPsi2sNonPrompt->Write();
	}
	/// Writing Histograms to Disk (Output will histAnalysis.root file) --> For check Corelations
	TFile foutput6("histSgnAnalysis.root","RECREATE");

	// before cuts
	if(kPromptCheck == kTRUE){
	histPromptJpsiPhiVersusEta->Write();
	histPromptJpsiEtaVersusPt->Write();
	histPromptJpsiEtaVersusPtPair->Write();
	histPromptJpsiMassVersusPt->Write();
	histPromptJpsiMassVersusPtPair->Write();
	histPromptJpsiMassVersusY->Write();
	histPromptJpsiPtVersusY->Write();
	histPromptJpsiPtPairVersusY->Write();
	histPromptJpsiPtPairVersusPt->Write();

	histPromptPsi2sPhiVersusEta->Write();
	histPromptPsi2sEtaVersusPt->Write();
	histPromptPsi2sEtaVersusPtPair->Write();
	histPromptPsi2sMassVersusPt->Write();
	histPromptPsi2sMassVersusPtPair->Write();
	histPromptPsi2sMassVersusY->Write();
	histPromptPsi2sPtVersusY->Write();
	histPromptPsi2sPtPairVersusY->Write();
	histPromptPsi2sPtPairVersusPt->Write();
  }
	if(kNonPromptCheck == kTRUE){
	histNonPromptJpsiPhiVersusEta->Write();
	histNonPromptJpsiEtaVersusPt->Write();
	histNonPromptJpsiEtaVersusPtPair->Write();
	histNonPromptJpsiMassVersusPt->Write();
	histNonPromptJpsiMassVersusPtPair->Write();
	histNonPromptJpsiMassVersusY->Write();
	histNonPromptJpsiPtVersusY->Write();
	histNonPromptJpsiPtPairVersusY->Write();
	histNonPromptJpsiPtPairVersusPt->Write();

	histNonPromptPsi2sPhiVersusEta->Write();
	histNonPromptPsi2sEtaVersusPt->Write();
	histNonPromptPsi2sEtaVersusPtPair->Write();
	histNonPromptPsi2sMassVersusPt->Write();
	histNonPromptPsi2sMassVersusPtPair->Write();
	histNonPromptPsi2sMassVersusY->Write();
	histNonPromptPsi2sPtVersusY->Write();
	histNonPromptPsi2sPtPairVersusY->Write();
	histNonPromptPsi2sPtPairVersusPt->Write();
  }
	/// Writing Histograms to Disk (Output will histCutAnalysis.root file) --> For check Corelations After Cuts
	TFile foutput7("histSgnCutAnalysis.root","RECREATE");

	// after cuts
	if(kPromptCheck == kTRUE){
	histPromptJpsiPhiVersusEtaAfterCuts->Write();
	histPromptJpsiEtaVersusPtAfterCuts->Write();
	histPromptJpsiEtaVersusPtPairAfterCuts->Write();
	histPromptJpsiMassVersusPtAfterCuts->Write();
	histPromptJpsiMassVersusPtPairAfterCuts->Write();
	histPromptJpsiMassVersusYAfterCuts->Write();
	histPromptJpsiPtVersusYAfterCuts->Write();
	histPromptJpsiPtPairVersusYAfterCuts->Write();
	histPromptJpsiPtPairVersusPtAfterCuts->Write();

	histPromptPsi2sPhiVersusEtaAfterCuts->Write();
	histPromptPsi2sEtaVersusPtAfterCuts->Write();
	histPromptPsi2sEtaVersusPtPairAfterCuts->Write();
	histPromptPsi2sMassVersusPtAfterCuts->Write();
	histPromptPsi2sMassVersusPtPairAfterCuts->Write();
	histPromptPsi2sMassVersusYAfterCuts->Write();
	histPromptPsi2sPtVersusYAfterCuts->Write();
	histPromptPsi2sPtPairVersusYAfterCuts->Write();
	histPromptPsi2sPtPairVersusPtAfterCuts->Write();
  }

	if(kNonPromptCheck == kTRUE){
	histNonPromptJpsiPhiVersusEtaAfterCuts->Write();
	histNonPromptJpsiEtaVersusPtAfterCuts->Write();
	histNonPromptJpsiEtaVersusPtPairAfterCuts->Write();
	histNonPromptJpsiMassVersusPtAfterCuts->Write();
	histNonPromptJpsiMassVersusPtPairAfterCuts->Write();
	histNonPromptJpsiMassVersusYAfterCuts->Write();
	histNonPromptJpsiPtVersusYAfterCuts->Write();
	histNonPromptJpsiPtPairVersusYAfterCuts->Write();
	histNonPromptJpsiPtPairVersusPtAfterCuts->Write();

	histNonPromptPsi2sPhiVersusEtaAfterCuts->Write();
	histNonPromptPsi2sEtaVersusPtAfterCuts->Write();
	histNonPromptPsi2sEtaVersusPtPairAfterCuts->Write();
	histNonPromptPsi2sMassVersusPtAfterCuts->Write();
	histNonPromptPsi2sMassVersusPtPairAfterCuts->Write();
	histNonPromptPsi2sMassVersusYAfterCuts->Write();
	histNonPromptPsi2sPtVersusYAfterCuts->Write();
	histNonPromptPsi2sPtPairVersusYAfterCuts->Write();
	histNonPromptPsi2sPtPairVersusPtAfterCuts->Write();
  }

	// Writing Histograms to Disk (Output will histCutAnalysis.root file) --> For Check Detector Parameters and Correlations
	TFile foutput8("histSgnDetectorAnalysis.root","RECREATE");

	histTPCHits->Write();
	histITSHits->Write();
	histTPCVersusITS->Write();

	// TODO: AFTER CUTS İÇİN DE YAP.
	TFile foutput9("histSgnBeautyAnalysis.root","RECREATE");

	histPseudoProperDecayLengthB0->Write();
	histPseudoProperDecayLengthBplus->Write();
	histPseudoProperDecayLengthBszero->Write();
	histPseudoProperDecayLengthBcplus->Write();
	histPseudoProperDecayLengthSigmaB->Write();
	histPseudoProperDecayLengthLambdaBzero->Write();
	histPseudoProperDecayLengthXiBminus->Write();
	histPseudoProperDecayLengthXiBzero->Write();
	histPseudoProperDecayLengthOmegaBminus->Write();

	histPseudoProperDecayTimeB0->Write();

	histPseudoProperDecayLengthJpsi->Write();
	histPseudoProperDecayTimeJpsi->Write();

	TFile foutput10("histSgnStatistics.root","RECREATE");

	histOrigin->Write();
	histJpsiPairStatistics->Write();
	histPsi2sPairStatistics->Write();
	histPromptJpsiPairStatistics->Write();
	histPromptPsi2sPairStatistics->Write();
	histNonPromptJpsiPairStatistics->Write();
	histNonPromptPsi2sPairStatistics->Write();
	histBeautyMesonsFromJpsi->Write();
	histBeautyBaryonsFromJpsi->Write();
	histBeautyMesonsFromPsi2s->Write();
	histBeautyBaryonsFromPsi2s->Write();

	TFile foutput11("histSgnCutStatistics.root","RECREATE");
	histOriginAfterCuts->Write();
	histJpsiPairStatisticsAfterCuts->Write();
	histPsi2sPairStatisticsAfterCuts->Write();
	histPromptJpsiPairStatisticsAfterCuts->Write();
	histPromptPsi2sPairStatisticsAfterCuts->Write();
	histNonPromptJpsiPairStatisticsAfterCuts->Write();
	histNonPromptPsi2sPairStatisticsAfterCuts->Write();
	histBeautyMesonsFromJpsiAfterCuts->Write();
	histBeautyBaryonsFromJpsiAfterCuts->Write();
	histBeautyMesonsFromPsi2sAfterCuts->Write();
	histBeautyBaryonsFromPsi2sAfterCuts->Write();

	LOG(info) << "STACK TEST SUCCESSFULL\n";
	return 0;
}
