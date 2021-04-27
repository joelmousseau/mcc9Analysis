////////////////////////////////////////////////////////////////////////
// Class:       SidebandTree
// Module Type: analyzer
// File:        SidebandTree_module.cc
//
// Generated at Thu Jun 18 11:14:25 2020 by Supraja Balasubramanian using artmod
// from cetpkgsupport v1_14_01.
////////////////////////////////////////////////////////////////////////

#include "SidebandTree_module.h"

using namespace std;
/*
#define NEventTypes 9 //0 = all (black), 1 = numuCC 1 pi0 (red), 2 = numuCC 0 pi0 (green), 3 = numuCC N pi0 (blue), 4 = numuNC 0 pi0 (yellow), 5 = numuNC N pi0 (pink), 6 = nue (light blue), 7 = outside FV (dark green), 8 = all bgd (purple)
#define NIntTypes 6 //0 = all, 1 = RES, 2 = QE, 3 = DIS, 4 = COHERENT, 5 = MEC
#define NMultiplicities 5 //All, 1, 2, 3, >3
#define NTrackTypes 7 //1 = muon, 2 = proton, 3 = charged pion, 4 = photon, 5 = other em, 6 = overlay, 7 = other
#define NShowerTypes 9 //1 = gamma1, 2 = gamma2, 3 = descendent, 4 = muon, 5 = proton, 6 = pion, 7 = other em, 8 = overlay, 9 = other
#define NPlanes 3
#define NFilterStages 4 //0 = no filter, 1 >= 1 candidate, 2 = 2 candidates, 3 >= 2 candidates
#define NSelectionStages 6 //0 = cc inclusive, 1 = charged pi veto, 2 = 2 showers, 3 = 1 lead, 4 = at least 1 sublead, 5 = selected ccpi0 event
#define NNewEventTypes 11
#define NCuts 10 //0 = all, 1 = ccinc, 2 = veto, 3 = 2 showers, 4 = lead cl, 5 = lead en, 6 = lead ra, 7 = sublead cl, 8 = sublead energy, 9 = npairs
#define nbins_pi0momentum 8
#define nbins_pi0costheta 12
#define nbins_mumomentum 6
#define nbins_mucostheta 9

#define naltbins_pi0momentum 7
#define naltbins_pi0costheta 9
#define naltbins_mumomentum 6
#define naltbins_mucostheta 9*/

//////////////////// FHICL PARAMETER TAGS //////////////////////////////
SidebandTree::SidebandTree(fhicl::ParameterSet const & p) : EDAnalyzer(p){
  fEventWeightTag = p.get<std::string>("EventWeightProducer");
  fMCTruthTag = p.get<std::string>("MCTruthProducer");
  fGTruthTag = p.get<std::string>("GTruthProducer");
  fMCParticleTag = p.get<std::string>("MCParticleProducer");
  fPFParticleTag = p.get<std::string>("PFParticleProducer");
  fPFParticleMetadataTag = p.get<std::string>("PFParticleMetadataProducer");
  fShowerTag = p.get<std::string>("ShowerProducer");
  fMCShowerTag = p.get<std::string>("MCShowerProducer");
  fClusterTag = p.get<std::string>("ClusterProducer");
  fTrackTag = p.get<std::string>("TrackProducer");
  fTrackFitterTag = p.get<std::string>("TrackFitterProducer");
  fCaloTag = p.get<std::string>("CaloProducer");
  fHitTag = p.get<std::string>("HitProducer");
  fHitPartAssnTag = p.get<std::string>("HitPartAssnProducer");
  fMCSFitResultTag = p.get<std::string>("MCSFitResultProducer");
  fIsMC = p.get<bool>("IsMC",true);
  fIsDirt = p.get<bool>("IsDirt",true);
  theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
  detClocks   = lar::providerFrom<detinfo::DetectorClocksService>();
  // SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  geom = lar::providerFrom<geo::Geometry>();
  // _min_track_len = p.get<double>("MinTrackLength", 0.1);

}//Fhicl parameter tags


void SidebandTree::beginJob()  {
  art::ServiceHandle<art::TFileService> tfs;

  trunklength = 3.;
  trunklength2 = 2.;
  trunkwidth = 1.;

  nchargedpi = 0;
  nproton = 0;
  nneutron = 0;
  nmuon = 0;
  nneutralpi = 0;
  nelectron = 0;
  nother = 0;

  nexiting = 0;
  ncontained = 0;

  NNEUTRINO = 0;
  NALL = 0;

  m_wire_spacing = 0.3;
  m_width_dqdx_box = 1.0;
  m_length_dqdx_box = 4.0;
  m_planes = 3;

  dist3d_min = 4.;
  dist3d_max = 68.;
  ang3d_min = 0.7;
  energy2_min = 20.;

  overlay_gain[0] = 245; //235.5; //245;
  overlay_gain[1] = 252; //249.7; //252;
  overlay_gain[2] = 237.6; //248.2;

  data_gain[0] = 232; //230.3; //232;
  data_gain[1] = 249; //237.6; //249;
  data_gain[2] = 243.7; //238.4;  

  bias_2 = 0.877156; //0.76625;
  altbias_2 = 0.806; //0.71;

  n_0_0 = 0;
  n_0_sh = 0;
  n_sh_0 = 0;
  n_sh_sh = 0;
  n_0_t = 0;
  n_t_0 = 0;
  n_t_t = 0;
  n_t_sh = 0;
  n_sh_t = 0;
  
  Gamma1_E1_id = -999;
  Gamma1_E2_id = -999;
  Gamma2_E1_id = -999;
  Gamma2_E2_id = -999;

  signalMassLow    = 60.0;
  signalMassHigh   = 180.0;
  sidebandMassHigh = 100000.0;
  sidebandMassLow  = 180.0;

  minConversionDist = 6.0; 
  maxConversionDist = 82.0;
  minYPlaneEnergy   = 50.0;
  maxYPlaneEnergy   = 300;
  minRadialAngle    = 0.96;

  _tree = 0;
  if(!_tree){
    _tree = tfs->make<TTree>("tree","tree");
    _tree->Branch("run",&_run,"run/I");
    _tree->Branch("subRun",&_subRun,"subRun/I");
    _tree->Branch("pot",&_pot,"pot/D");
  }
  _run = 0;
  _subRun = 0;
  _pot = 0;

  _eventtree = 0;
  if(!_eventtree){
    _eventtree = tfs->make<TTree>("eventtree","eventtree");
    
    _eventtree->Branch("_f_All_Genie_Weights", "std::vector<float>", &_f_All_Genie_Weights);
    _eventtree->Branch("_f_RPA_CCQE_Genie_Weights", "std::vector<float>", &_f_RPA_CCQE_Genie_Weights);
    _eventtree->Branch("_f_XSecShape_CCMEC_Genie_Weights", "std::vector<float>", &_f_XSecShape_CCMEC_Genie_Weights);
    _eventtree->Branch("_f_AxFFCCQEshape_Genie_Weights", "std::vector<float>", &_f_AxFFCCQEshape_Genie_Weights);
    _eventtree->Branch("_f_VecFFCCQEshape_Genie_Weights", "std::vector<float>", &_f_VecFFCCQEshape_Genie_Weights);
    _eventtree->Branch("_f_DecayAngMEC_Genie_Weights", "std::vector<float>", &_f_DecayAngMEC_Genie_Weights);
    _eventtree->Branch("_f_Theta_Delta2Npi_Genie_Weights", "std::vector<float>", &_f_Theta_Delta2Npi_Genie_Weights);
    _eventtree->Branch("_f_expskin_FluxUnisim_Weights", "std::vector<float>", &_f_expskin_FluxUnisim_Weights);
    _eventtree->Branch("_f_horncurrent_FluxUnisim_Weights", "std::vector<float>", &_f_horncurrent_FluxUnisim_Weights);
    _eventtree->Branch("_f_kminus_PrimaryHadronNormalization_Weights", "std::vector<float>", &_f_kminus_PrimaryHadronNormalization_Weights);
    _eventtree->Branch("_f_kplus_PrimaryHadronFeynmanScaling_Weights", "std::vector<float>", &_f_kplus_PrimaryHadronFeynmanScaling_Weights);
    _eventtree->Branch("_f_kzero_PrimaryHadronSanfordWang_Weights", "std::vector<float>", &_f_kzero_PrimaryHadronSanfordWang_Weights);
    _eventtree->Branch("_f_nucleoninexsec_FluxUnisim_Weights", "std::vector<float>", &_f_nucleoninexsec_FluxUnisim_Weights);
    _eventtree->Branch("_f_nucleonqexsec_FluxUnisim_Weights", "std::vector<float>", &_f_nucleonqexsec_FluxUnisim_Weights);
    _eventtree->Branch("_f_nucleontotxsec_FluxUnisim_Weights", "std::vector<float>", &_f_nucleontotxsec_FluxUnisim_Weights);
    _eventtree->Branch("_f_piminus_PrimaryHadronSWCentralSplineVariation_Weights", "std::vector<float>", &_f_piminus_PrimaryHadronSWCentralSplineVariation_Weights);
    _eventtree->Branch("_f_pioninexsec_FluxUnisim_Weights", "std::vector<float>", &_f_pioninexsec_FluxUnisim_Weights);
    _eventtree->Branch("_f_pionqexsec_FluxUnisim_Weights", "std::vector<float>", &_f_pionqexsec_FluxUnisim_Weights);
    _eventtree->Branch("_f_piontotxsec_FluxUnisim_Weights", "std::vector<float>", &_f_piontotxsec_FluxUnisim_Weights);
    _eventtree->Branch("_f_piplus_PrimaryHadronSWCentralSplineVariation_Weights", "std::vector<float>", &_f_piplus_PrimaryHadronSWCentralSplineVariation_Weights);
    _eventtree->Branch("_fCVWeight",&_fCVWeight,"_fCVWeight/F");
    _eventtree->Branch("_fNuCCNC",&_fNuCCNC,"_fNuCCNC/I");
    _eventtree->Branch("_fNuMode",&_fNuMode,"_fNuMode/I");
    _eventtree->Branch("_fNuEnergy",&_fNuEnergy,"_fNuEnergy/F");
    _eventtree->Branch("_fNuVtxX",&_fNuVtxX,"_fNuVtxX/F");
    _eventtree->Branch("_fNuVtxY",&_fNuVtxY,"_fNuVtxY/F");
    _eventtree->Branch("_fNuVtxZ",&_fNuVtxZ,"_fNuVtxZ/F");
    _eventtree->Branch("_fNuPDG",&_fNuPDG,"_fNuPDG/I");
    _eventtree->Branch("_fNuInFV",&_fNuInFV,"_fNuInFV/I");
    _eventtree->Branch("_fResNum",&_fResNum,"_fResNum/I");

    _eventtree->Branch("_fLepPdg",&_fLepPdg,"_fLepPdg/I");
    _eventtree->Branch("_fLepP",&_fLepP,"_fLepP/F");
    _eventtree->Branch("_fLepPx",&_fLepPx,"_fLepPx/F");
    _eventtree->Branch("_fLepPy",&_fLepPy,"_fLepPy/F");
    _eventtree->Branch("_fLepPz",&_fLepPz,"_fLepPz/F");
    _eventtree->Branch("_fLepCosTheta",&_fLepCosTheta,"_fLepCosTheta/F");
    _eventtree->Branch("_fLepPhi",&_fLepPhi,"_fLepPhi/F");

    _eventtree->Branch("_fPi0E",&_fPi0E,"_fPi0E/F");
    _eventtree->Branch("_fPi0P",&_fPi0P,"_fPi0P/F");
    _eventtree->Branch("_fPi0M",&_fPi0M,"_fPi0M/F");
    _eventtree->Branch("_fPi0Px",&_fPi0Px,"_fPi0Px/F");
    _eventtree->Branch("_fPi0Py",&_fPi0Py,"_fPi0Py/F");
    _eventtree->Branch("_fPi0Pz",&_fPi0Pz,"_fPi0Pz/F");
    _eventtree->Branch("_fPi0CosTheta",&_fPi0CosTheta,"_fPi0CosTheta/F");
    _eventtree->Branch("_fPi0Phi",&_fPi0Phi,"_fPi0Phi/F");

    _eventtree->Branch("_fNpi0",&_fNpi0,"_fNpi0/I");
    _eventtree->Branch("_fNmuon",&_fNmuon,"_fNmuon/I");
    _eventtree->Branch("_fNproton",&_fNproton,"_fNproton/I");
    _eventtree->Branch("_fNpiplus",&_fNpiplus,"_fNpiplus/I");
    _eventtree->Branch("_fNneutron",&_fNneutron,"_fNneutron/I");

    _eventtree->Branch("_fGamma1_id",&_fGamma1_id,"_fGamma1_id/I");
    _eventtree->Branch("_fGamma2_id",&_fGamma2_id,"_fGamma2_id/I");
    _eventtree->Branch("_fGamma1_E",&_fGamma1_E,"_fGamma1_E/F");
    _eventtree->Branch("_fGamma2_E",&_fGamma2_E,"_fGamma2_E/F");
    _eventtree->Branch("_fGamma12_CosTheta",&_fGamma12_CosTheta,"_fGamma12_CosTheta/F");
    _eventtree->Branch("_fGamma12_Angle",&_fGamma12_Angle,"_fGamma12_Angle/F");
    _eventtree->Branch("_fGamma1_DirX",&_fGamma1_DirX,"_fGamma1_DirX/F");
    _eventtree->Branch("_fGamma1_DirY",&_fGamma1_DirY,"_fGamma1_DirY/F");
    _eventtree->Branch("_fGamma1_DirZ",&_fGamma1_DirZ,"_fGamma1_DirZ/F");
    _eventtree->Branch("_fGamma2_DirX",&_fGamma2_DirX,"_fGamma2_DirX/F");
    _eventtree->Branch("_fGamma2_DirY",&_fGamma2_DirY,"_fGamma2_DirY/F");
    _eventtree->Branch("_fGamma2_DirZ",&_fGamma2_DirZ,"_fGamma2_DirZ/F");
    _eventtree->Branch("_fGamma1E1_id",&_fGamma1E1_id,"_fGamma1E1_id/I");
    _eventtree->Branch("_fGamma1E2_id",&_fGamma1E2_id,"_fGamma1E2_id/I");
    _eventtree->Branch("_fGamma2E1_id",&_fGamma2E1_id,"_fGamma2E1_id/I");
    _eventtree->Branch("_fGamma2E2_id",&_fGamma2E2_id,"_fGamma2E2_id/I");


    _eventtree->Branch("_fFlashChi2",&_fFlashChi2,"_fFlashChi2/F");
    _eventtree->Branch("_fNTracks",&_fNTracks,"_fNTracks/I");
    _eventtree->Branch("_fNShowers",&_fNShowers,"_fNShowers/I");
    _eventtree->Branch("_fNHitsU",&_fNHitsU,"_fNHitsU/I");
    _eventtree->Branch("_fNHitsV",&_fNHitsV,"_fNHitsV/I");
    _eventtree->Branch("_fNHitsY",&_fNHitsY,"_fNHitsY/I");
   
    _eventtree->Branch("_fCandidateVertexX",&_fCandidateVertexX,"_fCandidateVertexX/F");
    _eventtree->Branch("_fCandidateVertexY",&_fCandidateVertexY,"_fCandidateVertexY/F");
    _eventtree->Branch("_fCandidateVertexZ",&_fCandidateVertexZ,"_fCandidateVertexZ/F");

    //_fTwoPhotonInvariantMass
    _eventtree->Branch("_fTwoPhotonInvariantMass","std::vector<float>",&_fTwoPhotonInvariantMass);

    _eventtree->Branch("_fCandidateMuonTrueComposition",&_fCandidateMuonTrueComposition,"_fCandidateMuonTrueComposition/I");
    _eventtree->Branch("_fCandidateMuonIsContained",&_fCandidateMuonIsContained,"_fCandidateMuonIsContained/I");
    _eventtree->Branch("_fCandidateMuonDist3d",&_fCandidateMuonDist3d,"_fCandidateMuonDist3d/F");
    _eventtree->Branch("_fCandidateMuonAng3d",&_fCandidateMuonAng3d,"_fCandidateMuonAng3d/F");
    _eventtree->Branch("_fCandidateMuonPID",&_fCandidateMuonPID,"_fCandidateMuonPID/F");
    _eventtree->Branch("_fCandidateMuonStartX",&_fCandidateMuonStartX,"_fCandidateMuonStartX/F");
    _eventtree->Branch("_fCandidateMuonStartY",&_fCandidateMuonStartY,"_fCandidateMuonStartY/F");
    _eventtree->Branch("_fCandidateMuonStartZ",&_fCandidateMuonStartZ,"_fCandidateMuonStartZ/F");
    _eventtree->Branch("_fCandidateMuonEndX",&_fCandidateMuonEndX,"_fCandidateMuonEndX/F");
    _eventtree->Branch("_fCandidateMuonEndY",&_fCandidateMuonEndY,"_fCandidateMuonEndY/F");
    _eventtree->Branch("_fCandidateMuonEndZ",&_fCandidateMuonEndZ,"_fCandidateMuonEndZ/F");
    _eventtree->Branch("_fCandidateMuonMomentum",&_fCandidateMuonMomentum,"_fCandidateMuonMomentum/F");
    _eventtree->Branch("_fCandidateMuonMCSMomentum",&_fCandidateMuonMCSMomentum,"_fCandidateMuonMCSMomentum/F");
    _eventtree->Branch("_fCandidateMuonRangeMomentum",&_fCandidateMuonRangeMomentum,"_fCandidateMuonRangeMomentum/F");
    _eventtree->Branch("_fCandidateMuonLength",&_fCandidateMuonLength,"_fCandidateMuonLength/F");
    _eventtree->Branch("_fCandidateMuonCostheta",&_fCandidateMuonCostheta,"_fCandidateMuonCostheta/F");
    _eventtree->Branch("_fCandidateMuonPhi",&_fCandidateMuonPhi,"_fCandidateMuonPhi/F");
    _eventtree->Branch("_fCandidateMuonThetaXZ",&_fCandidateMuonThetaXZ,"_fCandidateMuonThetaXZ/F");
    _eventtree->Branch("_fCandidateMuonThetaYZ",&_fCandidateMuonThetaYZ,"_fCandidateMuonThetaYZ/F");
    _eventtree->Branch("_fCandidateMuonMomentumResolution",&_fCandidateMuonMomentumResolution,"_fCandidateMuonMomentumResolution/F");
    _eventtree->Branch("_fCandidateMuonRangeMomentumResolution",&_fCandidateMuonRangeMomentumResolution,"_fCandidateMuonRangeMomentumResolution/F");
    _eventtree->Branch("_fCandidateMuonMCSMomentumResolution",&_fCandidateMuonMCSMomentumResolution,"_fCandidateMuonMCSMomentumResolution/F");
    _eventtree->Branch("_fCandidateMuonAngleResolution",&_fCandidateMuonAngleResolution,"_fCandidateMuonAngleResolution/F");

    _eventtree->Branch("_fCandidatePi0Momentum",&_fCandidatePi0Momentum,"_fCandidatePi0Momentum/F");
    _eventtree->Branch("_fCandidatePi0Costheta",&_fCandidatePi0Costheta,"_fCandidatePi0Costheta/F");
    _eventtree->Branch("_fCandidatePi0Phi",&_fCandidatePi0Phi,"_fCandidatePi0Phi/F");
    _eventtree->Branch("_fCandidatePi0Mass",&_fCandidatePi0Mass, "_fCandidatePi0Mass/F");
    _eventtree->Branch("_fCandidatePi0Energy",&_fCandidatePi0Energy,"_fCandidatePi0Energy/F");
    _eventtree->Branch("_fCandidatePi0Angle12",&_fCandidatePi0Angle12,"_fCandidatePi0Angle12/F");
    _eventtree->Branch("_fCandidatePi0MomentumResolution",&_fCandidatePi0MomentumResolution,"_fCandidatePi0MomentumResolution/F");
    _eventtree->Branch("_fCandidatePi0AngleResolution",&_fCandidatePi0AngleResolution,"_fCandidatePi0AngleResolution/F");

    _eventtree->Branch("_fCandidateLeadTrueComposition",&_fCandidateLeadTrueComposition,"_fCandidateLeadTrueComposition/I");
    _eventtree->Branch("_fCandidateLeadScore",&_fCandidateLeadScore,"_fCandidateLeadScore/F");
    _eventtree->Branch("_fCandidateLeadLength",&_fCandidateLeadLength,"_fCandidateLeadLength/F");
    _eventtree->Branch("_fCandidateLeadOpenAngle",&_fCandidateLeadOpenAngle,"_fCandidateLeadOpenAngle/F");
    _eventtree->Branch("_fCandidateLeadCostheta",&_fCandidateLeadCostheta,"_fCandidateLeadCostheta/F");
    _eventtree->Branch("_fCandidateLeadPhi",&_fCandidateLeadPhi,"_fCandidateLeadPhi/F");
    _eventtree->Branch("_fCandidateLeadThetaXZ",&_fCandidateLeadThetaXZ,"_fCandidateLeadThetaXZ/F");
    _eventtree->Branch("_fCandidateLeadThetaYZ",&_fCandidateLeadThetaYZ,"_fCandidateLeadThetaYZ/F");
    _eventtree->Branch("_fCandidateLeadDist3d",&_fCandidateLeadDist3d,"_fCandidateLeadDist3d/F");
    _eventtree->Branch("_fCandidateLeadAng3d",&_fCandidateLeadAng3d,"_fCandidateLeadAng3d/F");
    _eventtree->Branch("_fCandidateLeadStartX",&_fCandidateLeadStartX,"_fCandidateLeadStartX/F");
    _eventtree->Branch("_fCandidateLeadStartY",&_fCandidateLeadStartY,"_fCandidateLeadStartY/F");
    _eventtree->Branch("_fCandidateLeadStartZ",&_fCandidateLeadStartZ,"_fCandidateLeadStartZ/F");
    _eventtree->Branch("_fCandidateLeadEnergy0",&_fCandidateLeadEnergy0,"_fCandidateLeadEnergy0/F");
    _eventtree->Branch("_fCandidateLeadEnergy1",&_fCandidateLeadEnergy1,"_fCandidateLeadEnergy1/F");
    _eventtree->Branch("_fCandidateLeadEnergy2",&_fCandidateLeadEnergy2,"_fCandidateLeadEnergy2/F");
    _eventtree->Branch("_fCandidateLeaddEdx0",&_fCandidateLeaddEdx0,"_fCandidateLeaddEdx0/F");
    _eventtree->Branch("_fCandidateLeaddEdx1",&_fCandidateLeaddEdx1,"_fCandidateLeaddEdx1/F");
    _eventtree->Branch("_fCandidateLeaddEdx2",&_fCandidateLeaddEdx2,"_fCandidateLeaddEdx2/F");
    _eventtree->Branch("_fCandidateLeaddEdx3",&_fCandidateLeaddEdx3,"_fCandidateLeaddEdx3/F");
    _eventtree->Branch("_fCandidateLeadEnergy0Resolution",&_fCandidateLeadEnergy0Resolution,"_fCandidateLeadEnergy0Resolution/F");
    _eventtree->Branch("_fCandidateLeadEnergy1Resolution",&_fCandidateLeadEnergy1Resolution,"_fCandidateLeadEnergy1Resolution/F");
    _eventtree->Branch("_fCandidateLeadEnergy2Resolution",&_fCandidateLeadEnergy2Resolution,"_fCandidateLeadEnergy2Resolution/F");
    _eventtree->Branch("_fCandidateLeadAngleResolution",&_fCandidateLeadAngleResolution,"_fCandidateLeadAngleResolution/F");
    _eventtree->Branch("_fCandidateLeadPurity0",&_fCandidateLeadPurity0,"_fCandidateLeadPurity0/F");
    _eventtree->Branch("_fCandidateLeadPurity1",&_fCandidateLeadPurity1,"_fCandidateLeadPurity1/F");
    _eventtree->Branch("_fCandidateLeadPurity2",&_fCandidateLeadPurity2,"_fCandidateLeadPurity2/F");
    _eventtree->Branch("_fCandidateLeadCompleteness0",&_fCandidateLeadCompleteness0,"_fCandidateLeadCompleteness0/F");
    _eventtree->Branch("_fCandidateLeadCompleteness1",&_fCandidateLeadCompleteness1,"_fCandidateLeadCompleteness1/F");
    _eventtree->Branch("_fCandidateLeadCompleteness2",&_fCandidateLeadCompleteness2,"_fCandidateLeadCompleteness2/F");

    _eventtree->Branch("_fCandidateSubleadTrueComposition",&_fCandidateSubleadTrueComposition,"_fCandidateSubleadTrueComposition/I");
    _eventtree->Branch("_fCandidateSubleadScore",&_fCandidateSubleadScore,"_fCandidateSubleadScore/F");
    _eventtree->Branch("_fCandidateSubleadLength",&_fCandidateSubleadLength,"_fCandidateSubleadLength/F");
    _eventtree->Branch("_fCandidateSubleadOpenAngle",&_fCandidateSubleadOpenAngle,"_fCandidateSubleadOpenAngle/F");
    _eventtree->Branch("_fCandidateSubleadCostheta",&_fCandidateSubleadCostheta,"_fCandidateSubleadCostheta/F");
    _eventtree->Branch("_fCandidateSubleadPhi",&_fCandidateSubleadPhi,"_fCandidateSubleadPhi/F");
    _eventtree->Branch("_fCandidateSubleadThetaXZ",&_fCandidateSubleadThetaXZ,"_fCandidateSubleadThetaXZ/F");
    _eventtree->Branch("_fCandidateSubleadThetaYZ",&_fCandidateSubleadThetaYZ,"_fCandidateSubleadThetaYZ/F");
    _eventtree->Branch("_fCandidateSubleadDist3d",&_fCandidateSubleadDist3d,"_fCandidateSubleadDist3d/F");
    _eventtree->Branch("_fCandidateSubleadAng3d",&_fCandidateSubleadAng3d,"_fCandidateSubleadAng3d/F");
    _eventtree->Branch("_fCandidateSubleadStartX",&_fCandidateSubleadStartX,"_fCandidateSubleadStartX/F");
    _eventtree->Branch("_fCandidateSubleadStartY",&_fCandidateSubleadStartY,"_fCandidateSubleadStartY/F");
    _eventtree->Branch("_fCandidateSubleadStartZ",&_fCandidateSubleadStartZ,"_fCandidateSubleadStartZ/F");
    _eventtree->Branch("_fCandidateSubleadEnergy0",&_fCandidateSubleadEnergy0,"_fCandidateSubleadEnergy0/F");
    _eventtree->Branch("_fCandidateSubleadEnergy1",&_fCandidateSubleadEnergy1,"_fCandidateSubleadEnergy1/F");
    _eventtree->Branch("_fCandidateSubleadEnergy2",&_fCandidateSubleadEnergy2,"_fCandidateSubleadEnergy2/F");
    _eventtree->Branch("_fCandidateSubleaddEdx0",&_fCandidateSubleaddEdx0,"_fCandidateSubleaddEdx0/F");
    _eventtree->Branch("_fCandidateSubleaddEdx1",&_fCandidateSubleaddEdx1,"_fCandidateSubleaddEdx1/F");
    _eventtree->Branch("_fCandidateSubleaddEdx2",&_fCandidateSubleaddEdx2,"_fCandidateSubleaddEdx2/F");
    _eventtree->Branch("_fCandidateSubleaddEdx3",&_fCandidateSubleaddEdx3,"_fCandidateSubleaddEdx3/F");
    _eventtree->Branch("_fCandidateSubleadEnergy0Resolution",&_fCandidateSubleadEnergy0Resolution,"_fCandidateSubleadEnergy0Resolution/F");
    _eventtree->Branch("_fCandidateSubleadEnergy1Resolution",&_fCandidateSubleadEnergy1Resolution,"_fCandidateSubleadEnergy1Resolution/F");
    _eventtree->Branch("_fCandidateSubleadEnergy2Resolution",&_fCandidateSubleadEnergy2Resolution,"_fCandidateSubleadEnergy2Resolution/F");
    _eventtree->Branch("_fCandidateSubleadAngleResolution",&_fCandidateSubleadAngleResolution,"_fCandidateSubleadAngleResolution/F");
    _eventtree->Branch("_fCandidateSubleadPurity0",&_fCandidateSubleadPurity0,"_fCandidateSubleadPurity0/F");
    _eventtree->Branch("_fCandidateSubleadPurity1",&_fCandidateSubleadPurity1,"_fCandidateSubleadPurity1/F");
    _eventtree->Branch("_fCandidateSubleadPurity2",&_fCandidateSubleadPurity2,"_fCandidateSubleadPurity2/F");
    _eventtree->Branch("_fCandidateSubleadCompleteness0",&_fCandidateSubleadCompleteness0,"_fCandidateSubleadCompleteness0/F");
    _eventtree->Branch("_fCandidateSubleadCompleteness1",&_fCandidateSubleadCompleteness1,"_fCandidateSubleadCompleteness1/F");
    _eventtree->Branch("_fCandidateSubleadCompleteness2",&_fCandidateSubleadCompleteness2,"_fCandidateSubleadCompleteness2/F");

    _eventtree->Branch("_fTwoMIPPi0Momentum",&_fTwoMIPPi0Momentum,"_fTwoMIPPi0Momentum/F");
    _eventtree->Branch("_fTwoMIPPi0Costheta",&_fTwoMIPPi0Costheta,"_fTwoMIPPi0Costheta/F");
    _eventtree->Branch("_fTwoMIPPi0Phi",&_fTwoMIPPi0Phi,"_fTwoMIPPi0Phi/F");
    _eventtree->Branch("_fTwoMIPPi0Mass", &_fTwoMIPPi0Mass, "&_fTwoMIPPi0Mass/F");
    _eventtree->Branch("_fTwoMIPPi0Energy",&_fTwoMIPPi0Energy,"_fTwoMIPPi0Energy/F");
    _eventtree->Branch("_fTwoMIPPi0Angle12",&_fTwoMIPPi0Angle12,"_fTwoMIPPi0Angle12/F");
    _eventtree->Branch("_fTwoMIPPi0MomentumResolution",&_fTwoMIPPi0MomentumResolution,"_fTwoMIPPi0MomentumResolution/F");
    _eventtree->Branch("_fTwoMIPPi0AngleResolution",&_fTwoMIPPi0AngleResolution,"_fTwoMIPPi0AngleResolution/F");

    _eventtree->Branch("_fTwoMIPLeadTrueComposition",&_fTwoMIPLeadTrueComposition,"_fTwoMIPLeadTrueComposition/I");
    _eventtree->Branch("_fTwoMIPLeadScore",&_fTwoMIPLeadScore,"_fTwoMIPLeadScore/F");
    _eventtree->Branch("_fTwoMIPLeadLength",&_fTwoMIPLeadLength,"_fTwoMIPLeadLength/F");
    _eventtree->Branch("_fTwoMIPLeadOpenAngle",&_fTwoMIPLeadOpenAngle,"_fTwoMIPLeadOpenAngle/F");
    _eventtree->Branch("_fTwoMIPLeadCostheta",&_fTwoMIPLeadCostheta,"_fTwoMIPLeadCostheta/F");
    _eventtree->Branch("_fTwoMIPLeadPhi",&_fTwoMIPLeadPhi,"_fTwoMIPLeadPhi/F");
    _eventtree->Branch("_fTwoMIPLeadThetaXZ",&_fTwoMIPLeadThetaXZ,"_fTwoMIPLeadThetaXZ/F");
    _eventtree->Branch("_fTwoMIPLeadThetaYZ",&_fTwoMIPLeadThetaYZ,"_fTwoMIPLeadThetaYZ/F");
    _eventtree->Branch("_fTwoMIPLeadDist3d",&_fTwoMIPLeadDist3d,"_fTwoMIPLeadDist3d/F");
    _eventtree->Branch("_fTwoMIPLeadAng3d",&_fTwoMIPLeadAng3d,"_fTwoMIPLeadAng3d/F");
    _eventtree->Branch("_fTwoMIPLeadStartX",&_fTwoMIPLeadStartX,"_fTwoMIPLeadStartX/F");
    _eventtree->Branch("_fTwoMIPLeadStartY",&_fTwoMIPLeadStartY,"_fTwoMIPLeadStartY/F");
    _eventtree->Branch("_fTwoMIPLeadStartZ",&_fTwoMIPLeadStartZ,"_fTwoMIPLeadStartZ/F");
    _eventtree->Branch("_fTwoMIPLeadEnergy0",&_fTwoMIPLeadEnergy0,"_fTwoMIPLeadEnergy0/F");
    _eventtree->Branch("_fTwoMIPLeadEnergy1",&_fTwoMIPLeadEnergy1,"_fTwoMIPLeadEnergy1/F");
    _eventtree->Branch("_fTwoMIPLeadEnergy2",&_fTwoMIPLeadEnergy2,"_fTwoMIPLeadEnergy2/F");
    _eventtree->Branch("_fTwoMIPLeaddEdx0",&_fTwoMIPLeaddEdx0,"_fTwoMIPLeaddEdx0/F");
    _eventtree->Branch("_fTwoMIPLeaddEdx1",&_fTwoMIPLeaddEdx1,"_fTwoMIPLeaddEdx1/F");
    _eventtree->Branch("_fTwoMIPLeaddEdx2",&_fTwoMIPLeaddEdx2,"_fTwoMIPLeaddEdx2/F");
    _eventtree->Branch("_fTwoMIPLeaddEdx3",&_fTwoMIPLeaddEdx3,"_fTwoMIPLeaddEdx3/F");
    _eventtree->Branch("_fTwoMIPLeadEnergy0Resolution",&_fTwoMIPLeadEnergy0Resolution,"_fTwoMIPLeadEnergy0Resolution/F");
    _eventtree->Branch("_fTwoMIPLeadEnergy1Resolution",&_fTwoMIPLeadEnergy1Resolution,"_fTwoMIPLeadEnergy1Resolution/F");
    _eventtree->Branch("_fTwoMIPLeadEnergy2Resolution",&_fTwoMIPLeadEnergy2Resolution,"_fTwoMIPLeadEnergy2Resolution/F");
    _eventtree->Branch("_fTwoMIPLeadAngleResolution",&_fTwoMIPLeadAngleResolution,"_fTwoMIPLeadAngleResolution/F");
    _eventtree->Branch("_fTwoMIPLeadPurity0",&_fTwoMIPLeadPurity0,"_fTwoMIPLeadPurity0/F");
    _eventtree->Branch("_fTwoMIPLeadPurity1",&_fTwoMIPLeadPurity1,"_fTwoMIPLeadPurity1/F");
    _eventtree->Branch("_fTwoMIPLeadPurity2",&_fTwoMIPLeadPurity2,"_fTwoMIPLeadPurity2/F");
    _eventtree->Branch("_fTwoMIPLeadCompleteness0",&_fTwoMIPLeadCompleteness0,"_fTwoMIPLeadCompleteness0/F");
    _eventtree->Branch("_fTwoMIPLeadCompleteness1",&_fTwoMIPLeadCompleteness1,"_fTwoMIPLeadCompleteness1/F");
    _eventtree->Branch("_fTwoMIPLeadCompleteness2",&_fTwoMIPLeadCompleteness2,"_fTwoMIPLeadCompleteness2/F");

    _eventtree->Branch("_fTwoMIPSubleadTrueComposition",&_fTwoMIPSubleadTrueComposition,"_fTwoMIPSubleadTrueComposition/I");
    _eventtree->Branch("_fTwoMIPSubleadScore",&_fTwoMIPSubleadScore,"_fTwoMIPSubleadScore/F");
    _eventtree->Branch("_fTwoMIPSubleadLength",&_fTwoMIPSubleadLength,"_fTwoMIPSubleadLength/F");
    _eventtree->Branch("_fTwoMIPSubleadOpenAngle",&_fTwoMIPSubleadOpenAngle,"_fTwoMIPSubleadOpenAngle/F");
    _eventtree->Branch("_fTwoMIPSubleadCostheta",&_fTwoMIPSubleadCostheta,"_fTwoMIPSubleadCostheta/F");
    _eventtree->Branch("_fTwoMIPSubleadPhi",&_fTwoMIPSubleadPhi,"_fTwoMIPSubleadPhi/F");
    _eventtree->Branch("_fTwoMIPSubleadThetaXZ",&_fTwoMIPSubleadThetaXZ,"_fTwoMIPSubleadThetaXZ/F");
    _eventtree->Branch("_fTwoMIPSubleadThetaYZ",&_fTwoMIPSubleadThetaYZ,"_fTwoMIPSubleadThetaYZ/F");
    _eventtree->Branch("_fTwoMIPSubleadDist3d",&_fTwoMIPSubleadDist3d,"_fTwoMIPSubleadDist3d/F");
    _eventtree->Branch("_fTwoMIPSubleadAng3d",&_fTwoMIPSubleadAng3d,"_fTwoMIPSubleadAng3d/F");
    _eventtree->Branch("_fTwoMIPSubleadStartX",&_fTwoMIPSubleadStartX,"_fTwoMIPSubleadStartX/F");
    _eventtree->Branch("_fTwoMIPSubleadStartY",&_fTwoMIPSubleadStartY,"_fTwoMIPSubleadStartY/F");
    _eventtree->Branch("_fTwoMIPSubleadStartZ",&_fTwoMIPSubleadStartZ,"_fTwoMIPSubleadStartZ/F");
    _eventtree->Branch("_fTwoMIPSubleadEnergy0",&_fTwoMIPSubleadEnergy0,"_fTwoMIPSubleadEnergy0/F");
    _eventtree->Branch("_fTwoMIPSubleadEnergy1",&_fTwoMIPSubleadEnergy1,"_fTwoMIPSubleadEnergy1/F");
    _eventtree->Branch("_fTwoMIPSubleadEnergy2",&_fTwoMIPSubleadEnergy2,"_fTwoMIPSubleadEnergy2/F");
    _eventtree->Branch("_fTwoMIPSubleaddEdx0",&_fTwoMIPSubleaddEdx0,"_fTwoMIPSubleaddEdx0/F");
    _eventtree->Branch("_fTwoMIPSubleaddEdx1",&_fTwoMIPSubleaddEdx1,"_fTwoMIPSubleaddEdx1/F");
    _eventtree->Branch("_fTwoMIPSubleaddEdx2",&_fTwoMIPSubleaddEdx2,"_fTwoMIPSubleaddEdx2/F");
    _eventtree->Branch("_fTwoMIPSubleaddEdx3",&_fTwoMIPSubleaddEdx3,"_fTwoMIPSubleaddEdx3/F");
    _eventtree->Branch("_fTwoMIPSubleadEnergy0Resolution",&_fTwoMIPSubleadEnergy0Resolution,"_fTwoMIPSubleadEnergy0Resolution/F");
    _eventtree->Branch("_fTwoMIPSubleadEnergy1Resolution",&_fTwoMIPSubleadEnergy1Resolution,"_fTwoMIPSubleadEnergy1Resolution/F");
    _eventtree->Branch("_fTwoMIPSubleadEnergy2Resolution",&_fTwoMIPSubleadEnergy2Resolution,"_fTwoMIPSubleadEnergy2Resolution/F");
    _eventtree->Branch("_fTwoMIPSubleadAngleResolution",&_fTwoMIPSubleadAngleResolution,"_fTwoMIPSubleadAngleResolution/F");
    _eventtree->Branch("_fTwoMIPSubleadPurity0",&_fTwoMIPSubleadPurity0,"_fTwoMIPSubleadPurity0/F");
    _eventtree->Branch("_fTwoMIPSubleadPurity1",&_fTwoMIPSubleadPurity1,"_fTwoMIPSubleadPurity1/F");
    _eventtree->Branch("_fTwoMIPSubleadPurity2",&_fTwoMIPSubleadPurity2,"_fTwoMIPSubleadPurity2/F");
    _eventtree->Branch("_fTwoMIPSubleadCompleteness0",&_fTwoMIPSubleadCompleteness0,"_fTwoMIPSubleadCompleteness0/F");
    _eventtree->Branch("_fTwoMIPSubleadCompleteness1",&_fTwoMIPSubleadCompleteness1,"_fTwoMIPSubleadCompleteness1/F");
    _eventtree->Branch("_fTwoMIPSubleadCompleteness2",&_fTwoMIPSubleadCompleteness2,"_fTwoMIPSubleadCompleteness2/F");

    _eventtree->Branch("_fMultiPairPi0Momentum","std::vector<float>",&_fMultiPairPi0Momentum);
    _eventtree->Branch("_fMultiPairPi0Costheta","std::vector<float>",&_fMultiPairPi0Costheta);
    _eventtree->Branch("_fMultiPairPi0Phi","std::vector<float>",&_fMultiPairPi0Phi);
    _eventtree->Branch("_fMultiPairPi0Mass","std::vector<float>",&_fMultiPairPi0Mass);
    _eventtree->Branch("_fMultiPairPi0Energy","std::vector<float>",&_fMultiPairPi0Energy);
    _eventtree->Branch("_fMultiPairPi0Angle12","std::vector<float>",&_fMultiPairPi0Angle12);
    _eventtree->Branch("_fMultiPairPi0MomentumResolution","std::vector<float>",&_fMultiPairPi0MomentumResolution);

    _eventtree->Branch("_fMultiPairLeadTrueComposition",&_fMultiPairLeadTrueComposition,"_fMultiPairLeadTrueComposition/I");
    _eventtree->Branch("_fMultiPairLeadScore",&_fMultiPairLeadScore,"_fMultiPairLeadScore/F");
    _eventtree->Branch("_fMultiPairLeadLength",&_fMultiPairLeadLength,"_fMultiPairLeadLength/F");
    _eventtree->Branch("_fMultiPairLeadOpenAngle",&_fMultiPairLeadOpenAngle,"_fMultiPairLeadOpenAngle/F");
    _eventtree->Branch("_fMultiPairLeadCostheta",&_fMultiPairLeadCostheta,"_fMultiPairLeadCostheta/F");
    _eventtree->Branch("_fMultiPairLeadPhi",&_fMultiPairLeadPhi,"_fMultiPairLeadPhi/F");
    _eventtree->Branch("_fMultiPairLeadThetaXZ",&_fMultiPairLeadThetaXZ,"_fMultiPairLeadThetaXZ/F");
    _eventtree->Branch("_fMultiPairLeadThetaYZ",&_fMultiPairLeadThetaYZ,"_fMultiPairLeadThetaYZ/F");
    _eventtree->Branch("_fMultiPairLeadDist3d",&_fMultiPairLeadDist3d,"_fMultiPairLeadDist3d/F");
    _eventtree->Branch("_fMultiPairLeadAng3d",&_fMultiPairLeadAng3d,"_fMultiPairLeadAng3d/F");
    _eventtree->Branch("_fMultiPairLeadStartX",&_fMultiPairLeadStartX,"_fMultiPairLeadStartX/F");
    _eventtree->Branch("_fMultiPairLeadStartY",&_fMultiPairLeadStartY,"_fMultiPairLeadStartY/F");
    _eventtree->Branch("_fMultiPairLeadStartZ",&_fMultiPairLeadStartZ,"_fMultiPairLeadStartZ/F");
    _eventtree->Branch("_fMultiPairLeadEnergy0",&_fMultiPairLeadEnergy0,"_fMultiPairLeadEnergy0/F");
    _eventtree->Branch("_fMultiPairLeadEnergy1",&_fMultiPairLeadEnergy1,"_fMultiPairLeadEnergy1/F");
    _eventtree->Branch("_fMultiPairLeadEnergy2",&_fMultiPairLeadEnergy2,"_fMultiPairLeadEnergy2/F");
    _eventtree->Branch("_fMultiPairLeaddEdx0",&_fMultiPairLeaddEdx0,"_fMultiPairLeaddEdx0/F");
    _eventtree->Branch("_fMultiPairLeaddEdx1",&_fMultiPairLeaddEdx1,"_fMultiPairLeaddEdx1/F");
    _eventtree->Branch("_fMultiPairLeaddEdx2",&_fMultiPairLeaddEdx2,"_fMultiPairLeaddEdx2/F");
    _eventtree->Branch("_fMultiPairLeaddEdx3",&_fMultiPairLeaddEdx3,"_fMultiPairLeaddEdx3/F");
    _eventtree->Branch("_fMultiPairLeadEnergy0Resolution",&_fMultiPairLeadEnergy0Resolution,"_fMultiPairLeadEnergy0Resolution/F");
    _eventtree->Branch("_fMultiPairLeadEnergy1Resolution",&_fMultiPairLeadEnergy1Resolution,"_fMultiPairLeadEnergy1Resolution/F");
    _eventtree->Branch("_fMultiPairLeadEnergy2Resolution",&_fMultiPairLeadEnergy2Resolution,"_fMultiPairLeadEnergy2Resolution/F");
    _eventtree->Branch("_fMultiPairLeadAngleResolution",&_fMultiPairLeadAngleResolution,"_fMultiPairLeadAngleResolution/F");
    _eventtree->Branch("_fMultiPairLeadPurity0",&_fMultiPairLeadPurity0,"_fMultiPairLeadPurity0/F");
    _eventtree->Branch("_fMultiPairLeadPurity1",&_fMultiPairLeadPurity1,"_fMultiPairLeadPurity1/F");
    _eventtree->Branch("_fMultiPairLeadPurity2",&_fMultiPairLeadPurity2,"_fMultiPairLeadPurity2/F");
    _eventtree->Branch("_fMultiPairLeadCompleteness0",&_fMultiPairLeadCompleteness0,"_fMultiPairLeadCompleteness0/F");
    _eventtree->Branch("_fMultiPairLeadCompleteness1",&_fMultiPairLeadCompleteness1,"_fMultiPairLeadCompleteness1/F");
    _eventtree->Branch("_fMultiPairLeadCompleteness2",&_fMultiPairLeadCompleteness2,"_fMultiPairLeadCompleteness2/F");
    
    //These should be vectors
    _eventtree->Branch("_fMultiPairSubleadTrueComposition", "std::vector<int>", &_fMultiPairSubleadTrueComposition);
    _eventtree->Branch("_fMultiPairSubleadScore", "std::vector<float>", &_fMultiPairSubleadScore);
    _eventtree->Branch("_fMultiPairSubleadLength", "std::vector<float>",&_fMultiPairSubleadLength);
    _eventtree->Branch("_fMultiPairSubleadOpenAngle", "std::vector<float>",&_fMultiPairSubleadOpenAngle);
    _eventtree->Branch("_fMultiPairSubleadCostheta", "std::vector<float>",&_fMultiPairSubleadCostheta);
    _eventtree->Branch("_fMultiPairSubleadPhi", "std::vector<float>",&_fMultiPairSubleadPhi);
    _eventtree->Branch("_fMultiPairSubleadThetaXZ", "std::vector<float>",&_fMultiPairSubleadThetaXZ);
    _eventtree->Branch("_fMultiPairSubleadThetaYZ", "std::vector<float>",&_fMultiPairSubleadThetaYZ);
    _eventtree->Branch("_fMultiPairSubleadDist3d", "std::vector<float>",&_fMultiPairSubleadDist3d);
    _eventtree->Branch("_fMultiPairSubleadAng3d", "std::vector<float>",&_fMultiPairSubleadAng3d);
    _eventtree->Branch("_fMultiPairSubleadStartX","std::vector<float>", &_fMultiPairSubleadStartX);
    _eventtree->Branch("_fMultiPairSubleadStartY","std::vector<float>", &_fMultiPairSubleadStartY);
    _eventtree->Branch("_fMultiPairSubleadStartZ","std::vector<float>", &_fMultiPairSubleadStartZ);
    _eventtree->Branch("_fMultiPairSubleadEnergy0","std::vector<float>", &_fMultiPairSubleadEnergy0);
    _eventtree->Branch("_fMultiPairSubleadEnergy1","std::vector<float>", &_fMultiPairSubleadEnergy1);
    _eventtree->Branch("_fMultiPairSubleadEnergy2","std::vector<float>", &_fMultiPairSubleadEnergy2);
    _eventtree->Branch("_fMultiPairSubleaddEdx0", "std::vector<float>", &_fMultiPairSubleaddEdx0);
    _eventtree->Branch("_fMultiPairSubleaddEdx1", "std::vector<float>", &_fMultiPairSubleaddEdx1);
    _eventtree->Branch("_fMultiPairSubleaddEdx2", "std::vector<float>", &_fMultiPairSubleaddEdx2);
    _eventtree->Branch("_fMultiPairSubleaddEdx3", "std::vector<float>", &_fMultiPairSubleaddEdx3);
    _eventtree->Branch("_fMultiPairSubleadEnergy0Resolution", "std::vector<float>", &_fMultiPairSubleadEnergy0Resolution);
    _eventtree->Branch("_fMultiPairSubleadEnergy1Resolution", "std::vector<float>", &_fMultiPairSubleadEnergy1Resolution);
    _eventtree->Branch("_fMultiPairSubleadEnergy2Resolution", "std::vector<float>", &_fMultiPairSubleadEnergy2Resolution);
    _eventtree->Branch("_fMultiPairSubleadAngleResolution", "std::vector<float>", &_fMultiPairSubleadAngleResolution);
    _eventtree->Branch("_fMultiPairSubleadPurity0", "std::vector<float>", &_fMultiPairSubleadPurity0);
    _eventtree->Branch("_fMultiPairSubleadPurity1", "std::vector<float>", &_fMultiPairSubleadPurity1);
    _eventtree->Branch("_fMultiPairSubleadPurity2", "std::vector<float>", &_fMultiPairSubleadPurity2);
    _eventtree->Branch("_fMultiPairSubleadCompleteness0", "std::vector<float>", &_fMultiPairSubleadCompleteness0);
    _eventtree->Branch("_fMultiPairSubleadCompleteness1", "std::vector<float>", &_fMultiPairSubleadCompleteness1);
    _eventtree->Branch("_fMultiPairSubleadCompleteness2", "std::vector<float>", &_fMultiPairSubleadCompleteness2);

    _eventtree->Branch("_fHiMassPi0Momentum",&_fHiMassPi0Momentum,"_fHiMassPi0Momentum/F");
    _eventtree->Branch("_fHiMassPi0Costheta",&_fHiMassPi0Costheta,"_fHiMassPi0Costheta/F");
    _eventtree->Branch("_fHiMassPi0Phi",&_fHiMassPi0Phi,"_fHiMassPi0Phi/F");
    _eventtree->Branch("_fHiMassPi0Mass",&_fHiMassPi0Mass, "_fHiMassPi0Mass/F");
    _eventtree->Branch("_fHiMassPi0Energy",&_fHiMassPi0Energy,"_fHiMassPi0Energy/F");
    _eventtree->Branch("_fHiMassPi0Angle12",&_fHiMassPi0Angle12,"_fHiMassPi0Angle12/F");
    _eventtree->Branch("_fHiMassPi0MomentumResolution",&_fHiMassPi0MomentumResolution,"_fHiMassPi0MomentumResolution/F");

    _eventtree->Branch("_fHiMassLeadTrueComposition",&_fHiMassLeadTrueComposition,"_fHiMassLeadTrueComposition/I");
    _eventtree->Branch("_fHiMassLeadScore",&_fHiMassLeadScore,"_fHiMassLeadScore/F");
    _eventtree->Branch("_fHiMassLeadLength",&_fHiMassLeadLength,"_fHiMassLeadLength/F");
    _eventtree->Branch("_fHiMassLeadOpenAngle",&_fHiMassLeadOpenAngle,"_fHiMassLeadOpenAngle/F");
    _eventtree->Branch("_fHiMassLeadCostheta",&_fHiMassLeadCostheta,"_fHiMassLeadCostheta/F");
    _eventtree->Branch("_fHiMassLeadPhi",&_fHiMassLeadPhi,"_fHiMassLeadPhi/F");
    _eventtree->Branch("_fHiMassLeadThetaXZ",&_fHiMassLeadThetaXZ,"_fHiMassLeadThetaXZ/F");
    _eventtree->Branch("_fHiMassLeadThetaYZ",&_fHiMassLeadThetaYZ,"_fHiMassLeadThetaYZ/F");
    _eventtree->Branch("_fHiMassLeadDist3d",&_fHiMassLeadDist3d,"_fHiMassLeadDist3d/F");
    _eventtree->Branch("_fHiMassLeadAng3d",&_fHiMassLeadAng3d,"_fHiMassLeadAng3d/F");
    _eventtree->Branch("_fHiMassLeadStartX",&_fHiMassLeadStartX,"_fHiMassLeadStartX/F");
    _eventtree->Branch("_fHiMassLeadStartY",&_fHiMassLeadStartY,"_fHiMassLeadStartY/F");
    _eventtree->Branch("_fHiMassLeadStartZ",&_fHiMassLeadStartZ,"_fHiMassLeadStartZ/F");
    _eventtree->Branch("_fHiMassLeadEnergy0",&_fHiMassLeadEnergy0,"_fHiMassLeadEnergy0/F");
    _eventtree->Branch("_fHiMassLeadEnergy1",&_fHiMassLeadEnergy1,"_fHiMassLeadEnergy1/F");
    _eventtree->Branch("_fHiMassLeadEnergy2",&_fHiMassLeadEnergy2,"_fHiMassLeadEnergy2/F");
    _eventtree->Branch("_fHiMassLeaddEdx0",&_fHiMassLeaddEdx0,"_fHiMassLeaddEdx0/F");
    _eventtree->Branch("_fHiMassLeaddEdx1",&_fHiMassLeaddEdx1,"_fHiMassLeaddEdx1/F");
    _eventtree->Branch("_fHiMassLeaddEdx2",&_fHiMassLeaddEdx2,"_fHiMassLeaddEdx2/F");
    _eventtree->Branch("_fHiMassLeaddEdx3",&_fHiMassLeaddEdx3,"_fHiMassLeaddEdx3/F");
    _eventtree->Branch("_fHiMassLeadEnergy0Resolution",&_fHiMassLeadEnergy0Resolution,"_fHiMassLeadEnergy0Resolution/F");
    _eventtree->Branch("_fHiMassLeadEnergy1Resolution",&_fHiMassLeadEnergy1Resolution,"_fHiMassLeadEnergy1Resolution/F");
    _eventtree->Branch("_fHiMassLeadEnergy2Resolution",&_fHiMassLeadEnergy2Resolution,"_fHiMassLeadEnergy2Resolution/F");
    _eventtree->Branch("_fHiMassLeadAngleResolution",&_fHiMassLeadAngleResolution,"_fHiMassLeadAngleResolution/F");
    _eventtree->Branch("_fHiMassLeadPurity0",&_fHiMassLeadPurity0,"_fHiMassLeadPurity0/F");
    _eventtree->Branch("_fHiMassLeadPurity1",&_fHiMassLeadPurity1,"_fHiMassLeadPurity1/F");
    _eventtree->Branch("_fHiMassLeadPurity2",&_fHiMassLeadPurity2,"_fHiMassLeadPurity2/F");
    _eventtree->Branch("_fHiMassLeadCompleteness0",&_fHiMassLeadCompleteness0,"_fHiMassLeadCompleteness0/F");
    _eventtree->Branch("_fHiMassLeadCompleteness1",&_fHiMassLeadCompleteness1,"_fHiMassLeadCompleteness1/F");
    _eventtree->Branch("_fHiMassLeadCompleteness2",&_fHiMassLeadCompleteness2,"_fHiMassLeadCompleteness2/F");
    

    _eventtree->Branch("_fHiMassSubleadTrueComposition",&_fHiMassSubleadTrueComposition,"_fHiMassSubleadTrueComposition/I");
    _eventtree->Branch("_fHiMassSubleadScore",&_fHiMassSubleadScore,"_fHiMassSubleadScore/F");
    _eventtree->Branch("_fHiMassSubleadLength",&_fHiMassSubleadLength,"_fHiMassSubleadLength/F");
    _eventtree->Branch("_fHiMassSubleadOpenAngle",&_fHiMassSubleadOpenAngle,"_fHiMassSubleadOpenAngle/F");
    _eventtree->Branch("_fHiMassSubleadCostheta",&_fHiMassSubleadCostheta,"_fHiMassSubleadCostheta/F");
    _eventtree->Branch("_fHiMassSubleadPhi",&_fHiMassSubleadPhi,"_fHiMassSubleadPhi/F");
    _eventtree->Branch("_fHiMassSubleadThetaXZ",&_fHiMassSubleadThetaXZ,"_fHiMassSubleadThetaXZ/F");
    _eventtree->Branch("_fHiMassSubleadThetaYZ",&_fHiMassSubleadThetaYZ,"_fHiMassSubleadThetaYZ/F");
    _eventtree->Branch("_fHiMassSubleadDist3d",&_fHiMassSubleadDist3d,"_fHiMassSubleadDist3d/F");
    _eventtree->Branch("_fHiMassSubleadAng3d",&_fHiMassSubleadAng3d,"_fHiMassSubleadAng3d/F");
    _eventtree->Branch("_fHiMassSubleadStartX",&_fHiMassSubleadStartX,"_fHiMassSubleadStartX/F");
    _eventtree->Branch("_fHiMassSubleadStartY",&_fHiMassSubleadStartY,"_fHiMassSubleadStartY/F");
    _eventtree->Branch("_fHiMassSubleadStartZ",&_fHiMassSubleadStartZ,"_fHiMassSubleadStartZ/F");
    _eventtree->Branch("_fHiMassSubleadEnergy0",&_fHiMassSubleadEnergy0,"_fHiMassSubleadEnergy0/F");
    _eventtree->Branch("_fHiMassSubleadEnergy1",&_fHiMassSubleadEnergy1,"_fHiMassSubleadEnergy1/F");
    _eventtree->Branch("_fHiMassSubleadEnergy2",&_fHiMassSubleadEnergy2,"_fHiMassSubleadEnergy2/F");
    _eventtree->Branch("_fHiMassSubleaddEdx0",&_fHiMassSubleaddEdx0,"_fHiMassSubleaddEdx0/F");
    _eventtree->Branch("_fHiMassSubleaddEdx1",&_fHiMassSubleaddEdx1,"_fHiMassSubleaddEdx1/F");
    _eventtree->Branch("_fHiMassSubleaddEdx2",&_fHiMassSubleaddEdx2,"_fHiMassSubleaddEdx2/F");
    _eventtree->Branch("_fHiMassSubleaddEdx3",&_fHiMassSubleaddEdx3,"_fHiMassSubleaddEdx3/F");
    _eventtree->Branch("_fHiMassSubleadEnergy0Resolution",&_fHiMassSubleadEnergy0Resolution,"_fHiMassSubleadEnergy0Resolution/F");
    _eventtree->Branch("_fHiMassSubleadEnergy1Resolution",&_fHiMassSubleadEnergy1Resolution,"_fHiMassSubleadEnergy1Resolution/F");
    _eventtree->Branch("_fHiMassSubleadEnergy2Resolution",&_fHiMassSubleadEnergy2Resolution,"_fHiMassSubleadEnergy2Resolution/F");
    _eventtree->Branch("_fHiMassSubleadAngleResolution",&_fHiMassSubleadAngleResolution,"_fHiMassSubleadAngleResolution/F");
    _eventtree->Branch("_fHiMassSubleadPurity0",&_fHiMassSubleadPurity0,"_fHiMassSubleadPurity0/F");
    _eventtree->Branch("_fHiMassSubleadPurity1",&_fHiMassSubleadPurity1,"_fHiMassSubleadPurity1/F");
    _eventtree->Branch("_fHiMassSubleadPurity2",&_fHiMassSubleadPurity2,"_fHiMassSubleadPurity2/F");
    _eventtree->Branch("_fHiMassSubleadCompleteness0",&_fHiMassSubleadCompleteness0,"_fHiMassSubleadCompleteness0/F");
    _eventtree->Branch("_fHiMassSubleadCompleteness1",&_fHiMassSubleadCompleteness1,"_fHiMassSubleadCompleteness1/F");
    _eventtree->Branch("_fHiMassSubleadCompleteness2",&_fHiMassSubleadCompleteness2,"_fHiMassSubleadCompleteness2/F");


    _eventtree->Branch("_fTrackTrueComposition", "std::vector<int>", &_fTrackTrueComposition);
    _eventtree->Branch("_fTrackIsContained", "std::vector<int>", &_fTrackIsContained);
    _eventtree->Branch("_fTrackDist3d", "std::vector<float>", &_fTrackDist3d);
    _eventtree->Branch("_fTrackAng3d", "std::vector<float>", &_fTrackAng3d);
    _eventtree->Branch("_fTrackPID", "std::vector<float>", &_fTrackPID);
    _eventtree->Branch("_fTrackStartX", "std::vector<float>", &_fTrackStartX);
    _eventtree->Branch("_fTrackStartY", "std::vector<float>", &_fTrackStartY);
    _eventtree->Branch("_fTrackStartZ", "std::vector<float>", &_fTrackStartZ);
    _eventtree->Branch("_fTrackEndX", "std::vector<float>", &_fTrackEndX);
    _eventtree->Branch("_fTrackEndY", "std::vector<float>", &_fTrackEndY);
    _eventtree->Branch("_fTrackEndZ", "std::vector<float>", &_fTrackEndZ);
    _eventtree->Branch("_fTrackMomentum", "std::vector<float>", &_fTrackMomentum);
    _eventtree->Branch("_fTrackMCSMomentum", "std::vector<float>", &_fTrackMCSMomentum);
    _eventtree->Branch("_fTrackRangeMomentum", "std::vector<float>", &_fTrackRangeMomentum);
    _eventtree->Branch("_fTrackLength", "std::vector<float>", &_fTrackLength);
    _eventtree->Branch("_fTrackCostheta", "std::vector<float>", &_fTrackCostheta);
    _eventtree->Branch("_fTrackPhi", "std::vector<float>", &_fTrackPhi);
    _eventtree->Branch("_fTrackThetaXZ", "std::vector<float>", &_fTrackThetaXZ);
    _eventtree->Branch("_fTrackThetaYZ", "std::vector<float>", &_fTrackThetaYZ);

    _eventtree->Branch("_fShowerTrueComposition", "std::vector<int>", &_fShowerTrueComposition);
    _eventtree->Branch("_fTrackShowerScore", "std::vector<float>", &_fTrackShowerScore);
    _eventtree->Branch("_fShowerLength", "std::vector<float>", &_fShowerLength);
    _eventtree->Branch("_fShowerOpenAngle", "std::vector<float>", &_fShowerOpenAngle);
    _eventtree->Branch("_fShowerCostheta", "std::vector<float>", &_fShowerCostheta);
    _eventtree->Branch("_fShowerPhi", "std::vector<float>", &_fShowerPhi);
    _eventtree->Branch("_fShowerThetaXZ", "std::vector<float>", &_fShowerThetaXZ);
    _eventtree->Branch("_fShowerThetaYZ", "std::vector<float>", &_fShowerThetaYZ);
    _eventtree->Branch("_fShowerDist3d", "std::vector<float>", &_fShowerDist3d);
    _eventtree->Branch("_fShowerAng3d", "std::vector<float>", &_fShowerAng3d);
    _eventtree->Branch("_fShowerStartX", "std::vector<float>", &_fShowerStartX);
    _eventtree->Branch("_fShowerStartY", "std::vector<float>", &_fShowerStartY);
    _eventtree->Branch("_fShowerStartZ", "std::vector<float>", &_fShowerStartZ);
    _eventtree->Branch("_fShowerEnergy0", "std::vector<float>", &_fShowerEnergy0);
    _eventtree->Branch("_fShowerEnergy1", "std::vector<float>", &_fShowerEnergy1);
    _eventtree->Branch("_fShowerEnergy2", "std::vector<float>", &_fShowerEnergy2);
    _eventtree->Branch("_fShowerdEdx0", "std::vector<float>", &_fShowerdEdx0);
    _eventtree->Branch("_fShowerdEdx1", "std::vector<float>", &_fShowerdEdx1);
    _eventtree->Branch("_fShowerdEdx2", "std::vector<float>", &_fShowerdEdx2);
    _eventtree->Branch("_fShowerdEdx3", "std::vector<float>", &_fShowerdEdx3);
    _eventtree->Branch("_fShowerdDirX", "std::vector<float>", &_fShowerdDirX);
    _eventtree->Branch("_fShowerdDirY", "std::vector<float>", &_fShowerdDirY);
    _eventtree->Branch("_fShowerdDirZ", "std::vector<float>", &_fShowerdDirZ);
    _eventtree->Branch("_fShowerEnergy0Resolution","std::vector<float>",&_fShowerEnergy0Resolution);
    _eventtree->Branch("_fShowerEnergy1Resolution","std::vector<float>",&_fShowerEnergy1Resolution);
    _eventtree->Branch("_fShowerEnergy2Resolution","std::vector<float>",&_fShowerEnergy2Resolution);
    _eventtree->Branch("_fShowerAngleResolution","std::vector<float>",&_fShowerAngleResolution);
    _eventtree->Branch("_fShowerPurity0","std::vector<float>",&_fShowerPurity0);
    _eventtree->Branch("_fShowerPurity1","std::vector<float>",&_fShowerPurity1);
    _eventtree->Branch("_fShowerPurity2","std::vector<float>",&_fShowerPurity2);
    _eventtree->Branch("_fShowerCompleteness0","std::vector<float>",&_fShowerCompleteness0);
    _eventtree->Branch("_fShowerCompleteness1","std::vector<float>",&_fShowerCompleteness1);
    _eventtree->Branch("_fShowerCompleteness2","std::vector<float>",&_fShowerCompleteness2);
    _eventtree->Branch("_fLeadShowerIDX", &_fLeadShowerIDX, "_fLeadShowerIDX/I");
    _eventtree->Branch("_fSubLeadShowerIDX", &_fSubLeadShowerIDX, "_fSubLeadShowerIDX/I");
    
    _eventtree->Branch("_fHasCandidateNeutrino",&_fHasCandidateNeutrino,"_fHasCandidateNeutrino/I");
    _eventtree->Branch("_fHasCandidateMuon",&_fHasCandidateMuon,"_fHasCandidateMuon/I");
    _eventtree->Branch("_fNChargedPiCandidates",&_fNChargedPiCandidates,"_fNChargedPiCandidates/I");
    _eventtree->Branch("_fNProtonCandidates",&_fNProtonCandidates,"_fNProtonCandidates/I");
    _eventtree->Branch("_fNLeadCandidates",&_fNLeadCandidates,"_fNLeadCandidates/I");
    _eventtree->Branch("_fNSubleadCandidates",&_fNSubleadCandidates,"_fNSubleadCandidates/I");
    _eventtree->Branch("_fNPi0Candidates",&_fNPi0Candidates,"_fNPi0Candidates/I");
    _eventtree->Branch("_fNPairCandidates",&_fNPairCandidates,"_fNPairCandidates/I");
    _eventtree->Branch("_fNOtherFancyPairs",&_fNOtherFancyPairs,"_fNOtherFancyPairs/I");
    _eventtree->Branch("_fNHiMassPairs", &_fNHiMassPairs, "_fNHiMassPairs/I");
    
  }
  
}//Begin job

void SidebandTree::analyze(art::Event const & e){

  // Get the MC information    
  art::Handle<std::vector<simb::MCTruth>> MCTruthHandle;
  std::vector<art::Ptr<simb::MCTruth>> MCTruthVector;
  if (e.getByLabel(fMCTruthTag, MCTruthHandle)){art::fill_ptr_vector(MCTruthVector, MCTruthHandle);}

  art::Handle<std::vector<simb::GTruth>> GTruthHandle;
  std::vector<art::Ptr<simb::GTruth>> GTruthVector;
  if (e.getByLabel(fGTruthTag, GTruthHandle)){art::fill_ptr_vector(GTruthVector, GTruthHandle);}
    
  art::Handle<std::vector<simb::MCParticle>> MCParticleHandle;
  std::vector<art::Ptr<simb::MCParticle>> MCParticleVector;
  if (e.getByLabel(fMCParticleTag, MCParticleHandle)){art::fill_ptr_vector(MCParticleVector, MCParticleHandle);}

  //Get the other reco products
  art::Handle<std::vector<recob::PFParticle>> PFParticleHandle;
  std::vector<art::Ptr<recob::PFParticle>> PFParticleVector;
  if (e.getByLabel(fPFParticleTag, PFParticleHandle)){art::fill_ptr_vector(PFParticleVector, PFParticleHandle);}

  art::Handle<std::vector<sim::MCShower>> MCShowerHandle;
  std::vector<art::Ptr<sim::MCShower>> MCShowerVector;
  if (e.getByLabel(fMCShowerTag, MCShowerHandle)){art::fill_ptr_vector(MCShowerVector, MCShowerHandle);}

  art::Handle<std::vector<recob::Shower>> ShowerHandle;
  std::vector<art::Ptr<recob::Shower>> ShowerVector;
  if (e.getByLabel(fShowerTag, ShowerHandle)){art::fill_ptr_vector(ShowerVector, ShowerHandle);}

  art::Handle<std::vector<recob::Cluster>> ClusterHandle;
  std::vector<art::Ptr<recob::Cluster>> ClusterVector;
  if (e.getByLabel(fClusterTag, ClusterHandle)){art::fill_ptr_vector(ClusterVector, ClusterHandle);}

  art::Handle<std::vector<recob::Track>> TrackHandle;
  std::vector<art::Ptr<recob::Track>> TrackVector;
  if (e.getByLabel(fTrackTag, TrackHandle)){art::fill_ptr_vector(TrackVector, TrackHandle);}

  art::Handle<std::vector<recob::Track>> TrackFitterHandle;
  std::vector<art::Ptr<recob::Track>> TrackFitterVector;
  if (e.getByLabel(fTrackFitterTag, TrackFitterHandle)){art::fill_ptr_vector(TrackFitterVector, TrackFitterHandle);}

  art::Handle<std::vector<recob::Hit>> HitHandle;
  std::vector<art::Ptr<recob::Hit>> HitVector;
  if (e.getByLabel(fHitTag, HitHandle)){art::fill_ptr_vector(HitVector, HitHandle);}

  art::Handle<std::vector<recob::MCSFitResult> > MCSFitResultHandle;
  std::vector<art::Ptr<recob::MCSFitResult>> MCSFitResultVector;
  if(e.getByLabel(fMCSFitResultTag, MCSFitResultHandle)){art::fill_ptr_vector(MCSFitResultVector, MCSFitResultHandle);}
    
  art::FindManyP<recob::Hit>  fmhs(ShowerHandle, e, fShowerTag);
  art::FindManyP<recob::Hit>  fmht(TrackHandle, e, fTrackTag);
  art::FindManyP<recob::Hit>  fmhc(ClusterHandle, e, fClusterTag);
  art::FindManyP<recob::Cluster>  fmcs(ShowerHandle, e, fShowerTag);
  art::FindManyP<recob::Shower> pfPartToShowerAssoc(PFParticleHandle, e, "pandora");
  art::FindManyP<recob::Cluster> pfPartToClusterAssoc(PFParticleHandle, e, "pandora");
  art::FindManyP<recob::Track> trackFitterAssoc(PFParticleHandle, e, "pandoraKalmanShower");
  std::vector<art::Ptr<recob::Track> > nuTracks;
  std::map< art::Ptr<recob::Shower> , float > nuShowers_ScoreMap; //map from shower to score
  std::map< art::Ptr<recob::Track> , float > nuTracks_ScoreMap; //map from track to score
  std::map< art::Ptr<recob::Shower> , int > nuShowerTrackMap;
  art::FindManyP<anab::T0> nuFlashScoreAssoc(PFParticleHandle, e, "flashmatch");

  //CUT VARIABLES
  //cut 1: cc inclusive filter
  std::vector<art::Ptr<recob::Vertex> > nuVertex;
  std::vector<art::Ptr<recob::Track> > nuMuon;
  //cut 2: charged pi veto
  int nchargedpicand = 0;
  //cut 3: nshowers >= 2
  std::vector<art::Ptr<recob::Shower> > nuShowers;
  //cut 4: >= 1 lead
  std::vector<art::Ptr<recob::Shower>> TheLeadVector;
  std::vector<art::Ptr<recob::Shower>> LeadCandidates; //all showers that pass lead requirements
  std::vector<art::Ptr<recob::Shower>> LeadCandidatesCL;
  std::vector<art::Ptr<recob::Shower>> LeadCandidatesEN;
  std::vector<art::Ptr<recob::Shower>> LeadCandidatesRA;
  //cut 5: >= 1 sublead
  std::vector<art::Ptr<recob::Shower>> TheSubleadVector;
  std::vector<art::Ptr<recob::Shower>> SubleadVectorHighMass;
  std::vector<art::Ptr<recob::Shower>> SubleadCandidates; //all showers that pass sublead requirements
  std::vector<art::Ptr<recob::Shower>> SubleadCandidatesCL;
  std::vector<art::Ptr<recob::Shower>> SubleadCandidatesEN;
  //cut 6: exactly one pair passes mass cut
  int npairs = 0;
  int npairsHiMass = 0;
  int notherfancypairs = 0;

  //CLEAR TREE VARIABLES
  _f_All_Genie_Weights.clear();
  _f_RPA_CCQE_Genie_Weights.clear();
  _f_XSecShape_CCMEC_Genie_Weights.clear();
  _f_AxFFCCQEshape_Genie_Weights.clear();
  _f_VecFFCCQEshape_Genie_Weights.clear();
  _f_DecayAngMEC_Genie_Weights.clear();
  _f_Theta_Delta2Npi_Genie_Weights.clear();
  _f_expskin_FluxUnisim_Weights.clear();
  _f_horncurrent_FluxUnisim_Weights.clear();
  _f_kminus_PrimaryHadronNormalization_Weights.clear();
  _f_kplus_PrimaryHadronFeynmanScaling_Weights.clear();
  _f_kzero_PrimaryHadronSanfordWang_Weights.clear();
  _f_nucleoninexsec_FluxUnisim_Weights.clear();
  _f_nucleonqexsec_FluxUnisim_Weights.clear();
  _f_nucleontotxsec_FluxUnisim_Weights.clear();
  _f_piminus_PrimaryHadronSWCentralSplineVariation_Weights.clear();
  _f_pioninexsec_FluxUnisim_Weights.clear();
  _f_pionqexsec_FluxUnisim_Weights.clear();
  _f_piontotxsec_FluxUnisim_Weights.clear();
  _f_piplus_PrimaryHadronSWCentralSplineVariation_Weights.clear();
  
  _fCVWeight = 1.;

  //All tracks
  _fTrackTrueComposition.clear();
  _fTrackIsContained.clear();
  _fTrackDist3d.clear();
  _fTrackAng3d.clear();
  _fTrackPID.clear();
  _fTrackStartX.clear();
  _fTrackStartY.clear();
  _fTrackStartZ.clear();
  _fTrackEndX.clear();
  _fTrackEndY.clear();
  _fTrackEndZ.clear();
  _fTrackMomentum.clear();
  _fTrackMCSMomentum.clear();
  _fTrackRangeMomentum.clear();
  _fTrackLength.clear();
  _fTrackCostheta.clear();
  _fTrackPhi.clear();
  _fTrackThetaXZ.clear();
  _fTrackThetaYZ.clear();

  //All showers
  _fShowerTrueComposition.clear();
  _fTrackShowerScore.clear();
  _fShowerLength.clear();
  _fShowerOpenAngle.clear();
  _fShowerCostheta.clear();
  _fShowerPhi.clear();
  _fShowerThetaXZ.clear();
  _fShowerThetaYZ.clear();
  _fShowerDist3d.clear();
  _fShowerAng3d.clear();
  _fShowerStartX.clear();
  _fShowerStartY.clear();
  _fShowerStartZ.clear();
  _fShowerEnergy0.clear();
  _fShowerEnergy1.clear();
  _fShowerEnergy2.clear();
  _fShowerdEdx0.clear();
  _fShowerdEdx1.clear();
  _fShowerdEdx2.clear();
  _fShowerdEdx3.clear();
  _fShowerdDirX.clear();
  _fShowerdDirY.clear();
  _fShowerdDirZ.clear();
  _fShowerEnergy0Resolution.clear();
  _fShowerEnergy1Resolution.clear();
  _fShowerEnergy2Resolution.clear();
  _fShowerAngleResolution.clear();
  _fShowerPurity0.clear();
  _fShowerPurity1.clear();
  _fShowerPurity2.clear();
  _fShowerCompleteness0.clear();
  _fShowerCompleteness1.clear();
  _fShowerCompleteness2.clear();
  _fTwoPhotonInvariantMass.clear();
  _fLeadShowerIDX = -999;
  _fSubLeadShowerIDX = -999;

  //////
  int MuonTrueComposition = -999;
  int MuonIsContained = -999;
  float MuonDist3d = -999;
  //float MuonAng3d = -999;
  float VertexX = -999;
  float VertexY = -999;
  float VertexZ = -999;
  float MuonStartX = -999;
  float MuonStartY = -999;
  float MuonStartZ = -999;
  float MuonEndX = -999;
  float MuonEndY = -999;
  float MuonEndZ = -999;
  float MuonMomentum = -999;
  float MuonMCSMomentum = -999;
  float MuonRangeMomentum = -999;
  float MuonLength = -999;
  float MuonCostheta = -999;
  float MuonPhi = -999;
  int MuonID = -999;
  float MuonPID = -999999;

  //Signal Pion Vars
  _fCandidatePi0Momentum = -999.9;
  _fCandidatePi0Costheta = -999.9;
  _fCandidatePi0Phi = -999.9;
  _fCandidatePi0Mass = -999.9;
  _fCandidatePi0Energy = -999.9;;
  _fCandidatePi0Angle12 = -999.9;;
  _fCandidatePi0MomentumResolution = -999.9;;
  _fCandidatePi0AngleResolution = -999.9;

   //Candidate photons
  _fCandidateLeadTrueComposition = -999;
  _fCandidateLeadScore = -999.9;
  _fCandidateLeadLength = -999.9;
  _fCandidateLeadOpenAngle = -999.9;
  _fCandidateLeadCostheta = -999.9;
  _fCandidateLeadPhi = -999.9;
  _fCandidateLeadThetaXZ = -999.9;
  _fCandidateLeadThetaYZ = -999.9;
  _fCandidateLeadDist3d = -999.9;
  _fCandidateLeadAng3d = -999.9;
  _fCandidateLeadStartX = -999.9;
  _fCandidateLeadStartY = -999.9;
  _fCandidateLeadStartZ = -999.9;
  _fCandidateLeadEnergy0 = -999.9;
  _fCandidateLeadEnergy1 = -999.9;
  _fCandidateLeadEnergy2 = -999.9;
  _fCandidateLeaddEdx0 = -999.9;
  _fCandidateLeaddEdx1 = -999.9;
  _fCandidateLeaddEdx2 = -999.9;
  _fCandidateLeaddEdx3 = -999.9;
  _fCandidateLeadEnergy0Resolution = -999.9;
  _fCandidateLeadEnergy1Resolution = -999.9;
  _fCandidateLeadEnergy2Resolution = -999.9;
  _fCandidateLeadAngleResolution = -999.9;
  _fCandidateLeadPurity0 = -999.9;
  _fCandidateLeadPurity1 = -999.9;
  _fCandidateLeadPurity2 = -999.9;
  _fCandidateLeadCompleteness0 = -999.9;
  _fCandidateLeadCompleteness1 = -999.9;
  _fCandidateLeadCompleteness2 = -999.9;

  _fCandidateSubleadTrueComposition = -999;
  _fCandidateSubleadScore = -999.9;
  _fCandidateSubleadLength = -999.9;
  _fCandidateSubleadOpenAngle = -999.9;
  _fCandidateSubleadCostheta = -999.9;
  _fCandidateSubleadPhi = -999.9;
  _fCandidateSubleadThetaXZ = -999.9;
  _fCandidateSubleadThetaYZ = -999.9;
  _fCandidateSubleadDist3d = -999.9;
  _fCandidateSubleadAng3d = -999.9;
  _fCandidateSubleadStartX = -999.9;
  _fCandidateSubleadStartY = -999.9;
  _fCandidateSubleadStartZ = -999.9;
  _fCandidateSubleadEnergy0 = -999.9;
  _fCandidateSubleadEnergy1 = -999.9;
  _fCandidateSubleadEnergy2 = -999.9;
  _fCandidateSubleaddEdx0 = -999.9;
  _fCandidateSubleaddEdx1 = -999.9;
  _fCandidateSubleaddEdx2 = -999.9;
  _fCandidateSubleaddEdx3 = -999.9;
  _fCandidateSubleadEnergy0Resolution = -999.9;
  _fCandidateSubleadEnergy1Resolution = -999.9;
  _fCandidateSubleadEnergy2Resolution = -999.9;
  _fCandidateSubleadAngleResolution = -999.9;
  _fCandidateSubleadPurity0 = -999.9;
  _fCandidateSubleadPurity1 = -999.9;
  _fCandidateSubleadPurity2 = -999.9;
  _fCandidateSubleadCompleteness0 = -999.9;
  _fCandidateSubleadCompleteness1 = -999.9;
  _fCandidateSubleadCompleteness2 = -999.9;
 
  //2 Mip Pion vars
  _fTwoMIPPi0Momentum = -999.9;
  _fTwoMIPPi0Costheta = -999.9;
  _fTwoMIPPi0Phi = -999.9;
  _fTwoMIPPi0Mass = -999.9;
  _fTwoMIPPi0Energy = -999.9;
  _fTwoMIPPi0Angle12 = -999.9;
  _fTwoMIPPi0MomentumResolution = -999.9;
  _fTwoMIPPi0AngleResolution = -999.9;

  //2 MIP Sideband photons
  _fTwoMIPLeadTrueComposition = -999;
  _fTwoMIPLeadScore = -999.9;
  _fTwoMIPLeadLength = -999.9;
  _fTwoMIPLeadOpenAngle = -999.9;
  _fTwoMIPLeadCostheta = -999.9;
  _fTwoMIPLeadPhi = -999.9;
  _fTwoMIPLeadThetaXZ = -999.9;
  _fTwoMIPLeadThetaYZ = -999.9;
  _fTwoMIPLeadDist3d = -999.9;
  _fTwoMIPLeadAng3d = -999.9;
  _fTwoMIPLeadStartX = -999.9;
  _fTwoMIPLeadStartY = -999.9;
  _fTwoMIPLeadStartZ = -999.9;
  _fTwoMIPLeadEnergy0 = -999.9;
  _fTwoMIPLeadEnergy1 = -999.9;
  _fTwoMIPLeadEnergy2 = -999.9;
  _fTwoMIPLeaddEdx0 = -999.9;
  _fTwoMIPLeaddEdx1 = -999.9;
  _fTwoMIPLeaddEdx2 = -999.9;
  _fTwoMIPLeaddEdx3 = -999.9;
  _fTwoMIPLeadEnergy0Resolution = -999.9;
  _fTwoMIPLeadEnergy1Resolution = -999.9;
  _fTwoMIPLeadEnergy2Resolution = -999.9;
  _fTwoMIPLeadAngleResolution = -999.9;
  _fTwoMIPLeadPurity0 = -999.9;
  _fTwoMIPLeadPurity1 = -999.9;
  _fTwoMIPLeadPurity2 = -999.9;
  _fTwoMIPLeadCompleteness0 = -999.9;
  _fTwoMIPLeadCompleteness1 = -999.9;
  _fTwoMIPLeadCompleteness2 = -999.9;

  _fTwoMIPSubleadTrueComposition = -999;
  _fTwoMIPSubleadScore = -999.9;
  _fTwoMIPSubleadLength = -999.9;
  _fTwoMIPSubleadOpenAngle = -999.9;
  _fTwoMIPSubleadCostheta = -999.9;
  _fTwoMIPSubleadPhi = -999.9;
  _fTwoMIPSubleadThetaXZ = -999.9;
  _fTwoMIPSubleadThetaYZ = -999.9;
  _fTwoMIPSubleadDist3d = -999.9;
  _fTwoMIPSubleadAng3d = -999.9;
  _fTwoMIPSubleadStartX = -999.9;
  _fTwoMIPSubleadStartY = -999.9;
  _fTwoMIPSubleadStartZ = -999.9;
  _fTwoMIPSubleadEnergy0 = -999.9;
  _fTwoMIPSubleadEnergy1 = -999.9;
  _fTwoMIPSubleadEnergy2 = -999.9;
  _fTwoMIPSubleaddEdx0 = -999.9;
  _fTwoMIPSubleaddEdx1 = -999.9;
  _fTwoMIPSubleaddEdx2 = -999.9;
  _fTwoMIPSubleaddEdx3 = -999.9;
  _fTwoMIPSubleadEnergy0Resolution = -999.9;
  _fTwoMIPSubleadEnergy1Resolution = -999.9;
  _fTwoMIPSubleadEnergy2Resolution = -999.9;
  _fTwoMIPSubleadAngleResolution = -999.9;
  _fTwoMIPSubleadPurity0 = -999.9;
  _fTwoMIPSubleadPurity1 = -999.9;
  _fTwoMIPSubleadPurity2 = -999.9;
  _fTwoMIPSubleadCompleteness0 = -999.9;
  _fTwoMIPSubleadCompleteness1 = -999.9;
  _fTwoMIPSubleadCompleteness2 = -999.9;

  //Multi pairs
  _fMultiPairPi0Momentum.clear();
  _fMultiPairPi0Costheta.clear();
  _fMultiPairPi0Phi.clear();
  _fMultiPairPi0Mass.clear();
  _fMultiPairPi0Energy.clear();
  _fMultiPairPi0Angle12.clear();
  _fMultiPairPi0MomentumResolution.clear();
  _fMultiPairPi0AngleResolution.clear();

  //Multi pairs leading photons
  _fMultiPairLeadTrueComposition = -999;
  _fMultiPairLeadScore = -999.9;
  _fMultiPairLeadLength = -999.9;
  _fMultiPairLeadOpenAngle = -999.9;
  _fMultiPairLeadCostheta = -999.9;
  _fMultiPairLeadPhi = -999.9;
  _fMultiPairLeadThetaXZ = -999.9;
  _fMultiPairLeadThetaYZ = -999.9;
  _fMultiPairLeadDist3d = -999.9;
  _fMultiPairLeadAng3d = -999.9;
  _fMultiPairLeadStartX = -999.9;
  _fMultiPairLeadStartY = -999.9;
  _fMultiPairLeadStartZ = -999.9;
  _fMultiPairLeadEnergy0 = -999.9;
  _fMultiPairLeadEnergy1 = -999.9;
  _fMultiPairLeadEnergy2 = -999.9;
  _fMultiPairLeaddEdx0 = -999.9;
  _fMultiPairLeaddEdx1 = -999.9;
  _fMultiPairLeaddEdx2 = -999.9;
  _fMultiPairLeaddEdx3 = -999.9;
  _fMultiPairLeadEnergy0Resolution = -999.9;
  _fMultiPairLeadEnergy1Resolution = -999.9;
  _fMultiPairLeadEnergy2Resolution = -999.9;
  _fMultiPairLeadAngleResolution = -999.9;
  _fMultiPairLeadPurity0 = -999.9;
  _fMultiPairLeadPurity1 = -999.9;
  _fMultiPairLeadPurity2 = -999.9;
  _fMultiPairLeadCompleteness0 = -999.9;
  _fMultiPairLeadCompleteness1 = -999.9;
  _fMultiPairLeadCompleteness2 = -999.9;

  _fMultiPairSubleadTrueComposition.clear();
  _fMultiPairSubleadScore.clear();
  _fMultiPairSubleadLength.clear();
  _fMultiPairSubleadOpenAngle.clear();
  _fMultiPairSubleadCostheta.clear();
  _fMultiPairSubleadPhi.clear();
  _fMultiPairSubleadThetaXZ.clear();
  _fMultiPairSubleadThetaYZ.clear();
  _fMultiPairSubleadDist3d.clear();
  _fMultiPairSubleadAng3d.clear();
  _fMultiPairSubleadStartX.clear();
  _fMultiPairSubleadStartY.clear();
  _fMultiPairSubleadStartZ.clear();
  _fMultiPairSubleadEnergy0.clear();
  _fMultiPairSubleadEnergy1.clear();
  _fMultiPairSubleadEnergy2.clear();
  _fMultiPairSubleaddEdx0.clear();
  _fMultiPairSubleaddEdx1.clear();
  _fMultiPairSubleaddEdx2.clear();
  _fMultiPairSubleaddEdx3.clear();
  _fMultiPairSubleadEnergy0Resolution.clear();
  _fMultiPairSubleadEnergy1Resolution.clear();
  _fMultiPairSubleadEnergy2Resolution.clear();
  _fMultiPairSubleadAngleResolution.clear();
  _fMultiPairSubleadPurity0.clear();
  _fMultiPairSubleadPurity1.clear();
  _fMultiPairSubleadPurity2.clear();
  _fMultiPairSubleadCompleteness0.clear();
  _fMultiPairSubleadCompleteness1.clear();
  _fMultiPairSubleadCompleteness2.clear();

  //Hi Mass Pion vars
  _fHiMassPi0Momentum = -999.9;
  _fHiMassPi0Costheta = -999.9;
  _fHiMassPi0Phi = -999.9;
  _fHiMassPi0Mass = -999.9;
  _fHiMassPi0Energy = -999.9;
  _fHiMassPi0Angle12 = -999.9;
  _fHiMassPi0MomentumResolution = -999.9;
  _fHiMassPi0AngleResolution = -999.9;

  //Hi Mass Sideband photons
  _fHiMassLeadTrueComposition = -999;
  _fHiMassLeadScore = -999.9;
  _fHiMassLeadLength = -999.9;
  _fHiMassLeadOpenAngle = -999.9;
  _fHiMassLeadCostheta = -999.9;
  _fHiMassLeadPhi = -999.9;
  _fHiMassLeadThetaXZ = -999.9;
  _fHiMassLeadThetaYZ = -999.9;
  _fHiMassLeadDist3d = -999.9;
  _fHiMassLeadAng3d = -999.9;
  _fHiMassLeadStartX = -999.9;
  _fHiMassLeadStartY = -999.9;
  _fHiMassLeadStartZ = -999.9;
  _fHiMassLeadEnergy0 = -999.9;
  _fHiMassLeadEnergy1 = -999.9;
  _fHiMassLeadEnergy2 = -999.9;
  _fHiMassLeaddEdx0 = -999.9;
  _fHiMassLeaddEdx1 = -999.9;
  _fHiMassLeaddEdx2 = -999.9;
  _fHiMassLeaddEdx3 = -999.9;
  _fHiMassLeadEnergy0Resolution = -999.9;
  _fHiMassLeadEnergy1Resolution = -999.9;
  _fHiMassLeadEnergy2Resolution = -999.9;
  _fHiMassLeadAngleResolution = -999.9;
  _fHiMassLeadPurity0 = -999.9;
  _fHiMassLeadPurity1 = -999.9;
  _fHiMassLeadPurity2 = -999.9;
  _fHiMassLeadCompleteness0 = -999.9;
  _fHiMassLeadCompleteness1 = -999.9;
  _fHiMassLeadCompleteness2 = -999.9;

  _fHiMassSubleadTrueComposition = -999;
  _fHiMassSubleadScore = -999.9;
  _fHiMassSubleadLength = -999.9;
  _fHiMassSubleadOpenAngle = -999.9;
  _fHiMassSubleadCostheta = -999.9;
  _fHiMassSubleadPhi = -999.9;
  _fHiMassSubleadThetaXZ = -999.9;
  _fHiMassSubleadThetaYZ = -999.9;
  _fHiMassSubleadDist3d = -999.9;
  _fHiMassSubleadAng3d = -999.9;
  _fHiMassSubleadStartX = -999.9;
  _fHiMassSubleadStartY = -999.9;
  _fHiMassSubleadStartZ = -999.9;
  _fHiMassSubleadEnergy0 = -999.9;
  _fHiMassSubleadEnergy1 = -999.9;
  _fHiMassSubleadEnergy2 = -999.9;
  _fHiMassSubleaddEdx0 = -999.9;
  _fHiMassSubleaddEdx1 = -999.9;
  _fHiMassSubleaddEdx2 = -999.9;
  _fHiMassSubleaddEdx3 = -999.9;
  _fHiMassSubleadEnergy0Resolution = -999.9;
  _fHiMassSubleadEnergy1Resolution = -999.9;
  _fHiMassSubleadEnergy2Resolution = -999.9;
  _fHiMassSubleadAngleResolution = -999.9;
  _fHiMassSubleadPurity0 = -999.9;
  _fHiMassSubleadPurity1 = -999.9;
  _fHiMassSubleadPurity2 = -999.9;
  _fHiMassSubleadCompleteness0 = -999.9;
  _fHiMassSubleadCompleteness1 = -999.9;
  _fHiMassSubleadCompleteness2 = -999.9;

  _fFlashChi2 = -999.9;

  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// EVENT TRUTH INFO ///////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  int Npi0 = 0;
  int Nmuon = 0;
  int Npiplus = 0;
  int Nproton = 0;
  int Nneutron = 0;
  
  std::vector<art::Ptr<simb::MCParticle>> NuPi0;
  std::vector<art::Ptr<simb::MCParticle>> TwoGammaV;
  std::vector<art::Ptr<simb::MCParticle>> Gamma1_descendents;
  std::vector<art::Ptr<simb::MCParticle>> Gamma2_descendents;
    
  if(fIsMC){
      
    if(MCTruthVector[0]->NeutrinoSet()){ //neutrino information
      _fNuPDG = MCTruthVector[0]->GetNeutrino().Nu().PdgCode();
      _fNuCCNC = MCTruthVector[0]->GetNeutrino().CCNC();
      _fNuVtxX = MCTruthVector[0]->GetNeutrino().Nu().Vx();
      _fNuVtxY = MCTruthVector[0]->GetNeutrino().Nu().Vy();
      _fNuVtxZ = MCTruthVector[0]->GetNeutrino().Nu().Vz();
      _fNuEnergy =  1000*(MCTruthVector[0]->GetNeutrino().Nu().E());
      _fNuMode = MCTruthVector[0]->GetNeutrino().Mode();
      if(inFV(_fNuVtxX,_fNuVtxY,_fNuVtxZ) == true){_fNuInFV = 1;}
      if(inFV(_fNuVtxX,_fNuVtxY,_fNuVtxZ) == false){_fNuInFV = 0;}
  
      _fLepPdg = MCTruthVector[0]->GetNeutrino().Lepton().PdgCode();
      if(_fLepPdg == 13 || _fLepPdg == -13){
  _fLepP = (MCTruthVector[0]->GetNeutrino().Lepton().P())*1000;
  _fLepPx = (MCTruthVector[0]->GetNeutrino().Lepton().Px())*1000;
  _fLepPy = (MCTruthVector[0]->GetNeutrino().Lepton().Py())*1000;
  _fLepPz = (MCTruthVector[0]->GetNeutrino().Lepton().Pz())*1000;
  _fLepCosTheta = _fLepPz/_fLepP;
  _fLepPhi = getPhi(_fLepPx, _fLepPy, _fLepPz);       
      }

      if(GTruthVector.size() == 1){
  _fResNum = GTruthVector[0]->fResNum;
      }
  
    }//neutrino set
      
    ////// HOW MANY TRUE PI0'S? //////
           
    for(auto const& mcpart : MCParticleVector){ //list of mc particles

      if(std::abs(mcpart->PdgCode()) == 111 && mcpart->Process() == "primary" && mcpart->P()*1000 < 800.){ 
  Npi0++;
  NuPi0.push_back(mcpart);
      }
      if(std::abs(mcpart->PdgCode()) == 13 && mcpart->Process() == "primary" && mcpart->P()*1000 > 150.){ 
  Nmuon++;
      }
      if(std::abs(mcpart->PdgCode()) == 211 && mcpart->Process() == "primary" && mcpart->P()*1000 > 70.){ 
  Npiplus++;
      }
      if(std::abs(mcpart->PdgCode()) == 2212 && mcpart->Process() == "primary"){ 
  Nproton++;
      }
      if(std::abs(mcpart->PdgCode()) == 2112 && mcpart->Process() == "primary"){ 
  Nneutron++;
      }
    }//mcparticles

    ////// TRUE SIGNAL EVENTS //////
    if(_fNuInFV == 1 && std::abs(_fNuPDG) == 14 && _fNuCCNC == 0 && Npi0 == 1 && Nmuon == 1 && Npiplus == 0){
   
      if(NuPi0.size() == 1){
  _fPi0P = NuPi0.at(0)->P()*1000;
  _fPi0E = NuPi0.at(0)->E()*1000;
  _fPi0Px = NuPi0.at(0)->Px();
  _fPi0Py = NuPi0.at(0)->Py();
  _fPi0Pz = NuPi0.at(0)->Pz();
  _fPi0CosTheta = (_fPi0Pz*1000)/_fPi0P;
  _fPi0Phi = getPhi(_fPi0Px, _fPi0Py, _fPi0Pz);
      }
  
      //get pi0 daughters (gammas)
      auto const nupi0 = NuPi0.at(0);
      for(int i = 0; i < nupi0->NumberDaughters(); i++){
  for(auto const& mcpart : MCParticleVector){
    if(mcpart->TrackId() != nupi0->Daughter(i)) continue; 
    TwoGammaV.push_back(mcpart);
  }
      }

      //get higher and lower energy gammas
      float Gamma1_tempe = -9999.;
      float Gamma2_tempe = -9999.;
  
      for(auto const& gamma : TwoGammaV){
  if(gamma->E() > Gamma1_tempe){
    Gamma1_tempe = gamma->E();
    _fGamma1_id = gamma->TrackId();
  }
      }
      for(auto const& gamma : TwoGammaV){
  if(gamma->TrackId() == _fGamma1_id) continue;
  _fGamma2_id = gamma->TrackId();
  Gamma2_tempe = gamma->E();
      }

      _fGamma1_E = Gamma1_tempe*1000;
      _fGamma2_E = Gamma2_tempe*1000;
  
      for(auto const& gamma : TwoGammaV){

  //Gamma1
  if(gamma->TrackId() == _fGamma1_id){
    if(gamma->NumberDaughters() == 2){
      _fGamma1E1_id = gamma->Daughter(0);
      _fGamma1E2_id = gamma->Daughter(1);
    }//2 daughters
  }//gamma1

  //Gamma2
  if(gamma->TrackId() == _fGamma2_id){
    if(gamma->NumberDaughters() == 2){
      _fGamma2E1_id = gamma->Daughter(0);
      _fGamma2E2_id = gamma->Daughter(1);
    }//2 daughters
  }//gamma1
    
      }
      //////////////// ALL DESCENDENTS //////////////
  
      getKids(MCParticleVector,_fGamma1_id,Gamma1_descendents);
      getKids(MCParticleVector,_fGamma2_id,Gamma2_descendents);

      for(auto const& mcpart : TwoGammaV){
    
  if(mcpart->TrackId() == _fGamma1_id){
    _fGamma1_DirX = mcpart->Px()/mcpart->P();
    _fGamma1_DirY = mcpart->Py()/mcpart->P();
    _fGamma1_DirZ = mcpart->Pz()/mcpart->P();
  }
    
  if(mcpart->TrackId() == _fGamma2_id){
    _fGamma2_DirX = mcpart->Px()/mcpart->P();
    _fGamma2_DirY = mcpart->Py()/mcpart->P();
    _fGamma2_DirZ = mcpart->Pz()/mcpart->P();
  }
    
      }

      _fGamma12_CosTheta = _fGamma1_DirX*_fGamma2_DirX + _fGamma1_DirY*_fGamma2_DirY + _fGamma1_DirZ*_fGamma2_DirZ;
      _fGamma12_Angle = radToDeg(std::acos(_fGamma12_CosTheta));

      _fPi0M = pi0Mass(_fGamma1_E,_fGamma2_E,_fGamma12_CosTheta);  

    }// end of is true signal
      
  }//fIsMC

  _fNpi0 = Npi0;
  _fNmuon = Nmuon;
  _fNpiplus = Npiplus;
  _fNproton = Nproton;
  _fNneutron = Nneutron;
      
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////// EVENT WEIGHTS //////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double cvweight_only = 1.;
  double splineweight = 1.;
  double cvweight = 1.;
  //double avgdirtweight = 1.;
      
  if(fIsMC){

    art::Handle<std::vector<evwgh::MCEventWeight>> EventWeightHandle;
    e.getByLabel(fEventWeightTag, EventWeightHandle);
    std::vector<art::Ptr<evwgh::MCEventWeight>> EventWeightVector;
    art::fill_ptr_vector(EventWeightVector, EventWeightHandle);
     
    if (EventWeightVector.size() > 0) {
      art::Ptr<evwgh::MCEventWeight> evt_wgt = EventWeightVector.at(0);

      for (auto entry : evt_wgt->fWeight) {

  std::cout<<"WEIGHT MAPS! "<<entry.first<<std::endl;
    
  if(entry.first == "TunedCentralValue_Genie"){
    if(entry.second.at(0) >= 0 && entry.second.at(0) <= 100000.){cvweight_only = entry.second.at(0);}
  }//cv

  if(entry.first == "splines_general_Spline"){
    if(entry.second.at(0) >= 0 && entry.second.at(0) <= 100000.){splineweight = entry.second.at(0);}
  }//spline

      }//eventweight

      cvweight = splineweight*cvweight_only;
      if(_fResNum == 17 || _fResNum == 9){cvweight = 0;}
      _fCVWeight = cvweight;

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////// STORE EVENT WEIGHTS FOR SYSTEMATIC UNCERTAINTIES /////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      
      for (auto entry : evt_wgt->fWeight) {

  if(_fResNum == 17 || _fResNum == 9) continue;

  ///////////////// cross-section knobs /////////////////////////////

  //1
  if(entry.first == "All_Genie"){
    for(int i = 0; i < 100; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_All_Genie_Weights.push_back(sysweight);
    }
  }

  //2
  if(entry.first == "RPA_CCQE_Genie"){
    for(int i = 0; i < 2; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_RPA_CCQE_Genie_Weights.push_back(sysweight);
    }
  }

  //3
  if(entry.first == "XSecShape_CCMEC_Genie"){
    for(int i = 0; i < 2; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_XSecShape_CCMEC_Genie_Weights.push_back(sysweight);
    }
  }

  //4
  if(entry.first == "AxFFCCQEshape_Genie"){
    for(int i = 0; i < 2; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_AxFFCCQEshape_Genie_Weights.push_back(sysweight);
    }
  }

  //5
  if(entry.first == "VecFFCCQEshape_Genie"){
    for(int i = 0; i < 2; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_VecFFCCQEshape_Genie_Weights.push_back(sysweight);
    }
  }

  //6
  if(entry.first == "DecayAngMEC_Genie"){
    for(int i = 0; i < 2; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_DecayAngMEC_Genie_Weights.push_back(sysweight);
    }
  }

  //7
  if(entry.first == "Theta_Delta2Npi_Genie"){
    for(int i = 0; i < 2; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_Theta_Delta2Npi_Genie_Weights.push_back(sysweight);
    }
  }

  ///////////////// flux knobs /////////////////////////////

  //1
  if(entry.first == "expskin_FluxUnisim"){
    for(int i = 0; i < 100; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_expskin_FluxUnisim_Weights.push_back(sysweight);
    }
  }//flux1

  //2
  if(entry.first == "horncurrent_FluxUnisim"){
    for(int i = 0; i < 100; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_horncurrent_FluxUnisim_Weights.push_back(sysweight);
    }
  }//flux2

  //3
  if(entry.first == "kminus_PrimaryHadronNormalization"){
    for(int i = 0; i < 100; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_kminus_PrimaryHadronNormalization_Weights.push_back(sysweight);
    }
  }//flux3

  //4
  if(entry.first == "kplus_PrimaryHadronFeynmanScaling"){
    for(int i = 0; i < 100; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_kplus_PrimaryHadronFeynmanScaling_Weights.push_back(sysweight);
    }
  }//flux4

  //5
  if(entry.first == "kzero_PrimaryHadronSanfordWang"){
    for(int i = 0; i < 100; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_kzero_PrimaryHadronSanfordWang_Weights.push_back(sysweight);
    }
  }//flux5

  //6
  if(entry.first == "nucleoninexsec_FluxUnisim"){
    for(int i = 0; i < 100; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_nucleoninexsec_FluxUnisim_Weights.push_back(sysweight);
    }
  }//flux6

  //7
  if(entry.first == "nucleonqexsec_FluxUnisim"){
    for(int i = 0; i < 100; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_nucleonqexsec_FluxUnisim_Weights.push_back(sysweight);
    }
  }//flux7

  //8
  if(entry.first == "nucleontotxsec_FluxUnisim"){
    for(int i = 0; i < 100; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_nucleontotxsec_FluxUnisim_Weights.push_back(sysweight);
    }
  }//flux8

  //9
  if(entry.first == "piminus_PrimaryHadronSWCentralSplineVariation"){
    for(int i = 0; i < 100; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_piminus_PrimaryHadronSWCentralSplineVariation_Weights.push_back(sysweight);
    }
  }//flux9

  //10
  if(entry.first == "pioninexsec_FluxUnisim"){
    for(int i = 0; i < 100; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_pioninexsec_FluxUnisim_Weights.push_back(sysweight);
    }
  }//flux10

  //11
  if(entry.first == "pionqexsec_FluxUnisim"){
    for(int i = 0; i < 100; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_pionqexsec_FluxUnisim_Weights.push_back(sysweight);
    }
  }//flux11

  //12
  if(entry.first == "piontotxsec_FluxUnisim"){
    for(int i = 0; i < 100; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_piontotxsec_FluxUnisim_Weights.push_back(sysweight);
    }
  }//flux12

  //13
  if(entry.first == "piplus_PrimaryHadronSWCentralSplineVariation"){
    for(int i = 0; i < 100; i++){
      double sysweight = 1.;
      if(entry.second.at(i) >= 0 && entry.second.at(i) <= 100000.){sysweight = entry.second.at(i)*splineweight;}
      _f_piplus_PrimaryHadronSWCentralSplineVariation_Weights.push_back(sysweight);
    }
  }//flux13
  
      }//systematics event weights loop

    }//event has weights
  }//fIsMC

   //DETECTOR
  auto const TPC = (*geom).begin_TPC();
  auto ID = TPC.ID();
  m_Cryostat = ID.Cryostat;
  m_TPC = ID.TPC;
  _time2cm = theDetector->SamplingRate() / 1000.0 * theDetector->DriftVelocity( theDetector->Efield(), theDetector->Temperature() );//found in ProtoShowerPandora_tool.cc

  std::cout<<"Test event loop with m_TPC = "<<m_TPC<<" m_Cryostat = "<<m_Cryostat<<" _time2cm = "<<_time2cm<<std::endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////// NEUTRINO PFPARTICLES /////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //NEUTRINO-INDUCED PFPARTICLES
  std::cout<<"FOUND "<<PFParticleVector.size()<<" PFPARTICLES! "<<std::endl;
  std::vector<art::Ptr<recob::PFParticle> > nuPFParticles;
  std::vector<art::Ptr<recob::PFParticle> > nuVertexPFP;

  //Get neutrino slice, and children of primary neutrino pfp
  for(auto const& pfp : PFParticleVector){
    
    //GET NEUTRINO PFP
    if(pfp->IsPrimary() && std::abs(pfp->PdgCode()) == 14){
      nuVertexPFP.push_back(pfp);
    }//neutrino slice
    
    bool isNuChild = false;
    if(pfp->IsPrimary())
      continue;
    
    auto parent = pfp->Parent();
    
    for(auto const& pfp_parent : PFParticleVector){
      if(pfp_parent->Self() != parent)
    continue;
      if(!pfp_parent->IsPrimary())
    continue;
      if(std::abs(pfp_parent->PdgCode()) != 12 && std::abs(pfp_parent->PdgCode()) != 14 && std::abs(pfp_parent->PdgCode()) != 16)
    continue;
      isNuChild = true;
    }
    
    if(isNuChild == true){
      nuPFParticles.push_back(pfp);
    }
    
  }//pfparticle loop

  //Pass only neutrino events
  if(nuVertexPFP.size() == 1){

    /////////////////////////////////////////////////
    ///////// CANDIDATE NEUTRINO VERTEX /////////////
    /////////////////////////////////////////////////
    
    art::FindManyP<recob::Vertex> pfPartToVertexAssoc(PFParticleHandle, e, "pandora");

    std::cout<<"HEY! EVENT HAS "<<nuVertexPFP.size()<<" CANDIDATE NEUTRINO PFP!"<<std::endl;
    
    for(auto const& pfp : nuVertexPFP){
      const auto associatedVertex = pfPartToVertexAssoc.at(pfp.key());
      
      std::cout<<"HEY! CANDIDATE NEUTRINO PFP HAS "<<associatedVertex.size()<<" ASSOCIATED VERTICES!"<<std::endl;

      const std::vector<art::Ptr<anab::T0>> T0_flashchi_v = nuFlashScoreAssoc.at(pfp.key());
    
      if (T0_flashchi_v.size() == 1){
         _fFlashChi2 = T0_flashchi_v.at(0)->TriggerConfidence();

      }
      
      if (associatedVertex.size() > 0){
  for(size_t i = 0; i < associatedVertex.size(); i++){
    auto const vtx = associatedVertex[i];
    nuVertex.push_back(vtx);
  }
      }//associated vertices
      
    }//nupfp
    
    std::cout<<"HEY! FINALLY FOUND "<<nuVertex.size()<<" CANDIDATE VERTICES!"<<std::endl;

    /////////////////////////////////////////////////
    ///////// CANDIDATE MUON TRACK ///////////////////
    /////////////////////////////////////////////////

    //Find tracks & showers associated with neutrino pfp
    art::FindManyP<larpandoraobj::PFParticleMetadata> pfPartToMetaDataAssoc(PFParticleHandle, e, "pandora");
    art::FindManyP<recob::Track> pfPartToTrackAssoc(PFParticleHandle, e, "pandora");
    art::FindManyP<anab::T0> pfPartToT0Assoc(PFParticleHandle, e, "NuCCproducer");
    
    for(auto const& pfp : nuPFParticles){
      const std::vector<art::Ptr<anab::T0>> T0_muon = pfPartToT0Assoc.at(pfp.key());
      std::cout<<"HEY! THIS DAUGHTER PFPARTICLE HAS "<<T0_muon.size()<<" ASSOCIATED MUON T0 INFO! "<<std::endl;

      if(T0_muon.size() == 1){
  const auto associatedTracks = pfPartToTrackAssoc.at(pfp.key()); //tracks associated with candidate muon pfp

  if(associatedTracks.size() == 0){std::cout<<"HEY! THE CANDIDATE MUON PFP HAS NO ASSOCIATED TRACKS!"<<std::endl;}
  if(associatedTracks.size() > 1){std::cout<<"HEY! THE CANDIDATE MUON PFP HAS MANY ASSOCIATED TRACKS!"<<std::endl;}
  if(associatedTracks.size() == 1){
    for(size_t i = 0; i < associatedTracks.size(); i++){
      auto const tr = associatedTracks[i];
      nuMuon.push_back(tr);
      
    }//associated track
  }//1 associated track
      }//is candidate muon pfp
    }//neutrino daughters

    std::cout<<"HEY! THIS EVENT HAS "<<nuMuon.size()<<" CANDIDATE MUON TRACKS!"<<std::endl;

    /////////// PASS ONLY CC INCLUSIVE FILTERED EVENTS
    if(nuVertex.size() == 1 && nuMuon.size() == 1){

      for(auto const& vtx : nuVertex){
  const recob::Vertex::Point_t &neutrino_vtx = vtx->position(); //nuVertex.at(0)->front()->position();
  VertexX = neutrino_vtx.X();
  VertexY = neutrino_vtx.Y();
  VertexZ = neutrino_vtx.Z();
      }//nuvtx
      
      auto muonTrack = nuMuon.front();
      MuonID = muonTrack->ID();
      MuonLength = muonTrack->Length();
      MuonCostheta = std::cos(muonTrack->Theta());
      MuonPhi = muonTrack->Phi();

      //1 = muon, 2 = proton, 3 = charged pion, 4 = photon, 5 = other em, 6 = overlay, 7 = other
      if(fIsMC){MuonTrueComposition = getTrueTrackComposition(e, fHitTag, HitVector, fmht.at(muonTrack.key()), MCParticleVector, fHitPartAssnTag, 2);}

      std::cout<<MuonID<<std::endl;

      trkf::TrackMomentumCalculator tmc;
      MuonMCSMomentum = (MCSFitResultVector.at(muonTrack.key())->bestMomentum())*1000;
      MuonRangeMomentum = 1000*tmc.GetTrackMomentum(MuonLength,13);
    
      float startdist = std::sqrt(std::pow(muonTrack->Vertex().X() - VertexX,2) + std::pow(muonTrack->Vertex().Y() - VertexY,2) + std::pow(muonTrack->Vertex().Z() - VertexZ,2));
      float enddist = std::sqrt(std::pow(muonTrack->End().X() - VertexX,2) + std::pow(muonTrack->End().Y() - VertexY,2) + std::pow(muonTrack->End().Z() - VertexZ,2));
      if(startdist <= enddist){
  MuonStartX = muonTrack->Vertex().X();
  MuonStartY = muonTrack->Vertex().Y();
  MuonStartZ = muonTrack->Vertex().Z();
  MuonEndX = muonTrack->End().X();
  MuonEndY = muonTrack->End().Y();
  MuonEndZ = muonTrack->End().Z();
      }
      if(startdist > enddist){
  MuonStartX = muonTrack->End().X();
  MuonStartY = muonTrack->End().Y();
  MuonStartZ = muonTrack->End().Z();
  MuonEndX = muonTrack->Vertex().X();
  MuonEndY = muonTrack->Vertex().Y();
  MuonEndZ = muonTrack->Vertex().Z();
      }

      MuonDist3d = std::min(startdist,enddist);
  
      if(inFV(muonTrack->Vertex().X(),muonTrack->Vertex().Y(),muonTrack->Vertex().Z()) == true && inFV(muonTrack->End().X(),muonTrack->End().Y(),muonTrack->End().Z()) == true){
  MuonIsContained = 1;
  MuonMomentum = MuonRangeMomentum;
      }
  
      if(inFV(muonTrack->Vertex().X(),muonTrack->Vertex().Y(),muonTrack->Vertex().Z()) == false || inFV(muonTrack->End().X(),muonTrack->End().Y(),muonTrack->End().Z()) == false){
  MuonIsContained = 0;
  MuonMomentum = MuonMCSMomentum;
      }
    }//ccinc

    _fCandidateVertexX = VertexX;
    _fCandidateVertexY = VertexY;
    _fCandidateVertexZ = VertexZ;
    _fCandidateMuonStartX = MuonStartX;
    _fCandidateMuonStartY = MuonStartY;
    _fCandidateMuonStartZ = MuonStartZ;
    _fCandidateMuonEndX = MuonEndX;
    _fCandidateMuonEndY = MuonEndY;
    _fCandidateMuonEndZ = MuonEndZ;
    _fCandidateMuonMomentum = MuonMomentum;
    _fCandidateMuonMCSMomentum = MuonMCSMomentum;
    _fCandidateMuonRangeMomentum = MuonRangeMomentum;
    _fCandidateMuonLength = MuonLength;
    _fCandidateMuonCostheta = MuonCostheta;
    _fCandidateMuonPhi = MuonPhi;
    _fCandidateMuonTrueComposition = MuonTrueComposition;
    _fCandidateMuonIsContained = MuonIsContained;
    _fCandidateMuonDist3d = MuonDist3d;
    if(fIsMC){
      if(_fNuInFV == 1 && std::abs(_fNuPDG) == 14 && _fNuCCNC == 0 && _fNpi0 == 1 && _fNmuon == 1 && _fNpiplus == 0){
  _fCandidateMuonMomentumResolution = (_fLepP - _fCandidateMuonMomentum)/_fLepP;
  _fCandidateMuonRangeMomentumResolution = (_fLepP - _fCandidateMuonRangeMomentum)/_fLepP;
  _fCandidateMuonMCSMomentumResolution = (_fLepP - _fCandidateMuonMCSMomentum)/_fLepP;
      }
    }//fIsMC
      
      for(auto const& pfp : nuPFParticles){
  const auto associatedTracks = pfPartToTrackAssoc.at(pfp.key());
  const auto associatedShowers = pfPartToShowerAssoc.at(pfp.key());
  const auto associatedMetaData = pfPartToMetaDataAssoc.at(pfp.key());
  const auto associatedTrackFitter = trackFitterAssoc.at(pfp.key());     
  
  float score = -999;
  
  //Associated metadata
  if(associatedMetaData.empty()){
    std::cout<<" PFParticle has no metadata!!!"<<std::endl;
  }
  if(associatedMetaData.size() > 1){
    std::cout<<" PFParticle has >1 metadata!!!"<<std::endl;
  }
  if(!associatedMetaData.empty()){
    if(associatedMetaData.size() == 1){
      const auto metaMap = associatedMetaData.at(0)->GetPropertiesMap();
      const auto scoreMap = metaMap.find("TrackScore");
      score = scoreMap->second;
      
    }
  }
  
  //Associated tracks
  if (associatedTracks.size() > 1){
    std::cout << "PFParticle has >1 track!" << std::endl;
  }
  if (associatedTracks.size() > 0){
    for(size_t i = 0; i < associatedTracks.size(); i++){
      auto const tr = associatedTracks[i];
      nuTracks.push_back(tr);
      nuTracks_ScoreMap[tr] = score;
      std::cout<<"PFP MATCHED TO TRACK!! SCORE = "<<score<<std::endl;
    }
  }
  
  //Associated showers
  if (associatedShowers.size() > 1){
    std::cout << "PFParticle has >1 shower!" << std::endl;
  }
  if (associatedShowers.size() > 0){
    for(size_t i = 0; i < associatedShowers.size(); i++){
      auto const sh = associatedShowers[i];
      nuShowers.push_back(sh);
      nuShowers_ScoreMap[sh] = score;
      std::cout<<"PFP MATCHED TO SHOWER!! SCORE = "<<score<<std::endl;
    }
  }
  
  //Associate track fitter to shower
  if(associatedShowers.size() == 1 && associatedTrackFitter.size() == 1){
    
    std::cout<<"FOUND EXACTLY ONE SHOWER AND ONE TRACK FITTER ASSOCIATED WITH THIS PFPARTICLE!"<<std::endl;
    
    for(size_t i = 0; i < associatedShowers.size(); i++){
      auto const sh = associatedShowers[i];
      
      for(size_t i = 0; i < associatedTrackFitter.size(); i++){
        auto const trk = associatedTrackFitter[i];
        int ID = trk->ID();
        
        nuShowerTrackMap[sh] = ID;
        
      }//track fitter
      
    }//shower
    
  }
  
      }//nupfp

  }//has neutrino
  
  _fNTracks = nuTracks.size();
  _fNShowers = nuShowers.size();
  if(nuVertex.size() == 1){_fHasCandidateNeutrino = 1;} else{_fHasCandidateNeutrino = 0;}
  if(nuMuon.size() == 1){_fHasCandidateMuon = 1;} else{_fHasCandidateMuon = 0;}

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////// CC PI0 SELECTION CUTS /////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////// CHARGED PION VETO /////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////

  art::FindMany<anab::ParticleID> fmpid(TrackHandle, e, "pandoracalipidSCE");
  
  for(auto const& track : nuTracks){

    auto length = track->Length();
    std::cout<<length<<std::endl;

    //Find 3-plane PID
    float chi2 = 99999999999;
    if(fmpid.isValid()){
      std::vector<const anab::ParticleID*> pids = fmpid.at(track.key());
      for (size_t ipid=0; ipid < pids.size(); ipid++){
  std::vector<anab::sParticleIDAlgScores> AlgScoresVec = pids[ipid]->ParticleIDAlgScores();
  for (size_t i_algscore=0; i_algscore < AlgScoresVec.size(); i_algscore++){
    anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
    //Get 3-plane pid
    if (AlgScore.fAlgName == "ThreePlaneProtonPID"){
      if (anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood){
        chi2 = std::log(AlgScore.fValue);   
      }//has chi2
    }//3-plane pid
  }//end loop though AlgScoresVec
      }//end loop over pid[ipid]
    }// fmpid.isValid()

    float startdist = std::sqrt(std::pow(track->Vertex().X() - VertexX,2) + std::pow(track->Vertex().Y() - VertexY,2) + std::pow(track->Vertex().Z() - VertexZ,2));
    float enddist = std::sqrt(std::pow(track->End().X() - VertexX,2) + std::pow(track->End().Y() - VertexY,2) + std::pow(track->End().Z() - VertexZ,2));
    float dist_from_vertex = std::min(startdist,enddist);
  
    //Avoid candidate muon
    if(track->ID() == MuonID){ MuonPID = chi2; continue;}
    //Avoid tracks starting away from vertex, i.e. cosmics or EM stuff
    if(dist_from_vertex > 4.) continue;
    //Avoid protons
    if(chi2 > -0.9 && inFV(track->Vertex().X(),track->Vertex().Y(),track->Vertex().Z()) == true && inFV(track->End().X(),track->End().Y(),track->End().Z()) == true && length < MuonLength) continue;
    //Absolute PID cut
    if(chi2 > -0.9) continue;
    nchargedpicand++;

  }//track loop
 
  _fCandidateMuonPID = MuonPID;

  std::cout<<nchargedpicand<<std::endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////// CHOOSE PHOTON CANDIDATES //////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  vector<int> leadingIdx;
  vector<int> subleadIdx;
  int showerNumber = 0;
  for(auto const& shower : nuShowers) {
    double dist3d = conversionDistance(shower, VertexX, VertexY, VertexZ);
    double ang3d = radialAngle(shower, VertexX, VertexY, VertexZ);
    double energy2 = (showerEnergy(2, fmhs.at(shower.key())));

    //lead candidate selection
    if((dist3d > minConversionDist && dist3d < maxConversionDist) || (energy2 > minYPlaneEnergy && energy2 < maxYPlaneEnergy && dist3d < minConversionDist)){
      LeadCandidatesCL.push_back(shower);
      if(energy2 > minYPlaneEnergy){
  LeadCandidatesEN.push_back(shower);
  if(ang3d > minRadialAngle){
    LeadCandidatesRA.push_back(shower);
    LeadCandidates.push_back(shower);
    leadingIdx.push_back(showerNumber);
  }
      }
    }
  ++showerNumber;  
  }//showers

  //FIND LEAD CANDIDATE
  int lead_id = 0;
  float lead_energy = -999.;
  for(auto const& shower : LeadCandidates){
    int id = shower->ID();
    float energy = showerEnergy(2, fmhs.at(shower.key()));
    if(energy > lead_energy){
      lead_energy = energy;
      lead_id = id;
    }//max energy
  }//lead loop
  
  vector<int>::const_iterator itShower = leadingIdx.begin();
  for(auto const& shower : LeadCandidates){
    if(shower->ID() == lead_id){
      TheLeadVector.push_back(shower);
      _fLeadShowerIDX = *itShower;
    }
    ++itShower;
  }
  
  //FIND SUBLEAD CANDIDATES
  showerNumber = 0;
  for(auto const& shower : nuShowers){
    int id = shower->ID();
    double dist3d = conversionDistance(shower, VertexX, VertexY, VertexZ);
    float energy2 = showerEnergy(2, fmhs.at(shower.key()));
      
    if(id == lead_id){
      ++showerNumber; //still need to increment this
      continue;
    } 

    if(dist3d > 1){
      SubleadCandidatesCL.push_back(shower);
      if(energy2 > 25){
  SubleadCandidatesEN.push_back(shower);
  SubleadCandidates.push_back(shower);
  subleadIdx.push_back(showerNumber);
      }//energy
    }//dist3d
  ++showerNumber;
  }//showers
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////// MASS CUT TO MATCH PHOTONS /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////
  
  if(TheLeadVector.size() == 1){
    auto thelead = TheLeadVector.front();
    double leadenergy = (showerEnergy(2, fmhs.at(thelead.key())))/altbias_2;
    itShower = subleadIdx.begin();
    //Sublead
    for(auto const& sublead : SubleadCandidates){
      double subleadenergy = (showerEnergy(2, fmhs.at(sublead.key())))/altbias_2;
      double angle12 = angleBetweenTwoShowers(thelead,sublead);
      double pi0mass = pi0Mass(leadenergy,subleadenergy,angle12);
      _fTwoPhotonInvariantMass.push_back(pi0mass);

      if(pi0mass > signalMassLow && pi0mass < signalMassHigh){
         npairs++;
         TheSubleadVector.push_back(sublead);
         _fSubLeadShowerIDX = *itShower;
      }
      //high pi0 mass sideband
      if(pi0mass >sidebandMassLow && pi0mass < sidebandMassHigh){
         SubleadVectorHighMass.push_back(sublead);
         ++npairsHiMass;
      }
     ++itShower;
    
    }

    /////////// AVOID OTHER PHOTON PAIRS FROM PI0 /////////////
  
    if(npairs == 1){

      auto sublead1 = TheSubleadVector.front();
      auto subleadid1 = sublead1->ID();

      /////// PATH 2: Check for other actual lead-sublead pairs ///////////

      for(auto const& shower : SubleadCandidates){
  auto showerid = shower->ID();
  double dist3d = conversionDistance(shower, VertexX, VertexY, VertexZ);
  double ang3d = radialAngle(shower, VertexX, VertexY, VertexZ);
  double energy2 = (showerEnergy(2, fmhs.at(shower.key())));
  double showerenergy = (showerEnergy(2, fmhs.at(shower.key())))/altbias_2;
      
  if(showerid == subleadid1) continue;

  for(auto const& shower2 : SubleadCandidates){
    auto showerid2 = shower2->ID();
    double dist3d2 = conversionDistance(shower2, VertexX, VertexY, VertexZ);
    double ang3d2 = radialAngle(shower2, VertexX, VertexY, VertexZ);
    double energy22 = (showerEnergy(2, fmhs.at(shower2.key())));
    double showerenergy2 = (showerEnergy(2, fmhs.at(shower2.key())))/altbias_2;

    if(showerid2 == showerid || showerid2 == subleadid1 || showerenergy2 > showerenergy) continue;
     
    /*
    (((dist3d > 6. && dist3d < 82.) || (energy2 > 50 && energy2 < 300 && dist3d < 6)) && (energy2 > 50) && (ang3d > 0.96)) ||
    (((dist3d2 > 6. && dist3d2 < 82.) || (energy22 > 50 && energy22 < 300 && dist3d2 < 6)) && (energy22 > 50) && (ang3d2 > 0.96))
    */
    
    //Cuts on shower energy, angle and conversion dist. 
    if(passesShowerCuts(dist3d, energy2, ang3d) || passesShowerCuts(dist3d2, energy22, ang3d2)){

      double angle12 = angleBetweenTwoShowers(shower,shower2);
      double pi0mass = pi0Mass(showerenergy,showerenergy2,angle12);
    
      if(pi0mass > signalMassLow && pi0mass < signalMassHigh){
        notherfancypairs++; //this checks subleading - subleading pairs
      }
    }//one is a lead

        
  }//loop2

      }//loop1

    }//1 pair

  }//1 lead vector

  _fNChargedPiCandidates = nchargedpicand;
  _fNLeadCandidates = LeadCandidates.size();
  _fNSubleadCandidates = SubleadCandidates.size();
  _fNPairCandidates = npairs;
  _fNOtherFancyPairs = notherfancypairs;
  _fNHiMassPairs     = npairsHiMass;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if(nuVertexPFP.size() == 1 && nuVertex.size() == 1 && nuMuon.size() == 1){

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////// TRACK DISTRIBUTIONS //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    for(auto const& track : nuTracks){
      auto score = nuTracks_ScoreMap[track];
      float length = track->Length();
      float costheta = std::cos(track->Theta());
      float phi = track->Phi();
      float momentum = (MCSFitResultVector.at(track.key())->bestMomentum())*1000;
      float pid = -99999.9;
      float dist3d = -9999.9;
      int iscontained = -9999;
      int composition = -999;
      float startx = -99999.9; float starty = -99999.9; float startz = -99999.9;
      float endx = -99999.9; float endy = -99999.9; float endz = -99999.9;
      float startdist = std::sqrt(std::pow(track->Vertex().X() - VertexX,2) + std::pow(track->Vertex().Y() - VertexY,2) + std::pow(track->Vertex().Z() - VertexZ,2));
      float enddist = std::sqrt(std::pow(track->End().X() - VertexX,2) + std::pow(track->End().Y() - VertexY,2) + std::pow(track->End().Z() - VertexZ,2));
      if(startdist <= enddist){
         startx = track->Vertex().X();
         starty = track->Vertex().Y();
         startz = track->Vertex().Z();
         endx = track->End().X();
         endy = track->End().Y();
         endz = track->End().Z();
      }
      if(startdist > enddist){
         startx = track->End().X();
         starty = track->End().Y();
         startz = track->End().Z();
         endx = track->Vertex().X();
         endy = track->Vertex().Y();
         endz = track->Vertex().Z();
      }
      dist3d = std::min(startdist,enddist);

      if(inFV(track->Vertex().X(),track->Vertex().Y(),track->Vertex().Z()) == true && inFV(track->End().X(),track->End().Y(),track->End().Z()) == true){
  iscontained = 1;
      }

      if(inFV(track->Vertex().X(),track->Vertex().Y(),track->Vertex().Z()) == false || inFV(track->End().X(),track->End().Y(),track->End().Z()) == false){
  iscontained = 0;
      }

      //1 = muon, 2 = proton, 3 = charged pion, 4 = photon, 5 = other em, 6 = overlay, 7 = other
      if(fIsMC){composition = getTrueTrackComposition(e, fHitTag, HitVector, fmht.at(track.key()), MCParticleVector, fHitPartAssnTag, 2);}

      //Find 3-plane PID
      if(fmpid.isValid()){
  std::vector<const anab::ParticleID*> pids = fmpid.at(track.key());
  for (size_t ipid=0; ipid < pids.size(); ipid++){
    std::vector<anab::sParticleIDAlgScores> AlgScoresVec = pids[ipid]->ParticleIDAlgScores();
    for (size_t i_algscore=0; i_algscore < AlgScoresVec.size(); i_algscore++){
      anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
      //Get 3-plane pid
      if (AlgScore.fAlgName == "ThreePlaneProtonPID"){
        if (anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood){
    pid = std::log(AlgScore.fValue);    
        }//has chi2
      }//3-plane pid
    }//end loop though AlgScoresVec
  }//end loop over pid[ipid]
      }// fmpid.isValid()

      _fTrackShowerScore.push_back(score);
      _fTrackTrueComposition.push_back(composition);
      _fTrackIsContained.push_back(iscontained);
      _fTrackDist3d.push_back(dist3d);
      _fTrackPID.push_back(pid);
      _fTrackStartX.push_back(startx);
      _fTrackStartY.push_back(starty);
      _fTrackStartZ.push_back(startz);
      _fTrackEndX.push_back(endx);
      _fTrackEndY.push_back(endy);
      _fTrackEndZ.push_back(endz);
      _fTrackMomentum.push_back(momentum);
      _fTrackLength.push_back(length);
      _fTrackCostheta.push_back(costheta);
      _fTrackPhi.push_back(phi);

    }//track loop

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////// SHOWER DISTRIBUTIONS /////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
    for(auto const& shower : nuShowers){
      int trackID = nuShowerTrackMap[shower];
      double score = nuShowers_ScoreMap[shower];
      double length = shower->Length();
      double openangle = radToDeg(shower->OpenAngle());
      double startx = shower->ShowerStart().X();
      double starty = shower->ShowerStart().Y();
      double startz = shower->ShowerStart().Z();
      double dirx = shower->Direction().X(); //std::acos(shower->Direction().X());
      double diry = shower->Direction().Y();
      double dirz = shower->Direction().Z();
      double dist3d = conversionDistance(shower, VertexX, VertexY, VertexZ);
      double ang3d = radialAngle(shower, VertexX, VertexY, VertexZ);
      double energy0 = (showerEnergy(0, fmhs.at(shower.key())));
      double energy1 = (showerEnergy(1, fmhs.at(shower.key())));
      double energy2 = (showerEnergy(2, fmhs.at(shower.key())));
      double dedx0 = -999; double dedx1 = -999; double dedx2 = -999; double dedxamalg = -999; 
      getdEdx(e, fCaloTag, fShowerTag, fPFParticleTag, fTrackTag, fTrackFitterTag, fClusterTag, shower, nuPFParticles, TrackFitterVector, trackID, trunklength, trunkwidth, dedx0, dedx1, dedx2,dedxamalg);

      int bestplane = bestPlane(fmhs.at(shower.key()));
      float bestdedx = -999.;
      if(bestplane == 0){
  bestdedx = dedx0;
      }
      if(bestplane == 1){
  bestdedx = dedx1;
      }
      if(bestplane == 2){
  bestdedx = dedx2;
      }
  
      if(fIsMC){
  float purity0 = getPurity(e, fHitTag, HitVector, fmhs.at(shower.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 0);
  float complete0 = getComplete(e, fHitTag, HitVector, fmhs.at(shower.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 0);
  float energyres0 = getEnergyRes(e, fHitTag, HitVector, fmhs.at(shower.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, _fGamma1_E, _fGamma2_E, energy2, fHitPartAssnTag, 0);

  float purity1 = getPurity(e, fHitTag, HitVector, fmhs.at(shower.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 1);
  float complete1 = getComplete(e, fHitTag, HitVector, fmhs.at(shower.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 1);
  float energyres1 = getEnergyRes(e, fHitTag, HitVector, fmhs.at(shower.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, _fGamma1_E, _fGamma2_E, energy2, fHitPartAssnTag, 1);

  float purity2 = getPurity(e, fHitTag, HitVector, fmhs.at(shower.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 2);
  float complete2 = getComplete(e, fHitTag, HitVector, fmhs.at(shower.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 2);
  float energyres2 = getEnergyRes(e, fHitTag, HitVector, fmhs.at(shower.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, _fGamma1_E, _fGamma2_E, energy2, fHitPartAssnTag, 2);
  
  //float trueenergy = getTrueEnergy(e, fHitTag, HitVector, fmhs.at(shower.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, _fGamma1_E, _fGamma2_E, fHitPartAssnTag, 2);
  float comp = getTrueComposition(e, fHitTag, HitVector, fmhs.at(shower.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, Gamma1_descendents, Gamma2_descendents, fHitPartAssnTag, 2);
  
  _fShowerTrueComposition.push_back(comp);
  _fShowerEnergy0Resolution.push_back(energyres0);
  _fShowerEnergy1Resolution.push_back(energyres1);
  _fShowerEnergy2Resolution.push_back(energyres2);
  _fShowerPurity0.push_back(purity0);
  _fShowerPurity1.push_back(purity1);
  _fShowerPurity2.push_back(purity2);
  _fShowerCompleteness0.push_back(complete0);
  _fShowerCompleteness1.push_back(complete1);
  _fShowerCompleteness2.push_back(complete2);
   
      }//fIsMC

      _fTrackShowerScore.push_back(score);
      _fShowerLength.push_back(length);
      _fShowerOpenAngle.push_back(openangle);
      _fShowerCostheta.push_back(dirz);
      _fShowerPhi.push_back(std::acos(dirx));
      _fShowerDist3d.push_back(dist3d);
      _fShowerAng3d.push_back(ang3d);
      _fShowerStartX.push_back(startx);
      _fShowerStartY.push_back(starty);
      _fShowerStartZ.push_back(startz);
      _fShowerEnergy0.push_back(energy0);
      _fShowerEnergy1.push_back(energy1);
      _fShowerEnergy2.push_back(energy2);
      _fShowerdEdx0.push_back(dedx0);
      _fShowerdEdx1.push_back(dedx1);
      _fShowerdEdx2.push_back(dedx2);
      _fShowerdEdx3.push_back(bestdedx);
      _fShowerdDirX.push_back(dirx);
      _fShowerdDirY.push_back(diry);
      _fShowerdDirZ.push_back(dirz);


    }//shower loop

  }//ccinc

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// SELECTED PI0 & GAMMA DISTRIBUTIONS ////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //IMPORTANT
  //All events should pass these cuts, even in the sidebands
  if(nuVertexPFP.size() == 1 && nuVertex.size() == 1 && nuMuon.size() == 1  && TheLeadVector.size() == 1) {
  
    //da selection nchargedpicand == 0 && nuShowers.size() >= 2 && TheSubleadVector.size() == 1 && npairs == 1 && notherfancypairs == 0)
    auto leadcandidate = TheLeadVector.front();
    
      
    int leadtrackID = nuShowerTrackMap[leadcandidate];
    double leadscore = nuShowers_ScoreMap[leadcandidate];
    double leadlength = leadcandidate->Length();
    double leadopenangle = radToDeg(leadcandidate->OpenAngle());
    double leadstartx = leadcandidate->ShowerStart().X();
    double leadstarty = leadcandidate->ShowerStart().Y();
    double leadstartz = leadcandidate->ShowerStart().Z();
    double leaddirx = leadcandidate->Direction().X();
    double leaddiry = leadcandidate->Direction().Y();
    double leaddirz = leadcandidate->Direction().Z();
    double leaddist3d = conversionDistance(leadcandidate, VertexX, VertexY, VertexZ);
    double leadang3d = radialAngle(leadcandidate, VertexX, VertexY, VertexZ);
    double leadenergy0 = (showerEnergy(0, fmhs.at(leadcandidate.key())));
    double leadenergy1 = (showerEnergy(1, fmhs.at(leadcandidate.key())));
    double leadenergy2 = (showerEnergy(2, fmhs.at(leadcandidate.key())))/altbias_2;
    double leaddedx0 = -999; double leaddedx1 = -999; double leaddedx2 = -999; double leaddedxamalg = -999; 
    getdEdx(e, fCaloTag, fShowerTag, fPFParticleTag, fTrackTag, fTrackFitterTag, fClusterTag, leadcandidate, nuPFParticles, TrackFitterVector, leadtrackID, trunklength, trunkwidth, leaddedx0, leaddedx1, leaddedx2,leaddedxamalg);

    int leadbestplane = bestPlane(fmhs.at(leadcandidate.key()));
    float leadbestdedx = -999.;
    if(leadbestplane == 0){
      leadbestdedx = leaddedx0;
    }
    if(leadbestplane == 1){
      leadbestdedx = leaddedx1;
    }
    if(leadbestplane == 2){
      leadbestdedx = leaddedx2;
    }

    double  leadpurity0 = -999.9;
    double  leadcomplete0 = -999.9;
    double  leadenergyres0 = -999.9;
    double  leadpurity1 = -999.9;
    double  leadcomplete1 = -999.9;
    double  leadenergyres1 = -999.9;
    double  leadpurity2 = -999.9;
    double  leadcomplete2 = -999.9;
    double  leadenergyres2 = -999.9;
    double  leadcomp = -999.9;

    if(fIsMC){
          leadpurity0 = getPurity(e, fHitTag, HitVector, fmhs.at(leadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 0);
          leadcomplete0 = getComplete(e, fHitTag, HitVector, fmhs.at(leadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 0);
          leadenergyres0 = getEnergyRes(e, fHitTag, HitVector, fmhs.at(leadcandidate.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, _fGamma1_E, _fGamma2_E, leadenergy2*altbias_2, fHitPartAssnTag, 0);
          leadpurity1 = getPurity(e, fHitTag, HitVector, fmhs.at(leadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 1);
          leadcomplete1 = getComplete(e, fHitTag, HitVector, fmhs.at(leadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 1);
          leadenergyres1 = getEnergyRes(e, fHitTag, HitVector, fmhs.at(leadcandidate.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, _fGamma1_E, _fGamma2_E, leadenergy2*altbias_2, fHitPartAssnTag, 1);
          leadpurity2 = getPurity(e, fHitTag, HitVector, fmhs.at(leadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 2);
          leadcomplete2 = getComplete(e, fHitTag, HitVector, fmhs.at(leadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 2);
          leadenergyres2 = getEnergyRes(e, fHitTag, HitVector, fmhs.at(leadcandidate.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, _fGamma1_E, _fGamma2_E, leadenergy2*altbias_2, fHitPartAssnTag, 2);
          leadcomp = getTrueComposition(e, fHitTag, HitVector, fmhs.at(leadcandidate.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, Gamma1_descendents, Gamma2_descendents, fHitPartAssnTag, 2);
    }      

    //Sidebands and signal events require at least one subleading candidate
    if(TheSubleadVector.size() > 0){
     //Intialize all this stuff
     int subleadtrackID = -999;
     double subleadscore = -999.9;
     double subleadlength = -999.9;
     double subleadopenangle = -999.9;
     double subleadstartx = -999.9;
     double subleadstarty = -999.9;
     double subleadstartz = -999.9;
     double subleaddirx =  -999.9;
     double subleaddiry = -999.9;
     double subleaddirz = -999.9;
     double subleaddist3d = -999.9;
     double subleadang3d = -999.9;
     double subleadenergy0 = -999.9;
     double subleadenergy1 = -999.9;
     double subleadenergy2 = -999.9;
     double subleaddedx0 = -999.9; 
     double subleaddedx1 = -999.9; 
     double subleaddedx2 = -999.9; 
     double subleaddedxamalg = -999.9;
     float subleadbestdedx = -999.9;
     int subleadbestplane = -999;

     //First we look at distributions that require exactly 1 subleading shower:

     if(TheSubleadVector.size() == 1){
       auto subleadcandidate = TheSubleadVector.front();
       subleadtrackID = nuShowerTrackMap[subleadcandidate];
       subleadscore = nuShowers_ScoreMap[subleadcandidate];
       subleadlength = subleadcandidate->Length();
       subleadopenangle = radToDeg(subleadcandidate->OpenAngle());
       subleadstartx = subleadcandidate->ShowerStart().X();
       subleadstarty = subleadcandidate->ShowerStart().Y();
       subleadstartz = subleadcandidate->ShowerStart().Z();
       subleaddirx = subleadcandidate->Direction().X();
       subleaddiry = subleadcandidate->Direction().Y();
       subleaddirz = subleadcandidate->Direction().Z();
       subleaddist3d = conversionDistance(subleadcandidate, VertexX, VertexY, VertexZ);
       subleadang3d = radialAngle(subleadcandidate, VertexX, VertexY, VertexZ);
       subleadenergy0 = (showerEnergy(0, fmhs.at(subleadcandidate.key())));
       subleadenergy1 = (showerEnergy(1, fmhs.at(subleadcandidate.key())));
       subleadenergy2 = (showerEnergy(2, fmhs.at(subleadcandidate.key())))/altbias_2;
       getdEdx(e, fCaloTag, fShowerTag, fPFParticleTag, fTrackTag, fTrackFitterTag, fClusterTag, subleadcandidate, nuPFParticles, TrackFitterVector, subleadtrackID, trunklength2, trunkwidth, subleaddedx0, subleaddedx1, subleaddedx2,subleaddedxamalg);
      
       subleadbestplane = bestPlane(fmhs.at(subleadcandidate.key()));
       
       if(subleadbestplane == 0){
         subleadbestdedx = subleaddedx0;
       }
       
       if(subleadbestplane == 1){
         subleadbestdedx = subleaddedx1;
       }
       
       if(subleadbestplane == 2){
         subleadbestdedx = subleaddedx2;
       }

       /*float leadpurity0    = -999.9;
       float leadcomplete0  = -999.9;
       float leadenergyres0 = -999.9;
       float leadpurity1    = -999.9;
       float leadcomplete1  = -999.9;
       float leadenergyres1 = -999.9;
       float leadpurity2    = -999.9;
       float leadcomplete2  = -999.9;
       float leadenergyres2 = -999.9;
       float leadcomp       = -999.9;
       */
       float subleadpurity0    = -999.9;
       float subleadcomplete0  = -999.9;
       float subleadenergyres0 = -999.9;
       float subleadpurity1    = -999.9;
       float subleadcomplete1  = -999.9;
       float subleadenergyres1 = -999.9;
       float subleadpurity2    = -999.9;
       float subleadcomplete2  = -999.9;
       float subleadenergyres2 = -999.9;
       float subleadcomp       = -999.9;

       double angle12 = angleBetweenTwoShowers(leadcandidate,subleadcandidate);
       double acosangle12 = radToDeg(std::acos(angle12));
       double pi0mass = -99999.;
       double pi0momentum = -99999.;
       double pi0costheta = -99999.;
       double pi0energy = -99999.;
      
       pi0mass = pi0Mass(leadenergy2,subleadenergy2,angle12);
       pi0momentum = Altpi0Momentum(leadenergy2,subleadenergy2,angle12);
       pi0energy = std::sqrt(std::pow(pi0mass,2) + std::pow(pi0momentum,2));
       pi0costheta = pi0Angle(pi0momentum, leadenergy2, subleadenergy2, leadcandidate, subleadcandidate, angle12);

       if(fIsMC){
          subleadpurity0 = getPurity(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 0);
          subleadcomplete0 = getComplete(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 0);
          subleadenergyres0 = getEnergyRes(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, _fGamma1_E, _fGamma2_E, subleadenergy2*altbias_2, fHitPartAssnTag, 0);
          subleadpurity1 = getPurity(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 1);
          subleadcomplete1 = getComplete(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 1);
          subleadenergyres1 = getEnergyRes(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, _fGamma1_E, _fGamma2_E, subleadenergy2*altbias_2, fHitPartAssnTag, 1);
          subleadpurity2 = getPurity(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 2);
          subleadcomplete2 = getComplete(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 2);
          subleadenergyres2 = getEnergyRes(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, _fGamma1_E, _fGamma2_E, subleadenergy2*altbias_2, fHitPartAssnTag, 2);
          subleadcomp = getTrueComposition(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, Gamma1_descendents, Gamma2_descendents, fHitPartAssnTag, 2);
        }//fIsMC

       //Signal Events
       if(nchargedpicand == 0 && nuShowers.size() >= 2 && npairs == 1 && notherfancypairs == 0){
        _fCandidateLeadScore= leadscore;
        _fCandidateLeadLength   = leadlength;
        _fCandidateLeadOpenAngle= leadopenangle;
        _fCandidateLeadCostheta = leaddirz;
        _fCandidateLeadPhi      = std::atan2(leaddiry, leaddirx);
        _fCandidateLeadDirZ     = leaddirz;
        _fCandidateLeadDirX     = leaddirx;
        _fCandidateLeadDirY     = leaddiry;


        _fCandidateLeadDist3d   = leaddist3d;
        _fCandidateLeadAng3d    = leadang3d;
        _fCandidateLeadStartX   = leadstartx;
        _fCandidateLeadStartY   = leadstarty;
        _fCandidateLeadStartZ   = leadstartz;
        _fCandidateLeadEnergy0  = leadenergy0;
        _fCandidateLeadEnergy1  = leadenergy1;
        _fCandidateLeadEnergy2  = leadenergy2;
        _fCandidateLeaddEdx0    = leaddedx0;
        _fCandidateLeaddEdx1    = leaddedx1;
        _fCandidateLeaddEdx2    = leaddedx2;
        _fCandidateLeaddEdx3    = leadbestdedx;

        _fCandidateSubleadScore = subleadscore;
        _fCandidateSubleadLength = subleadlength;
        _fCandidateSubleadOpenAngle = subleadopenangle;
        _fCandidateSubleadCostheta = subleaddirz;
        _fCandidateSubleadPhi      = std::atan2(subleaddiry, subleaddirx);
        _fCandidateSubleadDirZ     = subleaddirz;
        _fCandidateSubleadDirX     = subleaddirx;
        _fCandidateSubleadDirY     = subleaddiry;
        _fCandidateSubleadDist3d = subleaddist3d;
        _fCandidateSubleadAng3d = subleadang3d;
        _fCandidateSubleadStartX = subleadstartx;
        _fCandidateSubleadStartY = subleadstarty;
        _fCandidateSubleadStartZ = subleadstartz;
        _fCandidateSubleadEnergy0 = subleadenergy0;
        _fCandidateSubleadEnergy1 = subleadenergy1;
        _fCandidateSubleadEnergy2 = subleadenergy2;
        _fCandidateSubleaddEdx0 = subleaddedx0;
        _fCandidateSubleaddEdx1 = subleaddedx1;
        _fCandidateSubleaddEdx2 = subleaddedx2;
        _fCandidateSubleaddEdx3 = subleadbestdedx;

        //Candidate pi0
        _fCandidatePi0Momentum = pi0momentum;
        _fCandidatePi0Energy   = pi0energy;
        _fCandidatePi0Costheta = pi0costheta;
        _fCandidatePi0Angle12  = acosangle12;
        _fCandidatePi0Mass     = pi0mass;

        _fCandidateLeadTrueComposition = leadcomp;
        _fCandidateLeadEnergy0Resolution = leadenergyres0;
        _fCandidateLeadEnergy1Resolution = leadenergyres1;
        _fCandidateLeadEnergy2Resolution = leadenergyres2;
        _fCandidateLeadPurity0 = leadpurity0;
        _fCandidateLeadPurity1 = leadpurity1;
        _fCandidateLeadPurity2 = leadpurity2;
        _fCandidateLeadCompleteness0 = leadcomplete0;
        _fCandidateLeadCompleteness1 = leadcomplete1;
        _fCandidateLeadCompleteness2 = leadcomplete2;
        _fCandidateSubleadTrueComposition = subleadcomp;
        _fCandidateSubleadEnergy0Resolution = subleadenergyres0;
        _fCandidateSubleadEnergy1Resolution = subleadenergyres1;
        _fCandidateSubleadEnergy2Resolution = subleadenergyres2;
        _fCandidateSubleadPurity0 = subleadpurity0;
        _fCandidateSubleadPurity1 = subleadpurity1;
        _fCandidateSubleadPurity2 = subleadpurity2;
        _fCandidateSubleadCompleteness0 = subleadcomplete0;
        _fCandidateSubleadCompleteness1 = subleadcomplete1;
        _fCandidateSubleadCompleteness2 = subleadcomplete2;

        if(_fNuInFV == 1 && std::abs(_fNuPDG) == 14 && _fNuCCNC == 0 && _fNpi0 == 1 && _fNmuon == 1 && _fNpiplus == 0 && fIsMC){ //signal: CC 1pi0 0 pi+
               _fCandidatePi0MomentumResolution = (_fPi0P - pi0momentum)/_fPi0P;
           }//signal



       }//end of is reco signal

       //Pi+ Sideband
       if(nchargedpicand > 0 && nuShowers.size() >= 2 && npairs == 1 && notherfancypairs == 0){
          _fTwoMIPLeadScore= leadscore;
          _fTwoMIPLeadLength   = leadlength;
          _fTwoMIPLeadOpenAngle= leadopenangle;
          _fTwoMIPLeadCostheta = leaddirz;
          _fTwoMIPLeadPhi      = std::atan2(leaddiry, leaddirx);
          _fTwoMIPLeadDirZ     = leaddirz;
          _fTwoMIPLeadDirX     = leaddirx;
          _fTwoMIPLeadDirY     = leaddiry;
          _fTwoMIPLeadDist3d   = leaddist3d;
          _fTwoMIPLeadAng3d    = leadang3d;
          _fTwoMIPLeadStartX   = leadstartx;
          _fTwoMIPLeadStartY   = leadstarty;
          _fTwoMIPLeadStartZ   = leadstartz;
          _fTwoMIPLeadEnergy0  = leadenergy0;
          _fTwoMIPLeadEnergy1  = leadenergy1;
          _fTwoMIPLeadEnergy2  = leadenergy2;
          _fTwoMIPLeaddEdx0    = leaddedx0;
          _fTwoMIPLeaddEdx1    = leaddedx1;
          _fTwoMIPLeaddEdx2    = leaddedx2;
          _fTwoMIPLeaddEdx3    = leadbestdedx;

          _fTwoMIPSubleadScore = subleadscore;
          _fTwoMIPSubleadLength = subleadlength;
          _fTwoMIPSubleadOpenAngle = subleadopenangle;
          _fTwoMIPSubleadCostheta = subleaddirz;
          _fTwoMIPSubleadPhi      = std::atan2(subleaddiry, subleaddirx);
          _fTwoMIPSubleadDirZ     = subleaddirz;
          _fTwoMIPSubleadDirX     = subleaddirx;
          _fTwoMIPSubleadDirY     = subleaddiry;
          _fTwoMIPSubleadDist3d = subleaddist3d;
          _fTwoMIPSubleadAng3d = subleadang3d;
          _fTwoMIPSubleadStartX = subleadstartx;
          _fTwoMIPSubleadStartY = subleadstarty;
          _fTwoMIPSubleadStartZ = subleadstartz;
          _fTwoMIPSubleadEnergy0 = subleadenergy0;
          _fTwoMIPSubleadEnergy1 = subleadenergy1;
          _fTwoMIPSubleadEnergy2 = subleadenergy2;
          _fTwoMIPSubleaddEdx0 = subleaddedx0;
          _fTwoMIPSubleaddEdx1 = subleaddedx1;
          _fTwoMIPSubleaddEdx2 = subleaddedx2;
          _fTwoMIPSubleaddEdx3 = subleadbestdedx;

          //Candidate pi0
          _fTwoMIPPi0Momentum = pi0momentum;
          _fTwoMIPPi0Energy = pi0energy;
          _fTwoMIPPi0Costheta = pi0costheta;
          _fTwoMIPPi0Angle12 = acosangle12;
          _fTwoMIPPi0Mass     = pi0mass;

          _fTwoMIPLeadTrueComposition = leadcomp;
          _fTwoMIPLeadEnergy0Resolution = leadenergyres0;
          _fTwoMIPLeadEnergy1Resolution = leadenergyres1;
          _fTwoMIPLeadEnergy2Resolution = leadenergyres2;
          _fTwoMIPLeadPurity0 = leadpurity0;
          _fTwoMIPLeadPurity1 = leadpurity1;
          _fTwoMIPLeadPurity2 = leadpurity2;
          _fTwoMIPLeadCompleteness0 = leadcomplete0;
          _fTwoMIPLeadCompleteness1 = leadcomplete1;
          _fTwoMIPLeadCompleteness2 = leadcomplete2;

          _fTwoMIPSubleadTrueComposition = subleadcomp;
          _fTwoMIPSubleadEnergy0Resolution = subleadenergyres0;
          _fTwoMIPSubleadEnergy1Resolution = subleadenergyres1;
          _fTwoMIPSubleadEnergy2Resolution = subleadenergyres2;
          _fTwoMIPSubleadPurity0 = subleadpurity0;
          _fTwoMIPSubleadPurity1 = subleadpurity1;
          _fTwoMIPSubleadPurity2 = subleadpurity2;
          _fTwoMIPSubleadCompleteness0 = subleadcomplete0;
          _fTwoMIPSubleadCompleteness1 = subleadcomplete1;
          _fTwoMIPSubleadCompleteness2 = subleadcomplete2;          
       }

    }//subleading vector is size 1, one subleading candidate

    //More than one pair sideband
    if(TheSubleadVector.size() > 1 && nchargedpicand == 0 && nuShowers.size() >= 2 && npairs > 1  && notherfancypairs == 0){
      _fMultiPairLeadScore= leadscore;
      _fMultiPairLeadLength   = leadlength;
      _fMultiPairLeadOpenAngle= leadopenangle;
      _fMultiPairLeadCostheta = leaddirz;
      _fMultiPairLeadPhi      = std::atan2(leaddiry, leaddirx);
      _fMultiPairLeadDirZ     = leaddirz;
      _fMultiPairLeadDirX     = leaddirx;
      _fMultiPairLeadDirY     = leaddiry;
      _fMultiPairLeadDist3d   = leaddist3d;
      _fMultiPairLeadAng3d    = leadang3d;
      _fMultiPairLeadStartX   = leadstartx;
      _fMultiPairLeadStartY   = leadstarty;
      _fMultiPairLeadStartZ   = leadstartz;
      _fMultiPairLeadEnergy0  = leadenergy0;
      _fMultiPairLeadEnergy1  = leadenergy1;
      _fMultiPairLeadEnergy2  = leadenergy2;
      _fMultiPairLeaddEdx0    = leaddedx0;
      _fMultiPairLeaddEdx1    = leaddedx1;
      _fMultiPairLeaddEdx2    = leaddedx2;
      _fMultiPairLeaddEdx3    = leadbestdedx;      
      
      for(auto const& subleadcandidate : TheSubleadVector){
         int subleadtrackID = nuShowerTrackMap[subleadcandidate];
         double subleaddedx0 = -999.9; 
         double subleaddedx1 = -999.9; 
         double subleaddedx2 = -999.9; 
         double subleaddedxamalg = -999.9;
         float subleadbestdedx = -999.9;
         int subleadbestplane = -999;
         double angle12 = angleBetweenTwoShowers(leadcandidate,subleadcandidate);

         _fMultiPairSubleadScore.push_back(nuShowers_ScoreMap[subleadcandidate]);
         _fMultiPairSubleadLength.push_back(subleadcandidate->Length() );
         _fMultiPairSubleadOpenAngle.push_back(radToDeg(subleadcandidate->OpenAngle()) );
         _fMultiPairSubleadStartX.push_back(subleadcandidate->ShowerStart().X() );
         _fMultiPairSubleadStartY.push_back(subleadcandidate->ShowerStart().Y() );
         _fMultiPairSubleadStartZ.push_back(subleadcandidate->ShowerStart().Z() );
         _fMultiPairSubleadPhi.push_back(std::atan2(subleadcandidate->Direction().X(), subleadcandidate->Direction().Y()) );
         _fMultiPairSubleadDirX.push_back(subleadcandidate->Direction().X() );
         _fMultiPairSubleadDirY.push_back(subleadcandidate->Direction().Y() );
         _fMultiPairSubleadDirZ.push_back(subleadcandidate->Direction().Z() );
         _fMultiPairSubleadDist3d.push_back(conversionDistance(subleadcandidate, VertexX, VertexY, VertexZ) );
         _fMultiPairSubleadAng3d.push_back(radialAngle(subleadcandidate, VertexX, VertexY, VertexZ) );
         _fMultiPairSubleadEnergy0.push_back(showerEnergy(0, fmhs.at(subleadcandidate.key())));
         _fMultiPairSubleadEnergy1.push_back(showerEnergy(1, fmhs.at(subleadcandidate.key())));
         double subleadenergy2 = showerEnergy(2, fmhs.at(subleadcandidate.key()))/altbias_2;
         _fMultiPairSubleadEnergy2.push_back(subleadenergy2);
         getdEdx(e, fCaloTag, fShowerTag, fPFParticleTag, fTrackTag, fTrackFitterTag, fClusterTag, subleadcandidate, nuPFParticles, TrackFitterVector, subleadtrackID, trunklength2, trunkwidth, subleaddedx0, subleaddedx1, subleaddedx2,subleaddedxamalg);
        
         subleadbestplane = bestPlane(fmhs.at(subleadcandidate.key()));
         
         if(subleadbestplane == 0){
           subleadbestdedx = subleaddedx0;
         }
         
         if(subleadbestplane == 1){
           subleadbestdedx = subleaddedx1;
         }
         
         if(subleadbestplane == 2){
           subleadbestdedx = subleaddedx2;
         }

         _fMultiPairSubleaddEdx0.push_back(subleaddedx0);
         _fMultiPairSubleaddEdx1.push_back(subleaddedx1);
         _fMultiPairSubleaddEdx2.push_back(subleaddedx2);
         _fMultiPairSubleaddEdx3.push_back(subleadbestdedx);

         _fMultiPairPi0Angle12.push_back(angle12);
         float pi0mass = pi0Mass(leadenergy2,subleadenergy2,angle12);
         _fMultiPairPi0Mass.push_back(pi0mass);
         float pi0Momentum = Altpi0Momentum(leadenergy2,subleadenergy2,angle12); 
         _fMultiPairPi0Momentum.push_back(pi0Momentum);
         _fMultiPairPi0Energy.push_back( std::sqrt(std::pow(pi0mass,2) + std::pow(pi0Momentum,2) ) );
         _fMultiPairPi0Costheta.push_back( pi0Angle(pi0Momentum, leadenergy2, subleadenergy2, leadcandidate, subleadcandidate, angle12) );

         if(fIsMC){

          _fMultiPairSubleadPurity0.push_back(getPurity(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 0)  );
          _fMultiPairSubleadCompleteness0.push_back(getComplete(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 0) );
          _fMultiPairSubleadEnergy0Resolution.push_back(getEnergyRes(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, _fGamma1_E, _fGamma2_E, subleadenergy2*altbias_2, fHitPartAssnTag, 0) );
          
          _fMultiPairSubleadPurity1.push_back(getPurity(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 1)  );
          _fMultiPairSubleadCompleteness1.push_back(getComplete(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 1) );
          _fMultiPairSubleadEnergy1Resolution.push_back(getEnergyRes(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, _fGamma1_E, _fGamma2_E, subleadenergy2*altbias_2, fHitPartAssnTag, 1) );

          _fMultiPairSubleadPurity2.push_back(getPurity(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 2)  );
          _fMultiPairSubleadCompleteness2.push_back(getComplete(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 2) );
          _fMultiPairSubleadEnergy2Resolution.push_back(getEnergyRes(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, _fGamma1_E, _fGamma2_E, subleadenergy2*altbias_2, fHitPartAssnTag, 2) );
          
           _fMultiPairSubleadTrueComposition.push_back( getTrueComposition(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, Gamma1_descendents, Gamma2_descendents, fHitPartAssnTag, 2) );
        }//fIsMC



         

        }//end of loop over SubLeading Showers


    }//end of sub lead vector > 1  

   }//sub lead vector > 0
   
   //Finally, look at the other pairs which have a high invariant mass
   if(SubleadVectorHighMass.size() == 1 && nchargedpicand == 0 && nuShowers.size() >= 2 && npairs == 0  && notherfancypairs == 0){
     int subleadtrackID = -999;
     double subleadscore = -999.9;
     double subleadlength = -999.9;
     double subleadopenangle = -999.9;
     double subleadstartx = -999.9;
     double subleadstarty = -999.9;
     double subleadstartz = -999.9;
     double subleaddirx =  -999.9;
     double subleaddiry = -999.9;
     double subleaddirz = -999.9;
     double subleaddist3d = -999.9;
     double subleadang3d = -999.9;
     double subleadenergy0 = -999.9;
     double subleadenergy1 = -999.9;
     double subleadenergy2 = -999.9;
     double subleaddedx0 = -999.9; 
     double subleaddedx1 = -999.9; 
     double subleaddedx2 = -999.9; 
     double subleaddedxamalg = -999.9;
     float subleadbestdedx = -999.9;
     float subleadphi        = -999.9;
     int subleadbestplane = -999;

     auto subleadcandidate = SubleadVectorHighMass.front();
     subleadtrackID = nuShowerTrackMap[subleadcandidate];
     subleadscore = nuShowers_ScoreMap[subleadcandidate];
     subleadlength = subleadcandidate->Length();
     subleadopenangle = radToDeg(subleadcandidate->OpenAngle());

     subleadstartx = subleadcandidate->ShowerStart().X();
     subleadstarty = subleadcandidate->ShowerStart().Y();
     subleadstartz = subleadcandidate->ShowerStart().Z();
     subleadphi    = std::atan2(subleadstartx, subleadstarty);
     subleaddirx =  subleadcandidate->Direction().X();
     subleaddiry = subleadcandidate->Direction().Y();
     subleaddirz = subleadcandidate->Direction().Z();
     subleaddist3d = conversionDistance(subleadcandidate, VertexX, VertexY, VertexZ);
     subleadang3d = radialAngle(subleadcandidate, VertexX, VertexY, VertexZ);
     subleadenergy0 = (showerEnergy(0, fmhs.at(subleadcandidate.key())));
     subleadenergy1 = (showerEnergy(1, fmhs.at(subleadcandidate.key())));
     subleadenergy2 = (showerEnergy(2, fmhs.at(subleadcandidate.key())))/altbias_2;
     getdEdx(e, fCaloTag, fShowerTag, fPFParticleTag, fTrackTag, fTrackFitterTag, fClusterTag, subleadcandidate, nuPFParticles, TrackFitterVector, subleadtrackID, trunklength2, trunkwidth, subleaddedx0, subleaddedx1, subleaddedx2,subleaddedxamalg);
    
     subleadbestplane = bestPlane(fmhs.at(subleadcandidate.key()));
     
     if(subleadbestplane == 0){
       subleadbestdedx = subleaddedx0;
     }
     
     if(subleadbestplane == 1){
       subleadbestdedx = subleaddedx1;
     }
     
     if(subleadbestplane == 2){
       subleadbestdedx = subleaddedx2;
     }

     float leadpurity0    = -999.9;
     float leadcomplete0  = -999.9;
     float leadenergyres0 = -999.9;
     float leadpurity1    = -999.9;
     float leadcomplete1  = -999.9;
     float leadenergyres1 = -999.9;
     float leadpurity2    = -999.9;
     float leadcomplete2  = -999.9;
     float leadenergyres2 = -999.9;
     float leadcomp       = -999.9;
     
     float subleadpurity0    = -999.9;
     float subleadcomplete0  = -999.9;
     float subleadenergyres0 = -999.9;
     float subleadpurity1    = -999.9;
     float subleadcomplete1  = -999.9;
     float subleadenergyres1 = -999.9;
     float subleadpurity2    = -999.9;
     float subleadcomplete2  = -999.9;
     float subleadenergyres2 = -999.9;
     float subleadcomp       = -999.9;


     double angle12 = angleBetweenTwoShowers(leadcandidate,subleadcandidate);
     double acosangle12 = radToDeg(std::acos(angle12));
     double pi0mass = -99999.;
     double pi0momentum = -99999.;
     double pi0costheta = -99999.;
     double pi0energy = -99999.;
    
     pi0mass = pi0Mass(leadenergy2,subleadenergy2,angle12);
     pi0momentum = Altpi0Momentum(leadenergy2,subleadenergy2,angle12);
     pi0energy = std::sqrt(std::pow(pi0mass,2) + std::pow(pi0momentum,2));
     pi0costheta = pi0Angle(pi0momentum, leadenergy2, subleadenergy2, leadcandidate, subleadcandidate, angle12);

     if(fIsMC){
        subleadpurity0 = getPurity(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 0);
        subleadcomplete0 = getComplete(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 0);
        subleadenergyres0 = getEnergyRes(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, _fGamma1_E, _fGamma2_E, subleadenergy2*altbias_2, fHitPartAssnTag, 0);
        subleadpurity1 = getPurity(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 1);
        subleadcomplete1 = getComplete(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 1);
        subleadenergyres1 = getEnergyRes(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, _fGamma1_E, _fGamma2_E, subleadenergy2*altbias_2, fHitPartAssnTag, 1);
        subleadpurity2 = getPurity(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 2);
        subleadcomplete2 = getComplete(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, Gamma1_descendents, Gamma2_descendents, _fGamma1_id, _fGamma2_id, fHitPartAssnTag, 2);
        subleadenergyres2 = getEnergyRes(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, _fGamma1_E, _fGamma2_E, subleadenergy2*altbias_2, fHitPartAssnTag, 2);
        subleadcomp = getTrueComposition(e, fHitTag, HitVector, fmhs.at(subleadcandidate.key()), MCParticleVector, _fGamma1_id, _fGamma2_id, Gamma1_descendents, Gamma2_descendents, fHitPartAssnTag, 2);
      }//fIsMC

       _fHiMassLeadScore= leadscore;
       _fHiMassLeadLength   = leadlength;
       _fHiMassLeadOpenAngle= leadopenangle;
       _fHiMassLeadCostheta = leaddirz;
       _fHiMassLeadPhi      = std::atan2(leaddiry, leaddirx);
       _fHiMassLeadDirZ     = leaddirz;
       _fHiMassLeadDirX     = leaddirx;
       _fHiMassLeadDirY     = leaddiry;
       _fHiMassLeadDist3d   = leaddist3d;
       _fHiMassLeadAng3d    = leadang3d;
       _fHiMassLeadStartX   = leadstartx;
       _fHiMassLeadStartY   = leadstarty;
       _fHiMassLeadStartZ   = leadstartz;
       _fHiMassLeadEnergy0  = leadenergy0;
       _fHiMassLeadEnergy1  = leadenergy1;
       _fHiMassLeadEnergy2  = leadenergy2;
       _fHiMassLeaddEdx0    = leaddedx0;
       _fHiMassLeaddEdx1    = leaddedx1;
       _fHiMassLeaddEdx2    = leaddedx2;
       _fHiMassLeaddEdx3    = leadbestdedx;

       _fHiMassSubleadScore = subleadscore;
       _fHiMassSubleadLength = subleadlength;
       _fHiMassSubleadOpenAngle = subleadopenangle;
       _fHiMassSubleadCostheta = subleaddirz;
       _fHiMassSubleadPhi = subleadphi;
       _fHiMassSubleadDirZ     = subleaddirz;
       _fHiMassSubleadDirX     = subleaddirx;
       _fHiMassSubleadDirY     = subleaddiry;
       _fHiMassSubleadDist3d = subleaddist3d;
       _fHiMassSubleadAng3d = subleadang3d;
       _fHiMassSubleadStartX = subleadstartx;
       _fHiMassSubleadStartY = subleadstarty;
       _fHiMassSubleadStartZ = subleadstartz;
       _fHiMassSubleadEnergy0 = subleadenergy0;
       _fHiMassSubleadEnergy1 = subleadenergy1;
       _fHiMassSubleadEnergy2 = subleadenergy2;
       _fHiMassSubleaddEdx0 = subleaddedx0;
       _fHiMassSubleaddEdx1 = subleaddedx1;
       _fHiMassSubleaddEdx2 = subleaddedx2;
       _fHiMassSubleaddEdx3 = subleadbestdedx;

       //Candidate pi0
       _fHiMassPi0Momentum = pi0momentum;
       _fHiMassPi0Energy = pi0energy;
       _fHiMassPi0Costheta = pi0costheta;
       _fHiMassPi0Angle12 = acosangle12;
       _fHiMassPi0Mass     = pi0mass;

       _fHiMassLeadTrueComposition = leadcomp;
       _fHiMassLeadEnergy0Resolution = leadenergyres0;
       _fHiMassLeadEnergy1Resolution = leadenergyres1;
       _fHiMassLeadEnergy2Resolution = leadenergyres2;
       _fHiMassLeadPurity0 = leadpurity0;
       _fHiMassLeadPurity1 = leadpurity1;
       _fHiMassLeadPurity2 = leadpurity2;
       _fHiMassLeadCompleteness0 = leadcomplete0;
       _fHiMassLeadCompleteness1 = leadcomplete1;
       _fHiMassLeadCompleteness2 = leadcomplete2;

       _fHiMassSubleadTrueComposition = subleadcomp;
       _fHiMassSubleadEnergy0Resolution = subleadenergyres0;
       _fHiMassSubleadEnergy1Resolution = subleadenergyres1;
       _fHiMassSubleadEnergy2Resolution = subleadenergyres2;
       _fHiMassSubleadPurity0 = subleadpurity0;
       _fHiMassSubleadPurity1 = subleadpurity1;
       _fHiMassSubleadPurity2 = subleadpurity2;
       _fHiMassSubleadCompleteness0 = subleadcomplete0;
       _fHiMassSubleadCompleteness1 = subleadcomplete1;
       _fHiMassSubleadCompleteness2 = subleadcomplete2;   
   }

  }//passes nucc precuts

  _eventtree->Fill();
 
}//analyze

//////////////////////////////////////////////////////////////////////
/////////////////////// ALL FUNCTIONS HERE ///////////////////////////
//////////////////////////////////////////////////////////////////////

void SidebandTree::checkMatchedGammas(art::Event const & e, art::InputTag allhits_tag, std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> shower_hits, art::InputTag hitparticleassns_tag, int plane, int Gamma_id, std::vector<art::Ptr<simb::MCParticle>> TwoGammaV, int& samedaughters, int& otherdaughters){

  std::vector<int> thisgammadaughterids;
  std::vector<int> othergammadaughterids;
  for(auto const& mcpart: TwoGammaV){
    if(mcpart->TrackId() == Gamma_id){
      for(int i = 0; i < mcpart->NumberDaughters(); i++){
  thisgammadaughterids.push_back(mcpart->Daughter(i));
      }
    }//got right gamma
    if(mcpart->TrackId() != Gamma_id){
      for(int i = 0; i < mcpart->NumberDaughters(); i++){
  othergammadaughterids.push_back(mcpart->Daughter(i));
      }
    }//other gamma
  }//both gammas
  
  art::ServiceHandle<cheat::ParticleInventoryService> particleinventory;
  auto const& allhits_handle = e.getValidHandle<std::vector<recob::Hit>>(allhits_tag);

  //////////////////////////////////////////////////////
  ////////////// SHOWER HITS ///////////////////////////
  //////////////////////////////////////////////////////
  std::map<int,double> trkID_E; //map of track ID's and their deposited energies
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(allhits_handle,e,hitparticleassns_tag); //mc particles associated with all hits
  std::vector<simb::MCParticle const*> particle_vec; //vector of mc particles 
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec; //vector of associations
  
  for(size_t i = 0; i < shower_hits.size(); i++){
    art::Ptr<recob::Hit> hit = shower_hits[i];
    auto hitplane = hit->View();
    if(hitplane != plane) continue;

    particle_vec.clear(); match_vec.clear();
    particles_per_hit.get(hit.key(),particle_vec,match_vec); //mc particles matched to each shower hit

    //ALL MCPARTICLES MATCHED TO CURRENT SHOWER HIT
    for(size_t i_p=0; i_p < particle_vec.size(); ++i_p){
      trkID_E[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id for all shower-matched mc particles
    }//mcparticles per hit

  }//shower hit loop

  int matchestothis = 0;
  int matchestoother = 0;

  for(std::map<int,double>::iterator it = trkID_E.begin(); it!=trkID_E.end(); ++it){
    for(size_t i = 0; i < thisgammadaughterids.size(); i++){
      if(it->first == thisgammadaughterids.at(i)){
  matchestothis++;
      }
    }
    for(size_t i = 0; i < othergammadaughterids.size(); i++){
      if(it->first == othergammadaughterids.at(i)){
  matchestoother++;
      }
    }
  }

  samedaughters = matchestothis;
  otherdaughters = matchestoother;

}//checkForElectronPair

double SidebandTree::getTrueComposition(art::Event const & e, art::InputTag allhits_tag, std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> shower_hits, std::vector<art::Ptr<simb::MCParticle>> all_mcparticles, int Gamma1_id, int Gamma2_id, std::vector<art::Ptr<simb::MCParticle>> Gamma1_descendents, std::vector<art::Ptr<simb::MCParticle>> Gamma2_descendents, art::InputTag hitparticleassns_tag, int plane){
  
  int truecomp = -99999; //1 = gamma1, 2 = gamma2, 3 = descendent, 4 = muon, 5 = proton, 6 = pion, 7 = other em, 8 = overlay, 9 = the rest

  int mcparticleid = -999;
  float mcenergy = -999.;
  float overlayenergy = -999.;
  int pdgcode = -999;
  int momid = -999;
  bool matchGamma1 = false;
  bool matchGamma2 = false;
  bool descendentGamma1 = false;
  bool descendentGamma2 = false;
  
  if(fIsMC){

    getTrueParticle(e, allhits_tag, all_hits, shower_hits, hitparticleassns_tag, 2, mcparticleid, mcenergy);
    getOverlayEnergy(e, allhits_tag, all_hits, shower_hits, hitparticleassns_tag, 2, overlayenergy);
  
    for(auto const& mcpart : all_mcparticles){
      if(mcparticleid == mcpart->TrackId()){
  momid = mcpart->Mother();
  pdgcode = mcpart->PdgCode();
      }
    }//mcparticle loop

    if(momid == Gamma1_id){
      matchGamma1 = true;
    }
    if(momid == Gamma2_id){
      matchGamma2 = true;
    }

    //mostly overlay
    if(overlayenergy > mcenergy){
      truecomp = 8;
    }//overlay
  
    if(overlayenergy <= mcenergy){
    
      if(matchGamma1 == true){
  truecomp = 1;
      }//gamma1
    
      if(matchGamma2 == true){
  truecomp = 2;
      }//gamma2

      if(matchGamma1 == false && matchGamma2 == false){
  for(auto const & mcpart : Gamma1_descendents){
    if(mcparticleid == mcpart->TrackId()){descendentGamma1 = true;}
  }
  for(auto const & mcpart : Gamma2_descendents){
    if(mcparticleid == mcpart->TrackId()){descendentGamma2 = true;}
  }
      }

      if(descendentGamma1 == true || descendentGamma2 == true){
  truecomp = 3;
      }
    
      if(std::abs(pdgcode) == 13){
  truecomp = 4;
      }
    
      if(std::abs(pdgcode) == 2212){
  truecomp = 5;
      }
    
      if(std::abs(pdgcode) == 211){
  truecomp = 6;
      }
    
      if((std::abs(pdgcode) == 11) && matchGamma1 == false && matchGamma2 == false && descendentGamma1 == false && descendentGamma2 == false){
  truecomp = 7;
      }

      if(matchGamma1 == false && matchGamma2 == false && descendentGamma1 == false && descendentGamma2 == false && std::abs(pdgcode) != 13 && std::abs(pdgcode) != 2212 && std::abs(pdgcode) != 211 && std::abs(pdgcode) != 11){
  truecomp = 9;
      }
    
    }//not overlay


  }//fIsMC

  return truecomp;

}//getTrueComposition

int SidebandTree::bestPlane(std::vector<art::Ptr<recob::Hit>> hits){

  if(hits.size() == 0){std::cout<<"SHOWER HAS NO HITS!"<<std::endl;}

  int plane = -999;
  
  int nhitsU = 0;
  int nhitsV = 0;
  int nhitsY = 0;
  int nhitsMax = -999;
  
  for(auto const& ihit : hits){
    auto hitplane = ihit->View();
    if(hitplane == 0){nhitsU++;}
    if(hitplane == 1){nhitsV++;}
    if(hitplane == 2){nhitsY++;}
  }//shower hits

  nhitsMax = std::max(nhitsU,std::max(nhitsV,nhitsY));

  if(nhitsMax == nhitsU){plane = 0;}
  if(nhitsMax == nhitsV){plane = 1;}
  if(nhitsMax == nhitsY){plane = 2;}

  return plane;
  
}//bestPlane

///////////////////// CHECK TRUTH //////////////////////////////

double SidebandTree::getTrueTrackComposition(art::Event const & e, art::InputTag allhits_tag, std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> track_hits, std::vector<art::Ptr<simb::MCParticle>> all_mcparticles, art::InputTag hitparticleassns_tag, int plane){
  
  int truecomp = -99999; //1 = muon, 2 = proton, 3 = charged pion, 4 = photon, 5 = other em, 6 = overlay, 7 = other

  int mcparticleid = -999;
  float mcenergy = -999.;
  float overlayenergy = -999.;
  int pdgcode = -999;
  int mompdgcode = -999;
  int momid = -999;
 
  if(fIsMC){

    getTrueParticle(e, allhits_tag, all_hits, track_hits, hitparticleassns_tag, 2, mcparticleid, mcenergy);
    getOverlayEnergy(e, allhits_tag, all_hits, track_hits, hitparticleassns_tag, 2, overlayenergy);
  
    for(auto const& mcpart : all_mcparticles){
      if(mcparticleid == mcpart->TrackId()){
  momid = mcpart->Mother();
  pdgcode = mcpart->PdgCode();
  for(auto const& mcpart2 : all_mcparticles){
    if(mcpart2->TrackId() == momid){
      mompdgcode = mcpart2->PdgCode();
    }
  }
      }
    }//mcparticle loop

    //mostly overlay
    if(overlayenergy > mcenergy){
      truecomp = 6;
    }//overlay
  
    if(overlayenergy <= mcenergy){

      if(std::abs(pdgcode) == 13){
  truecomp = 1;
      }

      if(std::abs(pdgcode) == 2212){
  truecomp = 2;
      }

      if(std::abs(pdgcode) == 211){
  truecomp = 3;
      }
    
      if(std::abs(mompdgcode) == 22){
  truecomp = 4;
      }
    
      if(std::abs(pdgcode) == 11){
  truecomp = 5;
      }
    
      if(std::abs(pdgcode) != 13 && std::abs(pdgcode) != 11 && std::abs(pdgcode) != 2212 && std::abs(pdgcode) != 211 && std::abs(pdgcode) != 22){
  truecomp = 7;
      }
    
    }//not overlay


  }//fIsMC

  return truecomp;

}//getTrueTrackComposition



double SidebandTree::getEnergyRes(art::Event const & e, art::InputTag allhits_tag, std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> shower_hits, std::vector<art::Ptr<simb::MCParticle>> all_mcparticles, int Gamma1_id, int Gamma2_id, float Gamma1_E, float Gamma2_E, float energy, art::InputTag hitparticleassns_tag, int plane){

  art::ServiceHandle<cheat::ParticleInventoryService> particleinventory;
  auto const& allhits_handle = e.getValidHandle<std::vector<recob::Hit>>(allhits_tag);

  double energyres = -99999;

  int shower_mcparticleid = -999; //mc particle id that contributes most energy to the shower
  int mcparticle_id = -999; //mcparticle that contributes most energy to shower
  float max_e = -999.; //max true energy deposition in the shower
  
  //////////////////////////////////////////////////////
  ////////////// SHOWER HITS ///////////////////////////
  //////////////////////////////////////////////////////
  std::map<int,double> trkID_E; //map of track ID's and their deposited energies
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(allhits_handle,e,hitparticleassns_tag); //mc particles associated with all hits
  std::vector<simb::MCParticle const*> particle_vec; //vector of mc particles 
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec; //vector of associations
  
  for(size_t i = 0; i < shower_hits.size(); i++){
    art::Ptr<recob::Hit> hit = shower_hits[i];
    auto hitplane = hit->View();
    if(hitplane != plane) continue;

    particle_vec.clear(); match_vec.clear();
    particles_per_hit.get(hit.key(),particle_vec,match_vec); //mc particles matched to each shower hit

    //ALL MCPARTICLES MATCHED TO CURRENT SHOWER HIT
    for(size_t i_p=0; i_p < particle_vec.size(); ++i_p){

      trkID_E[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id for all shower-matched mc particles
     
      //Particle with max energy contribution to hit
      if( trkID_E[particle_vec[i_p]->TrackId()] > max_e ){ 
  max_e = trkID_E[ particle_vec[i_p]->TrackId() ];
  mcparticle_id = particle_vec[i_p]->TrackId();
      }//max energy particle
      
    }//mcparticles per hit
    
  }//shower hit loop

  shower_mcparticleid = mcparticle_id;

  ///////////////////////////////////////////////////////////////////////////////////////////////
  int momid = -999;
  for(auto const& mcpart : all_mcparticles){
    if(shower_mcparticleid == mcpart->TrackId()){
      momid = mcpart->Mother();
    }
  }//mcparticle loop

  //matches to gamma1
  if(momid == Gamma1_id){
    if(Gamma1_E != 0.){energyres = ((energy/altbias_2) - Gamma1_E)/Gamma1_E;}
  }//g1

  //matches to gamma2
  if(momid == Gamma2_id){
    if(Gamma2_E != 0.){energyres = ((energy/altbias_2) - Gamma2_E)/Gamma2_E;}
  }//g2

  return energyres;
  
}//getEnergyRes


double SidebandTree::getTrueEnergy(art::Event const & e, art::InputTag allhits_tag, std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> shower_hits, std::vector<art::Ptr<simb::MCParticle>> all_mcparticles, int Gamma1_id, int Gamma2_id, float Gamma1_E, float Gamma2_E, art::InputTag hitparticleassns_tag, int plane){

  art::ServiceHandle<cheat::ParticleInventoryService> particleinventory;
  auto const& allhits_handle = e.getValidHandle<std::vector<recob::Hit>>(allhits_tag);

  double trueenergy = -99999;

  int shower_mcparticleid = -999; //mc particle id that contributes most energy to the shower
  int mcparticle_id = -999; //mcparticle that contributes most energy to shower
  float max_e = -999.; //max true energy deposition in the shower
  
  //////////////////////////////////////////////////////
  ////////////// SHOWER HITS ///////////////////////////
  //////////////////////////////////////////////////////
  std::map<int,double> trkID_E; //map of track ID's and their deposited energies
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(allhits_handle,e,hitparticleassns_tag); //mc particles associated with all hits
  std::vector<simb::MCParticle const*> particle_vec; //vector of mc particles 
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec; //vector of associations
  
  for(size_t i = 0; i < shower_hits.size(); i++){
    art::Ptr<recob::Hit> hit = shower_hits[i];
    auto hitplane = hit->View();
    if(hitplane != plane) continue;

    particle_vec.clear(); match_vec.clear();
    particles_per_hit.get(hit.key(),particle_vec,match_vec); //mc particles matched to each shower hit

    //ALL MCPARTICLES MATCHED TO CURRENT SHOWER HIT
    for(size_t i_p=0; i_p < particle_vec.size(); ++i_p){

      trkID_E[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id for all shower-matched mc particles
     
      //Particle with max energy contribution to hit
      if( trkID_E[particle_vec[i_p]->TrackId()] > max_e ){ 
  max_e = trkID_E[ particle_vec[i_p]->TrackId() ];
  mcparticle_id = particle_vec[i_p]->TrackId();
      }//max energy particle
      
    }//mcparticles per hit
    
  }//shower hit loop

  shower_mcparticleid = mcparticle_id;

  ///////////////////////////////////////////////////////////////////////////////////////////////
  int momid = -999;
  for(auto const& mcpart : all_mcparticles){
    if(shower_mcparticleid == mcpart->TrackId()){
      momid = mcpart->Mother();
    }
  }//mcparticle loop

  //matches to gamma1
  if(momid == Gamma1_id){
    trueenergy = Gamma1_E;
  }//g1

  //matches to gamma2
  if(momid == Gamma2_id){
    trueenergy = Gamma2_E;
  }//g2

  return trueenergy;
  
}//getTrueEnergy
 
double SidebandTree::getPurity(art::Event const & e, art::InputTag allhits_tag, std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> shower_hits, std::vector<art::Ptr<simb::MCParticle>> all_mcparticles, std::vector<art::Ptr<simb::MCParticle>> Gamma1_descendents, std::vector<art::Ptr<simb::MCParticle>> Gamma2_descendents, int Gamma1_id, int Gamma2_id, art::InputTag hitparticleassns_tag, int plane){

  art::ServiceHandle<cheat::ParticleInventoryService> particleinventory;
  auto const& allhits_handle = e.getValidHandle<std::vector<recob::Hit>>(allhits_tag);

  double purity = -99999;

  int shower_mcparticleid = -999; //mc particle id that contributes most energy to the shower
  float allmc_shower = -999.;
  float gamma1_shower = -999.;
  float gamma2_shower = -999.;
 
  int mcparticle_id = -999; //mcparticle that contributes most energy to shower
  float max_e = -999.; //max true energy deposition in the shower
  
  float mconly_shower = 0.; 
  float g1_shower = 0.; //energy deposited by gamma1 in the shower
  float g2_shower = 0.; //energy deposited by gamma2 in the shower
  
  std::vector<art::Ptr<recob::Hit>> overlayHits;
    
  //////////////////////////////////////////////////////
  ////////////// SHOWER HITS ///////////////////////////
  //////////////////////////////////////////////////////
  std::map<int,double> trkID_E; //map of track ID's and their deposited energies
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(allhits_handle,e,hitparticleassns_tag); //mc particles associated with all hits
  std::vector<simb::MCParticle const*> particle_vec; //vector of mc particles 
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec; //vector of associations
  
  for(size_t i = 0; i < shower_hits.size(); i++){
    art::Ptr<recob::Hit> hit = shower_hits[i];
    auto hitplane = hit->View();
    if(hitplane != plane) continue;

    particle_vec.clear(); match_vec.clear();
    particles_per_hit.get(hit.key(),particle_vec,match_vec); //mc particles matched to each shower hit

    if(particle_vec.size() == 0){
      overlayHits.push_back(hit);
    }

    //ALL MCPARTICLES MATCHED TO CURRENT SHOWER HIT
    for(size_t i_p=0; i_p < particle_vec.size(); ++i_p){

      trkID_E[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id for all shower-matched mc particles
     
      //Particle with max energy contribution to hit
      if( trkID_E[particle_vec[i_p]->TrackId()] > max_e ){ 
  max_e = trkID_E[ particle_vec[i_p]->TrackId() ];
  mcparticle_id = particle_vec[i_p]->TrackId();
      }//max energy particle
      
    }//mcparticles per hit
    
  }//shower hit loop

  //////////////////////////////////////////////////////
  ////////// ALL HITS IN THE EVENT /////////////////////
  //////////////////////////////////////////////////////
  std::map<int,double> trkID_Eall; //map of track ID's and their deposited energies
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hitall(allhits_handle,e,hitparticleassns_tag); //mc particles associated with all hits
  std::vector<simb::MCParticle const*> particle_vecall; //vector of mc particles 
  std::vector<anab::BackTrackerHitMatchingData const*> match_vecall; //vector of associations
  
  for(size_t i = 0; i < all_hits.size(); i++){
    art::Ptr<recob::Hit> hit = all_hits[i];
    auto hitplane = hit->View();
    if(hitplane != plane) continue;
 
    particle_vecall.clear(); match_vecall.clear();
    particles_per_hitall.get(hit.key(),particle_vecall,match_vecall); //mc particles matched to each event hit

    //ALL MCPARTICLES MATCHED TO CURRENT EVENT HIT
    for(size_t i_p=0; i_p < particle_vecall.size(); ++i_p){

      trkID_Eall[ particle_vecall[i_p]->TrackId() ] += match_vecall[i_p]->energy; //store energy per track id for all matched mc particles
    
    }//mcparticles per hit

  }//all hit loop

  ///////////////////////////////////////////////////
  ////////// GET CONTRIBUTIONS //////////////////////
  ///////////////////////////////////////////////////

  /////////////// SHOWER PARTICLE MAP //////////
  for(std::map<int,double>::iterator it = trkID_E.begin(); it!=trkID_E.end(); ++it){

    mconly_shower += it->second;
    
    //Check how much energy gamma1 deposited in the shower
    for(auto const& mcpart : Gamma1_descendents){
      if(mcpart->TrackId() == it->first){
  g1_shower += it->second;
      }
    }//g1

    //Check how much energy gamma2 deposited in the shower
    for(auto const& mcpart : Gamma2_descendents){
      if(mcpart->TrackId() == it->first){
  g2_shower += it->second;
      }
    }//g2
    
  }//trkide

  float overlay_shower = showerEnergy(2,overlayHits);

  allmc_shower = overlay_shower + mconly_shower;
  gamma1_shower = g1_shower;
  gamma2_shower = g2_shower;
  shower_mcparticleid = mcparticle_id;

  ///////////////////////////////////////////////////////////////////////////////////////////////
  int momid = -999;
  for(auto const& mcpart : all_mcparticles){
    if(shower_mcparticleid == mcpart->TrackId()){
      momid = mcpart->Mother();
    }
  }//mcparticle loop

  //matches to gamma1
  if(momid == Gamma1_id){
    if(allmc_shower != 0.){purity = gamma1_shower/allmc_shower;}
  }//g1

  //matches to gamma2
  if(momid == Gamma2_id){
    if(allmc_shower != 0.){purity = gamma2_shower/allmc_shower;}
  }//g2

  return purity;
  
}//getPurity

double SidebandTree::getComplete(art::Event const & e, art::InputTag allhits_tag, std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> shower_hits, std::vector<art::Ptr<simb::MCParticle>> all_mcparticles, std::vector<art::Ptr<simb::MCParticle>> Gamma1_descendents, std::vector<art::Ptr<simb::MCParticle>> Gamma2_descendents, int Gamma1_id, int Gamma2_id, art::InputTag hitparticleassns_tag, int plane){

  art::ServiceHandle<cheat::ParticleInventoryService> particleinventory;
  auto const& allhits_handle = e.getValidHandle<std::vector<recob::Hit>>(allhits_tag);

  double complete = -99999;

  int shower_mcparticleid = -999; //mc particle id that contributes most energy to the shower
  float gamma1_shower = -999.;
  float gamma2_shower = -999.;
  float gamma1_event = -999.;
  float gamma2_event = -999.;
  
  int mcparticle_id = -999; //mcparticle that contributes most energy to shower
  float max_e = -999.; //max true energy deposition in the shower
  
  float g1_shower = 0.; //energy deposited by gamma1 in the shower
  float g2_shower = 0.; //energy deposited by gamma2 in the shower
  float g1_event = 0.; //energy deposited by gamma1 in the event
  float g2_event = 0.; //energy deposited by gamma2 in the event

  std::vector<art::Ptr<recob::Hit>> overlayHits;
    
  //////////////////////////////////////////////////////
  ////////////// SHOWER HITS ///////////////////////////
  //////////////////////////////////////////////////////
  std::map<int,double> trkID_E; //map of track ID's and their deposited energies
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(allhits_handle,e,hitparticleassns_tag); //mc particles associated with all hits
  std::vector<simb::MCParticle const*> particle_vec; //vector of mc particles 
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec; //vector of associations
  
  for(size_t i = 0; i < shower_hits.size(); i++){
    art::Ptr<recob::Hit> hit = shower_hits[i];
    auto hitplane = hit->View();
    if(hitplane != plane) continue;

    particle_vec.clear(); match_vec.clear();
    particles_per_hit.get(hit.key(),particle_vec,match_vec); //mc particles matched to each shower hit

    //ALL MCPARTICLES MATCHED TO CURRENT SHOWER HIT
    for(size_t i_p=0; i_p < particle_vec.size(); ++i_p){

      trkID_E[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id for all shower-matched mc particles
     
      //Particle with max energy contribution to hit
      if( trkID_E[particle_vec[i_p]->TrackId()] > max_e ){ 
  max_e = trkID_E[ particle_vec[i_p]->TrackId() ];
  mcparticle_id = particle_vec[i_p]->TrackId();
      }//max energy particle
      
    }//mcparticles per hit
    
  }//shower hit loop

  //////////////////////////////////////////////////////
  ////////// ALL HITS IN THE EVENT /////////////////////
  //////////////////////////////////////////////////////
  std::map<int,double> trkID_Eall; //map of track ID's and their deposited energies
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hitall(allhits_handle,e,hitparticleassns_tag); //mc particles associated with all hits
  std::vector<simb::MCParticle const*> particle_vecall; //vector of mc particles 
  std::vector<anab::BackTrackerHitMatchingData const*> match_vecall; //vector of associations
  
  for(size_t i = 0; i < all_hits.size(); i++){
    art::Ptr<recob::Hit> hit = all_hits[i];
    auto hitplane = hit->View();
    if(hitplane != plane) continue;
 
    particle_vecall.clear(); match_vecall.clear();
    particles_per_hitall.get(hit.key(),particle_vecall,match_vecall); //mc particles matched to each event hit

    //ALL MCPARTICLES MATCHED TO CURRENT EVENT HIT
    for(size_t i_p=0; i_p < particle_vecall.size(); ++i_p){

      trkID_Eall[ particle_vecall[i_p]->TrackId() ] += match_vecall[i_p]->energy; //store energy per track id for all matched mc particles
    
    }//mcparticles per hit

  }//all hit loop

  ///////////////////////////////////////////////////
  ////////// GET CONTRIBUTIONS //////////////////////
  ///////////////////////////////////////////////////

  /////////////// SHOWER PARTICLE MAP //////////
  for(std::map<int,double>::iterator it = trkID_E.begin(); it!=trkID_E.end(); ++it){

    //Check how much energy gamma1 deposited in the shower
    for(auto const& mcpart : Gamma1_descendents){
      if(mcpart->TrackId() == it->first){
  g1_shower += it->second;
      }
    }//g1

    //Check how much energy gamma2 deposited in the shower
    for(auto const& mcpart : Gamma2_descendents){
      if(mcpart->TrackId() == it->first){
  g2_shower += it->second;
      }
    }//g2
    
  }//trkide

   /////////////// ALL PARTICLE MAP //////////
  for(std::map<int,double>::iterator it = trkID_Eall.begin(); it!=trkID_Eall.end(); ++it){

    //Check how much energy gamma1 deposited in the shower
    for(auto const& mcpart : Gamma1_descendents){
      if(mcpart->TrackId() == it->first){
  g1_event += it->second;
      }
    }//g1

    //Check how much energy gamma2 deposited in the shower
    for(auto const& mcpart : Gamma2_descendents){
      if(mcpart->TrackId() == it->first){
  g2_event += it->second;
      }
    }//g2
    
  }//trkide

  gamma1_shower = g1_shower;
  gamma2_shower = g2_shower;
  gamma1_event = g1_event;
  gamma2_event = g2_event;
  shower_mcparticleid = mcparticle_id;

  ///////////////////////////////////////////////////////////////////////////////////////////////
  int momid = -999;
  for(auto const& mcpart : all_mcparticles){
    if(shower_mcparticleid == mcpart->TrackId()){
      momid = mcpart->Mother();
    }
  }//mcparticle loop

  //matches to gamma1
  if(momid == Gamma1_id){
    if(gamma1_event != 0.){complete = gamma1_shower/gamma1_event;}
  }//g1

  //matches to gamma2
  if(momid == Gamma2_id){
    if(gamma2_event != 0.){complete = gamma2_shower/gamma2_event;}
  }//g2

  return complete;
  
}//getComplete


void SidebandTree::getdEdx(art::Event const & e, art::InputTag calo_tag, art::InputTag shower_tag, art::InputTag pfp_tag, art::InputTag track_tag, art::InputTag trackfitter_tag, art::InputTag cluster_tag, art::Ptr<recob::Shower> shower, std::vector<art::Ptr<recob::PFParticle>> nuPFParticles, std::vector<art::Ptr<recob::Track>> TrackFitterVector, int trackID, float trunklength, float trunkwidth, double& dedx0, double& dedx1, double& dedx2, double& amalgdedx){

  auto const& trackfitter_handle = e.getValidHandle<std::vector<recob::Track>>(trackfitter_tag);
  auto const& track_handle = e.getValidHandle<std::vector<recob::Track>>(track_tag);
  auto const& shower_handle = e.getValidHandle<std::vector<recob::Shower>>(shower_tag);
  auto const& cluster_handle = e.getValidHandle<std::vector<recob::Cluster>>(cluster_tag);
  auto const& pfp_handle = e.getValidHandle<std::vector<recob::PFParticle>>(pfp_tag);
  
  art::FindManyP<anab::Calorimetry> TrackCaloAssoc(trackfitter_handle,e,calo_tag);
  art::FindManyP<recob::Shower> pfPartToShowerAssoc(pfp_handle,e,"pandora");
  art::FindManyP<recob::Cluster> pfPartToClusterAssoc(pfp_handle,e,"pandora");
  art::FindManyP<recob::Hit> fmhs(shower_handle,e,shower_tag);
  art::FindManyP<recob::Hit> fmhc(cluster_handle,e,cluster_tag);
  art::FindManyP<recob::Hit> fmht(track_handle,e,track_tag);

  std::cout<<trunklength<<trunkwidth<<std::endl;
  
  dedx0 = -999;
  dedx1 = -999;
  dedx2 = -999;
  amalgdedx = -999;
  
  float track_dedx0 = -99999.;
  float track_dedx1 = -99999.;
  float track_dedx2 = -99999.;

  //track fitter dedx
  bool matchedToTrackFitter = false;
  //int trackID = -9999;
  //trackID = nuShowerTrackMap[shower];
  std::cout<<trackID<<std::endl;
  for(auto const& track : TrackFitterVector){
    if(track->ID() != trackID) continue;
    matchedToTrackFitter = true;
    TrackFitdEdx(TrackCaloAssoc.at(track.key()), trunklength, track_dedx0, track_dedx1, track_dedx2);
  }//track fitter loop
  
  //shower dedx
  double amalgamated_dedx = -999.;
  std::vector<double> shower_dQdx_0;
  std::vector<double> shower_dQdx_1;
  std::vector<double> shower_dQdx_2;
  std::vector<double> shower_dEdx_0;
  std::vector<double> shower_dEdx_1;
  std::vector<double> shower_dEdx_2;
  double median_dedx_0 = -999.;
  double median_dedx_1 = -999.;    
  double median_dedx_2 = -999.;
  
  TVector3 shr_dir = shower->Direction();
  float angle_0 = getAnglewrtWires(shr_dir, 0);
  float angle_1 = getAnglewrtWires(shr_dir, 1);
  float angle_2 = getAnglewrtWires(shr_dir, 2);
  float dirz = shower->Direction().Z();
  float diru = std::cos(angle_0);
  float dirv =  std::cos(angle_1);
  float acosdirz = std::acos(dirz);
  
  std::cout<<"CHECK ANGLES!!! Theta wrt z = "<<acosdirz<<" wrt u = "<<angle_0<<" wrt v = "<<angle_1<<" cosine theta wrt z = "<<dirz<<" wrt u = "<<diru<<" wrt v = "<<dirv<<std::endl;
  
  //find matched clusters
  std::vector<art::Ptr<recob::Cluster> > matchedClusters;
  for(auto const& pfp : nuPFParticles){
    const auto associatedShowers = pfPartToShowerAssoc.at(pfp.key());
    const auto associatedClusters = pfPartToClusterAssoc.at(pfp.key());
    bool matched_to_shower = false;
    
    if (associatedShowers.size() > 0){
      for(size_t i = 0; i < associatedShowers.size(); i++){
    auto const sh = associatedShowers[i];
    if(sh == shower){
      matched_to_shower = true;
    }//matched
      }//associated showers
    }//has associated showers
    
    if(matched_to_shower == true){
      if (associatedClusters.size() > 0){
    for(size_t i = 0; i < associatedClusters.size(); i++){
      auto const cl = associatedClusters[i];
      matchedClusters.push_back(cl);
    }//associated clusters
      }//has associated clusters
    }//matched
  }//nupfp
  
  //Clusters associated with current shower
  for(auto const& cluster : matchedClusters){
    std::vector<double> dQdx_0;
    std::vector<double> dQdx_1;
    std::vector<double> dQdx_2;
    
    if(cluster->View() == 0){dQdx_0 = CalcdQdxShower(shower, cluster, fmhc.at(cluster.key()), 0, trunklength);}
    if(cluster->View() == 1){dQdx_1 = CalcdQdxShower(shower, cluster, fmhc.at(cluster.key()), 1, trunklength);}
    if(cluster->View() == 2){dQdx_2 = CalcdQdxShower(shower, cluster, fmhc.at(cluster.key()), 2, trunklength);}
    
    //plane 0
    for(size_t i = 0; i < dQdx_0.size(); i++){
      double thisdQdx = dQdx_0.at(i);
      double thisdEdx = QtoEConversion(dQdx_0.at(i));
      shower_dQdx_0.push_back(thisdQdx);
      shower_dEdx_0.push_back(thisdEdx);
    }//0
    
    //plane 1
    for(size_t i = 0; i < dQdx_1.size(); i++){
      double thisdQdx = dQdx_1.at(i);
      double thisdEdx = QtoEConversion(dQdx_1.at(i));
      shower_dQdx_1.push_back(thisdQdx);
      shower_dEdx_1.push_back(thisdEdx);
    }//1
    
    //plane 2
    for(size_t i = 0; i < dQdx_2.size(); i++){
      double thisdQdx = dQdx_2.at(i);
      double thisdEdx = QtoEConversion(dQdx_2.at(i));
      shower_dQdx_2.push_back(thisdQdx);
      shower_dEdx_2.push_back(thisdEdx);
    }//2
    
  }//cluster loop
  
  //calculate output variables
  median_dedx_0 = getMedian(shower_dEdx_0);
  median_dedx_1 = getMedian(shower_dEdx_1);
  median_dedx_2 = getMedian(shower_dEdx_2);
  amalgamated_dedx = getAmalgamateddEdx(angle_0, angle_1, angle_2, median_dedx_0, median_dedx_1, median_dedx_2, shower_dEdx_0.size(), shower_dEdx_1.size(), shower_dEdx_2.size());

  if(matchedToTrackFitter == true){
    dedx0 = track_dedx0;
    dedx1 = track_dedx1;
    dedx2 = track_dedx2;
  }
  if(matchedToTrackFitter == false){
    dedx0 = median_dedx_0;
    dedx1 = median_dedx_1;
    dedx2 = median_dedx_2;
  }
  amalgdedx = amalgamated_dedx;
  
}//getdEdx


////////////// TRACK DEDX //////////////////////////////////////
void SidebandTree::TrackFitdEdx(std::vector<art::Ptr<anab::Calorimetry>> trackcaloobjects, float trunklength, float &dedxU, float &dedxV, float &dedxY){

  dedxU = -1;
  dedxV = -1;
  dedxY = -1;

  float work_function = 23.6e-6;
  float recombination_factor = 0.62;

  for (const auto &tkcalo : trackcaloobjects){

    if (tkcalo->ResidualRange().size() == 0)
      continue;

    auto pl = tkcalo->PlaneID().Plane;
    
    std::vector<float> dedx4cm;
    
    for (size_t ic = 0; ic < tkcalo->ResidualRange().size(); ++ic){
      
      if ((tkcalo->ResidualRange().back() - tkcalo->ResidualRange()[ic]) < trunklength){
  float this_dqdx = tkcalo->dQdx()[ic];
  float this_dedx = -999;
    if(!fIsMC){this_dedx = this_dqdx*data_gain[pl]*work_function*(1/recombination_factor);}
  if(fIsMC){this_dedx = this_dqdx*overlay_gain[pl]*work_function*(1/recombination_factor);}
        dedx4cm.push_back(this_dedx);
      }
    }
    
    float dedx4cm_med = -1.;
    if (dedx4cm.size() > 0)
      {
  std::sort(dedx4cm.begin(), dedx4cm.end());
  if (dedx4cm.size() % 2 == 1)
    dedx4cm_med = dedx4cm[dedx4cm.size() / 2];
  else
    dedx4cm_med = 0.5 * (dedx4cm[dedx4cm.size() / 2] + dedx4cm[dedx4cm.size() / 2 - 1]);
      }

    
    if (pl == 0)
      {
  dedxU = dedx4cm_med;
      }
    if (pl == 1)
      {
  dedxV = dedx4cm_med;
      }
    if (pl == 2)
      {
  dedxY = dedx4cm_med;
      }

  } // for all calorimetry objects associated to the track

  return;
}//TrackFitdEdx


std::vector<double> SidebandTree::CalcdQdxShower(const art::Ptr<recob::Shower>& shower, const art::Ptr<recob::Cluster>& thiscluster, std::vector<art::Ptr<recob::Hit>> hits,  int plane, float trunklength){

  std::vector<double> dqdx;
  
  //get the 3D shower direction
  //note: in previous versions of the code there was a check to fix cases where the shower direction was inverted - this hasn't been implemented
  TVector3 shower_dir(shower->Direction().X(), shower->Direction().Y(),shower->Direction().Z());
  
  //calculate the pitch for this plane
  double pitch = getPitch(shower_dir, plane); 
  
  //keep only clusters on the plane
  if(thiscluster->View() == plane) {
  
    //calculate the cluster direction
    std::vector<double> cluster_axis = {cos(thiscluster->StartAngle()), sin(thiscluster->StartAngle())};    
    
    //get the cluster start and end in CM
    std::cout<<"for plane/tpc/cryo:"<<plane<<"/"<<m_TPC<<"/"<<m_Cryostat<<", fXTicksOffset: "<<theDetector->GetXTicksOffset(plane, m_TPC, m_Cryostat)<<" fXTicksCoefficient: "<<theDetector->GetXTicksCoefficient(m_TPC, m_Cryostat)<<std::endl;
    
    //convert the cluster start and end positions to time and wire coordinates
    std::vector<double> cluster_start = {thiscluster->StartWire() * m_wire_spacing,(thiscluster->StartTick() - theDetector->TriggerOffset())* _time2cm};
    std::vector<double> cluster_end = {thiscluster->EndWire() * m_wire_spacing,(thiscluster->EndTick() - theDetector->TriggerOffset())* _time2cm };
    
    //check that the cluster has non-zero length
    double length = sqrt(pow(cluster_end[0] - cluster_start[0], 2) + pow(cluster_end[1] - cluster_start[1], 2));
    if (length <= 0){ 
      std::cout<<"skipping cluster on plane "<<plane<<", length = "<<length<<std::endl;
      //continue;
    }

    if(length > 0){
    
      //draw a rectangle around the cluster axis 
      std::vector<std::vector<double>> rectangle = buildRectangle(cluster_start, cluster_axis, m_width_dqdx_box, trunklength);  

      //for each hit in the cluster
      for (art::Ptr<recob::Hit> &thishit: hits){  
  //get the hit position in cm from the wire and time
  std::vector<double> thishit_pos = {thishit->WireID().Wire * m_wire_spacing, (thishit->PeakTime() - theDetector->TriggerOffset())* _time2cm};
  
  //check if inside the box
  bool v2 = isInsidev2(thishit_pos, rectangle, trunklength);
  if (v2){
    double q = GetQHit(thishit, plane); 
    double this_dqdx = q/pitch; 
    dqdx.push_back(this_dqdx);
  }//if hit falls inside the box
  
      }//for each hit inthe cluster
    }//cluster has length
  }//cluster on plane
  return dqdx;
  
}




void SidebandTree::getOverlayEnergy(art::Event const & e, art::InputTag allhits_tag, std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> shower_hits, art::InputTag hitparticleassns_tag, int plane, float& overlay_energy){

  art::ServiceHandle<cheat::ParticleInventoryService> particleinventory;
  auto const& allhits_handle = e.getValidHandle<std::vector<recob::Hit>>(allhits_tag);

  overlay_energy = 0;
   
  std::map<int,double> trkID_E; //map of track ID's and their deposited energies
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(allhits_handle,e,hitparticleassns_tag); //mc particles associated with all hits
  std::vector<simb::MCParticle const*> particle_vec; //vector of mc particles 
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec; //vector of associations

  std::vector<art::Ptr<recob::Hit>> overlay_hits;

  //Shower hits loop
  for(size_t i = 0; i < shower_hits.size(); i++){
    art::Ptr<recob::Hit> hit = shower_hits[i];
    auto hitplane = hit->View();
    if(hitplane != plane) continue;

    particle_vec.clear(); match_vec.clear();
    particles_per_hit.get(hit.key(),particle_vec,match_vec); //mc particles matched to each shower hit

    if(particle_vec.size() == 0){
      std::cout<<"FOUND OVERLAY HIT IN SHOWER!"<<std::endl;
      overlay_hits.push_back(hit);
    }

  }//shower hit loop

  overlay_energy = showerEnergy(plane,overlay_hits);

}


float SidebandTree::XOffset(float t){
  auto const &detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  auto const &detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
  double g4Ticks = detClocks->TPCG4Time2Tick(t) + detProperties->GetXTicksOffset(0, 0, 0) - detProperties->TriggerOffset();
  float xoffset = detProperties->ConvertTicksToX(g4Ticks, 0, 0, 0);
  xoffset += 0.6;
  return xoffset;
}

 ///////////////////// CHECK TRUTH //////////////////////////////
void SidebandTree::getTrueParticle(art::Event const & e, art::InputTag allhits_tag, std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> shower_hits, art::InputTag hitparticleassns_tag, int plane, int& shower_mcparticleid, float& shower_mcenergy){

  art::ServiceHandle<cheat::ParticleInventoryService> particleinventory;
  auto const& allhits_handle = e.getValidHandle<std::vector<recob::Hit>>(allhits_tag);

  shower_mcparticleid = -999; //mc particle id that contributes most energy to the shower
  shower_mcenergy = -999;
  
  int mcparticle_id = -999; //mcparticle that contributes most energy to shower
  float max_e = -999.; //max true energy deposition in the shower
  
  //////////////////////////////////////////////////////
  ////////////// SHOWER HITS ///////////////////////////
  //////////////////////////////////////////////////////
  std::map<int,double> trkID_E; //map of track ID's and their deposited energies
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(allhits_handle,e,hitparticleassns_tag); //mc particles associated with all hits
  std::vector<simb::MCParticle const*> particle_vec; //vector of mc particles 
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec; //vector of associations
  
  for(size_t i = 0; i < shower_hits.size(); i++){
    art::Ptr<recob::Hit> hit = shower_hits[i];
    auto hitplane = hit->View();
    if(hitplane != plane) continue;

    particle_vec.clear(); match_vec.clear();
    particles_per_hit.get(hit.key(),particle_vec,match_vec); //mc particles matched to each shower hit

    //ALL MCPARTICLES MATCHED TO CURRENT SHOWER HIT
    for(size_t i_p=0; i_p < particle_vec.size(); ++i_p){

      trkID_E[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id for all shower-matched mc particles
     
      //Particle with max energy contribution to hit
      if( trkID_E[particle_vec[i_p]->TrackId()] > max_e ){ 
  max_e = trkID_E[ particle_vec[i_p]->TrackId() ];
  mcparticle_id = particle_vec[i_p]->TrackId();
      }//max energy particle
      
    }//mcparticles per hit

  }//shower hit loop

  shower_mcparticleid = mcparticle_id;
  shower_mcenergy = max_e;

}
 

///////////////////// CHECK TRUTH //////////////////////////////
void SidebandTree::newTruth(art::Event const & e, art::InputTag allhits_tag, std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> shower_hits, std::vector<art::Ptr<simb::MCParticle>> Gamma1_descendents, std::vector<art::Ptr<simb::MCParticle>> Gamma2_descendents, art::InputTag hitparticleassns_tag, int plane, int& shower_mcparticleid, int& shower_mcparticleid2, float& gamma1_shower, float& gamma2_shower, float& gamma1_event, float& gamma2_event, float& allmc_shower){

  art::ServiceHandle<cheat::ParticleInventoryService> particleinventory;
  auto const& allhits_handle = e.getValidHandle<std::vector<recob::Hit>>(allhits_tag);

  shower_mcparticleid = -999; //mc particle id that contributes most energy to the shower
  shower_mcparticleid2 = -999; //mc particle id that contributes second most energy to the shower

  allmc_shower = -999.;
  gamma1_shower = -999.;
  gamma2_shower = -999.;
  gamma1_event = -999.;
  gamma2_event = -999.;
  
  int mcparticle_id = -999; //mcparticle that contributes most energy to shower
  int mcparticle_id2 = -999; //mcparticle that contributes second most energy to shower
  float max_e = -999.; //max true energy deposition in the shower
  float max_e2 = -999.; //second max true energy deposition in the shower

  float mconly_shower = 0.; //energy deposited by gamma1 in the shower
  float g1_shower = 0.; //energy deposited by gamma1 in the shower
  float g2_shower = 0.; //energy deposited by gamma2 in the shower
  float g1_event = 0.; //energy deposited by gamma1 in the event
  float g2_event = 0.; //energy deposited by gamma2 in the event

  std::vector<art::Ptr<recob::Hit>> overlayHits;
    
  //////////////////////////////////////////////////////
  ////////////// SHOWER HITS ///////////////////////////
  //////////////////////////////////////////////////////
  std::map<int,double> trkID_E; //map of track ID's and their deposited energies
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(allhits_handle,e,hitparticleassns_tag); //mc particles associated with all hits
  std::vector<simb::MCParticle const*> particle_vec; //vector of mc particles 
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec; //vector of associations
  
  for(size_t i = 0; i < shower_hits.size(); i++){
    art::Ptr<recob::Hit> hit = shower_hits[i];
    auto hitplane = hit->View();
    if(hitplane != plane) continue;

    particle_vec.clear(); match_vec.clear();
    particles_per_hit.get(hit.key(),particle_vec,match_vec); //mc particles matched to each shower hit

    if(particle_vec.size() == 0){
      overlayHits.push_back(hit);
    }

    //ALL MCPARTICLES MATCHED TO CURRENT SHOWER HIT
    for(size_t i_p=0; i_p < particle_vec.size(); ++i_p){

      trkID_E[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id for all shower-matched mc particles
     
      //Particle with max energy contribution to hit
      if( trkID_E[particle_vec[i_p]->TrackId()] > max_e ){ 
  max_e = trkID_E[ particle_vec[i_p]->TrackId() ];
  mcparticle_id = particle_vec[i_p]->TrackId();
      }//max energy particle
      
    }//mcparticles per hit

    //loop over mc particles matched to shower hit
    for(size_t i_p=0; i_p < particle_vec.size(); ++i_p){
      
      if(particle_vec[i_p]->TrackId() == mcparticle_id) continue;
   
      //Particle with max energy contribution to hit
      if( trkID_E[particle_vec[i_p]->TrackId()] > max_e2 ){ 
  max_e2 = trkID_E[ particle_vec[i_p]->TrackId() ];
  mcparticle_id2 = particle_vec[i_p]->TrackId();
      }//second max energy particle
      
    }//end loop over particles per hit
  
  }//shower hit loop

  //////////////////////////////////////////////////////
  ////////// ALL HITS IN THE EVENT /////////////////////
  //////////////////////////////////////////////////////
  std::map<int,double> trkID_Eall; //map of track ID's and their deposited energies
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hitall(allhits_handle,e,hitparticleassns_tag); //mc particles associated with all hits
  std::vector<simb::MCParticle const*> particle_vecall; //vector of mc particles 
  std::vector<anab::BackTrackerHitMatchingData const*> match_vecall; //vector of associations
  
  for(size_t i = 0; i < all_hits.size(); i++){
    art::Ptr<recob::Hit> hit = all_hits[i];
    auto hitplane = hit->View();
    if(hitplane != plane) continue;
 
    particle_vecall.clear(); match_vecall.clear();
    particles_per_hitall.get(hit.key(),particle_vecall,match_vecall); //mc particles matched to each event hit

    //ALL MCPARTICLES MATCHED TO CURRENT EVENT HIT
    for(size_t i_p=0; i_p < particle_vecall.size(); ++i_p){

      trkID_Eall[ particle_vecall[i_p]->TrackId() ] += match_vecall[i_p]->energy; //store energy per track id for all matched mc particles
    
    }//mcparticles per hit

  }//all hit loop

  std::cout<<"TURTLE!!! THE SHOWER MATCHES TO "<<trkID_E.size()<<" MC PARTICLES, and THE EVENT MATCHES TO "<<trkID_Eall.size()<<" MC PARTICLES! "<<std::endl;

  ///////////////////////////////////////////////////
  ////////// GET CONTRIBUTIONS //////////////////////
  ///////////////////////////////////////////////////

  /////////////// SHOWER PARTICLE MAP //////////
  for(std::map<int,double>::iterator it = trkID_E.begin(); it!=trkID_E.end(); ++it){

    mconly_shower += it->second;
    
    //Check how much energy gamma1 deposited in the shower
    for(auto const& mcpart : Gamma1_descendents){
      if(mcpart->TrackId() == it->first){
  g1_shower += it->second;
      }
    }//g1

    //Check how much energy gamma2 deposited in the shower
    for(auto const& mcpart : Gamma2_descendents){
      if(mcpart->TrackId() == it->first){
  g2_shower += it->second;
      }
    }//g2
    
  }//trkide

   /////////////// ALL PARTICLE MAP //////////
  for(std::map<int,double>::iterator it = trkID_Eall.begin(); it!=trkID_Eall.end(); ++it){

    //Check how much energy gamma1 deposited in the shower
    for(auto const& mcpart : Gamma1_descendents){
      if(mcpart->TrackId() == it->first){
  g1_event += it->second;
      }
    }//g1

    //Check how much energy gamma2 deposited in the shower
    for(auto const& mcpart : Gamma2_descendents){
      if(mcpart->TrackId() == it->first){
  g2_event += it->second;
      }
    }//g2
    
  }//trkide

  float overlay_shower = showerEnergy(2,overlayHits);

  allmc_shower = overlay_shower + mconly_shower;
  gamma1_shower = g1_shower;
  gamma2_shower = g2_shower;
  gamma1_event = g1_event;
  gamma2_event = g2_event;
  shower_mcparticleid = mcparticle_id;
  shower_mcparticleid2 = mcparticle_id2;
  
}

///////////////////// END SUBRUN //////////////////////////////
void SidebandTree::endSubRun(art::SubRun const & sr) {

  std::cout<<"MOROCCO! enters subrun loop"<<std::endl;
  
  art::Handle <sumdata::POTSummary> potsum_h ;
  if(sr.getByLabel("generator",potsum_h)){

    std::cout<<"ARGENTINA! got subrun from generator"<<std::endl;
    
    if(fIsMC){
      _pot = potsum_h->totpot;
      std::cout<<"NIGER! got pot = "<<_pot<<std::endl;
    }
  }
  
  if(!fIsMC){
    _run = sr.run();
    _subRun = sr.subRun();
  }
  _tree->Fill();
  
}//endsubrun


///////////////////// IN FV FUNCTION //////////////////////////////
bool SidebandTree::inFV(double x, double y, double z){
  // geo::GeometryCore const* fGeometry(lar::providerFrom<geo::Geometry>());
  //   double fDistToEdgeX             = fGeometry->DetHalfWidth()   - 20.;
  //   double fDistToEdgeY             = fGeometry->DetHalfHeight()  - 20.;
  //   double fDistToEdgeZ             = fGeometry->DetLength() / 2. - 10.;
  //   double distInX = x - fGeometry->DetHalfWidth();
  //   double distInY = y;
  //   double distInZ = z - 0.5 * fGeometry->DetLength();
    
  //   if (std::abs(distInX) < fDistToEdgeX && std::abs(distInY) < fDistToEdgeY && std::abs(distInZ) < fDistToEdgeZ)
  //     return true;

  //   return false;

   // if(x >= 5 && x <= 251){
    //   if(y >= -110 && y <= 110){
    //  if((z >= 20 && z < 675) || (z > 775 && z <= 986)){
    //    return true;
    //  }//z
    //   }//y
    // }//x

    if(x >= 10 && x <= 246){
      if(y >= -140 && y <= 140){
  if(z >= 10 && z <= 986.5){
    return true;
  }//z
      }//y
    }//x
  
  return false;
    
}//inFV function

///////////////////// GET PHI //////////////////////////////
double SidebandTree::getPhi(double px, double py, double pz){ //from M. Del Tutto's UBXSec_module.cc

  TVector3 dir(px,py,pz);
  
  // We are in the plane Z = 0 
  dir.SetZ(0);
  TVector3 phi_versor (1, 0, 0);  
  
  double phi = phi_versor.Angle(dir);
  
  // Just convention
  if (dir.Y() < 0)
    phi = -phi;
  
  return phi;
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////   CALORIMETRY FUNCTIONS   ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////// SHOWER RECONSTRUCTED ENERGY //////////////////////////////
float SidebandTree::showerEnergy(int plane, std::vector<art::Ptr<recob::Hit>> hits){

  float energy = 0.;

  //Calibration: converting charge to energy
  //Hit energy = hit integral [ADC] * gain [e-/ADC] * work function [MeV/e] * 1/recombination factor

 
  float work_function = 23.6e-6;
  float recombination_factor = 0.62;

  if(hits.size() == 0){std::cout<<"SHOWER HAS NO HITS!"<<std::endl;}
  
  for(auto const& ihit : hits){
    auto hitplane = ihit->View();
    if(hitplane != plane) continue; //get hits from current plane
    float current_hit_charge = ihit->Integral(); //ADC
    float current_hit_energy = -9999.;

    if(!fIsMC){current_hit_energy = current_hit_charge*data_gain[hitplane]*work_function*(1/recombination_factor);}
    if(fIsMC){current_hit_energy = current_hit_charge*overlay_gain[hitplane]*work_function*(1/recombination_factor);}
    energy += current_hit_energy;
  }//shower hits

  return energy;
  
}

float SidebandTree::showerEnergyDataGain(int plane, std::vector<art::Ptr<recob::Hit>> hits){

  float energy = 0.;

  //Calibration: converting charge to energy
  //Hit energy = hit integral [ADC] * gain [e-/ADC] * work function [MeV/e] * 1/recombination factor

 
  float work_function = 23.6e-6;
  float recombination_factor = 0.62;

  if(hits.size() == 0){std::cout<<"SHOWER HAS NO HITS!"<<std::endl;}

  for(auto const& ihit : hits){
    auto hitplane = ihit->View();
    if(hitplane != plane) continue; //get hits from current plane
    float current_hit_charge = ihit->Integral(); //ADC
    float current_hit_energy = -9999.;

    if(!fIsMC){current_hit_energy = current_hit_charge*data_gain[hitplane]*work_function*(1/recombination_factor);}
    if(fIsMC){current_hit_energy = current_hit_charge*data_gain[hitplane]*work_function*(1/recombination_factor);}
    energy += current_hit_energy;
  }//shower hits

  return energy;
  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

double SidebandTree::QtoEConversion(double Q){
  float work_function = 23.6e-6;
  float recombination_factor = 0.62;
  
  double E = Q*work_function*(1/recombination_factor);
  return E;
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////

double SidebandTree::getAnglewrtWires(TVector3 shower_dir,int plane){
  
  TVector3 wire_dir = getWireVec(plane);
  double cos_theta =  getCoswrtWires(shower_dir, wire_dir);
  
  double theta = acos(cos_theta);
  // return abs(theta);
  return abs(3.14/2 - theta);
  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////

double SidebandTree::getMedian(std::vector<double> thisvector){
  //here the size corresponds to the max index
  int size = thisvector.size() - 1;
  //if no entries, return nonsense value
  if (size <= 0) return NAN;
  
  //find index of median location
  int ind;
  if (size%2 == 0) ind = size/2;
  else ind = size/2 + 1;
  //std::cout<<"the median index in vector with size "<<size+1<<" and  max index "<<size<<" is "<<ind<<std::endl;
  
  double median = thisvector[ind];
  //std::cout<<"returning median value "<< median<<std::endl;
  //return the value at median index
  return median;    
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double SidebandTree::getAmalgamateddEdx(double angle_wrt_plane0, double angle_wrt_plane1, double angle_wrt_plane2, double median_plane0, double median_plane1, double median_plane2, int plane0_nhits, int plane1_nhits, int plane2_nhits){
  
  //if the shower is within 10 degrees of the wires on plane 2, consider planes 1 and 0
  if(angle_wrt_plane2 < degToRad(10)){
    
    //if it's too close to the wires on either of the planes, then stick with plane 2
    if (angle_wrt_plane1> degToRad(20)|| angle_wrt_plane0>degToRad(20) ){
      
      //but if it's outside of the range on plane 1, choose that
      if(angle_wrt_plane1> angle_wrt_plane0){
  return median_plane1;
      } else{
  return median_plane0;
      }
    }
  }
  if (plane2_nhits< 2){
    if (plane1_nhits >=2 ){           
      return median_plane1;
    } else if (plane0_nhits >=2 ){
      return median_plane0;
    }
  }
  
  return median_plane2;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double SidebandTree::degToRad(double deg){
  return (deg*3.14)/180;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double SidebandTree::radToDeg(double rad){
  return (rad*180)/3.14;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double SidebandTree::GetQHit(art::Ptr<recob::Hit> thishitptr, int plane){

  float charge = -999;
  if(!fIsMC){charge = thishitptr->Integral() * data_gain[plane];}
  if(fIsMC){charge = thishitptr->Integral() * overlay_gain[plane];}
  
  return charge;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double SidebandTree::getPitch(TVector3 shower_dir, int plane){
  //get the wire direction for this plane - values are hardcoded which isn't great but the TPC geom object gave weird values
  TVector3 wire_dir = getWireVec(plane);

  //take dot product of shower and wire dir
  double cos = getCoswrtWires(shower_dir, wire_dir);
  
  //want only positive values so take abs, normalize by the lengths of the shower and wire
  cos = abs(cos)/(wire_dir.Mag() * shower_dir.Mag()); 
  
  //If the cos is 0 shower is perpendicular and therefore get infinite distance 
  if (cos == 0){ return std::numeric_limits<double>::max(); }
  
  //output is always >= the wire spacing
  return m_wire_spacing/cos;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TVector3 SidebandTree::getWireVec(int plane){
  TVector3 wire_dir;
  if (plane == 0){
    wire_dir = {0., -sqrt(3) / 2., 1 / 2.};
  } else if (plane == 1){
    wire_dir = {0., sqrt(3) / 2., 1 / 2.};
  } else if (plane == 2) {
    wire_dir = {0., 0., 1.};
  }
  return wire_dir;
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double SidebandTree::getCoswrtWires(TVector3 shower_dir, TVector3 wire_dir){
  //take the dot product between the wire direction and the shower direction
  double cos = wire_dir.Dot(shower_dir);
  return cos;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> SidebandTree::buildRectangle(std::vector<double> cluster_start, std::vector<double> cluster_axis, double width, double length){
  std::vector<std::vector<double>> corners;
  
  //get the axis perpedicular to the cluster axis
  double perp_axis[2] = {-cluster_axis[1], cluster_axis[0]};
  
  //create a vector for each corner of the rectangle on the plane
  //c1 = bottom left corner
  std::vector<double> c1 = {cluster_start[0] + perp_axis[0] * width / 2,  cluster_start[1] + perp_axis[1] * width / 2};
  //c2 = top left corner
  std::vector<double> c2 = {c1[0] + cluster_axis[0] * length, c1[1] + cluster_axis[1] * length};
  //c3 = bottom right corner
  std::vector<double> c3 = {cluster_start[0] - perp_axis[0] * width / 2, cluster_start[1] - perp_axis[1] * width / 2};
  //c4 = top right corner
  std::vector<double> c4 ={c3[0] + cluster_axis[0] * length, c3[1] + cluster_axis[1] * length}; 
  
  //save each of the vectors
  corners.push_back(c1);
  corners.push_back(c2);
  corners.push_back(c4);
  corners.push_back(c3);
  return corners;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//determines if a point is inside the rectangle by summing the areas of the four triangles made by 
//if the point is inside, the sum of the triangles should exactly equal the area of the rectangle
//also returns true if the point is on the boundary
bool SidebandTree::isInsidev2(std::vector<double> thishit_pos, std::vector<std::vector<double >> rectangle, float trunklength){
  int n_vertices = (int)rectangle.size();
  //bool inside = false;
  int i, j = 0;
  double areas = 0;
  //for each pair of vertices
  for (i = 0, j = n_vertices-1; i < n_vertices; j = i++) {
    //calculate the area of a triangle with the point and two vertices
    double this_area = areaTriangle(rectangle[i][0], rectangle[i][1], rectangle[j][0], rectangle[j][1], thishit_pos[0], thishit_pos[1]);
    areas += this_area;
  }        
  //calc area of the rectangle
  double area_rectangle = m_width_dqdx_box*trunklength;
  
  //check the sum of areas match
  if (abs(areas - area_rectangle) <= 0.001 ){
    return true;
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//area of a triangle given three vertices
double SidebandTree::areaTriangle(double x1, double y1, double x2, double y2, double x3, double y3){
  double num = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2);
  return abs(num)/2;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////   TRUTH MATCHING FUNCTIONS   /////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////// GET DESCENDENTS //////////////////////////////
void SidebandTree::getKids(std::vector<art::Ptr<simb::MCParticle>> all_mcpart_vector, int mcpart_id, std::vector<art::Ptr<simb::MCParticle>>& mcpart_kids){

  std::vector< std::vector< TVector2 > > particleV; //vector of vector of (particle id, number of daughter) pairs
  std::vector< TVector2 > mapV; //pairs of (particle id, number of daughter) at each stage
  std::map<int, art::Ptr<simb::MCParticle> > mcpart_id_map; //map from mcparticle to mcparticle id
  
  for(auto const& mcpart: all_mcpart_vector){ 
    mcpart_id_map[mcpart->TrackId()] = mcpart;

    //Step n = 0, get our particle's immediate kids
    if(mcpart->TrackId() == mcpart_id){
      for(int i = 0; i < mcpart->NumberDaughters(); i++){
  for(auto const& mcpart2: all_mcpart_vector){
    int id = mcpart2->TrackId();
    int nd = mcpart2->NumberDaughters();
    if(mcpart2->TrackId() == mcpart->Daughter(i)){
      mcpart_kids.push_back(mcpart2);
      mapV.push_back( TVector2(id,nd) );
    }
  }
      }
    }//our particle
  }//mcparticle loop
  particleV.push_back(mapV);

  //Now get the rest of the descendents
  int n = 1;
  while(n){
    
    if(particleV.at(n-1).size() == 0) break;
    mapV.clear();
    
    for(size_t i = 0; i < particleV.at(n-1).size(); i++){ //loop over (particle id, number of daughter) pairs at previous stage
      int parent_id = particleV.at(n-1).at(i).X(); //mc particle id
      int parent_nd = particleV.at(n-1).at(i).Y(); //number of daughters
      auto const parent_mcpart = mcpart_id_map[parent_id]; 
      
      for(int j = 0; j < parent_nd; j++){ //daughters 
  
  for(auto const & mcpart : all_mcpart_vector){
    int id = mcpart->TrackId();
    int nd = mcpart->NumberDaughters();
    
    if(id == parent_mcpart->Daughter(j)){ 
      mapV.push_back( TVector2(id,nd) );
      mcpart_kids.push_back(mcpart);
      
    }//daughter mcpart
    
  }//mcpart loop
  
      }//ndaughters
      
    }//particle vector
    
    particleV.push_back(mapV);
    
    n++;
    
  }//while
  
}


///////////////////// NEW SHOWER EVALUATION //////////////////////////////
void SidebandTree::newEvaluation(art::Event const & e, art::InputTag allhits_tag, std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> shower_hits, std::vector<art::Ptr<simb::MCParticle>> Gamma1_descendents, std::vector<art::Ptr<simb::MCParticle>> Gamma2_descendents, art::InputTag hitparticleassns_tag, int plane, int& shower_mcparticleid, int& shower_mcparticleid2, float& shower_purity, float& shower_complete, float& reco_e, float& true_e, float& energy_res){

  art::ServiceHandle<cheat::ParticleInventoryService> particleinventory;
  auto const& allhits_handle = e.getValidHandle<std::vector<recob::Hit>>(allhits_tag);
  
  shower_mcparticleid = -999; //mc particle id that contributes most energy to the shower
  shower_mcparticleid2 = -999; //mc particle id that contributes second most energy to the shower
  
  shower_purity = -999.;
  shower_complete = -999.;

  true_e = -999.; //energy deposited by the true gamma in the event
  reco_e = -999.; //energy deposited by the reco shower in the event
  energy_res = -999.; //energy resolution = (true_e - reco_e)/true_e

  int mcparticle_id = -999; //mcparticle that contributes most energy to shower
  int mcparticle_id2 = -999; //mcparticle that contributes second most energy to shower
  float max_e = -999.; //max true energy deposition in the shower
  float max_e2 = -999.; //second max true energy deposition in the shower
  float gamma_e = -999.; //energy deposited by the gamma in the shower
  float g1_depe_shower = 0.; //energy deposited by gamma1 in the shower
  float g2_depe_shower = 0.; //energy deposited by gamma2 in the shower
  float g1_depe_event = 0.; //energy deposited by gamma1 in the event
  float g2_depe_event = 0.; //energy deposited by gamma2 in the event

  std::map<int,double> trkID_E; //map of track ID's and their deposited energies
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(allhits_handle,e,hitparticleassns_tag); //mc particles associated with all hits
  std::vector<simb::MCParticle const*> particle_vec; //vector of mc particles 
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec; //vector of associations

  reco_e = showerEnergy(plane,shower_hits);

  //Shower hits loop
  for(size_t i = 0; i < shower_hits.size(); i++){
    art::Ptr<recob::Hit> hit = shower_hits[i];
    auto hitplane = hit->View();
    if(hitplane != plane) continue;

    particle_vec.clear(); match_vec.clear();
    particles_per_hit.get(hit.key(),particle_vec,match_vec); //mc particles matched to each shower hit

    //loop over mc particles matched to shower hit
    for(size_t i_p=0; i_p < particle_vec.size(); ++i_p){
      
      trkID_E[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id
      //reco_e += match_vec[i_p]->energy; //calculate total energy deposited by shower
    
      //Particle with max energy contribution to hit
      if( trkID_E[particle_vec[i_p]->TrackId()] > max_e ){ 
  max_e = trkID_E[ particle_vec[i_p]->TrackId() ];
  mcparticle_id = particle_vec[i_p]->TrackId();
      }//max energy particle
      
    }//end loop over particles per hit

    //loop over mc particles matched to shower hit
    for(size_t i_p=0; i_p < particle_vec.size(); ++i_p){
      
      if(particle_vec[i_p]->TrackId() == mcparticle_id) continue;
   
      //Particle with max energy contribution to hit
      if( trkID_E[particle_vec[i_p]->TrackId()] > max_e2 ){ 
  max_e2 = trkID_E[ particle_vec[i_p]->TrackId() ];
  mcparticle_id2 = particle_vec[i_p]->TrackId();
      }//second max energy particle
      
    }//end loop over particles per hit
  
  }//shower hit loop

  
  for(std::map<int,double>::iterator it = trkID_E.begin(); it!=trkID_E.end(); ++it){

    //Check how much energy gamma1 deposited in the shower
    for(auto const& mcpart : Gamma1_descendents){
      if(mcpart->TrackId() == it->first){
  g1_depe_shower += it->second;
      }
    }//g1

    //Check how much energy gamma2 deposited in the shower
    for(auto const& mcpart : Gamma2_descendents){
      if(mcpart->TrackId() == it->first){
  g2_depe_shower += it->second;
      }
    }//g2
    
  }//trkide

  std::map<int,double> trkID_E_all;
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit_all(allhits_handle,e,hitparticleassns_tag); //mc particles associated with all hits
  std::vector<simb::MCParticle const*> particle_vec_all; //vector of mc particles 
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec_all; //vector of associations
  
  //All hits loop
  for(size_t i = 0; i < all_hits.size(); i++){
    art::Ptr<recob::Hit> hit = all_hits[i];
   auto hitplane = hit->View();
    if(hitplane != plane) continue;

    particle_vec_all.clear(); match_vec_all.clear();
    particles_per_hit_all.get(hit.key(),particle_vec_all,match_vec_all); //mc particles matched to each hit

    //loop over mc particles matched to hit
    for(size_t i_p=0; i_p < particle_vec_all.size(); ++i_p){
      trkID_E_all[ particle_vec_all[i_p]->TrackId() ] += match_vec_all[i_p]->energy; //store energy per track id
    }//end loop over particles per hit
  }//all hit loop

  for(std::map<int,double>::iterator it = trkID_E_all.begin(); it!=trkID_E_all.end(); ++it){

    //Check how much energy gamma1 deposited in the event
    for(auto const& mcpart : Gamma1_descendents){
      if(mcpart->TrackId() == it->first){
  g1_depe_event += it->second;
      }
    }//g1

    //Check how much energy gamma2 deposited in the event
    for(auto const& mcpart : Gamma2_descendents){
      if(mcpart->TrackId() == it->first){
  g2_depe_event += it->second;
      }
    }//g2
    
  }//trkide_all
  
  ////// Check whether the shower has more energy from g1 or g2
  if(g1_depe_shower == 0 && g2_depe_shower == 0){
    shower_purity = 0.;
    shower_complete = 0.;
    gamma_e = 0.;
  }

  if(g1_depe_shower != 0 || g2_depe_shower != 0){
    if(g1_depe_shower > g2_depe_shower){
      true_e = g1_depe_event;
      gamma_e = g1_depe_shower;
    }
    if(g1_depe_shower <= g2_depe_shower){
      true_e = g2_depe_event;
      gamma_e = g2_depe_shower;
    }
  }
  
  shower_mcparticleid = mcparticle_id;
  shower_mcparticleid2 = mcparticle_id2;
  
   //purity & completeness
  if(reco_e != 0){shower_purity = gamma_e/reco_e;}
  if(true_e != 0){shower_complete = gamma_e/true_e;}
  if(reco_e > 0. && true_e > 0.){energy_res = (true_e - reco_e)/true_e;}
 
}

///////////////////// NEW SHOWER EVALUATION //////////////////////////////
void SidebandTree::getOverlayFraction(art::Event const & e, art::InputTag allhits_tag, std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> shower_hits, art::InputTag hitparticleassns_tag, int plane, float& overlay_fraction){

  art::ServiceHandle<cheat::ParticleInventoryService> particleinventory;
  auto const& allhits_handle = e.getValidHandle<std::vector<recob::Hit>>(allhits_tag);

  float overlay_energy = 0;
  float total_energy = 0;
 
  std::map<int,double> trkID_E; //map of track ID's and their deposited energies
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(allhits_handle,e,hitparticleassns_tag); //mc particles associated with all hits
  std::vector<simb::MCParticle const*> particle_vec; //vector of mc particles 
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec; //vector of associations

  std::vector<art::Ptr<recob::Hit>> overlay_hits;

  //Shower hits loop
  for(size_t i = 0; i < shower_hits.size(); i++){
    art::Ptr<recob::Hit> hit = shower_hits[i];
    auto hitplane = hit->View();
    if(hitplane != plane) continue;

    particle_vec.clear(); match_vec.clear();
    particles_per_hit.get(hit.key(),particle_vec,match_vec); //mc particles matched to each shower hit

    if(particle_vec.size() == 0){
      std::cout<<"FOUND OVERLAY HIT IN SHOWER!"<<std::endl;
      overlay_hits.push_back(hit);
    }

  }//shower hit loop

  total_energy = showerEnergy(plane,shower_hits);
  overlay_energy = showerEnergy(plane,overlay_hits);

  if(total_energy != 0){overlay_fraction = overlay_energy/total_energy;}
}



///////////////////// GET TRUE DEPOSITED ENERGY //////////////////////////////
double SidebandTree::depEnergy(art::Event const & e, art::InputTag allhits_tag, std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<simb::MCParticle>> MCdescendents, art::InputTag hitparticleassns_tag, int plane){

  double deposited_energy = 0.;

  art::ServiceHandle<cheat::ParticleInventoryService> particleinventory;
  auto const& allhits_handle = e.getValidHandle<std::vector<recob::Hit>>(allhits_tag);
  std::map<int,double> trkID_E; //map of track ID's and their deposited energies
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(allhits_handle,e,hitparticleassns_tag); //mc particles associated with all hits
  std::vector<simb::MCParticle const*> particle_vec; //vector of mc particles 
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec; //vector of associations

  //Hit loop
  for(size_t i = 0; i < all_hits.size(); i++){
    
    art::Ptr<recob::Hit> hit = all_hits[i];
    auto hitplane = hit->View();
    if(hitplane != plane) continue;

    particle_vec.clear(); match_vec.clear();
    particles_per_hit.get(hit.key(),particle_vec,match_vec); //mc particles matched to each hit

    //loop over mc particles matched to hit
    for(size_t i_p=0; i_p < particle_vec.size(); ++i_p){
      trkID_E[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id
    }//end loop over particles per hit
  }//all hit loop

  for(std::map<int,double>::iterator it = trkID_E.begin(); it!=trkID_E.end(); ++it){

    //Check how much energy our mcparticle deposited in the event
    for(auto const& mcpart : MCdescendents){
      
      if(mcpart->TrackId() == it->first){
  
  deposited_energy += it->second;
      }
    }

  }//trkide map loop

  return deposited_energy;
  
 
}


///////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// PI0 FUNCTIONS ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


void SidebandTree::altVertexAngles(float VtxX, float VtxY, float VtxZ, art::Ptr<recob::Shower> shower, double& DirX, double& DirY, double& DirZ){

  // DirX = -99999;
  // DirY = -99999;
  // DirZ = -99999;
  auto startx = shower->ShowerStart().X();
  auto starty = shower->ShowerStart().Y();
  auto startz = shower->ShowerStart().Z();

  float vectorX = startx - VtxX;
  float vectorY = starty - VtxY;
  float vectorZ = startz - VtxZ;
  float vectorMagnitude = std::sqrt(std::pow(vectorX,2) + std::pow(vectorY,2) + std::pow(vectorZ,2));
  
  DirX = vectorX/vectorMagnitude;
  DirY = vectorY/vectorMagnitude;
  DirZ = vectorZ/vectorMagnitude;
  
}

void SidebandTree::altTrackAngles(float TrackX, float TrackY, float TrackZ, art::Ptr<recob::Shower> shower, double& DirX, double& DirY, double& DirZ){

  // DirX = -99999;
  // DirY = -99999;
  // DirZ = -99999;
  auto startx = shower->ShowerStart().X();
  auto starty = shower->ShowerStart().Y();
  auto startz = shower->ShowerStart().Z();

  float vectorX = startx - TrackX;
  float vectorY = starty - TrackY;
  float vectorZ = startz - TrackZ;
  float vectorMagnitude = std::sqrt(std::pow(vectorX,2) + std::pow(vectorY,2) + std::pow(vectorZ,2));
  
  DirX = vectorX/vectorMagnitude;
  DirY = vectorY/vectorMagnitude;
  DirZ = vectorZ/vectorMagnitude;
  
}


double SidebandTree::conversionDistance(art::Ptr<recob::Shower> shower, float vx, float vy, float vz){

  float dist3d = -99999.;

  auto startx = shower->ShowerStart().X();
  auto starty = shower->ShowerStart().Y();
  auto startz = shower->ShowerStart().Z();
  dist3d = std::sqrt(std::pow(vx - startx,2) + std::pow(vy - starty,2) + std::pow(vz - startz,2));

  return dist3d;

}


double SidebandTree::radialAngle(art::Ptr<recob::Shower> shower, float vx, float vy, float vz){

  float ang3d = -9999.;
  
  auto startx = shower->ShowerStart().X();
  auto starty = shower->ShowerStart().Y();
  auto startz = shower->ShowerStart().Z();
  auto dirx = shower->Direction().X();
  auto diry = shower->Direction().Y();
  auto dirz = shower->Direction().Z();
  float dist3d = std::sqrt(std::pow(vx - startx,2) + std::pow(vy - starty,2) + std::pow(vz - startz,2));

  ang3d = (((startx - vx)*dirx) + ((starty - vy)*diry) + ((startz - vz)*dirz))/dist3d;

  return ang3d;
}

///////////////////// IMPACT PARAMETER //////////////////////////////
double SidebandTree::impactParameter(art::Ptr<recob::Shower> shower1, art::Ptr<recob::Shower> shower2){

  double ip = -9999999.;
  
  auto startx1 = shower1->ShowerStart().X();
  auto starty1 = shower1->ShowerStart().Y();
  auto startz1 = shower1->ShowerStart().Z();
  auto dirx1 = shower1->Direction().X();
  auto diry1 = shower1->Direction().Y();
  auto dirz1 = shower1->Direction().Z();

  auto startx2 = shower2->ShowerStart().X();
  auto starty2 = shower2->ShowerStart().Y();
  auto startz2 = shower2->ShowerStart().Z();
  auto dirx2 = shower2->Direction().X();
  auto diry2 = shower2->Direction().Y();
  auto dirz2 = shower2->Direction().Z();

  double numerator = (diry1*dirz2 - dirz1*diry2)*(startx1 - startx2) + (dirz1*dirx2 - dirx1*dirz2)*(starty1 - starty2) + (dirx1*diry2 - diry1*dirx2)*(startz1 - startz2);
  double denominator = std::sqrt(std::pow(diry1*dirz2 - dirz1*diry2,2) + std::pow(dirz1*dirx2 - dirx1*dirz2,2) + std::pow(dirx1*diry2 - diry1*dirx2,2));

  if(denominator != 0){ip = numerator/denominator;}

  return ip;
}

///////////////////// IMPACT PARAMETER //////////////////////////////
double SidebandTree::impactParameterWithTrack(art::Ptr<recob::Track> track, art::Ptr<recob::Shower> shower, bool flipped){

  double ip = -9999999.;

  auto startx1 = -999; auto starty1 = -999; auto startz1 = -999; auto dirx1 = -999; auto diry1 = -999; auto dirz1 = -999;
  if(flipped == false){
    startx1 = track->Vertex().X();
    starty1 = track->Vertex().Y();
    startz1 = track->Vertex().Z();
    dirx1 = track->VertexDirection().X();
    diry1 = track->VertexDirection().Y();
    dirz1 = track->VertexDirection().Z();
  }
  if(flipped == true){
    startx1 = track->End().X();
    starty1 = track->End().Y();
    startz1 = track->End().Z();
    dirx1 = track->EndDirection().X();
    diry1 = track->EndDirection().Y();
    dirz1 = track->EndDirection().Z();
  }

  auto startx2 = shower->ShowerStart().X();
  auto starty2 = shower->ShowerStart().Y();
  auto startz2 = shower->ShowerStart().Z();
  auto dirx2 = shower->Direction().X();
  auto diry2 = shower->Direction().Y();
  auto dirz2 = shower->Direction().Z();

  double numerator = (diry1*dirz2 - dirz1*diry2)*(startx1 - startx2) + (dirz1*dirx2 - dirx1*dirz2)*(starty1 - starty2) + (dirx1*diry2 - diry1*dirx2)*(startz1 - startz2);
  double denominator = std::sqrt(std::pow(diry1*dirz2 - dirz1*diry2,2) + std::pow(dirz1*dirx2 - dirx1*dirz2,2) + std::pow(dirx1*diry2 - diry1*dirx2,2));

  if(denominator != 0){ip = numerator/denominator;}

  return ip;
}

///////////////////// ANGLE BETWEEN TWO SHOWERS //////////////////////////////
double SidebandTree::angleBetweenTwoShowers(art::Ptr<recob::Shower> shower1, art::Ptr<recob::Shower> shower2){

  //double angle = -9999999.;
  double cosangle = -9999999.;
  
  auto dirx1 = shower1->Direction().X();
  auto diry1 = shower1->Direction().Y();
  auto dirz1 = shower1->Direction().Z();

  auto dirx2 = shower2->Direction().X();
  auto diry2 = shower2->Direction().Y();
  auto dirz2 = shower2->Direction().Z();

  cosangle = (dirx1*dirx2) + (diry1*diry2) + (dirz1*dirz2);
  //angle = std::acos(cosangle);

  return cosangle;
}

///////////////////// ANGLE BETWEEN TWO SHOWERS //////////////////////////////
double SidebandTree::angleBetweenTrackShower(art::Ptr<recob::Track> track, art::Ptr<recob::Shower> shower, bool flipped){

  //double angle = -9999999.;
  double cosangle = -9999999.;

  auto dirx1 = -999; auto diry1 = -999; auto dirz1 = -999;
  if(flipped == false){
    dirx1 = track->VertexDirection().X();
    diry1 = track->VertexDirection().Y();
    dirz1 = track->VertexDirection().Z();
  }
  if(flipped == true){
    dirx1 = track->EndDirection().X();
    diry1 = track->EndDirection().Y();
    dirz1 = track->EndDirection().Z();
  }
  
  auto dirx2 = shower->Direction().X();
  auto diry2 = shower->Direction().Y();
  auto dirz2 = shower->Direction().Z();

  cosangle = (dirx1*dirx2) + (diry1*diry2) + (dirz1*dirz2);
  //angle = std::acos(cosangle);

  return cosangle;
}

///////////////////// CONVERSION LENGTH //////////////////////////////
double SidebandTree::conversionLength(float vx, float vy, float vz, float x, float y, float z){

  double cl = std::sqrt(std::pow(vx - x,2) + std::pow(vy - y,2) + std::pow(vz - z,2));

  return cl;
}

////////////////////////// PI0 MASS //////////////////////////////////
double SidebandTree::pi0Mass(double energy1, double energy2, double angle){

  double mass = -9999;

  //double cosine_angle = std::cos(angle);
  //mass = std::sqrt(2*energy1*energy2*(1 - cosine_angle));

  mass = std::sqrt(2*energy1*energy2*(1 - angle));

  return mass;
 
}

////////////////////////// PI0 MOMENTUM //////////////////////////////////
double SidebandTree::pi0Momentum(double mass, double energy1, double energy2, double angle){

  double momentum = -9999;

  //double cosine_angle = std::cos(angle);
  double alpha = (std::abs(energy1 - energy2))/(energy1 + energy2);

  //momentum = mass*std::sqrt(2/((1 - std::pow(alpha,2))*(1 - cosine_angle)));
  momentum = mass*std::sqrt(2/((1 - std::pow(alpha,2))*(1 - angle)));

  return momentum;
 
}

////////////////////////// ALTERNATE PI0 MOMENTUM //////////////////////////////////
double SidebandTree::Altpi0Momentum(double energy1, double energy2, double angle){

  double momentum = -9999;

  momentum = energy1 + (energy2*angle) + ((std::pow(energy2,2)*std::pow(1 - std::pow(angle,2),2))/(energy1 + (energy2*angle)));
  return momentum;
 
}

////////////////////////// PI0 ANGLE //////////////////////////////////
double SidebandTree::pi0Angle(double p, double energy1, double energy2, art::Ptr<recob::Shower> gamma1, art::Ptr<recob::Shower> gamma2, double angle){

  auto g1cosx = gamma1->Direction().X();
  auto g1cosy = gamma1->Direction().Y();
  auto g1cosz = gamma1->Direction().Z();

  auto g2cosx = gamma2->Direction().X();
  auto g2cosy = gamma2->Direction().Y();
  auto g2cosz = gamma2->Direction().Z();

  double theta1 = (energy2*std::sqrt(1 - std::pow(angle,2)))/(energy1 + (energy2*angle));
  double theta2 = std::acos(angle) - theta1;

  double a = energy1*g1cosx;
  double b = energy1*g1cosy;
  double c = energy1*g1cosz;

  double d = energy2*g2cosx;
  double e = energy2*g2cosy;
  double f = energy2*g2cosz;

  double g = energy1*p*std::cos(theta1);
  double h = energy2*p*std::cos(theta2);

  //double px = -(a*b*e*h + a*c*f*h - a*(std::pow(e,2))*g - a*(std::pow(f,2))*g + (std::pow(b,2))*(-d)*h + b*d*e*g - (std::pow(c,2))*d*h + c*d*f*g)/((std::pow(a,2))*(std::pow(e,2)) + (std::pow(a,2))*(std::pow(f,2)) - 2*a*b*d*e - 2*a*c*d*f + (std::pow(b,2))*(std::pow(d,2)) + (std::pow(b,2))*(std::pow(f,2)) - 2*b*c*e*f + (std::pow(c,2))*(std::pow(d,2)) + (std::pow(c,2))*(std::pow(e,2)));
  
  //double py = -((std::pow(a,2))*(-e)*h + a*b*d*h + a*d*e*g + b*c*f*h - b*(std::pow(d,2))*g - b*(std::pow(f,2))*g - (std::pow(c,2))*e*h + c*e*f*g)/((std::pow(a,2))*(std::pow(e,2)) + (std::pow(a,2))*(std::pow(f,2)) - 2*a*b*d*e - 2*a*c*d*f + (std::pow(b,2))*(std::pow(d,2)) + (std::pow(b,2))*(std::pow(f,2)) - 2*b*c*e*f + (std::pow(c,2))*(std::pow(d,2)) + (std::pow(c,2))*(std::pow(e,2)));
  
  double pz = -((std::pow(a,2))*(-f)*h + a*c*d*h + a*d*f*g - (std::pow(b,2))*f*h + b*c*e*h + b*e*f*g - c*(std::pow(d,2))*g - c*(std::pow(e,2))*g)/((std::pow(a,2))*(std::pow(e,2)) + (std::pow(a,2))*(std::pow(f,2)) - 2*a*b*d*e - 2*a*c*d*f + (std::pow(b,2))*(std::pow(d,2)) + (std::pow(b,2))*(std::pow(f,2)) - 2*b*c*e*f + (std::pow(c,2))*(std::pow(d,2)) + (std::pow(c,2))*(std::pow(e,2)));
  
  double pcosz = pz/p;

  return pcosz;
 
}

bool SidebandTree::passesShowerCuts(double dist3d, double energy2, double ang3d){
     return (((dist3d > minConversionDist && dist3d < maxConversionDist) || (energy2 > minYPlaneEnergy && energy2 < maxYPlaneEnergy && dist3d < maxConversionDist)) && (energy2 > minYPlaneEnergy) && (ang3d > minRadialAngle) );

  }



/////////////////////////////// END JOB ////////////////////////////
void SidebandTree::endJob(){


  mf::LogVerbatim("SidebandTree")<< "SidebandTree finished job";
  std::cout<<"THIS OVERLAY SAMPLE HAS "<<NALL<<" EVENTS AND "<<NNEUTRINO<<" NEUTRINO EVENTS"<<std::endl;
  std::cout<<"THIS EVENT HAS "<<ncontained<<" CONTAINED MUONS AND "<<nexiting<<" EXITING MUONS"<<std::endl;
  std::cout<<"MATT EVENTS HAVE "<<nchargedpi<<" CHARGED PI, "<<nproton<<" PROTON, "<<nneutron<<" NEUTRON, "<<nmuon<<" MUON, "<<nneutralpi<<" PI0, "<<nelectron<<" ELECTRON, "<<nother<<" OTHER"<<std::endl;
                                                 
                                                 
}//End job
