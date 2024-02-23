// -*- C++ -*-
//
// Package:    Analyser/ClusterAnalyser
// Class:      ClusterAnalyser
//
/**\class ClusterAnalyser ClusterAnalyser.cc Analyser/ClusterAnalyser/plugins/ClusterAnalyser.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Tue, 20 Feb 2024 12:56:38 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ForwardDetId/interface/HFNoseDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "TTree.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class ClusterAnalyser : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ClusterAnalyser(const edm::ParameterSet&);
  ~ClusterAnalyser();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  //  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file

  edm::EDGetTokenT<std::vector<reco::CaloCluster>>           layerClusCollection_;

  edm::ESHandle<CaloGeometry>          pG_;

  edm::EDGetTokenT<HGCRecHitCollection> hgcEEHitCollection_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcHEFHitCollection_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcHEBHitCollection_;


  ///Tree variables
  ////Tree variables
  /////EE
  Int_t nheeRH_;
  std::vector<Float_t> heerhE_;
  std::vector<Float_t> heerhEta_;
  std::vector<Float_t> heerhPhi_;
  std::vector<Float_t> heerhX_;
  std::vector<Float_t> heerhY_;
  std::vector<Float_t> heerhZ_;

  Int_t nhhefRH_;
  std::vector<Float_t> hhefrhE_;
  std::vector<Float_t> hhefrhEta_;
  std::vector<Float_t> hhefrhPhi_;
  std::vector<Float_t> hhefrhX_;
  std::vector<Float_t> hhefrhY_;
  std::vector<Float_t> hhefrhZ_;

  Int_t nhhebRH_;
  std::vector<Float_t> hhebrhE_;
  std::vector<Float_t> hhebrhEta_;
  std::vector<Float_t> hhebrhPhi_;
  std::vector<Float_t> hhebrhX_;
  std::vector<Float_t> hhebrhY_;
  std::vector<Float_t> hhebrhZ_;
  
  Int_t nhclus_;
  std::vector<Float_t> hclusE_;
  std::vector<Float_t> hclusEta_;
  std::vector<Float_t> hclusPhi_;
  std::vector<Float_t> hclusX_;
  std::vector<Float_t> hclusY_;
  std::vector<Float_t> hclusZ_;

  std::vector<Int_t> nrhclus_;
  std::vector<std::vector<Float_t>> hclusrhE_;
  std::vector<std::vector<Float_t>> hclusrhEta_;
  std::vector<std::vector<Float_t>> hclusrhPhi_;
  std::vector<std::vector<Float_t>> hclusrhX_;
  std::vector<std::vector<Float_t>> hclusrhY_;
  std::vector<std::vector<Float_t>> hclusrhZ_;

  TTree   *tree;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ClusterAnalyser::ClusterAnalyser(const edm::ParameterSet& iConfig){
  layerClusCollection_   = consumes<std::vector<reco::CaloCluster>>(iConfig.getParameter<edm::InputTag>("layerClusCollection"));

  hgcEEHitCollection_  = consumes<HGCRecHitCollection> (iConfig.getParameter<edm::InputTag>("hgcEEHitCollection"));
  hgcHEFHitCollection_  = consumes<HGCRecHitCollection> (iConfig.getParameter<edm::InputTag>("hgcHEFHitCollection"));
  hgcHEBHitCollection_  = consumes<HGCRecHitCollection> (iConfig.getParameter<edm::InputTag>("hgcHEBHitCollection"));


  edm::Service<TFileService> fs;
  tree    = fs->make<TTree>("EventTree", "Event data");
  
  ///EE
  tree->Branch("nheeRH",         &nheeRH_);
  tree->Branch("heerhE",         &heerhE_);
  tree->Branch("heerhEta",         &heerhEta_);
  tree->Branch("heerhPhi",         &heerhPhi_);
  tree->Branch("heerhX",         &heerhX_);
  tree->Branch("heerhY",         &heerhY_);
  tree->Branch("heerhZ",         &heerhZ_);

  tree->Branch("nhhebRH",         &nhhebRH_);
  tree->Branch("hhebrhE",         &hhebrhE_);
  tree->Branch("hhebrhEta",         &hhebrhEta_);
  tree->Branch("hhebrhPhi",         &hhebrhPhi_);
  tree->Branch("hhebrhX",         &hhebrhX_);
  tree->Branch("hhebrhY",         &hhebrhY_);
  tree->Branch("hhebrhZ",         &hhebrhZ_);

  tree->Branch("nhhefRH",         &nhhefRH_);
  tree->Branch("hhefrhE",         &hhefrhE_);
  tree->Branch("hhefrhEta",         &hhefrhEta_);
  tree->Branch("hhefrhPhi",         &hhefrhPhi_);
  tree->Branch("hhefrhX",         &hhefrhX_);
  tree->Branch("hhefrhY",         &hhefrhY_);
  tree->Branch("hhefrhZ",         &hhefrhZ_);

  tree->Branch("nhclus",         &nhclus_);
  tree->Branch("hclusE",         &hclusE_);
  tree->Branch("hclusEta",         &hclusEta_);
  tree->Branch("hclusPhi",         &hclusPhi_);
  tree->Branch("hclusX",         &hclusX_);
  tree->Branch("hclusY",         &hclusY_);
  tree->Branch("hclusZ",         &hclusZ_);

  tree->Branch("nrhclus",         &nrhclus_);
  tree->Branch("hclusrhE",         &hclusrhE_);
  tree->Branch("hclusrhEta",         &hclusrhEta_);
  tree->Branch("hclusrhPhi",         &hclusrhPhi_);
  tree->Branch("hclusrhX",         &hclusrhX_);
  tree->Branch("hclusrhY",         &hclusrhY_);
  tree->Branch("hclusrhZ",         &hclusrhZ_);


}

ClusterAnalyser::~ClusterAnalyser() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void ClusterAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;


  edm::Handle<HGCRecHitCollection> hgcEEHitsHandle;
  edm::Handle<HGCRecHitCollection> hgcHEFHitsHandle;
  edm::Handle<HGCRecHitCollection> hgcHEBHitsHandle;
  iEvent.getByToken(hgcEEHitCollection_,hgcEEHitsHandle);
  iEvent.getByToken(hgcHEFHitCollection_,hgcHEFHitsHandle);
  iEvent.getByToken(hgcHEBHitCollection_,hgcHEBHitsHandle);



  edm::Handle<std::vector<reco::CaloCluster>> layerCluster;
  iEvent.getByToken(layerClusCollection_, layerCluster);


  iSetup.get<CaloGeometryRecord>().get(pG_);
  const CaloGeometry* geo = pG_.product();


  ////Tree variables
  /////EE
  nheeRH_ = 0;
  heerhE_.clear();
  heerhEta_.clear();
  heerhPhi_.clear();
  heerhX_.clear();
  heerhY_.clear();
  heerhZ_.clear();

  nhhefRH_ = 0;
  hhefrhE_.clear();
  hhefrhEta_.clear();
  hhefrhPhi_.clear();
  hhefrhX_.clear();
  hhefrhY_.clear();
  hhefrhZ_.clear();

  nhhebRH_ = 0;
  hhebrhE_.clear();
  hhebrhEta_.clear();
  hhebrhPhi_.clear();
  hhebrhX_.clear();
  hhebrhY_.clear();
  hhebrhZ_.clear();
  
  nhclus_ = 0;
  hclusE_.clear();
  hclusEta_.clear();
  hclusPhi_.clear();
  hclusX_.clear();
  hclusY_.clear();
  hclusZ_.clear();

  hclusrhE_.clear();
  hclusrhEta_.clear();
  hclusrhPhi_.clear();
  hclusrhX_.clear();
  hclusrhY_.clear();
  hclusrhZ_.clear();

  
  for (std::vector<reco::CaloCluster>::const_iterator iclus = layerCluster->begin(); iclus != layerCluster->end(); ++iclus){
    double en = iclus->energy();
    double x = iclus->x();
    double y = iclus->y();
    double z = iclus->z();
    double eta = iclus->eta();
    double phi = iclus->phi();

    nhclus_++;
    //std::cout<<" en : eta : phi "<<en<<" "<<eta<<" "<<phi<<std::endl;

    hclusE_.push_back(en);
    hclusX_.push_back(x);
    hclusY_.push_back(y);
    hclusZ_.push_back(z);
    hclusEta_.push_back(eta);
    hclusPhi_.push_back(phi);

    ///hits and fractions
    std::vector<Float_t> rhen;
    std::vector<Float_t> rheta;
    std::vector<Float_t> rhphi;
    std::vector<Float_t> rhx;
    std::vector<Float_t> rhy;
    std::vector<Float_t> rhz;
    nrhclus_.clear();
    double tmphrh = 0;
    //all rechits and the energy fractions in this cluster
    //// https://cmssdt.cern.ch/lxr/source/Validation/RecoParticleFlow/plugins/PFAnalysisNtuplizer.cc#1137
    const auto& rechit_fracs = iclus->hitsAndFractions();
    for (const auto& rh : rechit_fracs) {
      /*      if (detids.find(rh.first.rawId()) != detids.end()) {
	continue;
      }
      */
      //detids[rh.first.rawId()] += cluster.energy() * rh.second;
      const auto id = rh.first;
      const GlobalPoint & rechitPoint = geo->getPosition(id);

      std::cout<<""<<std::endl;
      std::cout<<rh.second<<" "<<rechitPoint.x()<<" "<<rechitPoint.y()<<" "<<rechitPoint.z()<<std::endl;

      rhen.push_back(rh.second * en);
      rhx.push_back(rechitPoint.x());
      rhy.push_back(rechitPoint.y());
      rhz.push_back(rechitPoint.z());
      rheta.push_back(rechitPoint.eta());
      rhphi.push_back(rechitPoint.phi());

      /*hbherhEta_.push_back(rechitPoint.eta());
      hbherhPhi_.push_back(rechitPoint.phi());
      hbherhX_.push_back(rechitPoint.x());
      hbherhY_.push_back(rechitPoint.y());
      hbherhZ_.push_back(rechitPoint.z());
      */
      tmphrh++;

    }///for (const auto& rh : rechit_fracs)
    nrhclus_.push_back(tmphrh);
    hclusrhE_.push_back(rhen);
    hclusrhX_.push_back(rhx);
    hclusrhY_.push_back(rhy);
    hclusrhZ_.push_back(rhz);
    hclusrhEta_.push_back(rheta);
    hclusrhPhi_.push_back(rhphi);

  }///for (std::vector<reco::CaloCluster>::const_iterator iclus = layerCluster->

  ////Loop on rechit collection
  const HGCRecHitCollection* hgcEERecHits = nullptr;
  const HGCRecHitCollection* hgcHEBRecHits = nullptr;
  const HGCRecHitCollection* hgcHEFRecHits = nullptr;
  
  if ( !hgcEEHitsHandle.isValid() ){
    LogDebug("") << "Error! EE rechits can't get product!" << std::endl;
  } else{
    hgcEERecHits = hgcEEHitsHandle.product();
  }

  if ( !hgcHEFHitsHandle.isValid() ){
    LogDebug("") << "Error! HEB rechits can't get product!" << std::endl;
  } else{
    hgcHEBRecHits = hgcHEFHitsHandle.product();
  }
  
  if ( !hgcHEBHitsHandle.isValid() ){
    LogDebug("") << "Error! HEF rechits can't get product!" << std::endl;
  } else{
    hgcHEFRecHits = hgcHEBHitsHandle.product();
  }

  
  //EE
  for( HGCRecHitCollection::const_iterator heerechit = hgcEERecHits->begin(); heerechit != hgcEERecHits->end(); heerechit++ ){

    double Energy = heerechit->energy();

    const auto id = heerechit->id();
    const GlobalPoint & rechitPoint = geo->getPosition(id);
    
    double x = rechitPoint.x();
    double y = rechitPoint.y();
    double z = rechitPoint.z();
    double eta = rechitPoint.eta();
    double phi = rechitPoint.phi();

    heerhE_.push_back(Energy);
    heerhX_.push_back(x);
    heerhY_.push_back(y);
    heerhZ_.push_back(z);
    heerhEta_.push_back(eta);
    heerhPhi_.push_back(phi);

    nheeRH_++;
  }/// Loop over EE


  //HEB
  for( HGCRecHitCollection::const_iterator hhebrechit = hgcHEBRecHits->begin(); hhebrechit != hgcHEBRecHits->end(); hhebrechit++ ){

    double Energy = hhebrechit->energy();

    const auto id = hhebrechit->id();
    const GlobalPoint & rechitPoint = geo->getPosition(id);
    
    double x = rechitPoint.x();
    double y = rechitPoint.y();
    double z = rechitPoint.z();
    double eta = rechitPoint.eta();
    double phi = rechitPoint.phi();

    hhebrhE_.push_back(Energy);
    hhebrhX_.push_back(x);
    hhebrhY_.push_back(y);
    hhebrhZ_.push_back(z);
    hhebrhEta_.push_back(eta);
    hhebrhPhi_.push_back(phi);
    
    nhhebRH_++;
  }/// Loop over HEB


  //HEF
  for( HGCRecHitCollection::const_iterator hhefrechit = hgcHEFRecHits->begin(); hhefrechit != hgcHEFRecHits->end(); hhefrechit++ ){

    double Energy = hhefrechit->energy();

    const auto id = hhefrechit->id();
    const GlobalPoint & rechitPoint = geo->getPosition(id);
    
    double x = rechitPoint.x();
    double y = rechitPoint.y();
    double z = rechitPoint.z();
    double eta = rechitPoint.eta();
    double phi = rechitPoint.phi();

    hhefrhE_.push_back(Energy);
    hhefrhX_.push_back(x);
    hhefrhY_.push_back(y);
    hhefrhZ_.push_back(z);
    hhefrhEta_.push_back(eta);
    hhefrhPhi_.push_back(phi);
    
    nhhefRH_++;
  }/// Loop over HEF

  tree->Fill();
  
  

}

// ------------ method called once each job just before starting event loop  ------------
void ClusterAnalyser::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void ClusterAnalyser::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ClusterAnalyser::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ClusterAnalyser);
