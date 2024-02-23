import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
process = cms.Process('Analyse',Phase2C9)

process.load("FWCore.MessageService.MessageLogger_cfi")
#######################
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#############################

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )


process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
                                'file:/home/shilpi/work/CMSSW_12_1_0_pre2/src/23293.0_CloseByParticleGun+2026D49+CloseByParticle_Photon_ERZRanges_GenSimHLBeamSpotHGCALCloseBy+DigiTrigger+RecoGlobal+HARVESTGlobal/step3.root' 
                    )
                            
                            #secondaryFileNames = cms.untracked.vstring(
                            #)
                    
                            
                    )


process.InterimOutput = cms.OutputModule("PoolOutputModule",
                                         fileName = cms.untracked.string('myOutputFile.root'), 
                                          SelectEvents = cms.untracked.PSet(
                                              SelectEvents = cms.vstring("p")
                                              ),
                                         outputCommands = cms.untracked.vstring('keep *')
                                         
                                     )

                                                        



process.TFileService = cms.Service("TFileService", fileName = cms.string('hgc.root'))

process.analyse = cms.EDAnalyzer('ClusterAnalyser',
                                 layerClusCollection = cms.InputTag("hgcalLayerClusters","", "RECO"),
                                 hgcEEHitCollection  = cms.InputTag("HGCalRecHit","HGCEERecHits", "RECO"),  
                                 hgcHEFHitCollection = cms.InputTag("HGCalRecHit","HGCHEFRecHits", "RECO"),  
                                 hgcHEBHitCollection = cms.InputTag("HGCalRecHit","HGCHEBRecHits", "RECO")
)


process.p = cms.Path(
    process.analyse
)

#process.e = cms.EndPath(process.InterimOutput)  


#print(process.dumpPython())
