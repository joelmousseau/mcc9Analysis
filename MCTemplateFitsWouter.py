import ROOT, math, sys, os
import uproot
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import csv
from mpl_toolkits.mplot3d import Axes3D
from ROOT import TH1, TAxis, gROOT, TCanvas
####################################################################################################
def chanToHistogram(channel):
    if channel == "QE":
       return 0
    elif channel == "RES":
       return 1
    elif channel == "DIS":
       return 2
    elif channel == "MEC":
       return 3
    elif channel == "Other":
       return 4
    elif channel == "OFFV":
       return 5
    else:
       return -1

def dataify(array, bins, limits):
   counts, bin_edges = np.histogram(array, bins, range=limits)
   bin_centers       = (bin_edges[:-1] + bin_edges[1:]) / 2
   errs = np.sqrt(counts)
   return [bin_centers, counts, errs]
   
def getChan(interaction, isNC):
    if (isNC):
      return "NC / Other"
    elif(interaction == 0):
      return "QE"
    elif(interaction == 1):
      return "RES"
    elif(interaction == 2):
      return "DIS"
    elif(interaction == 10):
      return "2p2h"

def getPhi(pY, pX):
    return (np.arctan2(pY, pX) / np.pi)

def getQ2(nuE, nuMu, thetaMu):
    return 4*nuE*nuMu*math.pow(math.sin(thetaMu/2), 2)

def getW2(Ehad, Q2):
    targetMass = 0.989
    return 2*targetMass*Ehad + math.pow(targetMass, 2) - Q2

def getW(Ehad, Q2):
    return math.sqrt(getW2(Ehad, Q2))

def getXbj(Ehad, Q2):
    targetMass = 0.989
    if(Ehad > 0.0):
      return Q2 / (2*targetMass*Ehad)
    else:
      return -1.0

def getInel(Ehad, Enu):
    if(Enu > 0.0):
      return (Ehad / Enu)
    else:
      return -1.0

class TemplateFitter( ROOT.TPyMultiGenFunction ):
    def __init__( self, mcDict, extVector, dataVector) :
        print "CREATED"
        ROOT.TPyMultiGenFunction.__init__( self, self )
        self.MC_dict = mcDict
        self.Data    = dataVector
        self.Ext     = extVector
        self.ScaleFactors = []

    def setFitRange( self, xlo, xhi ) :
        self.xlo = xlo
        self.xhi = xhi
    
    def NDim( self ):
        print 'PYTHON NDim called'
        return 2

    def DoEval( self, args ) :
        #mcTemplatesP     = {"QE" : empty, "RES" : empty, "DIS" : empty, "MEC" : empty, "Other" : empty, "OFFV" : empty}
        scale = 1.0
        mcVector = [0.0]*len(self.Data)
        for itype in self.MC_dict:
            
            if itype   == 'QE'    :   scale = scale*args[0]
        
            elif itype == 'RES'   :   scale = scale*args[0]
        
            elif itype == 'DIS'   :   scale = scale*args[0]
        
            elif itype == 'MEC'   :   scale = scale*args[0]
        
            elif itype == 'Other' :   scale = scale*args[0]
            
            elif itype == 'OFFV'  :   scale = scale*args[2]
            
            else: break
            #print len(self.Data)
            #print len(self.Ext)
            #print range(len(self.MC_dict[itype]))
            for bin in range(len(self.Data)):
              mcVector[bin] += scale*self.MC_dict[itype][bin]
        
        
        chi2 = 0.0
        
        for bin in range(len(self.Data)):
            #if bin < self.xlo or bin > self.xhi: continue
           error = math.sqrt(self.Data[bin])
           data = self.Data[bin] - args[1]*self.Ext[bin]
           mc   = mcVector[bin]
           if not error : error = 1.0
           chi2 += math.pow( ( ( mc - data ) / error ), 2 )

        return chi2

#print(sys.argv[1])
#study_dir = "Track_Plots/"
InputFiles = ["/uboone/app/users/wvdp/RootTrees/v20/run1/nucc_nu_overlay_run1_mcc9.root", "/uboone/app/users/wvdp/RootTrees/v20/run1/nucc_on_data_run1_mcc9.root", "/uboone/app/users/wvdp/RootTrees/v20/run1/nucc_off_data_run1_mcc9.root"]

#OverlayScale  = 1.0
ExtScale     = 0.97
numMCTemplates = 6
empty = []

maxProtonChi2 = 88.0
minTrackScore = 0.5

#plt.figure()

overlayEvents = uproot.open(InputFiles[0])["NuCCanalyzer"]["Event"]
bnbEvents     = uproot.open(InputFiles[1])["NuCCanalyzer"]["Event"]
extEvents     = uproot.open(InputFiles[2])["NuCCanalyzer"]["Event"]

overlayPOT    = uproot.open(InputFiles[0])["NuCCanalyzer"]["subruns"]

mcPOT         = pd.Series(overlayPOT.array("pot")).sum()
dataPOT       = 4.206e+19
bnbSpills     = 9932159.0
extTriggers   = 14675888.0

OverlayScale = dataPOT / mcPOT
ExtScale     = bnbSpills / extTriggers

print "MC POT: %e Overlay Scale: %.3f Ext Scale: %.3f" % (mcPOT, OverlayScale, ExtScale)

overlayDaughters = uproot.open(InputFiles[0])["NuCCanalyzer"]["Daughters"]
trackDaughters   = pd.DataFrame(overlayDaughters.arrays(["track_range_mom_mu", "track_mcs_mom", "track_range_mom_p", "track_is_muon_candidate", "track_score", "track_chi2_proton", "track_chi2_muon", "track_dirx", "track_diry", "track_dirz", "run", "subrun", "event"] ) )
filteredEvents   = pd.DataFrame(overlayEvents.arrays(["run", "subrun", "event", "mc_nu_interaction_type", "mc_nu_ccnc", "nu_mu_cc_selected", "mc_nu_lepton_energy", "mc_nu_energy", "mc_nu_lepton_theta"]) )

dataDaughters = uproot.open(InputFiles[1])["NuCCanalyzer"]["Daughters"]
trackData     = pd.DataFrame(dataDaughters.arrays(["track_range_mom_mu", "track_mcs_mom", "track_range_mom_p", "track_is_muon_candidate", "track_score", "track_chi2_proton", "track_chi2_muon", "track_dirx", "track_diry", "track_dirz", "run", "subrun", "event"] ) )
filteredData  = pd.DataFrame(bnbEvents.arrays(["run", "subrun", "event", "nu_mu_cc_selected"]) )

extDaughters = uproot.open(InputFiles[2])["NuCCanalyzer"]["Daughters"]
trackExt     = pd.DataFrame(extDaughters.arrays(["track_range_mom_mu", "track_mcs_mom", "track_range_mom_p", "track_is_muon_candidate", "track_score", "track_chi2_proton", "track_chi2_muon", "track_dirx", "track_diry", "track_dirz", "run", "subrun", "event"] ) )
filteredExt  = pd.DataFrame(extEvents.arrays(["run", "subrun", "event", "nu_mu_cc_selected"]) )

overlayWeights = np.full(filteredEvents.shape[0], OverlayScale )

filteredEvents.insert(filteredEvents.shape[1], "DuplicatedEvent", filteredEvents.duplicated() ) #Tag the events which are duplicated
filteredEvents.insert(filteredEvents.shape[1], "mc_channel", [getChan(x, y) for x, y in zip(filteredEvents['mc_nu_interaction_type'], filteredEvents['mc_nu_ccnc'])] ) #Classify neutrino events based on CC / NC and event Type
filteredEvents.eval('mc_Ehad = mc_nu_energy - mc_nu_lepton_energy', inplace=True)#Insert the true energy transfer (nu)
filteredEvents.insert(filteredEvents.shape[1], "mc_expQ2", [getQ2(x, y, z) for x, y, z in zip(filteredEvents['mc_nu_energy'], filteredEvents['mc_nu_lepton_energy'], filteredEvents['mc_nu_lepton_theta'])] )
filteredEvents.insert(filteredEvents.shape[1], "mc_expW", [getW(x, y) for x, y in zip(filteredEvents['mc_Ehad'], filteredEvents['mc_expQ2'] ) ] )
filteredEvents.insert(filteredEvents.shape[1], "mc_expXbj", [getXbj(x, y) for x, y in zip(filteredEvents['mc_Ehad'], filteredEvents['mc_expQ2'] ) ] )
filteredEvents.insert(filteredEvents.shape[1], "mc_expY", [getInel(x, y) for x, y in zip(filteredEvents['mc_Ehad'], filteredEvents['mc_nu_energy'] ) ] )
filteredEvents.insert(filteredEvents.shape[1], "wgt", overlayWeights )

trackDaughters.insert(trackDaughters.shape[1], "phi", [getPhi(x, y) for x, y in zip(trackDaughters['track_diry'], trackDaughters['track_dirx'] ) ] )


extWeights      = np.full(filteredExt.shape[0],  ExtScale)

#filteredEvents.insert(filteredEvents.shape[1], "wgt", [getInel(x, y) for x, y in zip(filteredEvents['mc_Ehad'], filteredEvents['mc_nu_energy'] ) ] )

filteredData.insert(filteredData.shape[1], "DuplicatedEvent", filteredData.duplicated() ) #Tag the events which are duplicated
trackData.insert(trackData.shape[1], "phi", [getPhi(x, y) for x, y in zip(trackData['track_diry'], trackData['track_dirx'] ) ] )

filteredExt.insert(filteredExt.shape[1], "DuplicatedEvent", filteredExt.duplicated() ) #Tag the events which are duplicated
filteredExt.insert(filteredExt.shape[1], "wgt", extWeights )

trackExt.insert(trackExt.shape[1], "phi", [getPhi(x, y) for x, y in zip(trackExt['track_diry'], trackExt['track_dirx'] ) ] )

trackDaughters =  trackDaughters.set_index(['run', 'subrun', 'event'])
filteredEvents =  filteredEvents.set_index(['run', 'subrun', 'event'])

filteredData   =  filteredData.set_index(['run', 'subrun', 'event'])
trackData      =  trackData.set_index(['run', 'subrun', 'event'])

filteredExt    = filteredExt.set_index(['run', 'subrun', 'event'])
trackExt       = trackExt.set_index(['run', 'subrun', 'event'])

#create a dict of event info we want to associate with each daughter.
interactionInfo = {"DuplicatedEvent" : [], "mc_channel" : [], "nu_mu_cc_selected" : [], "mc_Ehad" : [], "mc_expQ2" : [], "mc_expXbj" : [], "mc_expY" : [], "mc_expW" : [], "wgt" : [] }

for index, row in trackDaughters.iterrows():

    #This gets around duplicate events. It's a bit clunky
    if type(filteredEvents.at[index, "DuplicatedEvent"]) is np.ndarray:
       for itype in interactionInfo:
           interactionInfo[itype].append( filteredEvents.at[index, itype][1] )
       
    else:
        for itype in interactionInfo:
           interactionInfo[itype].append( filteredEvents.at[index, itype] )
#print ((trackDaughters.loc[index]).query('track_is_muon_candidate == False & track_chi2_proton < @maxProtonChi2 & track_score > @minTrackScore')).count(axis=0)['phi']

dataInfo = {"DuplicatedEvent" : [], "nu_mu_cc_selected" : []}

for index, row in trackData.iterrows():
    
    #This gets around duplicate events. It's a bit clunky
    if type(filteredData.at[index, "DuplicatedEvent"]) is np.ndarray:
        for itype in dataInfo:
            dataInfo[itype].append( filteredData.at[index, itype][1] )

    else:
        for itype in dataInfo:
            dataInfo[itype].append( filteredData.at[index, itype] )

extInfo = {"DuplicatedEvent" : [], "nu_mu_cc_selected" : [], "wgt" : []}

for index, row in trackExt.iterrows():
    
    #This gets around duplicate events. It's a bit clunky
    if type(filteredExt.at[index, "DuplicatedEvent"]) is np.ndarray:
        for itype in extInfo:
            extInfo[itype].append( filteredExt.at[index, itype][1] )

    else:
        for itype in extInfo:
            extInfo[itype].append( filteredExt.at[index, itype] )

'''
for index, row in filteredEvents.iterrows():
    eventID = (int(row['run']), int(row['subrun']), int(row['event'])) #build a tuple of the run, subrun event number
    eventList.append(eventID)
    eventInt = pd.Series(row['mc_nu_interaction_type'])
    
    #print eventInt.repeat(row['num_matched_daughters'])
    interacionTypeDaughters = pd.concat( [interacionTypeDaughters, eventInt.repeat(row['num_daughters']) ], ignore_index=True )
#interacionTypeDaughters = pd.concat( [interacionTypeDaughters, eventInt], ignore_index=True )

print nEvents_passing
print interacionTypeDaughters.size
print len(trackDaughters.index)
'''
#interacionTypeDaughters = pd.Series(interactionList)

#filteredDaughters = trackDaughters.loc[ eventList ] #Only keep the daughters from the run, subrun, event tuple. This should now contain all the daughters of the neutrino events which pass the selection

#print interacionTypeDaughters

#print interacionTypeDaughters

muonMomentumRange   = (0.0, 2.0)
protonMomentumRange = (0.0, 1.5)
phiRange = (-1.0, 1.0)
for itype in interactionInfo:
   trackDaughters.insert(trackDaughters.shape[1], itype, interactionInfo[itype])

for itype in dataInfo:
   trackData.insert(trackData.shape[1], itype, dataInfo[itype])

for itype in extInfo:
   trackExt.insert(trackExt.shape[1], itype, extInfo[itype])

extMuons     = trackExt.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == True & track_score > @minTrackScore')
muonTracks   = trackDaughters.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == True & track_score > @minTrackScore')

#Do we see any primary tracks which fail the min score cut?
print "Selected MC events failing track score cut: %d" % len(trackDaughters.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == True & track_score < @minTrackScore').index)
print "Selected BNB events failing track score cut: %d" % len(trackData.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == True & track_score < @minTrackScore').index)

protonTracks = trackDaughters.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == False & track_chi2_proton < @maxProtonChi2 & track_score > @minTrackScore')
print protonTracks
print protonTracks.count(level="event", axis=0)
extProtons   = trackExt.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == False & track_chi2_proton < @maxProtonChi2 & track_score > @minTrackScore')


mcsMomentumPassed   = muonTracks['track_mcs_mom']
momentumStack = [muonTracks.query('mc_channel == "QE"')['track_mcs_mom'].to_numpy(), muonTracks.query('mc_channel == "RES"')['track_mcs_mom'].to_numpy(), muonTracks.query('mc_channel == "DIS"')['track_mcs_mom'].to_numpy(), muonTracks.query('mc_channel == "2p2h"')['track_mcs_mom'].to_numpy(), muonTracks.query('mc_channel == "NC / Other"')['track_mcs_mom'].to_numpy(), extMuons['track_mcs_mom'].to_numpy() ]
muonPhiStack = [muonTracks.query('mc_channel == "QE"')['phi'].to_numpy(), muonTracks.query('mc_channel == "RES"')['phi'].to_numpy(), muonTracks.query('mc_channel == "DIS"')['phi'].to_numpy(), muonTracks.query('mc_channel == "2p2h"')['phi'].to_numpy(), muonTracks.query('mc_channel == "NC / Other"')['phi'].to_numpy(), extMuons['phi'].to_numpy() ]
protonPhiStack = [protonTracks.query('mc_channel == "QE"')['phi'].to_numpy(), protonTracks.query('mc_channel == "RES"')['phi'].to_numpy(), protonTracks.query('mc_channel == "DIS"')['phi'].to_numpy(), protonTracks.query('mc_channel == "2p2h"')['phi'].to_numpy(), protonTracks.query('mc_channel == "NC / Other"')['phi'].to_numpy(), extProtons['phi'].to_numpy() ]
muonWeights   = [muonTracks.query('mc_channel == "QE"')['wgt'].to_numpy(), muonTracks.query('mc_channel == "RES"')['wgt'].to_numpy(), muonTracks.query('mc_channel == "DIS"')['wgt'].to_numpy(), muonTracks.query('mc_channel == "2p2h"')['wgt'].to_numpy(), muonTracks.query('mc_channel == "NC / Other"')['wgt'].to_numpy(), extMuons['wgt'].to_numpy() ]
protonWeights = [protonTracks.query('mc_channel == "QE"')['wgt'].to_numpy(), protonTracks.query('mc_channel == "RES"')['wgt'].to_numpy(), protonTracks.query('mc_channel == "DIS"')['wgt'].to_numpy(), protonTracks.query('mc_channel == "2p2h"')['wgt'].to_numpy(), protonTracks.query('mc_channel == "NC / Other"')['wgt'].to_numpy(), extProtons['wgt'].to_numpy() ]
protonStack = [protonTracks.query('mc_channel == "QE"')['track_range_mom_p'].to_numpy(), protonTracks.query('mc_channel == "RES"')['track_range_mom_p'].to_numpy(), protonTracks.query('mc_channel == "DIS"')['track_range_mom_p'].to_numpy(), protonTracks.query('mc_channel == "2p2h"')['track_range_mom_p'].to_numpy(), protonTracks.query('mc_channel == "NC / Other"')['track_range_mom_p'].to_numpy(), extProtons['track_range_mom_p'].to_numpy() ]

rangeMomentumPassed = muonTracks['track_range_mom_mu']

dataMuons   = trackData.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == True & track_score > @minTrackScore')['track_mcs_mom']
dataProtons = trackData.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == False & track_chi2_proton < @maxProtonChi2 & track_score > @minTrackScore')['track_range_mom_p']

dataMuonPhi   = trackData.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == True & track_score > @minTrackScore')['phi']
dataProtonPhi = trackData.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == False & track_chi2_proton < @maxProtonChi2 & track_score > @minTrackScore')['phi']


fig = plt.figure()
plt.hist(momentumStack, bins=18, stacked=True, range=muonMomentumRange, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights = muonWeights )
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Ext'])
plt.title("MCS Momentum")
plt.xlabel("Momentum (GeV/c)")
plt.ylabel("Number of Primary Muons")
data_hist = dataify(dataMuons.to_numpy(), 18, muonMomentumRange)
#print data_hist[0]
#print data_hist[1]
plt.errorbar(data_hist[0], data_hist[1], yerr=data_hist[2], fmt='o', color='black')
plt.xlim(0.0, 2.0)
plt.show()

plt.hist(protonStack, bins=50, stacked=True, range=protonMomentumRange, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights = protonWeights )
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Ext'])
plt.title("Proton Range Momentum")
plt.xlabel("Momentum (GeV/c)")
plt.ylabel("Number of Protons")
data_proton_hist = dataify(dataProtons.to_numpy(), 50, protonMomentumRange)
plt.errorbar(data_proton_hist[0], data_proton_hist[1], yerr=data_proton_hist[2], fmt='o', color='black')
plt.show()

plt.hist(muonPhiStack, bins=50, stacked=True, range=phiRange, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights = muonWeights )
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Ext'])
plt.title("Muon Phi Angle")
plt.xlabel("Angle / pi (radians)")
plt.ylabel("Number of Primary Muons")
data_phi_hist = dataify(dataMuonPhi.to_numpy(), 50, phiRange)
plt.errorbar(data_phi_hist[0], data_phi_hist[1], yerr=data_phi_hist[2], fmt='o', color='black')
plt.show()

plt.hist(protonPhiStack, bins=25, stacked=True, range=phiRange, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights = protonWeights )
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Ext'])
plt.title("Proton Phi Angle")
plt.xlabel("Angle / pi (radians)")
plt.ylabel("Number of Primary Muons")
data_proton_phi_hist = dataify(dataProtonPhi.to_numpy(), 25, phiRange)
plt.errorbar(data_proton_phi_hist[0], data_proton_phi_hist[1], yerr=data_proton_phi_hist[2], fmt='o', color='black')
plt.show()

'''
fig = plt.figure()
trueEnergy = pd.DataFrame(overlayEvents.arrays(["mc_nu_energy", "nu_mu_cc_selected", "mc_nu_interaction_type", "mc_nu_ccnc"]) )
data_hist = dataify(trueEnergy['mc_nu_energy'].to_numpy(), 50)
result = trueEnergy.query('nu_mu_cc_selected == True')

trueEnergy['mc_nu_energy'].hist(bins=50, grid=False)
result['mc_nu_energy'].hist(bins=50, grid=False)
plt.errorbar(data_hist[0], data_hist[1], yerr=data_hist[2], fmt='o', color='black')
plt.xlabel("Neturino Energy (GeV)")
plt.legend(['No Cuts', 'Selected'])
#df = trueEnergyDict)

#print array

plt.show()
'''
sys.exit()
