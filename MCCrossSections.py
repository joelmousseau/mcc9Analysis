import ROOT, math, sys, os
import uproot
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import csv
import re
from mpl_toolkits.mplot3d import Axes3D
from ROOT import TH1, TAxis, gROOT, TCanvas
from scipy.constants import Avogadro as avosNumber
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



argonDensity         = 1.3954 #g/cm3
nProtonsPerNucleus   = 18
nNeutronsPerNucleus = 40 - nProtonsPerNucleus
molMass              = 39.95 #g/mol
TPC_X                = 256.4 #cm
TPC_Y                = 300   #cm
TPC_Z                = 1000  #cm
TPC_Volume           = TPC_X*TPC_Y*TPC_Z #cm3
nTargetProtons       = (nProtonsPerNucleus*argonDensity*TPC_Volume*avosNumber / molMass)
nTargetNeutrons      = (nNeutronsPerNucleus*argonDensity*TPC_Volume*avosNumber / molMass)

neutrinoEnergyRange = (0.0, 10.0)
invariantMassRange = (0.0, 10.0)
phiRange           = (-1.0, 1.0)
nBins              = 100

InputFiles      = ["/uboone/app/users/wvdp/RootTrees/v20/run1/nucc_nu_overlay_run1_mcc9.root", "/uboone/app/users/wvdp/RootTrees/v20/run1/nucc_on_data_run1_mcc9.root", "/uboone/app/users/wvdp/RootTrees/v20/run1/nucc_off_data_run1_mcc9.root"]
fluxFileName    = "/uboone/app/users/afurmans/CCNproton_ccinc/analysis/UBAna/Flux/MCC8_FluxHistograms_Uncertainties.root"
splineFileName  = "/cvmfs/uboone.opensciencegrid.org/products/genie_xsec/v3_00_04a/NULL/G1810a0211a-k250-e1000/data/xsec_graphs.root"
splineDir       = "nu_mu_Ar40"
fluxFile        = ROOT.TFile.Open(fluxFileName)
splineFile      = ROOT.TFile.Open(splineFileName)
fluxHistogram   = fluxFile.Get("numu/numu_CV_AV_TPC")
splineGraphs    = ["%s/qel_cc_n" % splineDir, "%s/mec_cc" % splineDir]


fluxNeutrinos = np.array([fluxHistogram.GetBinContent(bin+1) for bin in range(fluxHistogram.GetNbinsX())], dtype=float)
fluxBins      = np.array([fluxHistogram.GetBinCenter(bin+1) for bin in range(fluxHistogram.GetNbinsX())], dtype=float)
splineList    = []

for splineName in splineGraphs:
    graph = splineFile.Get(splineName)
    splineList.append(np.array([graph.Eval(x) for x in range(graph.GetN())]) )

print "Number of Neutrinos %e Array Sum %e" % (fluxHistogram.Integral(), np.sum(fluxNeutrinos, dtype=float))
print "Number of protons %e. Number of neutrons %e" % (nTargetProtons, nTargetNeutrons)

binWidth = (neutrinoEnergyRange[1] - neutrinoEnergyRange[0])/nBins
print binWidth

crossSectionDenom = (nTargetNeutrons+nTargetProtons)*np.sum(fluxNeutrinos, dtype=float)*binWidth

#OverlayScale  = 1.0
ExtScale     = 0.97
numMCTemplates = 6
empty = []

#plt.figure()

overlayEvents = uproot.open(InputFiles[0])["NuCCanalyzer"]["Event"]

overlayPOT    = uproot.open(InputFiles[0])["NuCCanalyzer"]["subruns"]

mcPOT         = pd.Series(overlayPOT.array("pot")).sum()

#Divide by area
fluxPOT = 1.6e20
fluxScale   = 2.43e11

OverlayScale = fluxPOT/(mcPOT*crossSectionDenom)
print "MC POT: %e Overlay Scale: %.3f" % (mcPOT, OverlayScale)

overlayEvents   = pd.DataFrame(overlayEvents.arrays(["run", "subrun", "event", "mc_nu_interaction_type", "mc_nu_ccnc", "nu_mu_cc_selected", "mc_nu_lepton_energy", "mc_nu_energy", "mc_nu_lepton_theta"]) )

overlayWeights = np.full(overlayEvents.shape[0], OverlayScale )


overlayEvents.insert(overlayEvents.shape[1], "DuplicatedEvent", overlayEvents.duplicated() ) #Tag the events which are duplicated
overlayEvents.insert(overlayEvents.shape[1], "mc_channel", [getChan(x, y) for x, y in zip(overlayEvents['mc_nu_interaction_type'], overlayEvents['mc_nu_ccnc'])] ) #Classify neutrino events based on CC / NC and event Type
overlayEvents.eval('mc_Ehad = mc_nu_energy - mc_nu_lepton_energy', inplace=True)#Insert the true energy transfer (nu)
overlayEvents.insert(overlayEvents.shape[1], "mc_expQ2", [getQ2(x, y, z) for x, y, z in zip(overlayEvents['mc_nu_energy'], overlayEvents['mc_nu_lepton_energy'], overlayEvents['mc_nu_lepton_theta'])] )
overlayEvents.insert(overlayEvents.shape[1], "mc_expW", [getW(x, y) for x, y in zip(overlayEvents['mc_Ehad'], overlayEvents['mc_expQ2'] ) ] )
overlayEvents.insert(overlayEvents.shape[1], "mc_expXbj", [getXbj(x, y) for x, y in zip(overlayEvents['mc_Ehad'], overlayEvents['mc_expQ2'] ) ] )
overlayEvents.insert(overlayEvents.shape[1], "mc_expY", [getInel(x, y) for x, y in zip(overlayEvents['mc_Ehad'], overlayEvents['mc_nu_energy'] ) ] )
overlayEvents.insert(overlayEvents.shape[1], "wgt", overlayWeights )

#filteredEvents.insert(filteredEvents.shape[1], "wgt", [getInel(x, y) for x, y in zip(filteredEvents['mc_Ehad'], filteredEvents['mc_nu_energy'] ) ] )

stackCode = ''' [overlayEvents.query('mc_channel == "QE"')['VAR'].to_numpy(), overlayEvents.query('mc_channel == "RES"')['VAR'].to_numpy(), overlayEvents.query('mc_channel == "DIS"')['VAR'].to_numpy(), overlayEvents.query('mc_channel == "2p2h"')['VAR'].to_numpy()] '''

weightsArray = [overlayEvents.query('mc_channel == "QE"')['wgt'].to_numpy(), overlayEvents.query('mc_channel == "RES"')['wgt'], overlayEvents.query('mc_channel == "DIS"')['wgt'], overlayEvents.query('mc_channel == "2p2h"')['wgt'].to_numpy()]

exec( "energyStack = " + re.sub(r'VAR', 'mc_nu_energy', stackCode) )
exec( "wStack = " + re.sub(r'VAR', 'mc_expW', stackCode) )




fig = plt.figure()
plt.hist(energyStack, bins=nBins, stacked=True, range=neutrinoEnergyRange, color = ['b', 'g', 'y', 'r'], weights=weightsArray )
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other'])
plt.title("Neutrino Energy")
plt.xlabel("Energy (GeV)")
plt.ylabel("Number of Events / %e POT" % mcPOT)
plt.xlim(0.0, 7.0)
plt.show()

fig = plt.figure()
plt.hist(wStack, bins=100, stacked=True, range=invariantMassRange, color = ['b', 'g', 'y', 'r'], weights=weightsArray )
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other'])
plt.title("Invariant Mass")
plt.xlabel("W (GeV/c^2)")
plt.ylabel("Number of Events / %e POT" % mcPOT)
plt.xlim(0.0, 3.0)
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
