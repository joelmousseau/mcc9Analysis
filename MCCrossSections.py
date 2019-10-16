import ROOT, math, sys, os
import uproot
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import csv
import re
import scipy as sci
from mpl_toolkits.mplot3d import Axes3D
from ROOT import TH1, TAxis, gROOT, TCanvas
from scipy.constants import Avogadro as avosNumber
from scipy import stats
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
    W2 = getW2(Ehad, Q2)
    if(W2 >= 0.0):
      return math.sqrt(W2)
    else:
      return -1.0

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


neutrinoPDG          = 14
argonDensity         = 1.3954 #g/cm3
#argonNucleons        = 40
argonNucleons        = 1
#nProtonsPerNucleus   = 18
nProtonsPerNucleus   = 0.5
nNeutronsPerNucleus  = argonNucleons - nProtonsPerNucleus
molMass              = 39.95 #g/mol
TPC_X                = 256.35 #cm
TPC_Y                = 233   #cm
TPC_Z                = 1000  #cm
TPC_Volume           = TPC_X*TPC_Y*TPC_Z #cm3
TPC_Area             = TPC_X*TPC_Y #cm2
#TPC_Area             = 1.0
CryoR                = 191.61 #cm from gdml
CryoL                = 1086.49 #cm from gdml
CryoVolume           = np.pi*CryoL*CryoR**2
nTargetProtons       = (nProtonsPerNucleus*argonDensity*CryoVolume*avosNumber / molMass)
nTargetNeutrons      = (nNeutronsPerNucleus*argonDensity*CryoVolume*avosNumber / molMass)

neutrinoEnergyRange = (0.0, 10.0)
invariantMassRange = (0.0, 10.0)
phiRange           = (-1.0, 1.0)
nBins              = 200

fluxScale   = 2.43e11
detectorPositionScale = 1.029
#fluxScale   = 1.0
fluxPOT     = fluxScale
binWidth    = (neutrinoEnergyRange[1] - neutrinoEnergyRange[0])/nBins

InputFiles      = ["/uboone/app/users/wvdp/RootTrees/v20/run1/nucc_nu_overlay_run1_big_mcc9.root", "/uboone/app/users/wvdp/RootTrees/v20/run1/nucc_on_data_run1_mcc9.root", "/uboone/app/users/wvdp/RootTrees/v20/run1/nucc_off_data_run1_mcc9.root"]
fluxFileName    = "/uboone/app/users/afurmans/CCNproton_ccinc/analysis/UBAna/Flux/MCC8_FluxHistograms_Uncertainties.root"
splineFileName  = "/cvmfs/uboone.opensciencegrid.org/products/genie_xsec/v3_00_04a/NULL/G1810a0211a-k250-e1000/data/xsec_graphs.root"
splineDir       = "nu_mu_Ar40"
fluxFile        = ROOT.TFile.Open(fluxFileName)
splineFile      = ROOT.TFile.Open(splineFileName)
fluxHistogram   = fluxFile.Get("numu/numu_CV_AV_TPC")
splineGraphs    = ["%s/tot_cc" % splineDir, "%s/qel_cc_n" % splineDir, "%s/mec_cc" % splineDir]

overlayEvents = uproot.open(InputFiles[0])["NuCCanalyzer"]["Event"]

overlayPOT    = uproot.open(InputFiles[0])["NuCCanalyzer"]["subruns"]

mcPOT         = pd.Series(overlayPOT.array("pot")).sum()
#mcPOT = 1.0

fluxBins = np.empty(fluxHistogram.GetNbinsX()+1, dtype=float)
fluxHistogram.GetXaxis().GetLowEdge(fluxBins)
fluxNeutrinos = ((mcPOT*detectorPositionScale)/(TPC_Area*fluxScale))*np.array([fluxHistogram.GetBinContent(bin+1) for bin in range(fluxHistogram.GetNbinsX()+1)], dtype=float)
fluxBins = fluxBins[:-1]
fluxNeutrinos = fluxNeutrinos[:-1]
#print fluxBins
#print len(fluxNeutrinos)


#fluxBins      = np.array([fluxHistogram.GetBinCenter(bin+1) for bin in range(fluxHistogram.GetNbinsX())], dtype=float)
splineList    = []
splineBins    = []


for splineName in splineGraphs:
    
    graph = splineFile.Get(splineName)
    splineList.append(np.array(list(graph.GetY() ) ) )
    splineBins.append(np.array(list(graph.GetX() ) ) )

#totalNumus        = mcPOT*np.sum(fluxNeutrinos, dtype=float)/(TPC_Area*fluxScale)
totalNumus        = np.sum(fluxNeutrinos)
crossSectionDenom = (nTargetNeutrons+nTargetProtons)

#print splineList[0]
#print splineBins[0]

#print splineList[0]
print "Number of Neutrinos %e Array Sum %e" % (fluxHistogram.Integral(), np.sum(fluxNeutrinos, dtype=float))
print "Number of protons %e. Number of neutrons %e" % (nTargetProtons, nTargetNeutrons)
print "Number of Neutrinos per cm2 (scaled) %e" % totalNumus
print "Bin Width %.3f" % binWidth


#print binWidth



#OverlayScale  = 1.0
ExtScale     = 0.97
numMCTemplates = 6
empty = []

#plt.figure()



#Divide by area


OverlayScale = argonNucleons*1e38/(crossSectionDenom)
print "MC POT: %e Overlay Scale: %e" % (mcPOT, OverlayScale)

overlayEvents   = pd.DataFrame(overlayEvents.arrays(["run", "subrun", "event", "mc_nu_interaction_type", "mc_nu_ccnc", "nu_mu_cc_selected", "mc_nu_lepton_energy", "mc_nu_energy", "mc_nu_lepton_theta", "mc_nu_pdg"]) )

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

stackCode = ''' [overlayEvents.query('mc_channel == "QE" & mc_nu_pdg == @neutrinoPDG')['VAR'].to_numpy(), overlayEvents.query('mc_channel == "RES" & mc_nu_pdg == @neutrinoPDG')['VAR'].to_numpy(), overlayEvents.query('mc_channel == "DIS" & mc_nu_pdg == @neutrinoPDG')['VAR'].to_numpy(), overlayEvents.query('mc_channel == "2p2h" & mc_nu_pdg == @neutrinoPDG')['VAR'].to_numpy()] '''

#stackSumCode = ''' overlayEvents.query('mc_channel == "QE"')['VAR'].to_numpy() + overlayEvents.query('mc_channel == "RES"')['VAR'].to_numpy() + overlayEvents.query('mc_channel == "DIS"')['VAR'].to_numpy() + overlayEvents.query('mc_channel == "2p2h"')['VAR'].to_numpy() '''

weightsArray = [overlayEvents.query('mc_channel == "QE" & mc_nu_pdg == @neutrinoPDG')['wgt'].to_numpy(), overlayEvents.query('mc_channel == "RES" & mc_nu_pdg == @neutrinoPDG')['wgt'], overlayEvents.query('mc_channel == "DIS" & mc_nu_pdg == @neutrinoPDG')['wgt'], overlayEvents.query('mc_channel == "2p2h" & mc_nu_pdg == @neutrinoPDG')['wgt'].to_numpy()]

exec( "energyStack = " + re.sub(r'VAR', 'mc_nu_energy', stackCode) )
#exec( "energyStackSum = " + re.sub(r'VAR', 'mc_nu_energy', stackSumCode) )
exec( "wStack = " + re.sub(r'VAR', 'mc_expW', stackCode) )





crossSectionStack = []
fluxStack         = []
sum     = np.zeros(len(fluxBins[:-1]), dtype=float)
sumErrs = np.zeros(len(fluxBins[:-1]), dtype=float)
#print fluxBins[:-1]
bin_centers       = (fluxBins[:-1] + fluxBins[1:]) / 2
#print bin_centers
for i in range(len(energyStack)):
  counts, bin_edges, something = sci.stats.binned_statistic(energyStack[i], energyStack[i], "count", bins=fluxBins, range=neutrinoEnergyRange)
  crossSectionStack.append(np.nan_to_num(OverlayScale*counts/fluxNeutrinos[:-1]) )

  sumErrs =  np.add(sum, counts)
  sum = np.add(sum, np.nan_to_num(OverlayScale*counts/fluxNeutrinos[:-1]) )
  fluxStack.append(fluxBins)


#counts, bin_edges, something = plt.hist(crossSectionStack, bins=nBins, stacked=True, range=neutrinoEnergyRange, color = ['b', 'g', 'y', 'r']) #unscaled version to get the errors right
#plt.show()

#print sum
#plt.hist(energyStack, bins=nBins, stacked=True, range=neutrinoEnergyRange, color = ['b', 'g', 'y', 'r'], weights=weightsArray )
sumErrs = OverlayScale*np.sqrt(sumErrs)/fluxNeutrinos[:-1]

bin_centers       = (fluxBins[:-1] + fluxBins[1:]) / 2
splineValues = np.interp(bin_centers, splineBins[0], splineList[0], right=0.0)
fig = plt.figure(1)
#plt.scatter(fluxStack[0], sum, color = "black")
plt.errorbar(fluxBins[:-1], sum, yerr=sumErrs, fmt='o', color='black')
plt.plot(fluxBins[:-1], splineValues, color='red')
plt.legend(['Spline', 'CCInc'])
plt.title("Neutrino Energy")
plt.xlabel("Energy (GeV)")
plt.ylabel("Number of Events / %e POT" % mcPOT)
plt.xlim(0.0, 5.0)
plt.ylim(0.0, 200.0)
#plt.show()

fig = plt.figure(2)
plt.hist(energyStack, bins=fluxBins, stacked=True, range=neutrinoEnergyRange, weights=weightsArray, color = ['b', 'g', 'y', 'r'])
plt.plot(fluxBins[:-1], 33*fluxNeutrinos[:-1], color='black')
plt.legend(['MCC8 Flux', 'QE', 'RES', 'DIS', '2p2h'])
plt.title("Neutrino Energy")
plt.xlabel("E_{\nu} (GeV)")
plt.ylabel("Number of Events / %e POT" % mcPOT)
plt.xlim(0.0, 5.0)
#plt.show()

fig = plt.figure(3)
plt.hist(energyStack, bins=fluxBins[:-1], stacked=True, density=True, range=neutrinoEnergyRange, weights=weightsArray, color = ['b', 'g', 'y', 'r'])
fluxXSec = np.multiply(fluxNeutrinos[:-1], splineValues)
norm = np.trapz(fluxXSec, dx=0.05)
print "Normalization: %.2e" % norm
print fluxBins
print (fluxXSec/norm)

plt.plot(fluxBins[:-1], fluxXSec/norm, color='black')
plt.legend(['Flux X Xsection', 'QE', 'RES', 'DIS', '2p2h'])
plt.title("Neutrino Energy")
plt.xlabel("Energy (GeV)")
plt.ylabel("Abr.")
plt.xlim(0.0, 5.0)
plt.show()

sys.exit()
