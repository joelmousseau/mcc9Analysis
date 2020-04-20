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

def splitAndSort(inputArray, nDivisions):
  preSplit = np.sort(inputArray)
  dropEntries = preSplit.size % nDivisions
  print "Dropping %d entries prior to split" % dropEntries
  #preSplit = preSplit[dropEntries:]

  return np.array_split(preSplit, nDivisions)


def makeMCHistogram(mc, channel, binRange, nBins, filename, Titles):
  dir_name = "PlotDir"
  colors = {"QE":'b', "RES":'g', "DIS":'y', "2p2h":'r', "NC / Other":'grey', "Ext":'magenta'}

  plt.hist(mc, bins=nBins, stacked=False, range=binRange, color = colors[channel])
  plt.legend([channel])
  try:
    plotTitle, xAxisTitle, yAxisTitle = Titles
  except(ValueError):
    print "Pleast provide three titles for plot name, x and y axis names" 
    plotTitle  = ""
    xAxisTitle = ""
    yAxisTitle = "" 
  plt.title(plotTitle)
  plt.xlabel(xAxisTitle)
  plt.ylabel(yAxisTitle)
  plt.savefig("%s/%s%s.png" % ( dir_name, filename, channel.replace(" / Other", "")) )
  plt.close()

def make2DMCHistogram(mc, channel, binRange, nBins, filename, Titles):
  dir_name = "PlotDir"
  colors = {"QE":'b', "RES":'g', "DIS":'y', "2p2h":'r', "NC / Other":'grey', "Ext":'magenta'}
  zMin = 0.01

  try:
    xBins, yBins = binRange

  except(ValueError):
    print "Please provide a range of bins for each axis"
    return

  plt.hist2d(mc[0], mc[1], bins=nBins, range=binRange, cmin=zMin, normed=True)
  plt.legend([channel])

  try:
    plotTitle, xAxisTitle, yAxisTitle = Titles
  except(ValueError):
    print "Pleast provide three titles for plot name, x and y axis names" 
    plotTitle  = ""
    xAxisTitle = ""
    yAxisTitle = "" 
  plt.title(plotTitle)
  plt.xlabel(xAxisTitle)
  plt.ylabel(yAxisTitle)
  plt.savefig("%s/%s%s.png" % ( dir_name, filename, channel.replace(" / Other", "")) )
  plt.close()

def makeDataMCHistogram(mcList, mcWeights, dataList, binRange, nBins, filename, Titles):
  dir_name = "PlotDir"
  plt.hist(mcList, bins=nBins, stacked=True, range=binRange, color = ['b', 'g', 'y', 'r', 'grey', 'gold', 'magenta'], weights = mcWeights )
  plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Dirt', 'Ext'])
  #plotTitle, xAxisTitle, yAxisTitle =  Titles
  try:
    plotTitle, xAxisTitle, yAxisTitle = Titles
  except(ValueError):
    print "Pleast provide three titles for plot name, x and y axis names" 
    plotTitle  = ""
    xAxisTitle = ""
    yAxisTitle = "" 
  plt.title(plotTitle)
  plt.xlabel(xAxisTitle)
  plt.ylabel(yAxisTitle)
  
  data_hist = dataify(dataList, nBins, binRange)
  plt.errorbar(data_hist[0], data_hist[1], yerr=data_hist[2], fmt='o', color='black')
  plt.savefig("%s/%s.png" % ( dir_name, filename) )
  plt.close()

  makeDataMCRatioHistogram(mcList, mcWeights, dataList, binRange, nBins, filename, Titles)

def makeDataMCRatioHistogram(mcList, mcWeights, dataList, binRange, nBins, filename, Titles):
  dir_name  = "PlotDir"
  mcSum = np.full(nBins, 0.0 )
  for mc, weight in zip(mcList, mcWeights):
     mc_hist   = np.histogram(mc, bins=nBins, range=binRange, weights = weight )
     np.add(mc_hist[0], mcSum, out=mcSum)
  data_hist = dataify(dataList, nBins, binRange)
  print "MC sum: %e" % np.sum(mcSum)
  print "Data sum: %e" % np.sum(data_hist[1]) 
  ratio = np.divide(data_hist[1], mcSum)
  err   = np.multiply(ratio, np.divide(1.0, data_hist[2]))
  np.nan_to_num(ratio, copy=False)
  np.nan_to_num(err, copy=False)

  plt.errorbar(data_hist[0], ratio, yerr=err, fmt='o', color='black') #This ignores MC stats.
  try:
    plotTitle, xAxisTitle, yAxisTitle = Titles
  except(ValueError):
    print "Pleast provide three titles for plot name, x and y axis names" 
    plotTitle  = ""
    xAxisTitle = ""
    yAxisTitle = "" 
  plt.title(plotTitle)
  plt.xlabel(xAxisTitle)
  plt.ylabel(yAxisTitle)
  plt.savefig("%s/%sRatio.png" % ( dir_name, filename) )
  plt.close()
  

def KEtoMomentum(T, restMass):
    TotalE = T + restMass
    return math.sqrt(TotalE**2 - restMass**2)
    '''
    0.2 + 0.1 = 0.3 
    sqrt(0.3^2 - 0.1^2) = sqrt(0.09 - 0.01)
    sqrt(0.08) = 0.28
    '''

def getMultibin(dataSequence, weightsSequence):
  
  #Munerahs Bins in KE
  T_mu_bins         = (0.0, 0.12, 0.34, 0.46, 2.0)
  T_p_bins          = (0.0, 0.009, 0.15, 0.21, 1.5)

  #Convert these to momentum
  #P_mu_bins         = tuple(KEtoMomentum(T, 0.1) for T in T_mu_bins)
  #P_p_bins          = tuple(KEtoMomentum(T, 1.0) for T in T_p_bins)
  P_mu_bins          = (0.12, 0.34, 0.52, 0.80, 15.0)
  P_p_bins          = (0.1, 0.42, 0.56, 0.74, 3.0) 
  cos_theta_mu_bins = (-1.0, 0.48, 0.75, 0.89, 1.0)
  cos_theta_p_bins  = (-1.0, 0.48, 0.75, 0.89, 1.0)
  phi_mu_p_bins    = (0.0, 0.69, 0.95, 1.04, 2.0)
  

  bins         = [P_mu_bins, P_p_bins, cos_theta_mu_bins, cos_theta_p_bins, phi_mu_p_bins]
  #bins          = [T_mu_bins, T_p_bins]
  return sci.stats.binned_statistic_dd(dataSequence, weightsSequence, statistic="sum", bins=bins)

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

def getParticle(pdg):
    if (pdg == 13):
      return "muon"
    elif(pdg == 2212):
      return "proton"
    elif(pdg == 211):
      return "pion"
    elif(pdg == 11):
      return "electron"
    else:
      return "other"      

def getChanWeight(interaction, isNC):
    templateFactors = [6.47365e-01, 1.20327e-07, 6.02801e-08, 2.71038e-08, 3.42514e-01, 1.00000e-02]
    if (isNC):
      return (1.0 + templateFactors[4])
    elif(interaction == 0):
      return (1.0 + templateFactors[0])
    elif(interaction == 1):
      return (1.0 + templateFactors[1])
    elif(interaction == 2):
      return (1.0 + templateFactors[2])
    elif(interaction == 10):
      return (1.0 + templateFactors[3])



def dotProduct(df1, df2, dim1, dim2, dimList):
   returnSeries = pd.Series(0.0, index=df1.index.copy())
  # print df1['%s' % dim1]*df1['%s' % "track_dirx"]
  # print df2['%s' % dim2]*df2['%s' % "track_dirx"]
  # print df1['%s' % dim1]*df1['%s' % "track_dirx"] + df2['%s' % dim2]*df2['%s' % "track_dirx"]
   denomSeries = df1['%s' % dim1]*df2['%s' % dim2]
   #print denomSeries

   for dim in dimList:    
       returnSeries += (df1['%s' % dim1]*df1['%s' % dim]*df2['%s' % dim2]*df2['%s' % dim])/denomSeries
       #print returnSeries
   #print returnSeries
   return returnSeries.to_numpy()

def getPhi(pY, pX):
    return (np.arctan2(pY, pX) / np.pi)

def getTheta(pTotal, pZ):
    return (pZ / pTotal)

def getQ2(nuE, nuMu, thetaMu):
    return 4*nuE*nuMu*math.pow(math.sin(thetaMu/2), 2)

#Containment taken from mcc8 Nproton analysis, note 1065
def isContained(xStart, yStart, zStart, xEnd, yEnd, zEnd):
    return (checkContained(xStart, yStart, zStart) and checkContained(xEnd, yEnd, zEnd))

def checkContained(x, y, z):
    if(x > 245.0 or x < 10.0):
      return False
    elif(y > 100.0 or y < -100.0):  
      return False
    elif(z > 1030.0 or z <10.0):
      return False
    else:  
      return True

def getBestP(isContained, pOne, pTwo):
  return (pOne if isContained else pTwo)

def getW2(Ehad, Q2, interaction):
    targetMass = 0.989
    #if it's a 2p 2 h, you have 2ps as your target! Wakka Wakka Wakka!
    if(interaction == "2p2h"):
      targetMass = 2*targetMass
    return 2*targetMass*Ehad + math.pow(targetMass, 2) - Q2

def getW(Ehad, Q2, interaction):
    W2 = getW2(Ehad, Q2, interaction)
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

InputFiles = ["/uboone/data/users/joelam/stv-ntuples-new/numu_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/bnb_5e19_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/extC1_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/dirt_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/extC2_run1.root"]

#OverlayScale  = 1.0
ExtScale     = 0.97
numMCTemplates = 6
empty = []


maxProtonChi2 = 88.0
minTrackScore = 0.5

#Python library to read in ROOT ntuples files.
overlayEvents = uproot.open(InputFiles[0])["NuCCanalyzer"]["Event"]
bnbEvents     = uproot.open(InputFiles[1])["NuCCanalyzer"]["Event"]
extEventsC1    = uproot.open(InputFiles[2])["NuCCanalyzer"]["Event"]
extEventsC2   = uproot.open(InputFiles[4])["NuCCanalyzer"]["Event"]
dirtEvents    = uproot.open(InputFiles[3])["NuCCanalyzer"]["Event"]

overlayPOT    = uproot.open(InputFiles[0])["NuCCanalyzer"]["subruns"]
dirtPOT       = uproot.open(InputFiles[3])["NuCCanalyzer"]["subruns"]

#Scale factors, because we generate more simulation than data. We also do not take an equal ammount of on and off beam data (though it is close)
mcPOT         = pd.Series(overlayPOT.array("pot"))
sumPOT        = mcPOT.sum()

sumDirtPOT    = (pd.Series(dirtPOT.array("pot"))).sum()


useC2 = False 
dataPOT       = 4.08e+19
bnbSpills     = 9045263.0
extTriggersC1 = 33630174.0
extTriggersC2 = 31587147.0
extTriggers   = extTriggersC1
maxEvents     = 35000




ExtTemplateWeight = (1.0 + 1.00000e-02)

#Create frames of the event tree (which has information about the interaction) and the duaghters tree (which has information about the particles within the interaction).
#Do this for "overlay" simulation, beam data, and off-beam data
overlayDaughters = uproot.open(InputFiles[0])["NuCCanalyzer"]["Daughters"]
trackDaughters   = pd.DataFrame(overlayDaughters.arrays(["track_range_mom_mu", "track_mcs_mom", "track_range_mom_p", "track_is_muon_candidate", "track_score", "track_chi2_proton", "track_chi2_muon", "track_dirx", "track_diry", "track_dirz", "vx", "vy", "vz", "track_endx", "track_endy", "track_endz", "mc_pdg", "run", "subrun", "event"] ) )
filteredEvents   = pd.DataFrame(overlayEvents.arrays(["run", "subrun", "event", "mc_nu_interaction_type", "mc_nu_ccnc", "nu_mu_cc_selected", "mc_nu_lepton_energy", "mc_nu_energy", "mc_nu_lepton_theta"]) )
overlayCVWeights = pd.read_csv("AllBNBWeights.csv", names=["run", "subrun", "event", "wgt_tune"], dtype={"run" : int, "subrun" : int, "event" : int, "wgt_tune" : float})
dirtCVWeights    = pd.read_csv("DirtWeights1.csv", names=["run", "subrun", "event", "wgt_tune"], dtype={"run" : int, "subrun" : int, "event" : int, "wgt_tune" : float})
overlaySplineWeights = pd.read_csv("BNBSplineWeights.csv", names=["run", "subrun", "event", "wgt_spline"], dtype={"run" : int, "subrun" : int, "event" : int, "wgt_spline" : float})
dirtSplineWeights    = pd.read_csv("DirtSplineWeights.csv", names=["run", "subrun", "event", "wgt_spline"], dtype={"run" : int, "subrun" : int, "event" : int, "wgt_spline" : float})

#OverlayScale = dataPOT / (sumPOT*(float(maxEvents)/trackDaughters.shape[0]))

overlayCVWeights.insert(overlayCVWeights.shape[1], "DuplicatedEvent", overlayCVWeights.duplicated() ) #Tag the events which are duplicated
overlayCVWeights.replace([np.inf], 0.0, inplace=True)

dirtCVWeights.insert(dirtCVWeights.shape[1], "DuplicatedEvent", dirtCVWeights.duplicated() ) #Tag the events which are duplicated
dirtCVWeights.replace([np.inf], 0.0, inplace=True)

overlaySplineWeights.insert(overlaySplineWeights.shape[1], "DuplicatedEvent", overlaySplineWeights.duplicated() ) #Tag the events which are duplicated
overlaySplineWeights.replace([np.inf], 0.0, inplace=True)

dirtSplineWeights.insert(dirtSplineWeights.shape[1], "DuplicatedEvent", dirtSplineWeights.duplicated() ) #Tag the events which are duplicated
dirtSplineWeights.replace([np.inf], 0.0, inplace=True)

dataDaughters = uproot.open(InputFiles[1])["NuCCanalyzer"]["Daughters"]
trackData     = pd.DataFrame(dataDaughters.arrays(["track_range_mom_mu", "track_mcs_mom", "track_range_mom_p", "track_is_muon_candidate", "track_score", "track_chi2_proton", "track_chi2_muon", "track_dirx", "track_diry", "track_dirz","track_dirx", "track_diry", "track_dirz", "vx", "vy", "vz", "track_endx", "track_endy", "track_endz", "run", "subrun", "event"] ) )
filteredData  = pd.DataFrame(bnbEvents.arrays(["run", "subrun", "event", "nu_mu_cc_selected"]) )

extDaughtersC1 = uproot.open(InputFiles[2])["NuCCanalyzer"]["Daughters"]
trackExtC1     = pd.DataFrame(extDaughtersC1.arrays(["track_range_mom_mu", "track_mcs_mom", "track_range_mom_p", "track_is_muon_candidate", "track_score", "track_chi2_proton", "track_chi2_muon", "track_dirx", "track_diry", "track_dirz", "track_dirx", "track_diry", "track_dirz", "vx", "vy", "vz", "track_endx", "track_endy", "track_endz", "run", "subrun", "event"] ) )
filteredExtC1  = pd.DataFrame(extEventsC1.arrays(["run", "subrun", "event", "nu_mu_cc_selected"]) )

extDaughtersC2 = uproot.open(InputFiles[4])["NuCCanalyzer"]["Daughters"]
trackExtC2     = pd.DataFrame(extDaughtersC2.arrays(["track_range_mom_mu", "track_mcs_mom", "track_range_mom_p", "track_is_muon_candidate", "track_score", "track_chi2_proton", "track_chi2_muon", "track_dirx", "track_diry", "track_dirz", "track_dirx", "track_diry", "track_dirz", "vx", "vy", "vz", "track_endx", "track_endy", "track_endz", "run", "subrun", "event"] ) )
filteredExtC2  = pd.DataFrame(extEventsC2.arrays(["run", "subrun", "event", "nu_mu_cc_selected"]) )

dirtDaughters = uproot.open(InputFiles[3])["NuCCanalyzer"]["Daughters"]
trackDirt     = pd.DataFrame(dirtDaughters.arrays(["track_range_mom_mu", "track_mcs_mom", "track_range_mom_p", "track_is_muon_candidate", "track_score", "track_chi2_proton", "track_chi2_muon", "track_dirx", "track_diry", "track_dirz", "track_dirx", "track_diry", "track_dirz", "vx", "vy", "vz", "track_endx", "track_endy", "track_endz", "run", "subrun", "event"] ) )
filteredDirt  = pd.DataFrame(dirtEvents.arrays(["run", "subrun", "event", "nu_mu_cc_selected"]) )

if(useC2):
  filteredExt = pd.concat([filteredExtC1, filteredExtC2])
  trackExt    = pd.concat([trackExtC1, trackExtC2])
  extTriggers = extTriggers + extTriggersC2

else:
  filteredExt = filteredExtC1
  trackExt    = trackExtC1

ExtScale     = bnbSpills / extTriggers

#Here, we calculate some additional event information that isn't part of the input ROOT ntuple
#This is because the grad. student who created the files didn't include this information
filteredEvents.insert(filteredEvents.shape[1], "DuplicatedEvent", filteredEvents.duplicated() ) #Tag the events which are duplicated
filteredEvents.insert(filteredEvents.shape[1], "mc_channel", [getChan(x, y) for x, y in zip(filteredEvents['mc_nu_interaction_type'], filteredEvents['mc_nu_ccnc'])] ) #Classify neutrino events based on CC / NC and event Type
filteredEvents.eval('mc_Ehad = mc_nu_energy - mc_nu_lepton_energy', inplace=True) #Insert the true energy transfer (nu)
filteredEvents.insert(filteredEvents.shape[1], "mc_expQ2", [getQ2(x, y, z) for x, y, z in zip(filteredEvents['mc_nu_energy'], filteredEvents['mc_nu_lepton_energy'], filteredEvents['mc_nu_lepton_theta'])] )
filteredEvents.insert(filteredEvents.shape[1], "mc_expW", [getW(x, y, z) for x, y, z in zip(filteredEvents['mc_Ehad'], filteredEvents['mc_expQ2'], filteredEvents['mc_channel'] ) ] )
filteredEvents.insert(filteredEvents.shape[1], "mc_expXbj", [getXbj(x, y) for x, y in zip(filteredEvents['mc_Ehad'], filteredEvents['mc_expQ2'] ) ] )
filteredEvents.insert(filteredEvents.shape[1], "mc_expY", [getInel(x, y) for x, y in zip(filteredEvents['mc_Ehad'], filteredEvents['mc_nu_energy'] ) ] )
#filteredEvents.insert(filteredEvents.shape[1], "template_wgt", [getChanWeight(x, y) for x, y in zip(filteredEvents['mc_nu_interaction_type'], filteredEvents['mc_nu_ccnc'])] ) #Classify neutrino events based on CC / NC and event Type
filteredEvents.insert(filteredEvents.shape[1], "pot", mcPOT)
#filteredEvents.eval('combined_wgt = wgt*template_wgt', inplace=True)

trackDaughters.insert(trackDaughters.shape[1], "phi", [getPhi(x, y) for x, y in zip(trackDaughters['track_diry'], trackDaughters['track_dirx'] ) ] )
trackDaughters.insert(trackDaughters.shape[1], "DuplicatedEvent", trackDaughters.duplicated() ) #Tag the events which are duplicated
trackDaughters.insert(trackDaughters.shape[1], "isContained", [isContained(x, y, z, a, b, c) for x, y, z, a, b, c in zip(trackDaughters['vx'], trackDaughters['vy'], trackDaughters['vz'], trackDaughters['track_endx'], trackDaughters['track_endy'], trackDaughters['track_endz'])])
trackDaughters.insert(trackDaughters.shape[1], 'track_mom_best', [getBestP(x, y, z) for x, y, z in zip(trackDaughters['isContained'], trackDaughters['track_range_mom_mu'], trackDaughters['track_mcs_mom'])])
trackDaughters.insert(trackDaughters.shape[1], 'particle', [getParticle(x) for x in trackDaughters['mc_pdg'] ])

filteredDirt.insert(filteredDirt.shape[1], "DuplicatedEvent", filteredDirt.duplicated() ) #Tag the Dirt which are duplicated
'''
filteredDirt.insert(filteredDirt.shape[1], "mc_channel", [getChan(x, False) for x in zip(filteredDirt['mc_nu_interaction_type'], filtere)] ) #Classify neutrino Dirt based on CC / NC and event Type
filteredDirt.eval('mc_Ehad = mc_nu_energy - mc_nu_lepton_energy', inplace=True) #Insert the true energy transfer (nu)
filteredDirt.insert(filteredDirt.shape[1], "mc_expQ2", [getQ2(x, y, z) for x, y, z in zip(filteredDirt['mc_nu_energy'], filteredDirt['mc_nu_lepton_energy'], filteredDirt['mc_nu_lepton_theta'])] )
filteredDirt.insert(filteredDirt.shape[1], "mc_expW", [getW(x, y, z) for x, y, z in zip(filteredDirt['mc_Ehad'], filteredDirt['mc_expQ2'], filteredDirt['mc_channel'] ) ] )
filteredDirt.insert(filteredDirt.shape[1], "mc_expXbj", [getXbj(x, y) for x, y in zip(filteredDirt['mc_Ehad'], filteredDirt['mc_expQ2'] ) ] )
filteredDirt.insert(filteredDirt.shape[1], "mc_expY", [getInel(x, y) for x, y in zip(filteredDirt['mc_Ehad'], filteredDirt['mc_nu_energy'] ) ] )
'''
trackDirt.insert(trackDirt.shape[1], "phi", [getPhi(x, y) for x, y in zip(trackDirt['track_diry'], trackDirt['track_dirx'] ) ] )
trackDirt.insert(trackDirt.shape[1], "DuplicatedEvent", trackDirt.duplicated() ) #Tag the events which are duplicated
trackDirt.insert(trackDirt.shape[1], "isContained", [isContained(x, y, z, a, b, c) for x, y, z, a, b, c in zip(trackDirt['vx'], trackDirt['vy'], trackDirt['vz'], trackDirt['track_endx'], trackDirt['track_endy'], trackDirt['track_endz'])])
trackDirt.insert(trackDirt.shape[1], 'track_mom_best', [getBestP(x, y, z) for x, y, z in zip(trackDirt['isContained'], trackDirt['track_range_mom_mu'], trackDirt['track_mcs_mom'])])
#trackDirt.insert(trackDirt.shape[1], 'particle', [getParticle(x) for x in trackDirt['mc_pdg'] ])

#SAVE THIS
#print trackDaughters.loc[:100,['isContained', 'track_range_mom_mu', 'track_mcs_mom', 'track_mom_best']]

extWeights              = np.full(filteredExt.shape[0],  ExtScale)
extTemplateWeights      = np.full(filteredExt.shape[0],  ExtTemplateWeight)

filteredData.insert(filteredData.shape[1], "DuplicatedEvent", filteredData.duplicated() ) #Tag the events which are duplicated
trackData.insert(trackData.shape[1], "DuplicatedEvent", trackData.duplicated())
trackData.insert(trackData.shape[1], "phi", [getPhi(x, y) for x, y in zip(trackData['track_diry'], trackData['track_dirx'] ) ] )
trackData.insert(trackData.shape[1], "isContained", [isContained(x, y, z, a, b, c) for x, y, z, a, b, c in zip(trackData['vx'], trackData['vy'], trackData['vz'], trackData['track_endx'], trackData['track_endy'], trackData['track_endz'])])
trackData.insert(trackData.shape[1], 'track_mom_best', [getBestP(x, y, z) for x, y, z in zip(trackData['isContained'], trackData['track_range_mom_mu'], trackData['track_mcs_mom'])])
'''
print trackDaughters.loc[:100,['isContained', 'track_range_mom_mu', 'track_mcs_mom', 'track_mom_best']]
raw_input()

print trackData.loc[:100,['isContained', 'track_range_mom_mu', 'track_mcs_mom', 'track_mom_best']]
raw_input()
'''
filteredExt.insert(filteredExt.shape[1], "DuplicatedEvent", filteredExt.duplicated() ) #Tag the events which are duplicated
filteredExt.insert(filteredExt.shape[1], "wgt", extWeights )
filteredExt.insert(filteredExt.shape[1], "template_wgt", extTemplateWeights )
filteredExt.eval('combined_wgt = wgt*template_wgt', inplace=True)

trackExt.insert(trackExt.shape[1], "phi", [getPhi(x, y) for x, y in zip(trackExt['track_diry'], trackExt['track_dirx'] ) ] )
trackExt.insert(trackExt.shape[1], "DuplicatedEvent", trackExt.duplicated())
trackExt.insert(trackExt.shape[1], "isContained", [isContained(x, y, z, a, b, c) for x, y, z, a, b, c in zip(trackExt['vx'], trackExt['vy'], trackExt['vz'], trackExt['track_endx'], trackExt['track_endy'], trackExt['track_endz'])])
trackExt.insert(trackExt.shape[1], 'track_mom_best', [getBestP(x, y, z) for x, y, z in zip(trackExt['isContained'], trackExt['track_range_mom_mu'], trackExt['track_mcs_mom'])])

OverlayScale = dataPOT / sumPOT
DirtScale    = dataPOT / sumDirtPOT
print "MC POT: %e or %e Overlay Scale: %.3f Ext Scale: %.3f" % (sumPOT, sumPOT*(float(maxEvents)/trackDaughters.shape[0]), OverlayScale, ExtScale)
print "Total MC POT: %e total MC events: %d" % (mcPOT.sum(), trackDaughters.shape[0])
print "Total Dirt POT: %e" % sumDirtPOT

#Index the events and daugthers by the run, subrun, event tuple
#This is IMPORTANT. The only infomration we have to connect the two frames a priori is this set of 3 ints
#A single event can have multiple tracks (and often does!)
#Multiindexing makes our life much easier, cuz we can grab the event info for ANY track from it's multiindex

trackDaughters =  trackDaughters.set_index(['run', 'subrun', 'event'])
filteredEvents =  filteredEvents.set_index(['run', 'subrun', 'event'])

filteredData   =  filteredData.set_index(['run', 'subrun', 'event'])
trackData      =  trackData.set_index(['run', 'subrun', 'event'])

filteredExt    = filteredExt.set_index(['run', 'subrun', 'event'])
trackExt       = trackExt.set_index(['run', 'subrun', 'event'])

filteredDirt    = filteredDirt.set_index(['run', 'subrun', 'event'])
trackDirt       = trackDirt.set_index(['run', 'subrun', 'event'])

overlayCVWeights = overlayCVWeights.set_index(['run', 'subrun', 'event'])
dirtCVWeights    = dirtCVWeights.set_index(['run', 'subrun', 'event'])

overlaySplineWeights = overlaySplineWeights.set_index(['run', 'subrun', 'event'])
dirtSplineWeights    = dirtSplineWeights.set_index(['run', 'subrun', 'event'])

#Do this to make our loops and lookups a bit more efficienct

trackDaughters.sort_index()
filteredEvents.sort_index()


filteredEvents   =  filteredEvents[filteredEvents.DuplicatedEvent == False]
filteredData     =  filteredData[filteredData.DuplicatedEvent == False]
filteredExt      =  filteredExt[filteredExt.DuplicatedEvent == False]
filteredDirt     =  filteredDirt[filteredDirt.DuplicatedEvent == False]

trackDaughters   =  trackDaughters[trackDaughters.DuplicatedEvent == False]
trackData        =  trackData[trackData.DuplicatedEvent == False]
trackExt         =  trackExt[trackExt.DuplicatedEvent == False]
trackDirt        =  trackDirt[trackDirt.DuplicatedEvent == False]


overlayCVWeights =  overlayCVWeights[overlayCVWeights.DuplicatedEvent == False]
dirtCVWeights    =  dirtCVWeights[dirtCVWeights.DuplicatedEvent == False]
overlaySplineWeights =  overlaySplineWeights[overlaySplineWeights.DuplicatedEvent == False]
dirtSplineWeights    =  dirtSplineWeights[dirtSplineWeights.DuplicatedEvent == False]


#print filteredEvents.loc[(7001, 920, 46044)]
#raw_input()

#trackDaughters = trackDaughters.iloc[:maxEvents]

#print trackDaughters.loc[(7001, 920, 46044)]['track_is_muon_candidate']
#print filteredEvents.loc[(7001, 920, 46044)]['DuplicatedEvent']

numberFiltered = 0


#create a dict of event info we want to associate with each daughter.
#by doing this, we have the complete event information for each track.
#Really what we want is to look at the particles' properties as a funciton of the underlying event information
#This is extendible to any event varaible we want to associate to a particle
#interactionInfo = {"DuplicatedEvent" : [], "mc_channel" : [], "nu_mu_cc_selected" : [], "mc_Ehad" : [], "mc_expQ2" : [], "mc_expXbj" : [], "mc_expY" : [], "mc_expW" : [] } 
# "pot"
interactionInfo = ("mc_channel", "nu_mu_cc_selected", "mc_Ehad", "mc_expQ2", "mc_expXbj", "mc_expY", "mc_expW") 
#interactionInfo = ("DuplicatedEvent","mc_channel" ) 

for field in interactionInfo:
  trackDaughters   = trackDaughters.join(filteredEvents['%s' % field], on=["run", "subrun", "event"])

trackDaughters = trackDaughters.join(overlayCVWeights['wgt_tune'], on=["run", "subrun", "event"])
trackDaughters = trackDaughters.join(overlaySplineWeights['wgt_spline'], on=["run", "subrun", "event"])

'''
for index, row in trackDaughters.iterrows():
    if(numberFiltered % 10000 == 0):
      print "Filtered: %d" % numberFiltered  
    #This gets around duplicate events. It's a bit clunky
    #If we didn't have duplicate events, we could fairly easily join the daughters / event dataframes based on multindex
    if type(filteredEvents.at[index, "DuplicatedEvent"]) is np.ndarray:
       for itype in interactionInfo:
           interactionInfo[itype].append( filteredEvents.at[index, itype][1] )
       
    else:
        for itype in interactionInfo:
           interactionInfo[itype].append( filteredEvents.at[index, itype] )
    numberFiltered += 1
'''

#print trackDaughters.loc[(7001, 920, 46044)]
#raw_input()

#dataInfo = {"nu_mu_cc_selected" : []}

trackData   = trackData.join(filteredData["nu_mu_cc_selected"], on=["run", "subrun", "event"])
'''
for index, row in trackData.iterrows():
    
    #This gets around duplicate events. It's a bit clunky
    #If we didn't have duplicate events, we could fairly easily join the daughters / event dataframes based on multindex
    if type(filteredData.at[index, "DuplicatedEvent"]) is np.ndarray:
        for itype in dataInfo:
            dataInfo[itype].append( filteredData.at[index, itype][1] )

    else:
        for itype in dataInfo:
            dataInfo[itype].append( filteredData.at[index, itype] )
'''
extInfo = { "nu_mu_cc_selected" : [], "wgt" : [], "template_wgt" : [], "combined_wgt" : [] }

for field in extInfo:
  trackExt   = trackExt.join(filteredExt['%s' % field], on=["run", "subrun", "event"])
'''
for index, row in trackExt.iterrows():
    
    #This gets around duplicate events. It's a bit clunky
     #If we didn't have duplicate events, we could fairly easily join the daughters / event dataframes based on multindex
    if type(filteredExt.at[index, "DuplicatedEvent"]) is np.ndarray:
        for itype in extInfo:
            extInfo[itype].append( filteredExt.at[index, itype][1] )

    else:
        for itype in extInfo:
            extInfo[itype].append( filteredExt.at[index, itype] )

'''
trackDirt   = trackDirt.join(filteredDirt["nu_mu_cc_selected"], on=["run", "subrun", "event"])

trackDirt = trackDirt.join(dirtCVWeights['wgt_tune'], on=["run", "subrun", "event"])
trackDirt = trackDirt.join(dirtSplineWeights['wgt_spline'], on=["run", "subrun", "event"])



muonMomentumRange   = (0.0, 2.0)
protonMomentumRange = (0.0, 1.5)
phiRange = (-1.1, 1.1)
phiDiffRange = (0.0, 2.0)
cosineRange  = (-1.0, 1.0)
chi2Range    = (0.0, 400.0)


#associate all the event info with the particles (at last!)
'''
for itype in interactionInfo:
   trackDaughters.insert(trackDaughters.shape[1], itype, interactionInfo[itype])
'''
#sumPOT = trackDaughters['pot'].sum()
'''
print trackDaughters['pot']
print trackDaughters.groupby(level=["run", "subrun", "event"]).agg({"pot" : ["min"]})
sumPOT = trackDaughters.groupby(level=["run", "subrun", "event"]).agg({"pot" : ["min"]}).sum()
'''


overlayWeights = np.full(trackDaughters.shape[0], OverlayScale )
dirtWeights    = np.full(trackDirt.shape[0], DirtScale )


trackDaughters.insert(trackDaughters.shape[1], "pot_wgt", overlayWeights )
trackDaughters.eval('wgt = pot_wgt*wgt_tune*wgt_spline', inplace=True) 

trackDirt.insert(trackDirt.shape[1], "pot_wgt", dirtWeights )
trackDirt.eval('wgt = pot_wgt*wgt_tune*wgt_spline', inplace=True) 


#for itype in dataInfo:
#   trackData.insert(trackData.shape[1], itype, dataInfo[itype])
'''
for itype in extInfo:
   trackExt.insert(trackExt.shape[1], itype, extInfo[itype])
'''
extMuons     = trackExt.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == True & track_score > @minTrackScore')
dirtMuons    = trackDirt.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == True & track_score > @minTrackScore')
muonTracks   = trackDaughters.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == True & track_score > @minTrackScore')

#Sanity check: Do we see any primary tracks which fail the min score cut?
print "Selected MC events failing track score cut: %d" % len(trackDaughters.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == True & track_score < @minTrackScore').index)
print "Selected BNB events failing track score cut: %d" % len(trackData.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == True & track_score < @minTrackScore').index)

protonTracks = trackDaughters.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == False & track_chi2_proton < @maxProtonChi2 & track_score > @minTrackScore')
protonCandidateTracks = trackDaughters.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_score > @minTrackScore')
extProtonCandidates = trackExt.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_score > @minTrackScore')
#inclusiveProtonCandidates = [protonCandidateTracks.query('mc_channel == "QE"')['track_chi2_proton'].to_numpy(), protonCandidateTracks.query('mc_channel == "RES"')['track_chi2_proton'].to_numpy(), protonCandidateTracks.query('mc_channel == "DIS"')['track_chi2_proton'].to_numpy(), protonCandidateTracks.query('mc_channel == "2p2h"')['track_chi2_proton'].to_numpy(), protonCandidateTracks.query('mc_channel == "NC / Other"')['track_chi2_proton'].to_numpy(), extProtonCandidates['track_chi2_proton'].to_numpy()]
inclusiveProtonCandidates = [protonCandidateTracks.query('particle == "muon"')['track_chi2_proton'].to_numpy(), protonCandidateTracks.query('particle == "proton"')['track_chi2_proton'].to_numpy(), protonCandidateTracks.query('particle == "pion"')['track_chi2_proton'].to_numpy(), protonCandidateTracks.query('particle == "electron"')['track_chi2_proton'].to_numpy(), protonCandidateTracks.query('particle == "other"')['track_chi2_proton'].to_numpy(), extProtonCandidates['track_chi2_proton'].to_numpy()]
#protonCandidateWeights = [protonCandidateTracks.query('mc_channel == "QE"')['wgt'].to_numpy(), protonCandidateTracks.query('mc_channel == "RES"')['wgt'].to_numpy(), protonCandidateTracks.query('mc_channel == "DIS"')['wgt'].to_numpy(), protonCandidateTracks.query('mc_channel == "2p2h"')['wgt'].to_numpy(), protonCandidateTracks.query('mc_channel == "NC / Other"')['wgt'].to_numpy(), extProtonCandidates['wgt'].to_numpy()]
protonCandidateWeights = [protonCandidateTracks.query('particle == "muon"')['wgt'].to_numpy(), protonCandidateTracks.query('particle == "proton"')['wgt'].to_numpy(), protonCandidateTracks.query('particle == "pion"')['wgt'].to_numpy(), protonCandidateTracks.query('particle == "electron"')['wgt'].to_numpy(), protonCandidateTracks.query('particle == "other"')['wgt'].to_numpy(), extProtonCandidates['wgt'].to_numpy()]

 

leadingProtons = protonTracks.groupby(level=["run", "subrun", "event"]).agg({"track_range_mom_p" : ["max", "count"]})
leadingMuons   = muonTracks.groupby(level=["run", "subrun", "event"]).agg({"track_mcs_mom" : ["max", "count"]})
# Using ravel, and a string join, we can create better names for the columns:
leadingProtons.columns = ["_".join(x) for x in leadingProtons.columns.ravel()]
leadingMuons.columns = ["_".join(x) for x in leadingMuons.columns.ravel()]
#print leadingProtons.head()
#leadingProtons.eval('isLeadingP = (track_range_mom_p == track_range_mom_p_max)', inplace=True)
'''
for x in leadingProtons.columns.ravel():
   print x[1]
'''
#FOR THE FUTURE, ie pandas 0.25
#leadingProtons = protonTracks.groupby(level=["run", "subrun", "event"]).agg("leading_proton_p"=pd.NamedAgg('column = track_range_mom_p', aggfunc='max'), "number_of_protons"=pd.NamedAgg('column = track_range_mom_p', aggfunc='count') )
#leadingProtons
'''
for index, rows in leadingProtons.iterrows():
    #print rows
    print index
'''

protonTracks = protonTracks.join(leadingProtons, on=["run", "subrun", "event"])
protonTracks.eval('isLeadingP = (track_range_mom_p == track_range_mom_p_max)', inplace=True)

muonTracks = muonTracks.join(leadingMuons, on=["run", "subrun", "event"])
muonTracks.eval('isLeadingMu = (track_mcs_mom == track_mcs_mom_max)', inplace=True)

#print muonTracks.query('mc_channel == "QE" & track_mcs_mom_count > 2')['isLeadingMu']
#print muonTracks.loc[(7001, 920, 46044)]['track_chi2_muon']
#raw_input()


muonTracks   = muonTracks.join(leadingProtons, on=["run", "subrun", "event"])
protonTracks   = protonTracks.join(leadingMuons, on=["run", "subrun", "event"])

extProtons   = trackExt.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == False & track_chi2_proton < @maxProtonChi2 & track_score > @minTrackScore')
leadingExtProtons = extProtons.groupby(level=["run", "subrun", "event"]).agg({"track_range_mom_p" : ["max", "count"]})
leadingExtMuons   = extMuons.groupby(level=["run", "subrun", "event"]).agg({"track_mcs_mom" : ["max", "count"]})
# Using ravel, and a string join, we can create better names for the columns:
leadingExtProtons.columns = ["_".join(x) for x in leadingExtProtons.columns.ravel()]
leadingExtMuons.columns = ["_".join(x) for x in leadingExtMuons.columns.ravel()]

extProtons = extProtons.join(leadingExtProtons, on=["run", "subrun", "event"])
extProtons = extProtons.join(leadingExtMuons, on=["run", "subrun", "event"])
extProtons.eval('isLeadingP = (track_range_mom_p == track_range_mom_p_max)', inplace=True)

extMuons   = extMuons.join(leadingExtMuons, on=["run", "subrun", "event"])
extMuons   = extMuons.join(leadingExtProtons, on=["run", "subrun", "event"])
extMuons.eval('isLeadingMu = (track_mcs_mom == track_mcs_mom_max)', inplace=True)

dirtProtons   = trackDirt.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == False & track_chi2_proton < @maxProtonChi2 & track_score > @minTrackScore')
leadingDirtProtons = dirtProtons.groupby(level=["run", "subrun", "event"]).agg({"track_range_mom_p" : ["max", "count"]})
leadingDirtMuons   = dirtMuons.groupby(level=["run", "subrun", "event"]).agg({"track_mcs_mom" : ["max", "count"]})
# Using ravel, and a string join, we can create better names for the columns:
leadingDirtProtons.columns = ["_".join(x) for x in leadingDirtProtons.columns.ravel()]
leadingDirtMuons.columns = ["_".join(x) for x in leadingDirtMuons.columns.ravel()]

dirtProtons = dirtProtons.join(leadingDirtProtons, on=["run", "subrun", "event"])
dirtProtons = dirtProtons.join(leadingDirtMuons, on=["run", "subrun", "event"])
dirtProtons.eval('isLeadingP = (track_range_mom_p == track_range_mom_p_max)', inplace=True)

dirtMuons   = dirtMuons.join(leadingDirtMuons, on=["run", "subrun", "event"])
dirtMuons   = dirtMuons.join(leadingDirtProtons, on=["run", "subrun", "event"])
dirtMuons.eval('isLeadingMu = (track_mcs_mom == track_mcs_mom_max)', inplace=True)

#Define the general queries (n P, 1 P, etc.) here. Then, use regex to query any variable I choose
inclusiveMuons = '''[muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "QE"')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "RES"')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "DIS"')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "2p2h"')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "NC / Other"')['VAR'].to_numpy(), dirtMuons.query('track_mcs_mom_count == 1')['VAR'].to_numpy(), extMuons.query('track_mcs_mom_count == 1')['VAR'].to_numpy()]'''

inclusiveProtons = '''protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "QE"')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "RES"')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "DIS"')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "2p2h"')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "NC / Other"')['VAR'].to_numpy(), dirtProtons.query('track_mcs_mom_count == 1')['VAR'].to_numpy(), extProtons.query('track_mcs_mom_count == 1')['VAR'].to_numpy()'''

onePMuons = '''[muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "QE" & track_range_mom_p_count == 1')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "RES" & track_range_mom_p_count == 1')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "DIS" & track_range_mom_p_count == 1')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "2p2h" & track_range_mom_p_count == 1')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "NC / Other" & track_range_mom_p_count == 1')['VAR'].to_numpy(), dirtMuons.query('track_mcs_mom_count == 1 & track_range_mom_p_count == 1')['VAR'].to_numpy(), extMuons.query('track_mcs_mom_count == 1 & track_range_mom_p_count == 1')['VAR'].to_numpy()]'''

onePProtons = '''[protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "QE" & track_range_mom_p_count == 1')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "RES" & track_range_mom_p_count == 1')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "DIS" & track_range_mom_p_count == 1')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "2p2h" & track_range_mom_p_count == 1')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "NC / Other" & track_range_mom_p_count == 1')['VAR'].to_numpy(), dirtProtons.query('track_mcs_mom_count == 1 & track_range_mom_p_count == 1')['VAR'].to_numpy(), extProtons.query('track_mcs_mom_count == 1 & track_range_mom_p_count == 1')['VAR'].to_numpy()]'''

nPMuons = '''[muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "QE" & track_range_mom_p_count > 0')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "RES" & track_range_mom_p_count > 0')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "DIS" & track_range_mom_p_count > 0')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "2p2h" & track_range_mom_p_count > 0')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "NC / Other" & track_range_mom_p_count > 0')['VAR'].to_numpy(), dirtMuons.query('track_mcs_mom_count == 1 & track_range_mom_p_count > 0')['VAR'].to_numpy(), extMuons.query('track_mcs_mom_count == 1 & track_range_mom_p_count > 0')['VAR'].to_numpy()]'''

nPProtons = '''[protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "QE" & track_range_mom_p_count> 0 & isLeadingP == True')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "RES" & track_range_mom_p_count > 0 & isLeadingP == True')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "DIS" & track_range_mom_p_count > 0 & isLeadingP == True')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "2p2h" & track_range_mom_p_count > 0 & isLeadingP == True')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "NC / Other" & track_range_mom_p_count > 0 & isLeadingP == True')['VAR'].to_numpy(), dirtProtons.query('track_mcs_mom_count == 1 & track_range_mom_p_count > 0 & isLeadingP == True')['VAR'].to_numpy(), extProtons.query('track_mcs_mom_count == 1 & track_range_mom_p_count > 0 & isLeadingP == True')['VAR'].to_numpy()]'''

onePMuons = '''[muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "QE" & track_range_mom_p_count == 1')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "RES" & track_range_mom_p_count == 1')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "DIS" & track_range_mom_p_count == 1')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "2p2h" & track_range_mom_p_count == 1')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "NC / Other" & track_range_mom_p_count == 1')['VAR'].to_numpy(), dirtMuons.query('track_mcs_mom_count == 1 & track_range_mom_p_count == 1')['VAR'].to_numpy(), extMuons.query('track_mcs_mom_count == 1 & track_range_mom_p_count == 1')['VAR'].to_numpy()]'''

onePProtons = '''[protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "QE" & track_range_mom_p_count == 1')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "RES" & track_range_mom_p_count == 1')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "DIS" & track_range_mom_p_count == 1')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "2p2h" & track_range_mom_p_count == 1')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "NC / Other" & track_range_mom_p_count == 1')['VAR'].to_numpy(),
  dirtProtons.query('track_mcs_mom_count == 1 & track_range_mom_p_count == 1')['VAR'].to_numpy(),
  extProtons.query('track_mcs_mom_count == 1 & track_range_mom_p_count == 1')['VAR'].to_numpy()]'''

twoPMuons = '''[muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "QE" & track_range_mom_p_count == 2')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "RES" & track_range_mom_p_count == 2')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "DIS" & track_range_mom_p_count == 2')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "2p2h" & track_range_mom_p_count == 2')['VAR'].to_numpy(), muonTracks.query('track_mcs_mom_count == 1 & mc_channel == "NC / Other" & track_range_mom_p_count == 2')['VAR'].to_numpy(), dirtMuons.query('track_mcs_mom_count == 1 & track_range_mom_p_count == 2')['VAR'].to_numpy(), extMuons.query('track_mcs_mom_count == 1 & track_range_mom_p_count == 2')['VAR'].to_numpy()]'''

twoPProtons = '''[protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "QE" & track_range_mom_p_count == 2 & isLeadingP == True')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "RES" & track_range_mom_p_count == 2 & isLeadingP == True')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "DIS" & track_range_mom_p_count == 2 & isLeadingP == True')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "2p2h" & track_range_mom_p_count == 2 & isLeadingP == True')['VAR'].to_numpy(), protonTracks.query('track_mcs_mom_count == 1 & mc_channel == "NC / Other" & track_range_mom_p_count == 2 & isLeadingP == True')['VAR'].to_numpy(), dirtProtons.query('track_mcs_mom_count == 1 & track_range_mom_p_count == 2 & isLeadingP == True')['VAR'].to_numpy(), extProtons.query('track_mcs_mom_count == 1 & track_range_mom_p_count == 2 & isLeadingP == True')['VAR'].to_numpy()]'''

exec( "incMomentumStack   = "  + re.sub(r'VAR', 'track_mom_best', inclusiveMuons) )
exec( "incPhiStack        = "  + re.sub(r'VAR', 'phi', inclusiveMuons) )
exec( "incMuonWeights     = "  + re.sub(r'VAR', 'wgt', inclusiveMuons) )
exec( "incMuonTheta      =  "  + re.sub(r'VAR', 'track_dirz', inclusiveMuons)) 

exec( "nPMomentumStack    = "  + re.sub(r'VAR', 'track_mom_best', nPMuons) )
exec( "nPPhiStack         = "  + re.sub(r'VAR', 'phi', nPMuons) )
exec( "nPMuonWeights      = "  + re.sub(r'VAR', 'wgt', nPMuons) )
exec( "nPMuonTheta        = "  + re.sub(r'VAR', 'track_dirz', nPMuons)) 

exec( "onePMomentumStack  = "  + re.sub(r'VAR', 'track_mom_best', onePMuons) )
exec( "onePPhiStack       = "  + re.sub(r'VAR', 'phi', onePMuons) )
exec( "onePMuonWeights    = "  + re.sub(r'VAR', 'wgt', onePMuons) )
exec( "onePMuonTheta      = "  + re.sub(r'VAR', 'track_dirz', onePMuons) )

exec( "twoPMomentumStack  = "  + re.sub(r'VAR', 'track_mom_best', twoPMuons) )
exec( "twoPPhiStack       = "  + re.sub(r'VAR', 'phi', twoPMuons) )
exec( "twoPMuonWeights    = "  + re.sub(r'VAR', 'wgt', twoPMuons) )
exec( "twoPMuonTheta      = "  + re.sub(r'VAR', 'track_dirz', twoPMuons) )

exec( "incProtonStack     = "  + re.sub(r'VAR', 'track_range_mom_p', inclusiveProtons) )
exec( "incProtonPhiStack  = "  + re.sub(r'VAR', 'phi', inclusiveProtons) )
exec( "incProtonWeights   = "  + re.sub(r'VAR', 'wgt', inclusiveProtons) )
exec( "incProtonTheta     = "  + re.sub(r'VAR', 'track_dirz', inclusiveProtons) )

exec( "nPProtonStack      = "  + re.sub(r'VAR', 'track_range_mom_p_max', nPProtons) )
exec( "nPPhiProtonStack   = "  + re.sub(r'VAR', 'phi', nPProtons) )
exec( "nPProtonWeights    = "  + re.sub(r'VAR', 'wgt', nPProtons) )
exec( "nPProtonTheta      = "  + re.sub(r'VAR', 'track_dirz', nPProtons) )

exec( "onePProtonStack    = "  + re.sub(r'VAR', 'track_range_mom_p', onePProtons) )
exec( "onePPhiProtonStack = "  + re.sub(r'VAR', 'phi', onePProtons) )
exec( "onePProtonWeights  = "  + re.sub(r'VAR', 'wgt', onePProtons) )
exec( "onePProtonTheta    = "  + re.sub(r'VAR', 'track_dirz', onePProtons) )

exec( "twoPProtonStack    = "  + re.sub(r'VAR', 'track_range_mom_p_max', twoPProtons) )
exec( "twoPPhiProtonStack = "  + re.sub(r'VAR', 'phi', twoPProtons) )
exec( "twoPProtonWeights  = "  + re.sub(r'VAR', 'wgt', twoPProtons) )
exec( "twoPProtonTheta    = "  + re.sub(r'VAR', 'track_dirz', twoPProtons) )

#Variables for fitting

MCType = ("QE", "RES", "DIS", "2p2h", "NC / Other")
muonProtonTheta = []
muonProtonThetaWeights = []
muonProtonPhiDiff      = []
#print protonTracks.query('track_mcs_mom_count ==1 & isLeadingP==True & mc_channel == "QE"')
#print muonTracks.query('track_mcs_mom_count ==1 & track_range_mom_p_count > 0 & mc_channel == "QE"')
for chan in MCType:
  #muonProtonTheta.append(dotProduct(protonTracks.query('track_mcs_mom_count ==1 & track_range_mom_p_count == 1 & mc_channel == @chan'), muonTracks.query('track_mcs_mom_count ==1 & track_range_mom_p_count == 1 & mc_channel == @chan'), "track_mcs_mom", "track_range_mom_p", ["track_dirx", "track_diry", "track_dirz"]))
  muonProtonTheta.append(dotProduct(protonTracks.query('track_mcs_mom_count ==1 & isLeadingP==True & mc_channel == @chan'), muonTracks.query('track_mcs_mom_count ==1 & track_range_mom_p_count > 0 & mc_channel == @chan'), "track_mcs_mom", "track_range_mom_p", ["track_dirx", "track_diry", "track_dirz"]))

  muonProtonThetaWeights.append(protonTracks.query('track_mcs_mom_count ==1 & isLeadingP==True & mc_channel == @chan')['wgt'].to_numpy())
  muonProtonPhiDiff.append(np.absolute(muonTracks.query('track_mcs_mom_count ==1 & track_range_mom_p_count > 0 &  mc_channel == @chan')['phi'].to_numpy() - protonTracks.query('track_mcs_mom_count ==1  & isLeadingP==True & mc_channel == @chan')['phi'].to_numpy()))

muonProtonTheta.append(dotProduct(extProtons.query('track_mcs_mom_count ==1 & isLeadingP==True'), extMuons.query('track_mcs_mom_count ==1 & track_range_mom_p_count > 0'), "track_mcs_mom", "track_range_mom_p", ["track_dirx", "track_diry", "track_dirz"]))
muonProtonThetaWeights.append(extProtons.query('track_mcs_mom_count ==1 & isLeadingP==True')['wgt'].to_numpy())
muonProtonPhiDiff.append(extMuons.query('track_mcs_mom_count == 1 & track_range_mom_p_count > 0')['phi'].to_numpy() - extProtons.query('track_mcs_mom_count == 1 & isLeadingP==True')['phi'].to_numpy())

#Create the 5D MC templates here
mcTemplatesDict = {"QE" : [], "RES" : [], "DIS" : [], "2p2h" : [], "NC / Other" : [], "Ext" : []}

print "Size of phi diff: %d" % len(muonProtonPhiDiff[0])
print "Size of proton stack: %d" % len(nPProtonStack[0])
print "Size of Muon stack: %d" % len(nPMomentumStack[0]) 

dataMuons   = trackData.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == True & track_score > @minTrackScore')
dataProtons = trackData.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == False & track_chi2_proton < @maxProtonChi2 & track_score > @minTrackScore')
dataProtonCandidates = trackData.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_score > @minTrackScore')

leadingDataProtons = dataProtons.groupby(level=["run", "subrun", "event"]).agg({"track_range_mom_p" : ["max", "count"]})
leadingDataMuons   = dataMuons.groupby(level=["run", "subrun", "event"]).agg({"track_mcs_mom" : ["max", "count"]})

# Using ravel, and a string join, we can create better names for the columns:
leadingDataProtons.columns = ["_".join(x) for x in leadingDataProtons.columns.ravel()]
leadingDataMuons.columns = ["_".join(x) for x in leadingDataMuons.columns.ravel()]

dataProtons = dataProtons.join(leadingDataProtons, on=["run", "subrun", "event"])
dataProtons   = dataProtons.join(leadingDataMuons, on=["run", "subrun", "event"])
dataProtons.eval('isLeadingP = (track_range_mom_p == track_range_mom_p_max)', inplace=True)

dataMuons   = dataMuons.join(leadingDataMuons, on=["run", "subrun", "event"])
dataMuons   = dataMuons.join(leadingDataProtons, on=["run", "subrun", "event"])
dataMuons.eval('isLeadingMu = (track_mcs_mom == track_mcs_mom_max)', inplace=True)

dataMuonPhi   = trackData.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == True & track_score > @minTrackScore')
dataProtonPhi = trackData.query('DuplicatedEvent == False & nu_mu_cc_selected == True & track_is_muon_candidate == False & track_chi2_proton < @maxProtonChi2 & track_score > @minTrackScore')

File = ROOT.TFile.Open("TemplatesQEPlus20.root", "RECREATE")
#multibinData    = [dataMuons.query('track_range_mom_p_count > 0 & track_mcs_mom_count ==1')['track_mom_best'].to_numpy(), dataProtons.query('track_mcs_mom_count ==1 & isLeadingP == True')['track_range_mom_p_max'].to_numpy(), dataMuons.query('track_mcs_mom_count ==1 & track_range_mom_p_count > 0 ')['track_dirz'].to_numpy(), dataProtons.query('track_range_mom_p_count == 1 & track_mcs_mom_count ==1')['track_dirz'].to_numpy(), np.absolute(dataMuons.query('track_mcs_mom_count ==1 & track_range_mom_p_count > 0 ')['phi'].to_numpy() - dataProtons.query('track_mcs_mom_count ==1 & isLeadingP==True')['phi'].to_numpy())]
numDataEvents   = len(dataMuons.query('track_range_mom_p_count > 0 & track_mcs_mom_count ==1').index)
#print "Number of data Events: %d" % numDataEvents
#print "Number of proton data Events: %d" % len(np.full(numDataEvents, 1.0))
#This is a dummy, we weight the data to 1.0
T_mu          = dataMuons.query('track_range_mom_p_count > 0 & track_mcs_mom_count ==1')['track_mom_best'].to_numpy()
T_p           = dataProtons.query('track_mcs_mom_count ==1 & isLeadingP == True &  track_range_mom_p_count > 0')['track_range_mom_p_max'].to_numpy()
theta_mu      = dataMuons.query('track_mcs_mom_count ==1 & track_range_mom_p_count > 0 ')['track_dirz'].to_numpy()
theta_p       = dataProtons.query('track_mcs_mom_count ==1 & isLeadingP &  track_range_mom_p_count > 0')['track_dirz'].to_numpy()
phi_mu_p_diff = np.absolute(dataMuons.query('track_mcs_mom_count ==1 & track_range_mom_p_count > 0 ')['phi'].to_numpy() - dataProtons.query('track_mcs_mom_count ==1 & isLeadingP==True  & track_range_mom_p_count > 0')['phi'].to_numpy() )

multibinData = [T_mu, T_p, theta_mu, theta_p, phi_mu_p_diff]

weights = np.ones(numDataEvents, dtype=float)
multibinWeights = [weights, weights, weights, weights, weights]
projection = 0
#multibinWeights = [np.full(numDataEvents, 1.0, dtype=float), np.full(numDataEvents, 1.0, dtype=float), np.full(numDataEvents, 1.0, dtype=float), np.full(numDataEvents, 1.0, dtype=float), np.full(numDataEvents, 1.0, dtype=float)]
#print dataProtons.query('track_mcs_mom_count ==1 & isLeadingP == True')['track_range_mom_p_max']
print "Data array length: %d" % len(multibinData)
#raw_input()
#print multibinData[1]
#print np.absolute(dataMuons.query('track_mcs_mom_count ==1 & track_range_mom_p_count > 0 ')['phi'].to_numpy() - dataProtons.query('track_mcs_mom_count ==1 & isLeadingP==True')['phi'].to_numpy())
'''
stat, bins, counts = getMultibin(multibinData, multibinWeights)
dataHisto = stat[projection,].flatten()
dataHist = ROOT.TH1F("Data", "Data", len(dataHisto), 0, len(dataHisto))
for i in range(len(dataHisto)):
   dataHist.SetBinContent(i+1, dataHisto[i])
dataHist.Write()


#Add one for the ext
for i in range(len(MCType)+1):
  #[T_mu, T_p, theta_mu, theta_p, phi_mu_p_diff]
  
  data    = [nPMomentumStack[i], nPProtonStack[i], nPMuonTheta[i], nPProtonTheta[i], muonProtonPhiDiff[i]]
  weights = [nPMuonWeights[i], nPProtonWeights[i], nPMuonWeights[i], nPProtonWeights[i], muonProtonThetaWeights[i]]
  stat, bins, counts = getMultibin(data, weights)
#  print stat[projection,]
#  raw_input()
  histo = stat[projection,].flatten()
  templateHist    = ROOT.TH1F(("template%d" % i), ("Template %d" % i), len(histo), 0, len(histo))
  templateWeights = ROOT.TH1F(("templateWeight%d" % i), ("Template %d" % i),len(histo), 0, len(histo))
  weight = 1.0
  if (i == 0):
    weight = 1.2

  
  for j in range(len(histo)):
     templateHist.SetBinContent(j+1, histo[j])
     #print nPProtonWeights[i][0]
     templateWeights.SetBinContent(j+1,  weight)
  templateHist.Write()
  templateWeights.Write()     

File.Close()
'''
nProtonMomentumBins = 25
dir_name = "PlotDir"

figNumber = 1
fig = plt.figure(figNumber)

figNumber = figNumber + 1

#fig = plt.figure(figNumber)
#makeDataMCHistogram(nPMomentumStack, nPMuonWeights, dataMuons.query('track_range_mom_p_count > 0 & track_mcs_mom_count ==1')['track_mom_best'].to_numpy(), muonMomentumRange, nProtonMomentumBins, "nMuonP", ["Muon Momentum: (1 mu n p)", "Momentum (GeV/c)", "Number of Primary Muons"])
#figNumber = figNumber + 1

fig = plt.figure(figNumber)
#makeDataMCHistogram(onePMomentumStack, onePMuonWeights,dataMuons.query('track_range_mom_p_count == 1 & track_mcs_mom_count ==1')['track_mom_best'].to_numpy(), muonMomentumRange, 18, "OnePMuonP", ["Muon Momentum: (1 mu 1 p)", "Momentum (GeV/c)", "Number of Primary Muons"])
figNumber = figNumber + 1

#print splitAndSort(nPMuonTheta[0], 4)
#print splitAndSort(nPMuonTheta[1], 4)

disMuMomentum     = splitAndSort(nPMomentumStack[0], 4)
disProtonMomentum = splitAndSort(nPProtonStack[0], 4)

print "Muon momentum bins: %.2f, %.2f, %.2f, %.2f" % (disMuMomentum[0][0], disMuMomentum[1][0], disMuMomentum[2][0], disMuMomentum[3][0])
print "Proton momentum bins: %.2f, %.2f, %.2f, %.2f" % (disProtonMomentum[0][0], disProtonMomentum[1][0], disProtonMomentum[2][0], disProtonMomentum[3][0])


fig = plt.figure(figNumber)
for i in range(len(MCType)):
  channel = MCType[i]

  title = "nPMuonPMCOnly"
  fig = plt.figure(figNumber)
  makeMCHistogram(nPMomentumStack[i], channel, muonMomentumRange, 100, title, ["Muon Momentum: (1 mu n p)", "Momentum (GeV/c)", "Number of Primary Muons"])
  figNumber = figNumber + 1

  title = "nPProtonPMCOnly"
  fig = plt.figure(figNumber)
  makeMCHistogram(nPProtonStack[i], channel, protonMomentumRange, 100, title, ["Proton Momentum: (1 mu n p)", "Momentum (GeV/c)", "Number of Leading Protons"])
  figNumber = figNumber + 1

  title = "nPMuonPvsProtonPMCOnly"
  fig = plt.figure(figNumber)
  make2DMCHistogram([nPMomentumStack[i], nPProtonStack[i]], channel, [muonMomentumRange, protonMomentumRange], [100, 100], title, ["Proton vs. Muon Momentum  (1 mu n p)", "Muon Momentum: (GeV/c)", "Proton Momentum: (GeV/c)"])
  figNumber = figNumber + 1

  title = "nPMuonThetaMCOnly"
  fig = plt.figure(figNumber)
  makeMCHistogram(nPMuonTheta[i], channel, phiRange, 100, title, ["Muon Angle: (1 mu n p)", "Cosine (theta)", "Number of Primary Muons"])
  figNumber = figNumber + 1

  title = "nPProtonThetaMCOnly"
  fig = plt.figure(figNumber)
  makeMCHistogram(nPProtonTheta[i], channel, phiRange, 100, title, ["Proton Angle: (1 mu n p)", "Cosine (theta)", "Number of Primary Protons"])
  figNumber = figNumber + 1

  title = "nPMuonThetavsProtonThetaMCOnly"
  fig = plt.figure(figNumber)
  make2DMCHistogram([nPMuonTheta[i], nPProtonTheta[i]], channel, [phiRange, phiRange], [100, 100], title, ["Proton vs. Muon Angle (1 mu n p)", "Muon Angle: cos(theta)", "Proton Angle: cos(theta)"])
  figNumber = figNumber + 1

  fig = plt.figure(figNumber)
  title = "OnePMuPDeltaPhiMCOnly"
  makeMCHistogram(muonProtonPhiDiff[i], channel, phiDiffRange, 100, title, ["Muon-Proton Phi Difference (1 mu N p)", "Angle / pi (radians)", "Number of Muon-Proton Pairs"])
  figNumber = figNumber + 1


#plt.hist(muonProtonPhiDiff, bins=25, stacked=True, range=phiDiffRange, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights = muonProtonThetaWeights )
'''
plt.title("Muon-Proton Phi Difference (1 mu N p)")
plt.xlabel("Angle / pi (radians)")
plt.ylabel("Number of Muon-Proton Pairs")
'''


figNumber = figNumber + 1


#makeDataMCHistogram(mcList, mcWeights, dataList, binRange, nBins, filename, Titles)
makeDataMCHistogram(incPhiStack, incMuonWeights, dataMuons.query('track_mcs_mom_count ==1')['phi'].to_numpy(), phiRange, 50, "IncMuonPhi", ["Muon Phi Angle", "Angle / pi (radians)", "Number of Primary Muons"])

makeDataMCRatioHistogram(incPhiStack, incMuonWeights, dataMuons.query('track_mcs_mom_count ==1')['phi'].to_numpy(), phiRange, 50, "IncMuonPhi", ["Muon Phi Angle", "Angle / pi (radians)", "Number of Primary Muons"])

makeDataMCHistogram(incMomentumStack, incMuonWeights, dataMuons.query('track_mcs_mom_count ==1')['track_mom_best'].to_numpy(), muonMomentumRange, 50, "IncMuonP", ["Muon Momentum", "Momentum (GeV/c)", "Number of Primary Muons"])

makeDataMCHistogram(incProtonStack, incProtonWeights,  dataProtons.query('track_mcs_mom_count ==1')['track_range_mom_p'].to_numpy(), protonMomentumRange, 50, "IncProtonP", ["Proton Range Momentum", "Momentum (GeV/c)", "Number of Protons"])

'''
fig = plt.figure(figNumber)
plt.hist(twoPMomentumStack, bins=18, stacked=True, range=muonMomentumRange, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights = twoPMuonWeights )
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Ext'])
plt.title("Muon Momentum: (1 mu 2 p)")
plt.xlabel("Momentum (GeV/c)")
plt.ylabel("Number of Primary Muons")
data_hist = dataify(dataMuons.query('track_range_mom_p_count == 2 & track_mcs_mom_count ==1')['track_mom_best'].to_numpy(), 18, muonMomentumRange)
plt.errorbar(data_hist[0], data_hist[1], yerr=data_hist[2], fmt='o', color='black')
plt.xlim(0.0, 2.0)
plt.savefig(("%s/TwoPMuonP.png" % dir_name))
figNumber = figNumber + 1

fig = plt.figure(figNumber)
plt.hist(incProtonStack, bins=50, stacked=True, range=protonMomentumRange, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights = incProtonWeights )
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Ext'])
plt.title("Proton Range Momentum")
plt.xlabel("Momentum (GeV/c)")
plt.ylabel("Number of Protons")
data_proton_hist = dataify(dataProtons.query('track_mcs_mom_count ==1')['track_range_mom_p'].to_numpy(), 50, protonMomentumRange)
plt.errorbar(data_proton_hist[0], data_proton_hist[1], yerr=data_proton_hist[2], fmt='o', color='black')
plt.savefig(("%s/IncProtonP.png" % dir_name))
figNumber = figNumber + 1

fig = plt.figure(figNumber)
plt.hist(nPProtonStack, bins=nProtonMomentumBins, stacked=True, range=protonMomentumRange, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights = nPProtonWeights )
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Ext'])
plt.title("Leading Proton Range Momentum (1 mu n p)")
plt.xlabel("Momentum (GeV/c)")
plt.ylabel("Number of Protons")
data_proton_hist = dataify(dataProtons.query('track_mcs_mom_count ==1 & isLeadingP == True')['track_range_mom_p_max'].to_numpy(), nProtonMomentumBins, protonMomentumRange)
plt.errorbar(data_proton_hist[0], data_proton_hist[1], yerr=data_proton_hist[2], fmt='o', color='black')
plt.savefig(("%s/nPLeadingProtonP.png" % dir_name))
figNumber = figNumber + 1

#print crossSectionStack


fig = plt.figure(figNumber)
plt.hist(onePProtonStack, bins=50, stacked=True, range=protonMomentumRange, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights = onePProtonWeights )
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Ext'])
plt.title("Proton Range Momentum (1 mu 1 p)")
plt.xlabel("Momentum (GeV/c)")
plt.ylabel("Number of Protons")
data_proton_hist = dataify(dataProtons.query('track_range_mom_p_count == 1 & track_mcs_mom_count ==1')['track_range_mom_p'].to_numpy(), 50, protonMomentumRange)
plt.errorbar(data_proton_hist[0], data_proton_hist[1], yerr=data_proton_hist[2], fmt='o', color='black')
plt.savefig(("%s/OnePProtonP.png" % dir_name))
figNumber = figNumber + 1

fig = plt.figure(figNumber)
plt.hist(twoPProtonStack, bins=50, stacked=True, range=protonMomentumRange, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights = twoPProtonWeights )
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Ext'])
plt.title("Leading Proton Range Momentum (1 mu 2 p)")
plt.xlabel("Momentum (GeV/c)")
plt.ylabel("Number of Protons")
data_proton_hist = dataify(dataProtons.query('track_range_mom_p_count == 2 & track_mcs_mom_count ==1 & isLeadingP == True')['track_range_mom_p_max'].to_numpy(), 50, protonMomentumRange)
plt.errorbar(data_proton_hist[0], data_proton_hist[1], yerr=data_proton_hist[2], fmt='o', color='black')
plt.savefig(("%s/TwoPLeadingProtonP.png" % dir_name))
figNumber = figNumber + 1

fig = plt.figure(figNumber)
plt.hist(incPhiStack, bins=50, stacked=True, range=phiRange, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights = incMuonWeights )
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Ext'])
plt.title("Muon Phi Angle")
plt.xlabel("Angle / pi (radians)")
plt.ylabel("Number of Primary Muons")
data_phi_hist = dataify(dataMuons.query('track_mcs_mom_count ==1')['phi'].to_numpy(), 50, phiRange)
plt.errorbar(data_phi_hist[0], data_phi_hist[1], yerr=data_phi_hist[2], fmt='o', color='black')
plt.savefig(("%s/IncMuonPhi.png" % dir_name))
figNumber = figNumber + 1

fig = plt.figure(figNumber)
plt.hist(onePPhiStack , bins=50, stacked=True, range=phiRange, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights = onePMuonWeights )
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Ext'])
plt.title("Muon Phi Angle (1 mu 1 p)")
plt.xlabel("Angle / pi (radians)")
plt.ylabel("Number of Primary Muons (1 mu 1 p)")
data_phi_hist = dataify(dataMuons.query('track_range_mom_p_count == 1 & track_mcs_mom_count ==1')['phi'].to_numpy(), 50, phiRange)
plt.errorbar(data_phi_hist[0], data_phi_hist[1], yerr=data_phi_hist[2], fmt='o', color='black')
plt.savefig(("%s/OnePMuonPhi.png" % dir_name))
figNumber = figNumber + 1

fig = plt.figure(figNumber)
plt.hist(muonProtonTheta, bins=25, stacked=True, range=phiRange, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights = muonProtonThetaWeights )
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Ext'])
plt.title("Muon-Proton Opening Angle (1 mu N p)")
plt.xlabel("Cos(theta)")
plt.ylabel("Number of Muon-Proton Pairs")
data_opening_angle_hist = dataify(dotProduct(dataProtons.query('track_mcs_mom_count ==1 & isLeadingP==True'), dataMuons.query('track_mcs_mom_count ==1 & track_range_mom_p_count > 0'), "track_mcs_mom", "track_range_mom_p", ["track_dirx", "track_diry", "track_dirz"]), 25, phiRange)
plt.errorbar(data_opening_angle_hist[0], data_opening_angle_hist[1], yerr=data_opening_angle_hist[2], fmt='o', color='black')
plt.savefig(("%s/OnePMuPOpening.png" % dir_name))
figNumber = figNumber + 1

fig = plt.figure(figNumber)
plt.hist(muonProtonPhiDiff, bins=25, stacked=True, range=phiDiffRange, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights = muonProtonThetaWeights )
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Ext'])
plt.title("Muon-Proton Phi Difference (1 mu N p)")
plt.xlabel("Angle / pi (radians)")
plt.ylabel("Number of Muon-Proton Pairs")
data_phi_diff_hist = dataify(np.absolute(dataMuons.query('track_mcs_mom_count ==1 & track_range_mom_p_count > 0 ')['phi'].to_numpy() - dataProtons.query('track_mcs_mom_count ==1 & isLeadingP==True')['phi'].to_numpy()), 25, phiDiffRange)
plt.errorbar(data_phi_diff_hist[0], data_phi_diff_hist[1], yerr=data_phi_diff_hist[2], fmt='o', color='black')
plt.savefig(("%s/OnePMuPDeltaPhi.png" % dir_name))
figNumber = figNumber + 1

fig = plt.figure(figNumber)
plt.hist(nPMuonTheta, bins=40, stacked=True, range=cosineRange, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights =  nPMuonWeights)
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Ext'])
plt.title("Muon Cosine Theta (1 mu n p)")
plt.xlabel("Cosine (theta)")
plt.ylabel("Number of Muons")
data_cosine_angle_hist = dataify(dataMuons.query('track_mcs_mom_count ==1 & track_range_mom_p_count > 0 ')['track_dirz'].to_numpy(), 40, cosineRange)
plt.errorbar(data_cosine_angle_hist[0], data_cosine_angle_hist[1], yerr=data_cosine_angle_hist[2], fmt='o', color='black')
plt.savefig(("%s/nPMuonCosine.png" % dir_name))
figNumber = figNumber + 1

fig = plt.figure(figNumber)
plt.hist(nPProtonTheta, bins=40, stacked=True, range=cosineRange, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights = nPProtonWeights )
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Ext'])
plt.title("Proton Cosine Theta (1 mu n p)")
plt.xlabel("Cosine (theta)")
plt.ylabel("Number of Protons")
data_proton_cosine_angle_hist = dataify(dataProtons.query('track_mcs_mom_count ==1 & isLeadingP == True')['track_dirz'].to_numpy(), 40, cosineRange)
plt.errorbar(data_proton_cosine_angle_hist[0], data_proton_cosine_angle_hist[1], yerr=data_proton_cosine_angle_hist[2], fmt='o', color='black')
plt.savefig(("%s/nPLeadingProtonCosine.png" % dir_name))
figNumber = figNumber + 1

fig = plt.figure(figNumber)
plt.hist(onePMuonTheta, bins=20, stacked=True, range=cosineRange, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights = onePMuonWeights )
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Ext'])
plt.title("Muon Cosine Theta (1 mu 1 p)")
plt.xlabel("Cosine (theta)")
plt.ylabel("Number of Muons")
data_oneP_muon_cosine_angle_hist = dataify(dataMuons.query('track_range_mom_p_count == 1 & track_mcs_mom_count ==1')['track_dirz'].to_numpy(), 20, cosineRange)
plt.errorbar(data_oneP_muon_cosine_angle_hist[0], data_oneP_muon_cosine_angle_hist[1], yerr=data_oneP_muon_cosine_angle_hist[2], fmt='o', color='black')
plt.savefig(("%s/OnePMuonCosine.png" % dir_name))
figNumber = figNumber + 1

fig = plt.figure(figNumber)
plt.hist(onePProtonTheta, bins=20, stacked=True, range=cosineRange, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights = onePProtonWeights )
plt.legend(['QE', 'RES', 'DIS', '2p2h', 'NC / Other', 'Ext'])
plt.title("Proton Cosine Theta (1 mu 1 p)")
plt.xlabel("Cosine (theta)")
plt.ylabel("Number of Protons")
data_oneP_proton_cosine_angle_hist = dataify(dataProtons.query('track_range_mom_p_count == 1 & track_mcs_mom_count ==1')['track_dirz'].to_numpy(), 20, cosineRange)
plt.errorbar(data_oneP_proton_cosine_angle_hist[0], data_oneP_proton_cosine_angle_hist[1], yerr=data_oneP_proton_cosine_angle_hist[2], fmt='o', color='black')
plt.savefig(("%s/OnePProtonCosine.png" % dir_name))
figNumber = figNumber + 1

fig = plt.figure(figNumber)
plt.hist(inclusiveProtonCandidates, bins=40, stacked=True, range=chi2Range, color = ['b', 'g', 'y', 'r', 'grey', 'magenta'], weights = protonCandidateWeights )
plt.legend(['Muon', 'Proton', 'Pion', 'Electron', 'Other', 'Ext'])
plt.title("Chi2 Proton")
plt.xlabel("Chi2 p")
plt.ylabel("Number of Tracks")
data_oneP_proton_cosine_angle_hist = dataify(dataProtonCandidates['track_chi2_proton'].to_numpy(), 40, chi2Range)
plt.errorbar(data_oneP_proton_cosine_angle_hist[0], data_oneP_proton_cosine_angle_hist[1], yerr=data_oneP_proton_cosine_angle_hist[2], fmt='o', color='black')
plt.savefig(("%s/Chi2Proton.png" % dir_name))
figNumber = figNumber + 1

'''

sys.exit()
