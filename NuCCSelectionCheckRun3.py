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

def tagDuplicateEvents(inputFrame):
  inputFrame.insert(inputFrame.shape[1], "DuplicatedEvent", inputFrame.duplicated() )

def loadMCEventInfo(inputFrame):
  tagDuplicateEvents(inputFrame)
  inputFrame.insert(inputFrame.shape[1], "mc_channel", [getChan(x, y) for x, y in zip(inputFrame['mc_nu_interaction_type'], inputFrame['mc_nu_ccnc'])] ) #Classify neutrino events based on CC / NC and event Type
  #inputFrame.eval('mc_Ehad = mc_nu_energy - mc_nu_lepton_energy', inplace=True) #Insert the true energy transfer (nu)
  #inputFrame.insert(inputFrame.shape[1], "mc_expQ2", [getQ2(x, y, z) for x, y, z in zip(inputFrame['mc_nu_energy'], inputFrame['mc_nu_lepton_energy'], inputFrame['mc_nu_lepton_theta'])] )
  #inputFrame.insert(inputFrame.shape[1], "mc_expW", [getW(x, y, z) for x, y, z in zip(inputFrame['mc_Ehad'], inputFrame['mc_expQ2'], inputFrame['mc_channel'] ) ] )
  #inputFrame.insert(inputFrame.shape[1], "mc_expXbj", [getXbj(x, y) for x, y in zip(inputFrame['mc_Ehad'], inputFrame['mc_expQ2'] ) ] )
  #inputFrame.insert(inputFrame.shape[1], "mc_expY", [getInel(x, y) for x, y in zip(inputFrame['mc_Ehad'], inputFrame['mc_nu_energy'] ) ] )
  #inputFrame.insert(inputFrame.shape[1], "template_wgt", [getChanWeight(x, y) for x, y in zip(inputFrame['mc_nu_interaction_type'], inputFrame['mc_nu_ccnc'])] ) #Classify neutrino events based on CC / NC and event Type
  #inputFrame.insert(inputFrame.shape[1], "pot", mcPOT)
  #inputFrame.insert(inputFrame.shape[1], "isTrueCC", [isTrueCC(x, y) for x, y in zip(inputFrame['mc_pdg'], inputFrame['mc_nu_ccnc'])])
  inputFrame.insert(inputFrame.shape[1], 'flash_wgt', [getFlashWgt(x) for x in inputFrame['nu_flash_chi2'] ])


def loadTrackInfo(inputFrame, isMC=False):
  inputFrame.insert(inputFrame.shape[1], "DuplicatedEvent", inputFrame.duplicated())
  inputFrame.insert(inputFrame.shape[1], "phi", [getPhi(x, y) for x, y in zip(inputFrame['track_diry'], inputFrame['track_dirx'] ) ] )
  inputFrame.eval('track_chi2_ratio = track_chi2_proton / track_chi2_muon', inplace=True)
  inputFrame.insert(inputFrame.shape[1], "isContained", [isContained(x, y, z, a, b, c) for x, y, z, a, b, c in zip(inputFrame['vx'], inputFrame['vy'], inputFrame['vz'], inputFrame['track_endx'], inputFrame['track_endy'], inputFrame['track_endz'])])
  inputFrame.insert(inputFrame.shape[1], 'track_mom_best', [getBestP(x, y, z) for x, y, z in zip(inputFrame['isContained'], inputFrame['track_range_mom_mu'], inputFrame['track_mcs_mom'])])
  if(isMC):
    inputFrame.insert(inputFrame.shape[1], 'particle', [getParticle(x) for x in inputFrame['mc_pdg'] ])


def AggregateFrame(inputFrame, var, stat):
  if stat not in ["count", "max", "mean", "min"]:
    print "Cannot aggregate based on stat %s" % stat
    return
  
  statFrame = inputFrame.groupby(level=["run", "subrun", "event"]).agg({var: [stat]})

  statFrame.columns = ["_".join(x) for x in statFrame.columns.ravel()]

  inputFrame = inputFrame.join(statFrame['%s_%s' % (var, stat)], on=["run", "subrun", "event"])

  if(stat == "max"): 
    inputFrame.eval('isMax_%s = (%s == %s_max)' % (var, var, var), inplace=True)
  if(stat == "min"):
    inputFrame.eval('isMin_%s = (%s == %s_min)' % (var, var, var), inplace=True)  
  return inputFrame


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

def makeDataMCHistogram(mcList, mcWeights, dataList, binRange, nBins, filename, Titles, xRange = [], yRange = [], rRange = []):
  dir_name = "PlotDir"
  xLimits  = []
  yLimits  = []
  for x in xRange:
    xLimits.append(x)
  for y in yRange:
    yLimits.append(y)  
  if(len(xLimits) == 1):
    xLimits.insert(0, 0.0)
  if(len(yLimits) == 1):
    yLimits.insert(0, 0.0)
  '''
  print (mcList)
  raw_input()
  print (mcWeights)
  raw_input()
  '''
  #out = plt.hist(mcList, bins=nBins, stacked=True, range=binRange, weights = mcWeights )

  out = plt.hist(mcList, bins=nBins, stacked=True, range=binRange, color = ['b', 'g', 'y', 'r', 'grey', 'gold', 'magenta'], weights = mcWeights )

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
  
  if(len(xLimits) == 2):
    plt.xlim(xLimits)
  if(len(yLimits) == 2):
    plt.ylim(yLimits)
  
  data_hist = dataify(dataList, nBins, binRange)
  plt.errorbar(data_hist[0], data_hist[1], yerr=data_hist[2], fmt='o', color='black')
  plt.savefig("%s/%s.png" % ( dir_name, filename) )
  plt.close()

  makeDataMCRatioHistogram(mcList, mcWeights, dataList, binRange, nBins, filename, Titles, xLimits, rRange)

def makeDataMCRatioHistogram(mcList, mcWeights, dataList, binRange, nBins, filename, Titles, xRange = [], yRange = []):
  dir_name  = "PlotDir"
  xLimits  = []
  yLimits  = []
  for x in xRange:
    xLimits.append(x)
  for y in yRange:
    yLimits.append(y)  
  if(len(xLimits) == 1):
    xLimits.insert(0, 0.0)
  if(len(yLimits) == 1):
    yLimits.insert(0, -1.0)
  mcSum = np.full(nBins, 0.0 )
  for mc, weight in zip(mcList, mcWeights):
     mc_hist   = np.histogram(mc, bins=nBins, range=binRange, weights = weight )
     np.add(mc_hist[0], mcSum, out=mcSum)
  data_hist = dataify(dataList, nBins, binRange)
  MCScalarSum   = np.sum(mcSum)
  DataScalarSum = np.sum(data_hist[1])
  sumRatio = DataScalarSum / MCScalarSum 
  ratio = np.divide(data_hist[1], mcSum)
  err   = np.multiply(ratio, np.divide(1.0, data_hist[2]))
  np.nan_to_num(ratio, copy=False)
  np.nan_to_num(err, copy=False)

  fig, axi = plt.subplots() #create subplots so I can put a textbox in

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
  plt.ylabel("Data / MC")
  
  if(len(xLimits) == 2):
    plt.xlim(xLimits)
  if(len(yLimits) == 2):
    plt.ylim(yLimits)

  text = r'$\int \frac{data}{MC} = %.3f$' % sumRatio
  props = dict(boxstyle='round', facecolor='lightsteelblue', alpha=0.5)

  axi.text(0.75, 1.1, text, transform=axi.transAxes, fontsize=14, verticalalignment='top', bbox=props)
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
    else:
      return "Unkown"  

def getParticle(pdg):
    if (pdg == 13):
      return "muon"
    elif(pdg == 2212):
      return "proton"
    elif(pdg == 211 or pdg == -211):
      return "cpi"
    elif(pdg == 111):
      return "pi0"
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

def isFiducial(x, y, z):
  if(x > 245.0 or x < 10.0):
      return False
  elif(y > 100.0 or y < -100.0):  
      return False
  elif(z > 970.0 or z <10.0):
      return False
  else:  
      return True

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

def isTrueCC(pdg, isNC):
  if(pdg == 14 and not isNC):
    return True
  else:
    return False

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

def getFlashWgt(x):
  returnVal = expDegSix(x)
  if(returnVal < 0):
    return 1.0
  else:   
    return returnVal

def expDegSix(x):
   p = [ 1.99384917e-01, -1.05022405e-10,  6.34916308e-05, -3.44705817e-03,
         7.32059307e-02, -5.91696006e-01,  2.14667463e+00, -1.02545380e+00,
         3.65915734e-01]
   return(np.exp(-p[0]*x)*(p[1]*pow(x, 6) + p[2]*pow(x,5) + p[3]*pow(x,4) + p[4]*pow(x,3) + p[5]*pow(x,2) + p[6]*x + p[7]) + p[8])  

def tagCorrectSlice(nFlashes, isMinFlash, isMinNu):
   #If there is one flash == min, use the min flash score
   #otherwise use the min nu score
   return(isMinNu if nFlashes > 1 else isMinFlash )
   #return(nFlashes > 1 ? isMinFlash : isMinNu)

#InputFiles = ["/uboone/data/users/joelam/stv-ntuples-new/numu_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/bnb_5e19_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/extC1_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/dirt_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/extC2_run1.root"]

#InputFiles = ["/uboone/data/users/ametzler/DetVarNtuples/prodgenie_bnb_nu_overlay_DetVar_LYAttenuation.root", "/uboone/data/users/joelam/stv-ntuples-new/bnb_5e19_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/extC1_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/dirt_run1.root", "/uboone/data/users/joelam/stv-ntuples-new/extC2_run1.root"]


InputFiles = ["/uboone/data/users/joelam/stv-ntuples-new/numu_run3.root", "/uboone/data/users/joelam/stv-ntuples-new/bnb_1e19_run3.root", "/uboone/data/users/joelam/stv-ntuples-new/extG1_run3.root", "/uboone/data/users/joelam/stv-ntuples-new/dirt_run3.root"]


#OverlayScale  = 1.0
#ExtScale     = 0.97
numMCTemplates = 6
empty = []


minProtonChi2 = 60.0
maxMuonChi2   = 30.0
minRatioChi2  = 7.0
minTrackScore = 0.5
minNeutrinoScore = 0.1
minMuonTrackScore = 0.85
minTrackL = 20
maxVtxDist = 4
requiredGen = 2
numupdg = 14

#Python library to read in ROOT ntuples files.
overlayEvents = uproot.open(InputFiles[0])["NuCCanalyzer"]["Event"]
bnbEvents     = uproot.open(InputFiles[1])["NuCCanalyzer"]["Event"]
extEventsG1    = uproot.open(InputFiles[2])["NuCCanalyzer"]["Event"]
dirtEvents    = uproot.open(InputFiles[3])["NuCCanalyzer"]["Event"]

overlayPOT    = uproot.open(InputFiles[0])["NuCCanalyzer"]["subruns"]
dirtPOT       = uproot.open(InputFiles[3])["NuCCanalyzer"]["subruns"]

#dataPOT            = 9.534e+18
dataPOT           = 4.988e+18
bnbSpills         = 1193929.0                   
#bnbSpills          = 2299517.0
#extTriggers     = 137052818.0
extTriggers      = 73285067.0

#Scale factors, because we generate more simulation than data. We also do not take an equal ammount of on and off beam data (though it is close)
mcPOT         = pd.DataFrame(overlayPOT.arrays(["run", "subRun", "pot"] ) )
mcPOT.insert(mcPOT.shape[1], "DuplicatedEvent", mcPOT.duplicated() )
mcPOT = mcPOT[mcPOT.DuplicatedEvent == False]
sumPOT        = mcPOT['pot'].sum()
sumDirtPOT    = (pd.Series(dirtPOT.array("pot"))).sum()

#Create frames of the event tree (which has information about the interaction) and the duaghters tree (which has information about the particles within the interaction).
#Do this for "overlay" simulation, beam data, and off-beam data
#with uproot.open(InputFiles[0])["NuCCanalyzer"]["Event"] as overlayEvents:
filteredEvents   = pd.DataFrame(overlayEvents.arrays(["run", "subrun", "event", "mc_nu_interaction_type", "nu_score", "nu_flash_chi2", "obvious_cosmic_chi2", "nu_pdg", "daughters_start_contained", "nu_vx", "nu_vy", "nu_vz", "mc_nu_ccnc", "nu_mu_cc_selected", "mc_nu_lepton_energy", "mc_nu_energy", "mc_nu_lepton_theta"]) )

overlayDaughters = uproot.open(InputFiles[0])["NuCCanalyzer"]["Daughters"]
trackOverlay   = pd.DataFrame(overlayDaughters.arrays(["track_range_mom_mu", "track_mcs_mom", "track_range_mom_p", "track_is_muon_candidate", "track_score", "track_chi2_proton", "track_chi2_muon", "track_dirx", "track_diry", "track_dirz", "vx", "vy", "vz", "track_endx", "track_endy", "track_endz", "track_length", "vtx_distance", "generation", "mc_pdg", "run", "subrun", "event"] ) )
filteredEvents.insert(filteredEvents.shape[1], "isFiducial", [isFiducial(x, y, z) for x, y, z in zip(filteredEvents['nu_vx'], filteredEvents['nu_vy'], filteredEvents['nu_vz'])] )
filteredEvents.eval('flash_chi2_ratio = nu_flash_chi2 / obvious_cosmic_chi2', inplace=True)

overlayCVWeights = pd.read_csv("/uboone/data/users/joelam/MCWeights/Run3TuneWeights.csv", names=["run", "subrun", "event", "wgt_tune"], dtype={"run" : int, "subrun" : int, "event" : int, "wgt_tune" : float})
dirtCVWeights    = pd.read_csv("/uboone/data/users/joelam/MCWeights/DirtRun3CVWeights.csv", names=["run", "subrun", "event", "wgt_tune"], dtype={"run" : int, "subrun" : int, "event" : int, "wgt_tune" : float})

overlaySplineWeights = pd.read_csv("/uboone/data/users/joelam/MCWeights/Run3SplineWeights.csv", names=["run", "subrun", "event", "wgt_spline"], dtype={"run" : int, "subrun" : int, "event" : int, "wgt_spline" : float})
dirtSplineWeights    = pd.read_csv("/uboone/data/users/joelam/MCWeights/DirtRun3SplineWeights.csv", names=["run", "subrun", "event", "wgt_spline"], dtype={"run" : int, "subrun" : int, "event" : int, "wgt_spline" : float})

overlayCVWeights.insert(overlayCVWeights.shape[1], "DuplicatedEvent", overlayCVWeights.duplicated() ) #Tag the events which are duplicated
overlayCVWeights.replace([np.inf], 0.0, inplace=True)

dirtCVWeights.insert(dirtCVWeights.shape[1], "DuplicatedEvent", dirtCVWeights.duplicated() ) #Tag the events which are duplicated
dirtCVWeights.replace([np.inf], 0.0, inplace=True)

overlaySplineWeights.insert(overlaySplineWeights.shape[1], "DuplicatedEvent", overlaySplineWeights.duplicated() ) #Tag the events which are duplicated
overlaySplineWeights.replace([np.inf], 0.0, inplace=True)

dirtSplineWeights.insert(dirtSplineWeights.shape[1], "DuplicatedEvent", dirtSplineWeights.duplicated() ) #Tag the events which are duplicated
dirtSplineWeights.replace([np.inf], 0.0, inplace=True)

dataDaughters = uproot.open(InputFiles[1])["NuCCanalyzer"]["Daughters"]
trackData     = pd.DataFrame(dataDaughters.arrays(["track_range_mom_mu", "track_mcs_mom", "track_range_mom_p", "track_is_muon_candidate", "track_score", "track_chi2_proton", "track_chi2_muon", "track_dirx", "track_diry", "track_dirz","track_dirx", "track_diry", "track_dirz", "vx", "vy", "vz", "track_endx", "track_endy", "track_endz", "track_length", "vtx_distance", "is_track_daughter", "generation", "run", "subrun", "event"] ) )
filteredData  = pd.DataFrame(bnbEvents.arrays(["run", "subrun", "event", "nu_mu_cc_selected", "nu_score", "nu_flash_chi2", "obvious_cosmic_chi2", "daughters_start_contained", "nu_vx", "nu_vy", "nu_vz", "nu_pdg"]) )
filteredData.insert(filteredData.shape[1], "isFiducial", [isFiducial(x, y, z) for x, y, z in zip(filteredData['nu_vx'], filteredData['nu_vy'], filteredData['nu_vz'])] )
filteredData.eval('flash_chi2_ratio = nu_flash_chi2 / obvious_cosmic_chi2', inplace=True)

extDaughtersG1 = uproot.open(InputFiles[2])["NuCCanalyzer"]["Daughters"]
trackExtG1     = pd.DataFrame(extDaughtersG1.arrays(["track_range_mom_mu", "track_mcs_mom", "track_range_mom_p", "track_is_muon_candidate", "track_score", "track_chi2_proton", "track_chi2_muon", "track_dirx", "track_diry", "track_dirz", "track_dirx", "track_diry", "track_dirz", "vx", "vy", "vz", "track_endx", "track_endy", "track_endz", "track_length", "vtx_distance", "is_track_daughter", "generation", "run", "subrun", "event"] ) )
filteredExtG1  = pd.DataFrame(extEventsG1.arrays(["run", "subrun", "event", "nu_mu_cc_selected", "nu_score", "nu_flash_chi2", "obvious_cosmic_chi2", "daughters_start_contained", "nu_vx", "nu_vy", "nu_vz", "nu_pdg"]) )
filteredExtG1.insert(filteredExtG1.shape[1], "isFiducial", [isFiducial(x, y, z) for x, y, z in zip(filteredExtG1['nu_vx'], filteredExtG1['nu_vy'], filteredExtG1['nu_vz'])] )
filteredExtG1.eval('flash_chi2_ratio = nu_flash_chi2 / obvious_cosmic_chi2', inplace=True)
filteredExtG1.insert(filteredExtG1.shape[1], "wgt", (bnbSpills / (extTriggers) ) )

dirtDaughters = uproot.open(InputFiles[3])["NuCCanalyzer"]["Daughters"]
trackDirt     = pd.DataFrame(dirtDaughters.arrays(["track_range_mom_mu", "track_mcs_mom", "track_range_mom_p", "track_is_muon_candidate", "track_score", "track_chi2_proton", "track_chi2_muon", "track_dirx", "track_diry", "track_dirz", "track_dirx", "track_diry", "track_dirz", "vx", "vy", "vz", "track_endx", "track_endy", "track_endz", "track_length", "vtx_distance", "is_track_daughter", "generation", "run", "subrun", "event"] ) )
filteredDirt  = pd.DataFrame(dirtEvents.arrays(["run", "subrun", "event", "nu_mu_cc_selected", "nu_score", "nu_flash_chi2", "obvious_cosmic_chi2", "daughters_start_contained", "nu_vx", "nu_vy", "nu_vz", "nu_pdg"]) )
filteredDirt.insert(filteredDirt.shape[1], "isFiducial", [isFiducial(x, y, z) for x, y, z in zip(filteredDirt['nu_vx'], filteredDirt['nu_vy'], filteredDirt['nu_vz'])] )
filteredDirt.eval('flash_chi2_ratio = nu_flash_chi2 / obvious_cosmic_chi2', inplace=True)

filteredExt = filteredExtG1
trackExt    = trackExtG1

#Here, we calculate some additional event information that isn't part of the input ROOT ntuple
#This is because the grad. student who created the files didn't include this information
loadMCEventInfo(filteredEvents)
loadTrackInfo(trackOverlay, True)

tagDuplicateEvents(filteredDirt)
loadTrackInfo(trackDirt)

tagDuplicateEvents(filteredExt)
loadTrackInfo(trackExt)

tagDuplicateEvents(filteredData)
loadTrackInfo(trackData)

dataRuns = trackExt[['run', 'subrun']]
dataRuns.insert(dataRuns.shape[1], "DuplicatedEvent", dataRuns.duplicated() ) 

dataRuns   =  dataRuns[dataRuns.DuplicatedEvent == False]
dataRuns[['run', 'subrun']].to_csv("ExtRuns.txt", sep=" ", header=False, index_label=False, index=False) 

dataRuns = trackData[['run', 'subrun']]
dataRuns.insert(dataRuns.shape[1], "DuplicatedEvent", dataRuns.duplicated() ) 

dataRuns   =  dataRuns[dataRuns.DuplicatedEvent == False]
dataRuns[['run', 'subrun']].to_csv("OnBeamRuns.txt", sep=" ", header=False, index_label=False, index=False) 

#trackOverlay   = trackOverlay.query('run < @maxOverlayRun')
#filteredEvents = filteredEvents.query('run < @maxOverlayRun')

#Index the events and daugthers by the run, subrun, event tuple
#This is IMPORTANT. The only infomration we have to connect the two frames a priori is this set of 3 ints
#A single event can have multiple tracks (and often does!)
#Multiindexing makes our life much easier, cuz we can grab the event info for ANY track from it's multiindex
trackOverlay   =  trackOverlay.set_index(['run', 'subrun', 'event'])
filteredEvents =  filteredEvents.set_index(['run', 'subrun', 'event'])

filteredData   =  filteredData.set_index(['run', 'subrun', 'event'])
trackData      =  trackData.set_index(['run', 'subrun', 'event'])

filteredExt    = filteredExt.set_index(['run', 'subrun', 'event'])
trackExt       = trackExt.set_index(['run', 'subrun', 'event'])

filteredDirt    = filteredDirt.set_index(['run', 'subrun', 'event'])
trackDirt       = trackDirt.set_index(['run', 'subrun', 'event'])

#overlayCVWeights = overlayCVWeights.query('run < @maxOverlayRun')
overlayCVWeights = overlayCVWeights.set_index(['run', 'subrun', 'event'])
dirtCVWeights    = dirtCVWeights.set_index(['run', 'subrun', 'event'])

#overlaySplineWeights = overlaySplineWeights.query('run < @maxOverlayRun')
overlaySplineWeights = overlaySplineWeights.set_index(['run', 'subrun', 'event'])
dirtSplineWeights    = dirtSplineWeights.set_index(['run', 'subrun', 'event'])

filteredEvents   =  filteredEvents[filteredEvents.DuplicatedEvent == False]
filteredData     =  filteredData[filteredData.DuplicatedEvent == False]
filteredExt      =  filteredExt[filteredExt.DuplicatedEvent == False]
filteredDirt     =  filteredDirt[filteredDirt.DuplicatedEvent == False]

trackOverlay     =  trackOverlay[trackOverlay.DuplicatedEvent == False]
trackData        =  trackData[trackData.DuplicatedEvent == False]
trackExt         =  trackExt[trackExt.DuplicatedEvent == False]
trackDirt        =  trackDirt[trackDirt.DuplicatedEvent == False]

overlayCVWeights =  overlayCVWeights[overlayCVWeights.DuplicatedEvent == False]
dirtCVWeights    =  dirtCVWeights[dirtCVWeights.DuplicatedEvent == False]
overlaySplineWeights =  overlaySplineWeights[overlaySplineWeights.DuplicatedEvent == False]
dirtSplineWeights    =  dirtSplineWeights[dirtSplineWeights.DuplicatedEvent == False]

OverlayScale = dataPOT / sumPOT
DirtScale    = dataPOT / sumDirtPOT
print "MC POT: %e Data POT: %e Overlay Scale: %.3f Ext Scale: %.3f Dirt Scale: %.3f" % (sumPOT, dataPOT, OverlayScale, (bnbSpills / extTriggers), DirtScale)
print "Total MC POT: %e total MC events: %d" % (sumPOT, trackOverlay.shape[0])
print "Total Dirt POT: %e" % sumDirtPOT

#create a dict of event info we want to associate with each daughter.
#by doing this, we have the complete event information for each track.
#Really what we want is to look at the particles' properties as a funciton of the underlying event information
#This is extendible to any event varaible we want to associate to a particle
interactionInfo = ("mc_channel", "nu_mu_cc_selected","nu_score", "nu_pdg", "nu_flash_chi2", "obvious_cosmic_chi2", "flash_chi2_ratio", "nu_vx", "nu_vy", "nu_vz", "daughters_start_contained", "isFiducial", "flash_wgt") 

filteredEvents = AggregateFrame(filteredEvents, "nu_flash_chi2", "min")
filteredEvents = AggregateFrame(filteredEvents, "nu_score", "min")

#filteredEvents.eval('isCorrectSlice = (isMin_nu_score & isMin_nu_flash_chi2)', inplace=True)
#grouped = filteredEvents.query('isCorrectSlice == False').groupby(level=["run", "subrun", "event"]).agg( {"isMin_nu_score": "count"} )

grouped = filteredEvents.query('nu_flash_chi2 == nu_flash_chi2_min').groupby(level=["run", "subrun", "event"])
eventsFailed = grouped.size().to_frame(name='nu_flash_chi2_min_count') #How many tracks per event pass the min nu flash cut? This would be we have 2 flashes with an equal chi2
filteredEvents = filteredEvents.join(eventsFailed, on=["run", "subrun", "event"])
filteredEvents.insert(filteredEvents.shape[1], 'isCorrectSlice', [tagCorrectSlice(x, y, z) for x, y, z in zip(filteredEvents['nu_flash_chi2_min_count'], filteredEvents['isMin_nu_flash_chi2'], filteredEvents['isMin_nu_score'])])

filteredEvents = filteredEvents[filteredEvents.isCorrectSlice == True]

for field in interactionInfo:
  trackOverlay   = trackOverlay.join(filteredEvents['%s' % field], on=["run", "subrun", "event"])
  #nan = trackOverlay[['mc_channel', 'phi']]
  #print(nan[nan.isna().any(axis=1)])
  #raw_input() 


trackOverlay = trackOverlay.join(overlayCVWeights['wgt_tune'], on=["run", "subrun", "event"])
trackOverlay = trackOverlay.join(overlaySplineWeights['wgt_spline'], on=["run", "subrun", "event"])

extInfo = { "nu_mu_cc_selected", "nu_score", "nu_pdg", "nu_flash_chi2", "obvious_cosmic_chi2", "flash_chi2_ratio", "nu_vx", "nu_vy", "nu_vz", "daughters_start_contained", "isFiducial", "wgt" }

for field in extInfo:
  trackExt   = trackExt.join(filteredExt['%s' % field], on=["run", "subrun", "event"])
  

dirtInfo = ("nu_mu_cc_selected", "nu_score", "nu_flash_chi2", "obvious_cosmic_chi2", "flash_chi2_ratio", "nu_vx", "nu_vy", "nu_vz", "daughters_start_contained", "isFiducial", "nu_pdg")

for field in dirtInfo:
  trackDirt = trackDirt.join(filteredDirt['%s' % field], on=["run", "subrun", "event"])



weightsPreAverage = dirtCVWeights['wgt_tune'].to_numpy()
weightsPreAverageRMS = np.nanstd(weightsPreAverage)

plt.hist(weightsPreAverage, bins=150, stacked=False, range=(0.8, 1.4), color = 'black')
text = r'$\sigma = %.3f$' % weightsPreAverageRMS
ax = plt.gca()
ymax = ax.get_ylim()[1] 
xmax = ax.get_xlim()[1]
plt.text(0.7*xmax, 0.9*ymax, text, {'fontsize' : 18})
plt.title("Dirt Weights before Averaging")
plt.xlabel("Dirt Weight")
plt.ylabel("Number of Neutrino interactions")
plt.ticklabel_format(axis="y",style="sci", scilimits=(0, 4))
plt.savefig("PlotDir/DirtWeightsPreAverage.png")
plt.close()  

dirtCVWeightMeans     = dirtCVWeights.groupby(level=["run", "subrun", "event"]).agg({"wgt_tune" : ["mean"]})
dirtSplineWeightMeans = dirtSplineWeights.groupby(level=["run", "subrun", "event"]).agg({"wgt_spline" : ["mean"]})

dirtCVWeightMeans.columns = ["_".join(x) for x in dirtCVWeightMeans.columns.ravel()]
dirtSplineWeightMeans.columns = ["_".join(x) for x in dirtSplineWeightMeans.columns.ravel()]

weightsPostAverage = dirtCVWeightMeans['wgt_tune_mean'].to_numpy()
weightsPostAverageRMS = np.nanstd(weightsPostAverage)

plt.hist(weightsPostAverage, bins=150, stacked=False, range=(0.8, 1.4), color = 'black')
text = r'$\sigma = %.3f$' % weightsPostAverageRMS
ax = plt.gca()
ymax = ax.get_ylim()[1] 
xmax = ax.get_xlim()[1]
plt.text(0.7*xmax, 0.9*ymax, text, {'fontsize' : 18})
plt.title("Dirt Weights after Averaging")
plt.xlabel("Dirt Weight")
plt.ylabel("Number of Events")
plt.ticklabel_format(axis="y",style="sci", scilimits=(0, 4))
plt.savefig("PlotDir/DirtWeightsPostAverage.png")
plt.close()  

for field in dirtInfo:
  print ("Number of data tracks: %d" % len(trackData.index ) )
  trackData  = trackData.join(filteredData['%s' % field], on=["run", "subrun", "event"])
  

trackDirt = trackDirt.join(dirtCVWeightMeans['wgt_tune_mean'], on=["run", "subrun", "event"])
trackDirt = trackDirt.join(dirtSplineWeightMeans['wgt_spline_mean'], on=["run", "subrun", "event"])


muonMomentumRange   = (0.0, 2.0)
protonMomentumRange = (0.0, 1.5)
phiRange = (-1.1, 1.1)
isSelectedRange = (0.0, 1.0)
phiDiffRange = (0.0, 2.0)
trkScoreRange= (0.0, 1.0)
chi2Range    = (0.0, 50.0)
chi2PRange   = (0.0, 350.0)
vtxRange     = (0.0, 6.0)
lengthRange  = (0.0, 200.0)
flashChi2Range = (0.0, 300.0)
pdgRange     = (0, 30)


overlayWeights = np.full(trackOverlay.shape[0], OverlayScale )

dirtWeights    = np.full(trackDirt.shape[0], DirtScale )

trackOverlay.insert(trackOverlay.shape[1], "pot_wgt", overlayWeights )
trackOverlay.eval('wgt = pot_wgt*wgt_tune*wgt_spline', inplace=True)
trackOverlay.eval('wgt_opt = pot_wgt*wgt_tune*wgt_spline*flash_wgt', inplace=True)  

trackDirt.insert(trackDirt.shape[1], "pot_wgt", dirtWeights )
trackDirt.eval('wgt = pot_wgt*wgt_tune_mean*wgt_spline_mean', inplace=True)
trackDirt.eval('wgt_opt = pot_wgt*wgt_tune_mean*wgt_spline_mean', inplace=True)

trackOverlay.fillna({'track_chi2_ratio' : 0}, inplace=True)
#trackOverlay = trackOverlay[10:]

overlaySliceScoreStack = '''[trackOverlay.query('mc_channel == "QE"')['VAR'].to_numpy(), trackOverlay.query('mc_channel == "RES"')['VAR'].to_numpy(), trackOverlay.query('mc_channel == "DIS"')['VAR'].to_numpy(), trackOverlay.query('mc_channel == "2p2h"')['VAR'].to_numpy(), trackOverlay.query('mc_channel == "NC / Other" | mc_channel == "Unkown"')['VAR'].to_numpy(), trackDirt['VAR'].to_numpy(), trackExt['VAR'].to_numpy()]'''
exec( "incSliceScoreStack   = "  + re.sub(r'VAR', 'nu_score', overlaySliceScoreStack) )
exec( "incMuonPhiStack      = "  + re.sub(r'VAR', 'phi', overlaySliceScoreStack) )
exec( "incSliceScorekWeights     = "  + re.sub(r'VAR', 'wgt',    overlaySliceScoreStack) )
'''
print(incMuonPhiStack[4:7])
raw_input()

print(incSliceScorekWeights[4:7])
raw_input()
'''
#makeDataMCHistogram(incSliceScoreStack, incSliceScorekWeights, trackData['nu_score'].to_numpy(), trkScoreRange, 25, "IncSliceScore", ["Slice Score", "Score", "Number of Daughters"])
makeDataMCHistogram(incMuonPhiStack,    incSliceScorekWeights, trackData['phi'].to_numpy(), phiRange, 25, "IncMuonPhi", ["Muon Phi Angle", "Angle / pi (radians)", "Number of Daughters"])

exec( "incIsSelectedStack   = "  + re.sub(r'VAR', 'nu_mu_cc_selected', overlaySliceScoreStack) )

makeDataMCHistogram(incIsSelectedStack, incSliceScorekWeights, trackData['nu_mu_cc_selected'].to_numpy(), isSelectedRange, 2, "IncIsSelected", ["Selected", "Selected", "Number of Daughters"])

extNuScore     = trackExt.query('DuplicatedEvent == False')
dirtNuScore    = trackDirt.query('DuplicatedEvent == False')
overlayNuScore = trackOverlay.query('DuplicatedEvent == False')
dataNuScore    = trackData.query('DuplicatedEvent == False')

overlayTrackScoreStack = '''[overlayNuScore.query('mc_channel == "QE"')['VAR'].to_numpy(), overlayNuScore.query('mc_channel == "RES"')['VAR'].to_numpy(), overlayNuScore.query('mc_channel == "DIS"')['VAR'].to_numpy(), overlayNuScore.query('mc_channel == "2p2h"')['VAR'].to_numpy(), overlayNuScore.query('mc_channel == "NC / Other"')['VAR'].to_numpy(), dirtNuScore['VAR'].to_numpy(), extNuScore['VAR'].to_numpy()]'''

exec( "incTrkScoreStack   = "  + re.sub(r'VAR', 'track_score', overlayTrackScoreStack) )
exec( "incTrkScorekWeights     = "  + re.sub(r'VAR', 'wgt',    overlayTrackScoreStack) )

makeDataMCHistogram(incTrkScoreStack, incTrkScorekWeights, dataNuScore['track_score'].to_numpy(), trkScoreRange, 25, "IncTrkScore", ["Track Score", "Score", "Number of Tracks"])


extTrackScore     = trackExt.query('DuplicatedEvent == False & track_score > @minTrackScore')
dirtTrackScore    = trackDirt.query('DuplicatedEvent == False & track_score > @minTrackScore')
overlayTrackScore = trackOverlay.query('DuplicatedEvent == False & track_score > @minTrackScore')
dataTrackScore    = trackData.query('DuplicatedEvent == False & track_score > @minTrackScore')

overlayTrackScoreStack = '''[overlayTrackScore.query('mc_channel == "QE"')['VAR'].to_numpy(), overlayTrackScore.query('mc_channel == "RES"')['VAR'].to_numpy(), overlayTrackScore.query('mc_channel == "DIS"')['VAR'].to_numpy(), overlayTrackScore.query('mc_channel == "2p2h"')['VAR'].to_numpy(), overlayTrackScore.query('mc_channel == "NC / Other"')['VAR'].to_numpy(), dirtTrackScore['VAR'].to_numpy(), extTrackScore['VAR'].to_numpy()]'''

exec( "incVtxStack   = "  + re.sub(r'VAR', 'vtx_distance', overlayTrackScoreStack) )
exec( "incChiSqrStackWeights     = "  + re.sub(r'VAR', 'wgt',    overlayTrackScoreStack) )


#makeDataMCHistogram(incVtxStack, incChiSqrStackWeights, dataTrackScore['vtx_distance'].to_numpy(), vtxRange, 30, "VtxDistance", ["Vertex Distance", "Distance to Vertex (cm)", "Number of Tracks"])

exec( "incTrkLStack   = "  + re.sub(r'VAR', 'track_length', overlayTrackScoreStack) )

#makeDataMCHistogram(incTrkLStack, incChiSqrStackWeights, dataTrackScore['track_length'].to_numpy(), lengthRange, 20, "TrkL", ["Track Length", "Track Length (cm)", "Number of Tracks"])

extPIDScore     = trackExt.query('DuplicatedEvent == False & track_score > @minMuonTrackScore & vtx_distance < @maxVtxDist & track_length > @minTrackL & generation == @requiredGen')
dirtPIDScore    = trackDirt.query('DuplicatedEvent == False & track_score > @minMuonTrackScore  & vtx_distance < @maxVtxDist & track_length > @minTrackL & generation == @requiredGen')
overlayPIDScore = trackOverlay.query('DuplicatedEvent == False & track_score > @minMuonTrackScore  & vtx_distance < @maxVtxDist & track_length > @minTrackL & generation == @requiredGen')
dataPIDScore    = trackData.query('DuplicatedEvent == False & track_score > @minMuonTrackScore  & vtx_distance < @maxVtxDist & track_length > @minTrackL & generation == @requiredGen')

overlayPIDStack = '''[overlayPIDScore.query('mc_channel == "QE"')['VAR'].to_numpy(), overlayPIDScore.query('mc_channel == "RES"')['VAR'].to_numpy(), overlayPIDScore.query('mc_channel == "DIS"')['VAR'].to_numpy(), overlayPIDScore.query('mc_channel == "2p2h"')['VAR'].to_numpy(), overlayPIDScore.query('mc_channel == "NC / Other"')['VAR'].to_numpy(), dirtPIDScore['VAR'].to_numpy(), extPIDScore['VAR'].to_numpy()]'''


exec( "incChiSqrStack   = "  + re.sub(r'VAR', 'track_chi2_muon', overlayPIDStack) )
exec( "incChiSqrStackWeights     = "  + re.sub(r'VAR', 'wgt',    overlayPIDStack) )

makeDataMCHistogram(incChiSqrStack, incChiSqrStackWeights, dataPIDScore['track_chi2_muon'].to_numpy(), chi2Range, 50, "IncChi2Muon", ["Chi2 Muon", "Chi2", "Number of Tracks"])

exec( "incChiSqrPStack   = "  + re.sub(r'VAR', 'track_chi2_proton', overlayPIDStack) )
exec( "incChiSqrPStackWeights     = "  + re.sub(r'VAR', 'wgt',      overlayPIDStack) )

#makeDataMCHistogram(incChiSqrPStack, incChiSqrPStackWeights, dataPIDScore['track_chi2_proton'].to_numpy(), chi2PRange, 35, "IncChi2Proton", ["Chi2 Proton", "Chi2", "Number of Tracks"])

exec( "incChiSqrPMuStack   = "  + re.sub(r'VAR', 'track_chi2_ratio', overlayPIDStack) )
#exec( "incChiSqrPStackWeights     = "  + re.sub(r'VAR', 'wgt',      overlayPIDStack) )

makeDataMCHistogram(incChiSqrPMuStack, incChiSqrPStackWeights, dataPIDScore['track_chi2_ratio'].to_numpy(), chi2Range, 50, "IncChi2Ratio", ["Chi2 Ratio", "Chi2", "Number of Tracks"], [5.0, 25.0], [200.0, 3000.0])

extMuonCandidates      = trackExt.query('DuplicatedEvent == False & track_score > @minMuonTrackScore  & vtx_distance < @maxVtxDist & track_length > @minTrackL & generation == @requiredGen & track_chi2_proton > @minProtonChi2 & track_chi2_muon < @maxMuonChi2 & track_chi2_ratio > @minRatioChi2')
dirtMuonCandidates     = trackDirt.query('DuplicatedEvent == False & track_score > @minMuonTrackScore &  vtx_distance < @maxVtxDist & track_length > @minTrackL & generation == @requiredGen & track_chi2_proton > @minProtonChi2 & track_chi2_muon < @maxMuonChi2 & track_chi2_ratio > @minRatioChi2')
overlayMuonCandidates  = trackOverlay.query('DuplicatedEvent == False & track_score > @minMuonTrackScore & vtx_distance < @maxVtxDist & track_length > @minTrackL & generation == @requiredGen & track_chi2_proton > @minProtonChi2 & track_chi2_muon < @maxMuonChi2 & track_chi2_ratio > @minRatioChi2')
dataMuonCandidates     = trackData.query('DuplicatedEvent == False & track_score > @minMuonTrackScore &  vtx_distance < @maxVtxDist & track_length > @minTrackL & generation == @requiredGen & track_chi2_proton > @minProtonChi2 & track_chi2_muon < @maxMuonChi2 & track_chi2_ratio > @minRatioChi2')


#Select only the primary muon track based on the longest track
extMuonCandidates = AggregateFrame(extMuonCandidates, "track_length", "max")

#statFrame = extMuonCandidates.groupby(level=["run", "subrun", "event"]).agg({"track_length": ["max"]})
#statFrame.columns = ["_".join(x) for x in statFrame.columns.ravel()]
#extMuonCandidates = extMuonCandidates.join(statFrame['track_length_max'], on=["run", "subrun", "event"])  
#extMuonCandidates.eval('isLongestTrack = (track_length == track_length_max)', inplace=True)

dirtMuonCandidates = AggregateFrame(dirtMuonCandidates, "track_length", "max")

#statFrame = dirtMuonCandidates.groupby(level=["run", "subrun", "event"]).agg({"track_length": ["max"]})
#statFrame.columns = ["_".join(x) for x in statFrame.columns.ravel()]
#dirtMuonCandidates = dirtMuonCandidates.join(statFrame['track_length_max'], on=["run", "subrun", "event"])  
#dirtMuonCandidates.eval('isLongestTrack = (track_length == track_length_max)', inplace=True)

overlayMuonCandidates = AggregateFrame(overlayMuonCandidates, "track_length", "max")

#statFrame = overlayMuonCandidates.groupby(level=["run", "subrun", "event"]).agg({"track_length": ["max"]})
#statFrame.columns = ["_".join(x) for x in statFrame.columns.ravel()]
#overlayMuonCandidates = overlayMuonCandidates.join(statFrame['track_length_max'], on=["run", "subrun", "event"])  
#overlayMuonCandidates.eval('isLongestTrack = (track_length == track_length_max)', inplace=True)

dataMuonCandidates = AggregateFrame(dataMuonCandidates, "track_length", "max")

#statFrame = dataMuonCandidates.groupby(level=["run", "subrun", "event"]).agg({"track_length": ["max"]})
#statFrame.columns = ["_".join(x) for x in statFrame.columns.ravel()]
#dataMuonCandidates = dataMuonCandidates.join(statFrame['track_length_max'], on=["run", "subrun", "event"])  
#dataMuonCandidates.eval('isLongestTrack = (track_length == track_length_max)', inplace=True)

#isMax_track_length
overlayPrimMuonStack = '''[overlayMuonCandidates.query('mc_channel == "QE" & isMax_track_length == True')['VAR'].to_numpy(), overlayMuonCandidates.query('mc_channel == "RES" & isMax_track_length == True')['VAR'].to_numpy(), overlayMuonCandidates.query('mc_channel == "DIS" & isMax_track_length == True')['VAR'].to_numpy(), overlayMuonCandidates.query('mc_channel == "2p2h" & isMax_track_length == True')['VAR'].to_numpy(), overlayMuonCandidates.query('mc_channel == "NC / Other" & isMax_track_length == True')['VAR'].to_numpy(), dirtMuonCandidates.query('isMax_track_length == True')['VAR'].to_numpy(), extMuonCandidates.query('isMax_track_length== True')['VAR'].to_numpy()]'''

exec( "incPrimMuonStack   = "  + re.sub(r'VAR', 'track_length', overlayPrimMuonStack) )
exec( "incPrimMuonStackWeights     = "  + re.sub(r'VAR', 'wgt', overlayPrimMuonStack) )

makeDataMCHistogram(incPrimMuonStack, incPrimMuonStackWeights, dataMuonCandidates.query('isMax_track_length == True')['track_length'].to_numpy(), lengthRange, 20, "PrimMuonL", ["Track Length", "Track Length (cm)", "Number of Events"])

exec( "incPrimMuonChi2Mu  = "  + re.sub(r'VAR', 'track_chi2_muon', overlayPrimMuonStack) )

#makeDataMCHistogram(incPrimMuonChi2Mu, incPrimMuonStackWeights, dataMuonCandidates.query('isLongestTrack == True')['track_chi2_muon'].to_numpy(), chi2Range, 50, "PrimMuonChi2Muon", ["Chi2 Muon", "Chi2", "Number of Events"])

exec( "incPrimMuonChi2Proton  = "  + re.sub(r'VAR', 'track_chi2_proton', overlayPrimMuonStack) )

#makeDataMCHistogram(incPrimMuonChi2Proton, incPrimMuonStackWeights, dataMuonCandidates.query('isLongestTrack == True')['track_chi2_proton'].to_numpy(), chi2PRange, 35, "PrimMuonChi2Proton", ["Chi2 Proton", "Chi2", "Number of Events"])

exec( "incPrimMuonChi2Ratio  = "  + re.sub(r'VAR', 'track_chi2_ratio', overlayPrimMuonStack) )

#makeDataMCHistogram(incPrimMuonChi2Ratio, incPrimMuonStackWeights, dataMuonCandidates.query('isLongestTrack == True')['track_chi2_ratio'].to_numpy(), chi2Range, 50, "PrimMuonChi2Ratio", ["Chi2 Ratio", "Chi2", "Number of Events"])

exec( "incPrimMuonNuScoreStack   = "  + re.sub(r'VAR', 'nu_score', overlayPrimMuonStack) )

#makeDataMCHistogram(incPrimMuonNuScoreStack, incPrimMuonStackWeights, dataMuonCandidates.query('isLongestTrack == True')['nu_score'].to_numpy(), trkScoreRange, 50, "PrimMuonNuSocre", ["Topological Score", "Neutrino ID", "Number of Events"])

exec( "incPrimMuonChi2FlashStack   = "  + re.sub(r'VAR', 'nu_flash_chi2', overlayPrimMuonStack) )

makeDataMCHistogram(incPrimMuonChi2FlashStack, incPrimMuonStackWeights, dataMuonCandidates.query('isMax_track_length == True')['nu_flash_chi2'].to_numpy(), (0, 50), 50, "PrimMuonFlashChi2", ["Flash Chi2", "Chi2", "Number of Events"])

#makeDataMCHistogram(incPrimMuonChi2FlashStack, incPrimMuonStackWeights, dataMuonCandidates.query('isLongestTrack == True')['nu_flash_chi2'].to_numpy(), (0, 15), 60, "PrimMuonFlashChi2Zoom", ["Flash Chi2", "Chi2", "Number of Events"])


exec( "incPrimMuonDaughtersStack   = "  + re.sub(r'VAR', 'daughters_start_contained', overlayPrimMuonStack) )

#makeDataMCHistogram(incPrimMuonDaughtersStack, incPrimMuonStackWeights, dataMuonCandidates.query('isLongestTrack == True')['daughters_start_contained'].to_numpy(), isSelectedRange, 2, "PrimMuonDaugthersContained", ["Is Daughter Contained", "Daugthers Contained", "Number of Events"])

exec( "incPrimMuonDaughtersStack   = "  + re.sub(r'VAR', 'daughters_start_contained', overlayPrimMuonStack) )

#makeDataMCHistogram(incPrimMuonDaughtersStack, incPrimMuonStackWeights, dataMuonCandidates.query('isLongestTrack == True')['daughters_start_contained'].to_numpy(), isSelectedRange, 2, "PrimMuonDaugthersContained", ["Is Daughter Contained", "Daugthers Contained", "Number of Events"])

exec( "incPrimMuonPDGStack   = "  + re.sub(r'VAR', 'nu_pdg', overlayPrimMuonStack) )

#makeDataMCHistogram(incPrimMuonPDGStack, incPrimMuonStackWeights, dataMuonCandidates.query('isLongestTrack == True')['nu_pdg'].to_numpy(), pdgRange, 30, "PrimMuonPDG", ["Event PDG", "Pandora PDG", "Number of Events"])

exec( "incPrimMuonIsFiducialStack   = "  + re.sub(r'VAR', 'isFiducial', overlayPrimMuonStack) )

#makeDataMCHistogram(incPrimMuonIsFiducialStack, incPrimMuonStackWeights, dataMuonCandidates.query('isLongestTrack == True')['isFiducial'].to_numpy(), isSelectedRange, 2, "PrimMuonDaugthersIsFiducial", ["Is Fiducial", "Vertices in Fiducial Volume", "Number of Events"])
#exec( "incPrimMuonIsSelectedStack   = "  + re.sub(r'VAR', 'nu_mu_cc_selected', overlayPrimMuonStack) )

exec( "incPrimMuonFlashChi2Ratio   = "  + re.sub(r'VAR', 'flash_chi2_ratio', overlayPrimMuonStack) )

#makeDataMCHistogram(incPrimMuonFlashChi2Ratio, incPrimMuonStackWeights, dataMuonCandidates.query('isLongestTrack == True')['flash_chi2_ratio'].to_numpy(), (0,5.0), 50, "PrimMuonFlashChi2CosmicRatio", ["Flash Chi2", "Chi2 Ratio", "Number of Events"])
#exec( "incPrimMuonIsSelectedStack   = "  + re.sub(r'VAR', 'nu_mu_cc_selected', overlayPrimMuonStack) )


maxFlashChi2 = 10
minNeutrinoScoreFlashFails = 0.25
maxFlashChi2Ratio  = 5

extInclusiveEvents = extMuonCandidates.query('isMax_track_length == True & isFiducial == True & nu_pdg == @numupdg & daughters_start_contained == True & flash_chi2_ratio < @maxFlashChi2Ratio & nu_score > @minNeutrinoScore & (nu_flash_chi2 < @maxFlashChi2 | nu_score > @minNeutrinoScoreFlashFails)')
dirtInclusiveEvents = dirtMuonCandidates.query('isMax_track_length == True & isFiducial == True & nu_pdg == @numupdg & daughters_start_contained == True & flash_chi2_ratio < @maxFlashChi2Ratio & nu_score > @minNeutrinoScore & (nu_flash_chi2 < @maxFlashChi2 | nu_score > @minNeutrinoScoreFlashFails)')
overlayInclusiveEvents = overlayMuonCandidates.query('isMax_track_length == True & isFiducial == True & nu_pdg == @numupdg & daughters_start_contained == True & flash_chi2_ratio < @maxFlashChi2Ratio & nu_score > @minNeutrinoScore & (nu_flash_chi2 < @maxFlashChi2 | nu_score > @minNeutrinoScoreFlashFails)')
dataInclusiveEvents = dataMuonCandidates.query('isMax_track_length == True & isFiducial == True & nu_pdg == @numupdg & daughters_start_contained == True & flash_chi2_ratio < @maxFlashChi2Ratio & nu_score > @minNeutrinoScore & (nu_flash_chi2 < @maxFlashChi2 | nu_score > @minNeutrinoScoreFlashFails)')


extInclusiveEventsContained = extMuonCandidates.query('isMax_track_length == True & isFiducial == True & nu_pdg == @numupdg & daughters_start_contained == True & flash_chi2_ratio < @maxFlashChi2Ratio & nu_score > @minNeutrinoScore & isContained == True')
dirtInclusiveEventsContained = dirtMuonCandidates.query('isMax_track_length == True & isFiducial == True & nu_pdg == @numupdg & daughters_start_contained == True & flash_chi2_ratio < @maxFlashChi2Ratio & nu_score > @minNeutrinoScore & isContained == True')
overlayInclusiveEventsContained = overlayMuonCandidates.query('isMax_track_length== True & isFiducial == True & nu_pdg == @numupdg & daughters_start_contained == True & flash_chi2_ratio < @maxFlashChi2Ratio & nu_score > @minNeutrinoScore & isContained == True')
dataInclusiveEventsContained = dataMuonCandidates.query('isMax_track_length == True & isFiducial == True & nu_pdg == @numupdg & daughters_start_contained == True & flash_chi2_ratio < @maxFlashChi2Ratio & nu_score > @minNeutrinoScore & isContained == True')

overlayInclusiveStack = '''[overlayInclusiveEvents.query('mc_channel == "QE"')['VAR'].to_numpy(), overlayInclusiveEvents.query('mc_channel == "RES"')['VAR'].to_numpy(), overlayInclusiveEvents.query('mc_channel == "DIS"')['VAR'].to_numpy(), overlayInclusiveEvents.query('mc_channel == "2p2h"')['VAR'].to_numpy(), overlayInclusiveEvents.query('mc_channel == "NC / Other"')['VAR'].to_numpy(), dirtInclusiveEvents['VAR'].to_numpy(), extInclusiveEvents['VAR'].to_numpy()]'''
overlayInclusiveStackOptWgt = '''[overlayInclusiveEvents.query('mc_channel == "QE"')['VAR'].to_numpy(), overlayInclusiveEvents.query('mc_channel == "RES"')['VAR'].to_numpy(), overlayInclusiveEvents.query('mc_channel == "DIS"')['VAR'].to_numpy(), overlayInclusiveEvents.query('mc_channel == "2p2h"')['VAR'].to_numpy(), overlayInclusiveEvents.query('mc_channel == "NC / Other"')['VAR'].to_numpy(), dirtInclusiveEvents['wgt'].to_numpy(), extInclusiveEvents['wgt'].to_numpy()]'''


exec( "overlayIsSelectedInclusiveStack     = "  + re.sub(r'VAR', 'nu_mu_cc_selected', overlayInclusiveStack) )
exec( "overlayIsSelectedInclusiveWeights   = "  + re.sub(r'VAR', 'wgt', overlayInclusiveStack) )
exec( "overlayIsSelectedInclusiveOpticalWeights   = "  + re.sub(r'VAR', 'wgt_opt', overlayInclusiveStackOptWgt) )

makeDataMCHistogram(overlayIsSelectedInclusiveStack, overlayIsSelectedInclusiveWeights, dataInclusiveEvents['nu_mu_cc_selected'].to_numpy(), isSelectedRange, 2, "InclusiveEventsIsSelected", ["Passes Selection", "Is Selected", "Number of Events"])

exec( "overlayPrimMuonChi2FlashInclusiveStack     = "  + re.sub(r'VAR', 'nu_flash_chi2', overlayInclusiveStack) )

makeDataMCHistogram(overlayPrimMuonChi2FlashInclusiveStack, overlayIsSelectedInclusiveWeights, dataInclusiveEvents['nu_flash_chi2'].to_numpy(), (0, 50), 50, "InclusiveEventsPrimMuonFlashChi2", ["Flash Chi2", "Chi2", "Number of Events"], [], [], [0.0, 2.0])

makeDataMCHistogram(overlayPrimMuonChi2FlashInclusiveStack, overlayIsSelectedInclusiveOpticalWeights, dataInclusiveEvents['nu_flash_chi2'].to_numpy(), (0, 50), 50, "InclusiveEventsPrimMuonFlashChi2RW", ["Flash Chi2", "Chi2", "Number of Events"], [], [], [0.0, 2.0])

exec( "overlayPrimMuonPhiInclusiveStack     = "  + re.sub(r'VAR', 'phi', overlayInclusiveStack) )

makeDataMCHistogram(overlayPrimMuonPhiInclusiveStack, overlayIsSelectedInclusiveWeights, dataInclusiveEvents['phi'].to_numpy(), phiRange, 30, "InclusiveEventsPrimMuonPhi", ["Muon Phi Angle", "Angle / pi (radians)", "Number of Primary Muons"], [], [], [0.7, 1.4])
makeDataMCHistogram(overlayPrimMuonPhiInclusiveStack, overlayIsSelectedInclusiveOpticalWeights, dataInclusiveEvents['phi'].to_numpy(), phiRange, 30, "InclusiveEventsPrimMuonPhiRW", ["Muon Phi Angle", "Angle / pi (radians)", "Number of Primary Muons"], [], [], [0.7, 1.4])

exec( "overlayPrimMuonNuScoreStack     = "  + re.sub(r'VAR', 'nu_score', overlayInclusiveStack) )

makeDataMCHistogram(overlayPrimMuonNuScoreStack, overlayIsSelectedInclusiveWeights, dataInclusiveEvents['nu_score'].to_numpy(), trkScoreRange, 50, "InclusiveEventsPrimMuonNuScore", ["Neutrino Score (Before Weights)", "nu_score", "Number of Primary Muons"], [], [], [0.7, 1.4])
makeDataMCHistogram(overlayPrimMuonNuScoreStack, overlayIsSelectedInclusiveOpticalWeights, dataInclusiveEvents['nu_score'].to_numpy(), trkScoreRange, 50, "InclusiveEventsPrimMuonNuScoreRW", ["Neutrino Score (After Weights)", "nu_score", "Number of Primary Muons"], [], [], [0.7, 1.4])

exec( "overlayPrimMuonPStack     = "  + re.sub(r'VAR', 'track_mom_best', overlayInclusiveStack) )

makeDataMCHistogram(overlayPrimMuonPStack, overlayIsSelectedInclusiveWeights, dataInclusiveEvents['track_mom_best'].to_numpy(), muonMomentumRange, 50, "InclusiveEventsPrimMuonP", ["Best Muon Momentum (Before Weights)", "Momentum (GeV/c)", "Number of Primary Muons"], [], [], [0.0, 2.0])
makeDataMCHistogram(overlayPrimMuonPStack, overlayIsSelectedInclusiveOpticalWeights, dataInclusiveEvents['track_mom_best'].to_numpy(), muonMomentumRange, 50, "InclusiveEventsPrimMuonPRW", ["Best Muon Momentum  (After Weights)", "Momentum (GeV/c)", "Number of Primary Muons"], [], [], [0.0, 2.0])

exec( "overlayPrimMuonThetaStack     = "  + re.sub(r'VAR', 'track_dirz', overlayInclusiveStack) )

makeDataMCHistogram(overlayPrimMuonThetaStack, overlayIsSelectedInclusiveWeights, dataInclusiveEvents['track_dirz'].to_numpy(), phiRange, 30, "InclusiveEventsPrimMuonTheta", ["Outgoing Muon Angle (Before Weights)", "Cosine (theta)", "Number of Primary Muons"], [], [], [0.0, 2.0])
makeDataMCHistogram(overlayPrimMuonThetaStack, overlayIsSelectedInclusiveOpticalWeights, dataInclusiveEvents['track_dirz'].to_numpy(), phiRange, 30, "InclusiveEventsPrimMuonThetaRW", ["Outgoing Muon Angle   (After Weights)", "Cosine (theta)", "Number of Primary Muons"], [], [], [0.0, 2.0])


sys.exit()
