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
  inputFrame.eval('mc_Ehad = mc_nu_energy - mc_nu_lepton_energy', inplace=True) #Insert the true energy transfer (nu)
  inputFrame.insert(inputFrame.shape[1], "mc_expQ2", [getQ2(x, y, z) for x, y, z in zip(inputFrame['mc_nu_energy'], inputFrame['mc_nu_lepton_energy'], inputFrame['mc_nu_lepton_theta'])] )
  inputFrame.insert(inputFrame.shape[1], "mc_expW", [getW(x, y, z) for x, y, z in zip(inputFrame['mc_Ehad'], inputFrame['mc_expQ2'], inputFrame['mc_channel'] ) ] )
  inputFrame.insert(inputFrame.shape[1], "mc_expXbj", [getXbj(x, y) for x, y in zip(inputFrame['mc_Ehad'], inputFrame['mc_expQ2'] ) ] )
  inputFrame.insert(inputFrame.shape[1], "mc_expY", [getInel(x, y) for x, y in zip(inputFrame['mc_Ehad'], inputFrame['mc_nu_energy'] ) ] )
  #inputFrame.insert(inputFrame.shape[1], "template_wgt", [getChanWeight(x, y) for x, y in zip(inputFrame['mc_nu_interaction_type'], inputFrame['mc_nu_ccnc'])] ) #Classify neutrino events based on CC / NC and event Type
  inputFrame.insert(inputFrame.shape[1], "pot", mcPOT)

def loadMCTrackInfo(inputFrame):
  inputFrame.insert(inputFrame.shape[1], 'particle', [getParticle(x) for x in inputFrame['mc_pdg'] ])

def loadTrackInfo(inputFrame, isMC=False):
  inputFrame.insert(inputFrame.shape[1], "DuplicatedEvent", inputFrame.duplicated())
  inputFrame.insert(inputFrame.shape[1], "phi", [getPhi(x, y) for x, y in zip(inputFrame['track_diry'], inputFrame['track_dirx'] ) ] )
  inputFrame.eval('track_chi2_ratio = track_chi2_proton / track_chi2_muon', inplace=True)
  inputFrame.insert(inputFrame.shape[1], "isContained", [isContained(x, y, z, a, b, c) for x, y, z, a, b, c in zip(inputFrame['vx'], inputFrame['vy'], inputFrame['vz'], inputFrame['track_endx'], inputFrame['track_endy'], inputFrame['track_endz'])])
  inputFrame.insert(inputFrame.shape[1], 'track_mom_best', [getBestP(x, y, z) for x, y, z in zip(inputFrame['isContained'], inputFrame['track_range_mom_mu'], inputFrame['track_mcs_mom'])])
  if(isMC):
    inputFrame.insert(inputFrame.shape[1], 'particle', [getParticle(x) for x in inputFrame['mc_pdg'] ])

def AggregateFrame(inputFrame, var, stat):
  if stat not in ["count", "max", "mean"]:
    print "Cannot aggregate based on stat %s" % stat
    return
  
  statFrame = inputFrame.groupby(level=["run", "subrun", "event"]).agg({var: [stat]})

  statFrame.columns = ["_".join(x) for x in statFrame.columns.ravel()]

  inputFrame = inputFrame.join(statFrame['%s_%s' % (var, stat)], on=["run", "subrun", "event"])
  
  if(stat == "max"): 
    inputFrame.eval('isMax_%s = (%s == %s_max)' % (var, var, var), inplace=True)


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
  MCScalarSum   = np.sum(mcSum)
  DataScalarSum = np.sum(data_hist[1])
  sumRatio = DataScalarSum / MCScalarSum 
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
  plt.ylabel("Data / MC")
  text = r'$\int \frac{data}{MC} = %.3f$' % sumRatio
  ax = plt.gca()
  ymax = ax.get_ylim()[1] 
  xmax = ax.get_xlim()[1]
  #print "Min %.2f Max %.2f" % (ax.get_xlim()[0], ax.get_xlim()[1])
  plt.text(0.5*xmax, 0.9*ymax, text, {'fontsize' : 18})
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
    elif(pdg == -13):
      return "muon+"  
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

InputFiles = ["/uboone/data/users/suprajab/MCC9/March2020/CCpi0Trees/bnb_eventtree.root", "/uboone/data/users/suprajab/MCC9/March2020/CCpi0Trees/ext_eventtree.root", "/uboone/data/users/suprajab/MCC9/March2020/CCpi0Trees/dirt_eventtree.root"]
selectionVariables = ["_fHasCandidateNeutrino", "_fHasCandidateMuon", "_fNChargedPiCandidates", "_fNProtonCandidates", "_fNLeadCandidates", "_fNSubleadCandidates", "_fNPi0Candidates", "_fNPairCandidates", "_fNOtherFancyPairs"]
plottingVariables = ["_fShowerEnergy0"]
selectionVariables.extend(plottingVariables)

ExtScale     = 1.02767


#Python library to read in ROOT ntuples files.
bnbEvents     = uproot.open(InputFiles[0])["eventtree"]
eventsData    = pd.DataFrame(bnbEvents.arrays(selectionVariables) )

#SAVE THIS
#print trackOverlay.loc[:100,['isContained', 'track_range_mom_mu', 'track_mcs_mom', 'track_mom_best']]

#tagDuplicateEvents(eventsData)
print eventsData

eventsData   =  eventsData[eventsData.DuplicatedEvent == False]


muonMomentumRange   = (0.0, 2.0)
protonMomentumRange = (0.0, 1.5)
phiRange = (-1.0, 1.0)
isSelectedRange = (0.0, 1.0)
phiDiffRange = (0.0, 2.0)
trkScoreRange= (0.0, 1.0)
chi2Range    = (0.0, 50.0)
chi2PRange   = (0.0, 350.0)
vtxRange     = (0.0, 6.0)
lengthRange  = (0.0, 200.0)
flashChi2Range = (0.0, 300.0)
pdgRange     = (0, 30)

sys.exit()
