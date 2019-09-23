import ROOT, math, sys, os
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
InputFiles = ["mc_tempalte_histos.root", "bnb_tempalte_histos.root", "ext_tempalte_histos.root"]
OverlayScale = 0.116
#OverlayScale  = 1.0
ExtScale     = 0.97
numMCTemplates = 6

empty = []
colors = [634, 418, 401, 603, 922, 616]
channels = ["CCQE", "CCRES", "CCDIS", "2p2h", "Other", "OOFV"]
'''
ROOT.gROOT.Reset()
'''
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetStatX(0.9)
ROOT.gStyle.SetStatY(0.9)
ROOT.gStyle.SetStatW(0.15)
ROOT.gStyle.SetStatH(0.1)
ROOT.gStyle.SetLabelSize(0.05, "x")
ROOT.gStyle.SetLabelSize(0.05, "y")
ROOT.gStyle.SetLineWidth(3)
ROOT.gStyle.SetLineColor(1)

min = ROOT.Math.Factory.CreateMinimizer('Minuit2')
min.SetMaxFunctionCalls(1000000)
min.SetMaxIterations(100000)
min.SetTolerance(0.001)
min.SetPrintLevel(1)

mcTemplatesP     = {"QE" : empty, "RES" : empty, "DIS" : empty, "MEC" : empty, "Other" : empty, "OFFV" : empty}
mcTemplatesTheta = {"QE" : empty, "RES" : empty, "DIS" : empty, "MEC" : empty, "Other" : empty, "OFFV" : empty}

MC_file     = ROOT.TFile.Open(InputFiles[0])
OnBeam_file = ROOT.TFile.Open(InputFiles[1])
Ext_file    = ROOT.TFile.Open(InputFiles[2])

histoDataP     = OnBeam_file.Get("h_muon_mom0")
histoDataCos   = OnBeam_file.Get("h_muon_costheta0")
histoDataP.SetMarkerStyle(8)
histoDataCos.SetMarkerStyle(8)

histoExtP      = Ext_file.Get("h_muon_mom0")
histoExtCos    = Ext_file.Get("h_muon_costheta0")

histoExtP.Scale(ExtScale)
histoExtCos.Scale(ExtScale)

histoExtP.SetFillColor(1)
histoExtP.SetFillColor(1)

histoExtP.SetFillStyle(3005)
histoExtCos.SetFillStyle(3005)

for itype in mcTemplatesP:
    index = chanToHistogram(itype)
    mcInfoP     = []
    mcInfoTheta = []
    histoP = MC_file.Get("h_muon_mom%d" % index)
    for bin in range(histoP.GetNbinsX()):
       mcInfoP.append(OverlayScale*histoP.GetBinContent(bin+1))

    histoTheta = MC_file.Get("h_muon_costheta%d" % index)
    for bin in range(histoTheta.GetNbinsX()):
       mcInfoTheta.append(OverlayScale*histoTheta.GetBinContent(bin+1))
    mcTemplatesP[itype] = mcInfoP
    mcTemplatesTheta[itype] = mcInfoTheta

#Get the Data and Ext
MomentumData = []
AngleData    = []
MomentumExt  = []
AngleExt     = []

for i in range(histoDataP.GetNbinsX()):
    MomentumData.append(histoDataP.GetBinContent(i+1))
    MomentumExt.append(histoExtP.GetBinContent(i+1))

for i in range(histoDataCos.GetNbinsX()):
    AngleData.append(histoDataCos.GetBinContent(i+1))
    AngleExt.append(histoExtCos.GetBinContent(i+1))

#print MomentumData
#print MomentumExt
#AllMC = mcTemplatesP["QE"] + mcTemplatesP["RES"] + mcTemplatesP["DIS"] + mcTemplatesP["MEC"] + mcTemplatesP["Other"] + mcTemplatesP["OFFV"]
#print len(mcTemplatesP["QE"])
'''
bin = 3
print mcTemplatesP["QE"][bin]
print mcTemplatesP["RES"][bin]
print mcTemplatesP["DIS"][bin]
print mcTemplatesP["MEC"][bin]
print mcTemplatesP["Other"][bin]
print mcTemplatesP["OFFV"][bin]
print MomentumExt[bin]
print "Sum: %.2f" % ( mcTemplatesP["QE"][bin] + mcTemplatesP["RES"][bin] + mcTemplatesP["DIS"][bin] + mcTemplatesP["MEC"][bin] + mcTemplatesP["Other"][bin] + mcTemplatesP["OFFV"][bin] + MomentumExt[bin])
'''
#print len(mcTemplatesP["QE"])
#print AllMC


fitter = TemplateFitter(mcTemplatesP, MomentumExt, MomentumData)
fitter.setFitRange( 0, 50 )
min.SetFunction( fitter )
# parameter number, name, initial value, step size, minimum value, maximum value
min.SetLimitedVariable( 0, 'MC',  5.0, 0.01, 0.01, 10.0 )
min.SetLimitedVariable( 1, 'Ext',  1.0, 0.01, 0.01, 10.0 )
min.SetLimitedVariable( 2, 'OOFV',  1.0, 0.01, 0.01, 2.0 )
#min.SetLimitedVariable( 2, 'DIS',  1.0, 0.01, 0.01, 2.0 )
#min.SetLimitedVariable( 3, 'MEC',  1.0, 0.01, 0.01, 2.0 )
#min.SetLimitedVariable( 4, 'Other',  1.0, 0.01, 0.01, 2.0 )

#min.SetLimitedVariable( 6, 'Ext',  1.0, 0.01, 0.01, 10.0 )

min.Minimize()
momentumScale = []
momentumScale.append(min.X()[0])
momentumScale.append(min.X()[1])
momentumScale.append(min.X()[2])
#print "Wut %.2f" % momentumScale[0]


fitter = TemplateFitter(mcTemplatesTheta, AngleExt, AngleData)
fitter.setFitRange( 0, 50 )
min.SetFunction( fitter )
# parameter number, name, initial value, step size, minimum value, maximum value
min.SetLimitedVariable( 0, 'MC',  5.0, 0.01, 0.01, 10.0 )
min.SetLimitedVariable( 1, 'Ext',  1.0, 0.01, 0.01, 10.0 )
min.SetLimitedVariable( 2, 'OOFV',  1.0, 0.01, 0.01, 2.0 )
#min.SetLimitedVariable( 2, 'DIS',  1.0, 0.01, 0.01, 2.0 )
#min.SetLimitedVariable( 3, 'MEC',  1.0, 0.01, 0.01, 2.0 )
#min.SetLimitedVariable( 4, 'Other',  1.0, 0.01, 0.01, 2.0 )
#min.SetLimitedVariable( 5, 'OOFV',  1.0, 0.01, 0.01, 2.0 )



min.Minimize()
thetaScale = []
thetaScale.append(min.X()[0])
thetaScale.append(min.X()[1])
thetaScale.append(min.X()[2])

#print "Wut in the but %.2f" % momentumScale[0]

c1 = ROOT.TCanvas("Momentum", "mom", 950, 950)
c2 = ROOT.TCanvas("Angle",    "cos", 950, 950)

maximum = 2200.0

MomentumnStack = ROOT.THStack()
MomentumnStack.Add(histoExtP)
AngleStack     = ROOT.THStack()
AngleStack.Add(histoExtCos)

legend1         = ROOT.TLegend(0.17, 0.55, 0.3, 0.9)
legend2         = ROOT.TLegend(0.67, 0.55, 0.9, 0.9)

print "Read in stuff"
for i in range(numMCTemplates):
   
   print "Reading in template %d" % i
   
   histoP = MC_file.Get("h_muon_mom%d" % i)
   histoTheta = MC_file.Get("h_muon_costheta%d" % i)
   
   histoP.SetLineColor(colors[i])
   histoP.SetFillColor(colors[i])
   histoP.Scale(OverlayScale)
       #if(histoP.GetMaximum() > maximum):
       # maximum = histoP.GetMaximum()
   histoP.SetMaximum(maximum)
   MomentumnStack.Add(histoP)

   legend1.AddEntry(histoP, channels[i], "l")
   
#c2.cd()
   
   histoTheta.SetLineColor(colors[i])
   histoTheta.SetFillColor(colors[i])
   histoTheta.Scale(OverlayScale)
   if(histoTheta.GetMaximum() > maximum):
     maximum = histoTheta.GetMaximum()
   histoTheta.SetMaximum(maximum)
   AngleStack.Add(histoTheta)

   legend2.AddEntry(histoTheta, channels[i], "l")

c1.cd()

MomentumnStack.SetMaximum(maximum)
MomentumnStack.Draw("hist")
histoDataP.Draw("E0 same")
legend2.Draw("same")
c1.Print("PreFitMomentum.png", "png")

c2.cd()

AngleStack.SetMaximum(maximum)
AngleStack.Draw("hist")
histoDataCos.Draw("E0 same")
legend1.Draw("same")
c2.Print("PreFitAngle.png", "png")

c1 = ROOT.TCanvas("Momentum", "mom", 950, 950)
c2 = ROOT.TCanvas("Angle",    "cos", 950, 950)

maximum = 2200.0

print "Scaling %.2f" % momentumScale[1]
print "Scaling %.2f" % thetaScale[1]

MomentumnStack = ROOT.THStack()
histoExtP.Scale(momentumScale[1])
MomentumnStack.Add(histoExtP)
AngleStack     = ROOT.THStack()
histoExtCos.Scale(thetaScale[1])
AngleStack.Add(histoExtCos)

legend1         = ROOT.TLegend(0.17, 0.55, 0.3, 0.9)
legend2         = ROOT.TLegend(0.67, 0.55, 0.9, 0.9)

print "Read in stuff Again"

for i in range(numMCTemplates):
    
    print "Reading in template %d with scale factors %.2f and %.2f" % (i, momentumScale[0], thetaScale[0])
        
    histoP = MC_file.Get("h_muon_mom%d" % i)
    histoTheta = MC_file.Get("h_muon_costheta%d" % i)
    
    
    histoP.SetLineColor(colors[i])
    histoP.SetFillColor(colors[i])
    histoP.Scale(momentumScale[0])
    histoP.SetMaximum(maximum)
    MomentumnStack.Add(histoP)
            
    legend1.AddEntry(histoP, channels[i], "l")
                    
    histoTheta.SetLineColor(colors[i])
    histoTheta.SetFillColor(colors[i])
    histoTheta.Scale(thetaScale[0])
    AngleStack.Add(histoTheta)
                            
    legend2.AddEntry(histoTheta, channels[i], "l")

c1.cd()

MomentumnStack.SetMaximum(maximum)
MomentumnStack.Draw("hist")
histoDataP.Draw("E0 same")
legend2.Draw("same")
c1.Print("PostFitMomentum.png", "png")

c2.cd()

AngleStack.SetMaximum(maximum)
AngleStack.Draw("hist")
histoDataCos.Draw("E0 same")
legend1.Draw("same")
c2.Print("PostFitAngle.png", "png")



sys.exit()
