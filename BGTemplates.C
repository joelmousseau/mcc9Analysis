#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <math.h>
#include <TF1.h>
#include <TPaletteAxis.h>
#include <TAxis.h>
#include <TLatex.h>
#include <vector>
#include <TLegend.h>
#include <TColor.h>
#include <stdio.h>

void MakeTemplatePlots(){
  const int nTemplates = 6;
  const int extTemplateInt = nTemplates-1;
  const int pBin = 256;
  
  std::vector<int> good_colors;
  good_colors.clear();
  good_colors.push_back( kBlue+2 );
  good_colors.push_back( kGreen+3 ); 
  good_colors.push_back( kYellow+2 );
  good_colors.push_back( kRed+2 );
  good_colors.push_back( kGray+2 );
  good_colors.push_back( kMagenta+2 );
  good_colors.push_back( kOrange-3 );
  good_colors.push_back( kPink+2 );
  good_colors.push_back( kCyan+2 );
  
  std::vector<std::string> channels;
  channels.push_back("CCQE");
  channels.push_back("RES");
  channels.push_back("DIS");
  channels.push_back("2p2h");
  channels.push_back("NC / Other");
  channels.push_back("Ext");
  
  THStack *stack = new THStack("hs","Template Stack");
  
  TH1F *pBinOne = new TH1F("pBinOne", "", 1, 0, pBin);
  pBinOne->SetBinContent(1, 50.0);
  pBinOne->SetLineWidth(3);
  pBinOne->SetLineStyle(kDashed);
  TH1F *pBinTwo = new TH1F("pBinTwo", "", 1, pBin, 2*pBin);
  pBinTwo->SetBinContent(1, 50.0);
  pBinTwo->SetLineWidth(3);
  pBinTwo->SetLineStyle(kDashed);
  TH1F *pBinThree = new TH1F("pBinThree", "", 1, 2*pBin+2, 3*pBin);
  pBinThree->SetBinContent(1, 50.0);
  pBinThree->SetLineWidth(3);
  pBinThree->SetLineStyle(kDashed);
  TH1F *pBinFour = new TH1F("pBinFour", "", 1, 3*pBin, 4*pBin);
  pBinFour->SetBinContent(1, 50.0);
  pBinFour->SetLineWidth(3);
  pBinFour->SetLineStyle(kDashed);
  
  //Get the Histogram Files
  TFile *histoFile= new TFile("Templates.root", "READ");
  TH1F *Templates[nTemplates];
  TH1F *weights[nTemplates];
  gStyle->SetOptStat(0);
  TLegend *legend = new TLegend(0.8, 0.9, 0.9, 0.7);
  TH1F *data = (TH1F*) histoFile->Get("Data");
  data->SetMarkerColor(kBlack);
  data->SetLineColor(kBlack);
  data->SetMarkerStyle(kFullDotLarge);
  
  for(int i=0; i < nTemplates; ++i){
    Templates[i] = (TH1F*) histoFile->Get(Form("template%d", i));
    Templates[i]->SetLineColor(good_colors[i]);
     Templates[i]->SetLineWidth(2);
     Templates[i]->SetTitle("Flattened Templates");
     legend->AddEntry(Templates[i], channels[i].c_str(), "l");
    stack->Add(Templates[i]);
    //std::cout << "Template " << i << " Bin 7: " << Templates[i]->GetBinContent(7) << std::endl;
    weights[i] = (TH1F*) histoFile->Get(Form("templateWeight%d", i));
  }
  legend->AddEntry(data, "Data", "p");
  TCanvas *c1 = new TCanvas("c1", " ", 1600, 800);
  stack->Draw("");

  data->Draw("Ep same");
  pBinOne->Draw("hist same");
  pBinTwo->Draw("hist same");
  pBinThree->Draw("hist same");
  pBinFour->Draw("hist same");
  legend->Draw("same");

  c1->Print("StackedWithData.png", "png");
  c1->Update();
  c1->Clear();

  TLegend *legend2 = new TLegend(0.8, 0.9, 0.9, 0.7);
  for(int i=0; i < nTemplates; ++i){
     legend2->AddEntry(Templates[i], channels[i].c_str(), "l");
     if(i == 0)
       Templates[i]->Draw("hist");
     else
       Templates[i]->Draw("hist same");
  }

  pBinOne->Draw("hist same");
  pBinTwo->Draw("hist same");
  pBinThree->Draw("hist same");
  pBinFour->Draw("hist same");
  legend2->Draw("same");
  c1->Print("TemplatesOnly.png", "png");
  

}



void DoFits(){  
  const int nTemplates = 6;
  const int extTemplateInt = nTemplates-1;
  double intialFrac[nTemplates];
  double total  = 0.0;
  double extInt = 0.0;

  /*mc_hist.Scale(self.fractions["mc"] *
                          data_hist.Integral() /
                          mc_hist.Integral())
 */                         
  
  std::vector<int> good_colors;
  good_colors.clear();
  good_colors.push_back( kBlack );
  good_colors.push_back( kRed+2 );
  good_colors.push_back( kBlue+2 );
  good_colors.push_back( kYellow+2 );
  good_colors.push_back( kMagenta+2 );
  good_colors.push_back( kGreen+3 );
  good_colors.push_back( kOrange-3 );
  good_colors.push_back( kPink+2 );
  good_colors.push_back( kCyan+2 );
  
  //Get the Histogram Files
  //TFile *histoFile= new TFile("Templates.root", "READ");
  TFile *histoFile= new TFile("Templates.root", "READ");
  TH1F *Templates[nTemplates];
  TH1F *weights[nTemplates];
  TObjArray *stack = new TObjArray(nTemplates);
  TH1F *Data   = (TH1F*) histoFile->Get("Data");
  TH1F *PreFit = (TH1F*) Data->Clone("preFit");
  PreFit->Reset();
    
  for(int i=0; i < nTemplates; ++i){
    Templates[i] = (TH1F*) histoFile->Get(Form("template%d", i));
    weights[i]   = (TH1F*) histoFile->Get(Form("templateWeight%d", i));
    
    stack->Add(Templates[i]);
    total += Templates[i]->Integral();
    
    intialFrac[i] = 0.0;
    
    if(i == extTemplateInt)
      extInt = Templates[i]->Integral();
    
    
  }
    


/*
  for(int i=0; i < nTemplates; ++i){

    if(i == extTemplateInt)
     extInt = Templates[i]->Integral();
     
    if(i == 1){
     stack->Add(disResSum);
     total += Templates[i]->Integral();
    }
    
    else if (i == 2){
      total += Templates[i]->Integral();
    }
    
    else{ 
      
      stack->Add(Templates[i]);
      total += Templates[i]->Integral();      
    }

    intialFrac[i] = 0.0;
 }*/
  
 const int templateBins = PreFit->GetNbinsX();
 printf("Number of bins : %d \n", templateBins );

 double extFraction = extInt / total;
 printf("Ext fraction: %e \n", extFraction);
 
 for (int bin = 1; bin <= templateBins; ++bin){
  double sum = 0.0;

  for(int i=0; i < nTemplates; ++i){
    //std::cout << weights[i]->GetBinContent(bin) << std::endl;
    double weight = 1.0;
    //weight = weights[i]->GetBinContent(bin);
    sum += (Templates[i]->GetBinContent(bin)*weight);
  }  
  PreFit->SetBinContent(bin, sum);
   /* code */
 }


  TFractionFitter *fit = new TFractionFitter(Data, stack);
  TVirtualFitter *vFit = (TVirtualFitter*) fit->GetFitter();
 // vFit->SetParameter(0, "CCQE", intialFrac[0], 0.01, 0.0, 1.0);
  
  for(int i=0; i < nTemplates; ++i){
    fit->SetWeight(i, weights[i]);
    intialFrac[i] = Templates[i]->Integral() / total;
  }


const double minFrac = 0.5;
const double maxFrac = (1/minFrac);

printf("Constraint: %e - %e \n", minFrac*intialFrac[0], maxFrac*intialFrac[0]);

  fit->SetRangeX(1, 1024);
  
  fit->Constrain(0, minFrac*intialFrac[0], maxFrac*intialFrac[0]); //CCQE
  fit->Constrain(1, minFrac*intialFrac[1], maxFrac*intialFrac[1]); //CCRES
  fit->Constrain(2, minFrac*intialFrac[2], maxFrac*intialFrac[2]); //DIS
  fit->Constrain(3, minFrac*intialFrac[3], maxFrac*intialFrac[3]); //2p2h
  fit->Constrain(4, minFrac*intialFrac[4], maxFrac*intialFrac[4]); //NC
  fit->Constrain(5, minFrac*intialFrac[5], maxFrac*intialFrac[5]); //EXT
  
  //fit->Constrain(6, 0.0, 1.0);

  int status = fit->Fit();
  TH1F *result = (TH1F*) fit->GetPlot();
  TH1 *ccqe = fit->GetMCPrediction(4);
  ccqe->SetMarkerColor(kRed);
  ccqe->SetMarkerStyle(kFullDotLarge);
  Templates[5]->SetLineColor(kBlack);
  double QEScale = 0.0;
  double ScaleErr = 0.0;
  double fracTotal = 0.0;
  for(int param = 0; param < nTemplates; ++param){
     double scale    = 0.0;
     double scaleErr = 0.0;
     double dataMCRatio = 1.0;
     
     fit->GetResult(param, scale, scaleErr);
     fracTotal += scale;
     
     dataMCRatio = (Data->Integral() / Templates[param]->Integral());
     //dataMCRatio = (Data->Integral() / fit->GetMCPrediction(param)->Integral());
     printf("Intial Fraction %d : %e Fitted Fracion : %e Ratio : %.2f MC SF : %.2f \n", param,intialFrac[param], scale, (scale/intialFrac[param]), scale*dataMCRatio);
     
     /*
     if(param == 1){
       dataMCRatio = (Data->Integral() / (Templates[1]->Integral() + Templates[2]->Integral() ) );
       printf("Intial Fraction %d : %e Fitted Fracion : %e Ratio : %.2f MC SF : %.2f \n", param,intialFrac[1]+intialFrac[2], scale, (scale/(intialFrac[1]+intialFrac[2])), scale*dataMCRatio);
       
     }
     
     else{
       dataMCRatio = (Data->Integral() / Templates[param+1]->Integral());
       printf("Intial Fraction %d : %e Fitted Fracion : %e Ratio : %.2f MC SF : %.2f \n", param,intialFrac[param+1], scale, (scale/intialFrac[param+1]), scale*dataMCRatio);
     }*/

  }

  fit->GetResult(0, QEScale, ScaleErr);
  //ccqe->Draw("E0");
  //Templates[5]->Draw("hist same");

  Data->SetMarkerColor(kBlack);
  Data->SetLineColor(kBlack);
  Data->SetMarkerStyle(kFullDotLarge);
  Data->Rebin(16);
  Data->Draw("Ep");


  result->SetLineColor(kRed+1);
  result->SetLineWidth(2);
  result->Rebin(16);
  result->Draw("same");

  PreFit->SetLineColor(kBlue+1);
  PreFit->SetLineWidth(2);
  PreFit->Rebin(16);
  PreFit->Draw("same");


  printf(" Chi2: %.2f \n", fit->GetChisquare() );
  printf("Sum: %e \n", fracTotal);
  //std::cout << "Bin 7: " << result->GetBinContent(7) << " Scale: " << QEScale << " Chi2 " << fit->GetChisquare() << std::endl;

}//end of comparison plots
