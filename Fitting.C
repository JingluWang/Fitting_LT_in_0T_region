#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include "TMath.h"
#include "TH1D.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include "TMinuit.h"
#include <cstdlib>
#include <fstream>
#include <stdio.h>
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <TChain.h>
#include <TFile.h>
#include "TMath.h"
#include <TStyle.h>
#include "TLorentzVector.h"
#include <TApplication.h>
#include <TComplex.h>
#include "TEnv.h"
#include <TGraph.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include "TTree.h"
#include <TLegend.h>
#include <THStack.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TPad.h>
#include <TF1.h>

#include "GlobalV.h"

using namespace std;

double Bin_LLinclu[20];
double Bin_LTinclu[20];
double Bin_TLinclu[20];
double Bin_TTinclu[20];
double BinError_LLinclu[20];
double BinError_LTinclu[20];
double BinError_TLinclu[20];
double BinError_TTinclu[20];

double Bin_LL1[20];
double Bin_LL2[20];
double Bin_LL3[20];
double Bin_LL4[20];

double Bin_LT1[20];
double Bin_LT2[20];
double Bin_LT3[20];
double Bin_LT4[20]; 

double Bin_TL1[20];
double Bin_TL2[20];
double Bin_TL3[20];
double Bin_TL4[20];

double Bin_TT1[20];
double Bin_TT2[20];
double Bin_TT3[20];
double Bin_TT4[20];



double Alpha[5];
double FitError[5];

string ToString(double num)
{
  string res;
  stringstream ss;
  ss << num;
  ss >> res;
  return res;
}

/*****    Source file     *****/
TString file_inclu = "Sample/allYears_364253_Sherpa_222_NNPDF30NNLO_lllv_myOutput.root";

//emu + tau
TString SAMPLE = "emutau";

TString file_1 = "Sample/Pol_00_total.root"; 
TString file_2 = "Sample/Pol_0T_total.root";
TString file_3 = "Sample/Pol_T0_total.root";
TString file_4 = "Sample/Pol_TT_total.root";

/*
//emu
TString SAMPLE = "emu";

TString file_1 = "Sample/Pol_00.root";                                                  
TString file_2 = "Sample/Pol_0T.root";
TString file_3 = "Sample/Pol_T0.root";
TString file_4 = "Sample/Pol_TT.root";
*/
/******************************/

TH1D *Hist_LLInclusive = new TH1D("Hist_LLInclusive","",12,-3,3);
TH1D *Hist_LTInclusive = new TH1D("Hist_LTInclusive","",12,-3,3);
TH1D *Hist_TLInclusive = new TH1D("Hist_TLInclusive","",12,-3,3);
TH1D *Hist_TTInclusive = new TH1D("Hist_TTInclusive","",12,-3,3);

TH1D *Hist_LL1 = new TH1D("Hist_LL1","",12,-3,3);
TH1D *Hist_LL2 = new TH1D("Hist_LL2","",12,-3,3);
TH1D *Hist_LL3 = new TH1D("Hist_LL3","",12,-3,3);
TH1D *Hist_LL4 = new TH1D("Hist_LL4","",12,-3,3);

TH1D *Hist_LT1 = new TH1D("Hist_LT1","",12,-3,3);
TH1D *Hist_LT2 = new TH1D("Hist_LT2","",12,-3,3);
TH1D *Hist_LT3 = new TH1D("Hist_LT3","",12,-3,3);
TH1D *Hist_LT4 = new TH1D("Hist_LT4","",12,-3,3);

TH1D *Hist_TL1 = new TH1D("Hist_TL1","",12,-3,3);
TH1D *Hist_TL2 = new TH1D("Hist_TL2","",12,-3,3);
TH1D *Hist_TL3 = new TH1D("Hist_TL3","",12,-3,3);
TH1D *Hist_TL4 = new TH1D("Hist_TL4","",12,-3,3);

TH1D *Hist_TT1 = new TH1D("Hist_TT1","",12,-3,3);
TH1D *Hist_TT2 = new TH1D("Hist_TT2","",12,-3,3);
TH1D *Hist_TT3 = new TH1D("Hist_TT3","",12,-3,3);
TH1D *Hist_TT4 = new TH1D("Hist_TT4","",12,-3,3);

void FillHist(TString dir, TH1D *HIST1, TH1D *HIST2, TH1D *HIST3, TH1D *HIST4)
{
  Int_t passWZInclusive;
  Float_t Pt_Z;
  Float_t Pt_WZ;
  Float_t DY_WZ;
  Float_t CosThetaLepW;
  Float_t CosThetaLepZ;
  Float_t TotalWeightNoKFactor;

  TFile *f = new TFile(dir);
  TTree *tree = (TTree*)f->Get("nominal");
  tree->SetBranchAddress("passWZInclusive", &passWZInclusive);
  tree->SetBranchAddress("Pt_Z", &Pt_Z);
  tree->SetBranchAddress("Pt_WZ", &Pt_WZ);
  tree->SetBranchAddress("DY_WZ", &DY_WZ);
  tree->SetBranchAddress("CosThetaLepW", &CosThetaLepW);
  tree->SetBranchAddress("CosThetaLepZ", &CosThetaLepZ);
  tree->SetBranchAddress("TotalWeightNoKFactor", &TotalWeightNoKFactor);

  Int_t Nentries = (Int_t)tree->GetEntries();

  for(Int_t LoopNumber=0; LoopNumber<Nentries; LoopNumber++)
  {
    tree->GetEntry(LoopNumber);
    if(Pt_Z>200&&Pt_WZ<70)//Signal region 
    {
      if(abs(CosThetaLepW)<0.4&&abs(CosThetaLepZ)<0.4) HIST1->Fill(DY_WZ, TotalWeightNoKFactor*passWZInclusive);
      if(abs(CosThetaLepW)>0.4&&abs(CosThetaLepZ)<0.4) HIST2->Fill(DY_WZ, TotalWeightNoKFactor*passWZInclusive);
      if(abs(CosThetaLepW)<0.4&&abs(CosThetaLepZ)>0.4) HIST3->Fill(DY_WZ, TotalWeightNoKFactor*passWZInclusive);
      if(abs(CosThetaLepW)>0.4&&abs(CosThetaLepZ)>0.4) HIST4->Fill(DY_WZ, TotalWeightNoKFactor*passWZInclusive);
    }
  }
}


void GetFitValue()
{
  FillHist(file_inclu, Hist_LLInclusive, Hist_LTInclusive, Hist_TLInclusive, Hist_TTInclusive);
  FillHist(file_1, Hist_LL1, Hist_LT1, Hist_TL1, Hist_TT1); 
  FillHist(file_2, Hist_LL2, Hist_LT2, Hist_TL2, Hist_TT2);
  FillHist(file_3, Hist_LL3, Hist_LT3, Hist_TL3, Hist_TT3);
  FillHist(file_4, Hist_LL4, Hist_LT4, Hist_TL4, Hist_TT4);

  for(int i=1; i<13; i++)
  {
    Bin_LLinclu[i] = Hist_LLInclusive->GetBinContent(i);
    Bin_LTinclu[i] = Hist_LTInclusive->GetBinContent(i);
    Bin_TLinclu[i] = Hist_TLInclusive->GetBinContent(i);
    Bin_TTinclu[i] = Hist_TTInclusive->GetBinContent(i);
    BinError_LLinclu[i] = Hist_LLInclusive->GetBinError(i);
    BinError_LTinclu[i] = Hist_LTInclusive->GetBinError(i);
    BinError_TLinclu[i] = Hist_TLInclusive->GetBinError(i);
    BinError_TTinclu[i] = Hist_TTInclusive->GetBinError(i);

    Bin_LL1[i] = Hist_LL1->GetBinContent(i);
    Bin_LL2[i] = Hist_LL2->GetBinContent(i);
    Bin_LL3[i] = Hist_LL3->GetBinContent(i);
    Bin_LL4[i] = Hist_LL4->GetBinContent(i);
    
    Bin_LT1[i] = Hist_LT1->GetBinContent(i);
    Bin_LT2[i] = Hist_LT2->GetBinContent(i);
    Bin_LT3[i] = Hist_LT3->GetBinContent(i);
    Bin_LT4[i] = Hist_LT4->GetBinContent(i); 

    Bin_TL1[i] = Hist_TL1->GetBinContent(i);
    Bin_TL2[i] = Hist_TL2->GetBinContent(i);
    Bin_TL3[i] = Hist_TL3->GetBinContent(i);
    Bin_TL4[i] = Hist_TL4->GetBinContent(i);

    Bin_TT1[i] = Hist_TT1->GetBinContent(i);
    Bin_TT2[i] = Hist_TT2->GetBinContent(i);
    Bin_TT3[i] = Hist_TT3->GetBinContent(i);
    Bin_TT4[i] = Hist_TT4->GetBinContent(i);
  }
}


void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  double t0, t1, t2, t3, Sum = 0;

  for(int i=1; i<13; i++)
  {
    t0 = (par[1]*Bin_LL1[i]+par[2]*Bin_LL2[i]+par[3]*Bin_LL3[i]+par[4]*Bin_LL4[i]-Bin_LLinclu[i])/sqrt(Bin_LLinclu[i]);
    t0 *= t0;
    t1 = (par[1]*Bin_LT1[i]+par[2]*Bin_LT2[i]+par[3]*Bin_LT3[i]+par[4]*Bin_LT4[i]-Bin_LTinclu[i])/sqrt(Bin_LTinclu[i]);    
    t1 *= t1;
    t2 = (par[1]*Bin_TL1[i]+par[2]*Bin_TL2[i]+par[3]*Bin_TL3[i]+par[4]*Bin_TL4[i]-Bin_TLinclu[i])/sqrt(Bin_TLinclu[i]);
    t2 *= t2;
    t3 = (par[1]*Bin_TT1[i]+par[2]*Bin_TT2[i]+par[3]*Bin_TT3[i]+par[4]*Bin_TT4[i]-Bin_TTinclu[i])/sqrt(Bin_TTinclu[i]);
    t3 *= t3;
    Sum = Sum + t0 + t1 + t2 + t3;
  }

  f=Sum;

  cout<<f<<endl;
}
 

void ToFit()
{
  GetFitValue();

  TMinuit *fit;
  Double_t arglist[10];
  arglist[0] = 0.01;
  Int_t iflag = 0;
  fit = new TMinuit(4);
  fit->SetFCN(fcn);

  fit->mncler();  

  for(int i = 1; i < 5; i++)
  {
    TString str = ToString(double(i));
    fit->mnparm(i,"Alpha"+str,1.0,0.5,0.0,3.0,iflag);
  }
 
  fit->SetMaxIterations(10000);
  fit->SetPrintLevel(3);

  fit->SetErrorDef(1.0);
  arglist[0] = 1.0;
  fit->mnexcm("SET ERR", arglist, 1, iflag);
  

  arglist[0] = 1.0;
  arglist[1] = 0.01;

  fit->mnexcm("SIMPLEX", arglist, 1, iflag);
  fit->mnexcm("MIGRAD", arglist, 2, iflag);
  fit->mnexcm("HESSE", arglist, 1, iflag);

  fit->mnimpr(); 

  for(int i = 1;i < 5; i++)
  {
    fit->GetParameter(i,Alpha[i],FitError[i]);
  }

}


void DrawPlot(double *FIT_INCLU, double *BIN_INCLU, double *BINERROR_INCLU, double *POL1, double *POL2, double *POL3, double *POL4, TString LATEX, TString SAVE)
{
  TCanvas *c = new TCanvas("c","",1000,750);
  c->SetMargin(0.13,0.04,0.25,0.03);

  TH1D *hist_fit = new TH1D("hist_fit","",12,-3,3);
  TH1D *hist_data = new TH1D("hist_data","",12,-3,3);
  TH1D *hist_D = new TH1D("hist_D","",12,-3,3);
  TH1D *hist_1 = new TH1D("hist_1","",12,-3,3);
  TH1D *hist_2 = new TH1D("hist_2","",12,-3,3);
  TH1D *hist_3 = new TH1D("hist_3","",12,-3,3);
  TH1D *hist_4 = new TH1D("hist_4","",12,-3,3);

  double mmmax = FIT_INCLU[1];
  for(int i=1; i<13; i++)
  {
    hist_fit->SetBinContent(i,FIT_INCLU[i]);

    hist_data->SetBinContent(i,BIN_INCLU[i]);
    hist_data->SetBinError(i,BINERROR_INCLU[i]);

    hist_D->SetBinContent(i,BIN_INCLU[i]);
    hist_D->SetBinError(i,BINERROR_INCLU[i]);
    
    hist_1->SetBinContent(i,Alpha[1]*POL1[i]);
    hist_2->SetBinContent(i,Alpha[2]*POL2[i]);
    hist_3->SetBinContent(i,Alpha[3]*POL3[i]);
    hist_4->SetBinContent(i,Alpha[4]*POL4[i]);

    if(FIT_INCLU[i]>mmmax) mmmax = FIT_INCLU[i];
  }
  double mmax = 2*mmmax;

  hist_data->SetDirectory(0);
  hist_data->SetLineColor(kRed);
  hist_data->SetStats(0);
  hist_data->SetLineStyle(1);
  hist_data->SetLineWidth(2);
  hist_data->SetMarkerColor(2);
  hist_data->SetMarkerStyle(20);
  hist_data->SetMarkerSize(0.5);

  hist_fit->SetDirectory(0);
  hist_fit->SetLineColor(kBlue);
  hist_fit->SetStats(0);
  hist_fit->SetTitle(0);
  hist_fit->SetLineStyle(1);
  hist_fit->SetLineWidth(3);

  hist_1->SetDirectory(0);
  hist_1->SetFillColor(kViolet+2);
  hist_1->SetStats(0);
  hist_1->SetTitle(0);
  hist_1->SetLineWidth(0);

  hist_2->SetDirectory(0);
  hist_2->SetFillColor(kAzure+6);
  hist_2->SetStats(0);
  hist_2->SetTitle(0);
  hist_2->SetLineWidth(0);

  hist_3->SetDirectory(0);
  hist_3->SetFillColor(kOrange-7);
  hist_3->SetStats(0);
  hist_3->SetTitle(0);
  hist_3->SetLineWidth(0);

  hist_4->SetDirectory(0);
  hist_4->SetFillColor(kOrange-2);
  hist_4->SetStats(0);
  hist_4->SetTitle(0);
  hist_4->SetLineWidth(0);


  THStack *hs = new THStack("hs","");
  hs->Add(hist_1);
  hs->Add(hist_2);
  hs->Add(hist_3);
  hs->Add(hist_4);


  hist_fit->SetMaximum(mmax);
  hist_fit->SetMinimum(0);

  hist_fit->GetXaxis()->SetTitle("DY_WZ");
  hist_fit->GetXaxis()->SetLabelSize(0.03);
  hist_fit->GetXaxis()->SetLabelFont(72);
  hist_fit->GetXaxis()->SetTitleSize(30);
  hist_fit->GetXaxis()->SetTitleFont(73);

  hist_fit->GetYaxis()->SetTitle("Events");
  hist_fit->GetYaxis()->SetLabelSize(0.03);
  hist_fit->GetYaxis()->SetLabelFont(72);
  hist_fit->GetYaxis()->SetTitleSize(30);
  hist_fit->GetYaxis()->SetTitleFont(73);


  c->cd();
  hist_fit->Draw("HIST");
  hs->Draw("HISTsame");
  hist_fit->Draw("HISTsame");
  hist_data->Draw("E0same");

  TLegend *legend = new TLegend(0.98, 0.72, 0.68, 0.97);
  legend->AddEntry(hist_data,"Signal");
  legend->AddEntry(hist_fit,"Post-fit");
  legend->AddEntry(hist_1,"WLZL");
  legend->AddEntry(hist_2,"WLZT");
  legend->AddEntry(hist_3,"WTZL");
  legend->AddEntry(hist_4,"WTZT");
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetLineColor(0);
  legend->SetLineWidth(0);
  legend->SetTextSize(0.03);
  legend->SetTextFont(72);
  legend->Draw("same");

  TLatex *text = new TLatex(0.25,0.9,LATEX);
  text->SetNDC();
  text->SetTextSize(0.04);
  text->SetTextFont(72);
  text->Draw("same");

  TPad *pad = new TPad("pad", "pad", 0, 0, 1, 0.245);
  pad->SetTopMargin(1);
  pad->SetBottomMargin(0.38);
  pad->SetLeftMargin(0.13);
  pad->SetRightMargin(0.04);
  c->cd();
  pad->Draw("same");
  pad->cd();

  hist_D->Divide(hist_fit);

  hist_D->SetLineColor(kBlack);
  hist_D->SetLineStyle(1);
  hist_D->SetLineWidth(2);
  hist_D->SetStats(0);
  hist_D->SetMarkerColor(1);
  hist_D->SetMarkerStyle(20);
  hist_D->SetMarkerSize(0.5);

  hist_D->GetXaxis()->SetTitle("DY_WZ");
  hist_D->GetXaxis()->SetTitleOffset(3.3);
  hist_D->GetYaxis()->SetTitle("Signal/Post");

  hist_D->GetYaxis()->SetLabelSize(0.13);
  hist_D->GetYaxis()->SetLabelFont(72);
  hist_D->GetXaxis()->SetLabelSize(0.13);
  hist_D->GetXaxis()->SetLabelFont(72);

  hist_D->GetXaxis()->SetTitleSize(30);
  hist_D->GetXaxis()->SetTitleFont(73);
  hist_D->GetYaxis()->SetTitleSize(30);
  hist_D->GetYaxis()->SetTitleFont(73);

  hist_D->GetYaxis()->SetNdivisions(4);
  hist_D->GetYaxis()->SetRangeUser(0,2);

  hist_D->Draw("E0");

  TF1 *referline = new TF1("referline", "1", -3, 3);
  referline->Draw("same");

  c->SaveAs("Plots/"+SAMPLE+"/"+SAVE+".png");
}


int main()
{
  ToFit();
  cout << "ToFit clear" << endl;  

  /**********************/
  /*****   Output   *****/
  /**********************/
 
  for(int i = 1; i < 5; i++)
  {
    cout <<  i << " : " << setw(20) << Alpha[i] << "  Error :  " << FitError[i]  << endl;
  }
 
  /************************/
  /*****   Plotting   *****/
  /************************/

  double Fit00[20], Fit0T[20], FitT0[20], FitTT[20];
  for(int i=1; i<13; i++)
  {
    Fit00[i] = Alpha[1]*Bin_LL1[i]+Alpha[2]*Bin_LL2[i]+Alpha[3]*Bin_LL3[i]+Alpha[4]*Bin_LL4[i];
    Fit0T[i] = Alpha[1]*Bin_LT1[i]+Alpha[2]*Bin_LT2[i]+Alpha[3]*Bin_LT3[i]+Alpha[4]*Bin_LT4[i];
    FitT0[i] = Alpha[1]*Bin_TL1[i]+Alpha[2]*Bin_TL2[i]+Alpha[3]*Bin_TL3[i]+Alpha[4]*Bin_TL4[i];
    FitTT[i] = Alpha[1]*Bin_TT1[i]+Alpha[2]*Bin_TT2[i]+Alpha[3]*Bin_TT3[i]+Alpha[4]*Bin_TT4[i];
  }

  DrawPlot(Fit00, Bin_LLinclu, BinError_LLinclu, Bin_LL1, Bin_LL2, Bin_LL3, Bin_LL4, "region 00", "Fit00");
  DrawPlot(Fit0T, Bin_LTinclu, BinError_LTinclu, Bin_LT1, Bin_LT2, Bin_LT3, Bin_LT4, "region 0T", "Fit0T");
  DrawPlot(FitT0, Bin_TLinclu, BinError_TLinclu, Bin_TL1, Bin_TL2, Bin_TL3, Bin_TL4, "region T0", "FitT0");
  DrawPlot(FitTT, Bin_TTinclu, BinError_TTinclu, Bin_TT1, Bin_TT2, Bin_TT3, Bin_TT4, "region TT", "FitTT");

  return 0;
}
