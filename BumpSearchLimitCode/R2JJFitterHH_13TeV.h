/*
#include <cstring>
#include <cerrno>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <unistd.h>
#include <errno.h>
#include <iomanip>
// ROOT headers
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TIterator.h"

#include "TLatex.h"
#include "TString.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"

// RooFit headers
#include "RooAbsPdf.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooMsgService.h"
#include "RooNLLVar.h"
#include "RooPlot.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TStyle.h"

// RooStats headers
#include "RooStats/HLFactory.h"

#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAbsData.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooProduct.h"
#include "RooExtendPdf.h"
#include "RooBernstein.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooCmdArg.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
*/
//#include "HiggsCSandWidth.h"
//#include "HiggsCSandWidth.cc"
//#include "RooPower.h"

//#include "extrap.cc"

#include <string>
#include <map>

using namespace RooFit;
using namespace RooStats ;


static const Int_t NCAT = 3;
Double_t MMIN = 1000.;
Double_t MMAX = 2650;
std::string filePOSTfix="";
double analysisLumi = 2.69; // Luminosity you use in your analysis
double nEventsInSignalMC = 0.; //number of events in Signal MC sample
int iGraviton = 0;

bool inTheList = false;

double signalScaler=0;//analysisLumi/nEventsInSignalMC; // assume signal cross section on 1/fb
int iGraviton = 0;

bool inTheList = false;

const int np(6);
double masses[np] = { 1200.0,1400.0,1600.0,1800.0,2000.0,2500.0 } ;

double signalScaler=0;//analysisLumi/nEventsInSignalMC; // assume signal cross section on 1/fb
double scaleFactorHP=1;// already done on 1 GeV Histo level
double scaleFactorLP=1;// already done on 1 GeV Histo level

std::vector<string> cat_names;
std::vector<string> method_names;

TString mainCut("1");

RooArgSet* defineVariables()
{
  // define variables of the input ntuple
  RooRealVar* mgg  = new RooRealVar("mgg","M(jet-jet)",MMIN,MMAX,"GeV");
  RooRealVar* evWeight   = new RooRealVar("evWeight","Reweightings",0.,10000.,"");
  RooRealVar* normWeight  = new RooRealVar("normWeight","Additionnal Weight",0.,10000000.,"");
  RooCategory* categories = new RooCategory("categories","event category 0") ;
  categories->defineType("4btag_cat0",0);
  categories->defineType("3btag_cat1",1);
  categories->defineType("2btag_cat2",2);
  RooArgSet* ntplVars = new RooArgSet(*mgg, *categories, *evWeight, *normWeight);
 
  return ntplVars;
}





void getVal(double mjj, int ichannel, double& alpha, double& n, double& sigma, double& mean) {

  std::map<string, TGraph*> gr;

  if (!iGraviton) {
    if (ichannel == 0) {
      double alpha0[np] = { 2.99994285083,1.46588196934,1.32019487083,1.24987586528,1.166365396,1.13823997853 } ; 
      double n0[np] = { 126.460030265,132.776297381,132.018630891,121.374740807,144.337572861,144.868422011 } ; 
      double sigma0[np] = { 45.7033847362,51.1543513626,58.3791610108,62.8435240971,66.7409003893,82.8829569711 } ; 
      double mean0[np] = { 1141.97051876,1328.64384766,1516.84908692,1706.22166888,1892.21279472,2359.89977439 } ; 
      gr["alpha"+ichannel] = new TGraph(np, masses, alpha0); 
      gr["n"+ichannel] = new TGraph(np, masses, n0); 
      gr["sigma"+ichannel] = new TGraph(np, masses, sigma0); 
      gr["mean"+ichannel] = new TGraph(np, masses, mean0); 
    } 
    else if (ichannel == 1) {
      double alpha1[np] = { 2.99992869473,1.23572089917,1.38704596759,1.43790885737,1.47418666312,1.11961838575 } ; 
      double n1[np] = { 127.58280731,134.559662569,126.077137231,138.005844913,137.205261077,142.861915389 } ; 
      double sigma1[np] = { 49.6014843743,51.8530982105,59.9136530318,67.5958048863,75.9686650749,84.6657692385 } ; 
      double mean1[np] = { 1138.7518106,1329.27922894,1512.35492927,1702.76282505,1890.47142027,2363.48968532 } ; 
      gr["alpha"+ichannel] = new TGraph(np, masses, alpha1); 
      gr["n"+ichannel] = new TGraph(np, masses, n1); 
      gr["sigma"+ichannel] = new TGraph(np, masses, sigma1); 
      gr["mean"+ichannel] = new TGraph(np, masses, mean1); 
    }
    else if (ichannel == 2) {
      double alpha2[np] = { 2.99992869473,1.23572089917,1.38704596759,1.43790885737,1.47418666312,1.11961838575 } ; 
      double n2[np] = { 127.58280731,134.559662569,126.077137231,138.005844913,137.205261077,142.861915389 } ; 
      double sigma2[np] = { 49.6014843743,51.8530982105,59.9136530318,67.5958048863,75.9686650749,84.6657692385 } ; 
      double mean2[np] = { 1138.7518106,1329.27922894,1512.35492927,1702.76282505,1890.47142027,2363.48968532 } ; 
      gr["alpha"+ichannel] = new TGraph(np, masses, alpha2); 
      gr["n"+ichannel] = new TGraph(np, masses, n2); 
      gr["sigma"+ichannel] = new TGraph(np, masses, sigma2); 
      gr["mean"+ichannel] = new TGraph(np, masses, mean2); 
    }
    else std::cout << "Error: check channel\n"; 
  }
  else if (iGraviton) {
    if (ichannel == 0) {
      double alpha0[np] = { 2.9999986789,1.32495487888,1.23799526302,1.16954277446,1.21623722849,1.25469474452 } ; 
      double n0[np] = { 137.838630302,126.095947032,138.035633381,146.013900331,132.830584043,145.834851327 } ; 
      double sigma0[np] = { 45.6843979549,51.0787006165,57.3445271952,60.7986544462,68.0860896875,85.7909735223 } ; 
      double mean0[np] = { 1141.24023757,1330.25308396,1517.83383081,1703.26522902,1890.69263049,2359.9586752 } ; 
      gr["alpha"+ichannel] = new TGraph(np, masses, alpha0); 
      gr["n"+ichannel] = new TGraph(np, masses, n0); 
      gr["sigma"+ichannel] = new TGraph(np, masses, sigma0); 
      gr["mean"+ichannel] = new TGraph(np, masses, mean0); 
    } 
    else if (ichannel == 1) {
      double alpha1[np] = { 1.97782922403,1.41896083125,1.2617805834,1.34539670468,1.31276457424,1.32940541309 } ; 
      double n1[np] = { 141.139758152,117.887448804,129.516555916,128.842102297,128.751760353,131.820871684 } ; 
      double sigma1[np] = { 47.2594027067,53.1273966322,60.8086392999,63.8498759854,73.2388044796,87.1485163363 } ; 
      double mean1[np] = { 1139.44448398,1326.1356103,1515.97716079,1701.67066795,1889.79779805,2359.87112233 } ; 
      gr["alpha"+ichannel] = new TGraph(np, masses, alpha1); 
      gr["n"+ichannel] = new TGraph(np, masses, n1); 
      gr["sigma"+ichannel] = new TGraph(np, masses, sigma1); 
      gr["mean"+ichannel] = new TGraph(np, masses, mean1); 
    }
    else if (ichannel == 2) {
      double alpha2[np] = { 1.97782922403,1.41896083125,1.2617805834,1.34539670468,1.31276457424,1.32940541309 } ; 
      double n2[np] = { 141.139758152,117.887448804,129.516555916,128.842102297,128.751760353,131.820871684 } ; 
      double sigma2[np] = { 47.2594027067,53.1273966322,60.8086392999,63.8498759854,73.2388044796,87.1485163363 } ; 
      double mean2[np] = { 1139.44448398,1326.1356103,1515.97716079,1701.67066795,1889.79779805,2359.87112233 } ; 
      gr["alpha"+ichannel] = new TGraph(np, masses, alpha2); 
      gr["n"+ichannel] = new TGraph(np, masses, n2); 
      gr["sigma"+ichannel] = new TGraph(np, masses, sigma2); 
      gr["mean"+ichannel] = new TGraph(np, masses, mean2); 
    }

    else std::cout << "Error: check channel\n"; 
  }
  else std::cout << "Error: check sample\n"; 

  alpha = gr["alpha"+ichannel] -> Eval(mjj) ; 
  n = gr["n"+ichannel] -> Eval(mjj) ; 
  sigma = gr["sigma"+ichannel] -> Eval(mjj) ; 
  mean = gr["mean"+ichannel] -> Eval(mjj) ; 

  return;

}






void getEfficiency(double mjj, int ichannel, double& eff) {

  std::map<std::string, TGraph*> gr;

  if (!iGraviton) {
    if (ichannel == 0){
      double eff0[np] = {0.154, 0.1561, 0.1470, 0.1396, 0.1287, 0.1054};
      gr["eff"+ichannel] = new TGraph(np, masses, eff0); 
    } 
    else if (ichannel == 1) {
      double eff1[np] = {0.124, 0.1268, 0.1253, 0.1233, 0.1238, 0.1245};
      gr["eff"+ichannel] = new TGraph(np, masses, eff1); 
    } else if (ichannel == 2) {
      double eff2[np] = {0.03, 0.03, 0.035, 0.036, 0.038, 0.041};
      gr["eff"+ichannel] = new TGraph(np, masses, eff2); 
    }
    else std::cout << "Error: check channel\n"; 
  }
  else if (iGraviton) {
    if (ichannel == 0){
      double eff0[np] = {0.2384, 0.2286, 0.2244, 0.2171, 0.2033, 0.1882};
      gr["eff"+ichannel] = new TGraph(np, masses, eff0); 
    } 
    else if (ichannel == 1) {
      double eff1[np] = {0.1871, 0.1939, 0.1902, 0.2011, 0.1910, 0.1681};
      gr["eff"+ichannel] = new TGraph(np, masses, eff1); 
    } else if (ichannel == 2) {
      double eff2[np] = {0.05, 0.054, 0.053, 0.055, 0.06, 0.067};
      gr["eff"+ichannel] = new TGraph(np, masses, eff2); 
    }
    else std::cout << "Error: check channel\n"; 
  }
  else std::cout << "Error: check sample\n"; 

  eff = gr["eff"+ichannel] -> Eval(mjj) / signalScaler ; 

  return;

}




void SetConstantParams(const RooArgSet* params) {

  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
  }  

}






void style(){
  TStyle *defaultStyle = new TStyle("defaultStyle","Default Style");
  defaultStyle->SetOptStat(0000);
  defaultStyle->SetOptFit(000); 
  defaultStyle->SetPalette(1);
  /////// pad ////////////
  defaultStyle->SetPadBorderMode(1);
  defaultStyle->SetPadBorderSize(1);
  defaultStyle->SetPadColor(0);
  defaultStyle->SetPadTopMargin(0.05);
  defaultStyle->SetPadBottomMargin(0.13);
  defaultStyle->SetPadLeftMargin(0.15);
  defaultStyle->SetPadRightMargin(0.04);
  /////// canvas /////////
  defaultStyle->SetCanvasBorderMode(0);
  defaultStyle->SetCanvasColor(0);
  defaultStyle->SetCanvasDefH(600);
  defaultStyle->SetCanvasDefW(600);
  /////// frame //////////
  defaultStyle->SetFrameBorderMode(0);
  defaultStyle->SetFrameBorderSize(1);
  defaultStyle->SetFrameFillColor(0); 
  defaultStyle->SetFrameLineColor(1);
  /////// label //////////
  defaultStyle->SetLabelOffset(0.005,"XY");
  defaultStyle->SetLabelSize(0.05,"XY");
  defaultStyle->SetLabelFont(42,"XY");
  /////// title //////////
  defaultStyle->SetTitleOffset(1.1,"X");
  defaultStyle->SetTitleSize(0.01,"X");
  defaultStyle->SetTitleOffset(1.10,"Y");
  defaultStyle->SetTitleSize(0.05,"Y");
  defaultStyle->SetTitleFont(42, "XYZ");
  /////// various ////////
  // defaultStyle->SetNdivisions(505,"Y");
  defaultStyle->SetLegendBorderSize(0); // For the axis titles:

  defaultStyle->SetTitleColor(1, "XYZ");
  defaultStyle->SetTitleFont(42, "XYZ");
  defaultStyle->SetTitleSize(0.06, "XYZ");
  defaultStyle->SetTitleXOffset(0.9);
  defaultStyle->SetTitleYOffset(1.05);
  
  // For the axis labels:
  defaultStyle->SetLabelColor(1, "XYZ");
  defaultStyle->SetLabelFont(42, "XYZ");
  defaultStyle->SetLabelOffset(0.007, "XYZ");
  defaultStyle->SetLabelSize(0.04, "XYZ");

  // For the axis:
    defaultStyle->SetAxisColor(1, "XYZ");
    defaultStyle->SetStripDecimals(kTRUE);
    defaultStyle->SetTickLength(0.03, "XYZ");
    defaultStyle->SetNdivisions(310, "XYZ");
    defaultStyle->SetPadTickX(1);
    defaultStyle->SetPadTickY(1);
    defaultStyle->cd();
  return;
}

