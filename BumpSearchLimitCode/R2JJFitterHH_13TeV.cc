/** \macro H2GGFitter.cc
 *
 * $Id: R2JJFitter.cc,v 1.14 2013/06/27 10:18:27 hinzmann Exp $
 *
 * Software developed for the CMS Detector at LHC
 *
 *
 *  \author Serguei Ganjour - CEA/IRFU/SPP, Saclay
 *  \modified by Maxime Gouzevitch for the Dijet Bump Search - IPNL, Lyon 
 *
 * Macro is implementing the unbinned maximum-likelihood model for 
 * the Higgs to gamma gamma analysis. PDF model and RooDataSets 
 * are stored in the workspace which is feeded to  HiggsAnalysis/CombinedLimit tools:
 * 
 * http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/HiggsAnalysis/CombinedLimit
 * http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/HiggsAnalysis/CombinedLimit/data/lhc-hcg/cms-jj-1fb/
 * 
 * The analysis root trees produced in a simple format 
 *
 *     TFile file(filename,"RECREATE", "X->jj input tree for unbinned maximum-likelihood fit");
*     TTree* outTree  = new TTree("XTojj","X->jj input tree for unbinned maximum-likelihood fit");
 *     Float_t mass;
 *     Int_t CAT3;
 *     Float_t weight;
 *
 *     outTree->Branch("mass",&mass,"mass/F");
 *     outTree->Branch("weight",&weight,"weight/F");
 *     outTree->Branch("CAT4",&CAT4,"CAT4/I");
 *     {
 *       .............
 *       outTree->Fill();
 *     }
 *
 *     file.Write();
 *     file.Close();
 *     delete outTree;
 *
 * are used as input files. They have to be produced for 
 * data and Monte Carlo signal and background data sets 
 * after all analysis selections to be applied. It is recommended to put   
 * loose kinematical cuts on pt1 and pt2 (20 GeV) since further selections 
 * are possible based on RooDataSets. 
 * It is recommended to use Root 5.28/00 (CMSSW_4_1_3).
 *
 *
 */
// Loading:  .L H2GGFitter.cc
// Running:  runfits("jj120-shapes-combined-Unbinned.root")  
//                

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

const int np(6);
double masses[np] = { 1200.0,1400.0,1600.0,1800.0,2000.0,2500.0 } ;

double signalScaler=0;//analysisLumi/nEventsInSignalMC; // assume signal cross section on 1/fb
double scaleFactorHP=1;// already done on 1 GeV Histo level
double scaleFactorLP=1;// already done on 1 GeV Histo level

void AddSigData(RooWorkspace*, Float_t, int, std::vector<string>);
void AddBkgData(RooWorkspace*, std::vector<string>);
void SigModelFit(RooWorkspace*, Float_t, TString signalname, std::vector<string>);
void BkgModelFit(RooWorkspace*, Bool_t, std::vector<string>, RooFitResult** fitresults);
void MakePlots(RooWorkspace*, Float_t, RooFitResult** , TString signalname, std::vector<string>);
void MakeSigWS(RooWorkspace* w, const char* filename, TString signalname, std::vector<string>);
void MakeBkgWS(RooWorkspace* w, const char* filename, std::vector<string>);
void SetConstantParams(const RooArgSet* params);
void MakeDataCard_1Channel(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName, int iChan, TString signalname, int signalsample, std::vector<string> cat_names, double mass);

TString mainCut("1");

RooArgSet* defineVariables()
{
  // define variables of the input ntuple
  RooRealVar* mgg  = new RooRealVar("mgg","M(jet-jet)",MMIN,MMAX,"GeV");
  RooRealVar* evWeight   = new RooRealVar("evWeight","Reweightings",0.,10000.,"");
  RooRealVar* normWeight  = new RooRealVar("normWeight","Additionnal Weight",0.,10000000.,"");
  RooCategory* categories = new RooCategory("categories","event category 0") ;
  categories->defineType("4btag_cat0",0);
  categories->defineType("3btag_HPHP_cat1",1);
  categories->defineType("3btag_HPLP_cat2",2);
  RooArgSet* ntplVars = new RooArgSet(*mgg, *categories, *evWeight, *normWeight);
 
  return ntplVars;
}


void runfits(const Float_t mass=1600, int signalsample = 0, Bool_t dobands = false)
{

//******************************************************************//
//  Running mode  corresponds to the following cases
//         - full run set:
//         - create signal and background data sets 
//         - make and fit signal and background  models 
//         - write signal and background workspaces in root files
//         - write data card

//*******************************************************************//


  TString signalname;
  if (signalsample==0)
  { signalname="HH";
  }

  TString fileBaseName("CMS_jj_"+signalname+TString::Format("_%.0f_13TeV", mass));

  vector<string> cat_names;
  cat_names.push_back("CMS_jj_4btag_cat0");
  cat_names.push_back("CMS_jj_3btag_HPHP_cat1");
  cat_names.push_back("CMS_jj_3btag_HPLP_cat2");


  TString fileBkgName("CMS_jj_bkg_HH_13TeV");
  TString card_name("Xvv_reduced_models_Bkg_HH_13TeV.rs");
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();
  RooFitResult* fitresults[NCAT];

  w->var("mgg")->setMin(MMIN);
  w->var("mgg")->setMax(MMAX);

// Add data to the workspace

  cout << "CREATE SIGNAL" << endl;

  if (inTheList) AddSigData(w, mass, signalsample,cat_names);

  SigModelFit(w, mass, signalname, cat_names);
  MakeSigWS(w, fileBaseName, signalname,cat_names);

  // ===============================

  AddBkgData(w,cat_names);
  BkgModelFit(w, dobands,cat_names, fitresults, signalname, mass);
  MakeBkgWS(w, fileBkgName,cat_names);


  // =======================

  cout << "CREATE DATACARD" << endl;
  MakeDataCard_1Channel(w, fileBaseName, fileBkgName, 0, signalname, signalsample, cat_names, mass);
  MakeDataCard_1Channel(w, fileBaseName, fileBkgName, 1, signalname, signalsample, cat_names, mass);
  MakeDataCard_1Channel(w, fileBaseName, fileBkgName, 2, signalname, signalsample, cat_names, mass);

  cout << "MAKE PLOTS" << endl;
  
// Make plots for data and fit results
  MakePlots(w, mass, fitresults, signalname, cat_names);


  return;



}



void AddSigData(RooWorkspace* w, Float_t mass, int signalsample, std::vector<string> cat_names) {

  Int_t ncat = NCAT;
  TString inDir   = "./MiniTrees/Signal_HH_13TeV/";


//****************************//
// Signal Data Set
//****************************//

  // Variables
  RooArgSet* ntplVars = defineVariables();


//signal300_tree_radcut.root
  int iMass = abs(mass);       
  string signal(inDir.Data());
  signal = signal +"dijetHH_" + filePOSTfix + "" + string(Form("HHOUT%d_miniTree.root", iMass));
  
  cout << " ================================================================================= signal_c_str() = " << signal.c_str() << endl;

  TFile sigFile1(signal.c_str());

  TTree* sigTree1 = (TTree*) sigFile1.Get("TCVARS");

// common preselection cut

  sigTree1->Print();

//****************************//
// Signal  Data Set
//****************************//
// Create non scaled signal dataset composed with  different productions 
// according to their cross sections

  RooDataSet sigScaled("sigScaled","Signal",sigTree1,*ntplVars,mainCut, "normWeight");

  cout << "Print Signal Scaled" << endl;
  //  sigScaled.Print("v");

  RooDataSet* sigToFit[NCAT];
  for (int c = 0; c < ncat; ++c) {
    TString cut(TString::Format("categories==%d",c));
    TString name = TString::Format("Sig_%s",cat_names.at(c).c_str());

    sigToFit[c] =  (RooDataSet*) sigScaled.reduce(*w->var("mgg"),mainCut+TString::Format(" && categories==%d",c));
    w->import(*sigToFit[c],Rename(TString::Format("Sig_%s",cat_names.at(c).c_str())));
    cout << "Sum Entries = " << sigToFit[c]->sumEntries() << " isWeighted ? = " << sigToFit[c]->isWeighted() << endl;
    
  }

  cout << "End Making Sig Data" << endl;
          
}


void AddBkgData(RooWorkspace* w, std::vector<string> cat_names) {

  Int_t ncat = NCAT;
  TString inDir   = "./MiniTrees/Background_HH_13TeV/";

// common preselection cut


  Float_t minMassFit(MMIN),maxMassFit(MMAX); 


//****************************//
// CMS Data Set
//****************************//
// retrieve the data tree;
// no common preselection cut applied yet; 



  TString infile = inDir + "dijetHH_data_miniTree.root"; 
  if (filePOSTfix.find("subtr") != string::npos) 
    infile = inDir + "dijetHH_data_subtr_miniTree.root"; 

  TFile dataFile(infile.Data());

  TTree* dataTree     = (TTree*) dataFile.Get("TCVARS");

  // Variables
  RooArgSet* ntplVars = defineVariables();

  RooDataSet Data("Data","dataset",dataTree,*ntplVars,"","normWeight");

  
  RooDataSet* dataToFit[NCAT];
  for (int c = 0; c < ncat; ++c) {

    dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("mgg"),mainCut+TString::Format(" && categories==%d",c));
    w->import(*dataToFit[c],Rename(TString::Format("Data_%s",cat_names.at(c).c_str())));

    cout << "Sum Entries data = " << dataToFit[c]->sumEntries() << " isWeighted ? = " << dataToFit[c]->isWeighted() << endl;
  }

// Create full data set without categorization
  RooDataSet* data    = (RooDataSet*) Data.reduce(*w->var("mgg"),mainCut);
  w->import(*data, Rename("Data"));
  data->Print("v");

}



void SigModelFit(RooWorkspace* w, Float_t mass, TString signalname, std::vector<string> cat_names) {

  Int_t ncat = NCAT;
  Float_t MASS = mass - 50.;

//******************************************//
// Fit signal with model pdfs
//******************************************//
// retrieve pdfs and datasets from workspace to fit with pdf models

  cout << "Fit signal" << endl;

  RooDataSet* sigToFit[NCAT];
  RooAbsPdf* jjSig[NCAT];

  Float_t minMassFit(MMIN),maxMassFit(MMAX); 


// Fit Signal 


  for (int c = 0; c < ncat; ++c) {
    cout << "---------- category = " << c << endl;
    if (inTheList)
	sigToFit[c]   = (RooDataSet*) w->data(TString::Format("Sig_%s",cat_names.at(c).c_str()));
    
    jjSig[c]     = (RooAbsPdf*)  w->pdf(signalname+"_jj"+TString::Format("_%s",cat_names.at(c).c_str()));

    cerr << ("jj_"+signalname+TString::Format("_sig_m0_%s",cat_names.at(c).c_str())) << endl;
    ((RooRealVar*) w->var("jj_"+signalname+TString::Format("_sig_m0_%s",cat_names.at(c).c_str())))->setVal(MASS);
  
    cout << " Mass = " << MASS << endl;
      
    if (inTheList)
      jjSig[c]     ->fitTo(*sigToFit[c],Range(mass*0.8,mass*1.3),SumW2Error(kTRUE),RooFit::PrintEvalErrors(-1));
    else {
      double val_alpha, val_n, val_sigma, val_m0;
      getVal(mass, c, val_alpha, val_n, val_sigma, val_m0);
      ((RooRealVar*) w->var("jj_"+signalname+TString::Format("_sig_m0_%s",cat_names.at(c).c_str())))->setVal(val_m0);
      ((RooRealVar*) w->var("jj_"+signalname+TString::Format("_sig_sigma_%s",cat_names.at(c).c_str())))->setVal(val_sigma);
      ((RooRealVar*) w->var("jj_"+signalname+TString::Format("_sig_alpha_%s",cat_names.at(c).c_str())))->setVal(val_alpha);
      ((RooRealVar*) w->var("jj_"+signalname+TString::Format("_sig_n_%s",cat_names.at(c).c_str())))->setVal(val_n);
    }

    

    cout << " fitted " << endl;

// IMPORTANT: fix all pdf parameters to constant
    w->defineSet(TString::Format("SigPdfParam_%s",cat_names.at(c).c_str()), RooArgSet(*w->var("jj_"+signalname+TString::Format("_sig_m0_%s",cat_names.at(c).c_str())),
								   *w->var("jj_"+signalname+TString::Format("_sig_sigma_%s",cat_names.at(c).c_str())),
								   *w->var("jj_"+signalname+TString::Format("_sig_alpha_%s",cat_names.at(c).c_str())),
								   *w->var("jj_"+signalname+TString::Format("_sig_n_%s",cat_names.at(c).c_str())) 
										      //								   *w->var("jj_"+signalname+TString::Format("_sig_gsigma_%s",cat_names.at(c).c_str())),
										      //								   *w->var("jj_"+signalname+TString::Format("_sig_frac_%s",cat_names.at(c).c_str()))) 
										      ));

    cout << " defined " << endl;

    SetConstantParams(w->set(TString::Format("SigPdfParam_%s",cat_names.at(c).c_str())));
  }


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









void BkgModelFit(RooWorkspace* w, Bool_t dobands, std::vector<string> cat_names, RooFitResult** fitresult, TString signalname, double mass) {

  Int_t ncat = NCAT;

//******************************************//
// Fit background with model pdfs
//******************************************//

// retrieve pdfs and datasets from workspace to fit with pdf models

  cout << "Start background model fit" << endl;
  RooDataSet* data[NCAT];
  RooPlot* plotbkg_fit[NCAT];

// dobands and dosignal
  RooDataSet* signal[NCAT];
  RooAbsPdf*  jjSig[NCAT];


  Float_t minMassFit(MMIN),maxMassFit(MMAX); 

// Fit data with background pdf for data limit

  RooRealVar* mgg     = w->var("mgg");  
  mgg->setUnit("GeV");
  
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);


  for (int c = 0; c < ncat; ++c) {
    data[c]   = (RooDataSet*) w->data(TString::Format("Data_%s",cat_names.at(c).c_str()));
                    
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("p1mod_%s",cat_names.at(c).c_str()),"","@0",*w->var(TString::Format("bkg_fit_slope1_%s",cat_names.at(c).c_str())));
    RooFormulaVar *p1mod_clone = new RooFormulaVar(TString::Format("p1mod_clone_%s",cat_names.at(c).c_str()),"","@0",*w->var(TString::Format("bkg_fit_slope1_clone_%s",cat_names.at(c).c_str())));
    RooFormulaVar *p1mod_power = new RooFormulaVar(TString::Format("p1mod_power_%s",cat_names.at(c).c_str()),"","@0",*w->var(TString::Format("bkg_fit_slope1_power_%s",cat_names.at(c).c_str())));

     


    RooFormulaVar *sqrtS = new RooFormulaVar(TString::Format("sqrtS_%s",cat_names.at(c).c_str()),"","@0",*w->var("sqrtS"));
    RooFormulaVar *x = new RooFormulaVar(TString::Format("x_%s",cat_names.at(c).c_str()),"","@0/@1",RooArgList(*mgg, *sqrtS));


    // EXO-12-053 1-parameter function
    RooAbsPdf* bkg_fitTmp_2par = new RooGenericPdf(TString::Format("bkg_fit_2par_%s",cat_names.at(c).c_str()), "exp(-1*@1*@1*@0)", RooArgList(*x, *p1mod_clone));
    bkg_fitTmp_2par->fitTo(*data[c], Strategy(1),Minos(kFALSE), Range(minMassFit,maxMassFit),SumW2Error(kTRUE), Save(kTRUE),RooFit::PrintEvalErrors(-1));

    // EXO-12-053 1-parameter function
    RooAbsPdf* bkg_fitTmp_power = new RooGenericPdf(TString::Format("bkg_fit_power_%s",cat_names.at(c).c_str()), "TMath::Power(@0, -1*@1*@1)", RooArgList(*x, *p1mod_power));
    bkg_fitTmp_power->fitTo(*data[c], Strategy(1),Minos(kFALSE), Range(minMassFit,maxMassFit),SumW2Error(kTRUE), Save(kTRUE),RooFit::PrintEvalErrors(-1));

    RooAbsPdf* bkg_fitTmp = new RooGenericPdf(TString::Format("bkg_fit_%s",cat_names.at(c).c_str()), "exp(-1*@1*@1*@0)", RooArgList(*x, *p1mod));
    if (c > 0)  {
      RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("p2mod_%s",cat_names.at(c).c_str()),"","@0",*w->var(TString::Format("bkg_fit_slope2_%s",cat_names.at(c).c_str())));
      bkg_fitTmp = new RooGenericPdf(TString::Format("bkg_fit_%s",cat_names.at(c).c_str()), "exp(-1*@1*@1*@0/(1+@1*@1*@2*@0))", RooArgList(*x, *p1mod, *p2mod));
    }

    fitresult[c] =  bkg_fitTmp->fitTo(*data[c], Strategy(1),Minos(kFALSE), Range(minMassFit,maxMassFit),SumW2Error(kTRUE), Save(kTRUE),RooFit::PrintEvalErrors(-1));
    Double_t norm = data[c]->sumEntries();
    cout << "====================== norm = " << norm << endl;


    RooAbsReal* bkg_fitTmp2  = new RooRealVar(TString::Format("bkg_fit_%s_norm",cat_names.at(c).c_str()),"",norm,0.0,10000000);
    w->import(*bkg_fitTmp);
    w->import(*bkg_fitTmp2);


//************************************************//
// Plot jj background fit results per categories 
//************************************************//
// Plot Background Categories 
//****************************//

    TCanvas* ctmp = new TCanvas("ctmp","jj Background Categories",0,0,800,600);
    Int_t nBinsMass(33);
    plotbkg_fit[c] = mgg->frame(nBinsMass);
    plotbkg_fit[c]->SetXTitle("M^{red}_{jj} (GeV)");
    plotbkg_fit[c]->SetYTitle(Form("Events / %i GeV", (int) plotbkg_fit[c]->getFitRangeBinW()));

    //    plotbkg_fit[c]->SetAxisRange(0.1,70,"Y");

    data[c]->plotOn(plotbkg_fit[c],Invisible());//LineColor(kWhite),MarkerColor(kWhite));    

    bkg_fitTmp->plotOn(plotbkg_fit[c],LineColor(kRed),Range(minMassFit+25., maxMassFit-25.),NormRange("fitrange"),RooFit::PrintEvalErrors(-1)); 
    //  bkg_fitTmp_2par->plotOn(plotbkg_fit[c],LineColor(kMagenta),Range("fitrange"),NormRange("fitrange"),RooFit::PrintEvalErrors(-1)); 
    bkg_fitTmp_power->plotOn(plotbkg_fit[c],LineColor(kBlue),Range(minMassFit+25., maxMassFit-25.),NormRange("fitrange"),RooFit::PrintEvalErrors(-1)); 
    data[c]->plotOn(plotbkg_fit[c]);    


    plotbkg_fit[c]->Draw();  

    TLatex *latexLabel = new TLatex();
    latexLabel->SetTextSize(0.75 * ctmp->GetTopMargin());
    latexLabel->SetNDC();
    latexLabel->SetTextFont(42); // helvetica
    latexLabel->DrawLatex(0.78, 0.96, "2.7 fb^{-1} (13 TeV)");
    latexLabel->SetTextFont(61); // helvetica bold face
    latexLabel->DrawLatex(0.19, 0.88, "CMS");

    //    TLatex *latexLabelPrel = new TLatex();
    latexLabel->SetTextFont(52); // helvetica bold face
    latexLabel->DrawLatex(0.19, 0.83, "Preliminary");


    


    jjSig[c]       = (RooAbsPdf*)  w->pdf(signalname+"_jj"+TString::Format("_%s",cat_names.at(c).c_str()));

    double signalNorm =  0;
    if (inTheList) {
      signal[c]       = (RooDataSet*) w->data(TString::Format("Sig_%s",cat_names.at(c).c_str()));
      signalNorm = signal[c]->sumEntries()*signalScaler*20;
    } else {
      double val_eff;
      getEfficiency(mass, c, val_eff);
       signalNorm = val_eff*20;
    }



    jjSig[c]->plotOn(plotbkg_fit[c], LineColor(kMagenta), LineStyle(2), Normalization(signalNorm ,RooAbsReal::NumEvent)); 










    RooArgSet* set = new RooArgSet(*mgg);
    /*
    cout << " ================================ PRINTING ===================================" << endl;

  

      // ================================================================== data_obs_CMS
      //        cat_names.push_back("CMS_jj_4btag_cat0");

    data[c]->plotOn(plotbkg_fit[c],LineColor(kWhite),MarkerColor(kWhite));    

    bkg_fitTmp_2par->plotOn(plotbkg_fit[c],LineColor(kBlue),Range("fitrange"),NormRange("fitrange"),RooFit::PrintEvalErrors(-1)); 
    */

    double normalisation =  data[c]->sumEntries();
    const double alpha = 1 - 0.6827;

    double chi2 = 0, chi2_3par = 0;
    double rss0=0, rss_3par=0;
    for (int i = 0; i < (plotbkg_fit[c]->getHist(TString::Format("h_Data_%s",cat_names.at(c).c_str())))->GetN() ; i++)
      {
	Double_t x0,y0;
	(plotbkg_fit[c]->getHist(TString::Format("h_Data_%s",cat_names.at(c).c_str())))->GetPoint(i, x0, y0);
	int N = y0;
	double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
	double U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1);

	double xc, yc;
	double xlowErr = (plotbkg_fit[c]->getHist(TString::Format("h_Data_%s",cat_names.at(c).c_str())))->GetErrorXlow(i);
	double xhighErr = (plotbkg_fit[c]->getHist(TString::Format("h_Data_%s",cat_names.at(c).c_str())))->GetErrorXhigh(i);
	(plotbkg_fit[c]->getHist(TString::Format("h_Data_%s",cat_names.at(c).c_str())))->GetPoint(i, xc, yc);


	if (N > -1) (plotbkg_fit[c]->getHist(TString::Format("h_Data_%s",cat_names.at(c).c_str())))->SetPointError(i, 0, 0, N-L, U-N);
	else (plotbkg_fit[c]->getHist(TString::Format("h_Data_%s",cat_names.at(c).c_str())))->SetPointError(i, 0, 0, 0, 0);
	
	cout << "total integral = "  <<  data[c]->sumEntries();


	double xmin = xc-fabs(xlowErr), xmax = xc+fabs(xhighErr);
	
	mgg->setRange("A",xmin,xmax);

	RooAbsReal* intBin0 = bkg_fitTmp_2par->createIntegral(*set,*set,"A") ;

	RooAbsReal* intBin_3par = bkg_fitTmp->createIntegral(*set,*set,"A") ;
    
	double dintBin0 = intBin0->getVal();
	double dintBin_3par = intBin_3par->getVal();
	cout << "=================== Bin = " << i << " xmin = " << xmin << " xmax = " << xmax << " xlowErr = " << xlowErr << " xhighErr = " << xhighErr << " N-L = " << N-L << " U-N = " << U-N << " bin content = " << y0 << " intBin = " << dintBin0 << " unnormalised integral = " << dintBin0*normalisation << " yc = " << yc << endl;

	if (dintBin0*normalisation >= yc) chi2 += TMath::Power((dintBin0*normalisation - yc)/(U-N),2);
	else  chi2 += TMath::Power((dintBin0*normalisation - yc)/(N-L),2);

	if (dintBin_3par*normalisation >= yc) chi2_3par += TMath::Power((dintBin_3par*normalisation - yc)/(U-N),2);
	else  chi2_3par += TMath::Power((dintBin_3par*normalisation - yc)/(N-L),2);

	rss0 += TMath::Power((dintBin0*normalisation - yc),2);
	rss_3par += TMath::Power((dintBin_3par*normalisation - yc),2);

      }

    double p1_10 = 1;
    double p2_10 = (plotbkg_fit[c]->getHist(TString::Format("h_Data_%s",cat_names.at(c).c_str())))->GetN() - 3;
    double Ftest_10 = (rss0-rss_3par)/p1_10 / (rss_3par/p2_10);
    double good_CL23 =  1.-TMath::FDistI(Ftest_10,p1_10,p2_10);
    
    cout <<  " p1_10 = " << p1_10 << " p2_10 = " << p2_10 << " Ftest_23 = " << Ftest_10 << endl;

    cout << "chi2 = " << chi2 << endl;
    cout << "chi2_3par = " << chi2_3par << endl;
    cout << "F Test CL = " << good_CL23 << endl;




//********************************************************************************//

    if (dobands) {

      RooAbsPdf *cpdf; cpdf = bkg_fitTmp_2par;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotbkg_fit[c]->getObject(1));
      
      for (int i=1; i<(plotbkg_fit[c]->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotbkg_fit[c]->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotbkg_fit[c]->GetXaxis()->GetBinUpEdge(i);
	double center  = plotbkg_fit[c]->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);
	nlim->setVal(nombkg);
	mgg->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data[c]),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	// printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
	onesigma->SetPoint(i-1,center,nombkg);
	if (fabs(nlim->getErrorLo())>1e-5) onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	else onesigma->SetPointError(i-1,0.,0.,nlim->getErrorHi(),nlim->getErrorHi());

	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	// eventually if cl = 0.95 this is the usual 1.92!      
	
	
	minim.migrad();
	minim.minos(*nlim);
	
	twosigma->SetPoint(i-1,center,nombkg);
	if (fabs(nlim->getErrorLo())>1e-5) twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	else twosigma->SetPointError(i-1,0.,0.,nlim->getErrorHi(),nlim->getErrorHi());
	
	
	delete nll;
	delete epdf;
	
      }
      mgg->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma->SetLineColor(kYellow+1);
      twosigma->SetFillColor(kYellow+1);
      twosigma->SetMarkerColor(kYellow+1);
      twosigma->Draw("L3 SAME");
      
      onesigma->SetLineColor(kGreen+1);
      onesigma->SetFillColor(kGreen+1);
      onesigma->SetMarkerColor(kGreen+1);
      onesigma->Draw("L3 SAME");
      

      TLatex *lat  = new TLatex(MMIN+1000,10.,Form("#scale[1.0]{Exp. #chi^{2} = %.1f}",chi2));
      lat->SetTextSize(0.04);
      //    lat->Draw();
      TLatex *lat_3par  = new TLatex(MMIN+1000,6.,Form("#scale[1.0]{Lev. Exp. #chi^{2} = %.1f}",chi2_3par));
      lat_3par->SetTextSize(0.04);
      //    lat_3par->Draw();

      TLatex *lat_ftest  = new TLatex(MMIN+1000,20.,Form("#scale[1.0]{F Test 1 vs 2 par. CL = %.2f}",good_CL23));
      lat_ftest->SetTextSize(0.03);
      //   lat_ftest->Draw();

      plotbkg_fit[c]->SetTitle("");
      plotbkg_fit[c]->Draw("SAME"); 
      plotbkg_fit[c]->SetTitle("");      

 
 


      plotbkg_fit[c]->GetYaxis()->SetRangeUser(1.001e-1,70);

      ctmp->SetLogy();

          
      TLegend *legmc = new TLegend(0.47,0.60,0.92,0.8);
      TLegend *legmc2 = new TLegend(0.35,0.70,0.49,0.8);
      legmc->SetTextFont(62);
      legmc2->SetTextFont(62);
      

      legmc->AddEntry(plotbkg_fit[c]->getObject(3),"Data ","EP");
      legmc->AddEntry(plotbkg_fit[c]->getObject(1),"Background model","L");
      legmc->AddEntry(plotbkg_fit[c]->getObject(2),"Alternative background model","L");
      legmc->AddEntry(plotbkg_fit[c]->getObject(4),"G_{Bulk}, #sigma(M_{X}=1.6 TeV) = 20 fb","L");

      if(dobands)legmc2->AddEntry(onesigma,"68% CL","F");
      if(dobands)legmc2->AddEntry(twosigma,"95% CL","F"); // not...

      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc2->SetBorderSize(0);
      legmc2->SetFillStyle(0);
      legmc->Draw();
      legmc2->Draw();

      
      latexLabel->SetTextFont(42); // helvetica
      //    latexLabel->SetTextSize(0.6 * ctmp->GetTopMargin());
      latexLabel->SetTextColor(kBlue);
      latexLabel->DrawLatex(0.50, 0.88, "pp #rightarrow X #rightarrow HH  #rightarrow b#bar{b}b#bar{b}");
    
      if (c == 0)    latexLabel->DrawLatex(0.54, 0.83, " 4 b tag category");
      else if (c == 1)    latexLabel->DrawLatex(0.54, 0.83, " 3 b tag category");
      else if (c == 2)    latexLabel->DrawLatex(0.54, 0.83, " 2 b tag category");
      


      ctmp->Update();

      string out("plots/backgrounds");
      out = out + "" + filePOSTfix.c_str() + Form("_channel%d", c) + "_withband.pdf";
      ctmp->SaveAs(out.c_str());


      out = string("plots/backgrounds");
      out = out + "" + filePOSTfix.c_str() + Form("_channel%d", c) + "_withband_log.pdf";
      ctmp->SaveAs(out.c_str());

      out = string("plots/backgrounds");
      out = out + "" + filePOSTfix.c_str() + Form("_channel%d", c) + "_withband_log.root";
      ctmp->SaveAs(out.c_str());

      out = string("plots/backgrounds");
      out = out + "" + filePOSTfix.c_str() + Form("_channel%d", c) + "_withband_log.C";
      ctmp->SaveAs(out.c_str());

      out = string("plots/backgrounds");
      out = out + "" + filePOSTfix.c_str() + Form("_channel%d", c) + "_withband_parameters.root";


      TFile * output = new TFile(out.c_str(), "RECREATE");
      onesigma->Write("onesigma");
      twosigma->Write("twosigma");
      p1mod->Write();
      output->Close();
    }


 

  }


}


void SetConstantParams(const RooArgSet* params) {

  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
  }  

}

void MakePlots(RooWorkspace* w, Float_t mass, RooFitResult** fitresults, TString signalname, std::vector<string> cat_names) {

  cout << "Start plotting" << endl; 

  Int_t ncat = NCAT;

// retrieve data sets from the workspace
  RooDataSet* dataAll         = (RooDataSet*) w->data("Data");
  RooDataSet* signalAll       = (RooDataSet*) w->data("Sig");

  RooDataSet* data[9];  
  RooDataSet* signal[9];
  //  RooAbsPdf*  jjGaussSig[9];
  //  RooAbsPdf*  jjCBSig[9];
  RooAbsPdf*  jjSig[9];
  RooAbsPdf*  bkg_fit[9];  
//  RooAbsPdf*  bkg_fit2[9];  

  for (int c = 0; c < ncat; ++c) {
    data[c]         = (RooDataSet*) w->data(TString::Format("Data_%s",cat_names.at(c).c_str()));
//    signal[c]       = (RooDataSet*) w->data(TString::Format("Sig_%s",cat_names.at(c).c_str()));
    if (inTheList) signal[c]       = (RooDataSet*) w->data(TString::Format("Sig_%s",cat_names.at(c).c_str()));
    //    jjGaussSig[c]  = (RooAbsPdf*)  w->pdf(TString::Format("jjGaussSig_%s",cat_names.at(c).c_str()));
    //    jjCBSig[c]     = (RooAbsPdf*)  w->pdf(TString::Format("jjCBSig_%s",cat_names.at(c).c_str()));
    jjSig[c]       = (RooAbsPdf*)  w->pdf(signalname+"_jj"+TString::Format("_%s",cat_names.at(c).c_str()));
    bkg_fit[c]       = (RooAbsPdf*)  w->pdf(TString::Format("bkg_fit_%s",cat_names.at(c).c_str()));
//    bkg_fit2[c]      = (RooAbsPdf*)  w->pdf(TString::Format("bkg_fit2_%s",cat_names.at(c).c_str()));
  }

// retrieve mass observable from the workspace
  RooRealVar* mgg     = w->var("mgg");  
  mgg->setUnit("GeV");

// retrieve pdfs after the fits
// Signal Model

//  RooAbsPdf* jjGaussSigAll  = w->pdf("jjGaussSig"+signalname);
//  RooAbsPdf* jjCBSigAll     = w->pdf("jjCBSig"+signalname);
  RooAbsPdf* jjSigAll       = w->pdf(signalname+"_jj");

//  RooAbsPdf* bkg_fitAll       = w->pdf("bkg_fit");
  RooAbsPdf* bkg_fitAll       = w->pdf("bkg_fitAll");
  
  cout << "Progress plotting 1" << endl;
 
//****************************//
// Plot jj Fit results
//****************************//


  Float_t minMassFit(MMIN),maxMassFit(MMAX); 
  Float_t MASS(mass);

  Int_t nBinsMass(100);



  cout << "Progress plotting 2" << endl;

//********************************************//
// Plot jj signal fit results per categories 
//********************************************//
// Plot Signal Categories 
//****************************//

  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);

//  TCanvas* c2 = new TCanvas("c2","jj Categories",0,0,1000,1000);

//  c2->Divide(3,3);
  RooPlot* plotjj[9];
  for (int c = 0; c < ncat; ++c) {
    plotjj[c] = mgg->frame(Range(minMassFit,maxMassFit),Bins(nBinsMass));
    if (inTheList) signal[c]->plotOn(plotjj[c],LineColor(kWhite),MarkerColor(kWhite));    

    jjSig[c]  ->plotOn(plotjj[c]);
    //    jjSig[c]  ->plotOn(plotjj[c],Components("jjGaussSig"+signalname+TString::Format("_%s",cat_names.at(c).c_str())),LineStyle(kDashed),LineColor(kGreen),RooFit::PrintEvalErrors(-1));
    //   jjSig[c]  ->plotOn(plotjj[c],Components("jjCBSig"+signalname+TString::Format("_%s",cat_names.at(c).c_str())),LineStyle(kDashed),LineColor(kRed),RooFit::PrintEvalErrors(-1));
    

    //    jjSig[c]  ->paramOn(plotjj[c]);
    if (inTheList) signal[c]  ->plotOn(plotjj[c]);

        
    TCanvas* dummy = new TCanvas(Form("dummy_%d",c), "dummy",0, 0, 400, 400);
    TH1F *hist = new TH1F(Form("hist_%d",c), "hist", 400, minMassFit, maxMassFit);
 
    plotjj[c]->SetTitle("");      
    plotjj[c]->SetMinimum(0.0);
    plotjj[c]->SetMaximum(1.40*plotjj[c]->GetMaximum());
    plotjj[c]->GetXaxis()->SetTitle("m_{jj} [GeV]");

    TCanvas* ctmp_sig = new TCanvas(Form("ctmp_sig_%d",c),"jj Background Categories",0,0,500,500);
    plotjj[c]->Draw();  
//    hist->Draw("same");
    plotjj[c]->Print();

    plotjj[c]->Draw("SAME");  
    TLegend *legmc = new TLegend(0.62,0.75,0.92,0.9);
    legmc->AddEntry(plotjj[c]->getObject(2),"Simulation","LPE");
    legmc->AddEntry(plotjj[c]->getObject(1),"Parametric Model","L");
    //    legmc->AddEntry(plotjj[c]->getObject(3),"Crystal Ball component","L");
    //    legmc->AddEntry(plotjj[c]->getObject(2),"Gaussian Outliers","L");
    
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();    
    

    TLatex *lat  = new TLatex(minMassFit+1.5,0.85*plotjj[c]->GetMaximum(),"#scale[1.0]{CMS Preliminary}");
    lat->Draw();

    int iMass = abs(mass);

    //ctmp_sig->SaveAs("plots/sigmodel_"+signalname+TString::Format("%d_%s.png", iMass, cat_names.at(c).c_str()));
    string pdfout("plots/sigmodel_");
    cout << pdfout.c_str() << endl;
    pdfout = pdfout + "" + filePOSTfix.c_str() + "" + signalname.Data() + "" + Form("%d_%s.pdf", iMass, cat_names.at(c).c_str());
    cout << pdfout.c_str() << endl;

    string pngout("plots/sigmodel_");
    pngout = pngout + "" + filePOSTfix.c_str() + "" + signalname.Data() + "" + Form("%d_%s.pdf", iMass, cat_names.at(c).c_str());
    cout << pngout.c_str() << endl;

    ctmp_sig->SaveAs(pngout.c_str());
    ctmp_sig->SaveAs(pdfout.c_str());

  }


}


void MakeSigWS(RooWorkspace* w, const char* fileBaseName, TString signalname, std::vector<string> cat_names) {

  TString wsDir   = "workspaces/"+filePOSTfix;
  Int_t ncat = NCAT;


//********************************//
// Retrieve P.D.F.s
//********************************//

  RooAbsPdf* jjSigPdf[6];

// (1) import signal P.D.F.s

  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");


  for (int c = 0; c < ncat; ++c) {
    jjSigPdf[c] = (RooAbsPdf*)  w->pdf(signalname+"_jj"+TString::Format("_%s",cat_names.at(c).c_str()));
    wAll->import(*w->pdf(signalname+"_jj"+TString::Format("_%s",cat_names.at(c).c_str())));
  }

// (2) Systematics on energy scale and resolution

  wAll->factory("CMS_sig_p1_jes_13TeV[0.0,-5.0,5.0]");
  wAll->factory("CMS_jj_sig_p1_jes_13TeV[0.012,0.012,0.012]");
  wAll->factory("sum::CMS_sig_p1_jes_sum(1.0,prod::CMS_sig_p1_jes_prod(CMS_sig_p1_jes_13TeV, CMS_jj_sig_p1_jes_13TeV))");
    for (int c = 0; c < ncat; ++c) {
      wAll->factory("prod::CMS_jj_"+signalname+"_sig_m0_"+TString::Format("%s",cat_names.at(c).c_str())+"(jj_"+signalname+"_sig_m0_"+TString::Format("%s",cat_names.at(c).c_str())+", CMS_sig_p1_jes_sum)");
    }

// (3) Systematics on resolution: create new sigmas

    // apply JER resolution smearing from 
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
    // Assume ~7% smearing +- 3% for each jet. Divide by sqrt(2) for both

  wAll->factory("CMS_sig_p2_jer_13TeV[0.0,-5.0,5.0]");
  wAll->factory("CMS_jj_sig_p2_jer_13TeV[0.02,0.02,0.02]");
  wAll->factory("sum::CMS_sig_p2_jer_sum(1.00,prod::CMS_sig_p2_jer_prod(CMS_sig_p2_jer_13TeV, CMS_jj_sig_p2_jer_13TeV))");

    for (int c = 0; c < ncat; ++c) {
      wAll->factory("prod::CMS_jj_"+signalname+"_sig_sigma_"+TString::Format("%s",cat_names.at(c).c_str())+"(jj_"+signalname+"_sig_sigma_"+TString::Format("%s",cat_names.at(c).c_str())+", CMS_sig_p2_jer_sum)");
    }

    //    for (int c = 0; c < ncat; ++c) {
    //wAll->factory("prod::CMS_jj_"+signalname+"_sig_gsigma_"+TString::Format("%s",cat_names.at(c).c_str())+"(jj_"+signalname+"_sig_gsigma_"+TString::Format("%s",cat_names.at(c).c_str())+", CMS_sig_p2_jer_sum)");
    //    }

// (4) do reparametrization of signal
  for (int c = 0; c < ncat; ++c) {
    wAll->factory(
		  "EDIT::"+signalname+"_jj"+TString::Format("_sig_%s(",cat_names.at(c).c_str())+signalname+"_jj"+TString::Format("_%s,",cat_names.at(c).c_str()) +
		  " jj_"+signalname+TString::Format("_sig_m0_%s=CMS_jj_",cat_names.at(c).c_str())+signalname+TString::Format("_sig_m0_%s, ", cat_names.at(c).c_str()) +
		  " jj_"+signalname+TString::Format("_sig_sigma_%s=CMS_jj_",cat_names.at(c).c_str())+signalname+TString::Format("_sig_sigma_%s)", cat_names.at(c).c_str()));
  }

  TString filename(wsDir+TString(fileBaseName)+".root");
  wAll->writeToFile(filename);
  cout << "Write signal workspace in: " << filename << " file" << endl;

  return;
}


void MakeBkgWS(RooWorkspace* w, const char* fileBaseName, std::vector<string> cat_names) {

  TString wsDir   = "workspaces/"+filePOSTfix;
  Int_t ncat = NCAT;  


//********************************//
// Retrieve the datasets and PDFs
//********************************//

  RooDataSet* data[NCAT];
  RooExtendPdf* bkg_fitPdf[NCAT];

// (1) import everything

  cout << "Start importing everything" << endl;

  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");

  for (int c = 0; c < ncat; ++c) {
 
    cout << "For category " << c << endl;
    data[c]      = (RooDataSet*) w->data(TString::Format("Data_%s",cat_names.at(c).c_str()));
    ((RooRealVar*) data[c]->get()->find("mgg"))->setBins(MMAX-MMIN) ;
    RooDataHist* dataBinned = data[c]->binnedClone();
    bkg_fitPdf[c] = (RooExtendPdf*)  w->pdf(TString::Format("bkg_fit_%s",cat_names.at(c).c_str()));
    //   wAll->import(*data[c], Rename(TString::Format("data_obs_%s",cat_names.at(c).c_str())));
    wAll->import(*dataBinned, Rename(TString::Format("data_obs_%s",cat_names.at(c).c_str())));
   wAll->import(*w->pdf(TString::Format("bkg_fit_%s",cat_names.at(c).c_str())));
   wAll->import(*w->function(TString::Format("bkg_fit_%s_norm",cat_names.at(c).c_str())));

   double mean = (wAll->var(TString::Format("bkg_fit_%s_norm",cat_names.at(c).c_str())))->getVal();
   double min = (wAll->var(TString::Format("bkg_fit_%s_norm",cat_names.at(c).c_str())))->getMin();
   double max = (wAll->var(TString::Format("bkg_fit_%s_norm",cat_names.at(c).c_str())))->getMax();
   wAll->factory(TString::Format("CMS_bkg_fit_%s_norm[%g,%g,%g]", cat_names.at(c).c_str(), mean, min, max));

    mean = (wAll->var(TString::Format("bkg_fit_slope1_%s",cat_names.at(c).c_str())))->getVal();
    min = (wAll->var(TString::Format("bkg_fit_slope1_%s",cat_names.at(c).c_str())))->getMin();
    max = (wAll->var(TString::Format("bkg_fit_slope1_%s",cat_names.at(c).c_str())))->getMax();

   wAll->factory(TString::Format("CMS_bkg_fit_slope1_%s[%g,%g,%g]", cat_names.at(c).c_str(), mean, min, max));


   if (c > 0){
     mean = (wAll->var(TString::Format("bkg_fit_slope2_%s",cat_names.at(c).c_str())))->getVal();
     min = (wAll->var(TString::Format("bkg_fit_slope2_%s",cat_names.at(c).c_str())))->getMin();
     max = (wAll->var(TString::Format("bkg_fit_slope2_%s",cat_names.at(c).c_str())))->getMax();

     wAll->factory(TString::Format("CMS_bkg_fit_slope2_%s[%g,%g,%g]", cat_names.at(c).c_str(), mean, min, max));
   }

    cout << "Done For category " << c << endl;    
  }
  
  
  cout << "Imported" << endl;

// (2) do reparametrization of background

  for (int c = 0; c < ncat; ++c) {

    TString sFormat = 	    TString::Format("EDIT::CMS_bkg_fit_%s(bkg_fit_%s,",cat_names.at(c).c_str(),cat_names.at(c).c_str()) +
      TString::Format(" bkg_fit_%s_norm=CMS_bkg_fit_%s_norm,", cat_names.at(c).c_str(),cat_names.at(c).c_str())+
      TString::Format(" bkg_fit_slope1_%s=CMS_bkg_fit_slope1_%s,", cat_names.at(c).c_str(),cat_names.at(c).c_str())+
      TString::Format(" bkg_fit_slope2_%s=CMS_bkg_fit_slope2_%s)", cat_names.at(c).c_str(),cat_names.at(c).c_str());

    if (c == 0) 
      TString sFormat = 	    TString::Format("EDIT::CMS_bkg_fit_%s(bkg_fit_%s,",cat_names.at(c).c_str(),cat_names.at(c).c_str()) +
	TString::Format(" bkg_fit_%s_norm=CMS_bkg_fit_%s_norm,", cat_names.at(c).c_str(),cat_names.at(c).c_str())+
	TString::Format(" bkg_fit_slope1_%s=CMS_bkg_fit_slope1_%s)", cat_names.at(c).c_str(),cat_names.at(c).c_str());
    


    cout << sFormat.Data() << endl;

    wAll->factory(sFormat.Data());
  } 


  TString filename(wsDir+TString(fileBaseName)+".root");
  wAll->writeToFile(filename);
  cout << "Write background workspace in: " << filename << " file" << endl;

  std::cout << "observation ";
  for (int c = 0; c < ncat; ++c) {
    std::cout << "  " << (wAll->data(TString::Format("data_obs_%s",cat_names.at(c).c_str())))->sumEntries();
    (wAll->data(TString::Format("data_obs_%s",cat_names.at(c).c_str())))->Print();
  }
  std::cout << std::endl;
  
  for (int c = 0; c < ncat; ++c) {
    printf("CMS_bkg_fit_slope1_%s  param  %.4f  %.3f   # Mean and absolute uncertainty on background slope\n",
	   cat_names.at(c).c_str(), (wAll->var(TString::Format("CMS_bkg_fit_slope1_%s",cat_names.at(c).c_str())))->getVal(), 10.);
    if (c > 0) printf("CMS_bkg_fit_slope2_%s  param  %.4f  %.3f   # Mean and absolute uncertainty on background slope\n",
	   cat_names.at(c).c_str(), (wAll->var(TString::Format("CMS_bkg_fit_slope2_%s",cat_names.at(c).c_str())))->getVal(), 10.);
  }

  cout << "BKG WS DONE" << endl;

  return;
}















void MakeDataCard_1Channel(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName, int iChan, TString signalname, int signalsample, std::vector<string> cat_names, double mass) {

  TString cardDir = "datacards/"+filePOSTfix;
  Int_t ncat = NCAT;
  TString wsDir   = "../workspaces/"+filePOSTfix;
//**********************//
// Retrieve the datasets
//**********************//

  cout << "Start retrieving dataset" << endl;

  RooDataSet* data[NCAT];
  RooDataSet* signal[NCAT];
  for (int c = 0; c < NCAT; ++c) {
    data[c]        = (RooDataSet*) w->data(TString::Format("Data_%s",cat_names.at(c).c_str()));
    if (inTheList) signal[c] = (RooDataSet*) w->data(TString::Format("Sig_%s",cat_names.at(c).c_str()));
  }

//*****************************//
// Print Expected event yields
//*****************************//

  cout << "======== Expected Events Number =====================" << endl;  
  cout << "#Events data:        " <<  w->data("Data")->sumEntries()  << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data %s:   ",cat_names.at(c).c_str()) << data[c]->sumEntries()  << endl;
  }
  cout << ".........Expected Signal ............................" << endl;  

  for (int c = 0; c < ncat; ++c) {
    if (inTheList)
      cout << TString::Format("#Events Signal %s: ",cat_names.at(c).c_str()) << signal[c]->sumEntries() << endl;
    else {
      double val_eff;
      getEfficiency(mass, iChan, val_eff);
      cout << TString::Format("#Events Signal %s: ", cat_names.at(c).c_str()) << val_eff << endl;
    }
  }
  cout << "====================================================" << endl;  

//*************************//
// Print Data Crd int file
//*************************//


  TString filename(cardDir+TString(fileBaseName)+Form("_%s.txt",cat_names[iChan].c_str()));
  ofstream outFile(filename);

  cout << "================================================================ signalScaler = " << signalScaler << endl;

  // Pythia HP+HP
  //if(((signalsample==0))&&(iChan==0))
  //    scaleFactor*=(scaleFactorHP*scaleFactorHP);
  // Pythia HP+LP
  //if(((signalsample==0))&&(iChan==1))
  //    scaleFactor*=(scaleFactorHP*scaleFactorLP);

  outFile << "# Fully Hadronic HH analysis" << endl;
  outFile << "imax 1" << endl;
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;

  outFile << Form("shapes data_obs %s ", cat_names[iChan].c_str()) << wsDir+TString(fileBkgName)+".root" << Form(" w_all:data_obs_%s", cat_names[iChan].c_str()) << endl;
  outFile << Form("shapes bkg_fit_jj %s ", cat_names[iChan].c_str()) <<  wsDir+TString(fileBkgName)+".root" << Form(" w_all:CMS_bkg_fit_%s", cat_names[iChan].c_str()) << endl;
  outFile << Form("shapes HH_jj %s ", cat_names[iChan].c_str()) << wsDir+TString::Format("CMS_jj_HH_%.0f_13TeV.root", mass) << Form(" w_all:HH_jj_sig_%s", cat_names[iChan].c_str()) << endl;
  outFile << "---------------" << endl;
  outFile << Form("bin          %s", cat_names[iChan].c_str()) << endl;
  outFile <<  "observation   "  <<  Form("%.10lg",data[iChan]->sumEntries()) << endl;
  outFile << "------------------------------" << endl;

  outFile << "bin                      "<< Form("%s       %s      ", cat_names[iChan].c_str(), cat_names[iChan].c_str()) << endl;
  outFile << "process                 HH_jj     bkg_fit_jj     " << endl;
  outFile << "process                 0        1          " << endl;
  if(signalname=="HH"){
    if (inTheList) outFile <<  "rate                      " 
			   << " " << signal[iChan]->sumEntries()*signalScaler << " " << 1 << endl;
    else {
      double val_eff;
      getEfficiency(mass, iChan, val_eff)
      outFile <<  "rate                      " 
			   << " " << val_eff*signalScaler << " " << 1 << endl;
    }
  }

  outFile << "--------------------------------" << endl;
  outFile << "# signal scaled by " << signalScaler << " to a cross section of 1 fb and also scale factor of " << 1 << " are applied." << endl;
  
  outFile << "lumi_13TeV       lnN  1.027      - " << endl;
  outFile << "CMS_pu_13TeV              lnN  1.020      - # pileup impact of W mass tag" << endl;
  outFile << "CMS_mass_res_j_13TeV    lnN  1.02       - # JEs and JER uncertainty on H mass tag" << endl;
  outFile << "CMS_jet_mass_13TeV     lnN  1.10       - # differenec between H tag and W tag efficiencies" << endl;
  outFile << "CMS_PDF_Scales      lnN  1.02       - # selection efficiency" << endl;
  if(iChan==0)
    outFile << "CMS_htag_13TeV    lnN  1.17       - # btag efficiency" << endl;
  else if (iChan == 1 || iChan == 2)
    outFile << "CMS_htag_13TeV    lnN  0.95       - # btag efficiency" << endl;
    
  
  outFile << "CMS_eff_htag_tau21_sf_13TeV     lnN  1.056/0.946       - # tau21 efficiency" << endl;
  /*
  if (iChan == 0)
    outFile << "CMS_eff_htag_tau21_sf_13TeV     lnN  1.27/0.76       - # tau21 efficiency" << endl;
  else if (iChan == 1)
    outFile << "CMS_eff_htag_tau21_sf_13TeV   lnN  1.27/0.76       - # tau21 efficiency" << endl;
  else if(iChan == 2)
    outFile << "CMS_eff_htag_tau21_sf_13TeV     lnN  1.27/0.76       - # tau21 efficiency" << endl;
  */
  // HPHP 3btag: ((1.03+0.13)^2 - (1.03)^2) / 1.03^2 :  1.27/0.76
  // HPLP 3btag: ((1.03+0.13)(0.88+0.49) - (1.03)*0.88) / 1.03*0.88 :  1.75/0.38


  /*  
  if((iChan==0)||(iChan==3)){
  outFile << "CMS_eff_vtag_tau21_sf         lnN  1.15       - # tau21 efficiency" << endl;
  } else {
  // anti-correlated the high purity (1.076*1.076) and low purity (0.54*1.076) categories
  outFile << "CMS_eff_vtag_tau21_sf         lnN  0.58      - # tau21 efficiency" << endl;
  }
  */
  //  outFile << "CMS_scale_j         lnN  1.120 	   - # jet energy scale" << endl;
  //  outFile << "CMS_res_j         lnN  1.040	- # jet energy resolution" << endl;


  outFile << "# Parametric shape uncertainties, entered by hand." << endl;
  outFile << Form("CMS_sig_p1_jes_13TeV    param   0.0   1.0   # dijet mass shift due to JES uncertainty") << endl;
  outFile << Form("CMS_sig_p2_jer_13TeV     param   0.0   1.0   # dijet mass resolution shift due to JER uncertainty") << endl;
 
  outFile << Form("CMS_bkg_fit_%s_norm           flatParam  # Normalization uncertainty on background slope",cat_names[iChan].c_str()) << endl;

  outFile << Form("CMS_bkg_fit_slope1_%s         flatParam  # Mean and absolute uncertainty on background slope",cat_names[iChan].c_str()) << endl;

  if (iChan > 0) outFile << Form("CMS_bkg_fit_slope2_%s         flatParam  # Mean and absolute uncertainty on background levelled parameter",cat_names[iChan].c_str()) << endl;

  outFile.close();

  cout << "Write data card in: " << filename << " file" << endl;

  return;
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


void R2JJFitterHH_13TeV(double mass, std::string postfix="", bool dobands=false, int graviton = 0, double nEvents = 50000.)
{
    filePOSTfix=postfix;
    if (postfix.find("Graviton") != string::npos) iGraviton = true;
    nEventsInSignalMC = nEvents;
    signalScaler=analysisLumi/nEventsInSignalMC;
    style();

    for (int i = 0; i < np; i++){
      if ( fabs(mass - masses[i]) < 1e-5 ) {
	inTheList = true;
	break;
      }
    }



    runfits(mass, 0, dobands);

}
