using namespace RooFit;
using namespace RooStats;

#include <string>
#include <map>
#include "TGraph.h"
#include "TString.h"



static const Int_t NCAT = 3;
Double_t MMIN = 1000.;
Double_t MMAX = 2650;
std::string filePOSTfix="";
double analysisLumi = 12.9; // Luminosity you use in your analysis
double nEventsInSignalMC = 0.; //number of events in Signal MC sample
int iGraviton = 0;

bool inTheList = false;

const int np(6);
double masses[np] = { 1200.0,1400.0,1600.0,1800.0,2000.0,2500.0 } ;

double signalScaler=0;//analysisLumi/nEventsInSignalMC; // assume signal cross section on 1/fb

std::vector<string> cat_names;
std::vector<string> method_names;
TString mainCut("1");

void MakeDataCard_1Channel(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName, int iChan, TString signalname, int signalsample, double mass) {

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


  for (int imethod == 0; imethod < 3; imethod++){

    TString filename(cardDir+TString(fileBaseName)+Form("_%s",cat_names[iChan].c_str())+Form("_%s.txt",method_names[imethod].c_str()));
    ofstream outFile(filename);

    cout << "================================================================ signalScaler = " << signalScaler << endl;


    outFile << "# Fully Hadronic HH analysis" << endl;
    outFile << "imax 1" << endl;
    outFile << "kmax *" << endl;
    outFile << "---------------" << endl;

    outFile << Form("shapes data_obs bin_%s ", cat_names[iChan].c_str()) << wsDir+TString(fileBkgName)+".root" << Form(" w_all:data_obs_%s", cat_names[iChan].c_str()) << endl;

    if (imethod == 0) outFile << Form("shapes bkg_fit bin_%s ", cat_names[iChan].c_str()) <<  wsDir+TString(fileBkgName)+".root" << Form(" w_all:CMS_bkg_fit_%s", cat_names[iChan].c_str()) << endl;
    else outFile << Form("shapes bkg_fit bin_%s ", cat_names[iChan].c_str()) <<  wsDir+TString(fileBkgName)+".root" << Form(" w_all:CMS_bkg_fit_%s", cat_names[2].c_str()) << endl;
  

    outFile << Form("shapes HH bin_%s ", cat_names[iChan].c_str()) << wsDir+TString::Format("HH_%.0f_13TeV.root", mass) << Form(" w_all:HH_sig_%s", cat_names[iChan].c_str()) << endl;
    outFile << "---------------" << endl;
    outFile << Form("bin          bin_%s", cat_names[iChan].c_str()) << endl;
    outFile <<  "observation   "  <<  Form("%.10lg",data[iChan]->sumEntries()) << endl;
    outFile << "------------------------------" << endl;

    outFile << "bin                      "<< Form("bin_%s       bin_%s      ", cat_names[iChan].c_str(), cat_names[iChan].c_str()) << endl;
    outFile << "process                 HH     bkg_fit      " << endl;
    outFile << "process                 0        1          " << endl;
    if(signalname=="HH"){
      if (inTheList) outFile <<  "rate                      " 
			     << " " << signal[iChan]->sumEntries()*signalScaler << " " << 1 << endl;
      else {
	double val_eff;
	getEfficiency(mass, iChan, val_eff);
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


    outFile << "# Parametric shape uncertainties, entered by hand." << endl;
    outFile << Form("CMS_sig_p1_jes_13TeV    param   0.0   1.0   # dijet mass shift due to JES uncertainty") << endl;
    outFile << Form("CMS_sig_p2_jer_13TeV     param   0.0   1.0   # dijet mass resolution shift due to JER uncertainty") << endl;
 


    if (imethod == 0){
      outFile << Form("CMS_bkg_fit_%s_norm           flatParam  # Normalization uncertainty on background slope",cat_names[iChan].c_str()) << endl;
      outFile << Form("CMS_bkg_fit_slope1_eps_%s     flatParam  # Mean and absolute uncertainty on background slope",cat_names[iChan].c_str()) << endl;
      if (iChan > 0) outFile << Form("CMS_bkg_fit_slope2_eps_%s     flatParam  # Mean and absolute uncertainty on background slope",cat_names[iChan].c_str()) << endl;
    } else if (imethod == 1){
      outFile << Form("CMS_bkg_fit_%s_norm           flatParam  # Normalization uncertainty on background slope",cat_names[2].c_str()) << endl;
      outFile << Form("CMS_bkg_fit_slope1_eps_%s     param  0.0  1.0  # Mean and absolute uncertainty on background slope",cat_names[2].c_str()) << endl;
      outFile << Form("CMS_bkg_fit_slope2_eps_%s     param  0.0  1.0  # Mean and absolute uncertainty on background slope",cat_names[2].c_str()) << endl;
    } else if (imethod == 2){
      outFile << Form("CMS_bkg_fit_%s_norm           flatParam  # Normalization uncertainty on background slope",cat_names[iChan].c_str()) << endl;
      outFile << Form("CMS_bkg_fit_slope1_eps_%s     flatParam  # Mean and absolute uncertainty on background slope",cat_names[2].c_str()) << endl;
      outFile << Form("CMS_bkg_fit_slope2_eps_%s     flatParam  # Mean and absolute uncertainty on background slope",cat_names[2].c_str()) << endl;
    }


    outFile.close();

    cout << "Write data card in: " << filename << " file" << endl;

  }

  return;
}


void R2JJDatacardsMakerHH_13TeV(double mass, std::string postfix="", bool dobands=false, int graviton = 0, double nEvents = 50000.)
{
    filePOSTfix=postfix;
    if (postfix.find("Graviton") != string::npos) iGraviton = true;
    nEventsInSignalMC = nEvents;
    signalScaler=analysisLumi/nEventsInSignalMC;
    //style();

    for (int i = 0; i < np; i++){
      if ( fabs(mass - masses[i]) < 1e-5 ) {
	inTheList = true;
	break;
      }
    }

	TString signalname="HH";
	TString fileBaseName(TString::Format("%s_%.0f_13TeV", signalname.Data(), mass));

	
	
cat_names.push_back("4btag_cat0");
cat_names.push_back("3btag_cat1");
cat_names.push_back("2btag_cat2");
	
	// classical bump hunt
  method_names.push_back("bumphunt");
  // fit 2btag cateogry and use in 3,4 btag
  method_names.push_back("seqfrank");
  // fit at the same time 2 btag and 3 or 4 btag
  method_names.push_back("simfrank");
	
	TString fileBkgName("bkg_HH_13TeV");
	TString card_name("Xvv_reduced_models_Bkg_HH_13TeV.rs");

  
	TFile* tf1=TFile::Open(Form("makeDatacards/%s_%s.root",filePOSTfix.data(),fileBaseName.Data()));
	RooWorkspace* w= tf1->Get("HLFactory_ws");
	tf1->Close();
	
	int signalsample = 0;
	//cout << "CREATE DATACARD" << endl;
	MakeDataCard_1Channel(w, fileBaseName, fileBkgName, 0, signalname, 0, mass);
	MakeDataCard_1Channel(w, fileBaseName, fileBkgName, 1, signalname, 0, mass);
	MakeDataCard_1Channel(w, fileBaseName, fileBkgName, 2, signalname, 0, mass);


    //runfits(mass, 0, dobands);

}