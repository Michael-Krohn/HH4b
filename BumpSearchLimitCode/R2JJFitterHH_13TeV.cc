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
 *
 *
 */


#include "R2JJFitterHH_13TeV.h"

#include <string>
#include <map>



void R2JJFitterHH_13TeV(double mass, std::string postfix="", bool dobands=false, int graviton = 0, double nEvents = 50000.)
{
    filePOSTfix=postfix;
    if (postfix.find("Graviton") != string::npos) iGraviton = true;
    nEventsInSignalMC = nEvents;
    signalScaler=analysisLumi/nEventsInSignalMC;

    for (int i = 0; i < np; i++){
      if ( fabs(mass - masses[i]) < 1e-5 ) {
	inTheList = true;
	break;
      }
    }



    runfits(mass, 0, dobands);

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

  TString fileBaseName(TString::Format("%s_%.0f_13TeV",signalname.Data(), mass));

  vector<string> cat_names;
  cat_names.push_back("4btag_cat0");
  cat_names.push_back("3btag_cat1");
  cat_names.push_back("2btag_cat2");


  TString fileBkgName("bkg_HH_13TeV");
  TString card_name("Xvv_reduced_models_Bkg_HH_13TeV.rs");
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();

  w->var("mgg")->setMin(MMIN);
  w->var("mgg")->setMax(MMAX);

// Add data to the workspace

  cout << "CREATE SIGNAL" << endl;

  if (inTheList) AddSigData(w, mass, signalsample,cat_names);

  SigModelFit(w, mass, signalname, cat_names);
  MakeSigWS(w, fileBaseName, signalname,cat_names);

  // ===============================

  AddBkgData(w,cat_names);
  BkgModelFit(w, dobands,cat_names, signalname, mass);
  MakeBkgWS(w, fileBkgName,cat_names);


  // =======================
  TFile* tf1=new TFile(Form("makeDatacards/%s_%s.root",filePOSTfix.data(),fileBaseName.Data()),"recreate");
  w->Write();
  tf1->Close();



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
    
    jjSig[c]     = (RooAbsPdf*)  w->pdf(TString::Format("%s_sig_%s",signalname.Data(), cat_names.at(c).c_str()));

    ((RooRealVar*) w->var(signalname+TString::Format("_sig_m0_%s",cat_names.at(c).c_str())))->setVal(MASS);
  
    cout << " Mass = " << MASS << endl;
      
    if (inTheList)
      jjSig[c]     ->fitTo(*sigToFit[c],Range(mass*0.8,mass*1.3),SumW2Error(kTRUE),RooFit::PrintEvalErrors(-1));
    else {
      double val_alpha, val_n, val_sigma, val_m0;
      getVal(mass, c, val_alpha, val_n, val_sigma, val_m0);
      ((RooRealVar*) w->var(signalname+TString::Format("_sig_m0_%s",cat_names.at(c).c_str())))->setVal(val_m0);
      ((RooRealVar*) w->var(signalname+TString::Format("_sig_sigma_%s",cat_names.at(c).c_str())))->setVal(val_sigma);
      ((RooRealVar*) w->var(signalname+TString::Format("_sig_alpha_%s",cat_names.at(c).c_str())))->setVal(val_alpha);
      ((RooRealVar*) w->var(signalname+TString::Format("_sig_n_%s",cat_names.at(c).c_str())))->setVal(val_n);
    }

    

    cout << " fitted " << endl;

// IMPORTANT: fix all pdf parameters to constant
    w->defineSet(TString::Format("SigPdfParam_%s",cat_names.at(c).c_str()), 
		 RooArgSet(*w->var(TString::Format("%s_sig_m0_%s",signalname.Data(), cat_names.at(c).c_str())),
			   *w->var(TString::Format("%s_sig_sigma_%s",signalname.Data(), cat_names.at(c).c_str())),
			   *w->var(TString::Format("%s_sig_alpha_%s",signalname.Data(), cat_names.at(c).c_str())),
			   *w->var(TString::Format("%s_sig_n_%s",signalname.Data(), cat_names.at(c).c_str())) ));

    cout << " defined " << endl;

    SetConstantParams(w->set(TString::Format("SigPdfParam_%s",cat_names.at(c).c_str())));
  }

  cout << " ========== Parameters constant ============ " << endl;

}



void BkgModelFit(RooWorkspace* w, Bool_t dobands, std::vector<string> cat_names, TString signalname, double mass) {

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

    bkg_fitTmp->fitTo(*data[c], Strategy(1),Minos(kFALSE), Range(minMassFit,maxMassFit),SumW2Error(kTRUE), Save(kTRUE),RooFit::PrintEvalErrors(-1));
    Double_t norm = data[c]->sumEntries();
    cout << "====================== norm = " << norm << endl;


    RooAbsReal* bkg_fitTmp2  = new RooRealVar(TString::Format("bkg_fit_%s_norm",cat_names.at(c).c_str()),"",norm,0.0,10000000);
    w->import(*bkg_fitTmp);
    w->import(*bkg_fitTmp2);

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
    jjSigPdf[c] = (RooAbsPdf*)  w->pdf(signalname+TString::Format("_sig_%s",cat_names.at(c).c_str()));
    wAll->import(*jjSigPdf[c]);
	cout<<c<<endl;
  }

// (2) Systematics on energy scale and resolution

  wAll->factory("CMS_sig_p1_jes_13TeV[0.0,-5.0,5.0]");
  wAll->factory("CMS_jj_sig_p1_jes_13TeV[0.012,0.012,0.012]");
  wAll->factory("sum::CMS_sig_p1_jes_sum(1.0,prod::CMS_sig_p1_jes_prod(CMS_sig_p1_jes_13TeV, CMS_jj_sig_p1_jes_13TeV))");
    for (int c = 0; c < ncat; ++c) {
      wAll->factory("prod::CMS_"+signalname+"_sig_m0_"+TString::Format("%s",cat_names.at(c).c_str())+"("+signalname+"_sig_m0_"+TString::Format("%s",cat_names.at(c).c_str())+", CMS_sig_p1_jes_sum)");
    }

// (3) Systematics on resolution: create new sigmas

    // apply JER resolution smearing from 
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
    // Assume ~7% smearing +- 3% for each jet. Divide by sqrt(2) for both

  wAll->factory("CMS_sig_p2_jer_13TeV[0.0,-5.0,5.0]");
  wAll->factory("CMS_jj_sig_p2_jer_13TeV[0.02,0.02,0.02]");
  wAll->factory("sum::CMS_sig_p2_jer_sum(1.00,prod::CMS_sig_p2_jer_prod(CMS_sig_p2_jer_13TeV, CMS_jj_sig_p2_jer_13TeV))");

    for (int c = 0; c < ncat; ++c) {
      wAll->factory("prod::CMS_"+signalname+"_sig_sigma_"+TString::Format("%s",cat_names.at(c).c_str())+"("+signalname+"_sig_sigma_"+TString::Format("%s",cat_names.at(c).c_str())+", CMS_sig_p2_jer_sum)");
    }


// (4) do reparametrization of signal
  for (int c = 0; c < ncat; ++c) {
    wAll->factory(
		  TString::Format("EDIT::%s_sig_%s",signalname.Data(),cat_names.at(c).c_str()) +
		  TString::Format("(%s_sig_%s, ",signalname.Data(),cat_names.at(c).c_str()) + 
		  TString::Format("%s_sig_m0_%s=CMS_%s_sig_m0_%s, ",signalname.Data(),cat_names.at(c).c_str(),signalname.Data(),cat_names.at(c).c_str()) +
		  TString::Format("%s_sig_sigma_%s=CMS_%s_sig_sigma_%s)",signalname.Data(),cat_names.at(c).c_str(),signalname.Data(),cat_names.at(c).c_str()));
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
    wAll->import(*data[c], Rename(TString::Format("data_obs_unbinned_%s",cat_names.at(c).c_str())));

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
    double errMean = (wAll->var(TString::Format("bkg_fit_slope1_%s",cat_names.at(c).c_str())))->getError();
   
    cout << " ============================ cat " << c << " slope 1 = " << mean << " err1 = " << errMean << endl;

    (wAll->var(TString::Format("bkg_fit_slope1_%s",cat_names.at(c).c_str())))->setConstant(true);

    wAll->factory(TString::Format("CMS_bkg_fit_slope1_eps_%s[0.0,-5.0,5.0]", cat_names.at(c).c_str()));
    wAll->factory(TString::Format("CMS_bkg_fit_slope1_err_%s[%g, %g, %g]",  cat_names.at(c).c_str(), errMean, errMean, errMean));
    wAll->factory(TString::Format("sum::CMS_bkg_fit_slope1_sum_%s(1.0,prod::CMS_sig_slope1_jes_prod_%s(CMS_bkg_fit_slope1_eps_%s, CMS_bkg_fit_slope1_err_%s))", cat_names.at(c).c_str(), cat_names.at(c).c_str(), cat_names.at(c).c_str(), cat_names.at(c).c_str()));


    wAll->factory(TString::Format("prod::CMS_bkg_fit_slope1_%s(bkg_fit_slope1_%s,CMS_bkg_fit_slope1_sum_%s)", cat_names.at(c).c_str(), cat_names.at(c).c_str(), cat_names.at(c).c_str()));


   if (c > 0){
     mean = (wAll->var(TString::Format("bkg_fit_slope2_%s",cat_names.at(c).c_str())))->getVal();
     min = (wAll->var(TString::Format("bkg_fit_slope2_%s",cat_names.at(c).c_str())))->getMin();
     max = (wAll->var(TString::Format("bkg_fit_slope2_%s",cat_names.at(c).c_str())))->getMax();
     errMean = (wAll->var(TString::Format("bkg_fit_slope2_%s",cat_names.at(c).c_str())))->getError();

     cout << " ============================ cat " << c << " slope 2 = " << mean << " err2 = " << errMean << endl;

     (wAll->var(TString::Format("bkg_fit_slope2_%s",cat_names.at(c).c_str())))->setConstant(true);

     wAll->factory(TString::Format("CMS_bkg_fit_slope2_eps_%s[0.0,-5.0,5.0]", cat_names.at(c).c_str()));
     wAll->factory(TString::Format("CMS_bkg_fit_slope2_err_%s[%g, %g, %g]",  cat_names.at(c).c_str(), errMean, errMean, errMean));
     wAll->factory(TString::Format("sum::CMS_bkg_fit_slope2_sum_%s(1.0,prod::CMS_sig_slope2_jes_prod_%s(CMS_bkg_fit_slope2_eps_%s, CMS_bkg_fit_slope2_err_%s))", 
				 cat_names.at(c).c_str(), cat_names.at(c).c_str(), cat_names.at(c).c_str(), cat_names.at(c).c_str()));


     wAll->factory(TString::Format("prod::CMS_bkg_fit_slope2_%s(bkg_fit_slope2_%s,CMS_bkg_fit_slope2_sum_%s)", cat_names.at(c).c_str(), cat_names.at(c).c_str(), cat_names.at(c).c_str()));

   }

    cout << "Done For category " << c << endl;    
  }
  
  
  cout << "Imported" << endl;

// (2) do reparametrization of background

  for (int c = 0; c < ncat; ++c) {

        
   
    TString sFormat = 	    TString::Format("EDIT::CMS_bkg_fit_%s(bkg_fit_%s,",cat_names.at(c).c_str(),cat_names.at(c).c_str()) +
      TString::Format(" bkg_fit_%s_norm=CMS_bkg_fit_%s_norm,", cat_names.at(c).c_str(),cat_names.at(c).c_str())+
      TString::Format(" bkg_fit_slope1_%s=CMS_bkg_fit_slope1_%s)", cat_names.at(c).c_str(),cat_names.at(c).c_str());
  
    
    if (c > 0) {
      sFormat = 	    TString::Format("EDIT::CMS_bkg_fit_%s(bkg_fit_%s,",cat_names.at(c).c_str(),cat_names.at(c).c_str()) +
	TString::Format(" bkg_fit_%s_norm=CMS_bkg_fit_%s_norm,", cat_names.at(c).c_str(),cat_names.at(c).c_str())+
	TString::Format(" bkg_fit_slope1_%s=CMS_bkg_fit_slope1_%s,", cat_names.at(c).c_str(),cat_names.at(c).c_str())+
	TString::Format(" bkg_fit_slope2_%s=CMS_bkg_fit_slope2_%s)", cat_names.at(c).c_str(),cat_names.at(c).c_str());
    }
    

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

  /*   
  for (int c = 0; c < ncat; ++c) {
    printf("CMS_bkg_fit_slope1_%s  param  %.4f  %.3f   # Mean and absolute uncertainty on background slope\n",
	   cat_names.at(c).c_str(), (wAll->var(TString::Format("CMS_bkg_fit_slope1_%s",cat_names.at(c).c_str())))->getVal(), 10.);
    if (c > 0) printf("CMS_bkg_fit_slope2_%s  param  %.4f  %.3f   # Mean and absolute uncertainty on background slope\n",
	   cat_names.at(c).c_str(), (wAll->var(TString::Format("CMS_bkg_fit_slope2_%s",cat_names.at(c).c_str())))->getVal(), 10.);
  }
  */
  cout << "BKG WS DONE" << endl;

  return;
}

