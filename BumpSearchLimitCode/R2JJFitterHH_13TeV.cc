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
    style();

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

  TString fileBaseName(TString::Format("%s_%.0f_13TeV", signalname.Data(), mass));

  vector<string> cat_names;
  cat_names.push_back("4btag_cat0");
  cat_names.push_back("3btag_cat1");
  cat_names.push_back("2btag_cat2");


  TString fileBkgName("bkg_HH_13TeV");
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
  
// Make plots for signal shape
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


    


    jjSig[c]       = (RooAbsPdf*)  w->pdf(signalname+TString::Format("_sig_%s",cat_names.at(c).c_str()));

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
	
	//	cout << "total integral = "  <<  data[c]->sumEntries();


	double xmin = xc-fabs(xlowErr), xmax = xc+fabs(xhighErr);
	
	mgg->setRange("A",xmin,xmax);

	RooAbsReal* intBin0 = bkg_fitTmp_2par->createIntegral(*set,*set,"A") ;

	RooAbsReal* intBin_3par = bkg_fitTmp->createIntegral(*set,*set,"A") ;
    
	double dintBin0 = intBin0->getVal();
	double dintBin_3par = intBin_3par->getVal();
	//	cout << "=================== Bin = " << i << " xmin = " << xmin << " xmax = " << xmax << " xlowErr = " << xlowErr << " xhighErr = " << xhighErr << " N-L = " << N-L << " U-N = " << U-N << " bin content = " << y0 << " intBin = " << dintBin0 << " unnormalised integral = " << dintBin0*normalisation << " yc = " << yc << endl;

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
    
    //cout <<  " p1_10 = " << p1_10 << " p2_10 = " << p2_10 << " Ftest_23 = " << Ftest_10 << endl;

    cout << "chi2 = " << chi2 << endl;
    cout << "chi2_3par = " << chi2_3par << endl;
    //    cout << "F Test CL = " << good_CL23 << endl;


 


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


void MakePlots(RooWorkspace* w, Float_t mass, RooFitResult** fitresults, TString signalname, std::vector<string> cat_names) {

  cout << "Start plotting" << endl; 

  Int_t ncat = NCAT;

// retrieve data sets from the workspace

  RooDataSet* data[9];  
  RooDataSet* signal[9];
  RooAbsPdf*  jjSig[9];

  for (int c = 0; c < ncat; ++c) {
    if (inTheList) signal[c]       = (RooDataSet*) w->data(TString::Format("Sig_%s",cat_names.at(c).c_str()));
    jjSig[c]       = (RooAbsPdf*)  w->pdf(signalname+TString::Format("_sig_%s",cat_names.at(c).c_str()));
  }

// retrieve mass observable from the workspace
  RooRealVar* mgg     = w->var("mgg");  
  mgg->setUnit("GeV");

// retrieve pdfs after the fits
// Signal Model

  RooAbsPdf* jjSigAll       = w->pdf(signalname);

 
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

  RooPlot* plotjj[9];
  for (int c = 0; c < ncat; ++c) {
    plotjj[c] = mgg->frame(Range(minMassFit,maxMassFit),Bins(nBinsMass));
    if (inTheList) signal[c]->plotOn(plotjj[c],LineColor(kWhite),MarkerColor(kWhite));    

    jjSig[c]  ->plotOn(plotjj[c]);

    if (inTheList) signal[c]  ->plotOn(plotjj[c]);

        
    TCanvas* dummy = new TCanvas(Form("dummy_%d",c), "dummy",0, 0, 400, 400);
    TH1F *hist = new TH1F(Form("hist_%d",c), "hist", 400, minMassFit, maxMassFit);
 
    plotjj[c]->SetTitle("");      
    plotjj[c]->SetMinimum(0.0);
    plotjj[c]->SetMaximum(1.40*plotjj[c]->GetMaximum());
    plotjj[c]->GetXaxis()->SetTitle("m_{jj} [GeV]");

    TCanvas* ctmp_sig = new TCanvas(Form("ctmp_sig_%d",c),"jj Background Categories",0,0,500,500);
    plotjj[c]->Draw();  
    plotjj[c]->Print();

    plotjj[c]->Draw("SAME");  
    TLegend *legmc = new TLegend(0.62,0.75,0.92,0.9);
    legmc->AddEntry(plotjj[c]->getObject(2),"Simulation","LPE");
    legmc->AddEntry(plotjj[c]->getObject(1),"Parametric Model","L");
    
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();    
    

    TLatex *lat  = new TLatex(minMassFit+1.5,0.85*plotjj[c]->GetMaximum(),"#scale[1.0]{CMS Preliminary}");
    lat->Draw();

    int iMass = abs(mass);

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
    jjSigPdf[c] = (RooAbsPdf*)  w->pdf(signalname+TString::Format("_sig_%s",cat_names.at(c).c_str()));
    wAll->import(*jjSigPdf[c]);
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
      wAll->factory("prod::CMS_"+signalname+"_sig_sigma_"+TString::Format("%s",cat_names.at(c).c_str())+"(jj_"+signalname+"_sig_sigma_"+TString::Format("%s",cat_names.at(c).c_str())+", CMS_sig_p2_jer_sum)");
    }


// (4) do reparametrization of signal
  for (int c = 0; c < ncat; ++c) {
    wAll->factory(
		  TString::Format("EDIT::%s_sig_%s",signalname.Data(),cat_names.at(c).c_str()) +
		  TString::Format("(%s_%s, ",signalname.Data(),cat_names.at(c).c_str()) + 
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


  outFile << "# Fully Hadronic HH analysis" << endl;
  outFile << "imax 1" << endl;
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;

  outFile << Form("shapes data_obs bin_%s ", cat_names[iChan].c_str()) << wsDir+TString(fileBkgName)+".root" << Form(" w_all:data_obs_%s", cat_names[iChan].c_str()) << endl;
  outFile << Form("shapes bkg_fit bin_%s ", cat_names[iChan].c_str()) <<  wsDir+TString(fileBkgName)+".root" << Form(" w_all:CMS_bkg_fit_%s", cat_names[2].c_str()) << endl;
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


  outFile << "# Parametric shape uncertainties, entered by hand." << endl;
  outFile << Form("CMS_sig_p1_jes_13TeV    param   0.0   1.0   # dijet mass shift due to JES uncertainty") << endl;
  outFile << Form("CMS_sig_p2_jer_13TeV     param   0.0   1.0   # dijet mass resolution shift due to JER uncertainty") << endl;
 
  outFile << Form("CMS_bkg_fit_%s_norm           flatParam  # Normalization uncertainty on background slope",cat_names[2].c_str()) << endl;

  //  outFile << Form("CMS_bkg_fit_slope1_%s         flatParam  # Mean and absolute uncertainty on background slope",cat_names[iChan].c_str()) << endl;
  //   wAll->factory(TString::Format("CMS_bkg_fit_slope1_eps_%s[0.0,-5.0,5.0]", cat_names.at(c).c_str()));
  outFile << Form("CMS_bkg_fit_slope1_eps_%s     param  0.0  1.0  # Mean and absolute uncertainty on background slope",cat_names[2].c_str()) << endl;
  outFile << Form("CMS_bkg_fit_slope2_eps_%s     param  0.0  1.0  # Mean and absolute uncertainty on background slope",cat_names[2].c_str()) << endl;

    //outFile << Form("CMS_bkg_fit_slope2_%s         flatParam  # Mean and absolute uncertainty on background levelled parameter",cat_names[iChan].c_str()) << endl;

  outFile.close();

  cout << "Write data card in: " << filename << " file" << endl;

  return;
}


