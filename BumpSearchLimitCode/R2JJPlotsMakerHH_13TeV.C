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
double analysisLumi = 2.69; // Luminosity you use in your analysis
double nEventsInSignalMC = 0.; //number of events in Signal MC sample
int iGraviton = 0;

bool inTheList = false;

const int np(6);
double masses[np] = { 1200.0,1400.0,1600.0,1800.0,2000.0,2500.0 } ;

double signalScaler=0;//analysisLumi/nEventsInSignalMC; // assume signal cross section on 1/fb

TString mainCut("1");

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

void MakePlots(RooWorkspace* w, Float_t mass,TString signalname, std::vector<string> cat_names) {

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


void MakeBkgPlots(RooWorkspace* w, Bool_t dobands, std::vector<string> cat_names, TString signalname, double mass){
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
	w->import(*bkg_fitTmp_power);
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

void R2JJPlotsMakerHH_13TeV(double mass, std::string postfix="", bool dobands=false, int graviton = 0, double nEvents = 50000.)
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

	TString signalname="HH";
	TString fileBaseName(TString::Format("%s_%.0f_13TeV", signalname.Data(), mass));

	vector<string> cat_names;
	cat_names.push_back("4btag_cat0");
	cat_names.push_back("3btag_cat1");
	cat_names.push_back("2btag_cat2");
	
	TString fileBkgName("bkg_HH_13TeV");
	TString card_name("Xvv_reduced_models_Bkg_HH_13TeV.rs");

  
	TFile* tf1=TFile::Open(Form("makeDatacards/%s_%s.root",filePOSTfix.data(),fileBaseName.Data()));
	RooWorkspace* w= tf1->Get("HLFactory_ws");
	tf1->Close();
	
	cout << "MAKE PLOTS" << endl;
	// Make plots for signal shape
	MakePlots(w, mass, signalname, cat_names);
	MakeBkgPlots(w, dobands,cat_names,  signalname, mass);
	//MakeBkgPlots(w, mass, signalname, cat_names);
	
    //runfits(mass, 0, dobands);

}




