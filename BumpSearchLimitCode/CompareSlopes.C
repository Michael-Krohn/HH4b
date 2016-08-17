
using namespace RooFit;
using namespace RooStats ;

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

void CompareSlopes(){

   // Load the combine Library 
   gSystem->Load("libHiggsAnalysisCombinedLimit.so");
   style();




   // Open the dummy H->gg workspace 
   TFile *f_hgg = TFile::Open("workspaces/Graviton_subtrCMS_jj_bkg_HH_13TeV.root");
   RooWorkspace *w_hgg = (RooWorkspace*)f_hgg->Get("w_all");
   // w_all:CMS_bkg_fit_CMS_jj_4btag_cat0
   // The observable (CMS_hgg_mass in the workspace)
   RooRealVar *mass =  w_hgg->var("mgg");

   string file;

   // Get three of the functions inside, exponential, linear polynomial, power law
   file = string("CMS_bkg_fit_CMS_jj_4btag_cat0");
   RooAbsPdf *pdf_exp_cat0 = w_hgg->pdf(file.c_str());
   file = string("CMS_bkg_fit_CMS_jj_3btag_cat1");
   RooAbsPdf *pdf_exp_cat1 = w_hgg->pdf(file.c_str());
   file = string("CMS_bkg_fit_CMS_jj_2btag_cat2");
   RooAbsPdf *pdf_exp_cat2 = w_hgg->pdf(file.c_str());

   //   (w_hgg->var("CMS_bkg_fit_slope1_eps_CMS_jj_3btag_HPLP_cat2"))->setConstant(true);
   //   (w_hgg->var("CMS_bkg_fit_slope2_eps_CMS_jj_3btag_HPLP_cat2"))->setConstant(true);

   RooDataSet *data_cat0 = w_hgg->data("data_obs_unbinned_CMS_jj_4btag_cat0");
   RooDataSet *data_cat1 = w_hgg->data("data_obs_unbinned_CMS_jj_3btag_cat1");
   RooDataSet *data_cat2 = w_hgg->data("data_obs_unbinned_CMS_jj_2btag_cat2");  // THIS IS NOW 2 btag category

   /*
   pdf_exp_cat0->Print();
   pdf_exp_cat1->Print();
   pdf_exp_cat2->Print();
   */

   TFile *f_sigmas_cat0 = TFile::Open("plots/backgroundsGraviton_subtr_channel0_withband_parameters.root");

   TGraphAsymmErrors*	onesigma_cat0 = f_sigmas_cat0->Get("onesigma;1");
   TGraphAsymmErrors*	twosigma_cat0 = f_sigmas_cat0->Get("twosigma;1");

   TFile *f_sigmas_cat1 = TFile::Open("plots/backgroundsGraviton_subtr_channel1_withband_parameters.root");

   TGraphAsymmErrors*	onesigma_cat1 = f_sigmas_cat1->Get("onesigma;1");
   TGraphAsymmErrors*	twosigma_cat1 = f_sigmas_cat1->Get("twosigma;1");



   TCanvas* ctmp = new TCanvas("ctmp","Shapes",0,0,500,500);
   ctmp->SetLogy();

   // Make a plot (data is a toy dataset)
   RooPlot *plot = mass->frame(Range(1000, 2600), Bins(32)); //  data->plotOn(plot);
   data_cat0->plotOn(plot);
   plot->Draw();

   twosigma_cat0->SetLineColor(kYellow+1);
   twosigma_cat0->SetFillColor(kYellow+1);
   twosigma_cat0->SetMarkerColor(kYellow+1);

   onesigma_cat0->SetLineColor(kGreen+1);
   onesigma_cat0->SetFillColor(kGreen+1);
   onesigma_cat0->SetMarkerColor(kGreen+1);


   twosigma_cat0->Draw("L3 SAME");
   onesigma_cat0->Draw("L3 SAME");

 
   //   data_cat0->plotOn(plot);

   pdf_exp_cat2->plotOn(plot,RooFit::LineColor(kBlue));
   pdf_exp_cat0->plotOn(plot,RooFit::LineColor(kRed));
   plot->SetTitle("4b cat: PDFs extracted from different categories");
   plot->Draw("SAME");



   plot->SetMinimum(1e-1);
 
    TLegend *legmc = new TLegend(0.52,0.65,0.87,0.85);
    legmc->AddEntry(plot->getObject(0),"4 btag cat. ","EP");
    legmc->AddEntry(plot->getObject(1),"Shape from 2 btag cat.","L");
    legmc->AddEntry(plot->getObject(2),"Shape from 4 btag cat.","L");

    legmc->SetFillStyle(0);

    legmc->Draw();

    ctmp->SetLogy();
    out = string("plots/backgrounds_bump_hunt_and_alphabet_4btag_cat0.pdf");
    ctmp->SaveAs(out.c_str());


    ctmp->Clear();
    RooPlot *plot = mass->frame(Range(1000, 2600), Bins(32)); //  data->plotOn(plot);
    data_cat1->plotOn(plot);
   

    twosigma_cat1->SetLineColor(kYellow+1);
    twosigma_cat1->SetFillColor(kYellow+1);
    twosigma_cat1->SetMarkerColor(kYellow+1);

    onesigma_cat1->SetLineColor(kGreen+1);
    onesigma_cat1->SetFillColor(kGreen+1);
    onesigma_cat1->SetMarkerColor(kGreen+1);

    plot->Draw(); 

    twosigma_cat1->Draw("L3 SAME");
    onesigma_cat1->Draw("L3 SAME");
    
    //   data_cat1->plotOn(plot);

    pdf_exp_cat2->plotOn(plot,RooFit::LineColor(kBlue));
    pdf_exp_cat1->plotOn(plot,RooFit::LineColor(kRed));
    plot->SetTitle("3b cat: PDFs extracted from different categories");
    
    plot->SetMinimum(1e-1);
    plot->Draw("SAME");  


    TLegend *legmc = new TLegend(0.52,0.65,0.87,0.85);
    legmc->AddEntry(plot->getObject(0),"3 btag cat. ","EP");
    legmc->AddEntry(plot->getObject(1),"Shape from 2 btag cat.","L");
    legmc->AddEntry(plot->getObject(2),"Shape from 3 btag cat.","L");

    legmc->SetFillStyle(0);

    legmc->Draw();

    ctmp->SetLogy();
    out = string("plots/backgrounds_bump_hunt_and_alphabet_3btag_cat1.pdf");
    ctmp->SaveAs(out.c_str());

    //    twosigma_cat1->Draw("L3");
    //  onesigma_cat1->Draw("L3 SAME");

}
