double Search(string cfilename){

  ifstream fileInput;
  int offset;
  string line;
  char* search = "rate"; // test variable to search in file
  
  double output = 0;
 
  // open file to search
  fileInput.open(cfilename.c_str());
  if(fileInput.is_open()) {
    while(!fileInput.eof()) {
      getline(fileInput, line);
      if ((offset = line.find(search, 0)) != string::npos) {
	cout << "found: " << line << endl;
	char* ssizes = strtok(line.c_str()," ,\t");
	int iterator = -1;
	  while (ssizes != NULL)
	    {
	      iterator++;
	      if (iterator == 1) {
		output = atof(ssizes);
		cout << "atof = " << output << endl;
	      }
	      ssizes = strtok (NULL, " ");
	    }
      }
    }
    fileInput.close();
  }
  else cout << "Unable to open file.";

  return output;

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


void ShapesPlotCompar(int category){

  using namespace RooFit;
  using namespace RooStats ;

 string Particle = "Graviton";
// string Particle = "Radion";

 style();

 TCanvas* ctmp = new TCanvas("ctmp","jj Background Categories",0,0,600,600);

  cout << "category " << category << endl; 

   // Load the combine Library 
  //   gSystem->Load("libHiggsAnalysisCombinedLimit.so");

   //   string file;
   if (category == 0) file = string("4btag_cat0");
   if (category == 1) file = string("3btag_HPHP_cat1");

   string filein = "HH_jj_sig_CMS_jj_" + file; 

   int amass[2] = {1600, 1600};
   int colors[2] = {kRed, kBlue}
   RooPlot *plot;
   RooAbsPdf *pdf_exp;
   RooRealVar *mass; 
   TFile *f_hgg ;
   string datacard;

   TLegend *legmc = new TLegend(0.70,0.55,0.99,0.70);
   if (Particle.find("Graviton") != string::npos ) legmc->SetHeader("Bulk Graviton");
   if (Particle.find("Radion") != string::npos ) legmc->SetHeader("Radion #sigma = 10 fb");

   for (int i = 0; i <2; i++){

     cout << "i = " << i << endl;

     if (i == 0){ 

       f_hgg = TFile::Open(Form("workspaces/%s_subtrCMS_jj_HH_%d_13TeV.root", Particle.c_str(), amass[i]));  
       datacard = Form("datacards/%s_subtrCMS_jj_HH_%d_13TeV_CMS_jj_", Particle.c_str(), amass[i]) + file + ".txt";

     } else {
       
       f_hgg = TFile::Open(Form("workspaces/%sCMS_jj_HH_%d_13TeV.root", Particle.c_str(), amass[i]));  
       datacard = Form("datacards/%sCMS_jj_HH_%d_13TeV_CMS_jj_", Particle.c_str(), amass[i]) + file + ".txt";

     }

     cout << "i = " << i << endl;
     cout << datacard.c_str() << endl;

     RooWorkspace *w_hgg = (RooWorkspace*)f_hgg->Get("w_all");

     
     double Norm = Search(datacard);
     cout << "Normalisation = " << Norm << endl;



       // Make a plot (data is a toy dataset)
     if (i == 0) {
       mass = w_hgg->var("mgg");
       plot = mass->frame(); //  data->plotOn(plot);
  
     }


     cout << "i = " << i << endl;
     pdf_exp = w_hgg->pdf(filein.c_str());
  
     pdf_exp->plotOn(plot,RooFit::LineColor(colors[i]), RooFit::LineStyle(i+1), Range("fitrange"),NormRange("fitrange"), Normalization(10*Norm,RooAbsReal::NumEvent));
     //   pdf_pol->plotOn(plot,RooFit::LineColor(kBlue));
     //  pdf_pow->plotOn(plot,RooFit::LineColor(kRed));
     if (i == 0) legmc->AddEntry(plot->getObject(i),Form("M_{jj}^{red}", amass[i]),"L");
     else legmc->AddEntry(plot->getObject(i),Form("M_{jj}", amass[i]),"L");
     cout << "i = " << i << endl;
   }
   plot->SetTitle("; M_{jj} [GeV]; Events / GeV");
   plot->Draw();


    TLatex *latexLabel = new TLatex();
    latexLabel->SetTextSize(0.75 * ctmp->GetTopMargin());
    latexLabel->SetNDC();
    latexLabel->SetTextFont(42); // helvetica
    latexLabel->DrawLatex(0.73, 0.96, "2.7 fb^{-1} (13 TeV)");
    latexLabel->SetTextFont(61); // helvetica bold face
    latexLabel->SetTextAlign(32);
    latexLabel->DrawLatex(0.90, 0.88, "CMS");


    //    TLatex *latexLabelPrel = new TLatex();
    latexLabel->SetTextFont(52); // helvetica bold face
    latexLabel->SetTextAlign(32);
    latexLabel->DrawLatex(0.90, 0.83, "Simulation");

    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->SetTextFont(62);


    legmc->Draw();

   latexLabel->SetTextColor(kBlue);
    latexLabel->SetTextAlign(22);    
    latexLabel->DrawLatex(0.70, 0.35, "pp #rightarrow X #rightarrow HH  #rightarrow b#bar{b}b#bar{b}");


    if (category == 0)    latexLabel->DrawLatex(0.70, 0.30, " 4 b tag category");
    else if (category == 1)    latexLabel->DrawLatex(0.70, 0.30, " 3 b tag category");
        
    latexLabel->SetTextColor(kBlack);
    latexLabel->DrawLatex(0.70, 0.25, "#sigma (M_{X} = 1.6 TeV) = 10 fb");

    ctmp->SaveAs(Form("plots/ShapesComparison_cat%d_%s_compar.png", category, Particle.c_str()));
    ctmp->SaveAs(Form("plots/ShapesComparison_cat%d_%s_compar.C", category, Particle.c_str()));
    ctmp->SaveAs(Form("plots/ShapesComparison_cat%d_%s_compar.pdf", category, Particle.c_str()));

   // Get three of the functions inside, exponential, linear polynomial, power law
 

   // Fit the functions to the data to set the "prefit" state (note this can and should be redone with combine when doing 
   // bias studies as one typically throws toys from the "best-fit"
   //RooDataSet *data = w_hgg->data("roohist_data_mass_cat1_toy1_cutrange__CMS_hgg_mass");
   //pdf_exp->fitTo(*data);  // index 0
   //   pdf_pow->fitTo(*data); // index 1 
   //   pdf_pol->fitTo(*data);   // index 2



}
