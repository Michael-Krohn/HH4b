
void BiasStudies(int mass, int ievt)
{
  TFile *_file0 = TFile::Open(Form("toys/mlfit_%d_mu%d.root", mass, ievt));
 
  TCanvas* c1 = new TCanvas("c1", "");
  TTree* tree_fit_sb = (TTree*) _file0->Get("tree_fit_sb;1");

  //  tree_fit_sb->Draw(Form("(mu-%d)/(muLoErr+muHiErr)*2>>h(50,-4,4)", ievt));
  tree_fit_sb->Draw(Form("(mu-%d)/(muLoErr)>>hBUp(50,-4,4)", ievt), Form("mu-%d > 0 && muHiErr < 4 && muLoErr < 4", ievt));
  TH1D* hBUp = (TH1D*) _file0->Get("hBUp");
  tree_fit_sb->Draw(Form("(mu-%d)/(muHiErr)>>hBDo(50,-4,4)", ievt), Form("mu-%d < 0 && muLoErr < 4 && muLoErr < 4", ievt));
  TH1D* hBias = (TH1D*) _file0->Get("hBDo");
  hBias->SetTitle(Form("(mu-%d)/Err and mX = %d", ievt, mass));
  hBias->Add(hBUp);

  hBias->Draw();

  TF1* funG = new TF1("funG", "[0]*TMath::Gaus(x,[1],[2])");
  //hist->SetStats(0);

  gStyle->SetOptFit(0011);

  funG->SetParName(1, "bias");
  funG->SetParName(2, "sigma");

  funG->SetParameter(0, 40.);
  funG->SetParameter(1, 0.);
  funG->SetParameter(2, 1.);
  
  hBias->Fit("funG");


  c1->SaveAs(Form("toys/Bias_Graviton%d_mu%d.png",mass, ievt));
  c1->SaveAs(Form("toys/Bias_Graviton%d_mu%d.pdf",mass, ievt));

  tree_fit_sb->Draw(Form("mu>>hmu(20,-10.,10+%d*4.)", ievt));
  TH1D* hmu = (TH1D*) _file0->Get("hmu");
  hmu->Fit("gaus");

  hmu->SetStats(0);
  c1->SaveAs(Form("toys/Mu_Graviton%d_mu%d.png",mass, ievt));
  c1->SaveAs(Form("toys/Mu_Graviton%d_mu%d.pdf",mass, ievt));

  c1->Clear();

  tree_fit_sb->Draw("mu/muErr>>hmuSigma(20,0.,10)");
  TH1D* hmuSigma = (TH1D*) _file0->Get("hmuSigma");
  hmuSigma->SetTitle(Form("#{sigma} = %d/Err and mX = %d", ievt, mass));

  hmuSigma->Fit("gaus");

  hmuSigma->SetStats(0);
  c1->SaveAs(Form("toys/MuSigmas_Graviton%d_mu%d.png",mass, ievt));
  c1->SaveAs(Form("toys/MuSigmas_Graviton%d_mu%d.pdf",mass, ievt));

  c1->Clear();

  tree_fit_sb->Draw("muErr>>hErr(50,0,10)");
  TH1F *hErr = (TH1F*)c1->GetPrimitive("hErr");

  tree_fit_sb->Draw("muHiErr>>hErrUp(50,0,10)");
  TH1F *hErrUp = (TH1F*)c1->GetPrimitive("hErrUp");

  tree_fit_sb->Draw("muLoErr>>hErrDo(50,0,10)");
  TH1F *hErrDo = (TH1F*)c1->GetPrimitive("hErrDo");


  hErr->SetLineColor(kBlack);
  hErr->SetLineStyle(2);
  hErr->Draw();

  hErrUp->SetLineColor(kRed);
  hErrUp->Draw("SAME");

  hErrDo->SetLineColor(kBlue);
  hErrDo->Draw("SAME");


  c1->SaveAs(Form("toys/Err_Graviton%d_mu%d.png",mass, ievt));

}
