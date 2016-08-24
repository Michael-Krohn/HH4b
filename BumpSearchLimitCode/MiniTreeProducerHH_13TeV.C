{


  double mgg, mjj,evWeight, mtot, normWeight;
 Int_t categories;

 evWeight = 1.0;
 normWeight = 1;
 bool bSubstr = true;
 bool bData = false;

 // string sInFile = string("input/MassPlotFineBins_Moriond_Silver_tau21cut.root");
 // if (bSubstr) sInFile = string("input/MassPlotFineBins_subtr_Moriond_Silver_tau21cut.root");

 for (int itech= 1; itech < 2; itech++){

   bSubstr = itech;


   string sInFile = string("input/MassPlotFineBins_Moriond_Silver.root");
   if (bSubstr) sInFile = string("input/MassPlotFineBins_subtr_Moriond_Silver.root");
   cout << sInFile.c_str() << endl;
   TFile file0(sInFile.c_str(), "read");
   

   string sOutFile("MiniTrees/Background_HH_13TeV/dijetHH_miniTree.root");
   if (bSubstr && bData) sOutFile = string("MiniTrees/Background_HH_13TeV/dijetHH_data_subtr_miniTree.root");
   else if (bSubstr && !bData) sOutFile = string("MiniTrees/Background_HH_13TeV/dijetHH_qcd_subtr_miniTree.root");
   else if (!bSubstr && bData) sOutFile = string("MiniTrees/Background_HH_13TeV/dijetHH_data_miniTree.root");
   else if (!bSubstr && !bData) sOutFile = string("MiniTrees/Background_HH_13TeV/dijetHH_qcd_miniTree.root");
   
   TFile f1(sOutFile.c_str(), "recreate");
   f1.cd();
   
   TTree *TCVARS = new TTree("TCVARS", "HH selection");
   TCVARS->Branch("mgg",&mgg,"mgg/D");
   
   TCVARS->Branch("evWeight",&evWeight,"evWeight/D");
   TCVARS->Branch("normWeight",&normWeight,"normWeight/D");
   
   TCVARS->Branch("categories",&categories,"categories/I");
   
   
   TRandom3 R;
   
   
   
   for (Int_t iCat = 0; iCat < 3; iCat++){
     double sum = 0;
     
     TH1D* hMass = (TH1D*) file0.Get(Form("QCD_cat%d;1",iCat));
     if (bData)  hMass = (TH1D*) file0.Get(Form("Data_cat%d;1",iCat));
     
     TAxis* Axis =   hMass->GetXaxis();
     for (int i = 1 ; i < hMass->GetNbinsX()+1; i++){
       double N = hMass->GetBinContent(i);
       
     mgg = Axis->GetBinCenter(i);
     

     int intPart = TMath::Nint(N);
     double resid = intPart - N;
     double rnd = R.Uniform(1.);
     if (resid > 0) normWeight = rnd > resid ? intPart : intPart-1; 
     else normWeight = rnd > fabs(resid) ?  intPart+0. : intPart+1.; 
     
     if (bData) normWeight = N;
     //     normWeight = intPart*1.0;
     sum += normWeight;

     categories = iCat;
 
     if (N > 1e-10 && N <1) {
       //     cout << "i = " << i << " N = " << N << " normWeight = " << normWeight << " binCenter = " << hMass->GetBinCenter(i)<< " categories = " << categories << endl;
       // cout << "N = " << N << " intPart = " << intPart << " resid = " << resid << " rnd = " << rnd << " norm = " << normWeight << endl; 
     }

    if (N > 1e-10 && mgg > 1000 && mgg < 5000) TCVARS->Fill();
   }

     cout << "sum = " << sum << " integral = " << hMass->Integral() << endl;
     
   }
   
   cout << sInFile.c_str() << endl;
   cout << sOutFile.c_str() << endl;
   
   TCVARS->Write();
   f1.Close();
   file0.Close();
 } 

}


