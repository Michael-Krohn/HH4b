{


  double mgg, mjj,evWeight, mtot, normWeight, massrange;
  int categories;

  evWeight = 1.0;
  normWeight = 1;

  //bool bSubtr = true;//false;
  bool bSubtr = true;

  for (int itech= 1; itech < 2; itech++){

    bSubstr = itech;

    cout << "bla 0" << endl;

    string inFile("input/MassPlotFineBins_Moriond_Silver.root");
    if (bSubtr) inFile = string("input/MassPlotFineBins_subtr_Moriond_Silver.root");
    cout<<itech<<","<<bSubtr<<","<<inFile<<endl;

    TFile* file0 = new TFile(inFile.c_str(), "read");
 
    for (int sighyp = 1; sighyp <2; sighyp++){
   
      cout << "bla 1" << endl;

      massrange=21;
 
      for (int iMass = 1; iMass<massrange; iMass++){
   
	cout << "iMass = " << iMass << endl;
     
	if ( iMass == 1 || iMass == 3  || iMass == 5  || iMass == 7 || iMass == 9) continue;
	if (iMass > 10 && iMass < 15) continue;
	if (iMass > 15 && iMass < 20) continue;
  

	int Mass = 1000+iMass*100;   
	string outFile("dijetHH_RadionHH");
	if (bSubtr) outFile = string("dijetHH_Radion_subtrHH");
	if(sighyp == 1 && !bSubtr) outFile = string("dijetHH_GravitonHH");
	else if(sighyp == 1 && bSubtr) outFile = string("dijetHH_Graviton_subtrHH");
	cout << "here " <<  outFile.c_str()  << endl;

	string sOutFile = "MiniTrees/Signal_HH_13TeV/" + outFile + Form("OUT%d_miniTree.root", Mass);
 
	cout << sOutFile.c_str() << endl;
	TFile* f1 = new TFile(sOutFile.c_str(), "recreate");
       

	TTree *TCVARS = new TTree("TCVARS", "hh selection");
	TCVARS->Branch("mgg",&mgg,"mgg/D");
       
	TCVARS->Branch("evWeight",&evWeight,"evWeight/D");
	TCVARS->Branch("normWeight",&normWeight,"normWeight/D");
 
     
	TCVARS->Branch("categories",&categories,"categories/I");
       
   
	for (int iCat = 0; iCat < 3; iCat++){
	  TH1D* hMass = (TH1D*) file0->Get(Form("Radion%d_cat%d;1",Mass,iCat));
	  if (sighyp == 1) hMass = (TH1D*) file0->Get(Form("Graviton%d_cat%d;1",Mass,iCat));
	 
	  TAxis* Axis =   hMass->GetXaxis();
	  for (int i = 1 ; i < hMass->GetNbinsX()+1; i++){
	    int N = hMass->GetBinContent(i);
	   
	    if (i%10 == 0) cout << "i = " << i << "N = " << N << " binCenter = " << hMass->GetBinCenter(i) << endl;
	   
	    mgg = Axis->GetBinCenter(i);
	   
	    normWeight = N;
	    categories = iCat;
	    if (N > 1e-10 && mgg > 999 && mgg < 5000) TCVARS->Fill();
	  }
	}
   
	TCVARS->Write();
	f1->Close();

	cout << "bla 2" << endl;
   
     
      }
    }
 
 
    file0->Close();

  }

}


