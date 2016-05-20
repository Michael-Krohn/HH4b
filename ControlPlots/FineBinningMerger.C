{
  double rescale = 1.;//3368/4528.2; // Factor for data/MC matching
  double Lumi = 2.69*rescale;
  double HPHP = 1.;//0.979*0.979;
  double HPLP = 1;



  TFile *_file0 = TFile::Open("HH4b_subjetBTagged_76X_May/ControlPlots_1GeV.root");

  string stech[2] = {string(""), string("_subtr")};

  for (int itech = 0; itech<2; itech++){

    TFile *_fileOut = TFile::Open(Form("HH4b_subjetBTagged_76X_May/MassPlotFineBins%s_Moriond_Silver.root", stech[itech].c_str()), "RECREATE");
 



    //====================================================================================//
    
    int masses[8] = {1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000};
    string categories[3] = {string("4btag"), string("3btagExactly_HPHP"), string("2btagExactly_HPHP")};
    

    for (int icat = 0; icat < 3; icat++){

      TH1F* TotalMass1GeV_Data = (TH1F*) _file0->Get(Form("TotalMass1GeV%s_%s_Data;1", stech[itech].c_str(), categories[icat].c_str()));
      TotalMass1GeV_Data->Write(Form("Data_cat%d", icat));
    

      for (int i = 0; i<8; i++){
	TH1F* TotalMass1GeV_Radion_13TeV = (TH1F*)   _file0->Get(Form("TotalMass1GeV%s_%s_Radion_m%d_13TeV;1", stech[itech].c_str(), categories[icat].c_str(), masses[i]));
	TH1F* TotalMass1GeV_Graviton_13TeV = (TH1F*)   _file0->Get(Form("TotalMass1GeV%s_%s_Graviton_m%d_13TeV;1",  stech[itech].c_str(), categories[icat].c_str(), masses[i]));
	
	TotalMass1GeV_Radion_13TeV->Scale(HPHP);
	TotalMass1GeV_Graviton_13TeV->Scale(HPHP);
	
	TotalMass1GeV_Radion_13TeV->Write(Form("Radion%d_cat%d", masses[i], icat));
	TotalMass1GeV_Graviton_13TeV->Write(Form("Graviton%d_cat%d", masses[i], icat));
	
      }
    }
    
    _fileOut->Close();

  }
  _file0->Close();

}
