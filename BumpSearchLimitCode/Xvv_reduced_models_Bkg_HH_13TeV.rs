mgg[1000,3000];


HH_sig_m0_4btag_cat0[2000.0, 900.0, 3100.0];
HH_sig_sigma_4btag_cat0[40, 20.0, 1000.0];
HH_sig_alpha_4btag_cat0[ 0.8, 0.0, 3.0]; 
HH_sig_n_4btag_cat0[130, 0.00001, 1000.0]; 



HH_sig_4btag_cat0   = CBShape(mgg, HH_sig_m0_4btag_cat0, HH_sig_sigma_4btag_cat0, HH_sig_alpha_4btag_cat0, HH_sig_n_4btag_cat0);


bkg_fit_slope_4btag_cat0[1000.0,0, 10000000];
bkg_fit_slope1_4btag_cat0[10., -100.0, 100.0];
bkg_fit_slope1_clone_4btag_cat0[10., -100.0, 100.0];
bkg_fit_slope1_power_4btag_cat0[2., -8.0, 8.0];



HH_sig_m0_3btag_cat1[2000.0, 900.0, 3100.0];
HH_sig_sigma_3btag_cat1[100, 20.0, 1000.0];
HH_sig_alpha_3btag_cat1[ 0.8, 0.0, 3.0]; 
HH_sig_n_3btag_cat1[130, 0.00001, 1000.0]; 


HH_sig_3btag_cat1      = CBShape(mgg, HH_sig_m0_3btag_cat1, HH_sig_sigma_3btag_cat1, HH_sig_alpha_3btag_cat1, HH_sig_n_3btag_cat1);




bkg_fit_slope_3btag_cat1[1000.0,0, 10000000];
bkg_fit_slope1_3btag_cat1[10., -100.0, 100.0];
bkg_fit_slope1_clone_3btag_cat1[10., -100.0, 100.0];
bkg_fit_slope1_power_3btag_cat1[2., -8.0, 8.0];
bkg_fit_slope2_3btag_cat1[0., -500.0, 500.0];




HH_sig_m0_2btag_cat2[2000.0, 900.0, 3100.0];
HH_sig_sigma_2btag_cat2[100, 20.0, 1000.0];
HH_sig_alpha_2btag_cat2[ 0.8, 0.0, 3.0]; 
HH_sig_n_2btag_cat2[130, 0.00001, 1000.0]; 


HH_sig_2btag_cat2     = CBShape(mgg, HH_sig_m0_2btag_cat2, HH_sig_sigma_2btag_cat2, HH_sig_alpha_2btag_cat2, HH_sig_n_2btag_cat2);



bkg_fit_slope_2btag_cat2[1000.0,0, 10000000];
bkg_fit_slope1_2btag_cat2[10., -100.0, 100.0];
bkg_fit_slope1_clone_2btag_cat2[10., -100.0, 100.0];
bkg_fit_slope1_power_2btag_cat2[2., -8.0, 8.0];
bkg_fit_slope2_2btag_cat2[0., -500.0, 500.0];


wei[1,0,10];

sqrtS[13000., 13000., 13000.]
