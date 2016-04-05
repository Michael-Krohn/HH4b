mgg[1000,3000];


jj_HH_sig_m0_CMS_jj_4btag_cat0[2000.0, 900.0, 3100.0];
jj_HH_sig_sigma_CMS_jj_4btag_cat0[50, 20.0, 1000.0];
jj_HH_sig_alpha_CMS_jj_4btag_cat0[ 0.8, 0.0, 3.0]; 
jj_HH_sig_n_CMS_jj_4btag_cat0[130, 0.00001, 1000.0]; 



HH_jj_CMS_jj_4btag_cat0   = CBShape(mgg, jj_HH_sig_m0_CMS_jj_4btag_cat0, jj_HH_sig_sigma_CMS_jj_4btag_cat0, jj_HH_sig_alpha_CMS_jj_4btag_cat0, jj_HH_sig_n_CMS_jj_4btag_cat0);


bkg_fit_slope_CMS_jj_4btag_cat0[1000.0,0, 10000000];
bkg_fit_slope1_CMS_jj_4btag_cat0[-50., -200.0, 0.0];




jj_HH_sig_m0_CMS_jj_3btag_HPHP_cat1[2000.0, 900.0, 3100.0];
jj_HH_sig_sigma_CMS_jj_3btag_HPHP_cat1[50, 20.0, 1000.0];
jj_HH_sig_alpha_CMS_jj_3btag_HPHP_cat1[ 0.8, 0.0, 3.0]; 
jj_HH_sig_n_CMS_jj_3btag_HPHP_cat1[130, 0.00001, 1000.0]; 


HH_jj_CMS_jj_3btag_HPHP_cat1      = CBShape(mgg, jj_HH_sig_m0_CMS_jj_3btag_HPHP_cat1, jj_HH_sig_sigma_CMS_jj_3btag_HPHP_cat1, jj_HH_sig_alpha_CMS_jj_3btag_HPHP_cat1, jj_HH_sig_n_CMS_jj_3btag_HPHP_cat1);




bkg_fit_slope_CMS_jj_3btag_HPHP_cat1[1000.0,0, 10000000];
bkg_fit_slope1_CMS_jj_3btag_HPHP_cat1[-50., -200.0, 0.0];





jj_HH_sig_m0_CMS_jj_3btag_HPLP_cat2[2000.0, 900.0, 3100.0];
jj_HH_sig_sigma_CMS_jj_3btag_HPLP_cat2[50, 20.0, 1000.0];
jj_HH_sig_alpha_CMS_jj_3btag_HPLP_cat2[ 0.8, 0.0, 3.0]; 
jj_HH_sig_n_CMS_jj_3btag_HPLP_cat2[130, 0.00001, 1000.0]; 


HH_jj_CMS_jj_3btag_HPLP_cat2     = CBShape(mgg, jj_HH_sig_m0_CMS_jj_3btag_HPLP_cat2, jj_HH_sig_sigma_CMS_jj_3btag_HPLP_cat2, jj_HH_sig_alpha_CMS_jj_3btag_HPLP_cat2, jj_HH_sig_n_CMS_jj_3btag_HPLP_cat2);



bkg_fit_slope_CMS_jj_3btag_HPLP_cat2[1000.0,0, 10000000];
bkg_fit_slope1_CMS_jj_3btag_HPLP_cat2[-50., -200.0, 0.0];



wei[1,0,10];

sqrtS[13000., 13000., 13000.]
