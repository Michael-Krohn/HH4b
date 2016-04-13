#include </usr/local/root/root-6.04.02/include/TGraph.h>
#include <string>
#include <map>
#include <iostream>

using namespace std;

void getVal(double mjj, std::string sample, std::string channel, double& alpha, double& n, double& sigma, double& mean) {

  const int np(6);
  
  double masses[np] = { 1200.0,1400.0,1600.0,1800.0,2000.0,2500.0 } ;

  std::map<std::string, TGraph*> gr;

  if (sample == "Radion") {
    if (channel == "3b") {
      double alpha3b[np] = { 2.99992869473,1.23572089917,1.38704596759,1.43790885737,1.47418666312,1.11961838575 } ; 
      double n3b[np] = { 127.58280731,134.559662569,126.077137231,138.005844913,137.205261077,142.861915389 } ; 
      double sigma3b[np] = { 49.6014843743,51.8530982105,59.9136530318,67.5958048863,75.9686650749,84.6657692385 } ; 
      double mean3b[np] = { 1138.7518106,1329.27922894,1512.35492927,1702.76282505,1890.47142027,2363.48968532 } ; 
      gr[sample+"_"+"alpha"+channel] = new TGraph(np, masses, alpha3b); 
      gr[sample+"_"+"n"+channel] = new TGraph(np, masses, n3b); 
      gr[sample+"_"+"sigma"+channel] = new TGraph(np, masses, sigma3b); 
      gr[sample+"_"+"mean"+channel] = new TGraph(np, masses, mean3b); 
    }
    else if (channel == "4b") {
      double alpha4b[np] = { 2.99994285083,1.46588196934,1.32019487083,1.24987586528,1.166365396,1.13823997853 } ; 
      double n4b[np] = { 126.460030265,132.776297381,132.018630891,121.374740807,144.337572861,144.868422011 } ; 
      double sigma4b[np] = { 45.7033847362,51.1543513626,58.3791610108,62.8435240971,66.7409003893,82.8829569711 } ; 
      double mean4b[np] = { 1141.97051876,1328.64384766,1516.84908692,1706.22166888,1892.21279472,2359.89977439 } ; 
      gr[sample+"_"+"alpha"+channel] = new TGraph(np, masses, alpha4b); 
      gr[sample+"_"+"n"+channel] = new TGraph(np, masses, n4b); 
      gr[sample+"_"+"sigma"+channel] = new TGraph(np, masses, sigma4b); 
      gr[sample+"_"+"mean"+channel] = new TGraph(np, masses, mean4b); 
    } 
    else std::cout << "Error: check channel\n"; 
  }
  else if (sample == "Graviton") {
    if (channel == "3b") {
      double alpha3b[np] = { 1.97782922403,1.41896083125,1.2617805834,1.34539670468,1.31276457424,1.32940541309 } ; 
      double n3b[np] = { 141.139758152,117.887448804,129.516555916,128.842102297,128.751760353,131.820871684 } ; 
      double sigma3b[np] = { 47.2594027067,53.1273966322,60.8086392999,63.8498759854,73.2388044796,87.1485163363 } ; 
      double mean3b[np] = { 1139.44448398,1326.1356103,1515.97716079,1701.67066795,1889.79779805,2359.87112233 } ; 
      gr[sample+"_"+"alpha"+channel] = new TGraph(np, masses, alpha3b); 
      gr[sample+"_"+"n"+channel] = new TGraph(np, masses, n3b); 
      gr[sample+"_"+"sigma"+channel] = new TGraph(np, masses, sigma3b); 
      gr[sample+"_"+"mean"+channel] = new TGraph(np, masses, mean3b); 
    }
    else if (channel == "4b") {
      double alpha4b[np] = { 2.9999986789,1.32495487888,1.23799526302,1.16954277446,1.21623722849,1.25469474452 } ; 
      double n4b[np] = { 137.838630302,126.095947032,138.035633381,146.013900331,132.830584043,145.834851327 } ; 
      double sigma4b[np] = { 45.6843979549,51.0787006165,57.3445271952,60.7986544462,68.0860896875,85.7909735223 } ; 
      double mean4b[np] = { 1141.24023757,1330.25308396,1517.83383081,1703.26522902,1890.69263049,2359.9586752 } ; 
      gr[sample+"_"+"alpha"+channel] = new TGraph(np, masses, alpha4b); 
      gr[sample+"_"+"n"+channel] = new TGraph(np, masses, n4b); 
      gr[sample+"_"+"sigma"+channel] = new TGraph(np, masses, sigma4b); 
      gr[sample+"_"+"mean"+channel] = new TGraph(np, masses, mean4b); 
    } 
    else std::cout << "Error: check channel\n"; 
  }
  else std::cout << "Error: check sample\n"; 

  alpha = gr[sample+"_alpha"+channel] -> Eval(mjj) ; 
  n = gr[sample+"_n"+channel] -> Eval(mjj) ; 
  sigma = gr[sample+"_sigma"+channel] -> Eval(mjj) ; 
  mean = gr[sample+"_mean"+channel] -> Eval(mjj) ; 

  return;

}

int extrap () {

  double mjj(1200.), alpha(0), n(0), sigma(0), mean(0); 
  std::string sample("Graviton"), channel("3b") ;
  getVal(mjj, sample, channel, alpha, n, sigma, mean);
  std::cout << " mjj = " << mjj 
            << " alpha = " << alpha 
            << " n = " << n 
            << " sigma = " << sigma 
            << " mean = " << mean 
            << std::endl ; 

  return 0;
}
