#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TMultiGraph.h>
#include <TRandom.h>
#include <TLegend.h>
#include <TLegendEntry.h>

#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"
#include "ComponentAnalyticField.hh"
#include "GeometrySimple.hh"
#include "Sensor.hh"
#include "ViewField.hh"
#include "SolidTube.hh"
#include "MediumGas.hh"
#include "ViewCell.hh"
#include "SolidBox.hh"
#include "AvalancheMicroscopic.hh"
#include "ViewSignal.hh"
#include "ViewDrift.hh"
#include "ViewMedium.hh"
#include "TrackElectron.hh"
#include "Medium.hh"
#include "ViewCell.hh"
#include "Track.hh"
using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
 
  const double pressure = 1 * AtmosphericPressure;
  const double temperature = 273.15;
 
  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetTemperature(temperature);
  gas->SetPressure(pressure);
  gas->SetComposition("Ar", 66.4, "iC4H10", 33.6);

  //**ONLY NEED THIS SECTION IF YOU NEED TO SETUP A NEW GASTABLE**//


  //Set the field range to be covered by the gas table. 
  const int nFields = 20;
  const double emin =    100.;
  const double emax = 10000.;
  //Flag to request logarithmic spacing.
  const bool useLog = true;
  gas->SetFieldGrid(emin, emax, nFields, useLog); 

  const int ncoll = 10;
  // Switch on debugging to print the Magboltz output.
  gas->DisableDebugging();
  // Run Magboltz to generate the gas table.
  gas->GenerateGasTable(ncoll);
  //gas->DisableDebugging();
  //Save the table. 
  gas->WriteGasFile("ar_66_ic4h10_34.gas");


  //If using a different gas, will need to load the correct gasfile//
  gas->LoadGasFile("ar_66_ic4h10_34.gas");


  ComponentAnalyticField* cmp = new ComponentAnalyticField();

  //Define radius of sense wire, diameter of field wire, dimensions of tube//
  const double rWire= 0.00125, dfWire=0.015;
  const double rTube= 1; 
  const double lTube= 53.5; 

  //Setup geometry//
  GeometrySimple* geo = new GeometrySimple();
  SolidBox* tube = new SolidBox(0.,0.,0.,rTube,rTube,lTube);
  //Add tube and gas to geometry//
  geo->AddSolid(tube, gas);
  //Add geometry to analytic field//
  cmp->SetGeometry(geo);

  //Setup and add wires and tube to the field, may not need a tube//
  const double vWire= 2700.;
  //const double vTube = 0.;
  cmp->AddWire(0.,0., 2*rWire, vWire, "s");
  cmp->AddWire(rTube/2,0, dfWire, -1900, "g1");
  cmp->AddWire(-rTube/2,0, dfWire,-1900, "g1");
  cmp->AddWire(0,rTube/2, dfWire, -2400, "g4");
  cmp->AddWire(0,-rTube/2, dfWire, -2400, "g4");
  cmp->AddWire(.5,.5, dfWire, 0., "g5");
  cmp->AddWire(0.5,-.5, dfWire, 0., "g5");
  cmp->AddWire(-0.5,.5, dfWire, 0., "g5");
  cmp->AddWire(-0.5,-.5, dfWire, 0., "g5");
  //Repeat components every 1 unit (cm)//
  cmp->SetPeriodicityY(1);
 
  
  cmp->AddReadout("s");

  //Create sensor and add component field to it, not always necessary unless calculating signals//
  Sensor* sensor=new Sensor();
  sensor->AddComponent(cmp);
  sensor->AddElectrode(cmp, "s");

  //Visualize the field//
  ViewField* field= new ViewField();
  field->SetSensor(sensor);
  //May need to change the range//
  field->SetElectricFieldRange(0, 10000.);
  field->PlotContour("e");
  //Prepare to make a plot of drift velocity versus electric field//
  double efield[50], dvel[50];
  double vx,vy,vz;
  TGraph *drifefie= new TGraph();
  
  for (int i=0; i<50; i++){
    //loop over positions you are interested in//
    double pos[2]= {0.01*i,0.01*i};

    //  Not sure how what[] comes into play, but it is needed to tell 
    //  EvaluatePotential whether to get electric field components, e-field 
    //  mag. or get the potential. EvalutatePotential only reads what[0]. 
    //
    //     if (what[0] > 30.) {
    //        return Z-component of electric field;
    //      } else if (what[0] > 20.) {
    //        return Y-component of electric field;
    //      } else if (what[0] > 10.) {
    //        return X-component of electric field;
    //      } else if (what[0] > 0.) {
    //        return Magnitude of electric field;
    //      }
    //     return potential;

    double what[3]={11,1,2};
    double ex=(field->EvaluatePotential(pos, what));
    what[0]=what[0]*2;
    double ey= (field->EvaluatePotential(pos, what));
    what[0]=what[0]*3;
    double ez=(field->EvaluatePotential(pos, what));

//  ElectronVelocity only returns a bool. However, it does store the components 
//  of drift velocity at that point in vx,vy,vz
    gas->ElectronVelocity(ex, ey, ez, 0., 0., 0., vx,vy,vz);

    //Store the e-field and drift velocity magnitudes//
    efield[i]=sqrt(ex*ex+ey*ey+ez*ez);
    dvel[i]=sqrt(vx*vx+vy*vy+vz*vz);
  
    //Plot the points//
    if (efield[i]<10000.){
      drifefie->SetPoint(i, efield[i], dvel[i]);
    }
  }

  //NEXT COMMENTS NECESSARY ONLY IF YOU WANT TO COMPARE TO THE GASTABLE RESULTS
  //ALSO THESE MUST BE CHANGED IF USING A GAS OTHER THAN Ar-C2H6 50/50.

  // Double_t xgpoints[20]={ 1.31578947e-01, 1.67667761e-01, 2.13654834e-01, 2.72255011e-01, 3.46927750e-01, 4.42081353e-01, 5.63333210e-01, 7.17841419e-01, 9.14727363e-01, 1.16561420e+00, 1.48531302e+00, 1.89269722e+00, 2.41181672e+00, 3.07331772e+00, 3.91625190e+00, 4.99038183e+00, 6.35911873e+00, 8.10326594e+00, 1.03257891e+01, 1.31578947e+01};

  //  Double_t ygpoints[20]={ 2.46561830e+00, 2.98075892e+00, 3.48895177e+00, 3.95223049e+00, 4.34517668e+00, 4.66439922e+00, 4.91746082e+00, 5.11043863e+00, 5.24838191e+00, 5.32293867e+00, 5.32005531e+00, 5.23377813e+00, 5.07646060e+00, 4.89688261e+00, 4.74778419e+00, 4.66468487e+00, 4.63585997e+00, 4.65042164e+00, 4.69895847e+00, 4.80436817e+00}; 

  //  TGraph* garf=new TGraph();
  //  for (Int_t i=0; i<18; i++){
  //    garf->SetPoint(i, xgpoints[i]*1000, ygpoints[i]/1000);
  //  }



  TMultiGraph* mg= new TMultiGraph();
  mg->Add(drifefie);
  //mg->Add(garf);
  drifefie->SetMarkerStyle(20);
  drifefie->SetMarkerSize(1);
  //garf->SetMarkerStyle(21); garf->SetMarkerColor(2);
  //mg->Draw("AP");
 
  // ViewDrift is if you want to view the drift lines of the electron\s //
  
  ViewDrift* drift= new ViewDrift();

  // Setup time //
  const double tmin=0., tmax=1., tstep=.000001;
  const int ntimebins=int((tmax-tmin)/tstep);
  sensor->SetTimeWindow(0., tstep, ntimebins);
  
  // Use to simulate either an avalanche of electrons or just one single drifting
  // electron. For drift->Plot(bool, bool):
  //        Plot(true, true)= 2-D plot
  //        Plot(false, true)=3-D plot

  // AvalancheMicroscopic* aval= new AvalancheMicroscopic();
  // aval->EnablePlotting(drift);
  // aval->SetSensor(sensor);
  // aval->EnableSignalCalculation();
  // aval->SetTimeWindow(0.,1000.);
  // const double x0=0.4  , y0=0., z0=0., t0=0.;
  // aval->EnableAvalancheSizeLimit(150);
  // aval->AvalancheElectron(x0, y0, z0, t0,0.);
  // drift->Plot(true,true);
  

// Yet to be correctly implemented due to C2H6 not being in their table of cross
// sections. TrackElectron simulates a high-energy electron traversing the cell

  TRandom *rand=new TRandom(0);
  TrackElectron* elec=new TrackElectron();
  elec->SetSensor(sensor);
  elec->EnablePlotting(drift);
  elec->SetEnergy(4e9);
  TH1F *dist= new TH1F("distance", "distance", 20,0,0.5);
  AvalancheMicroscopic* aval= new AvalancheMicroscopic();
  aval->EnablePlotting(drift);
  aval->EnableDistanceHistogramming(1);
  aval->SetDistanceHistogram(dist, 'r');
  aval->SetSensor(sensor);
  aval->EnableSignalCalculation();
  aval->SetTimeWindow(0.,1000.);
  double xcls, ycls, zcls, tcls, ecls, extra;
  int ncls, count=0, place=0;
  double xi,yi,zi,ti,ei,xf,yf,zf,tf,ef;
  int state,spot=0;
  int min=0;
  TArrayD *time= new TArrayD(), *xpos= new TArrayD(), *smalt= new TArrayD(), *minx=new TArrayD();
  for (Int_t i=0; i<10000; i++){
    elec->NewTrack(-1, (rand->Uniform(-.49,.49)), 0., 0., 1,0,0);
    place++;
    std::cout<<place<<"->loop iteration \n";
    while (elec->GetCluster(xcls, ycls, zcls, tcls, ncls, ecls, extra)){
      double location=sqrt(xcls*xcls+ycls*ycls+zcls*zcls);
      //std::cout<<count<<std::endl;
      if (location<=rTube/2){
	aval->EnableAvalancheSizeLimit(20);
	aval->DriftElectron(xcls, ycls, zcls, tcls, ecls);
	aval->GetElectronEndpoint(0,xi,yi,zi,ti,ei,xf,yf,zf,tf,ef,state);
	if (sqrt(xf*xf+yf*yf)<=rWire){
	  time->Set(count+1);
	  time->AddAt((tf-ti),count);
	  xpos->Set(count+1);
	  xpos->AddAt(xi,count);
	  count++;
	}
	
      }
    }
    int a=spot;
    //std::cout<<(time->GetSize())-spot<<std::endl;
    for (Int_t p=0; p<((time->GetSize())-spot);p++){
      if ((time->GetAt(a))>(time->GetAt(spot+p))){
    	a=(spot+p);
      }
    }
    min=a;
    smalt->Set(i+1); minx->Set(i+1);
    smalt->AddAt(time->GetAt(min), i); minx->AddAt(xpos->GetAt(min), i);
  
    spot=count;
    //a=spot+1;
    //std::cout<<xpos->GetAt(min)<<std::endl;
    // std::cout<<smalt->GetAt(i)<<std::endl;
  }
  
  //drift->Plot(true,true);
  TCanvas *can=new TCanvas("DriftTime", "Drift Time", 600, 600);
  TCanvas *can2=new TCanvas("DriftPos", "Drift Position", 600, 600);
  dist->SetBinContent(1,0.);
  //dist->Draw();
  // TH1D *times=new TH1D("DriftsTime", "Drift Times", 75,0,150);
  // for (int k; k<time->GetSize(); k++){
  //   times->Fill(time->GetAt(k));
  // }
  // times->Draw();

  TH1D *mintimes=new TH1D("MinTime", "Minimum Drift Times", 75,0,150), *minxpos=new TH1D("minxpos", "X-pos of min. Drift Time",50 ,-.1, .1);
  for (int k; k<smalt->GetSize(); k++){
    mintimes->Fill(smalt->GetAt(k));
    minxpos->Fill(minx->GetAt(k));
  }
  can->cd();
  mintimes->Draw();
  can2->cd();
  minxpos->Draw();
  

  //** Makes plot of drift time vs x-position**// 

  // TGraph *tvsy=new TGraph();
  // for (Int_t j=0; j<time->GetSize(); j++){
  //   tvsy->SetPoint(j, xpos->GetAt(j), time->GetAt(j));
  // }
  // tvsy->SetMarkerSize(1);
  // tvsy->SetMarkerStyle(20);
  // tvsy->Draw("AP");

  //can->SaveSource("minimum_drift_times");

  //DON'T GET RID OF THIS OR YOUR STUFF WON'T WORK//
  app.Run(kTRUE);

}
