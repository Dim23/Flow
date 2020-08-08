
#include <iostream>
#include "TF1.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TTree.h"
#include <string>
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TStopwatch.h"

TH1F *hBimp;
TH1I *hNpart, *hNcoll;
TH2F *hBimpvsNpart, *hBimpvsNcoll;
TH1F *HphiGen = new TH1F("HphiGen", "HphiGen", 10000, -0.5, 6.5);
TRandom3 *ar;
TRandom3 *br;

TFile *fi;

double v1 = 0, v2 = 0.2, v3 = 0, v4 = 0, v5 = 0, v6 = 0;
static const float PI = 3.1415926535;
static const float pTmax = 4900.;
static const float ptMax = 5.0;

float calc_v2(double b, double eta, double pt)
{
  float a1, a2, a3, a4;
  a1 = 0.4397 * exp(-(b - 4.526) * (b - 4.526) / 72.0) + 0.636;
  a2 = 1.916 / (b + 2) + 0.1;
  a3 = 4.79 * 0.0001 * (b - 0.621) * (b - 10.172) * (b - 23) + 1.2; // this is >0 for b>0
  a4 = 0.135 * exp(-0.5 * (b - 10.855) * (b - 10.855) / 4.607 / 4.607) + 0.0120;

  float temp1 = pow(pt, a1) / (1 + exp((pt - 3.0) / a3));
  float temp2 = pow(pt + 0.1, -a2) / (1 + exp(-(pt - 4.5) / a3));
  float temp3 = 0.01 / (1 + exp(-(pt - 4.5) / a3));

  //v2 = (a4 * (temp1 + temp2) + temp3) * exp (-0.5 * eta * eta / 6.27 / 6.27);

  // Adjust flow rapidity dependence to better match PHOBOS 200 GeV Au+Au data
  // JGL 9/9/2019
  // See JS ToG talk at https://indico.bnl.gov/event/6764/

  v2 = (a4 * (temp1 + temp2) + temp3) * exp(-0.5 * eta * eta / 2.0 / 2.0);

  return v2;
}

// New parameterization for vn
void jjia_minbias_new(double b, double eta, double pt)
{
  v2 = calc_v2(b, eta, pt);

  float fb = 0.97 + 1.06 * exp(-0.5 * b * b / 3.2 / 3.2);
  v3 = pow(fb * sqrt(v2), 3);

  float gb = 1.096 + 1.36 * exp(-0.5 * b * b / 3.0 / 3.0);
  gb = gb * sqrt(v2);
  v4 = pow(gb, 4);
  v5 = pow(gb, 5);
  v6 = pow(gb, 6);
  v1 = 0;
}

double F(double a, double pt)
{
  return (1 + 2 * v1 * cos(a) + 2 * v2 * cos(2 * (a)) + 2 * v3 * cos(3 * (a)) + 2 * v4 * cos(4 * (a)) + 2 * v5 * cos(5 * (a)) + 2 * v6 * cos(6 * (a))) / (2 * PI);
}
double ra, rb;
double FGen(double a, double b, double pt)
{
  ra = 2 * PI * a;
  rb = (1 + 2 * v1 + 2 * v2 + 2 * v3 + 2 * v4 + 2 * v5 + 2 * v6) / (2 * PI) * b;
  while (F(ra, pt) < rb)
  {
    rb = (1 + 2 * v1 + 2 * v2 + 2 * v3 + 2 * v4 + 2 * v5 + 2 * v6) / (2 * PI) * br->Rndm();
    ra = 2 * PI * ar->Rndm();
  }
  return ra;
}

double GetProjectedRandom(double b, TH2F *const &hist)
// Projects Npart or Ncoll distribution for given impact parameter and gives you random Npart or Ncoll
// using the projected distribution as a probability function (i.e. gives you random number within RMS)
{
  TH1D *hist1D = hist->ProjectionY("tmp", hist->GetXaxis()->FindBin(b), hist->GetXaxis()->FindBin(b));
  double result = hist1D->GetRandom();
  return result;
}

double *DrawMult(TString glauberFileName = "/home/dim2/FLOW5/merge_hist_glaub_200gev.root")
{
  const double f = 0.1, MeanMult = 250.;
  //const int Nevents = 2e5;

  hBimp = (TH1F *)fi->Get("hBimp");
  hNpart = (TH1I *)fi->Get("hNpart");
  hNcoll = (TH1I *)fi->Get("hNcoll");
  hBimpvsNpart = (TH2F *)fi->Get("hBimpvsNpart");
  hBimpvsNcoll = (TH2F *)fi->Get("hBimpvsNcoll");

  double b = hBimp->GetRandom();
  double Npart = GetProjectedRandom(b, hBimpvsNpart);
  double Ncoll = GetProjectedRandom(b, hBimpvsNcoll);
  double Na = f * Npart + (1 - f) * Ncoll;
  double NaNorm = (f * hNpart->GetMean() + (1 - f) * hNcoll->GetMean()) / MeanMult;
  int mult = (Na / NaNorm);
  double *bmult = new double[2];
  bmult[0] = b;
  bmult[1] = mult;
  return bmult;
}

double dndpT(double pT)
{
  double temp;
  double pTmin = 100.;
  double pTmax = 4900.;
  double pT0 = 500.;
  double pT1 = 3000.;
  double T0 = 300.;
  if (pT < pTmin)
    temp = 0;
  else if (pT <= pT0)
    temp = 1;
  else if (pT <= pT1)
    temp = exp(-(pT - pT0) / T0);
  else
    temp = exp(-(pT1 - pT0) / T0) * pow((pT1 / pT), 7.);
  return temp;
}

void generator(const char *OUTfile_name = "~/FLOW5/OUT/Vin1nonflow.root", int run0 = 1000000, TString glauberFileName = "/home/dim2/FLOW5/merge_hist_glaub_200gev.root")
{


  TRandom3 rnd;
  ar = new TRandom3();
  br = new TRandom3();
  TRandom3 *TRpsi = new TRandom3();
  TRandom3 *TReta = new TRandom3();
  TRandom3 *TRpt = new TRandom3();

  fi = new TFile(glauberFileName.Data(), "read");

  rnd.SetSeed(0);
  gRandom->SetSeed(0);
  ar->SetSeed(0);
  br->SetSeed(0);

  TRpsi->SetSeed(0);
  TReta->SetSeed(0);
  TRpt->SetSeed(0);
cout<<"ar "<<ar->Rndm()<<" br "<<br->Rndm()<<" rnd "<<rnd.Rndm()<<" TRpsi "<<TRpsi->Rndm()<<" TReta "<<TReta->Rndm()<<" TRpt "<<TRpt->Rndm()<<endl;
cout<<"ar "<<ar->Rndm()<<" br "<<br->Rndm()<<" rnd "<<rnd.Rndm()<<" TRpsi "<<TRpsi->Rndm()<<" TReta "<<TReta->Rndm()<<" TRpt "<<TRpt->Rndm()<<endl;
cout<<"ar "<<ar->Rndm()<<" br "<<br->Rndm()<<" rnd "<<rnd.Rndm()<<" TRpsi "<<TRpsi->Rndm()<<" TReta "<<TReta->Rndm()<<" TRpt "<<TRpt->Rndm()<<endl;
cout<<"ar "<<ar->Rndm()<<" br "<<br->Rndm()<<" rnd "<<rnd.Rndm()<<" TRpsi "<<TRpsi->Rndm()<<" TReta "<<TReta->Rndm()<<" TRpt "<<TRpt->Rndm()<<endl;

  TFile *f0 = new TFile(OUTfile_name, "recreate");
  f0->cd();

  Int_t event;
  Int_t Nch = 0;
  Float_t b;
  Float_t psi_RP;
  Float_t phi[1500];
  Bool_t flow[1500];
  Float_t dPHI[1500];
  Float_t meta[1500];
  Float_t mpt[1500];
  Int_t nonflow = 2.;
  Float_t nonflowrate = 0.2;

  TTree *tree = new TTree("tree", "My Tree");
  //tree->Branch("EVENT", &event, "event/I");
  tree->Branch("nh", &Nch, "Nch/I");
  tree->Branch("b", &b, "b/F");
  tree->Branch("rp", &psi_RP, "psi_RP/F");
  tree->Branch("phi0", &phi, "phi[Nch]/F");
  //tree->Branch("dPHI", &dPHI, "dPHI[Nch]/F");
  tree->Branch("eta", &meta[0], "meta[Nch]/F");
  tree->Branch("pt", &mpt[0], "mpt[Nch]/F");
  tree->Branch("bFlow", &flow[0], "flow[Nch]/O");

  float A = PI / 3.;
  float B = PI / 2.;
  float C = PI;
  float D = 5. * PI / 4.;

  double *bmult = new double[2];

  for (int k = 0; k < run0; k++)
  {
    if (k % 10000 == 0)
    {
      cout << k << " ivents were simulated,the process is " << k * 100 / run0 << "% complete" << endl;
    }
    bmult = DrawMult(glauberFileName);
    b = (float)(bmult[0]);
    Nch = (int)(bmult[1]);
    /*while(Nch<20){bmult=DrawMult(glauberFileName);
    b=(float)(bmult[0]);
    Nch=(int)(bmult[1]);}*/

    psi_RP = (TRpsi->Rndm()) * 2 * PI;
    for (int j = 0; j < Nch; j++)
    {

    momentum:
      mpt[j] = pTmax * (TRpt->Rndm());
      if (rnd.Rndm() > dndpT(mpt[j]))
        goto momentum;
      mpt[j] = mpt[j] / 1000.;

      meta[j] = 4 * TReta->Rndm() - 2;
    }
    for (int d = 0; d < Nch; d++)
    {
      jjia_minbias_new(b, meta[d], mpt[d]);

      dPHI[d] = FGen(ar->Rndm(), br->Rndm(), mpt[d]);
      phi[d] = dPHI[d] + psi_RP;
      flow[d] = 1;
      if (phi[d] > 2 * PI)
      {
        phi[d] = phi[d] - 2 * PI;
      }

     //while((phi[d]>=A && phi[d]<=B) || (phi[d]>=C && phi[d]<=D)){dPHI[d]=FGen(ar->Rndm(),br->Rndm(),mpt[d]);phi[d]=dPHI[d]+psi_RP;if(phi[d]>2*PI){phi[d]=phi[d]-2*PI;} }

      if(nonflow>0 && (d+nonflow-1)<Nch && TRpt->Rndm()<nonflowrate){for(int kk=1;kk<nonflow;kk++){d++;mpt[d]=mpt[d-1];dPHI[d]=dPHI[d-1];phi[d]=phi[d-1];flow[d]=0;/*meta[d]=meta[d-1];*/ };};

      HphiGen->Fill(phi[d]);
    }
    event = k;
    tree->Fill();
  }
  tree->Write();
  HphiGen->Write();
  f0->Close();
  cout << "DONE generator" << endl;
}
