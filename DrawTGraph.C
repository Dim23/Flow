#include <TProfile.h>
#include <iostream>
#include <fstream>

// Draws 2 TGraphErrors (upper panel) with their gr1/gr2 ratio (lower pannel)
TCanvas *DrawTGraph(TGraphErrors *const &gr1, TGraphErrors *const &gr2, TString str="", 
                    Double_t yRatio_low=0.89, Double_t yRatio_high=1.11)
{
  // Setting up global variables for the plot
  gROOT->SetStyle("Pub");
	gROOT->ForceStyle();
	gStyle->SetPalette(kDarkRainBow);
	gStyle->SetErrorX(0);

  // Read points
  Double_t *vx_gr1 = gr1->GetX();
  Double_t *vy_gr1 = gr1->GetY();
  Double_t *vx_gr2 = gr2->GetX();
  Double_t *vy_gr2 = gr2->GetY();

  // Read errors
  Double_t *ex_gr1 = gr1->GetEX();
  Double_t *ey_gr1 = gr1->GetEY();
  Double_t *ex_gr2 = gr2->GetEX();
  Double_t *ey_gr2 = gr2->GetEY();

  Int_t n1bins = gr1->GetN();

  // Initialization of the canvas & pads
  TCanvas *canv = new TCanvas(Form("canv"),Form("Canvas"),900,800);
  canv->cd();
  TPad *padUp = new TPad(Form("padUp"),"v2 vs pt",0.,0.33,1.,1.,0,-1,0);
  TPad *padDown = new TPad(Form("padDown"),"Ratio v2",0.,0.,1.,0.33,0,-1,0);

  double padUW;
	double padUH;
	double padDW;
	double padDH;

  padUp->SetBorderSize(0);
  padDown->SetBorderSize(0);
  
  padUp->SetBottomMargin(0.);
  padDown->SetTopMargin(0.005);
  
  padUW = padUp->GetWw()*padUp->GetAbsWNDC();
  padUH = padUp->GetWh()*padUp->GetAbsHNDC();
  padDW = padDown->GetWw()*padDown->GetAbsWNDC();
  padDH = padDown->GetWh()*padDown->GetAbsHNDC();
  
  padUp->Draw();
  padDown->Draw();

  // Setting up markers & colors for TGraphErrors
  gr1->SetMarkerStyle(26);
  gr1->SetMarkerSize(1.6);
  gr1->SetLineColor(kRed);
  gr1->SetMarkerColor(kRed);

  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(1.6);
  gr2->SetLineColor(kBlack);
  gr2->SetMarkerColor(kBlack);

  // Draw TGraphErrors in the upper pad
  padUp->cd();

  gr1->GetXaxis()->SetLimits(0.95*vx_gr1[0],1.05*vx_gr1[n1bins-1]);

  gr1->GetXaxis()->SetLabelSize(0.06);
  gr1->GetYaxis()->SetLabelSize(0.06);
  gr1->GetXaxis()->SetTitleSize(0.07);
  gr1->GetYaxis()->SetTitleSize(0.07);
  gr1->GetYaxis()->SetTitleOffset(1.08);

  gr1->Draw("AP");
  gr2->Draw("P");

  TLegend *leg_pt = new TLegend(0.468,0.04,0.89,0.309);
  leg_pt->SetBorderSize(0);
  leg_pt->SetHeader(str.Data(),"C");
  // leg_pt->SetHeader(Form("Au+Au,#sqrt{s_{NN}}=200 GeV"),"C");
  leg_pt->AddEntry(gr1,Form("%s",gr1->GetTitle()),"p");
  leg_pt->AddEntry(gr2,Form("%s",gr2->GetTitle()),"p");

  leg_pt->Draw();

  //Draw gr1/gr2 ratio in the bottom pad
  padDown->cd();

  std::vector<Double_t> v1X;
	std::vector<Double_t> v1Y;
	std::vector<Double_t> v1Xerr;
	std::vector<Double_t> v1Yerr;
  std::vector<Double_t> v2X;
	std::vector<Double_t> v2Y;
	std::vector<Double_t> v2Xerr;
	std::vector<Double_t> v2Yerr;
  std::vector<Double_t> vRatioY;
	std::vector<Double_t> vRatioYerr;

  for (int i=0; i<gr1->GetN();i++)
  {
    v1X.push_back(vx_gr1[i]);
    v1Y.push_back(vy_gr1[i]);
    v1Xerr.push_back(ex_gr1[i]);
    v1Yerr.push_back(ey_gr1[i]);

    v2Y.push_back((Double_t) gr2->Eval(v1X.at(i),0,"S"));
    v2Yerr.push_back(ey_gr2[i]);

    vRatioY.push_back(v1Y.at(i)/v2Y.at(i));
    vRatioYerr.push_back(
      TMath::Sqrt(
        TMath::Power(v1Yerr.at(i)/v2Y.at(i),2) + 
        TMath::Power(v1Y.at(i)*v2Yerr.at(i)/(v2Y.at(i)*v2Y.at(i)),2)
      )
    );
  }
  TGraphErrors *grRatio = new TGraphErrors(v1X.size(),&v1X[0],&vRatioY[0],&v1Xerr[0],&vRatioYerr[0]);

  grRatio->GetXaxis()->SetLabelSize(0.11);
  grRatio->GetYaxis()->SetLabelSize(0.11);
  grRatio->GetXaxis()->SetTitleSize(0.12);
  grRatio->GetYaxis()->SetTitleSize(0.12);

  grRatio->GetYaxis()->SetTitle(Form("%s/%s",gr1->GetTitle(),gr2->GetTitle()));
  grRatio->GetYaxis()->SetTitleOffset(0.65);
  grRatio->GetXaxis()->SetTitle(Form("%s",gr1->GetXaxis()->GetTitle()));
  padDown->SetBottomMargin(0.3);
  grRatio->GetYaxis()->SetNdivisions(504);
  grRatio->GetXaxis()->SetTickLength(3*12/padUH);
  grRatio->GetYaxis()->SetTickLength(2.6*12/padUW);
  grRatio->GetYaxis()->SetRangeUser(yRatio_low,yRatio_high);
  // grRatio->GetYaxis()->SetRangeUser(0.81,1.19);

  grRatio->SetMarkerStyle(26);
  grRatio->SetMarkerSize(1.6);
  grRatio->SetLineColor(kRed);
  grRatio->SetMarkerColor(kRed);

  grRatio->GetXaxis()->SetLimits(0.95*vx_gr1[0],1.05*vx_gr1[n1bins-1]);

  grRatio->Draw("AP");

  TLine lineOne;
	lineOne.SetLineStyle(1);
	TLine line90;
	line90.SetLineWidth(2.);
	line90.SetLineStyle(2);	
	TLine line110;
	line110.SetLineWidth(2.);
	line110.SetLineStyle(2);

  // lineOne.SetLineColor(kRed);
  lineOne.DrawLine(0.95*vx_gr1[0],1.,  1.05*vx_gr1[n1bins-1],1.);
  line90.DrawLine( 0.95*vx_gr1[0],.95, 1.05*vx_gr1[n1bins-1],.95);
  line110.DrawLine(0.95*vx_gr1[0],1.05,1.05*vx_gr1[n1bins-1],1.05);

  return canv;
}

// Draws N TGraphErrors (upper panel) with their grN/gr1 ratio (lower pannel)
TCanvas *DrawTGraph(std::vector<TGraphErrors*> vgr, TString str, 
                    Double_t yRatio_low=0.89, Double_t yRatio_high=1.11,
                    Double_t x_low=0.0, Double_t x_high=1.0,
                    Double_t y_low=0.0, Double_t y_high=1.0,
                    Double_t leg_x_low=0.22, Double_t leg_y_low=0.55,
                    Double_t leg_x_high=0.55, Double_t leg_y_high=0.89)
{
  // Setting up global variables for the plot
  gROOT->SetStyle("Pub");
	gROOT->ForceStyle();
	gStyle->SetPalette(kDarkRainBow);
	gStyle->SetErrorX(0);

  std::vector<Double_t*> vx_gr, vy_gr, ex_gr, ey_gr;
  std::vector<Int_t> nbins;
  for (int i=0; i<vgr.size();i++)
  {
    // Read points
    vx_gr.push_back(vgr.at(i)->GetX());
    vy_gr.push_back(vgr.at(i)->GetY());

    // Read errors
    ex_gr.push_back(vgr.at(i)->GetEX());
    ey_gr.push_back(vgr.at(i)->GetEY());

    nbins.push_back(vgr.at(i)->GetN());
  }

  // Initialization of the canvas & pads
  TCanvas *canv = new TCanvas(Form("canv"),Form("Canvas"),900,800);
  if (vgr.size() < 2) return canv;
  canv->cd();
  TPad *padUp = new TPad(Form("padUp"),"v2 vs pt",0.,0.33,1.,1.,0,-1,0);
  TPad *padDown = new TPad(Form("padDown"),"Ratio v2",0.,0.,1.,0.33,0,-1,0);

  double padUW;
	double padUH;
	double padDW;
	double padDH;

  padUp->SetBorderSize(0);
  padDown->SetBorderSize(0);
  
  padUp->SetBottomMargin(0.);
  padDown->SetTopMargin(0.005);
  
  padUW = padUp->GetWw()*padUp->GetAbsWNDC();
  padUH = padUp->GetWh()*padUp->GetAbsHNDC();
  padDW = padDown->GetWw()*padDown->GetAbsWNDC();
  padDH = padDown->GetWh()*padDown->GetAbsHNDC();
  
  padUp->Draw();
  padDown->Draw();

  // Draw TGraphErrors in the upper pad
  padUp->cd();

  // gr1->GetXaxis()->SetLimits(0.95*vx_gr1[0],1.05*vx_gr1[n1bins-1]);
  vgr.at(0)->GetXaxis()->SetLimits(x_low,x_high);
  vgr.at(0)->GetYaxis()->SetRangeUser(y_low,y_high);

  vgr.at(0)->GetXaxis()->SetLabelSize(0.06);
  vgr.at(0)->GetYaxis()->SetLabelSize(0.06);
  vgr.at(0)->GetXaxis()->SetTitleSize(0.07);
  vgr.at(0)->GetYaxis()->SetTitleSize(0.07);
  vgr.at(0)->GetYaxis()->SetTitleOffset(1.08);

  vgr.at(0)->Draw("AP PLC PMC");
  for (int i=1; i<vgr.size();i++)
  {
    vgr.at(i)->Draw("P PLC PMC");
  }

  // TLegend *leg_pt = new TLegend(0.568,0.02,0.89,0.295);
  TLegend *leg_pt = new TLegend(leg_x_low,leg_y_low,leg_x_high,leg_y_high);
  leg_pt->SetBorderSize(0);
  leg_pt->SetHeader(str.Data(),"C");
  for (int i=0; i<vgr.size();i++)
  {
    leg_pt->AddEntry(vgr.at(i),Form("%s",vgr.at(i)->GetTitle()),"p");
  }

  leg_pt->Draw();

  //Draw grN/gr1 ratio in the bottom pad
  padDown->cd();

  std::vector<Double_t> v1X;
	std::vector<Double_t> v1Y;
	std::vector<Double_t> v1Xerr;
	std::vector<Double_t> v1Yerr;
  std::vector<Double_t> v2X;
	std::vector<Double_t> v2Y;
	std::vector<Double_t> v2Xerr;
	std::vector<Double_t> v2Yerr;
  std::vector<Double_t> vRatioY;
	std::vector<Double_t> vRatioYerr;

  std::vector<TGraphErrors*> vgrRatio;
  for (int igr=1; igr<vgr.size();igr++)
  {
    v1X.clear();
    v1Y.clear();
    v1Xerr.clear();
    v1Yerr.clear();
    v2X.clear();
    v2Y.clear();
    v2Xerr.clear();
    v2Yerr.clear();
    vRatioY.clear();
    vRatioYerr.clear();
    for (int i=0; i<vgr.at(igr)->GetN();i++)
    {
      v1X.push_back(vx_gr.at(igr)[i]);
      v1Y.push_back(abs(vy_gr.at(igr)[i]));
      v1Xerr.push_back(ex_gr.at(igr)[i]);
      v1Yerr.push_back(ey_gr.at(igr)[i]);

      v2Y.push_back((Double_t) abs(vgr.at(0)->Eval(v1X.at(i),0,"S")));
      v2Yerr.push_back(ey_gr.at(0)[i]);

      vRatioY.push_back(v1Y.at(i)/v2Y.at(i));
      vRatioYerr.push_back(
        TMath::Sqrt(
          TMath::Power(v1Yerr.at(i)/v2Y.at(i),2) + 
          TMath::Power(v1Y.at(i)*v2Yerr.at(i)/(v2Y.at(i)*v2Y.at(i)),2)
        )
      );
    }
    vgrRatio.push_back(new TGraphErrors(v1X.size(),&v1X[0],&vRatioY[0],&v1Xerr[0],&vRatioYerr[0]));
  }
  
  padDown->SetBottomMargin(0.3);

  for (int igr=0; igr<vgrRatio.size();igr++)
  {
    vgrRatio.at(igr)->GetXaxis()->SetLabelSize(0.11);
    vgrRatio.at(igr)->GetYaxis()->SetLabelSize(0.11);
    vgrRatio.at(igr)->GetXaxis()->SetTitleSize(0.12);
    vgrRatio.at(igr)->GetYaxis()->SetTitleSize(0.12);

    // vgrRatio.at(igr)->GetYaxis()->SetTitle(Form("%s/%s",vgr.at(igr+1)->GetTitle(),vgr.at(0)->GetTitle()));
    vgrRatio.at(igr)->GetYaxis()->SetTitle(Form("Ratio"));
    vgrRatio.at(igr)->GetYaxis()->SetTitleOffset(0.65);
    vgrRatio.at(igr)->GetXaxis()->SetTitle(Form("%s",vgr.at(0)->GetXaxis()->GetTitle()));
    vgrRatio.at(igr)->GetYaxis()->SetNdivisions(504);
    vgrRatio.at(igr)->GetXaxis()->SetTickLength(3*12/padUH);
    vgrRatio.at(igr)->GetYaxis()->SetTickLength(2.6*12/padUW);
    vgrRatio.at(igr)->GetYaxis()->SetRangeUser(yRatio_low,yRatio_high);

    vgrRatio.at(igr)->SetMarkerStyle(vgr.at(igr+1)->GetMarkerStyle());
    vgrRatio.at(igr)->SetMarkerSize(1.6);
    vgrRatio.at(igr)->SetLineColor(vgr.at(igr+1)->GetMarkerStyle());
    vgrRatio.at(igr)->SetMarkerColor(vgr.at(igr+1)->GetMarkerStyle());

    // grRatio->GetXaxis()->SetLimits(0.95*vx_gr1[0],1.05*vx_gr1[n1bins-1]);
    if (igr==0)
    {
      vgrRatio.at(igr)->GetXaxis()->SetLimits(x_low,x_high);
      vgrRatio.at(igr)->Draw("AP PLC PMC");
    }
    // else{
    vgrRatio.at(igr)->Draw("P PLC PMC");
    // }

    TLine lineOne;
    lineOne.SetLineStyle(1);
    lineOne.SetLineColor(1);
    TLine line90;
    line90.SetLineWidth(2.);
    line90.SetLineStyle(2);	
    TLine line110;
    line110.SetLineWidth(2.);
    line110.SetLineStyle(2);

    // lineOne.SetLineColor(kRed);
    lineOne.DrawLine(x_low,1.,  x_high,1.);
    line90.DrawLine( x_low,.95, x_high,.95);
    line110.DrawLine(x_low,1.05,x_high,1.05);
  }

  return canv;
}

// Draws N TGraphErrors (upper panel) with their grN/gr1 ratio (lower pannel)
TCanvas *DrawTGraphNoRatio(std::vector<TGraphErrors*> vgr, TString str,
                           Double_t x_low=0.0, Double_t x_high=1.0)
{
  // Setting up global variables for the plot
  gROOT->SetStyle("Pub");
	gROOT->ForceStyle();
	gStyle->SetPalette(kDarkRainBow);
	gStyle->SetErrorX(0);

  std::vector<Double_t*> vx_gr, vy_gr, ex_gr, ey_gr;
  std::vector<Int_t> nbins;
  for (int i=0; i<vgr.size();i++)
  {
    // Read points
    vx_gr.push_back(vgr.at(i)->GetX());
    vy_gr.push_back(vgr.at(i)->GetY());

    // Read errors
    ex_gr.push_back(vgr.at(i)->GetEX());
    ey_gr.push_back(vgr.at(i)->GetEY());

    nbins.push_back(vgr.at(i)->GetN());
  }

  // Initialization of the canvas & pads
  TCanvas *canv = new TCanvas(Form("canv"),Form("Canvas"),900,600);
  canv->cd();

  // gr1->GetXaxis()->SetLimits(0.95*vx_gr1[0],1.05*vx_gr1[n1bins-1]);
  vgr.at(0)->GetXaxis()->SetLimits(x_low,x_high);
  // vgr.at(0)->GetYaxis()->SetRangeUser(y_low,y_high);

  vgr.at(0)->GetXaxis()->SetLabelSize(0.06);
  vgr.at(0)->GetYaxis()->SetLabelSize(0.06);
  vgr.at(0)->GetXaxis()->SetTitleSize(0.07);
  vgr.at(0)->GetYaxis()->SetTitleSize(0.07);
  vgr.at(0)->GetYaxis()->SetTitleOffset(1.08);

  vgr.at(0)->Draw("AP PLC PMC");
  for (int i=1; i<vgr.size();i++)
  {
    vgr.at(i)->Draw("P PLC PMC");
  }

  // TLegend *leg_pt = new TLegend(0.568,0.02,0.89,0.295);
  TLegend *leg_pt = new TLegend(0.22,0.55,0.55,0.89);
  leg_pt->SetBorderSize(0);
  leg_pt->SetHeader(str.Data(),"C");
  for (int i=0; i<vgr.size();i++)
  {
    leg_pt->AddEntry(vgr.at(i),Form("%s",vgr.at(i)->GetTitle()),"p");
  }

  leg_pt->Draw();

  return canv;
}

// Save TGraphs in file
void SaveTGraph(TString outFileName, TGraphErrors *const &gr1, TGraphErrors *const &gr2)
{
  TFile *fo = new TFile(outFileName.Data(),"recreate");
  fo->cd();

  TCanvas *canv = (TCanvas*) DrawTGraph(gr1, gr2);

  gr1->Write();
  gr2->Write();
  canv->Write();
}

// Convert TProfile into TGraphErrors
TGraphErrors *ConvertProfile(TProfile *const& pr)
{
    std::vector<Double_t> vX, vY, vEx, vEy;
    for (int ibin=0; ibin<pr->GetNbinsX(); ibin++)
    {
        vX.push_back(pr->GetBinCenter(ibin+1));
        vEx.push_back(0.);
        vY.push_back(pr->GetBinContent(ibin+1));
        vEy.push_back(pr->GetBinError(ibin+1));
    }
    TGraphErrors *gr = new TGraphErrors(vX.size(),&vX[0],&vY[0],&vEx[0],&vEy[0]);
    gr->SetTitle(pr->GetTitle());
    gr->SetName(Form("%s_graph", pr->GetName()));
    gr->GetXaxis()->SetTitle(pr->GetXaxis()->GetTitle());
    gr->GetYaxis()->SetTitle(pr->GetYaxis()->GetTitle());
    gr->SetMarkerStyle(pr->GetMarkerStyle());
    gr->SetMarkerStyle(pr->GetMarkerStyle());
    gr->SetMarkerSize(pr->GetMarkerSize());
    gr->SetMarkerSize(pr->GetMarkerSize());
    gr->SetMarkerColor(pr->GetMarkerColor());
    gr->SetMarkerColor(pr->GetMarkerColor());
    gr->SetLineColor(pr->GetLineColor());
    gr->SetLineColor(pr->GetLineColor());

    return gr;
}

// Test function. Shows basic usage of DrawTGraph()
void Test()
{
  // cent 0-10%
  const std::vector<Double_t> vPt_cent1  = {0.247,0.348,0.448,0.548,0.648,0.748,0.848,0.948,1.092,1.292,1.492,1.692,1.892,2.200,2.703,3.343};
  const std::vector<Double_t> vV2_cent1  = {0.01342,0.01488,0.01231,0.02085,0.01557,0.02236,0.02656,0.03014,0.04275,0.03826,0.03859,0.04492,0.06318,0.06910,0.07798,0.07481};
  const std::vector<Double_t> ePt_cent1  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  const std::vector<Double_t> eV2_cent1  = {0.00114,0.00131,0.00133,0.00170,0.00166,0.00203,0.00233,0.00265,0.00268,0.00304,0.00367,0.00476,0.00654,0.00672,0.01123,0.01667};
  const std::vector<Double_t> sysl_cent1 = {0.00972,0.00130,0.00123,0.00187,0.00112,0.00129,0.00116,0.00121,0.00145,0.00110,0.00091,0.00120,0.00167,0.00168,0.00224,0.00294};
  const std::vector<Double_t> sysh_cent1 = {0.00972,0.00130,0.00123,0.00187,0.00112,0.00129,0.00116,0.00121,0.00145,0.00110,0.00091,0.00120,0.00167,0.00168,0.00224,0.00294};
  const std::vector<Double_t> eePt_cent1 = {0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05};

  // cent 10-20%
  const std::vector<Double_t> vPt_cent2  = {0.247,0.348,0.448,0.548,0.648,0.748,0.848,0.948,1.092,1.292,1.492,1.692,1.892,2.200,2.703,3.343};
  const std::vector<Double_t> vV2_cent2  = {0.02194,0.02987,0.03696,0.04342,0.05052,0.05556,0.06572,0.07064,0.07773,0.09169,0.10236,0.11847,0.13255,0.13748,0.15166,0.14679};
  const std::vector<Double_t> ePt_cent2  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  const std::vector<Double_t> eV2_cent2  = {0.00061,0.00070,0.00078,0.00088,0.00098,0.00110,0.00125,0.00141,0.00122,0.00155,0.00204,0.00275,0.00365,0.00363,0.00640,0.00945,0.02444,0.05846};
  const std::vector<Double_t> sysl_cent2 = {0.00506,0.00264,0.00283,0.00348,0.00269,0.00236,0.00224,0.00213,0.00207,0.00202,0.00227,0.00238,0.00308,0.00270,0.00313,0.00322,0.00505,0.00160};
  const std::vector<Double_t> sysh_cent2 = {0.00506,0.00264,0.00283,0.00348,0.00269,0.00236,0.00224,0.00213,0.00207,0.00202,0.00227,0.00238,0.00308,0.00270,0.00313,0.00322,0.00505,0.00160};
  const std::vector<Double_t> eePt_cent2 = {0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05};

  TGraphErrors *grPHENIX[2];
  grPHENIX[0] = new TGraphErrors(vPt_cent1.size(),&vPt_cent1[0],&vV2_cent1[0],&ePt_cent1[0],&eV2_cent1[0]);
  grPHENIX[1] = new TGraphErrors(vPt_cent2.size(),&vPt_cent2[0],&vV2_cent2[0],&ePt_cent2[0],&eV2_cent2[0]);

  // Make titile names for graph and axes.
  // Those names will be used in DrawTGraph() to plot legend and axis names
  grPHENIX[0]->SetTitle("0-10%");
  grPHENIX[0]->GetXaxis()->SetTitle("p_{T}, [GeV/c]");
  grPHENIX[0]->GetYaxis()->SetTitle("v_{2}");

  grPHENIX[1]->SetTitle("10-20%");
  grPHENIX[1]->GetXaxis()->SetTitle("p_{T}, [GeV/c]");
  grPHENIX[1]->GetYaxis()->SetTitle("v_{2}");

  TCanvas *canv = (TCanvas*) DrawTGraph(grPHENIX[0],grPHENIX[1]);
  canv->SetName("canv0");
  std::vector<TGraphErrors*> vgr;
  vgr.push_back(grPHENIX[0]);
  vgr.push_back(grPHENIX[1]);
  TCanvas *canv1 = (TCanvas*) DrawTGraph(vgr,"Ratio Plot Title",0.89, 1.11, 0.2, 3.5, 0., 0.15, 0.22, 0.55, 0.55, 0.89);
  canv1->SetName("canv1");
  SaveTGraph("outfile.root",grPHENIX[0],grPHENIX[1]);
}
