
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
    leg_pt->AddEntry(vgr.at(i),Form("%s",vgr.at(i)->GetTitle()),"pe");
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




void AcceptanceRatio(){
  TFile *inputFile = new TFile("~/GIT/NoneFlow30_40.root","read");
  TGraphErrors *gr[4][8];
  char name[400];
  for (int icent=0; icent<8; icent++){
    for (int i=0; i<4; i++){
      sprintf(name,"gr_cent%i_%i",icent,i);
      gr[i][icent] = (TGraphErrors*)inputFile->Get(name);
    }
  }
  std::vector<TGraphErrors*> vgr[8];
  for (int icent=0; icent<8; icent++){
    for (int i=0; i<4; i++){
      vgr[icent].push_back(gr[i][icent]);
    }  
  }
  TCanvas *can[8];
  TLatex l[8];
  for (int icent=0; icent<8; icent++){
    //                                                    yRatio_low    x_low     y_low    leg_x_low  leg_x_high
    can[icent] = (TCanvas*) DrawTGraph(vgr[icent],"v2 ratio",0.89, 1.11, -0.005, 3.5, -0.01, 0.25, 0.65, 0.11, 0.89, 0.35);
    //                                                          yRatio_high  x_high   y_high     leg_y_low   leg_y_high
    sprintf(name,"Cent%i-%i%%",icent*10,(icent+1)*10);
    can[icent] -> SetName(name);
    l[icent].SetNDC();
    l[icent].SetTextSize(0.15);
    l[icent].SetTextAlign(21);  
    l[icent].DrawLatex(0.5,0.1,name);
    sprintf(name,"~/GIT/Graphics/Cent%i-%i%%.png",icent*10,(icent+1)*10);
    can[icent] -> SaveAs(name);
  }
}
