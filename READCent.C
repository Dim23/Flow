#include <TF1.h>
#include <TLine.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TLegend.h>
//#include "/home/dim2/GIT/DrawTGraph.C"

TGraphErrors *grv4[Nb], *grv2[Nb], *grvMc[Nb], *grvEp[Nb];
TGraphErrors *grIntv4[Nb], *grIntv2[Nb], *grIntvMc[Nb], *grIntvEp[Nb];
double RMSSV2[Nb], RMSSV4[Nb], RMSVMC[Nb], RMSVEP[Nb], binCent[Nb], IntMc[Nb], IntV2[Nb], IntV4[Nb], IntEv[Nb], ErMc[Nb], ErV2[Nb], ErV4[Nb], ErEv[Nb], RMSbinCent[Nb]; //RMS для SV2 и SV4

double sigmaX(TH1F *X, double weight)
{
    double RMS = X->GetStdDev();
    double neff = X->GetEffectiveEntries();
    return weight * RMS * RMS / (neff - 1);
}

double sigmaXY(TH1F *X, TH1F *Y, TH1F *XY, double weight)
{
    double RMS = XY->GetMean() - (X->GetMean()) * (Y->GetMean());
    double Xsum = X->GetSum();
    double Ysum = Y->GetSum();
    double XYsum = XY->GetSum();
    double neff = Xsum * Ysum / XYsum;
    return weight * RMS * RMS / (neff - 1);
}

void read(const char *infile, const char *savefile = "~/GIT/NoneFlow30_40.root")
{
    // Setting up global variables for the plot
    gROOT->SetStyle("Pub");
    gROOT->ForceStyle();

    double RMSdn2[Nb][NN];
    double RMSv2[Nb][NN];
    double RMSdn4[Nb][NN];
    double RMSv4[Nb][NN];
    double Rv2[Nb][NN];
    double Rv4[Nb][NN];
    double mc[Nb][NN];
    double RMSmc[Nb][NN];
    double RMSvep[Nb][NN];
    double RMSbinPt[Nb][NN];

    double dsv4[Nb][NN], dsv2[Nb][NN], sv4[Nb], sv2[Nb], v4[Nb], v2[Nb], vMC[Nb], vEP[Nb], cn2[Nb], cn4[Nb];

    TH1F *SV2[Nb]; //2 частичная референсная корреляция
    TH1F *SV4[Nb]; //4 частичная референсная корреляция
    TH1F *COSf1[Nb];
    TH1F *SINf1[Nb];
    TH1F *COSf1f2[Nb];
    TH1F *SINf1f2[Nb];
    TH1F *COSf1f2f3[Nb];
    TH1F *SINf1f2f3[Nb];

    TH1F *SV2_SV4[Nb];
    TH1F *HMC[Nb], *HRES[Nb], *HVobs[Nb], *HDiffRES[Nb][NN], *HDiffVobs[Nb][NN];
    TH1F *DiffSV2[Nb][NN];
    TH1F *DiffSV4[Nb][NN];
    TH1F *DiffMC[Nb][NN];
    TH1F *SV2_DiffSV2[Nb][NN];
    TH1F *SV2_DiffSV4[Nb][NN];
    TH1F *SV4_DiffSV2[Nb][NN];
    TH1F *SV4_DiffSV4[Nb][NN];
    TH1F *DiffSV2_DiffSV4[Nb][NN];
    TH1F *hBin_Pt[Nb][NN];
    TH1F *COSp1[Nb][NN];
    TH1F *SINp1[Nb][NN];
    TH1F *COSp1f2[Nb][NN];
    TH1F *SINp1f2[Nb][NN];
    TH1F *COSp1f2mf3[Nb][NN];
    TH1F *SINp1f2mf3[Nb][NN];
    TH1F *COSp1mf2mf3[Nb][NN];
    TH1F *SINp1mf2mf3[Nb][NN];

    double binPt[Nb][NN];
    double QMC[Nb][NN];
    double dn2[Nb][NN];
    double dn4[Nb][NN];
    double Diffv2[Nb][NN];
    double Diffv4[Nb][NN];
    double DiffvEP[Nb][NN];
    double Diffsv2[Nb][NN];
    double Diffsv4[Nb][NN];

    char str[200];

    TFile *f = new TFile(infile);
    for (int k = 0; k < Nb; k++)
    {
        sprintf(str, "%sCENT%d_%d", "sv2", k * 10, (k + 1) * 10);
        SV2[k] = (TH1F *)f->Get(str);

        sprintf(str, "%sCENT%d_%d", "sv4", k * 10, (k + 1) * 10);
        SV4[k] = (TH1F *)f->Get(str);

        sprintf(str, "%sCENT%d_%d", "sv2_sv4", k * 10, (k + 1) * 10);
        SV2_SV4[k] = (TH1F *)f->Get(str);

        sprintf(str, "%sCENT%d_%d", "mc", k * 10, (k + 1) * 10);
        HMC[k] = (TH1F *)f->Get(str);

        sprintf(str, "%sCENT%d_%d", "HRES", k * 10, (k + 1) * 10);
        HRES[k] = (TH1F *)f->Get(str);

        sprintf(str, "%sCENT%d_%d", "Vobs", k * 10, (k + 1) * 10);
        HVobs[k] = (TH1F *)f->Get(str);

        RMSbinCent[k] = 0;
        for (int m = 0; m < NN; m++)
        {
            sprintf(str, "%s%d_CENT%d_%d", "sv2Diff", m, k * 10, (k + 1) * 10);
            DiffSV2[k][m] = (TH1F *)f->Get(str);
            sprintf(str, "%s%d_CENT%d_%d", "sv4Diff", m, k * 10, (k + 1) * 10);
            DiffSV4[k][m] = (TH1F *)f->Get(str);
            sprintf(str, "%s%d_CENT%d_%d", "DiffMC", m, k * 10, (k + 1) * 10);
            DiffMC[k][m] = (TH1F *)f->Get(str);
            sprintf(str, "%s%d_CENT%d_%d", "SV2_DiffSV2", m, k * 10, (k + 1) * 10);
            SV2_DiffSV2[k][m] = (TH1F *)f->Get(str);
            sprintf(str, "%s%d_CENT%d_%d", "SV2_DiffSV4", m, k * 10, (k + 1) * 10);
            SV2_DiffSV4[k][m] = (TH1F *)f->Get(str);
            sprintf(str, "%s%d_CENT%d_%d", "SV4_DiffSV2", m, k * 10, (k + 1) * 10);
            SV4_DiffSV2[k][m] = (TH1F *)f->Get(str);
            sprintf(str, "%s%d_CENT%d_%d", "SV4_DiffSV4", m, k * 10, (k + 1) * 10);
            SV4_DiffSV4[k][m] = (TH1F *)f->Get(str);
            sprintf(str, "%s%d_CENT%d_%d", "DiffSV2_DiffSV4", m, k * 10, (k + 1) * 10);
            DiffSV2_DiffSV4[k][m] = (TH1F *)f->Get(str);
            sprintf(hBinPt, "%s%d_CENT%d_%d", "BIN_pt_", m, k * 10, (k + 1) * 10);
            sprintf(str, "%s%d_CENT%d_%d", "diffVobs", m, k * 10, (k + 1) * 10);
            HDiffVobs[k][m] = (TH1F *)f->Get(str);
            sprintf(str, "%s%d_CENT%d_%d", "diffRES", m, k * 10, (k + 1) * 10);
            HDiffRES[k][m] = (TH1F *)f->Get(str);
            //hBin_Pt[k][m] = 0.5 * (pt_bin[m] + pt_bin[m + 1]); //hBin_Pt[k][m]=(TH1F*)f->Get(HBinPt);
            RMSbinPt[k][m] = 0;
        }
    }

    //Референсный поток
    for (int k = 0; k < Nb; k++)
    {
        sv2[k] = SV2[k]->GetMean();
        v2[k] = pow(fabs(sv2[k]), 0.5);
        sv4[k] = SV4[k]->GetMean();
        cn4[k] = sv4[k] - 2 * (sv2[k] * sv2[k]);
        v4[k] = pow(fabs(cn4[k]), 0.25);
        vMC[k] = HMC[k]->GetMean();
        vEP[k] = HVobs[k]->GetMean() / pow(HRES[k]->GetMean(), 0.5);
        cout << " v2{2} is " << v2[k] << " v2{4} is " << v4[k] << " v2{MC} is " << vMC[k] << endl;
        binCent[k] = 5 + k * 10;
        //Дифференциальный поток
        cout << "Different. flow for cent " << 10 * (k + 0.5) << endl;
        for (int m = 0; m < NN; m++)
        {
            DiffvEP[k][m] = HDiffVobs[k][m]->GetMean() / pow(HDiffRES[k][m]->GetMean(), 0.5);
            dn2[k][m] = DiffSV2[k][m]->GetMean();
            Diffv2[k][m] = dn2[k][m] / pow(fabs(sv2[k]), 0.5);
            dsv4[k][m] = DiffSV4[k][m]->GetMean();
            dn4[k][m] = dsv4[k][m] - 2 * dn2[k][m] * sv2[k];
            Diffv4[k][m] = -dn4[k][m] / pow(fabs(cn4[k]), 0.75);
            binPt[k][m] = 0.5 * (pt_bin[m] + pt_bin[m + 1]);

            //cout <<"Different. flow "<<" v2{2} "<<Diffv2[k][m]<"+- "<<Diffv2[k][m]<<", v2{4} "<<Diffv4[k][m]<<" Pt "<<binPt[m]<<endl ;}
        }

        //погрешность для дифференциального потока
        for (int m = 0; m < NN; m++)
        {
            mc[k][m] = DiffMC[k][m]->GetMean();
            RMSmc[k][m] = pow(sigmaX(DiffMC[k][m], 1), 0.5);
            RMSvep[k][m] = sigmaX(HDiffVobs[k][m], 1) / HDiffRES[k][m]->GetMean() + 0.25 * HDiffVobs[k][m]->GetMean() * HDiffVobs[k][m]->GetMean() * sigmaX(HDiffRES[k][m], 1) / pow(HDiffRES[k][m]->GetMean(), 3);
            RMSvep[k][m] = pow(RMSvep[k][m], 0.5);
            RMSdn2[k][m] = sigmaX(DiffSV2[k][m], 1);
            RMSdn4[k][m] = DiffSV4[k][m]->GetStdDev();

            Rv2[k][m] = 0.25 * pow(sv2[k], -3) * (pow(dn2[k][m], 2) * sigmaX(SV2[k], 1) + 4 * sv2[k] * sv2[k] * sigmaX(DiffSV2[k][m], 1) - 4 * dn2[k][m] * sv2[k] * sigmaXY(SV2[k], DiffSV2[k][m], SV2_DiffSV2[k][m], 1));
            RMSv2[k][m] = pow(Rv2[k][m], 0.5);

            Rv4[k][m] = (pow(fabs(cn4[k]), -3.5)) * (pow((2 * sv2[k] * sv2[k] * dn2[k][m] - 3 * sv2[k] * dsv4[k][m] + 2 * sv4[k] * dn2[k][m]), 2) * sigmaX(SV2[k], 1) +
                                                     9 / 16 * pow(dn4[k][m], 2) * sigmaX(SV4[k], 1) + 4 * sv2[k] * sv2[k] * pow(fabs(cn4[k]), 2) * sigmaX(DiffSV2[k][m], 1) +

                                                     pow(fabs(cn4[k]), 2) * sigmaX(DiffSV4[k][m], 1) -

                                                     1.5 * fabs(dn4[k][m]) * (2 * sv2[k] * sv2[k] * dn2[k][m] - 3 * sv2[k] * dsv4[k][m] + 2 * sv4[k] * dn2[k][m]) * sigmaXY(SV2[k], SV4[k], SV2_SV4[k], 1) -

                                                     4 * sv2[k] * fabs(cn4[k]) * (2 * sv2[k] * sv2[k] * dn2[k][m] - 3 * sv2[k] * dsv4[k][m] + 2 * sv4[k] * dn2[k][m]) * sigmaXY(SV2[k], DiffSV2[k][m], SV2_DiffSV2[k][m], 1) +

                                                     2 * fabs(cn4[k]) * (2 * sv2[k] * sv2[k] * dn2[k][m] - 3 * sv2[k] * dsv4[k][m] + 2 * sv4[k] * dn2[k][m]) * sigmaXY(SV2[k], DiffSV4[k][m], SV2_DiffSV4[k][m], 1) +

                                                     3 * sv2[k] * (fabs(cn4[k])) * dn4[k][m] * sigmaXY(SV4[k], DiffSV2[k][m], SV4_DiffSV2[k][m], 1) -

                                                     .5 * fabs(cn4[k]) * sigmaXY(SV4[k], DiffSV4[k][m], SV4_DiffSV4[k][m], 1) -

                                                     4 * sv2[k] * pow(fabs(cn4[k]), 2) * sigmaXY(DiffSV2[k][m], DiffSV4[k][m], DiffSV2_DiffSV4[k][m], 1));

            RMSv4[k][m] = pow(Rv4[k][m], 0.5);
            cout << "Different. flow "
                 << " v2{2} " << Diffv2[k][m] << "+- " << RMSv2[k][m] << ", v2{4} " << Diffv4[k][m] << "+- " << RMSv4[k][m] << " Pt " << binPt[k][m] << endl;
        }
    }

    cout << "Ok " << endl;

    TFile *d_outfile = new TFile(savefile, "recreate");
    d_outfile->cd();


cout<<""<<endl;






    char name[400];
    for (int kk = 0; kk < Nb; kk++)
    {
        //Дифференциальный поток с погрешностями

        grv4[kk] = new TGraphErrors(NN, binPt[kk], Diffv4[kk], RMSbinPt[kk], RMSv4[kk]);
        sprintf(name, "gr_cent%i_2", kk);
        grv4[kk]->SetName(name);
        grv4[kk]->SetTitle("v_{2}{4}");
        grv4[kk]->GetYaxis()->SetRangeUser(-0.01, 0.26);
        grv4[kk]->GetXaxis()->SetRangeUser(0.1, 3.6);
        grv4[kk]->SetMarkerStyle(21);
        grv4[kk]->SetMarkerSize(1);
        grv4[kk]->SetMarkerColorAlpha(kBlack, 1);
        grv4[kk]->SetLineColorAlpha(kBlack, 1);
        grv4[kk]->SetLineWidth(1);
        grv4[kk]->GetXaxis()->SetTitle("Pt,Gev/c");
        grv4[kk]->GetYaxis()->SetTitle("V_{2}");
        grv4[kk]->Write();

        grv2[kk] = new TGraphErrors(NN, binPt[kk], Diffv2[kk], RMSbinPt[kk], RMSv2[kk]);
        sprintf(name, "gr_cent%i_1", kk);
        grv2[kk]->SetName(name);
        grv2[kk]->SetTitle("v_{2}{2}");
        grv2[kk]->SetMarkerStyle(20);
        grv2[kk]->SetMarkerSize(1);
        grv2[kk]->SetMarkerColorAlpha(kBlue, 1);
        grv2[kk]->SetLineColorAlpha(kBlue, 1);
        grv2[kk]->GetYaxis()->SetTitle("V_{2}");
        grv2[kk]->GetXaxis()->SetTitle("Pt,Gev/c");
        grv2[kk]->SetLineWidth(2);
        grv2[kk]->Write();

        grvMc[kk] = new TGraphErrors(NN, binPt[kk], mc[kk], RMSbinPt[kk], RMSmc[kk]);
        sprintf(name, "gr_cent%i_0", kk);
        grvMc[kk]->SetName(name);
        grvMc[kk]->SetTitle("v_{2}{MC}");
        grvMc[kk]->SetMarkerStyle(22);
        grvMc[kk]->SetMarkerSize(1);
        grvMc[kk]->SetMarkerColorAlpha(kGreen, 1);
        grvMc[kk]->SetLineColorAlpha(kGreen, 1);
        grvMc[kk]->SetLineWidth(2);
        grvMc[kk]->Write();

        grvEp[kk] = new TGraphErrors(NN, binPt[kk], DiffvEP[kk], RMSbinPt[kk], RMSvep[kk]);
        sprintf(name, "gr_cent%i_3", kk);
        grvEp[kk]->SetName(name);
        grvEp[kk]->SetTitle("v_{2}{EP}");
        grvEp[kk]->SetMarkerStyle(22);
        grvEp[kk]->SetMarkerSize(1);
        grvEp[kk]->SetMarkerColorAlpha(kRed, 1);
        grvEp[kk]->SetLineColorAlpha(kRed, 1);
        grvEp[kk]->SetLineWidth(2);
        grvEp[kk]->Write();

        //Рефренсный поток с погрешностями
        IntMc[kk] = vMC[kk];
        ErMc[kk] = pow(sigmaX(HMC[kk], 1), 0.5);
        IntV2[kk] = v2[kk];
        ErV2[kk] = pow(sigmaX(SV2[kk], 1), 0.5) / (2 * pow(sv2[kk], 0.5));
        IntV4[kk] = v4[kk];
        double Is4 = pow(fabs(cn4[kk]), -1.5) * (sigmaX(SV2[kk], sv2[kk] * sv2[kk]) + sigmaX(SV4[kk], 0.125) - 0.5 * sv2[kk] * sigmaXY(SV2[kk], SV4[kk], SV2_SV4[kk], 1));
        ErV4[kk] = pow(Is4, 0.5);
        IntEv[kk] = vEP[kk];
        ErEv[kk] = pow(sigmaX(HVobs[kk], 1) / HRES[kk]->GetMean() + 0.25 * HVobs[kk]->GetMean() * HVobs[kk]->GetMean() * sigmaX(HRES[kk], 1) / pow(HRES[kk]->GetMean(), 3), 0.5);

        double x[1] = {0.5}, y[1] = {vMC[kk]}, ex[1] = {0.}, ey[1] = {ErMc[kk]};
        grIntvMc[kk] = new TGraphErrors(1, x, y, ex, ey);
        grIntvMc[kk]->SetName("vMC");
        grIntvMc[kk]->SetMarkerStyle(20);
        grIntvMc[kk]->GetYaxis()->SetTitle("V_{2}");

        double x2[1] = {1.5}, y2[1] = {v2[kk]}, ex2[1] = {0.}, ey2[1] = {ErV2[kk]};
        grIntv2[kk] = new TGraphErrors(1, x2, y2, ex2, ey2);
        grIntv2[kk]->SetName("v2");
        grIntv2[kk]->SetMarkerStyle(20);

        double x4[1] = {2.5}, y4[1] = {v4[kk]}, ex4[1] = {0.}, ey4[1] = {ErV4[kk]};
        grIntv4[kk] = new TGraphErrors(1, x4, y4, ex4, ey4);
        grIntv4[kk]->SetName("v4");
        grIntv4[kk]->SetMarkerStyle(20);

        double x5[1] = {3.5}, y5[1] = {vEP[kk]}, ex5[1] = {0.}, ey5[1] = {ErEv[kk]};
        grIntvEp[kk] = new TGraphErrors(1, x5, y5, ex5, ey5);
        grIntvEp[kk]->SetName("vEP");
        grIntvEp[kk]->SetMarkerStyle(20);
        grIntvEp[kk]->GetYaxis()->SetTitle("V_{2}");
    }
    d_outfile->Close();

    Int_t cent = 3;

    TCanvas *cdif = new TCanvas("cdif", "def Flow All", 200, 10, 1600, 900);
    cdif->Divide(3, 2, 0, 0);

    TLegend *legdif[Nb];
    for (int k = 0; k < 6; k++)
    {
        cdif->cd(k + 1);

        if (k % 3 == 0)
        {
            legdif[k] = new TLegend(0.16, 0.72, 0.35, 0.995);
        }
        else
        {
            legdif[k] = new TLegend(0., 0.72, 0.24, 0.995);
        }
        legdif[k]->AddEntry(grv4[k], "V_{2}{4}", "pe");
        legdif[k]->AddEntry(grv2[k], "V_{2}{2}", "pe");
        legdif[k]->AddEntry(grvMc[k], "V_{2}{MCk}", "pe");
        legdif[k]->AddEntry(grvEp[k], "V_{2}{EP}", "pe");
        char strleg[200];
        sprintf(strleg, " cent: %i-%i%%", k * 10, (k + 1) * 10);

        legdif[k]->SetHeader(strleg);
        grv4[k]->Draw("AP");
        grv2[k]->Draw("SAME P");
        grvMc[k]->Draw("SAME P");
        grvEp[k]->Draw("SAME P");
        legdif[k]->Draw();
    }

    TCanvas *c2 = new TCanvas("c2", "Ref Flow Cent", 0, 0, 600, 400);
    c2->SetLeftMargin(0.1);
    c2->SetRightMargin(0.1);
    c2->SetBottomMargin(0.1);
    c2->SetTopMargin(0.1);
    TH1F *h = new TH1F("h", "", 4, 0., 4.);
    const Int_t nx = 4;
    const char *month[nx] = {"V_{2}{MC}", "V_{2}{2}", "V_{2}{4}", "V_{2}{EP}"};
    h->SetLineColorAlpha(kBlack, 1);
    h->SetCanExtend(TH1::kAllAxes);
    h->SetStats(0);
    for (Int_t i = 0; i < 4; i++)
    {
        h->Fill(month[i], 1);
    }
    h->LabelsDeflate("X");
    double ey4[1] = {ErV4[cent]};
    //ey4[0]= grIntv4[kk]->GetErrorY(0);
    double ymin = min(vMC[cent] - ey4[0], v4[cent] - ey4[0]);
    double ymax = max(v4[cent] + ey4[0], v2[cent] + ey4[0]);
    h->GetYaxis()->SetRangeUser(0.95 * ymin, 1.05 * ymax);
    h->Draw();

    TLine *line = new TLine(0., vMC[cent], 4., vMC[cent]);
    line->SetLineWidth(1);
    line->SetLineStyle(5);

    TLegend *leg2 = new TLegend(0.1, 0.6, 0.3, 0.9);
    leg2->AddEntry(grIntv4[cent], "V_{2}{4}", "pe");
    leg2->AddEntry(grIntv2[cent], "V_{2}{2}", "pe");
    leg2->AddEntry(grIntvMc[cent], "V_{2}{MC}", "pe");
    leg2->AddEntry(grIntvEp[cent], "V_{2}{EP}", "pe");
    leg2->SetTextSize(0.04);
    char strleg[200];
    sprintf(strleg, " cent: %i-%i%%", cent * 10, (cent + 1) * 10);
    leg2->SetHeader(strleg);
    grIntv4[cent]->Draw("SAME P");
    grIntv2[cent]->Draw("SAME P");
    grIntvMc[cent]->Draw("SAME P");
    grIntvEp[cent]->Draw("SAME P");
    line->Draw("SAME");
    leg2->Draw();
    cout << "Q2Er " << ErV2[cent] << "Q4Er " << ErV4[cent] << "McEr " << ErMc[cent] << endl;

    //Референсный поток по центральсности
    TCanvas *c3 = new TCanvas("c3", "Ref Flow be Centrality", 10, 10, 600, 400);
    TGraphErrors *grIv4 = new TGraphErrors(Nb, binCent, IntV4, RMSbinCent, ErV4);
    grIv4->SetMarkerStyle(21);
    grIv4->SetMarkerSize(1);
    grIv4->SetMarkerColorAlpha(kBlack, 1);
    grIv4->SetLineColorAlpha(kBlack, 1);
    grIv4->SetLineWidth(2);
    grIv4->GetXaxis()->SetTitle("Cent, %");
    grIv4->GetYaxis()->SetTitle("V_{2}");

    TGraphErrors *grIv2 = new TGraphErrors(Nb, binCent, IntV2, RMSbinCent, ErV2);
    grIv2->SetMarkerStyle(20);
    grIv2->SetMarkerSize(1);
    grIv2->SetMarkerColorAlpha(kBlue, 1);
    grIv2->SetLineColorAlpha(kBlue, 1);
    grIv2->SetLineWidth(2);

    TGraphErrors *grIvMc = new TGraphErrors(Nb, binCent, IntMc, RMSbinCent, ErMc);
    grIvMc->GetYaxis()->SetRangeUser(0, 0.15);
    grIvMc->SetMarkerStyle(22);
    grIvMc->SetMarkerSize(1);
    grIvMc->SetMarkerColorAlpha(kGreen, 1);
    grIvMc->SetLineColorAlpha(kGreen, 1);
    grIvMc->SetLineWidth(1);
    grIvMc->GetYaxis()->SetTitle("V_{2}");

    TGraphErrors *grIvEp = new TGraphErrors(Nb, binCent, IntEv, RMSbinCent, ErEv);
    grIvEp->SetName("diff_vEVENTplane");
    grIvEp->SetMarkerStyle(22);
    grIvEp->SetMarkerSize(1);
    grIvEp->SetMarkerColorAlpha(kRed, 1);
    grIvEp->SetLineColorAlpha(kRed, 1);
    grIvEp->SetLineWidth(1);

    TLegend *leg3 = new TLegend(0.15, 0.65, 0.3, 0.9);
    leg3->AddEntry(grIv4, "V_{2}{4}", "pe");
    leg3->AddEntry(grIv2, "V_{2}{2}", "pe");
    leg3->AddEntry(grIvMc, "V_{2}{MC}", "pe");
    leg3->AddEntry(grIvEp, "V_{2}{EP}", "pe");
    leg3->SetTextSize(0.04);
    grIv4->Draw("AP");
    grIv2->Draw("SAME P");
    grIvMc->Draw("SAME P");
    grIvEp->Draw("SAME P");
    leg3->Draw();

    //Дифференциальный поток полученный всеми способоами на одном канвасе
    TCanvas *c1 = new TCanvas("c1", "Def Flow for cent", 0, 0, 600, 400);
    grv4[cent]->Draw("AP");
    grv2[cent]->Draw("SAME P");
    grvMc[cent]->Draw("SAME P");
    grvEp[cent]->Draw("SAME P");
    TLegend *dleg = new TLegend(0.15, 0.6, 0.32, 0.9);
    dleg->AddEntry(grv4[cent], "V_{2}{4}", "pe");
    dleg->AddEntry(grv2[cent], "V_{2}{2}", "pe");
    dleg->AddEntry(grvMc[cent], "V_{2}{MC}", "pe");
    dleg->AddEntry(grvEp[cent], "V_{2}{EP}", "pe");
    dleg->SetTextSize(0.04);
    char strlegd[200];
    sprintf(strlegd, "cent: %i-%i%%", cent * 10, (cent + 1) * 10);
    dleg->SetHeader(strlegd);
    dleg->Draw();

    //Ratio дифференциального потока для центральностей cent*10-(cent+1)*10
    //DrawTGraph(grvEp[cent], grvMc[cent]);
}
