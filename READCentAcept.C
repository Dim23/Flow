#include <TF1.h>
#include <TLine.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>

TGraphErrors *grv4[Nb], *grv2[Nb],*grv4IMP[Nb], *grv2IMP[Nb], *grvMc[Nb], *grvEp[Nb];
TGraphErrors *grIntv4[Nb], *grIntv2[Nb],*grIntv4IMP[Nb], *grIntv2IMP[Nb], *grIntvMc[Nb], *grIntvEp[Nb];

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

void read(const char *infile, const char *savefile = "~/FLOW5/PLOT/Acept30_40.root")
{
double b_bin[9] = {0.0, 4.18, 6.01, 7.37, 8.52, 9.57, 10.55, 11.46, 12.31};
static const int Nb = 8;
double pt_bin[13] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 2.6, 3.0, 3.5};
static const int NN = 12;

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

double IMPdn2[Nb][NN];double IMPdv2[Nb][NN];
double IMPdn4[Nb][NN];double IMPdv4[Nb][NN];
double IMPsv4[Nb],IMPsv2[Nb],IMPv4[Nb],IMPv2[Nb],IMPcn2[Nb],IMPcn4[Nb];
double dCOSp1[Nb][NN];double dSINp1[Nb][NN];
double dCOSp1f2[Nb][NN];double dSINp1f2[Nb][NN];
double dCOSp1f2mf3[Nb][NN];double dSINp1f2mf3[Nb][NN];
double dCOSp1mf2mf3[Nb][NN];double dSINp1mf2mf3[Nb][NN];
double dCOSf1[Nb],dCOSf1f2[Nb],dCOSf1mf2[Nb],dCOSf1f2f3[Nb],dSINf1[Nb],dSINf1f2[Nb],dSINf1f2f3[Nb];

double binPt[Nb][NN];
double QMC[Nb][NN];
double dn2[Nb][NN];
double dn4[Nb][NN];
double Diffv2[Nb][NN];
double Diffv4[Nb][NN];
double DiffvEP[Nb][NN];
double Diffsv2[Nb][NN];
double Diffsv4[Nb][NN];

char strvObs[200];
char strRES[200];
char strCOVsv2sv2[200];
char strCOVsv2sv4[200];
char strCOVsv4sv2[200];
char strCOVsv4sv4[200];
char DiffCOVsv2sv4[200];
char strE[200];
char strW[200];
char strMC[200];

char strCOSp1[200];
char strSINp1[200];
char strCOSp1f2[200];
char strSINp1f2[200];
char strCOSp1f2mf3[200];
char strSINp1f2mf3[200];
char strCOSp1mf2mf3[200];
char strSINp1mf2mf3[200];
char hBinPt[200];
char str[200];
const char *STR;
    TFile *f = new TFile(infile);
    for (int k = 0; k < Nb; k++)
    {
        sprintf(str, "%sCENT%d_%d", "sv2", k * 10, (k + 1) * 10);
        STR = (char *)str;
        SV2[k] = (TH1F *)f->Get(STR);

        sprintf(str, "%sCENT%d_%d", "sv4", k * 10, (k + 1) * 10);
        STR = (char *)str;
        SV4[k] = (TH1F *)f->Get(STR);

        sprintf(str, "%sCENT%d_%d", "sv2_sv4", k * 10, (k + 1) * 10);
        STR = (char *)str;
        SV2_SV4[k] = (TH1F *)f->Get(STR);

        sprintf(str, "%sCENT%d_%d", "mc", k * 10, (k + 1) * 10);
        STR = (char *)str;
        HMC[k] = (TH1F *)f->Get(STR);

        sprintf(str, "%sCENT%d_%d", "HRES", k * 10, (k + 1) * 10);
        STR = (char *)str;
        HRES[k] = (TH1F *)f->Get(STR);

        sprintf(str, "%sCENT%d_%d", "Vobs", k * 10, (k + 1) * 10);
        STR = (char *)str;
        HVobs[k] = (TH1F *)f->Get(STR);

        sprintf(str, "%sCENT%d_%d", "cosf1", k * 10, (k + 1) * 10);
        STR = (char *)str;
        COSf1[k] = (TH1F *)f->Get(STR);

        sprintf(str, "%sCENT%d_%d", "cosf1f2", k * 10, (k + 1) * 10);
        STR = (char *)str;
        COSf1f2[k] = (TH1F *)f->Get(STR);

        sprintf(str, "%sCENT%d_%d", "cosf1f2f3", k * 10, (k + 1) * 10);
        STR = (char *)str;
        COSf1f2f3[k] = (TH1F *)f->Get(STR);

        sprintf(str, "%sCENT%d_%d", "sinf1", k * 10, (k + 1) * 10);
        STR = (char *)str;
        SINf1[k] = (TH1F *)f->Get(STR);

        sprintf(str, "%sCENT%d_%d", "sinf1f2", k * 10, (k + 1) * 10);
        STR = (char *)str;
        SINf1f2[k] = (TH1F *)f->Get(STR);

        sprintf(str, "%sCENT%d_%d", "sinf1f2f3", k * 10, (k + 1) * 10);
        STR = (char *)str;
        SINf1f2f3[k] = (TH1F *)f->Get(STR);

        RMSbinCent[k] = 0;
        for (int m = 0; m < NN; m++)
        {
            sprintf(strE, "%s%d_CENT%d_%d", "sv2Diff", m, k * 10, (k + 1) * 10);
            sprintf(strW, "%s%d_CENT%d_%d", "sv4Diff", m, k * 10, (k + 1) * 10);
            sprintf(strMC, "%s%d_CENT%d_%d", "DiffMC", m, k * 10, (k + 1) * 10);
            sprintf(strCOVsv2sv2, "%s%d_CENT%d_%d", "SV2_DiffSV2", m, k * 10, (k + 1) * 10);
            sprintf(strCOVsv2sv4, "%s%d_CENT%d_%d", "SV2_DiffSV4", m, k * 10, (k + 1) * 10);
            sprintf(strCOVsv4sv2, "%s%d_CENT%d_%d", "SV4_DiffSV2", m, k * 10, (k + 1) * 10);
            sprintf(strCOVsv4sv4, "%s%d_CENT%d_%d", "SV4_DiffSV4", m, k * 10, (k + 1) * 10);
            sprintf(DiffCOVsv2sv4, "%s%d_CENT%d_%d", "DiffSV2_DiffSV4", m, k * 10, (k + 1) * 10);

            sprintf(strCOSp1, "%s%d_CENT%d_%d", "COSp1", m, k * 10, (k + 1) * 10);
            sprintf(strCOSp1f2, "%s%d_CENT%d_%d", "COSp1f2", m, k * 10, (k + 1) * 10);
            sprintf(strCOSp1f2mf3, "%s%d_CENT%d_%d", "COSp1f2mf3", m, k * 10, (k + 1) * 10);
            sprintf(strCOSp1mf2mf3, "%s%d_CENT%d_%d", "COSp1mf2mf3", m, k * 10, (k + 1) * 10);
            sprintf(strSINp1, "%s%d_CENT%d_%d", "SINp1", m, k * 10, (k + 1) * 10);
            sprintf(strSINp1f2, "%s%d_CENT%d_%d", "SINp1f2", m, k * 10, (k + 1) * 10);
            sprintf(strSINp1f2mf3, "%s%d_CENT%d_%d", "SINp1f2mf3", m, k * 10, (k + 1) * 10);
            sprintf(strSINp1mf2mf3, "%s%d_CENT%d_%d", "SINp1mf2mf3", m, k * 10, (k + 1) * 10);
            sprintf(hBinPt, "%s%d_CENT%d_%d", "BIN_pt_", m, k * 10, (k + 1) * 10);
            sprintf(strvObs, "%s%d_CENT%d_%d", "diffVobs", m, k * 10, (k + 1) * 10);
            sprintf(strRES, "%s%d_CENT%d_%d", "diffRES", m, k * 10, (k + 1) * 10);

            const char *rvObs = (char *)strvObs;
            const char *rRES = (char *)strRES;
            const char *EastH = (char *)strE;
            const char *WastH = (char *)strW;
            const char *SMC = (char *)strMC;
            const char *stCOVsv2sv2 = (char *)strCOVsv2sv2;
            const char *stCOVsv2sv4 = (char *)strCOVsv2sv4;
            const char *stCOVsv4sv2 = (char *)strCOVsv4sv2;
            const char *stCOVsv4sv4 = (char *)strCOVsv4sv4;
            const char *DifCOVsv2sv4 = (char *)DiffCOVsv2sv4;
const char *stCOSp1 = (char *)strCOSp1;
            const char *stCOSp1f2 = (char *)strCOSp1f2;
            const char *stCOSp1f2mf3 = (char *)strCOSp1f2mf3;
            const char *stCOSp1mf2mf3 = (char *)strCOSp1mf2mf3;
            const char *stSINp1 = (char *)strSINp1;
            const char *stSINp1f2 = (char *)strSINp1f2;
            const char *stSINp1f2mf3 = (char *)strSINp1f2mf3;
            const char *stSINp1mf2mf3 = (char *)strSINp1mf2mf3;
            const char *HBinPt = (char *)hBinPt;
            HDiffVobs[k][m] = (TH1F *)f->Get(rvObs);
            HDiffRES[k][m] = (TH1F *)f->Get(rRES);
            DiffSV2[k][m] = (TH1F *)f->Get(EastH);
            DiffSV4[k][m] = (TH1F *)f->Get(WastH);
            DiffMC[k][m] = (TH1F *)f->Get(SMC);
            SV2_DiffSV2[k][m] = (TH1F *)f->Get(stCOVsv2sv2);
            SV2_DiffSV4[k][m] = (TH1F *)f->Get(stCOVsv2sv4);
            SV4_DiffSV2[k][m] = (TH1F *)f->Get(stCOVsv4sv2);
            SV4_DiffSV4[k][m] = (TH1F *)f->Get(stCOVsv4sv4);
            DiffSV2_DiffSV4[k][m] = (TH1F *)f->Get(DifCOVsv2sv4);


COSp1[k][m]=(TH1F*)f->Get(stCOSp1);
COSp1f2[k][m]=(TH1F*)f->Get(stCOSp1f2);
COSp1f2mf3[k][m]=(TH1F*)f->Get(stCOSp1f2mf3);
COSp1mf2mf3[k][m]=(TH1F*)f->Get(stCOSp1mf2mf3);
SINp1[k][m]=(TH1F*)f->Get(stSINp1);
SINp1f2[k][m]=(TH1F*)f->Get(stSINp1f2);
SINp1f2mf3[k][m]=(TH1F*)f->Get(stSINp1f2mf3);
SINp1mf2mf3[k][m]=(TH1F*)f->Get(stSINp1mf2mf3);
//cout <<"Mean  "<<m<<endl;
dCOSp1[k][m]=COSp1[k][m]->GetMean();
dCOSp1f2[k][m]=COSp1f2[k][m]->GetMean();
dCOSp1f2mf3[k][m]=COSp1f2mf3[k][m]->GetMean();
dCOSp1mf2mf3[k][m]=COSp1mf2mf3[k][m]->GetMean();
dSINp1[k][m]=SINp1[k][m]->GetMean();
dSINp1f2[k][m]=SINp1f2[k][m]->GetMean();
dSINp1f2mf3[k][m]=SINp1f2mf3[k][m]->GetMean();
dSINp1mf2mf3[k][m]=SINp1mf2mf3[k][m]->GetMean();

            //hBin_Pt[k][m] = 0.5 * (pt_bin[m] + pt_bin[m + 1]); //hBin_Pt[k][m]=(TH1F*)f->Get(HBinPt);
            RMSbinPt[k][m] = 0;
        }
    }

    //Референсный поток
    for (int k = 0; k < Nb; k++)
    {

dCOSf1[k]=COSf1[k]->GetMean();
cout <<"Start read  "<<endl;
dCOSf1f2[k]=COSf1f2[k]->GetMean();
cout <<"Start read  "<<endl;
dCOSf1f2f3[k]=COSf1f2f3[k]->GetMean();
dSINf1[k]=SINf1[k]->GetMean();
dSINf1f2[k]=SINf1f2[k]->GetMean();
dSINf1f2f3[k]=SINf1f2f3[k]->GetMean();

        sv2[k] = SV2[k]->GetMean();v2[k] = pow(fabs(sv2[k]), 0.5);
        IMPsv2[k]=sv2[k]-(pow(COSf1[k]->GetMean(),2)+pow(SINf1[k]->GetMean(),2));IMPv2[k]=pow(fabs(IMPsv2[k]),0.5);

        sv4[k] = SV4[k]->GetMean();
        cn4[k] = sv4[k] - 2 * (sv2[k] * sv2[k]);
        v4[k] = pow(fabs(cn4[k]), 0.25);

IMPcn4[k]=sv4[k]-2*(sv2[k]*sv2[k])-4*dCOSf1[k]*dCOSf1f2f3[k]+4*dSINf1[k]*dSINf1f2f3[k]-pow(dCOSf1f2[k],2)-pow(dSINf1f2[k],2)+
4*dCOSf1f2[k]*(pow(dCOSf1[k],2)-pow(dSINf1[k],2))+8*dSINf1f2[k]*dSINf1[k]*dCOSf1[k]+
8*sv2[k]*(pow(dCOSf1[k],2)+pow(dSINf1[k],2))-6*pow((pow(dCOSf1[k],2)+pow(dSINf1[k],2)),2);
IMPv4[k]=pow(fabs(IMPcn4[k]),0.25);

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

IMPdn2[k][m]=DiffSV2[k][m]->GetMean()-dCOSp1[k][m]*dCOSf1[k]-dSINp1[k][m]*dSINf1[k];IMPdv2[k][m]=IMPdn2[k][m]/pow(fabs(IMPsv2[k]),0.5);

            dsv4[k][m] = DiffSV4[k][m]->GetMean();
            dn4[k][m] = dsv4[k][m] - 2 * dn2[k][m] * sv2[k];
            Diffv4[k][m] = -dn4[k][m] / pow(fabs(cn4[k]), 0.75);

IMPdn4[k][m]=dsv4[k][m]-2*dn2[k][m]*sv2[k]-dCOSp1[k][m]*dCOSf1f2f3[k]+dSINp1[k][m]*dSINf1f2f3[k]-dCOSf1[k]*dCOSp1mf2mf3[k][m]+dSINf1[k]*dSINp1mf2mf3[k][m]-2*dCOSf1[k]*dCOSp1f2mf3[k][m]-
2*dSINf1[k]*dSINp1f2mf3[k][m]-dCOSp1f2[k][m]*dCOSf1f2[k]-dSINp1f2[k][m]*dSINf1f2[k]+2*dCOSf1f2[k]*(dCOSp1[k][m]*dCOSf1[k]-dSINp1[k][m]*dSINf1[k])+2*dSINf1f2[k]*(dCOSp1[k][m]*dSINf1[k]+dSINp1[k][m]*dCOSf1[k])+
4*sv2[k]*(dCOSp1[k][m]*dCOSf1[k]+dSINp1[k][m]*dSINf1[k])+2*dCOSp1f2[k][m]*(pow(dCOSf1[k],2)-pow(dSINf1[k],2))+4*dSINp1f2[k][m]*dCOSf1[k]*dSINf1[k]+4*dn2[k][m]*(pow(dCOSf1[k],2)+pow(dSINf1[k],2))-
6*(pow(dCOSf1[k],2)-pow(dSINf1[k],2))*(dCOSp1[k][m]*dCOSf1[k]-dSINp1[k][m]*dSINf1[k])-12*dCOSf1[k]*dSINf1[k]*(dSINp1[k][m]*dCOSf1[k]+dCOSp1[k][m]*dSINf1[k]);

IMPdv4[k][m]=-IMPdn4[k][m]/pow(fabs(IMPcn4[k]),0.75);

            binPt[k][m] = 0.5 * (pt_bin[m] + pt_bin[m + 1]);

            //cout <<"Different. flow "<<" v2{2} "<<Diffv2[k][m]<"+- "<<Diffv2[k][m]<<", v2{4} "<<Diffv4[k][m]<<" Pt "<<binPt[m]<<endl ;}
        }

        //погрешность для дифференциального потока
        for (int m = 0; m < NN; m++)
        {
            mc[k][m] = DiffMC[k][m]->GetMean();
            RMSmc[k][m] = pow(sigmaX(DiffMC[k][m], 1), 0.5);
            RMSvep[k][m] = sigmaX(HDiffVobs[k][m], 1)/ HDiffRES[k][m]->GetMean()+0.25*HDiffVobs[k][m]->GetMean()*HDiffVobs[k][m]->GetMean()*sigmaX(HDiffRES[k][m], 1)/pow(HDiffRES[k][m]->GetMean(),3) ;
            RMSvep[k][m]=pow(RMSvep[k][m],0.5);
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

    /*TFile *d_outfile = new TFile(savefile, "recreate");
    d_outfile->cd();

    grv4->Write("grV4");
    grv2->Write("grV2");
    grv->Write("grMc");
    grvE->Write("grEP");
    grMC->Write("grIntMC");
    grEP->Write("grIntEP");
    grsv4->Write("grIntV2");
    grsv2->Write("grIntV4");
    d_outfile->Close();*/

double IntIMPv4[Nb],IntIMPv2[Nb];

    for (int kk = 0; kk < Nb; kk++)
    {
        //Рефренсный поток с погрешностями
        IntMc[kk] = vMC[kk]; ErMc[kk] = pow(sigmaX(HMC[kk], 1), 0.5);
        IntV2[kk] = v2[kk];ErV2[kk] = pow(sigmaX(SV2[kk], 1), 0.5)/(2*pow(sv2[kk],0.5));
IntIMPv2[kk] = IMPv2[kk];
        IntV4[kk] = v4[kk];
        double Is4 = pow(fabs(cn4[kk]), -1.5) * (sigmaX(SV2[kk], sv2[kk] * sv2[kk]) + sigmaX(SV4[kk], 0.125) - 0.5 * sv2[kk] * sigmaXY(SV2[kk], SV4[kk], SV2_SV4[kk], 1));
IntIMPv4[kk] = IMPv4[kk];
        ErV4[kk] = pow(Is4, 0.5);
        IntEv[kk] = vEP[kk]; 
ErEv[kk] = pow(sigmaX(HVobs[kk],1)/HRES[kk]->GetMean()+0.25*HVobs[kk]->GetMean()*HVobs[kk]->GetMean()*sigmaX(HRES[kk], 1)/pow(HRES[kk]->GetMean(),3), 0.5)      ;

    //Дифференциальный поток с погрешностями
    grv4[kk] = new TGraphErrors(NN, binPt[kk], Diffv4[kk], RMSbinPt[kk], RMSv4[kk]);
    grv4[kk]->SetName("diff_v4");
    grv4[kk]->GetYaxis()->SetRangeUser(0, 0.25);
    grv4[kk]->SetMarkerStyle(21);
    grv4[kk]->SetMarkerSize(1);
    grv4[kk]->SetMarkerColorAlpha(kBlack, 1);
    grv4[kk]->SetLineColorAlpha(kBlack, 1);
    grv4[kk]->SetLineWidth(1);
    grv4[kk]->GetXaxis()->SetTitle("Pt,Gev/c");
    grv4[kk]->GetYaxis()->SetTitle("V_{2}");
    grv4[kk]->SetTitle("");

    grv2[kk] = new TGraphErrors(NN, binPt[kk], Diffv2[kk], RMSbinPt[kk], RMSv2[kk]);
    grv2[kk]->SetName("diff_v2");
    grv2[kk]->SetMarkerStyle(20);
    grv2[kk]->SetMarkerSize(1);
    grv2[kk]->SetMarkerColorAlpha(kBlue, 1);
    grv2[kk]->SetLineColorAlpha(kBlue, 1);
    grv2[kk]->SetLineWidth(2);

    grv4IMP[kk] = new TGraphErrors(NN, binPt[kk], IMPdv4[kk], RMSbinPt[kk], RMSv4[kk]);
    grv4IMP[kk]->SetName("IMPdiff_v4");
grv4IMP[kk]->GetYaxis()->SetRangeUser(0, 0.25);
    grv4IMP[kk]->SetMarkerStyle(25);
    grv4IMP[kk]->SetMarkerSize(1);
    grv4IMP[kk]->SetMarkerColorAlpha(kBlack, 1);
    grv4IMP[kk]->SetLineColorAlpha(kBlack, 1);
    grv4IMP[kk]->SetLineWidth(1);
grv4IMP[kk]->SetTitle("");

    grv2IMP[kk] = new TGraphErrors(NN, binPt[kk], IMPdv2[kk], RMSbinPt[kk], RMSv2[kk]);
    grv2IMP[kk]->SetName("IMPdiff_v2");
    grv2IMP[kk]->SetMarkerStyle(24);
    grv2IMP[kk]->SetMarkerSize(1);
    grv2IMP[kk]->SetMarkerColorAlpha(kBlue, 1);
    grv2IMP[kk]->SetLineColorAlpha(kBlue, 1);
    grv2IMP[kk]->SetLineWidth(2);

    grvMc[kk] = new TGraphErrors(NN, binPt[kk], mc[kk], RMSbinPt[kk], RMSmc[kk]);
    grvMc[kk]->SetName("diff_vMC");
    grvMc[kk]->SetMarkerStyle(22);
    grvMc[kk]->SetMarkerSize(1);
    grvMc[kk]->SetMarkerColorAlpha(kGreen, 1);
    grvMc[kk]->SetLineColorAlpha(kGreen, 1);
    grvMc[kk]->SetLineWidth(2);
    //grvMc[kk]->Draw("SAME P");
    

    grvEp[kk] = new TGraphErrors(NN, binPt[kk], DiffvEP[kk], RMSbinPt[kk], RMSvep[kk]);
    grvEp[kk]->SetName("diff_vEVENTplane");
    grvEp[kk]->SetMarkerStyle(22);
    grvEp[kk]->SetMarkerSize(1);
    grvEp[kk]->SetMarkerColorAlpha(kRed, 1);
    grvEp[kk]->SetLineColorAlpha(kRed, 1);
    grvEp[kk]->SetLineWidth(2);
    //grvEp[kk]->Draw("SAME P");
    

    //Рефренсный поток с погрешностями
    double x[1] = {0.5}, y[1] = {vMC[kk]}, ex[1] = {0.}, ey[1] = {ErMc[kk]};
    grIntvMc[kk] = new TGraphErrors(1, x, y, ex, ey);
    grIntvMc[kk]->SetName("vMC");
    grIntvMc[kk]->SetMarkerStyle(20);
    grIntvMc[kk]->GetYaxis()->SetTitle("V_{2}");

    double x2[1] = {1.5}, y2[1] = {v2[kk]}, ex2[1] = {0.}, ey2[1] = {ErV2[kk]};
    grIntv2[kk] = new TGraphErrors(1, x2, y2, ex2, ey2);
    grIntv2[kk]->SetName("v2");
    grIntv2[kk]->SetMarkerStyle(20);
y2[0] = IMPv2[kk];
    grIntv2IMP[kk] = new TGraphErrors(1, x2, y2, ex2, ey2);
    grIntv2IMP[kk]->SetName("v2");
    grIntv2IMP[kk]->SetMarkerStyle(24);

    double x4[1] = {2.5}, y4[1] = {v4[kk]}, ex4[1] = {0.}, ey4[1] = {ErV4[kk]};
    
    grIntv4[kk] = new TGraphErrors(1, x4, y4, ex4, ey4);
    grIntv4[kk]->SetName("v4");
    grIntv4[kk]->SetMarkerStyle(20);
y4[0] = IMPv4[kk];
grIntv4IMP[kk] = new TGraphErrors(1, x4, y4, ex4, ey4);
    grIntv4IMP[kk]->SetName("v4");
    grIntv4IMP[kk]->SetMarkerStyle(24);
    //cout << "MC eror " << ey[0] << " SV2 eror " << ey2[0] << " SV4 eror " << ey4[0] << endl;

    double x5[1] = {3.5}, y5[1] = {vEP[kk]}, ex5[1] = {0.}, ey5[1] = {ErEv[kk]};
    grIntvEp[kk] = new TGraphErrors(1, x5, y5, ex5, ey5);
    grIntvEp[kk]->SetName("vEP");
    grIntvEp[kk]->SetMarkerStyle(20);
    grIntvEp[kk]->GetYaxis()->SetTitle("V_{2}");

    }


    /*gPad->BuildLegend();
    gROOT->SetStyle("Pub");
    gStyle->SetOptStat(1111);*/

static const int cent=3;
TCanvas *c1 = new TCanvas("c1", "defferen. Flow", 10, 10, 600, 400);
  c1->SetLeftMargin(0.12);c1->SetRightMargin(0.12);
    c1->SetBottomMargin(0.1);
    c1->SetTopMargin(0.03);
double ymin0 = 0;//min(Diffv4[cent][NN-1] - RMSv4[cent][NN-1],0.);
    double ymax0 = 0.25;//max(Diffv4[cent][NN-1] + RMSv4[cent][NN-1], Diffv4[cent][NN-2] + RMSv4[cent][NN-1]);
    grv4[cent]->GetYaxis()->SetRangeUser(0.99 * ymin0, 1.05 * ymax0);

//grv4[cent]->Draw("AP");
//grv2[cent]->Draw("SAME P");
grv4IMP[cent]->Draw("AP");
grv2IMP[cent]->Draw("SAME P");
grvMc[cent]->Draw("SAME P");
grvEp[cent]->Draw("SAME P");


TLegend *leg1 = new TLegend(0.19,0.75,0.32,0.95);
   leg1->AddEntry(grv4IMP[cent],"V_{2}{4}","pe");
   leg1->AddEntry(grv2IMP[cent],"V_{2}{2}","pe");
   leg1->AddEntry(grvMc[cent],"V_{2}{MC}","pe");
   leg1->AddEntry(grvEp[cent],"V_{2}{EP}","pe");
sprintf(str, "cent: %i-%i%",cent * 10, (cent + 1) * 10);
        STR = (char *)str;
leg1->SetHeader(STR);
   // leg3 -> SetTextFont(62);
    leg1 -> SetTextSize(0.03);
    leg1 -> SetTextAlign(22);
  leg1 -> SetBorderSize(1);
leg1 ->Draw();






TLegend *legdif[Nb];

TCanvas *cdif = new TCanvas("cdif", "def Flow Flow", 10, 10, 800, 600);
    cdif->SetLeftMargin(0.11);cdif->SetRightMargin(0.11);
    cdif->SetBottomMargin(0.1);
    cdif->SetTopMargin(0.03);
cdif->Divide(3,2,0.001,0.001);char text1[800];

for(int k=0;k<6;k++){cdif->cd(k+1);



legdif[k]= new TLegend(0.1,0.65,0.32,0.9);
   legdif[k]->AddEntry(grv4IMP[cent],"V_{2}{4}","pe");
   legdif[k]->AddEntry(grv2IMP[cent],"V_{2}{2}","pe");
   legdif[k]->AddEntry(grvMc[cent],"V_{2}{MC}","pe");
   legdif[k]->AddEntry(grvEp[cent],"V_{2}{EP}","pe");

    
sprintf(str, "cent: %d-%d%",k * 10, (k + 1) * 10);
        STR = (char *)str;

legdif[k]->SetHeader(STR);
grv4IMP[k]->Draw("AP");
grv2IMP[k]->Draw("SAME P");
grvMc[k]->Draw("SAME P");
grvEp[k]->Draw("SAME P");legdif[k]-> Draw();}

TFile *d_outfile = new TFile(savefile, "recreate");
    d_outfile->cd();


    grvEp[cent]->Write("grMC");
    grvEp[cent]->Write("grEP");
    grv2IMP[cent]->Write("grV2");
    grv4IMP[cent]->Write("grV4");
    d_outfile->Close();


TCanvas *c2 = new TCanvas("c2", "Reference Flow", 10, 10, 600, 400);
    c2->SetLeftMargin(0.11);c2->SetRightMargin(0.11);
    c2->SetBottomMargin(0.1);
    c2->SetTopMargin(0.03);
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
    double ey4[1] ={ErV4[cent]};
    //ey4[0]= grIntv4[kk]->GetErrorY(0);
    double ymin = min(vMC[cent] - ey4[0], v4[cent] - ey4[0]);
    double ymax = max(v4[cent] + ey4[0], v2[cent] + ey4[0]);
    h->GetYaxis()->SetRangeUser(0.95 * ymin, 1.05 * ymax);
    h->Draw();

    TLine *line = new TLine(0., vMC[cent], 4., vMC[cent]);
    line->SetLineWidth(1);
    line->SetLineStyle(5);

TLegend *leg2 = new TLegend(0.19,0.75,0.32,0.95);
   leg2->AddEntry(grIntv4[cent],"V_{2}{4}","pe");
   leg2->AddEntry(grIntv2[cent],"V_{2}{2}","pe");
   leg2->AddEntry(grIntv4IMP[cent],"V_{2}{4} acept","pe");
   leg2->AddEntry(grIntv2IMP[cent],"V_{2}{2} acept","pe");
   leg2->AddEntry(grIntvMc[cent],"V_{2}{MC}","pe");
   leg2->AddEntry(grIntvEp[cent],"V_{2}{EP}","pe");
sprintf(str, "cent: %i-%i%",cent * 10, (cent + 1) * 10);
        STR = (char *)str;
leg2->SetHeader(STR);
   // leg3 -> SetTextFont(62);
    leg2 -> SetTextSize(0.03);
    leg2 -> SetTextAlign(22);
  leg2 -> SetBorderSize(1);

    grIntvMc[cent]->Draw("SAME P");
    grIntvEp[cent]->Draw("SAME P");
    grIntv4[cent]->Draw("SAME P");
    grIntv2[cent]->Draw("SAME P");
grIntv4IMP[cent]->Draw("SAME P");
    grIntv2IMP[cent]->Draw("SAME P");
    line->Draw("SAME");
leg2 ->Draw();
cout << "Q2Er "<<ErV2[cent]<<"Q4Er "<<ErV4[cent]<<"McEr "<<ErMc[cent]<<endl;


    TCanvas *c3 = new TCanvas("c3", "Reference Flow by Centrality", 10, 10, 600, 400);
    c3->SetLeftMargin(0.1);c3->SetRightMargin(0.1);
    c3->SetBottomMargin(0.1);
    c3->SetTopMargin(0.03);
    TGraphErrors *grIv4 = new TGraphErrors(Nb, binCent, IntV4, RMSbinCent, ErV4);
    grIv4->SetMarkerStyle(21);
    grIv4->SetMarkerSize(1);
    grIv4->SetMarkerColorAlpha(kBlack, 1);
    grIv4->SetLineColorAlpha(kBlack, 1);
    grIv4->SetLineWidth(2);
    grIv4->GetXaxis()->SetTitle("Cent, %");
    grIv4->GetYaxis()->SetTitle("V_{2}");
    //grv4[kk]->Draw("AP");
    grIv4->SetTitle("");

TGraphErrors *grIIMPv4 = new TGraphErrors(Nb, binCent, IntIMPv4, RMSbinCent, ErV4);
    grIIMPv4->SetMarkerStyle(25);
    grIIMPv4->SetMarkerSize(1);
    grIIMPv4->SetMarkerColorAlpha(kBlack, 1);
    grIIMPv4->SetLineColorAlpha(kBlack, 1);
    grIIMPv4->SetLineWidth(2);
     grIIMPv4->SetTitle("");


TGraphErrors *grIv2 = new TGraphErrors(Nb, binCent, IntV2, RMSbinCent, ErV2);
    grIv2->SetMarkerStyle(20);
    grIv2->SetMarkerSize(1);
    grIv2->SetMarkerColorAlpha(kBlue, 1);
    grIv2->SetLineColorAlpha(kBlue, 1);
    grIv2->SetLineWidth(2);
    //grIv2->SetTitle("V_{2}{2}");
TGraphErrors *grIIMPv2 = new TGraphErrors(Nb, binCent, IntIMPv2, RMSbinCent, ErV2);
    grIIMPv2->SetMarkerStyle(24);
    grIIMPv2->SetMarkerSize(1);
    grIIMPv2->SetMarkerColorAlpha(kBlue, 1);
    grIIMPv2->SetLineColorAlpha(kBlue, 1);
    grIIMPv2->SetLineWidth(2);

TGraphErrors *grIvMc = new TGraphErrors(Nb, binCent, IntMc, RMSbinCent, ErMc);
    grIvMc->GetYaxis()->SetRangeUser(0, 0.15);
    grIvMc->SetMarkerStyle(21);
    grIvMc->SetMarkerSize(1);
    grIvMc->SetMarkerColorAlpha(kGreen, 1);
    grIvMc->SetLineColorAlpha(kGreen, 1);
    grIvMc->SetLineWidth(1);
    grIvMc->GetYaxis()->SetTitle("V_{2}");
    //grIvMc->SetTitle("V_{2}{MC}");

    TGraphErrors *grIvEp = new TGraphErrors(Nb, binCent, IntEv, RMSbinCent, ErEv);
    grIvEp->SetName("diff_vEVENTplane");
    grIvEp->SetMarkerStyle(22);
    grIvEp->SetMarkerSize(1);
    grIvEp->SetMarkerColorAlpha(kRed, 1);
    grIvEp->SetLineColorAlpha(kRed, 1);
    grIvEp->SetLineWidth(1);
    //grIvEp->SetTitle("V_{2}{Event Plane}");

TLegend *leg3 = new TLegend(0.19,0.75,0.29,0.95);
   leg3->AddEntry(grIIMPv4,"V_{2}{4}","pe");
   leg3->AddEntry(grIIMPv2,"V_{2}{2}","pe");
   leg3->AddEntry(grIvMc,"V_{2}{MC}","pe");
   leg3->AddEntry(grIvEp,"V_{2}{EP}","pe");
   // leg3 -> SetTextFont(62);
    leg3 -> SetTextSize(0.03);
    leg3 -> SetTextAlign(22);
  leg3 -> SetBorderSize(1);
//grIv4->Draw("AP");
//grIv2->Draw("SAME P");
grIIMPv4->Draw("AP");
grIIMPv2->Draw("SAME P");
    grIvMc->Draw("SAME P");
    grIvEp->Draw("SAME P");
leg3 -> Draw();



}

