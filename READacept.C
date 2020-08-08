#include <TF1.h>
#include <TLine.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGraph.h>

double sigmaX(TH1F *X,double weight){
double RMS=X->GetStdDev();
double neff=X->GetEffectiveEntries();
return weight*RMS*RMS/(neff-1);
}

double sigmaXY(TH1F *X,TH1F *Y,TH1F *XY,double weight){
double RMS=XY->GetMean()-(X->GetMean())*(Y->GetMean());
double Xsum=X->GetSum();
double Ysum=Y->GetSum();
double XYsum=XY->GetSum();
double neff=Xsum*Ysum/XYsum;
return weight*RMS*RMS/(neff-1);
}
void read(const char *outfile,const char *savefile="~/FLOW5/PLOT/plot30_40.root"){
//double pt_bin[11]={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};int NN=10;
//double pt_bin[15]={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2};int NN=14;//бины пот импульсу для дифференциального потока
double pt_bin[25]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.5}; static const int NN=24;
double RMSSV2,RMSSV4,RMSVMC;//RMS для SV2 и SV4


double *IMPdn2=new double[NN];double *IMPdv2=new double[NN];
double *IMPdn4=new double[NN];double *IMPdv4=new double[NN];
double *dCOSp1=new double[NN];double *dSINp1=new double[NN];
double *dCOSp1f2=new double[NN];double *dSINp1f2=new double[NN];
double *dCOSp1f2mf3=new double[NN];double *dSINp1f2mf3=new double[NN];
double *dCOSp1mf2mf3=new double[NN];double *dSINp1mf2mf3=new double[NN];
double *RMSdn2=new double[NN];double *RMSv2=new double[NN];
double *RMSdn4=new double[NN];double *RMSv4=new double[NN];
double *Rv2=new double[NN];double *Rv4=new double[NN];
double *mc=new double[NN];double *RMSmc=new double[NN];
double *RMSbinPt=new double[NN];

double sv4=0,sv2=0,v4,v2,vMC,cn2,cn4;
double IMPsv4=0,IMPsv2=0,IMPv4,IMPv2,IMPcn2,IMPcn4;

TH1F *SV2;//2 частичная референсная корреляция 
TH1F *SV4;//4 частичная референсная корреляция
TH1F *SV2_SV4;
TH1F *HMC;
TH1F *COSf1;
TH1F *SINf1;
TH1F *COSf1f2;
TH1F *SINf1f2;
TH1F *COSf1f2f3;
TH1F *SINf1f2f3;

TH1F *DiffSV2[NN];TH1F *DiffSV4[NN];TH1F *DiffMC[NN];
TH1F *SV2_DiffSV2[NN];TH1F *SV2_DiffSV4[NN];
TH1F *SV4_DiffSV2[NN];TH1F *SV4_DiffSV4[NN];
TH1F *DiffSV2_DiffSV4[NN];

double *binPt=new double[NN];
double *QMC=new double[NN];
double *dn2=new double[NN];double *dn4=new double[NN];double *dsv4=new double[NN];
double *Diffv2=new double[NN]; double *Diffv4=new double[NN];
double *Diffsv2=new double[NN]; double *Diffsv4=new double[NN];

char strCOVsv2sv2[20];char strCOVsv2sv4[20];
char strCOVsv4sv2[20];char strCOVsv4sv4[20];
char DiffCOVsv2sv4[20];
char strE[20];char strW[20];char strMC[20];

char strCOSp1[20];char strSINp1[20];
char strCOSp1f2[20];char strSINp1f2[20];
char strCOSp1f2mf3[20];char strSINp1f2mf3[20];
char strCOSp1mf2mf3[20];char strSINp1mf2mf3[20];

TFile *f=new TFile(outfile);
SV2=(TH1F*)f->Get("sv2");
SV4=(TH1F*)f->Get("sv4");
SV2_SV4=(TH1F*)f->Get("sv2_sv4");
HMC=(TH1F*)f->Get("mc");
COSf1=(TH1F*)f->Get("cosf1");
COSf1f2=(TH1F*)f->Get("cosf1f2");
COSf1f2f3=(TH1F*)f->Get("cosf1f2f3");
SINf1=(TH1F*)f->Get("sinf1");
SINf1f2=(TH1F*)f->Get("sinf1f2");
SINf1f2f3=(TH1F*)f->Get("sinf1f2f3");

for(int m=0;m<NN;m++){
sprintf(strE,"%s%d","sv2Diff",m);
sprintf(strW,"%s%d","sv4Diff",m);
sprintf(strMC,"%s%d","DiffMC",m);
sprintf(strCOVsv2sv2,"%s%d","SV2_DiffSV2",m);
sprintf(strCOVsv2sv4,"%s%d","SV2_DiffSV4",m);
sprintf(strCOVsv4sv2,"%s%d","SV4_DiffSV2",m);
sprintf(strCOVsv4sv4,"%s%d","SV4_DiffSV4",m);
sprintf(DiffCOVsv2sv4,"%s%d","DiffSV2_DiffSV4",m);
sprintf(strCOSp1,"%s%d","COSp1",m);
sprintf(strCOSp1f2,"%s%d","COSp1f2",m);
sprintf(strCOSp1f2mf3,"%s%d","COSp1f2mf3",m);
sprintf(strCOSp1mf2mf3,"%s%d","COSp1mf2mf3",m);
sprintf(strSINp1,"%s%d","SINp1",m);
sprintf(strSINp1f2,"%s%d","SINp1f2",m);
sprintf(strSINp1f2mf3,"%s%d","SINp1f2mf3",m);
sprintf(strSINp1mf2mf3,"%s%d","SINp1mf2mf3",m);

const char *EastH=(char*)strE;
const char *WastH=(char*)strW;
const char *SMC=(char*)strMC;
const char *stCOVsv2sv2=(char*)strCOVsv2sv2;
const char *stCOVsv2sv4=(char*)strCOVsv2sv4;
const char *stCOVsv4sv2=(char*)strCOVsv4sv2;
const char *stCOVsv4sv4=(char*)strCOVsv4sv4;
const char *DifCOVsv2sv4=(char*)DiffCOVsv2sv4;
const char *stCOSp1=(char*)strCOSp1;
const char *stCOSp1f2=(char*)strCOSp1f2;
const char *stCOSp1f2mf3=(char*)strCOSp1f2mf3;
const char *stCOSp1mf2mf3=(char*)strCOSp1mf2mf3;
const char *stSINp1=(char*)strSINp1;
const char *stSINp1f2=(char*)strSINp1f2;
const char *stSINp1f2mf3=(char*)strSINp1f2mf3;
const char *stSINp1mf2mf3=(char*)strSINp1mf2mf3;
DiffSV2[m]=(TH1F*)f->Get(EastH);
DiffSV4[m]=(TH1F*)f->Get(WastH);
DiffMC[m]=(TH1F*)f->Get(SMC);
SV2_DiffSV2[m]=(TH1F*)f->Get(stCOVsv2sv2);
SV2_DiffSV4[m]=(TH1F*)f->Get(stCOVsv2sv4);
SV4_DiffSV2[m]=(TH1F*)f->Get(stCOVsv4sv2);
SV4_DiffSV4[m]=(TH1F*)f->Get(stCOVsv4sv4);
DiffSV2_DiffSV4[m]=(TH1F*)f->Get(DifCOVsv2sv4);
COSp1[m]=(TH1F*)f->Get(stCOSp1);
COSp1f2[m]=(TH1F*)f->Get(stCOSp1f2);
COSp1f2mf3[m]=(TH1F*)f->Get(stCOSp1f2mf3);
COSp1mf2mf3[m]=(TH1F*)f->Get(stCOSp1mf2mf3);
SINp1[m]=(TH1F*)f->Get(stSINp1);
SINp1f2[m]=(TH1F*)f->Get(stSINp1f2);
SINp1f2mf3[m]=(TH1F*)f->Get(stSINp1f2mf3);
SINp1mf2mf3[m]=(TH1F*)f->Get(stSINp1mf2mf3);
//cout <<"Mean  "<<m<<endl;
dCOSp1[m]=COSp1[m]->GetMean();
dCOSp1f2[m]=COSp1f2[m]->GetMean();
dCOSp1f2mf3[m]=COSp1f2mf3[m]->GetMean();
dCOSp1mf2mf3[m]=COSp1mf2mf3[m]->GetMean();
dSINp1[m]=SINp1[m]->GetMean();
dSINp1f2[m]=SINp1f2[m]->GetMean();
dSINp1f2mf3[m]=SINp1f2mf3[m]->GetMean();
dSINp1mf2mf3[m]=SINp1mf2mf3[m]->GetMean();
RMSbinPt[m]=0;
}

double dCOSf1,dCOSf1f2,dCOSf1mf2,dCOSf1f2f3,dSINf1,dSINf1f2,dSINf1f2f3;
dCOSf1=COSf1->GetMean();
cout <<"Start read  "<<endl;
dCOSf1f2=COSf1f2->GetMean();
cout <<"Start read  "<<endl;
dCOSf1f2f3=COSf1f2f3->GetMean();
dSINf1=SINf1->GetMean();
dSINf1f2=SINf1f2->GetMean();
dSINf1f2f3=SINf1f2f3->GetMean();

//Референсный поток
SV2->Sumw2();SV4->Sumw2();
sv2=SV2->GetMean();v2=pow(fabs(sv2),0.5);
IMPsv2=sv2-(pow(COSf1->GetMean(),2)+pow(SINf1->GetMean(),2));
IMPv2=pow(fabs(IMPsv2),0.5);

sv4=SV4->GetMean();cn4=sv4-2*(sv2*sv2);v4=pow(fabs(cn4),0.25);
IMPcn4=sv4-2*(sv2*sv2)-4*dCOSf1*dCOSf1f2f3+4*dSINf1*dSINf1f2f3-pow(dCOSf1f2,2)-pow(dSINf1f2,2)+
4*dCOSf1f2*(pow(dCOSf1,2)-pow(dSINf1,2))+8*dSINf1f2*dSINf1*dCOSf1+
8*sv2*(pow(dCOSf1,2)+pow(dSINf1,2))-6*pow((pow(dCOSf1,2)+pow(dSINf1,2)),2);
IMPv4=pow(fabs(IMPcn4),0.25);

vMC=HMC->GetMean();
cout <<" v2{2} is "<<v2<<" v2{4} is "<<v4 <<endl; 
cout <<" v2{2} is "<<IMPv2<<" v2{4} is "<<IMPv4 <<endl;
//Дифференциальный поток

for(int m=0;m<NN;m++){
dn2[m]=DiffSV2[m]->GetMean();Diffv2[m]=dn2[m]/pow(fabs(sv2),0.5);
IMPdn2[m]=DiffSV2[m]->GetMean()-dCOSp1[m]*dCOSf1-dSINp1[m]*dSINf1;IMPdv2[m]=IMPdn2[m]/pow(fabs(IMPsv2),0.5);
dsv4[m]=DiffSV4[m]->GetMean();dn4[m]=dsv4[m]-2*dn2[m]*sv2;Diffv4[m]=-dn4[m]/pow(fabs(cn4),0.75);

IMPdn4[m]=dsv4[m]-2*dn2[m]*sv2-dCOSp1[m]*dCOSf1f2f3+dSINp1[m]*dSINf1f2f3-dCOSf1*dCOSp1mf2mf3[m]+dSINf1*dSINp1mf2mf3[m]-2*dCOSf1*dCOSp1f2mf3[m]-2*dSINf1*dSINp1f2mf3[m]-dCOSp1f2[m]*dCOSf1f2-dSINp1f2[m]*dSINf1f2+2*dCOSf1f2*(dCOSp1[m]*dCOSf1-dSINp1[m]*dSINf1)+2*dSINf1f2*(dCOSp1[m]*dSINf1+dSINp1[m]*dCOSf1)+4*sv2*(dCOSp1[m]*dCOSf1+dSINp1[m]*dSINf1)+2*dCOSp1f2[m]*(pow(dCOSf1,2)-pow(dSINf1,2))+4*dSINp1f2[m]*dCOSf1*dSINf1+4*dn2[m]*(pow(dCOSf1,2)+pow(dSINf1,2))-6*(pow(dCOSf1,2)-pow(dSINf1,2))*(dCOSp1[m]*dCOSf1-dSINp1[m]*dSINf1)-12*dCOSf1*dSINf1*(dSINp1[m]*dCOSf1+dCOSp1[m]*dSINf1);

IMPdv4[m]=-IMPdn4[m]/pow(fabs(IMPcn4),0.75);
binPt[m]=0.5*(pt_bin[m]+pt_bin[m+1]);
//cout <<"Different. flow "<<" v2{2} is "<<Diffv2[m]<<", v2{4} is "<<Diffv4[m]<<" Pt "<<binPt[m]<<endl ;
cout <<"Different. flow "<<" v2{2} is "<<IMPdv2[m]<<", v2{4} is "<<IMPdv4[m]<<" Pt "<<binPt[m]<<" m="<<m<<endl ;
//cout <<"Different. flow "<<" v2{2} is "<<dn2[m]<<", v2{4} is "<<dsv4[m]<<" Pt "<<binPt[m]<<endl ;
}
//погрешность для дифференциального потока

for(int m=0;m<NN;m++){
mc[m]=DiffMC[m]->GetMean();
RMSmc[m]=pow(sigmaX(DiffMC[m],1),0.5);RMSdn2[m]=sigmaX(DiffSV2[m],1);RMSdn4[m]=DiffSV4[m]->GetStdDev();

Rv2[m]=0.25*pow(sv2,-3)*(pow(dn2[m],2)*sigmaX(SV2,1)+4*sv2*sv2*sigmaX(DiffSV2[m],1)-4*dn2[m]*sv2*sigmaXY(SV2,DiffSV2[m],SV2_DiffSV2[m],1));
RMSv2[m]=pow(Rv2[m],0.5);

Rv4[m]=(pow(fabs(cn4),-3.5))*(pow((2*sv2*sv2*dn2[m]-3*sv2*dsv4[m]+2*sv4*dn2[m]),2)*sigmaX(SV2,1)+

9/16*pow(dn4[m],2)*sigmaX(SV4,1)+

4*sv2*sv2*pow(fabs(cn4),2)*sigmaX(DiffSV2[m],1)+

pow(fabs(cn4),2)*sigmaX(DiffSV4[m],1)-

1.5*fabs(dn4[m])*(2*sv2*sv2*dn2[m]-3*sv2*dsv4[m]+2*sv4*dn2[m])*sigmaXY(SV2,SV4,SV2_SV4,1)-

4*sv2*fabs(cn4)*(2*sv2*sv2*dn2[m]-3*sv2*dsv4[m]+2*sv4*dn2[m])*sigmaXY(SV2,DiffSV2[m],SV2_DiffSV2[m],1)+

2*fabs(cn4)*(2*sv2*sv2*dn2[m]-3*sv2*dsv4[m]+2*sv4*dn2[m])*sigmaXY(SV2,DiffSV4[m],SV2_DiffSV4[m],1)+

3*sv2*(fabs(cn4))*dn4[m]*sigmaXY(SV4,DiffSV2[m],SV4_DiffSV2[m],1)-

1.5*fabs(cn4)*sigmaXY(SV4,DiffSV4[m],SV4_DiffSV4[m],1)-

4*sv2*pow(fabs(cn4),2)*sigmaXY(DiffSV2[m],DiffSV4[m],DiffSV2_DiffSV4[m],1));

RMSv4[m]=pow(Rv4[m],0.5);
}


//TF1 *momentdist = new TF1("momentdist","[1]*x+[0]", 0., 2.1);
//TF1 *momentdist = new TF1("momentdist","[1]*(exp(2*x)-1)/(exp(2*x)+1)+[0]", 0., 2.3);
/*
momentdist->FixParameter(0, 0.02);
momentdist->FixParameter(1, 0.1);
momentdist->SetLineColorAlpha(kBlack, 0.8);
momentdist->SetLineWidth(2);momentdist->SetLineStyle(5);
momentdist->GetXaxis()->SetTitle("Pt,Gev/c");
momentdist->GetYaxis()->SetTitle("V_{2}");
momentdist->SetTitle("V_{2}(Pt)");
momentdist->Draw();*/

TFile *file = new TFile(savefile,"recreate");
file->cd();

TGraphErrors *impgrv4 = new TGraphErrors(NN,binPt,IMPdv4,RMSbinPt,RMSv4);
impgrv4->SetName("diff_v4");
impgrv4->GetYaxis()->SetRangeUser(0,0.15);
impgrv4->GetXaxis()->SetRangeUser(0,2.5);
impgrv4->SetMarkerStyle(21);
impgrv4->SetMarkerSize(1);
impgrv4->SetMarkerColorAlpha(kBlack, 1);
impgrv4->SetLineColorAlpha(kBlack, 1);
impgrv4->SetLineWidth(2);
impgrv4->GetXaxis()->SetTitle("Pt,Gev/c");
impgrv4->GetYaxis()->SetTitle("V_{2}");
impgrv4->Draw("AP");
impgrv4->SetTitle("V_{2}{4} accept. affects");


TGraphErrors *grv4 = new TGraphErrors(NN,binPt,Diffv4,RMSbinPt,RMSv4);
grv4->SetMarkerStyle(25);
grv4->SetMarkerSize(1);
grv4->SetMarkerColorAlpha(kBlack, 1);
grv4->SetLineColorAlpha(kBlack, 1);
grv4->SetLineWidth(2);
grv4->Draw("SAME P");
grv4->SetTitle("V_{2}{4}");

TGraphErrors *grv2 = new TGraphErrors(NN,binPt,Diffv2,RMSbinPt,RMSv2);
grv2->SetName("diff_v2");
grv2->SetMarkerStyle(24);
grv2->SetMarkerSize(1);
grv2->SetMarkerColorAlpha(kBlue, 1);
grv2->SetLineColorAlpha(kBlue, 1);
grv2->GetYaxis()->SetTitle("V_{2}");
grv2->GetXaxis()->SetTitle("Pt,Gev/c");
grv2->SetLineWidth(2);
grv2->Draw("SAME P");
grv2->SetTitle("V_{2}{2}");

TGraphErrors *impgrv2 = new TGraphErrors(NN,binPt,IMPdv2,RMSbinPt,RMSv2);
impgrv2->SetMarkerStyle(20);
impgrv2->SetMarkerSize(1);
impgrv2->SetMarkerColorAlpha(kBlue, 1);
impgrv2->SetLineColorAlpha(kBlue, 1);
impgrv2->GetYaxis()->SetTitle("V_{2}");
impgrv2->GetXaxis()->SetTitle("Pt,Gev/c");
impgrv2->SetLineWidth(2);
impgrv2->Draw("SAME P");
impgrv2->SetTitle("V_{2}{2} accept. affects");

TGraphErrors *grv = new TGraphErrors(NN,binPt,mc,RMSbinPt,RMSmc);
grv->SetName("diff_vMC");
grv->SetMarkerStyle(22);
grv->SetMarkerSize(1);
grv->SetMarkerColorAlpha(kGreen, 1);
grv->SetLineColorAlpha(kGreen, 1);
grv->SetLineWidth(1);
grv->Draw("SAME P");
grv->SetTitle("V_{2}{MC}");


gPad->BuildLegend();
gROOT->SetStyle("Pub");
gStyle->SetOptStat(1111);

//Рефренсный поток с погрешностями
TCanvas *c2=new TCanvas("c2","Graph Draw Options",650,500);
double x[1]={1},y[1]={vMC},ex[1]={0.},ey[1]={pow(sigmaX(HMC,1),0.5)};
TGraphErrors *grMC = new TGraphErrors(1,x,y,ex,ey);
grMC->SetName("vMC");
grMC->SetMarkerStyle(20);
grMC->GetYaxis()->SetTitle("V_{2}");

double x2[1]={1.4},y2[1]={v2},ex2[1]={0.},ey2[1]={pow(sigmaX(SV2,0.25/sv2),0.5)};
TGraphErrors *grsv2 = new TGraphErrors(1,x2,y2,ex2,ey2);
grsv2->SetName("v2");
grsv2->SetMarkerStyle(24);

double impx2[1]={1.4},impy2[1]={IMPv2},impex2[1]={0.},impey2[1]={pow(sigmaX(SV2,0.25/sv2),0.5)};
TGraphErrors *impgrsv2 = new TGraphErrors(1,impx2,impy2,impex2,impey2);
impgrsv2->SetMarkerStyle(20);

double x4[1]={1.8},y4[1]={v4},ex4[1]={0.},ey4[1]={0};
double s4=pow(fabs(cn4),-1.5)*(sigmaX(SV2,sv2*sv2)+sigmaX(SV4,0.125)-0.5*sv2*sigmaXY(SV2,SV4,SV2_SV4,1));
ey4[0]=pow(s4,0.5);
TGraphErrors *grsv4 = new TGraphErrors(1,x4,y4,ex4,ey4);
grsv4->SetName("v4");
grsv4->SetMarkerStyle(24);
cout <<"MC eror "<<ey[0]<<" SV2 eror "<<ey2[0]<<" SV4 eror "<<ey4[0]<<endl;
grMC->GetXaxis()->SetRangeUser(0,1.9);
grMC->GetYaxis()->SetRangeUser(y4[0]-0.003,y4[0]+0.003);

double impx4[1]={1.8},impy4[1]={IMPv4},impex4[1]={0.},impey4[1]={0};
double imps4=pow(fabs(cn4),-1.5)*(sigmaX(SV2,sv2*sv2)+sigmaX(SV4,0.125)-0.5*sv2*sigmaXY(SV2,SV4,SV2_SV4,1));
impey4[0]=pow(imps4,0.5);
TGraphErrors *impgrsv4 = new TGraphErrors(1,impx4,impy4,impex4,impey4);
impgrsv4->SetMarkerStyle(20);
grsv4->SetMarkerStyle(24);
cout <<"MC eror "<<ey[0]<<" SV2 eror "<<ey2[0]<<" SV4 eror "<<ey4[0]<<endl;
grMC->GetXaxis()->SetRangeUser(0,1.9);
grMC->GetYaxis()->SetRangeUser(impey4[0]-2*ey4[0],y2[0]+2*ey2[0]);

TLine *line=new TLine(0.9,vMC,1.9,vMC);
line->SetLineWidth(1);
line->SetLineStyle(5);
grMC->Draw("AP");
grsv4->Draw("SAME P");
grsv2->Draw("SAME P");
line->Draw("SAME");
impgrsv4->Draw("SAME P");
impgrsv2->Draw("SAME P");


grv2->Write();
grv->Write();
grMC->Write();
grsv4->Write();
grsv2->Write();
file->Close();
}
