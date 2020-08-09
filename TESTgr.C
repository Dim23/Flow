//double* test(){double *a=new double[2];a[0]=1;a[1]=3;return a;}

double start(double N0){TFile *f = TFile::Open("/home/dim2/FLOW5/All.root");
TH1I *hNpart, *hNcoll;hNpart = (TH1I *)f->Get("hNpart");double c;c=hNpart->Integral(N0,600)/hNpart->Integral();//hNpart->Draw();
return c;}

void start2(){TFile *fDim = TFile::Open("~/FLOW5/OUT/Vin1nonflowNEW.root");
TH1F *hNpart;hNpart = (TH1F *)fDim->Get("HphiGen");
TF1 *f1 = new TF1("f1","(1+2*[0]*cos(2*x)+2*[1]*cos(3*x))/(2*3.1415926)*[2]",0,6.28);
hNpart->Fit(f1,"R");
}
