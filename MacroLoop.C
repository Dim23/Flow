
#include "MyClass.C"
#include "READCent.C"
//#include "READCent.C"
void start(const char *file_adress="/home/dim2/FLOW5/OUT/Vin1nonflow.root",const char *outfile="~/FLOW5/OUT/OUT.root"){
MyClass t(file_adress);
t.Book();//Centering("");
//Flatting("~/FLOW5/OUT/OUTcent.root");
t.Loop();t.SaveData(outfile);
}
void MacroLoop(){//start();
//read("/home/dim2/FLOW5/Macro/claster/10MsumALLnonflow.root");
//read("/home/dim2/FLOW5/Macro/claster/10MsumALLwiegt.root");
read("/home/dim2/FLOW5/OUT/OUT.root");


}


