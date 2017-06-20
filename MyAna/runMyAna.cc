#include "MyAna.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"

#include "StdArg.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

void Usage(const char* exename) {
  cout << "Usage : runMyAna -debug -filelist <data.list> -out <output.root> -num <Nevent> -skim" << endl;
}

int main(int argc, char *argv[])
{

// parse executable options :
//___________________________

  StdArg arg(argc,argv);
 
  cout << "======================================================================" << endl;
  for (int i=0;i<argc;i++)
    cout << argv[i] << " ";
  cout << endl;
  cout << "======================================================================" << endl;
   
  // enter all possible flags:
  arg.flags << "-debug" << "-skim" << "-mc"<<"-sig"<<"-ttbar";
 
  // enter all possible keys:
  arg.keys << "-filelist" << "-out" << "-num" << "-mode";
 
  bool debug   = false;
  bool skim    = false;
  bool issig   = false;
  string input  = "data.list";
  string output = "output_dummy.root";
  int nevent = -1;
  string mode   = "loop";

  try {  // parse command line
    arg.Process();
    debug = arg.Flag("-debug") ? true: false;
    skim  = arg.Flag("-skim")  ? true: false;
    issig  = arg.Flag("-sig")  ? true: false;
    if ( arg.Key("-filelist")    ) input  = arg.Get<string>("-filelist");
    if ( arg.Key("-out")         ) output = arg.Get<string>("-out");    
    if ( arg.Key("-num")         ) nevent = arg.Get<int>("-num");
    if ( arg.Key("-mode")        ) mode   = arg.Get<string>("-mode");
 }
  catch (StdArg::BadInput) {
    if (argc > 1) cout<< "Input error" <<endl;
    // usage in case of error or no parameters
    Usage(argv[0]);
    return 1;
  }

//
// Load ntuples or TCHAIN :
//_________________________
 
  TChain *tree = new TChain("events");
  TString fNameList = input.c_str();
  cout << "=> Load Chain from file: " << fNameList << endl;
  ifstream fList(fNameList.Data());
  if (!fList) {
    cout << "=> Can't open file " << fNameList << endl;
    return 1;
  }
  std::string lineFromFile;
  do {
    std::getline(fList, lineFromFile);
    if (lineFromFile=="") break;
    if(tree->Add(lineFromFile.c_str())) cout << "=>File '"
                                          << lineFromFile
                                          << "' has been loaded" << endl;
    else cout << ">>Can't load file '" << lineFromFile << "'" << endl;
    // datas files does not contain MC branch
  } while (fList.good());
  
  
    
  fList.close();

  cout << "======================================================================" << endl;
  cout << "Total Number of events = "  << tree->GetEntries() << endl;
  cout << "======================================================================" << endl;
  MyAna ana(tree); //instantiate one MyAna object
  ana.SetNevent(nevent);// configure the object of type MyAna
  ana.SetRootName(output);
  ana.SetDebugMode(debug);
  ana.SetDoSkim(skim);
  ana.SetSIGmode(issig);
  cout << "======================================================================" << endl;
  
  cout << "Calling MyAna::Loop()..." << endl;
  ana.Loop(); //Call the Loop method on the object ana

  return 0;
}
