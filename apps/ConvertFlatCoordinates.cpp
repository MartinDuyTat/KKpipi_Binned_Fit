// Martin Duy Tat 19th December 2020
/**
 * ConvertFlatCoordinates is the program for converting the coordinates of flat phase space samples from a rectangular basis to invariant masses, which will be used for reweighting purposes
 * @param 1 Name of input ROOT file
 * @param 2 Name of output ROOt file
 */

#include<iostream>
#include<string>
#include<vector>
#include"TTree.h"
#include"TChain.h"
#include"TFile.h"
#include"TLorentzVector.h"
#include"KKpipiMath.h"

int main(int argc, char *argv[]) {
  if(argc != 3) {
    std::cout << "Incorrect number of inputs\n";
    return 0;
  }
  std::cout << "Loading input ROOT file...";
  TChain Chain("FlatEvents");
  Chain.Add(argv[1]);
  std::vector<double> X(5), S(5);
  for(int i = 0; i < 5; i++) {
    Chain.SetBranchAddress(("x" + std::to_string(i + 1)).c_str(), &X[i]);
  }
  TFile Outfile(argv[2], "RECREATE");
  std::cout << "Input ROOT file loaded\n";
  std::cout << "Setting up output ROOT file...\n";
  TTree Tree("FlatEvents", "FlatEvents");
  std::vector<std::string> ss{"s01", "s03", "s12", "s23", "s012"};
  for(int i = 0; i < 5; i++) {
    Tree.Branch(ss[i].c_str(), &S[i]);
  }
  std::cout << "Ready to process events\n";
  for(int i = 0; i < Chain.GetEntries(); i++) {
    Chain.GetEntry(i);
    std::vector<TLorentzVector> P = KKpipiMath::ConvertEventTo4Vectors(KKpipiMath::ConvertXToMomenta(X));
    S[0] = (P[0] + P[1]).M2();
    S[1] = (P[0] + P[3]).M2();
    S[2] = (P[1] + P[2]).M2();
    S[3] = (P[2] + P[3]).M2();
    S[4] = (P[0] + P[1] + P[2]).M2();
    Tree.Fill();
  }
  std::cout << "Conversion done\n";
  Tree.Write();
  Outfile.Close();
  return 0;
}
