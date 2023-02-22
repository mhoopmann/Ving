#include "MSReader.h"
#include "NeoPepXMLParser.h"
#include <algorithm>

#include "VingData.h"
//#include "vingStructs.h"
#include "VingParameters.h"

using namespace std;
using namespace MSToolkit;

#define BDATE "Feb 8 2023"
#define VERSION "0.7.6"

bool cmdLine(int argc, char* argv[]);
void exportResults2(const char* fn, vector<sMS2>& v);
void marquee();
void usage();

string g_param;
VingParameters params;

int main(int argc, char* argv[]){

  marquee();
  if(!cmdLine(argc,argv)){
    usage();
    return 1;
  }

  if(!params.parse(g_param.c_str())){
    cout << "Cannot read params file, or file contains bad parameters." << endl;
    return 2;
  }

  VingData data(&params);
  //data.parseMzML();
  data.importMS2SearchResults();
  data.importMS3SearchResults();
  data.assessXLType();

  data.exportResults2();
  data.exportProXL();
  //exportResults2(params.output.c_str(),data.groups);

  return 0;
}

bool cmdLine(int argc, char* argv[]){
  if(argc!=2) return false;
  g_param=argv[1];
  return true;
}

void exportResults2(const char* fn, vector<sMS2>& v){
  FILE* f = fopen(fn, "wt");

  //Massive Header
  fprintf(f, "GroupID\tDesignation\tProbability\tMS2ScanNum\tMS2mz\tMS2charge\tMS2monoMass\tCalcNeutMass\tPPM\tSequence\tProtein(s)\tMS2Peptide\tMS2Protein\tMS2CalcNeutMass\tMS2Prob");
  for(size_t a=0;a<4;a++){
    fprintf(f,"\tMS3ScanNum-%d\tMS3mz\tMS3charge\tMS3selectedMass\tMS3Peptide\tMS3Protein\tMS3CalcNeutMass\tMS3Prob",(int)a+1);
  }
  fprintf(f,"\n");

  for(size_t a=0;a<v.size();a++){
    fprintf(f,"%d",(int)a);
    switch(v[a].type){
    case xlDeadEnd: fprintf(f,"\tDeadEnd"); break;
    case xlSingle: fprintf(f, "\tSingle"); break;
    case xlLoop: fprintf(f, "\tLoop"); break;
    case xlXL: fprintf(f, "\tXL"); break;
    case xlIncomplete: fprintf(f, "\tIncomplete"); break;
    default: fprintf(f, "\tUnknown"); break;
    }
    fprintf(f,"\t%.4lf",v[a].probability);
    fprintf(f, "\t%d", v[a].scan);
    fprintf(f, "\t%.4lf", v[a].mz);
    fprintf(f, "\t%d", v[a].charge);
    fprintf(f, "\t%.4lf", v[a].monoMZ*v[a].charge-v[a].charge*1.007276466);
    fprintf(f,"\t%.4lf",v[a].calcNeutMassG);
    fprintf(f,"\t%.2lf",v[a].ppmG);
    fprintf(f, "\t%s", v[a].sequence.c_str());
    fprintf(f, "\t%s", v[a].proteinS.c_str());
    if(!v[a].peptide.empty()){
      fprintf(f, "\t%s", v[a].peptide.c_str());
      fprintf(f, "\t%s", v[a].protein.c_str());
      fprintf(f, "\t%.4lf", v[a].calcNeutMass);
      fprintf(f, "\t%.4lf", v[a].probabilityMS2);
    } else fprintf(f,"\tn/a\tn/a\t0\t0");

    //And every MS3
    size_t b=0;
    while(b<v[a].ms3.size()){
      fprintf(f, "\t%d", v[a].ms3[b].scan);
      fprintf(f, "\t%.4lf", v[a].ms3[b].mz);
      fprintf(f, "\t%d", v[a].ms3[b].charge);
      fprintf(f, "\t%.4lf", v[a].ms3[b].mz * v[a].ms3[b].charge - v[a].ms3[b].charge * 1.007276466);
      fprintf(f, "\t%s", v[a].ms3[b].peptide.c_str());
      fprintf(f,"\t%s",v[a].ms3[b].protein.c_str());
      fprintf(f, "\t%.4lf", v[a].ms3[b].pepMass);
      fprintf(f, "\t%.4lf", v[a].ms3[b].prob);
      b++;
    }
    while(b<4){ //blanks
      fprintf(f,"\tn/a\t0\t0\t0\tn/a\tn/a\t0\t0");
      b++;
    }
    fprintf(f,"\n");
  }

  fclose(f);

}

void marquee() {
  cout << "Ving: Some sort of cleavable crosslinker analysis program" << endl;
  cout << "Copyright Michael Hoopmann, Institute for Systems Biology" << endl;
  cout << "Version: " << VERSION << endl;
  cout << "Date: " << BDATE << endl;
}

void usage(){
  cout << "\nUSAGE: Ving.exe <mzML> <MS2 PepXML> <MS3 PepXML> <outfile>" << endl;
  cout << "  mzML = mzML file of the spectra analyzed." << endl;
  cout << "  PepXML = PepXML of the MS2 or MS3 scan database search results with PeptideProphet analysis." << endl;
  cout << "  outfile = desired output file name." << endl;
}
