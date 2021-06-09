#include "MSReader.h"
#include "NeoPepXMLParser.h"
#include <algorithm>

using namespace std;
using namespace MSToolkit;

typedef struct sMS3{
  int scan;
  string peptide;
  string protein;
  double xcorr;
  double eval;
  double mz;
  double prob;
  double pepMass;
  int charge;
  sMS3(){
    scan=0;
    xcorr=0;
    eval=0;
    mz=0;
    charge=0;
    pepMass=0;
  }
} sMS3;

typedef struct sMS2{
  int scan;
  double mz;
  int charge;
  double monoMZ;
  vector<sMS3> ms3;
  sMS2(){ clear();}
  void clear(){
    scan=0;
    mz=0;
    charge=0;
    monoMZ=0;
    ms3.clear();
  }
} sMS2;

typedef struct sXLPep{
  size_t index;
  double mass;
} sXLPep;

void addXML(char* fn, vector<sMS2>& v);
void exportResults(char* fn, vector<sMS2>& v);
void marquee();
void usage();

bool compareXL(sXLPep& a, sXLPep& b);

int main(int argc, char* argv[]){

  marquee();
  if(argc!=4){
    usage();
    return 1;
  }

  MSReader r;
  Spectrum s;

  r.setFilter(MS2);
  r.addFilter(MS3);
  r.addFilter(MS1);

  vector<sMS2> groups;
  vector<sMS2> cycle;
  

  r.readFile(argv[1],s);
  while(s.getScanNumber()>0){

    if(s.getMsLevel()==1){
      for(size_t a=0;a<cycle.size();a++){
        if(cycle[a].ms3.size()>0) groups.push_back(cycle[a]);
      }
      cycle.clear();

    } else if(s.getMsLevel()==2){
      sMS2 ms2;
      ms2.scan=s.getScanNumber();
      ms2.charge=s.getCharge();
      ms2.mz=s.getMZ();
      ms2.monoMZ=s.getMonoMZ();
      cycle.push_back(ms2);
    } else {
      sMS3 ms3;
      ms3.scan=s.getScanNumber();
      ms3.charge=s.getCharge();
      ms3.mz=s.getMZ();
      double pMZ=0;
      char header[256];
      s.getRawFilter(header,256);
      char* tok=strtok(header," ");
      while(tok!=NULL){
        string st=tok;
        size_t pos=st.find_first_of('@');
        if(pos!=string::npos){
          st=st.substr(0,pos);
          pMZ=atof(st.c_str());
          break;
        }
        tok = strtok(NULL, " ");
      }
      if(pMZ==0) cout << "Cannot extract MS2 selected ion" << endl;
      size_t a=0;
      for(a=0;a<cycle.size();a++){
        if(fabs(cycle[a].mz-pMZ)<0.01) break;
      }
      if(a==cycle.size()) cout << "Cannot find parent MS2 scan." << endl;
      else cycle[a].ms3.push_back(ms3);
    }

    r.readFile(NULL,s);
  }
  for (size_t a = 0; a < cycle.size(); a++) {
    if (cycle[a].ms3.size() > 0) groups.push_back(cycle[a]);
  }

  cout << groups.size() << " total groups." << endl;

  addXML(argv[2],groups);
  exportResults(argv[3],groups);

  return 0;
}

void addXML(char* fn, vector<sMS2>& v){
  NeoPepXMLParser p;
  p.read(fn);

  cout << "Done reading" << endl;

  size_t vp=0;
  size_t v3=0;

  for(size_t a=0;a<p.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query.size();a++){
    CnpxSpectrumQuery* sq=&p.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query[a];
   
    char* tok;
    char str[256];
    strcpy(str,sq->spectrum.c_str());
    tok=strtok(str,".\n\r");
    tok=strtok(NULL,".\n\r");
    int scan=atoi(tok);

    while(vp<v.size() && v[vp].ms3[v3].scan<scan){
      v3++;
      if(v3==v[vp].ms3.size()){
        v3=0;
        vp++;
      }
    }

    if(v[vp].ms3[v3].scan!=scan){
      cout << "No match to " << scan << endl;
      continue;
    }

    if(sq->search_result[0].search_hit.size()==0) continue;

    CnpxSearchHit* sh= &sq->search_result[0].search_hit[0];
    if(sh->modification_info.modified_peptide.empty()) v[vp].ms3[v3].peptide=sh->peptide;
    else v[vp].ms3[v3].peptide=sh->modification_info.modified_peptide;

    v[vp].ms3[v3].protein=sh->protein;

    v[vp].ms3[v3].charge=sq->assumed_charge;

    v[vp].ms3[v3].pepMass=sh->calc_neutral_pep_mass;

    for(size_t b=0;b<sh->search_score.size();b++){
      CnpxSearchScore* ss=&sh->search_score[b];
      if(ss->name.compare("xcorr")==0) v[vp].ms3[v3].xcorr=atof(ss->value.c_str());
      if (ss->name.compare("expect") == 0) v[vp].ms3[v3].eval = atof(ss->value.c_str());
    }

    for(size_t b=0;b<sh->analysis_result.size();b++){
      if(sh->analysis_result[b].analysis.compare("peptideprophet")==0) {
        v[vp].ms3[v3].prob=sh->analysis_result[b].peptide_prophet_result.probability;
      }
    }
  }

}

void exportResults(char* fn, vector<sMS2>& v){
  FILE* f = fopen(fn, "wt");

  fprintf(f,"GroupID\tMS2ScanNum\tMS2mz\tMS2charge\tMS2monoMass\tMS3ScanCount\n");;

  for (size_t a = 0; a < v.size(); a++) {
    int lowE = 0;
    for (size_t b = 0; b < v[a].ms3.size(); b++) {
      if (v[a].ms3[b].prob > 0.8) lowE++;
    }

    if (lowE < 1) continue;

    double mm;
    if(v[a].monoMZ>0) mm=v[a].monoMZ*v[a].charge-1.007276466*v[a].charge;
    else mm= v[a].mz * v[a].charge - 1.007276466 * v[a].charge;

    fprintf(f,"%d\t%d\t%.4lf\t%d\t%.4lf\t%d\n",(int)a, v[a].scan, v[a].mz, v[a].charge,mm,(int)v[a].ms3.size());

    vector<sXLPep> xl;
    for(size_t b=0;b<v[a].ms3.size();b++) {
      if(v[a].ms3[b].prob<0.8) continue;
      if(v[a].ms3[b].peptide.find("K[")!=string::npos){
        sXLPep p;
        p.index=b;
        p.mass=v[a].ms3[b].pepMass;
        xl.push_back(p);
      }
    }

    if(xl.size()<1) continue;
    sort(xl.begin(),xl.end(),compareXL);

    double pm=0;
    int res=0;
    for(size_t b=0;b<xl.size();b++){
      pm=xl[b].mass;
      double dif=mm-pm;
      if((dif>89.9 && dif<90.2) || (dif>121.9 && dif<122.2)) {
        fprintf(f,"  DE\t%s\t%s\t%.4lf\t%.4lf\t\t%.4lf\n",v[a].ms3[xl[b].index].peptide.c_str(), v[a].ms3[xl[b].index].protein.c_str(), v[a].ms3[xl[b].index].pepMass, v[a].ms3[xl[b].index].prob, mm - pm);
        res++;
      } else {
        for(size_t c=b+1;c<xl.size();c++){
          double xm=pm+xl[c].mass;
          double dif=mm-xm;
          if((dif>16.9 && dif<19.2) || (dif>48.9 && dif<51.2)) {
            fprintf(f, "  XL\t%s + %s\t%s & %s\t%.2lf + %.2lf\t%.2lf + %.2lf\t\t%.4lf\n", v[a].ms3[xl[b].index].peptide.c_str(), v[a].ms3[xl[c].index].peptide.c_str(), v[a].ms3[xl[b].index].protein.c_str(), v[a].ms3[xl[c].index].protein.c_str(),v[a].ms3[xl[b].index].pepMass, v[a].ms3[xl[c].index].pepMass,v[a].ms3[xl[b].index].prob, v[a].ms3[xl[c].index].prob, mm - xm);
            res++;
          }
        }
      }
    }
    if(res==0){
      for (size_t b = 0; b < v[a].ms3.size(); b++) {
        if (v[a].ms3[b].prob < 0.8) continue;
        if (v[a].ms3[b].peptide.find("K[") != string::npos) {
          fprintf(f, "  U\t%s\t%s\t%.4lf\t%.4lf\t\t%.4lf\n", v[a].ms3[b].peptide.c_str(), v[a].ms3[b].protein.c_str(), v[a].ms3[b].pepMass, v[a].ms3[b].prob, mm - v[a].ms3[b].pepMass);
        }
      }
    }


    //fprintf(f, "\nGroup %d, Scan: %d\tm/z: %.4lf\tz: %d\n", (int)a, v[a].scan, v[a].mz, v[a].charge);
    //for (size_t b = 0; b < v[a].ms3.size(); b++) {
    //  if (v[a].ms3[b].prob < 0.8) continue;
    //  fprintf(f, "MS3 Scan: %d\tm/z:%.4lf\tz:%d\n", v[a].ms3[b].scan, v[a].ms3[b].mz, v[a].ms3[b].charge);
    //  if (!v[a].ms3[b].peptide.empty()) {
    //    fprintf(f, "  %s\t%s\te-value:%e\txcorr:%.4lf\tprobability:%.4lf\n", v[a].ms3[b].peptide.c_str(), v[a].ms3[b].protein.c_str(), v[a].ms3[b].eval, v[a].ms3[b].xcorr, v[a].ms3[b].prob);
    //  }
    //}
  }
}

void marquee() {
  cout << "XL3: Some sort of cleavable crosslinker analysis program" << endl;
  cout << "Copyright Michael Hoopmann, Institute for Systems Biology" << endl;
  cout << "Version/Build: 2021.06.04" << endl;
}

void usage(){
  cout << "\nUSAGE: XL3.exe <mzML> <PepXML> <outfile>" << endl;
  cout << "  mzML = mzML file of the spectra analyzed." << endl;
  cout << "  PepXML = PepXML of the MS3 scan database search results with PeptideProphet analysis." << endl;
  cout << "  outfile = desired output file name." << endl;
}


bool compareXL(sXLPep& a, sXLPep& b){
  return a.mass>b.mass;
}