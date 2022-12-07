#include "MSReader.h"
#include "NeoPepXMLParser.h"
#include <algorithm>

#include "VingData.h"
//#include "vingStructs.h"
#include "VingParameters.h"

using namespace std;
using namespace MSToolkit;

#define BDATE "Dec 7 2022"
#define VERSION "0.7.2"

//enum eXLType{
//  xlUnknown=0,
//  xlSingle,
//  xlDeadEnd,
//  xlLoop,
//  xlXL,
//  xlIncomplete,
//  xlSize
//};
//
//typedef struct sMS3{
//  int scan;
//  string peptide;
//  string protein;
//  double xcorr;
//  double eval;
//  double mz;
//  double prob;
//  double pepMass;   //mass including stub
//  double stubMass;  //just the stub
//  bool stubA;       //true if reacting to target A, false if reacting to target B.
//  int charge;
//  sMS3(){
//    scan=0;
//    xcorr=0;
//    eval=0;
//    mz=0;
//    charge=0;
//    pepMass=0;
//    stubMass=0;
//    stubA=true;
//  }
//} sMS3;
//
//typedef struct sMod{
//  char aa;
//  int pos;
//  double mass;
//}sMod;
//
//
//typedef struct sMS2{
//  int scan;
//  double mz;
//  int charge;
//  double monoMZ;
//  double calcNeutMass;
//  double probability;
//  double probabilityMS2;
//  string peptide;
//  string sequence;
//  string protein;
//  string proteinS;
//  eXLType type;
//  vector<sMS3> ms3;
//  vector<sMod> mods;
//  sMS2(){ clear();}
//  void clear(){
//    scan=0;
//    mz=0;
//    charge=0;
//    monoMZ=0;
//    calcNeutMass=0;
//    probability=0;
//    probabilityMS2=0;
//    type=xlUnknown;
//    ms3.clear();
//    mods.clear();
//  }
//} sMS2;
//
//typedef struct sXLPep{
//  size_t index;
//  double mass;
//  double stub;
//} sXLPep;

//void addXML2(const char* fn, vector<sMS2>& v);
//void addXML3(const char* fn, vector<sMS2>& v);
//void assessXLType(vector<sMS2>& v);
bool cmdLine(int argc, char* argv[]);
//void exportResults(const char* fn, vector<sMS2>& v);
void exportResults2(const char* fn, vector<sMS2>& v);
void marquee();
void usage();

//bool compareXL(sXLPep& a, sXLPep& b);

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
  data.parseMzML();
  data.importMS2SearchResults();
  data.importMS3SearchResults();
  data.assessXLType();

  //MSReader r;
  //Spectrum s;

  //r.setFilter(MS2);
  //r.addFilter(MS3);
  //r.addFilter(MS1);

  //vector<sMS2> groups;
  //vector<sMS2> cycle;
  //bool bMS2=false;
  //

  //r.readFile(params.mzML.c_str(),s);
  //while(s.getScanNumber()>0){

  //  if(s.getMsLevel()==1){
  //    for(size_t a=0;a<cycle.size();a++){
  //      if(bMS2) groups.push_back(cycle[a]);
  //    }
  //    cycle.clear();
  //    bMS2=false;

  //  } else if(s.getMsLevel()==2){
  //    sMS2 ms2;
  //    ms2.scan=s.getScanNumber();
  //    ms2.charge=s.getCharge();
  //    ms2.mz=s.getMZ();
  //    ms2.monoMZ=s.getMonoMZ();
  //    cycle.push_back(ms2);
  //    bMS2=true;
  //  } else {
  //    sMS3 ms3;
  //    ms3.scan=s.getScanNumber();
  //    ms3.charge=s.getCharge();
  //    ms3.mz=s.getMZ();
  //    double pMZ=0;
  //    char header[256];
  //    s.getRawFilter(header,256);
  //    char* tok=strtok(header," ");
  //    while(tok!=NULL){
  //      string st=tok;
  //      size_t pos=st.find_first_of('@');
  //      if(pos!=string::npos){
  //        st=st.substr(0,pos);
  //        pMZ=atof(st.c_str());
  //        break;
  //      }
  //      tok = strtok(NULL, " ");
  //    }
  //    if(pMZ==0) cout << "Cannot extract MS2 selected ion" << endl;
  //    size_t a=0;
  //    for(a=0;a<cycle.size();a++){
  //      if(fabs(cycle[a].mz-pMZ)<0.01) break;
  //    }
  //    if(a==cycle.size()) cout << "Cannot find parent MS2 scan." << endl;
  //    else cycle[a].ms3.push_back(ms3);
  //  }

  //  r.readFile(NULL,s);
  //}
  //for (size_t a = 0; a < cycle.size(); a++) {
  //  if (bMS2) groups.push_back(cycle[a]);
  //}

  //cout << groups.size() << " total groups." << endl;

  //addXML2(params.ms2SearchResult.c_str(),groups);
  //addXML3(params.ms3SearchResult.c_str(),groups);
  //assessXLType(groups);
  exportResults2(params.output.c_str(),data.groups);

  return 0;
}

bool cmdLine(int argc, char* argv[]){
  if(argc!=2) return false;
  g_param=argv[1];
  return true;
}

////Reads in PepXML from MS2 search results and assigns information to sMS2 objects
//void addXML2(const char* fn, vector<sMS2>& v) {
//  NeoPepXMLParser p;
//  p.read(fn);
//
//  cout << "Done reading" << endl;
//
//  size_t vp = 0;
//
//  for (size_t a = 0; a < p.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query.size(); a++) {
//    CnpxSpectrumQuery* sq = &p.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query[a];
//
//    char* tok;
//    char str[256];
//    strcpy(str, sq->spectrum.c_str());
//    tok = strtok(str, ".\n\r");
//    tok = strtok(NULL, ".\n\r");
//    int scan = atoi(tok);
//
//    while (vp < v.size() && v[vp].scan < scan) vp++;
//    if (v[vp].scan != scan) {
//      cout << "No match to MS2: " << scan << endl;
//      continue;
//    }
//
//    if (sq->search_result[0].search_hit.size() == 0) continue;
//
//    CnpxSearchHit* sh = &sq->search_result[0].search_hit[0];
//    if (sh->modification_info.empty()) v[vp].peptide = sh->peptide;
//    else v[vp].peptide = sh->modification_info[0].modified_peptide;
//
//    v[vp].protein = sh->protein;
//    v[vp].charge = sq->assumed_charge;
//    v[vp].calcNeutMass = sh->calc_neutral_pep_mass;
//
//    /*for (size_t b = 0; b < sh->search_score.size(); b++) {
//      CnpxSearchScore* ss = &sh->search_score[b];
//      if (ss->name.compare("xcorr") == 0) v[vp].ms3[v3].xcorr = atof(ss->value.c_str());
//      if (ss->name.compare("expect") == 0) v[vp].ms3[v3].eval = atof(ss->value.c_str());
//    }*/
//
//    for (size_t b = 0; b < sh->analysis_result.size(); b++) {
//      if (sh->analysis_result[b].analysis.compare("peptideprophet") == 0) {
//        v[vp].probabilityMS2 = sh->analysis_result[b].peptide_prophet_result.probability;
//      }
//    }
//
//    //parse mods
//    if (!sh->modification_info.empty()) {
//      for (size_t b = 0; b < sh->modification_info[0].mod_aminoacid_mass.size(); b++) {
//        CnpxModAminoAcidMass* mam = &sh->modification_info[0].mod_aminoacid_mass[b];
//        sMod mod;
//        mod.aa = sh->peptide[mam->position - 1];
//        mod.pos = mam->position;
//        mod.mass = mam->variable;
//        v[vp].mods.push_back(mod);
//      }
//    }
//
//  }
//
//}
//
//void addXML3(const char* fn, vector<sMS2>& v){
//
//  NeoPepXMLParser p;
//  p.read(fn);
//
//  cout << "Done reading" << endl;
//
//  size_t vp=0;
//  size_t v3=0;
//
//  while(vp < v.size() && v[vp].ms3.size()==0) vp++;
//
//  for(size_t a=0;a<p.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query.size();a++){
//    CnpxSpectrumQuery* sq=&p.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query[a];
//   
//    char* tok;
//    char str[256];
//    strcpy(str,sq->spectrum.c_str());
//    tok=strtok(str,".\n\r");
//    tok=strtok(NULL,".\n\r");
//    int scan=atoi(tok);
//
//    while(vp<v.size() && v[vp].ms3[v3].scan<scan){
//      v3++;
//      if(v3==v[vp].ms3.size()){
//        v3=0;
//        vp++;
//        while (vp < v.size() && v[vp].ms3.size() == 0) vp++;
//      }
//    }
//
//    if(v[vp].ms3[v3].scan!=scan){
//      cout << "No match to " << scan << endl;
//      continue;
//    }
//
//    if(sq->search_result[0].search_hit.size()==0) continue;
//
//    CnpxSearchHit* sh= &sq->search_result[0].search_hit[0];
//    if(sh->modification_info.empty()) v[vp].ms3[v3].peptide=sh->peptide;
//    else v[vp].ms3[v3].peptide=sh->modification_info[0].modified_peptide;
//
//    v[vp].ms3[v3].protein=sh->protein;
//    v[vp].ms3[v3].charge=sq->assumed_charge;
//    v[vp].ms3[v3].pepMass=sh->calc_neutral_pep_mass;
//
//    for(size_t b=0;b<sh->search_score.size();b++){
//      CnpxSearchScore* ss=&sh->search_score[b];
//      if(ss->name.compare("xcorr")==0) v[vp].ms3[v3].xcorr=atof(ss->value.c_str());
//      if (ss->name.compare("expect") == 0) v[vp].ms3[v3].eval = atof(ss->value.c_str());
//    }
//
//    for(size_t b=0;b<sh->analysis_result.size();b++){
//      if(sh->analysis_result[b].analysis.compare("peptideprophet")==0) {
//        v[vp].ms3[v3].prob=sh->analysis_result[b].peptide_prophet_result.probability;
//      }
//    }
//
//    //check for stubs
//    if(!sh->modification_info.empty()){
//      for(size_t b=0;b<sh->modification_info[0].mod_aminoacid_mass.size();b++){
//        CnpxModAminoAcidMass* mam=&sh->modification_info[0].mod_aminoacid_mass[b];
//        char aa=sh->peptide[mam->position-1];
//        if(params.crosslinker.targetA.find(aa)!=string::npos) {
//          size_t c;
//          for(c=0;c<params.crosslinker.stubA.size();c++){
//            if(fabs(params.crosslinker.stubA[c].mass-mam->variable)<0.1){
//              v[vp].ms3[v3].stubMass = params.crosslinker.stubA[c].mass;
//              v[vp].ms3[v3].stubA=true;
//              break;
//            } else if(fabs(params.crosslinker.stubA[c].mass + params.crosslinker.stubA[c].reaction - mam->variable) < 0.1){
//              v[vp].ms3[v3].stubMass=params.crosslinker.stubA[c].mass + params.crosslinker.stubA[c].reaction;
//              v[vp].ms3[v3].stubA = true;
//              break;
//            }
//          }
//          if(c<params.crosslinker.stubA.size())break;
//        } 
//        if (params.crosslinker.targetB.find(aa) != string::npos) {
//          size_t c;
//          for (c = 0; c < params.crosslinker.stubB.size(); c++) {
//            if (fabs(params.crosslinker.stubB[c].mass - mam->variable) < 0.1) {
//              v[vp].ms3[v3].stubMass = params.crosslinker.stubB[c].mass;
//              v[vp].ms3[v3].stubA = false;
//              break;
//            } else if (fabs(params.crosslinker.stubB[c].mass + params.crosslinker.stubB[c].reaction - mam->variable) < 0.1) {
//              v[vp].ms3[v3].stubMass = params.crosslinker.stubB[c].mass + params.crosslinker.stubB[c].reaction;
//              v[vp].ms3[v3].stubA = false;
//              break;
//            }
//          }
//          if (c < params.crosslinker.stubB.size())break;
//        }
//      }
//    }
//  }
//
//}
//
//void assessXLType(vector<sMS2>& v){
//
//  for (size_t a = 0; a < v.size(); a++) {
//
//    //Single peptides are defined by high probability after MS2 searching.
//    if(v[a].probabilityMS2>0.8) {
//      v[a].type=xlSingle;
//      v[a].sequence=v[a].peptide;
//      v[a].proteinS=v[a].protein;
//      for(size_t b=0;b<v[a].mods.size();b++){
//        if(params.crosslinker.targetA.find(v[a].mods[b].aa)!=string::npos){
//          if (fabs(v[a].mods[b].mass - (params.crosslinker.xlMass + params.crosslinker.capMassA)) < 0.01) {
//            v[a].type = xlDeadEnd; //generally make deadends known here
//            break;
//          }
//        }
//        if (params.crosslinker.targetB.find(v[a].mods[b].aa) != string::npos) {
//          if(fabs(v[a].mods[b].mass-(params.crosslinker.xlMass+params.crosslinker.capMassB))<0.01) {
//            v[a].type=xlDeadEnd; //generally make deadends known here
//            break;
//          }
//        }
//      }
//      continue;
//    }
//
//    //If MS2 group has no MS3 peptides of good probability, leave type as Unknown.
//    int lowE = 0;
//    for (size_t b = 0; b < v[a].ms3.size(); b++) {
//      if (v[a].ms3[b].prob > 0.8) lowE++;
//    }
//    if (lowE < 1) continue;
//    
//    //Look for XL-related results
//    //First compute the neutral precursor mass. If a monoisotopic peak was not provided, then
//    //estimate the neutral precursor mass from the selected ion peak.
//    double mm;
//    if (v[a].monoMZ > 0) mm = v[a].monoMZ * v[a].charge - 1.007276466 * v[a].charge;
//    else mm = v[a].mz * v[a].charge - 1.007276466 * v[a].charge;
//
//
//    vector<sXLPep> xl; //sort through all combinations of MS3 peptides with crosslink mods
//    for (size_t b = 0; b < v[a].ms3.size(); b++) {
//      if (v[a].ms3[b].prob < 0.8) continue;
//      if(v[a].ms3[b].stubMass!=0){
//      //if (v[a].ms3[b].peptide.find("H[") != string::npos) {
//        sXLPep p;
//        p.index = b;
//        p.mass = v[a].ms3[b].pepMass-v[a].ms3[b].stubMass;
//        p.stub=v[a].ms3[b].stubMass;
//        xl.push_back(p);
//      }
//    }
//    //no crosslink mods, then not a crosslink that we can see....
//    //TODO: try to fall back on alternate means, such as looking for reporter ions, etc.
//    if (xl.size() < 1) continue; 
//    sort(xl.begin(), xl.end(), compareXL);
//
//    double pm = 0;
//    int res = 0;
//    for (size_t b = 0; b < xl.size(); b++) {
//      pm = xl[b].mass;
//      double dif;
//      if(v[a].ms3[xl[b].index].stubA) dif = mm - (pm+params.crosslinker.xlMass+params.crosslinker.capMassA); //xl[b].stub);
//      else dif = mm - (pm + params.crosslinker.xlMass + params.crosslinker.capMassB);
//
//      //change this to check for the remaining mass of a dead-end
//      bool bDE=false;
//      if(fabs(dif)<0.1) bDE=true;
//      //if(a==1634) cout << b << "," << xl[b].mass << "\t" << (pm + params.crosslinker.xlMass + xl[b].stub) << "\t" << dif << endl;
//
//      //Dead-end has a mass difference of a long or short arm.
//      if (bDE) {
//        v[a].type = xlDeadEnd;
//        v[a].sequence = v[a].ms3[xl[b].index].peptide;
//        v[a].proteinS = v[a].ms3[xl[b].index].protein;
//        res++;
//        break;
//
//        //XL has two peptides plus the spacer arm
//      } else {
//        for (size_t c = b + 1; c < xl.size(); c++) {
//          double xm = xl[b].mass + xl[c].mass + params.crosslinker.xlMass;
//          double dif = mm - xm;
//
//          //if(a==7356){
//          //  cout << b << "," <<xl[b].mass << "\t" << c << "," << xl[c].mass <<"\t" << xm  << "\t" << dif << endl;
//          //}
//
//          //if ((dif > 16.9 && dif < 19.2) || (dif > 48.9 && dif < 51.2)) {
//          if(fabs(dif)<0.1){
//            v[a].type = xlXL;
//            v[a].sequence= v[a].ms3[xl[b].index].peptide + "+" + v[a].ms3[xl[c].index].peptide;
//            v[a].proteinS = v[a].ms3[xl[b].index].protein + "+" + v[a].ms3[xl[c].index].protein;
//            res++;
//            break;
//          }
//        }
//      }
//    }
//
//    //check for isotope error
//    //TODO: Mark isotope error somehow
//    if(res==0){
//      for (size_t b = 0; b < xl.size(); b++) {
//        pm = xl[b].mass;
//        double dif;
//        if (v[a].ms3[xl[b].index].stubA) dif = mm - (pm + params.crosslinker.xlMass + params.crosslinker.capMassA) - 1.003354835;
//        else dif = mm - (pm + params.crosslinker.xlMass + params.crosslinker.capMassB) - 1.003354835;
//
//        //change this to check for the remaining mass of a dead-end
//        bool bDE = false;
//        if (fabs(dif) < 0.1) bDE = true;
//
//        //Dead-end has a mass difference of a long or short arm.
//        if (bDE) {
//          v[a].type = xlDeadEnd;
//          v[a].sequence = v[a].ms3[xl[b].index].peptide;
//          v[a].proteinS = v[a].ms3[xl[b].index].protein;
//          res++;
//          break;
//
//          //XL has two peptides plus the spacer arm
//        } else {
//          for (size_t c = b + 1; c < xl.size(); c++) {
//            double xm = xl[b].mass + xl[c].mass + params.crosslinker.xlMass;
//            double dif = mm - xm - 1.003354835;
//
//            //if (a == 7356) {
//            //  cout << "offset: " << b << "," << xl[b].mass << "\t" << c << "," << xl[c].mass << "\t" << xm << "\t" << dif << endl;
//            //}
//
//            //if ((dif > 16.9 && dif < 19.2) || (dif > 48.9 && dif < 51.2)) {
//            if (fabs(dif) < 0.1) {
//              v[a].type = xlXL;
//              v[a].sequence = v[a].ms3[xl[b].index].peptide + "+" + v[a].ms3[xl[c].index].peptide;
//              v[a].proteinS = v[a].ms3[xl[b].index].protein + "+" + v[a].ms3[xl[c].index].protein;
//              res++;
//              break;
//            }
//          }
//        }
//      }
//    }
//
//    if (res == 0) {
//      for (size_t b = 0; b < v[a].ms3.size(); b++) {
//        if (v[a].ms3[b].prob < 0.8) continue;
//
//        //Cross-link is incomplete if there is evidence of a stub, but nothing adds up.
//        if (v[a].ms3[b].stubMass!=0){
//        //if (v[a].ms3[b].peptide.find("H[") != string::npos) {
//          v[a].type = xlIncomplete;
//          break;
//        }
//      }
//    }
//
//  }
//
//  int uCount=0;
//  int iCount=0;
//  int xCount=0;
//  int dCount=0;
//  int sCount=0;
//  for(size_t a=0;a<v.size();a++){
//    switch(v[a].type){
//    case xlUnknown: uCount++; break;
//    case xlIncomplete: iCount++; break;
//    case xlSingle: sCount++; break;
//    case xlDeadEnd: dCount++; break;
//    case xlXL: xCount++; break;
//    default: break;
//    }
//  }
//
//  cout << "Groups without any results: " << uCount << "  (" << (double)uCount/v.size()*100 << "%)" << endl;
//  cout << "Groups with single (unlinked) peptides: " << sCount << "  (" << (double)sCount / v.size() * 100 << "%)" << endl;
//  cout << "Groups with incomplete results: " << iCount << "  (" << (double)iCount / v.size() * 100 << "%)" << endl;
//  cout << "Groups with deadend peptides: " << dCount << "  (" << (double)dCount / v.size() * 100 << "%)" << endl;
//  cout << "Groups with crosslinked peptides: " << xCount << "  (" << (double)xCount / v.size() * 100 << "%)" << endl;
//
//}

//void exportResults(const char* fn, vector<sMS2>& v){
//  FILE* f = fopen(fn, "wt");
//
//  fprintf(f,"GroupID\tMS2ScanNum\tMS2mz\tMS2charge\tMS2monoMass\tMS3ScanCount\n");
//
//  for (size_t a = 0; a < v.size(); a++) {
//    int lowE = 0;
//    for (size_t b = 0; b < v[a].ms3.size(); b++) {
//      if (v[a].ms3[b].prob > 0.8) lowE++;
//    }
//
//    if (lowE < 1) continue;
//
//    double mm;
//    if(v[a].monoMZ>0) mm=v[a].monoMZ*v[a].charge-1.007276466*v[a].charge;
//    else mm= v[a].mz * v[a].charge - 1.007276466 * v[a].charge;
//
//    fprintf(f,"%d\t%d\t%.4lf\t%d\t%.4lf\t%d\n",(int)a, v[a].scan, v[a].mz, v[a].charge,mm,(int)v[a].ms3.size());
//
//    vector<sXLPep> xl;
//    for(size_t b=0;b<v[a].ms3.size();b++) {
//      if(v[a].ms3[b].prob<0.8) continue;
//      if(v[a].ms3[b].peptide.find("K[")!=string::npos){
//        sXLPep p;
//        p.index=b;
//        p.mass=v[a].ms3[b].pepMass;
//        xl.push_back(p);
//      }
//    }
//
//    if(xl.size()<1) continue;
//    sort(xl.begin(),xl.end(),compareXL);
//
//    double pm = 0;
//    int res = 0;
//    for (size_t b = 0; b < xl.size(); b++) {
//      pm = xl[b].mass;
//      double dif = mm - pm;
//      if ((dif > 89.9 && dif < 90.2) || (dif > 121.9 && dif < 122.2)) {
//        fprintf(f, "  DE\t%s\t%s\t%.4lf\t%.4lf\t\t%.4lf\n", v[a].ms3[xl[b].index].peptide.c_str(), v[a].ms3[xl[b].index].protein.c_str(), v[a].ms3[xl[b].index].pepMass, v[a].ms3[xl[b].index].prob, mm - pm);
//        res++;
//      } else {
//        for (size_t c = b + 1; c < xl.size(); c++) {
//          double xm = pm + xl[c].mass;
//          double dif = mm - xm;
//          if ((dif > 16.9 && dif < 19.2) || (dif > 48.9 && dif < 51.2)) {
//            fprintf(f, "  XL\t%s + %s\t%s & %s\t%.2lf + %.2lf\t%.2lf + %.2lf\t\t%.4lf\n", v[a].ms3[xl[b].index].peptide.c_str(), v[a].ms3[xl[c].index].peptide.c_str(), v[a].ms3[xl[b].index].protein.c_str(), v[a].ms3[xl[c].index].protein.c_str(), v[a].ms3[xl[b].index].pepMass, v[a].ms3[xl[c].index].pepMass, v[a].ms3[xl[b].index].prob, v[a].ms3[xl[c].index].prob, mm - xm);
//            res++;
//          }
//        }
//      }
//    }
//
//
//    //fprintf(f, "\nGroup %d, Scan: %d\tm/z: %.4lf\tz: %d\n", (int)a, v[a].scan, v[a].mz, v[a].charge);
//    //for (size_t b = 0; b < v[a].ms3.size(); b++) {
//    //  if (v[a].ms3[b].prob < 0.8) continue;
//    //  fprintf(f, "MS3 Scan: %d\tm/z:%.4lf\tz:%d\n", v[a].ms3[b].scan, v[a].ms3[b].mz, v[a].ms3[b].charge);
//    //  if (!v[a].ms3[b].peptide.empty()) {
//    //    fprintf(f, "  %s\t%s\te-value:%e\txcorr:%.4lf\tprobability:%.4lf\n", v[a].ms3[b].peptide.c_str(), v[a].ms3[b].protein.c_str(), v[a].ms3[b].eval, v[a].ms3[b].xcorr, v[a].ms3[b].prob);
//    //  }
//    //}
//  }
//}

void exportResults2(const char* fn, vector<sMS2>& v){
  FILE* f = fopen(fn, "wt");

  //Massive Header
  fprintf(f, "GroupID\tDesignation\tMS2ScanNum\tMS2mz\tMS2charge\tMS2monoMass\tSequence\tProtein(s)\tMS2Peptide\tMS2Protein\tMS2CalcNeutMass\tMS2Prob");
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
    fprintf(f, "\t%d", v[a].scan);
    fprintf(f, "\t%.4lf", v[a].mz);
    fprintf(f, "\t%d", v[a].charge);
    fprintf(f, "\t%.4lf", v[a].monoMZ*v[a].charge-v[a].charge*1.007276466);
    fprintf(f, "\t%s", v[a].sequence.c_str());
    fprintf(f, "\t%s", v[a].proteinS.c_str());
    if(v[a].type==xlSingle){
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


//bool compareXL(sXLPep& a, sXLPep& b){
//  return a.mass>b.mass;
//}