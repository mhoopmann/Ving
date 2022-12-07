#include "VingData.h"

using namespace MSToolkit;
using namespace std;
//using namespace mzParser;

VingData::VingData(VingParameters* p){
  params=p;
}

VingData::~VingData(){
  params=NULL;
}

void VingData::assessXLType() {
  for (size_t a = 0; a < groups.size(); a++) {

    //Single peptides are defined by high probability after MS2 searching.
    if (groups[a].probabilityMS2 > 0.8) {
      groups[a].type = xlSingle;
      groups[a].sequence = groups[a].peptide;
      groups[a].proteinS = groups[a].protein;
      for (size_t b = 0; b < groups[a].mods.size(); b++) {
        if (params->crosslinker.targetA.find(groups[a].mods[b].aa) != string::npos) {
          if (fabs(groups[a].mods[b].mass - (params->crosslinker.xlMass + params->crosslinker.capMassA)) < 0.01) {
            groups[a].type = xlDeadEnd; //generally make deadends known here
            break;
          }
        }
        if (params->crosslinker.targetB.find(groups[a].mods[b].aa) != string::npos) {
          if (fabs(groups[a].mods[b].mass - (params->crosslinker.xlMass + params->crosslinker.capMassB)) < 0.01) {
            groups[a].type = xlDeadEnd; //generally make deadends known here
            break;
          }
        }
      }
      continue;
    }

    //If MS2 group has no MS3 peptides of good probability, leave type as Unknown.
    int lowE = 0;
    for (size_t b = 0; b < groups[a].ms3.size(); b++) {
      if (groups[a].ms3[b].prob > 0.8) lowE++;
    }
    if (lowE < 1) continue;

    //Look for XL-related results
    //First compute the neutral precursor mass. If a monoisotopic peak was not provided, then
    //estimate the neutral precursor mass from the selected ion peak.
    double mm;
    if (groups[a].monoMZ > 0) mm = groups[a].monoMZ * groups[a].charge - 1.007276466 * groups[a].charge;
    else mm = groups[a].mz * groups[a].charge - 1.007276466 * groups[a].charge;

    vector<sXLPep> xl; //sort through all combinations of MS3 peptides with crosslink mods
    for (size_t b = 0; b < groups[a].ms3.size(); b++) {
      if (groups[a].ms3[b].prob < 0.8) continue;
      if (groups[a].ms3[b].stubMass != 0) {
        //if (groups[a].ms3[b].peptide.find("H[") != string::npos) {
        sXLPep p;
        p.index = b;
        p.mass = groups[a].ms3[b].pepMass - groups[a].ms3[b].stubMass;
        p.stub = groups[a].ms3[b].stubMass;
        xl.push_back(p);
      }
    }
    //no crosslink mods, then not a crosslink that we can see....
    //TODO: try to fall back on alternate means, such as looking for reporter ions, etc.
    if (xl.size() < 1) continue;
    sort(xl.begin(), xl.end(), compareXL);

    double pm = 0;
    int res = 0;
    for (size_t b = 0; b < xl.size(); b++) {
      pm = xl[b].mass;
      double dif;
      if(groups[a].ms3[xl[b].index].stubA) dif = mm - (pm+params->crosslinker.xlMass+params->crosslinker.capMassA);
      else dif = mm - (pm + params->crosslinker.xlMass + params->crosslinker.capMassB);

      //change this to check for the remaining mass of a dead-end
      bool bDE = false;
      if (fabs(dif) < 0.1) bDE = true;

      //Dead-end has a mass difference of a long or short arm.
      if (bDE) {
        groups[a].type = xlDeadEnd;
        groups[a].sequence = groups[a].ms3[xl[b].index].peptide;
        groups[a].proteinS = groups[a].ms3[xl[b].index].protein;
        res++;
        break;

        //XL has two peptides plus the spacer arm
      } else {
        for (size_t c = b + 1; c < xl.size(); c++) {
          double xm = xl[b].mass + xl[c].mass + params->crosslinker.xlMass;
          double dif = mm - xm;

          //if ((dif > 16.9 && dif < 19.2) || (dif > 48.9 && dif < 51.2)) {
          if (fabs(dif) < 0.1) {
            groups[a].type = xlXL;
            groups[a].sequence = groups[a].ms3[xl[b].index].peptide + "+" + groups[a].ms3[xl[c].index].peptide;
            groups[a].proteinS = groups[a].ms3[xl[b].index].protein + "+" + groups[a].ms3[xl[c].index].protein;
            res++;
            break;
          }
        }
      }
    }

    //check for isotope error
    //TODO: Mark isotope error somehow
    if (res == 0) {
      for (size_t b = 0; b < xl.size(); b++) {
        pm = xl[b].mass;
        double dif;
        if (groups[a].ms3[xl[b].index].stubA) dif = mm - (pm + params->crosslinker.xlMass + params->crosslinker.capMassA) - 1.003354835;
        else dif = mm - (pm + params->crosslinker.xlMass + params->crosslinker.capMassB) - 1.003354835;

        //change this to check for the remaining mass of a dead-end
        bool bDE = false;
        if (fabs(dif) < 0.1) bDE = true;

        //Dead-end has a mass difference of a long or short arm.
        if (bDE) {
          groups[a].type = xlDeadEnd;
          groups[a].sequence = groups[a].ms3[xl[b].index].peptide;
          groups[a].proteinS = groups[a].ms3[xl[b].index].protein;
          res++;
          break;

          //XL has two peptides plus the spacer arm
        } else {
          for (size_t c = b + 1; c < xl.size(); c++) {
            double xm = xl[b].mass + xl[c].mass + params->crosslinker.xlMass;
            double dif = mm - xm - 1.003354835;

            //if ((dif > 16.9 && dif < 19.2) || (dif > 48.9 && dif < 51.2)) {
            if (fabs(dif) < 0.1) {
              groups[a].type = xlXL;
              groups[a].sequence = groups[a].ms3[xl[b].index].peptide + "+" + groups[a].ms3[xl[c].index].peptide;
              groups[a].proteinS = groups[a].ms3[xl[b].index].protein + "+" + groups[a].ms3[xl[c].index].protein;
              res++;
              break;
            }
          }
        }
      }
    }

    if (res == 0) {
      for (size_t b = 0; b < groups[a].ms3.size(); b++) {
        if (groups[a].ms3[b].prob < 0.8) continue;

        //Cross-link is incomplete if there is evidence of a stub, but nothing adds up.
        if (groups[a].ms3[b].stubMass != 0) {
          //if (groups[a].ms3[b].peptide.find("H[") != string::npos) {
          groups[a].type = xlIncomplete;
          break;
        }
      }
    }

  }

  int uCount = 0;
  int iCount = 0;
  int xCount = 0;
  int dCount = 0;
  int sCount = 0;
  for (size_t a = 0; a < groups.size(); a++) {
    switch (groups[a].type) {
    case xlUnknown: uCount++; break;
    case xlIncomplete: iCount++; break;
    case xlSingle: sCount++; break;
    case xlDeadEnd: dCount++; break;
    case xlXL: xCount++; break;
    default: break;
    }
  }

  cout << "Groups without any results: " << uCount << "  (" << (double)uCount / groups.size() * 100 << "%)" << endl;
  cout << "Groups with single (unlinked) peptides: " << sCount << "  (" << (double)sCount / groups.size() * 100 << "%)" << endl;
  cout << "Groups with incomplete results: " << iCount << "  (" << (double)iCount / groups.size() * 100 << "%)" << endl;
  cout << "Groups with deadend peptides: " << dCount << "  (" << (double)dCount / groups.size() * 100 << "%)" << endl;
  cout << "Groups with crosslinked peptides: " << xCount << "  (" << (double)xCount / groups.size() * 100 << "%)" << endl;

}

//Reads in PepXML from MS2 search results and assigns information to sMS2 objects
bool VingData::importMS2SearchResults() {
  NeoPepXMLParser p;
  p.read(params->ms2SearchResult.c_str());
  cout << "Done reading" << endl;

  size_t vp = 0;

  for (size_t a = 0; a < p.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query.size(); a++) {
    CnpxSpectrumQuery* sq = &p.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query[a];

    char* tok;
    char str[256];
    strcpy(str, sq->spectrum.c_str());
    tok = strtok(str, ".\n\r");
    tok = strtok(NULL, ".\n\r");
    int scan = atoi(tok);

    while (vp < groups.size() && groups[vp].scan < scan) vp++;
    if (groups[vp].scan != scan) {
      cout << "No match to MS2: " << scan << endl;
      continue;
    }

    if (sq->search_result[0].search_hit.size() == 0) continue;

    CnpxSearchHit* sh = &sq->search_result[0].search_hit[0];
    if (sh->modification_info.empty()) groups[vp].peptide = sh->peptide;
    else groups[vp].peptide = sh->modification_info[0].modified_peptide;

    groups[vp].protein = sh->protein;
    groups[vp].charge = sq->assumed_charge;
    groups[vp].calcNeutMass = sh->calc_neutral_pep_mass;

    /*for (size_t b = 0; b < sh->search_score.size(); b++) {
      CnpxSearchScore* ss = &sh->search_score[b];
      if (ss->name.compare("xcorr") == 0) groups[vp].ms3[v3].xcorr = atof(ss->value.c_str());
      if (ss->name.compare("expect") == 0) groups[vp].ms3[v3].eval = atof(ss->value.c_str());
    }*/

    for (size_t b = 0; b < sh->analysis_result.size(); b++) {
      if (sh->analysis_result[b].analysis.compare("peptideprophet") == 0) {
        groups[vp].probabilityMS2 = sh->analysis_result[b].peptide_prophet_result.probability;
      }
    }

    //parse mods
    if (!sh->modification_info.empty()) {
      for (size_t b = 0; b < sh->modification_info[0].mod_aminoacid_mass.size(); b++) {
        CnpxModAminoAcidMass* mam = &sh->modification_info[0].mod_aminoacid_mass[b];
        sMod mod;
        mod.aa = sh->peptide[mam->position - 1];
        mod.pos = mam->position;
        mod.mass = mam->variable;
        groups[vp].mods.push_back(mod);
      }
    }

  }
  return true;
}

bool VingData::importMS3SearchResults() {
  NeoPepXMLParser p;
  if(!p.read(params->ms3SearchResult.c_str())) return false;
  cout << "Done reading" << endl;

  size_t vp = 0;
  size_t v3 = 0;

  while (vp < groups.size() && groups[vp].ms3.size() == 0) vp++;

  for (size_t a = 0; a < p.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query.size(); a++) {
    CnpxSpectrumQuery* sq = &p.msms_pipeline_analysis[0].msms_run_summary[0].spectrum_query[a];

    char* tok;
    char str[256];
    strcpy(str, sq->spectrum.c_str());
    tok = strtok(str, ".\n\r");
    tok = strtok(NULL, ".\n\r");
    int scan = atoi(tok);

    while (vp < groups.size() && groups[vp].ms3[v3].scan < scan) {
      v3++;
      if (v3 == groups[vp].ms3.size()) {
        v3 = 0;
        vp++;
        while (vp < groups.size() && groups[vp].ms3.size() == 0) vp++;
      }
    }

    if (groups[vp].ms3[v3].scan != scan) {
      cout << "No match to " << scan << endl;
      continue;
    }

    if (sq->search_result[0].search_hit.size() == 0) continue;

    CnpxSearchHit* sh = &sq->search_result[0].search_hit[0];
    if (sh->modification_info.empty()) groups[vp].ms3[v3].peptide = sh->peptide;
    else groups[vp].ms3[v3].peptide = sh->modification_info[0].modified_peptide;

    groups[vp].ms3[v3].protein = sh->protein;
    groups[vp].ms3[v3].charge = sq->assumed_charge;
    groups[vp].ms3[v3].pepMass = sh->calc_neutral_pep_mass;

    for (size_t b = 0; b < sh->search_score.size(); b++) {
      CnpxSearchScore* ss = &sh->search_score[b];
      if (ss->name.compare("xcorr") == 0) groups[vp].ms3[v3].xcorr = atof(ss->value.c_str());
      if (ss->name.compare("expect") == 0) groups[vp].ms3[v3].eval = atof(ss->value.c_str());
    }

    for (size_t b = 0; b < sh->analysis_result.size(); b++) {
      if (sh->analysis_result[b].analysis.compare("peptideprophet") == 0) {
        groups[vp].ms3[v3].prob = sh->analysis_result[b].peptide_prophet_result.probability;
      }
    }

    //check for stubs
    if (!sh->modification_info.empty()) {
      for (size_t b = 0; b < sh->modification_info[0].mod_aminoacid_mass.size(); b++) {
        CnpxModAminoAcidMass* mam = &sh->modification_info[0].mod_aminoacid_mass[b];
        char aa = sh->peptide[mam->position - 1];
        if (params->crosslinker.targetA.find(aa) != string::npos) {
          size_t c;
          for (c = 0; c < params->crosslinker.stubA.size(); c++) {
            if (fabs(params->crosslinker.stubA[c].mass - mam->variable) < 0.1) {
              groups[vp].ms3[v3].stubMass = params->crosslinker.stubA[c].mass;
              groups[vp].ms3[v3].stubA = true;
              break;
            } else if (fabs(params->crosslinker.stubA[c].mass + params->crosslinker.stubA[c].reaction - mam->variable) < 0.1) {
              groups[vp].ms3[v3].stubMass = params->crosslinker.stubA[c].mass + params->crosslinker.stubA[c].reaction;
              groups[vp].ms3[v3].stubA = true;
              break;
            }
          }
          if (c < params->crosslinker.stubA.size()) break;
        }
        if (params->crosslinker.targetB.find(aa) != string::npos) {
          size_t c;
          for (c = 0; c < params->crosslinker.stubB.size(); c++) {
            if (fabs(params->crosslinker.stubB[c].mass - mam->variable) < 0.1) {
              groups[vp].ms3[v3].stubMass = params->crosslinker.stubB[c].mass;
              groups[vp].ms3[v3].stubA = false;
              break;
            } else if (fabs(params->crosslinker.stubB[c].mass + params->crosslinker.stubB[c].reaction - mam->variable) < 0.1) {
              groups[vp].ms3[v3].stubMass = params->crosslinker.stubB[c].mass + params->crosslinker.stubB[c].reaction;
              groups[vp].ms3[v3].stubA = false;
              break;
            }
          }
          if (c < params->crosslinker.stubB.size()) break;
        }
      }
    }
  }

}

bool VingData::parseMzML(){
  MSReader r;
  Spectrum s;

  r.setFilter(MS2);
  r.addFilter(MS3);
  r.addFilter(MS1);

  vector<sMS2> cycle;
  bool bMS2 = false;

  groups.clear();

  if(!r.readFile(params->mzML.c_str(), s)) return false;
  while (s.getScanNumber() > 0) {

    if (s.getMsLevel() == 1) {
      for (size_t a = 0; a < cycle.size(); a++) {
        if (bMS2) groups.push_back(cycle[a]);
      }
      cycle.clear();
      bMS2 = false;

    } else if (s.getMsLevel() == 2) {
      sMS2 ms2;
      ms2.scan = s.getScanNumber();
      ms2.charge = s.getCharge();
      ms2.mz = s.getMZ();
      ms2.monoMZ = s.getMonoMZ();
      cycle.push_back(ms2);
      bMS2 = true;
    } else {
      sMS3 ms3;
      ms3.scan = s.getScanNumber();
      ms3.charge = s.getCharge();
      ms3.mz = s.getMZ();
      double pMZ = 0;
      char header[256];
      s.getRawFilter(header, 256);
      char* tok = strtok(header, " ");
      while (tok != NULL) {
        string st = tok;
        size_t pos = st.find_first_of('@');
        if (pos != string::npos) {
          st = st.substr(0, pos);
          pMZ = atof(st.c_str());
          break;
        }
        tok = strtok(NULL, " ");
      }
      if (pMZ == 0) cout << "Cannot extract MS2 selected ion" << endl;
      size_t a = 0;
      for (a = 0; a < cycle.size(); a++) {
        if (fabs(cycle[a].mz - pMZ) < 0.01) break;
      }
      if (a == cycle.size()) cout << "Cannot find parent MS2 scan." << endl;
      else cycle[a].ms3.push_back(ms3);
    }

    r.readFile(NULL, s);
  }
  for (size_t a = 0; a < cycle.size(); a++) {
    if (bMS2) groups.push_back(cycle[a]);
  }

  cout << groups.size() << " total groups." << endl;
}

bool VingData::compareXL(sXLPep& a, sXLPep& b) {
  return a.mass > b.mass;
}
