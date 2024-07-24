#include "VingData.h"

using namespace MSToolkit;
using namespace std;

VingData::VingData(VingParameters* p){
  params=p;
  maxMS3Count=0;
}

VingData::~VingData(){
  params=NULL;
}

void VingData::assessIncomplete(string massList, string protein){
  vector<sPepMass> v;
  char str[1024];
  char* tok;
  FILE* f=fopen(massList.c_str(),"rt");
  while(!feof(f)){
    if(fgets(str,1024,f)==NULL) continue;
    if(strlen(str)<2) continue;
    sPepMass pm;
    tok=strtok(str," \t\n\r");
    pm.peptide=tok;
    tok=strtok(NULL, " \t\n\r");
    pm.z=atoi(tok);
    tok = strtok(NULL, " \t\n\r");
    pm.monoMass=atof(tok);
    tok = strtok(NULL, " \t\n\r");
    while(tok!=NULL){
      pm.mz=atof(tok);
      v.push_back(pm);
      tok = strtok(NULL, " \t\n\r");
    }
  }
  fclose(f);
  sort(v.begin(),v.end(),compareMZ);

  int iCount=0;
  int mCount=0;
  for(size_t a=0;a<groups.size();a++){
    if(groups[a].type!=xlIncomplete) continue;
    if(groups[a].probability<0.9) continue;
    if(!protein.empty()){  //skip groups that have no peptide IDs to our protein of interest
      size_t b;
      for(b=0;b<groups[a].ms3.size();b++){
        if(groups[a].ms3[b].protein[0].find(protein)!=string::npos) break;
      }
      if(b==groups[a].ms3.size()) continue;
    }

    //Output basic information
    iCount++;
    double mm = groups[a].monoMZ * groups[a].charge - groups[a].charge* PROTON;
    cout << files[groups[a].fileID].base << "\t" << groups[a].scan << "\t" << groups[a].monoMZ << "\t" << groups[a].charge << "\t" << mm << endl;

    //Read in spectrum
    MSReader r;
    Spectrum s;
    string file= files[groups[a].fileID].base+files[groups[a].fileID].ext;
    r.setFilter(MS2);
    r.readFile(file.c_str(),s,groups[a].scan);

    int sPos=0;
    size_t aPos=0;
    vector<sPepMass> match;
    while(sPos<s.size() && aPos<v.size()){
      double ppm=s[sPos].mz/1e5; //10 ppm tolerance
      while(aPos<v.size() && v[aPos].mz<s[sPos].mz-ppm) aPos++;
      if(aPos==v.size()) break;
      while(aPos<v.size() && v[aPos].mz<s[sPos].mz+ppm) {
        //cout << "  Peptide Match: " << v[aPos].peptide << "\t" << v[aPos].z << "\t" << v[aPos].mz << " vs. " << s[sPos].mz << endl;
        match.push_back(v[aPos]);
        match.back().matchMZ=s[sPos].mz;
        aPos++;
      }
      sPos++;
    }

    //Export any selected peaks
    for (size_t b = 0; b < groups[a].ms3.size(); b++) {
      cout << "  Selected Peak: " << groups[a].ms3[b].mz;
      if (!groups[a].ms3[b].peptide.empty()) cout << "\t" << groups[a].ms3[b].peptide << "\t" << groups[a].ms3[b].prob;
      cout << endl;
    }

    if(match.size()<2) continue;
    
    //get peak rank of matches
    s.sortIntensityRev();
    for(size_t b=0;b<match.size();b++){
      int c;
      for(c=0;c<s.size();c++){
        if(fabs(match[b].matchMZ-s[c].mz)<0.0001) {
          cout << "  Peptide Match: " << match[b].peptide << "\t" << match[b].z << "\t" << match[b].mz << " vs. " << match[b].matchMZ << "\t" << c+1 << endl;
          break;
        }
      }
      if(c==s.size()) cout << "  Error: cannot find matched peak." << endl;
    }

    //iterate matches for crosslink
    //cout << "  Checking matches" << endl;
    bool bMatch=false;
    for(size_t b=0;b<match.size()-1;b++){
      for(size_t c=b+1;c<match.size();c++){
        double dif=mm-match[b].monoMass-match[c].monoMass- params->crosslinker.xlMass;
        if(fabs(dif)<2){
          if(!bMatch){
            bMatch=true;
            mCount++;
          }
          cout << "    Candidate XL: " << match[b].peptide << " + " << match[c].peptide << "\t" << match[b].monoMass + match[c].monoMass + params->crosslinker.xlMass << "\t" << dif << endl;
        }
      }
    }
  }

  cout << mCount << " of " << iCount << " incomplete results are plausible crosslinks." << endl;

}

void VingData::assessXLType() {
  for (size_t a = 0; a < groups.size(); a++) {

    //First compute the neutral precursor mass. If a monoisotopic peak was not provided, then
    //estimate the neutral precursor mass from the selected ion peak.
    double mm;
    if (groups[a].monoMZ > 0) mm = groups[a].monoMZ * groups[a].charge - PROTON * groups[a].charge;
    else mm = groups[a].mz * groups[a].charge - PROTON * groups[a].charge;
    double ppm = mm / 1e6 * params->ppmTolerance;

    //Single peptides are defined by high probability after MS2 searching.
    if (groups[a].probabilityMS2 > groups[a].probability) {
      groups[a].type = xlSingle;
      groups[a].xlMass=0;
      groups[a].sequence = groups[a].peptide;
      groups[a].proteinS = processProteins(groups[a].protein);
      groups[a].probability= groups[a].probabilityMS2;
      groups[a].calcNeutMassG=groups[a].calcNeutMass;
      groups[a].ppmG=(mm-groups[a].calcNeutMass)/groups[a].calcNeutMass*1e6;
      for (size_t b = 0; b < groups[a].mods.size(); b++) {
        if (params->crosslinker.targetA.find(groups[a].mods[b].aa) != string::npos) {
          for(size_t c=0;c<params->crosslinker.capMassA.size();c++){
            if (fabs(groups[a].mods[b].mass - (params->crosslinker.xlMass + params->crosslinker.capMassA[c])) < ppm) {
              groups[a].type = xlDeadEnd; //generally make deadends known here
              groups[a].xlMass = params->crosslinker.xlMass + params->crosslinker.capMassA[c];
              break;
            }
          }
        }
        if (params->crosslinker.targetB.find(groups[a].mods[b].aa) != string::npos) {
          for (size_t c = 0; c < params->crosslinker.capMassB.size(); c++) {
              if (fabs(groups[a].mods[b].mass - (params->crosslinker.xlMass + params->crosslinker.capMassB[c])) < ppm) {
              groups[a].type = xlDeadEnd; //generally make deadends known here
              groups[a].xlMass = params->crosslinker.xlMass + params->crosslinker.capMassB[c];
              break;
            }
          }
        }
      }
    }

    //If MS2 group has no MS3 peptides of good probability, leave type as Unknown.
    int lowE = 0;
    for (size_t b = 0; b < groups[a].ms3.size(); b++) {
      if (groups[a].ms3[b].prob > 0.0) lowE++;
    }
    if (lowE < 1) continue;

    //Look for XL-related results
    vector<sXLPep> xl; //sort through all combinations of MS3 peptides with crosslink mods
    for (size_t b = 0; b < groups[a].ms3.size(); b++) {
      if (groups[a].ms3[b].stubMass != 0) {
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
    for(int io=-1;io<=1;io++){ //check all isotope errors
      for (size_t b = 0; b < xl.size(); b++) { //check all MS3 peptides
        pm = xl[b].mass;
        double nm;
        vector<double>* vc;
        if (groups[a].ms3[xl[b].index].stubA) vc= &params->crosslinker.capMassA;
        else vc=&params->crosslinker.capMassB;

        //TODO: Revise this, perhaps. Maybe check for XL first. Then alternately check for all possible dead-end possibilities.
        for(size_t d=0;d<vc->size();d++){
          nm=(pm+params->crosslinker.xlMass+vc->at(d));

          double dif = mm-nm - io * ISOTOPE_OFFSET;

          //change this to check for the remaining mass of a dead-end
          bool bDE = false;
          if (fabs(dif) < ppm) bDE = true;

          //if (a == 193) cout << a << "\t" << nm << "\t" << dif << "\t" << io << "\t" << bDE << "\tp:" << groups[a].probability << "\t" << ppm << endl;

          //Dead-end has a mass difference of a long or short arm.
          if (bDE && groups[a].ms3[xl[b].index].prob>groups[a].probability) {
            groups[a].type = xlDeadEnd;
            groups[a].sequence = groups[a].ms3[xl[b].index].peptide;
            groups[a].proteinS = processProteins(groups[a].ms3[xl[b].index].protein);
            groups[a].probability = groups[a].ms3[xl[b].index].prob;
            groups[a].calcNeutMassG= nm;
            groups[a].ppmG = (mm-nm)/nm*1e6;
            groups[a].offset = io;
            groups[a].posA=groups[a].ms3[xl[b].index].pos;
            groups[a].posB=0;
            groups[a].xlMass= params->crosslinker.xlMass + vc->at(d);
          }
        }

        //XL has two peptides plus the spacer arm
        for (size_t c = b + 1; c < xl.size(); c++) { //check all combinations with other MS3 peptides
          double xm = xl[b].mass + xl[c].mass + params->crosslinker.xlMass - io * ISOTOPE_OFFSET;
          double dif = mm - xm;

          //if(a==867) cout << c << "\t" << xm << "\t" << dif << "\t" << io << endl;
          
          //crosslinks that have the same probability as a different result take precedence.
          // TODO
          // alphabetize cross-linked peptides before storing them; this prevents confusion with duplications downstream
          if (fabs(dif) < ppm && groups[a].ms3[xl[b].index].prob* groups[a].ms3[xl[c].index].prob>=groups[a].probability*groups[a].probability) {
            groups[a].type = xlXL;
            groups[a].sequence = groups[a].ms3[xl[b].index].peptide + "+" + groups[a].ms3[xl[c].index].peptide;
            groups[a].proteinS = processProteins(groups[a].ms3[xl[b].index].protein) + "+" + processProteins(groups[a].ms3[xl[c].index].protein);
            groups[a].probability= groups[a].ms3[xl[b].index].prob * groups[a].ms3[xl[c].index].prob;
            groups[a].calcNeutMassG = xm;
            groups[a].ppmG = (mm - xm) / xm * 1e6;
            groups[a].offset=io;
            groups[a].posA = groups[a].ms3[xl[b].index].pos;
            groups[a].posB = groups[a].ms3[xl[c].index].pos;
            groups[a].xlMass = params->crosslinker.xlMass;
          }
        }
      }
    }

    //Check to see if there are partial results if nothing found by this point.
    if (groups[a].type==xlUnknown) {
      groups[a].xlMass=0;
      for (size_t b = 0; b < groups[a].ms3.size(); b++) {
        if (groups[a].ms3[b].prob < groups[a].probability) continue;

        //Cross-link is incomplete if there is evidence of a stub, but nothing adds up.
        if (groups[a].ms3[b].stubMass != 0) {
          groups[a].type = xlIncomplete;
          groups[a].probability=groups[a].ms3[b].prob;
          groups[a].calcNeutMassG=groups[a].ms3[b].pepMass;
          groups[a].ppmG=(mm- groups[a].ms3[b].pepMass)/ groups[a].ms3[b].pepMass*1e6;
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

void VingData::exportJSON(){
  string outFile = params->output + ".json";
  FILE* f = fopen(outFile.c_str(), "wt");

  fprintf(f, "{\n");
  fprintf(f, " \"ScanGroup\":[\n");
  for (size_t a = 0; a < groups.size(); a++) {
    fprintf(f, "  {\n");
    fprintf(f, "   \"GroupID\":%d,\n", (int)a);
    fprintf(f, "   \"Designation\":");
    switch (groups[a].type) {
    case xlDeadEnd: fprintf(f, "\"DeadEnd\",\n"); break;
    case xlSingle: fprintf(f, "\"Single\",\n"); break;
    case xlLoop: fprintf(f, "\"Loop\",\n"); break;
    case xlXL: fprintf(f, "\"XL\",\n"); break;
    case xlIncomplete: fprintf(f, "\"Incomplete\",\n"); break;
    default: fprintf(f, "\"Unknown\",\n"); break;
    }
    fprintf(f, "   \"Probability\":%.4lf,\n", groups[a].probability);
    fprintf(f, "   \"MSFile\":\"%s\",\n", files[groups[a].fileID].base.c_str());
    fprintf(f, "   \"Sequence\":\"%s\",\n", groups[a].sequence.c_str());
    fprintf(f, "   \"Proteins\":\"%s\",\n", groups[a].proteinS.c_str());

    fprintf(f, "   \"MS2_Scan\":{\n");
    fprintf(f, "    \"ScanNumber\":%d\n", groups[a].scan);
    fprintf(f, "    \"m/z\":%.4lf\n", groups[a].mz);
    fprintf(f, "    \"Charge\":%d\n", groups[a].charge);
    fprintf(f, "    \"Monoisotopic_Mass\":%.4lf\n", groups[a].monoMZ * groups[a].charge - groups[a].charge * 1.007276466);
    fprintf(f, "    \"Peptide\":\"%s\"\n", groups[a].peptide.c_str());
    fprintf(f, "    \"Proteins\":\"%s\"\n", processProteins(groups[a].protein).c_str());
    fprintf(f, "    \"Peptide_Neutral_Mass\":%.4lf\n", groups[a].calcNeutMass);
    fprintf(f, "    \"Probability\":%.4lf\n", groups[a].probabilityMS2);
    fprintf(f, "   }\n");

    fprintf(f, "   \"MS3_Scan\":[\n");
    for(size_t b=0;b < groups[a].ms3.size();b++) {
      fprintf(f,"    {\n");
      fprintf(f, "     \"ScanNumber\":%d\n", groups[a].ms3[b].scan);
      fprintf(f, "     \"m/z\":%.4lf\n", groups[a].ms3[b].mz);
      fprintf(f, "     \"Charge\":%d\n", groups[a].ms3[b].charge);
      fprintf(f, "     \"Monoisotopic_Mass\":%.4lf\n", groups[a].ms3[b].mz * groups[a].ms3[b].charge - groups[a].ms3[b].charge * 1.007276466);
      fprintf(f, "     \"Peptide\":\"%s\"\n", groups[a].ms3[b].peptide.c_str());
      fprintf(f, "     \"Proteins\":\"%s\"\n", processProteins(groups[a].ms3[b].protein).c_str());
      fprintf(f, "     \"Precursor_Neutral_Mass\":%.4lf\n", groups[a].ms3[b].pepMass);
      fprintf(f, "     \"Peptide_Neutral_Mass\":%.4lf\n", groups[a].ms3[b].pepMass - groups[a].ms3[b].stubMass);
      fprintf(f, "     \"Probability\":%.4lf\n", groups[a].ms3[b].prob);
      fprintf(f, "    }");
      if (b < groups[a].ms3.size() - 1) fprintf(f, ",");
      fprintf(f, "\n");
    }
    fprintf(f, "   ]\n");

    fprintf(f, "  }");
    if(a<groups.size()-1) fprintf(f,",");
    fprintf(f, "\n");
  }
  fprintf(f, " ]\n");
  fprintf(f, "}\n");
}

void VingData::exportProXL(){

  typedef struct sProXLPSM{
    string scanFile;
    int scanNumber;
    int precursorCharge;
    double probability;
  } sProXLPSM;

  typedef struct sProXL{
    eXLType type;
    string reportedPeptide;
    string peptideA;
    string peptideB;
    int posA;
    int posB;
    vector<sProXLPSM> psm;
  } sProXL;

  vector<sProXL> v;
  map<string,size_t> m;
  size_t index;

  for(size_t a=0;a<groups.size();a++){
    if(groups[a].type==xlUnknown) continue;
    if(groups[a].type==xlIncomplete) continue;
    
    sProXL p;
    map<string,size_t>::iterator i=m.find(groups[a].sequence);
    if(i==m.end()){
      p.reportedPeptide=groups[a].sequence;
      p.type=groups[a].type;
      bool bFirst=true;
      for(size_t b=0;b<p.reportedPeptide.size();b++){
        if(p.reportedPeptide[b]>='A' && p.reportedPeptide[b]<='Z'){
          if(bFirst) p.peptideA+=p.reportedPeptide[b];
          else p.peptideB += p.reportedPeptide[b];
        } else if(p.reportedPeptide[b]=='+') bFirst=false;
      }
      p.posA=groups[a].posA;
      p.posB=groups[a].posB;
      index=v.size();
      v.push_back(p);
      m.insert(pair<string,size_t>(p.reportedPeptide,index));
    } else index=i->second;

    sProXLPSM psm;
    psm.scanFile=files[groups[a].fileID].base+ files[groups[a].fileID].ext;
    psm.scanNumber=groups[a].scan;
    psm.probability=groups[a].probability;
    v[index].psm.push_back(psm);
    
  }

  string outFile=params->output+".proxl.xml";
  FILE* f = fopen(outFile.c_str(), "wt");

  fprintf(f,"<?xml version=\"1.0\" encoding=\"UTF - 8\" standalone=\"yes\"?>\n");
  fprintf(f,"<proxl_input fasta_filename = \"%s\">\n",dbName.c_str());

  exportProXLSearchProgramInfo(f);
  exportProXLLinkers(f);

  fprintf(f," <reported_peptides>\n");
  for(size_t a=0;a<v.size();a++){
    fprintf(f, "  <reported_peptide reported_peptide_string=\"%s\" ",v[a].reportedPeptide.c_str());
    if(v[a].type==xlXL) fprintf(f, " type=\"crosslink\">\n");
    else fprintf(f, " type=\"single\">\n");
    fprintf(f, "   <peptides>\n");
    fprintf(f, "    <peptide sequence=\"%s\">\n",v[a].peptideA.c_str());
    if(v[a].type==xlXL){
      fprintf(f, "     <linked_positions>\n");
      fprintf(f, "      <linked_position position=\"%d\">\n", v[a].posA);
      fprintf(f, "     </linked_positions>\n");
    }
    fprintf(f, "    </peptide>\n");
    if(!v[a].peptideB.empty()) {
      fprintf(f, "    <peptide sequence=\"%s\">\n", v[a].peptideB.c_str());
      fprintf(f, "     <linked_positions>\n");
      fprintf(f, "      <linked_position position=\"%d\">\n",v[a].posB);
      fprintf(f, "     </linked_positions>\n");
      fprintf(f, "    </peptide>\n");
    }
    fprintf(f, "   </peptides>\n");
    fprintf(f, "   <psms>\n");
    for(size_t b=0;b<v[a].psm.size();b++){
      fprintf(f,"    <psm scan_file_name=\"%s\" scan_number=\"%d\" precursor_charge=\"%d\">\n",v[a].psm[b].scanFile.c_str(),v[a].psm[b].scanNumber,v[a].psm[b].precursorCharge);
      fprintf(f,"     <filterable_psm_annotations>\n");
      fprintf(f,"      <filterable_psm_annotation search_program=\"iProphet\" annotation_name=\"IProphet Score\" value=\"%.4lf\"/>\n",v[a].psm[b].probability);
      fprintf(f,"     </filterable_psm_annotations>\n");
      fprintf(f,"    </psm>\n");
    }
    fprintf(f, "   </psms>\n");
    fprintf(f, "  </reported_peptide>\n");
  }
  fprintf(f, " </reported_peptide>\n");
  fprintf(f, "</proxl_input>\n");
  fclose(f);
}

void VingData::exportProXLLinkers(FILE* f){
  fprintf(f, " <linkers>\n");
  fprintf(f, "  <linker name = \"DSSO\">\n");
  fprintf(f, "   <monolink_masses>\n");
  fprintf(f, "    <monolink_mass mass=\"176.014330\" />\n");
  fprintf(f, "    <monolink_mass mass=\"279.077700\" />\n");
  fprintf(f, "   </monolink_masses>\n");
  fprintf(f, "   <crosslink_masses>\n");
  fprintf(f, "    <crosslink_mass mass=\"158.003765\" />\n");
  fprintf(f, "   </crosslink_masses>\n");
  fprintf(f, "   <linked_ends>\n");
  fprintf(f, "    <linked_end>\n");
  fprintf(f, "     <residues>\n");
  fprintf(f, "      <residue>K</residue>\n");
  fprintf(f, "     </residues>\n");
  fprintf(f, "     <protein_termini>\n");
  fprintf(f, "      <protein_terminus terminus_end=\"n\" distance_from_terminus=\"0\" />\n");
  fprintf(f, "      <protein_terminus terminus_end=\"n\" distance_from_terminus=\"1\" />\n");
  fprintf(f, "     </protein_termini>\n");
  fprintf(f, "    </linked_end>\n");
  fprintf(f, "   </linked_ends>\n");
  fprintf(f, "  </linker>\n");
  fprintf(f, " </linkers>\n");
}

void VingData::exportProXLSearchProgramInfo(FILE* f){
  fprintf(f, " <search_program_info>\n");
  fprintf(f, "  <search_programs>\n");
  fprintf(f, "   <search_program name=\"iProphet\" display_name=\"iProphet\" version=\"TPP v6.2.0 Parhelion - dev, Build 202208191617 - 8701 (Windows_NT - x86_64)\">\n");
  fprintf(f, "    <psm_annotation_types>\n");
  fprintf(f, "     <filterable_psm_annotation_types>\n");
  fprintf(f, "      <filterable_psm_annotation_type name=\"IProphet Score\" description=\"iProphet Probability Score\" filter_direction=\"above\" default_filter=\"true\" default_filter_value=\"0.9\" />\n");
  fprintf(f, "     </filterable_psm_annotation_types>\n");
  fprintf(f, "    </psm_annotation_types>\n");
  fprintf(f, "   </search_program>\n");
  //fprintf(f, "   <search_program name=\"Ving\" display_name=\"Ving\" version=\"pfft...whatever\">\n");
  //fprintf(f, "    <psm_annotation_types>\n");
  //fprintf(f, "     <filterable_psm_annotation_types>\n");
  //fprintf(f, "      <filterable_psm_annotation_type name=\"Score\" description=\"Ving Score\" filter_direction=\"above\" default_filter=\"false\" default_filter_value=\"0.9\" />\n");
  //fprintf(f, "     </filterable_psm_annotation_types>\n");
  //fprintf(f, "    </psm_annotation_types>\n");
  //fprintf(f, "   </search_program>\n");
  fprintf(f, "  </search_programs>\n");
  fprintf(f, "  <default_visible_annotations>\n");
  fprintf(f, "   <visible_psm_annotations>\n");
  fprintf(f, "    <search_annotation search_program=\"iProphet\" annotation_name=\"IProphet Score\" />\n");
  fprintf(f, "   </visible_psm_annotations>\n");
  fprintf(f, "  </default_visible_annotations>\n");
  fprintf(f, "  <annotation_sort_order>\n");
  fprintf(f, "   <psm_annotation_sort_order>\n");
  fprintf(f, "    <search_annotation search_program=\"iProphet\" annotation_name=\"IProphet Score\" />\n");
  fprintf(f, "   </psm_annotation_sort_order>\n");
  fprintf(f, "  </annotation_sort_order>\n");
  fprintf(f, " </search_program_info>\n");
}

void VingData::exportResults2() {
  FILE* f = fopen(params->output.c_str(), "wt");

  //Massive Header
  fprintf(f, "GroupID\tDesignation\tProbability\tMSFile\tMS2ScanNum\tMS2mz\tMS2charge\tMS2monoMass\tCalcNeutMass\tPPM\tIsotopeOffset\tXL_Mass\tSequence\tProtein(s)\tMS2Peptide\tMS2Protein\tMS2CalcNeutMass\tMS2Prob");
  for (size_t a = 0; a < maxMS3Count; a++) {
    fprintf(f, "\tMS3ScanNum-%d\tMS3mz\tMS3charge\tMS3selectedMass\tMS3Peptide\tMS3Protein\tMS3CalcNeutMass\tMS3PepMass\tMS3Prob", (int)a + 1);
  }
  fprintf(f, "\n");

  for (size_t a = 0; a < groups.size(); a++) {
    fprintf(f, "%d", (int)a);
    switch (groups[a].type) {
    case xlDeadEnd: fprintf(f, "\tDeadEnd"); break;
    case xlSingle: fprintf(f, "\tSingle"); break;
    case xlLoop: fprintf(f, "\tLoop"); break;
    case xlXL: fprintf(f, "\tXL"); break;
    case xlIncomplete: fprintf(f, "\tIncomplete"); break;
    default: fprintf(f, "\tUnknown"); break;
    }
    fprintf(f, "\t%.4lf", groups[a].probability);
    fprintf(f,"\t%s",files[groups[a].fileID].base.c_str());
    fprintf(f, "\t%d", groups[a].scan);
    fprintf(f, "\t%.4lf", groups[a].mz);
    fprintf(f, "\t%d", groups[a].charge);
    fprintf(f, "\t%.4lf", groups[a].monoMZ * groups[a].charge - groups[a].charge * 1.007276466);
    fprintf(f, "\t%.4lf", groups[a].calcNeutMassG);
    fprintf(f, "\t%.2lf", groups[a].ppmG);
    fprintf(f,"\t%d",groups[a].offset);
    fprintf(f,"\t%.4lf", groups[a].xlMass);
    fprintf(f, "\t%s", groups[a].sequence.c_str());
    fprintf(f, "\t%s", groups[a].proteinS.c_str());
    if (!groups[a].peptide.empty()) {
      fprintf(f, "\t%s", groups[a].peptide.c_str());
      fprintf(f, "\t%s", processProteins(groups[a].protein).c_str());
      fprintf(f, "\t%.4lf", groups[a].calcNeutMass);
      fprintf(f, "\t%.4lf", groups[a].probabilityMS2);
    } else fprintf(f, "\tn/a\tn/a\t0\t0");

    //And every MS3
    size_t b = 0;
    while (b < groups[a].ms3.size()) {
      fprintf(f, "\t%d", groups[a].ms3[b].scan);
      fprintf(f, "\t%.4lf", groups[a].ms3[b].mz);
      fprintf(f, "\t%d", groups[a].ms3[b].charge);
      fprintf(f, "\t%.4lf", groups[a].ms3[b].mz * groups[a].ms3[b].charge - groups[a].ms3[b].charge * 1.007276466);
      fprintf(f, "\t%s", groups[a].ms3[b].peptide.c_str());
      fprintf(f, "\t%s", processProteins(groups[a].ms3[b].protein).c_str());
      fprintf(f, "\t%.4lf", groups[a].ms3[b].pepMass);
      fprintf(f, "\t%.4lf", groups[a].ms3[b].pepMass-groups[a].ms3[b].stubMass);
      fprintf(f, "\t%.4lf", groups[a].ms3[b].prob);
      b++;
    }
    while (b < maxMS3Count) { //blanks
      fprintf(f, "\tn/a\t0\t0\t0\tn/a\tn/a\t0\t0\t0");
      b++;
    }
    fprintf(f, "\n");
  }

  fclose(f);

}

//Reads in PepXML from MS2 search results and assigns information to sMS2 objects
bool VingData::importMS2SearchResults() {
  NeoPepXMLParser p;
  p.read(params->ms2SearchResult.c_str());
  cout << "Done reading" << endl;

  bMS2Prophet[0]=false;
  bMS2Prophet[1]=false;
  for(size_t a=0;a<p.msms_pipeline_analysis[0].analysis_summary.size();a++){
    CnpxAnalysisSummary* as=&p.msms_pipeline_analysis[0].analysis_summary[a];
    if(as->analysis.compare("peptideprophet")==0) bMS2Prophet[0]=true;
    if(as->analysis.compare("interprophet")==0) bMS2Prophet[1]=true;
  }
  if(params->probabilityType && !bMS2Prophet[1]){
    cout << "Error: MS2 search results do not contain iProphet analysis." << endl;
    return false;
  }
  if (!params->probabilityType && !bMS2Prophet[0]) {
    cout << "Error: MS2 search results do not contain PeptideProphet analysis." << endl;
    return false;
  }

  for(size_t z=0;z<p.msms_pipeline_analysis[0].msms_run_summary.size();z++){
    CnpxMSMSRunSummary* rs = &p.msms_pipeline_analysis[0].msms_run_summary[z];
    sMzML mf;
    size_t fileID=files.size();
    size_t pos=rs->base_name.find_last_of('\\');
    if(pos==string::npos) pos = rs->base_name.find_last_of('/');
    if(pos==string::npos) pos=0;
    if(pos>0) pos++;
    mf.base=rs->base_name.substr(pos);
    mf.ext=rs->raw_data;
    mf.file=rs->base_name + rs->raw_data;
    mf.base_name=rs->base_name;
    files.push_back(mf);

    if(dbName.empty()){
      CnpxSearchSummary* ss=&rs->search_summary[0];
      CnpxSearchDatabase* sd=&ss->search_database[0];
      pos=sd->local_path.find_last_of('\\');
      if (pos == string::npos) pos = sd->local_path.find_last_of('/');
      if (pos == string::npos) pos = 0;
      if (pos > 0) pos++;
      dbName= sd->local_path.substr(pos);
    }

    size_t vp = parseMzML(mf,fileID);
  
    for (size_t a = 0; a < rs->spectrum_query.size(); a++) {
      CnpxSpectrumQuery* sq = &rs->spectrum_query[a];

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

      groups[vp].protein.push_back(sh->protein);
      for(size_t b=0;b<sh->alternative_protein.size();b++) groups[vp].protein.push_back(sh->alternative_protein[b].protein);
      groups[vp].charge = sq->assumed_charge;
      groups[vp].calcNeutMass = sh->calc_neutral_pep_mass;

    /*for (size_t b = 0; b < sh->search_score.size(); b++) {
      CnpxSearchScore* ss = &sh->search_score[b];
      if (ss->name.compare("xcorr") == 0) groups[vp].ms3[v3].xcorr = atof(ss->value.c_str());
      if (ss->name.compare("expect") == 0) groups[vp].ms3[v3].eval = atof(ss->value.c_str());
    }*/

      for (size_t b = 0; b < sh->analysis_result.size(); b++) {
        if (!params->probabilityType && sh->analysis_result[b].analysis.compare("peptideprophet") == 0) {
          groups[vp].probabilityMS2 = sh->analysis_result[b].peptide_prophet_result.probability;
        } else if (params->probabilityType && sh->analysis_result[b].analysis.compare("interprophet") == 0) {
          groups[vp].probabilityMS2 = sh->analysis_result[b].interprophet_result.probability;
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
  }
  return true;
}

bool VingData::importMS3SearchResults() {
  NeoPepXMLParser p;
  if(!p.read(params->ms3SearchResult.c_str())) return false;
  cout << "Done reading" << endl;

  bMS3Prophet[0] = false;
  bMS3Prophet[1] = false;
  for (size_t a = 0; a < p.msms_pipeline_analysis[0].analysis_summary.size(); a++) {
    CnpxAnalysisSummary* as = &p.msms_pipeline_analysis[0].analysis_summary[a];
    if (as->analysis.compare("peptideprophet") == 0) bMS3Prophet[0] = true;
    if (as->analysis.compare("interprophet") == 0) bMS3Prophet[1] = true;
  }
  if (params->probabilityType && !bMS3Prophet[1]) {
    cout << "Error: MS3 search results do not contain iProphet analysis." << endl;
    return false;
  }
  if (!params->probabilityType && !bMS3Prophet[0]) {
    cout << "Error: MS3 search results do not contain PeptideProphet analysis." << endl;
    return false;
  }

  size_t vp = 0;
  size_t v3 = 0;

  for (size_t z = 0; z < p.msms_pipeline_analysis[0].msms_run_summary.size(); z++) {
    CnpxMSMSRunSummary* rs = &p.msms_pipeline_analysis[0].msms_run_summary[z];

    size_t pos = rs->base_name.find_last_of('\\');
    if (pos == string::npos) pos = rs->base_name.find_last_of('/');
    if (pos == string::npos) pos = 0;
    if (pos > 0) pos++;
    string base = rs->base_name.substr(pos);

    size_t x;
    for(x=0;x<files.size();x++){
      if(files[x].base.compare(base)==0) break;
    }
    if(x==files.size()){
      cout << "WARNING: " << base << " not in list of files analyzed. Skipping." << endl;
      continue;
    }
    size_t fileID=x;
    vp=files[fileID].startPos;
    while (vp < groups.size() && groups[vp].ms3.size() == 0) vp++;

    for (size_t a = 0; a < rs->spectrum_query.size(); a++) {
      CnpxSpectrumQuery* sq = &rs->spectrum_query[a];

      char* tok;
      char str[256];
      strcpy(str, sq->spectrum.c_str());
      tok = strtok(str, ".\n\r");
      tok = strtok(NULL, ".\n\r");
      int scan = atoi(tok);

      v3=0;
      while (vp < groups.size() && groups[vp].ms3[v3].scan < scan) {
        v3++;
        if (v3 == groups[vp].ms3.size()) {
          v3 = 0;
          vp++;
          while (vp < groups.size() && groups[vp].ms3.size() == 0) vp++;
        }
      }

      if (groups[vp].ms3[v3].scan != scan) {
        cout << "No match to MS3 " << scan << " in " << files[fileID].base_name << "\t" << files[fileID].base << endl;
        cout << "v3=" << v3 << "\t" << "scan=" << groups[vp].ms3[v3].scan << "\t" << "size=" << groups[vp].ms3.size() << endl;
        cout << vp << " of " << groups.size() << " has fID: " << groups[vp].fileID << "\t" << fileID << "\t" << groups[vp].scan << endl;
        for(size_t b=0;b<groups[vp].ms3.size();b++){
          cout << " " << b << " " << groups[vp].ms3[b].scan << endl;
        }

        continue;
      }

      if (sq->search_result[0].search_hit.size() == 0) continue;

      CnpxSearchHit* sh = &sq->search_result[0].search_hit[0];
      if (sh->modification_info.empty()) groups[vp].ms3[v3].peptide = sh->peptide;
      else groups[vp].ms3[v3].peptide = sh->modification_info[0].modified_peptide;

      groups[vp].ms3[v3].protein.push_back(sh->protein);
      for(size_t b=0;b<sh->alternative_protein.size();b++) groups[vp].ms3[v3].protein.push_back(sh->alternative_protein[b].protein);

      groups[vp].ms3[v3].charge = sq->assumed_charge;
      groups[vp].ms3[v3].pepMass = sh->calc_neutral_pep_mass;

      for (size_t b = 0; b < sh->search_score.size(); b++) {
        CnpxSearchScore* ss = &sh->search_score[b];
        if (ss->name.compare("xcorr") == 0) groups[vp].ms3[v3].xcorr = atof(ss->value.c_str());
        if (ss->name.compare("expect") == 0) groups[vp].ms3[v3].eval = atof(ss->value.c_str());
      }

      for (size_t b = 0; b < sh->analysis_result.size(); b++) {
        if (!params->probabilityType && sh->analysis_result[b].analysis.compare("peptideprophet") == 0) {
          groups[vp].ms3[v3].prob = sh->analysis_result[b].peptide_prophet_result.probability;
        } else if(params->probabilityType && sh->analysis_result[b].analysis.compare("interprophet") == 0) {
          groups[vp].ms3[v3].prob = sh->analysis_result[b].interprophet_result.probability;
        }
      }

      //check for stubs
      //TODO: track site of stubs for inferring link positions.
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
                groups[vp].ms3[v3].pos=mam->position;
                break;
              } else if (fabs(params->crosslinker.stubA[c].mass + params->crosslinker.stubA[c].reaction - mam->variable) < 0.1) {
                groups[vp].ms3[v3].stubMass = params->crosslinker.stubA[c].mass + params->crosslinker.stubA[c].reaction;
                groups[vp].ms3[v3].stubA = true;
                groups[vp].ms3[v3].pos = mam->position;
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
                groups[vp].ms3[v3].pos = mam->position;
                break;
              } else if (fabs(params->crosslinker.stubB[c].mass + params->crosslinker.stubB[c].reaction - mam->variable) < 0.1) {
                groups[vp].ms3[v3].stubMass = params->crosslinker.stubB[c].mass + params->crosslinker.stubB[c].reaction;
                groups[vp].ms3[v3].stubA = false;
                groups[vp].ms3[v3].pos = mam->position;
                break;
              }
            }
            if (c < params->crosslinker.stubB.size()) break;
          }
        }
      }
    }
  }
  return true;
}

size_t VingData::parseMzML(sMzML& fn, size_t id){
  MSReader r;
  Spectrum s;

  r.setFilter(MS2);
  r.addFilter(MS3);
  r.addFilter(MS1);

  vector<sMS2> cycle;
  bool bMS2 = false;

  size_t start_pos=groups.size();

  cout << "Reading: " << fn.file << " ... ";
  if(!r.readFile(fn.file.c_str(), s)){
    string st=fn.base+fn.ext;
    cout << " Not found, attempting CWD: " << st << " ... ";
    if (!r.readFile(st.c_str(), s)) return SIZE_MAX;
  }
    
  //Set progress meter
  int lastScan=r.getLastScan();
  int iPercent=0;
  printf("%2d%%", iPercent);
  fflush(stdout);

  int count=0;
  while (s.getScanNumber() > 0) {

    if (s.getMsLevel() == 1) {
      for (size_t a = 0; a < cycle.size(); a++) {
        if (bMS2) {
          cycle[a].fileID=id;
          groups.push_back(cycle[a]);
          if(cycle[a].ms3.size()>maxMS3Count)maxMS3Count=cycle[a].ms3.size();
          count++;
        }
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

    //Update progress meter
    int iTmp = (int)((double)s.getScanNumber() / lastScan * 100);
    if (iTmp > iPercent) {
      iPercent = iTmp;
      printf("\b\b\b%2d%%", iPercent);
      fflush(stdout);
    }

    r.readFile(NULL, s);
  }

  //Finalize progress meter
  //printf("\b\b\b100%%");
  cout << endl;

  for (size_t a = 0; a < cycle.size(); a++) {
    if (bMS2) {
      cycle[a].fileID=id;
      groups.push_back(cycle[a]);
      if (cycle[a].ms3.size() > maxMS3Count)maxMS3Count = cycle[a].ms3.size();
      count++;
    }
  }

  cout << count << " new groups." << endl;
  cout << groups.size() << " total groups." << endl;
  files[id].startPos=start_pos;

  //Diagnostics only
  //for(size_t a=start_pos;a<groups.size(); a++ ){
  //  cout << "New " << id << "\t" << groups[a].scan;
  //  for(size_t b=0;b<groups[a].ms3.size();b++) cout << "\t" << groups[a].ms3[b].scan;
  //  cout << endl;
  //}
  return start_pos;
}

string VingData::processProteins(vector<string>& proteins){
  bool bDecoy = true;
  for (size_t e = 0; e < proteins.size(); e++) {
    if (proteins[e].find(params->decoy) != 0) {
      bDecoy = false;
      break;
    }
  }

  string p;
  for (size_t e = 0; e < proteins.size(); e++) {
    if (bDecoy) {
      if (!p.empty()) p += ";";
      p += proteins[e];
    } else {
      if (proteins[e].find(params->decoy) == 0) continue;
      if (!p.empty()) p += ";";
      p += proteins[e];
    }
  }

  return p;
}

bool VingData::compareXL(sXLPep& a, sXLPep& b) {
  return a.mass > b.mass;
}

bool VingData::compareMZ(sPepMass& a, sPepMass& b) {
  return a.mz < b.mz;
}
