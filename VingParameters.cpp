#include "VingParameters.h"

using namespace std;

bool VingParameters::checkToken(char* tok) {
  if (tok == NULL) {
    //if (log != NULL) log->addError("Error in [XL_PARAMS]");
    //else printf("Error in [XL_PARAMS]\n");
    return false;
  }
  return true;
}

void VingParameters::exportDefault(string ver){
  FILE* f = fopen("ving_default_params.conf", "wt");

  fprintf(f, "# Ving %s parameter file\n\n", ver.c_str());
  fprintf(f, "# This file was autogenerated. Please alter and rename before use.\n");
  fprintf(f, "# All parameters are separated from their values by an equals sign ('=')\n");
  fprintf(f, "# Anything after a '#' will be ignored for the remainder of the line.\n");
  fprintf(f, "# Remove the '#' at the beginning of a parameter line to activate that parameter.\n");
  //fprintf(f, "# Please see online documentation at:\n");
  //fprintf(f, "# http://www.kojak-ms.org/param");

  fprintf(f, "\n#\n");
  fprintf(f, "# XL_PARAMS defines pre - configured crosslinker settings.The table can be expanded to suit new crosslinkers.\n");
  fprintf(f, "# This table MUST precede the crosslink parameter.\n");
  fprintf(f, "# Columns are : [ID Num] [Name] [Quenching Reagent] [Target A] [Target B] [XL Mass] [MonoCap A] [MonoCap B] [Stub A, reaction] [Stub B, reaction] [Reporter Ion(s)]\n");
  fprintf(f, "#\n");
  fprintf(f, "[XL_PARAMS]\n");
  fprintf(f, "1   Hyp1    na     H    H    362.026364  -13.979265            -13.979265            13.979265,0;348.047099,0                               13.979265,0;348.047099,0                               334.067834\n");
  fprintf(f, "2   Hyp8    na     H    H    566.162524  -13.979265            -13.979265            13.979265,0;552.183259,0                               13.979265,0;552.183259,0                               538.203994\n");
  fprintf(f, "3   DSSO    NH2    nK   nK   158.003765  17.026549             17.026549             54.010565,0;103.993200,-18.010565                      54.010565,0;103.993200,-18.010565                      0\n");
  fprintf(f, "4   DSSO    H2O    nK   nK   158.003765  18.010565             18.010565             54.010565,0;103.993200,-18.010565                      54.010565,0;103.993200,-18.010565                      0\n");
  fprintf(f, "5   DHSO    na     H    H    250.0372    -13.979265            -13.979265            13.979265,0;100.0273,0;150.0099,-18.010565;236.0579,0  13.979265,0;100.0273,0;150.0099,-18.010565;236.0579,0  0\n");
  fprintf(f, "6   DSSO    Tris   nK   nK   158.003765  121.073965;18.010565  121.073965;18.010565  54.010565,0;103.993200,-18.010565                      54.010565,0;103.993200,-18.010565                      0\n");
  fprintf(f, "[END_XL_PARAMS]\n");
  fprintf(f, "\n");
  fprintf(f, "crosslinker = 6\n");
  fprintf(f, "ms2_search_result = ms2.ipro.pep.xml\n");
  fprintf(f, "ms3_search_result = ms3.ipro.pep.xml\n");
  fprintf(f, "output_name = ving_results.txt\n");
  fprintf(f, "\n");
  fprintf(f, "probability_type = 1  #0 = PeptideProphet, 1 = iProphet\n");
  fprintf(f, "ppm_tolerance = 20.0  #parent ion mass tolerance, in ppm, default=10.0\n");
  fprintf(f, "decoy_identifier = DECOY\n");

  fclose(f);
}

bool VingParameters::parse(const char* fn){
  FILE* f;
  char str[4096];
  xlTable.clear();

  f = fopen(fn, "rt");
  if (f == NULL) {
    //if (log != NULL) log->addError("Cannot open config file!");
    //else printf("Cannot open config file!\n");
    return false;
  }

  vector<string> vParamLines;
  bool bReadXLTable = false;
  while (!feof(f)) {
    if (!fgets(str, 4096, f)) continue;
    if (strlen(str) == 0) continue;
    if (str[0] == '#') continue;
    string s = str;

    if (s.find("[XL_PARAMS]") == 0) {
      bReadXLTable = true;
      continue;
    }
    if (s.find("[END_XL_PARAMS]") == 0) {
      bReadXLTable = false;
      continue;
    }

    if (bReadXLTable) {
      if(!parseXLTable(str)){
        return false;
      }
      continue;
    }

    parseParam(str);

  }
  fclose(f);

  return true;
}

void VingParameters::parseParam(const char* str){
  char* tok;
  char c_cmd[4096];
  char param[32];
  char tmpstr[256];

  string tstr;
  vector<string> values;
  strcpy(c_cmd, str);

  //Pre-process entire line to remove characters that should not be read
  //Replace first # with a terminator
  tok = strstr(c_cmd, "#");
  if (tok != NULL) strncpy(tok, "\0", 1);

  //if we have only white space, exit here
  strcpy(tmpstr, c_cmd);
  tok = strtok(tmpstr, " \t\n\r");
  if (tok == NULL) return;

  //Check if we have a parameter (has '=' in it) or lots of random text.
  tok = strstr(c_cmd, "=");
  if (tok == NULL) {
    //if (log != NULL) log->addParameterWarning("Unknown parameter line in config file: " + string(cmd));
    //else printf("Unknown parameter line in config file: %s\n", cmd);
    return;
  }

  //Process parameters
  //Read parameter into param name (before = sign) and value (after = sign)
  tok = strtok(c_cmd, " \t=\n\r");
  if (tok == NULL) return;
  strcpy(param, tok);
  tok = strtok(NULL, " \t=\n\r");
  if (tok == NULL) {
    //warn(param, 0);
    return;
  } else {
    while (tok != NULL) {
      tstr = tok;
      values.push_back(tstr);
      tok = strtok(NULL, " \t=\n\r");
    }
  }

  //Look up parameter and assign the value
  if (strcmp(param, "crosslinker") == 0) {
    int id=atoi(values[0].c_str());
    size_t a;
    for(a=0;a<xlTable.size();a++){
      if(xlTable[a].id==id) break;
    }
    if(a==xlTable.size()){
      //warn
    } else {
      crosslinker=xlTable[a];
    }
  } else if (strcmp(param, "decoy_identifier") == 0) {
    decoy = values[0];
  } else if (strcmp(param, "ms2_search_result") == 0) {
    ms2SearchResult=values[0];
  } else if (strcmp(param, "ms3_search_result") == 0) {
    ms3SearchResult = values[0];
  } else if (strcmp(param, "mzml") == 0) {
    mzML = values[0];
  } else if (strcmp(param, "output_name") == 0) {
    output = values[0];
  } else if (strcmp(param, "ppm_tolerance") == 0) {
    ppmTolerance = atof(values[0].c_str());
  } else if (strcmp(param, "probability_type") == 0) {
    if(atoi(values[0].c_str())==1) probabilityType=true;
    else probabilityType=false;
  }
}

bool VingParameters::parseXLTable(char* str){
  vingXL x;
  char* tok = strtok(str, " \t\n\r");
  vector<double> stubs;
  vector<string> words;
  x.id = atoi(tok);

  size_t a;
  for (a = 0; a < xlTable.size(); a++) if (xlTable[a].id == x.id) break;
  if (a < xlTable.size()) {
    //if (log != NULL) log->addError("Duplicate ID number in [XL_PARAMS]");
    //else printf("Duplicate ID number in [XL_PARAMS]\n");
    return false;
  }

  tok = strtok(NULL, " \t\n\r");
  if (checkToken(tok)) x.name = tok;
  else return false;

  tok = strtok(NULL, " \t\n\r");
  //if (checkToken(tok)) x.xlQuench = tok;
  //else return false;

  tok = strtok(NULL, " \t\n\r");
  if (checkToken(tok)) x.targetA = tok;
  else return false;

  tok = strtok(NULL, " \t\n\r");
  if (checkToken(tok)) x.targetB = tok;
  else return false;

  tok = strtok(NULL, " \t\n\r");
  if (checkToken(tok)) x.xlMass = atof(tok);
  else return false;

  tok = strtok(NULL, " \t\n\r");
  if (checkToken(tok)) splitMasses(tok, x.capMassA,';');
  else return false;

  tok = strtok(NULL, " \t\n\r");
  if (checkToken(tok)) splitMasses(tok, x.capMassB,';');
  else return false;

  tok = strtok(NULL, " \t\n\r");
  if(checkToken(tok)) splitWords(tok,words,';');
  else return false;
  if (!words.empty()) {
    for(size_t a=0;a<words.size();a++) {
      splitMasses(words[a], stubs,',');
      vingStub s;
      s.mass=stubs[0];
      if(stubs.size()>1) s.reaction=stubs[1];
      x.stubA.push_back(s);
    }
  } else return false;
  
  tok = strtok(NULL, " \t\n\r");
  if (checkToken(tok)) splitWords(tok, words, ';');
  else return false;
  if (!words.empty()) {
    for (size_t a = 0; a < words.size(); a++) {
      splitMasses(words[a], stubs, ',');
      vingStub s;
      s.mass = stubs[0];
      if (stubs.size() > 1) s.reaction = stubs[1];
      x.stubB.push_back(s);
    }
  } else return false;

  tok = strtok(NULL, " \t\n\r");
  if (checkToken(tok)) splitMasses(tok, x.reporterIons);
  else return false;

  xlTable.push_back(x);
  return true;
}

void VingParameters::splitMasses(char*& c, vector<double>& v) {
  string s;
  v.clear();
  for (size_t a = 0; a < strlen(c); a++) {
    if (c[a] == ',') {
      double m=atof(s.c_str());
      if(m!=0) v.push_back(m);
      s.clear();
    } else s += c[a];
  }
  v.push_back(atof(s.c_str()));
}

void VingParameters::splitMasses(string str, vector<double>& v, char token) {
  string s;
  v.clear();
  for (size_t a = 0; a < str.size(); a++) {
    if (str[a] == token) {
      double m = atof(s.c_str());
      if (m != 0) v.push_back(m);
      s.clear();
    } else s += str[a];
  }
  if(!s.empty()) v.push_back(atof(s.c_str()));
}

void VingParameters::splitWords(char*& c, std::vector<std::string>& v, char token){
  string s;
  v.clear();
  for (size_t a = 0; a < strlen(c); a++) {
    if (c[a] == token) {
      v.push_back(s);
      s.clear();
    } else s += c[a];
  }
  if(!s.empty()) v.push_back(s);
}
