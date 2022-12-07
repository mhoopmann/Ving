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
  } else if (strcmp(param, "ms2_search_result") == 0) {
    ms2SearchResult=values[0];
  } else if (strcmp(param, "ms3_search_result") == 0) {
    ms3SearchResult = values[0];
  } else if (strcmp(param, "mzml") == 0) {
    mzML = values[0];
  } else if (strcmp(param, "output_name") == 0) {
    output = values[0];
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
  if (checkToken(tok)) x.capMassA=atof(tok);
  else return false;

  tok = strtok(NULL, " \t\n\r");
  if (checkToken(tok)) x.capMassB = atof(tok);
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
