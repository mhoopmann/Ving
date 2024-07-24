#ifndef _VINGDATA_H
#define _VINGDATA_H


#include "VingParameters.h"
#include "MSReader.h"
#include "NeoPepXMLParser.h"

#include <algorithm>
#include <cmath>

#define ISOTOPE_OFFSET 1.003354835
#define PROTON 1.007276466

enum eXLType {
  xlUnknown = 0,
  xlSingle,
  xlDeadEnd,
  xlLoop,
  xlXL,
  xlIncomplete,
  xlSize
};

typedef struct sMS3 {
  int scan=0;
  std::string peptide;
  std::vector<std::string> protein;
  double xcorr=0;
  double eval=0;
  double mz=0;
  double prob=0;
  double pepMass=0;   //mass including stub
  double stubMass=0;  //just the stub
  bool stubA=true;
  int charge=0;
  int pos=0;          //site of crosslinker, relative to peptide (1-based)
  std::string processProteins(std::string decoy) {
    bool bDecoy = true;
    for (size_t e = 0; e < protein.size(); e++) {
      if (protein[e].find(decoy) != 0) {
        bDecoy = false;
        break;
      }
    }

    std::string p;
    for (size_t e = 0; e < protein.size(); e++) {
      if (bDecoy) {
        if (!p.empty()) p += ";";
        p += protein[e];
      } else {
        if (protein[e].find(decoy) == 0) continue;
        if (!p.empty()) p += ";";
        p += protein[e];
      }
    }

    return p;
  }
} sMS3;

typedef struct sMod {
  char aa;
  int pos;
  double mass;
}sMod;


typedef struct sMS2 {
  size_t fileID=0;
  int scan=0;
  double mz=0;
  int charge=0;
  int offset=0;
  double monoMZ=0;
  double calcNeutMass=0;
  double calcNeutMassG=0;
  double ppmG=0;
  int posA=0;
  int posB=0;
  double probability=0;
  double probabilityMS2=0;
  double xlMass=0; //additional mass from the crosslinker. For deadends, includes the partially reacted reagent
  eXLType type=xlUnknown;
  std::string peptide;
  std::string sequence;
  std::vector<std::string> protein;
  std::string proteinS;
  std::vector<sMS3> ms3;
  std::vector<sMod> mods;
  void clear() {
    fileID = 0;
    scan = 0;
    mz = 0;
    charge = 0;
    offset = 0;
    monoMZ = 0;
    calcNeutMass = 0;
    calcNeutMassG = 0;
    ppmG = 0;
    posA = 0;
    posB = 0;
    probability = 0;
    probabilityMS2 = 0;
    xlMass=0;
    type = xlUnknown;
    ms3.clear();
    mods.clear();
    peptide.clear();
    sequence.clear();
    protein.clear();
    proteinS.clear();
  }
  std::string processProteins(std::string decoy) {
    bool bDecoy = true;
    for (size_t e = 0; e < protein.size(); e++) {
      if (protein[e].find(decoy) != 0) {
        bDecoy = false;
        break;
      }
    }

    std::string p;
    for (size_t e = 0; e < protein.size(); e++) {
      if (bDecoy) {
        if (!p.empty()) p += ";";
        p += protein[e];
      } else {
        if (protein[e].find(decoy) == 0) continue;
        if (!p.empty()) p += ";";
        p += protein[e];
      }
    }

    return p;
  }
} sMS2;

typedef struct sXLPep {
  size_t index;
  double mass;
  double stub;
} sXLPep;

typedef struct sMzML {
  std::string file;
  std::string base;
  std::string base_name;
  std::string ext;
  size_t startPos;
} sMzML;

typedef struct sPepMass{
  std::string peptide;
  double monoMass=0;
  double mz=0;
  int z=0;
  double matchMZ=0;
} sPepMass;

class VingData {
public:
  VingData(VingParameters* p);
  ~VingData();

  std::vector<sMS2> groups;
  std::vector<sMzML> files;

  void assessIncomplete(std::string massList, std::string protein="");
  void assessXLType();
  void exportJSON();
  void exportProXL();
  void exportResults2();
  bool importMS2SearchResults();
  bool importMS3SearchResults();
  size_t parseMzML(sMzML& fn, size_t id);

private:
  VingParameters* params = NULL;

  bool bMS2Prophet[2]={false};  //0=PeptideProphet, 1=iProphet
  bool bMS3Prophet[2]={false};

  size_t maxMS3Count;
  std::string dbName;

  static bool compareXL(sXLPep& a, sXLPep& b);
  static bool compareMZ(sPepMass& a, sPepMass& b);

  void exportProXLLinkers(FILE* f);
  void exportProXLSearchProgramInfo(FILE* f);
  std::string processProteins(std::vector<std::string>& proteins);
};

#endif

