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
  int scan;
  std::string peptide;
  std::string protein;
  double xcorr;
  double eval;
  double mz;
  double prob;
  double pepMass;   //mass including stub
  double stubMass;  //just the stub
  bool stubA;
  int charge;
  int pos;          //site of crosslinker, relative to peptide (1-based)
  sMS3() {
    scan = 0;
    xcorr = 0;
    eval = 0;
    mz = 0;
    charge = 0;
    pepMass = 0;
    stubMass = 0;
    stubA=true;
    pos=0;
    prob=0;
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
  eXLType type=xlUnknown;
  std::string peptide;
  std::string sequence;
  std::string protein;
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
    type = xlUnknown;
    ms3.clear();
    mods.clear();
    peptide.clear();
    sequence.clear();
    protein.clear();
    proteinS.clear();
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

class VingData {
public:
  VingData(VingParameters* p);
  ~VingData();

  std::vector<sMS2> groups;
  std::vector<sMzML> files;

  void assessXLType();
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

  void exportProXLLinkers(FILE* f);
  void exportProXLSearchProgramInfo(FILE* f);
};

#endif

