#ifndef _VINGDATA_H
#define _VINGDATA_H


#include "VingParameters.h"
#include "MSReader.h"
#include "NeoPepXMLParser.h"

#include <algorithm>

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
  sMS3() {
    scan = 0;
    xcorr = 0;
    eval = 0;
    mz = 0;
    charge = 0;
    pepMass = 0;
    stubMass = 0;
    stubA=true;
  }
} sMS3;

typedef struct sMod {
  char aa;
  int pos;
  double mass;
}sMod;


typedef struct sMS2 {
  int scan;
  double mz;
  int charge;
  double monoMZ;
  double calcNeutMass;
  double calcNeutMassG;
  double ppmG;
  double probability;
  double probabilityMS2;
  std::string peptide;
  std::string sequence;
  std::string protein;
  std::string proteinS;
  eXLType type;
  std::vector<sMS3> ms3;
  std::vector<sMod> mods;
  sMS2() { clear(); }
  void clear() {
    scan = 0;
    mz = 0;
    charge = 0;
    monoMZ = 0;
    calcNeutMass = 0;
    calcNeutMassG=0;
    ppmG=0;
    probability = 0;
    probabilityMS2 = 0;
    type = xlUnknown;
    ms3.clear();
    mods.clear();
  }
} sMS2;

typedef struct sXLPep {
  size_t index;
  double mass;
  double stub;
} sXLPep;

class VingData {
public:
  VingData(VingParameters* p);
  ~VingData();

  std::vector<sMS2> groups;

  void assessXLType();
  bool importMS2SearchResults();
  bool importMS3SearchResults();
  bool parseMzML();

private:
  VingParameters* params = NULL;

  bool bMS2Prophet[2]={false};  //0=PeptideProphet, 1=iProphet
  bool bMS3Prophet[2]={false};

  static bool compareXL(sXLPep& a, sXLPep& b);
};

#endif

