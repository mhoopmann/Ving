#ifndef _VINGPARAMETERS_H
#define _VINGPARAMETERS_H

#include <string>
#include <vector>
#include <cstring>

typedef struct vingStub{
  double mass=0;
  double reaction=0;
} vingStub;

typedef struct vingXL {
  std::string name;
  std::string targetA;
  std::string targetB;
  int id=0;
  double xlMass = 0;
  std::vector<double> capMassA;
  std::vector<double> capMassB;
  std::vector<vingStub> stubA;
  std::vector<vingStub> stubB;
  std::vector<double> reporterIons;
} vingXL;

class VingParameters {
public:

  void exportDefault(std::string ver);
  bool parse(const char* fn);

  vingXL crosslinker;
  std::string mzML;
  std::string ms2SearchResult;
  std::string ms3SearchResult;
  std::string output;
  double ppmTolerance=10.0;
  bool probabilityType=false;

private:
  //double* log=NULL;
  std::vector<vingXL> xlTable;

  bool checkToken(char* tok);
  void parseParam(const char* str);
  bool parseXLTable(char* str);
  void splitMasses(char*& c, std::vector<double>& v);
  void splitMasses(std::string c, std::vector<double>& v, char token);
  void splitWords(char*& c, std::vector<std::string>& v, char token);

};

#endif

