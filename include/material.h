#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>
#include <string>
#include <map>

#include "inputparser.h"

// Notes:
// Should DN data be private?
// Are all of these getters/setters necesseary?

/// Material object
/**
 *  Cross sections
 */
class Material
{
 public:
  Material(unsigned int G, unsigned int N);
  ~Material();

  /// Get the total cross section
  double* getSigma_t(unsigned int group)
  { return &_sigma_t[group-1]; }

  /// Get an element of the scattering matrix
  double* getSigma_s(unsigned int group_in, unsigned int group_out, unsigned int order = 0)
  { return &_sigma_s_ker[ si(group_in, group_out,order) ]; }

  /// Get the total scattering cross section
  double* getSigma_s(unsigned int group_in)
  { return &_sigma_s[ group_in-1 ]; }

  /// Get the production cross section
  double* getNuSigma_f(unsigned int group)
  { return &_nu_sigma_f[group-1]; }

  /// Get the fission spectrum
  double* getFissionSpectrum(unsigned int group)
  { return &_fissionSpectrum[group-1]; }


  /// Get the neutron speed
  double* getNeutronSpeed(unsigned int group)
  { return &_speed[group-1]; };

  void putSigma_t(unsigned int group, double crossSection);
  void putSigma_t(std::vector<double>& crossSection);
  
  void putNuSigma_f(unsigned int group, double crossSection);
  void putNuSigma_f(std::vector<double>& crossSection);

  void putSigma_s(unsigned int groupIn, unsigned int groupOut, unsigned int order, double crossSection);
  void putSigma_s(std::vector<double>& crossSection);

  void putFissionSpectrum(unsigned int group, double fraction);
  void putFissionSpectrum(std::vector<double>& spectrum);
  
  void putNeutronSpeed(unsigned int group, double speed);
  void putNeutronSpeed(std::vector<double>& speed);

  int scatteringOrder;
  int numGroups;

  bool hasDelayedNeutrons;
  int numDelayedGroups;
  double lambda[6];
  double beta[6];
  double betaEff;

  double* _sigma_t;
  double* _sigma_s_ker;
  double* _sigma_s;
  double* _nu_sigma_f;

  double* _sigma_t_adder;
  double* _nu_sigma_f_scaling;

  double* _fissionSpectrum;

  double* _speed;

 private:
  unsigned int si(unsigned int groupIn, unsigned int groupOut, unsigned int order)
  { return order*numGroups*numGroups + (groupOut-1)*numGroups + groupIn - 1; }

};

/// Factory for creating materials and keeping a list of available materials
class MaterialFactory
{
 public:
  ~MaterialFactory();
  int createMaterialsFromInput(InputParser& input);
  void readNextMaterial(InputParser& input, int matNum);
  Material* createNewMaterial(std::string& materialName, unsigned int block, unsigned int G, unsigned int N);
  Material* getMaterial(std::string& materialName);
  Material* getMaterial(unsigned int& materialID);
  void list();

  int numGroups;
  
  static std::map<std::string, Material*> _materialMap;
  static std::map<unsigned int, Material*> _materialID;
};

#endif
