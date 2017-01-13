#include "material.h"
#include "log.h"

#include <iostream>

// Instantiaion of static material list
std::map<std::string, Material*> MaterialFactory::_materialMap;
std::map<unsigned int, Material*> MaterialFactory::_materialID;

/// Multigroup constructor
Material::Material(unsigned int G, unsigned int N) :
  scatteringOrder(N), numGroups(G)
{
  _sigma_t = new double [numGroups];
  _sigma_s_ker = new double [numGroups*numGroups*(scatteringOrder+1)];
  _sigma_s = new double [numGroups];
  _nu_sigma_f = new double [numGroups];
  _fissionSpectrum = new double [numGroups];
  _speed = new double [numGroups];

  _sigma_t_adder = new double [numGroups];
  _nu_sigma_f_scaling = new double [numGroups];

  for (unsigned int g=0; g<numGroups; g++) {
    _sigma_t[g] = 0.0;
    _sigma_s[g] = 0.0;
    _nu_sigma_f[g] = 0.0;
    _fissionSpectrum[g] = 0.0;
    _speed[g] = 1.0;
    _sigma_t_adder[g] = 0.0;
    _nu_sigma_f_scaling[g] = 1.0;
  }

  hasDelayedNeutrons = false;
  for (int d=0; d<6; d++) {
    lambda[d] = 0.0;
    beta[d] = 0.0;
  }
  betaEff = 0.0;
}


Material::~Material()
{
  if (_sigma_t) delete [] _sigma_t;
  if (_sigma_s_ker) delete [] _sigma_s_ker;
  if (_sigma_s) delete [] _sigma_s;
  if (_nu_sigma_f) delete [] _nu_sigma_f;
  if (_sigma_t_adder) delete [] _sigma_t_adder;
  if (_nu_sigma_f_scaling) delete [] _nu_sigma_f_scaling;
  if (_fissionSpectrum) delete [] _fissionSpectrum;
  if (_speed) delete [] _speed;
}

void
Material::putSigma_t(unsigned int group, double crossSection)
{
  _sigma_t[group-1] = crossSection;
}

void
Material::putSigma_t(std::vector<double>& crossSection)
{
  for (int g=0; g<numGroups; g++) {
    _sigma_t[g] = crossSection[g];
  }
}

void
Material::putNuSigma_f(unsigned int group, double crossSection)
{
  _nu_sigma_f[group-1] = crossSection;
}

void
Material::putNuSigma_f(std::vector<double>& crossSection)
{
  for (int g=0; g<numGroups; g++) {
    _nu_sigma_f[g] = crossSection[g];
  }
}

void
Material::putSigma_s(unsigned int groupIn, unsigned int groupOut, unsigned int order, double crossSection)
{
  _sigma_s_ker[ si(groupIn, groupOut, order) ] = crossSection;
  if (order==0) _sigma_s[groupIn-1] += crossSection;
}

void
Material::putSigma_s(std::vector<double>& crossSection)
{
  if (crossSection.size()/(numGroups*numGroups) != scatteringOrder+1) {
    std::cerr << "Size error in scattering vector" << std::endl;
    return;
  }
  for (unsigned int order=0; order<=scatteringOrder; order++) {
    for (unsigned int groupOut=1; groupOut<=numGroups; groupOut++) {
      for (unsigned int groupIn=1; groupIn<=numGroups; groupIn++) {
        _sigma_s_ker[ si(groupIn, groupOut, order) ]
          = crossSection[ (groupOut-1)*numGroups + groupIn - 1];
        if (order==0) _sigma_s[groupIn-1] += crossSection[ (groupOut-1)*numGroups + groupIn - 1];
      }
    }
  }
}

void
Material::putFissionSpectrum(unsigned int group, double fraction)
{
  _fissionSpectrum[group-1] = fraction;
}

void
Material::putFissionSpectrum(std::vector<double>& spectrum)
{
  for (int g=0; g<numGroups; g++) {
    _fissionSpectrum[g] = spectrum[g];
  }
}


void
Material::putNeutronSpeed(unsigned int group, double speed)
{
  _speed[group-1] = speed;
}

void
Material::putNeutronSpeed(std::vector<double>& speed)
{
  for (int g=0; g<numGroups; g++) {
    _speed[g] = speed[g];
  }
}


MaterialFactory::~MaterialFactory()
{
  for (std::map<std::string, Material*>::iterator it=_materialMap.begin(); it!=_materialMap.end(); ++it) {
    delete it->second;
    it->second = NULL;
  }
  _materialMap.clear();
  _materialID.clear();

}

int
MaterialFactory::createMaterialsFromInput(InputParser& input)
{
  // Make materials
  std::vector<std::string> path(1,"Material");
  for (int i=1; i<=input.countSets(path); i++) {
    LOG_DBG("Reading material ", i);
    readNextMaterial(input, i);
  }

  return numGroups;
}

void
MaterialFactory::readNextMaterial(InputParser& input, int matNum)
{

  std::vector<std::string> path(1,"Material");
  std::vector<double> v;

  v = input.getVector(path, "numGroups", matNum);
  if (v.size() != 1) {
    LOG_ERR("Must specify a scalar number of groups.");
  }
  int tempNumGroups = int (v[0]);
  numGroups = tempNumGroups;

  std::string matName = input.getString(path, "name", matNum);
  if (matName == "empty") {
    LOG_ERR("Must specify the name of all materials; use the name parameter.");
  }

  v = input.getVector(path, "block", matNum);
  if (v.size() == 0) {
    LOG_ERR("Must specify a block number for each material.  Use the block parameter.");
  }
  int block = int (v[0]);

  Material* mat1 = createNewMaterial(matName, (unsigned int)block,  (unsigned int)numGroups, 0 );

  
  v = input.getVector(path, "sigmaT", matNum);
  if (v.size() == 0) {
    LOG_ERR("Must specify a total cross section for each material; use the sigmaT parameter.");
  }
  mat1->putSigma_t( v );

  v = input.getVector(path, "sigmaS", matNum);
  if (v.size() == 0) {
    LOG_ERR("Must specify a scattering cross section for each material; use the sigmaS parameter.");
  }
  mat1->putSigma_s( v );

  // Read optional fission XS
  v = input.getVector(path, "nu_sigmaF", matNum);
  if (v.size() > 0) {
    mat1->putNuSigma_f( v );

    v = input.getVector(path, "fissionSpectrum", matNum);
    if (v.size() == 0) LOG_ERR(matName, " is fissile; must specify fission spectrum.");
    mat1->putFissionSpectrum( v );

    // Read optional delayed neutron data
    v = input.getVector(path, "dnFractionEff", matNum);
    if ( v.size()==1 ) {
      mat1->betaEff = v[0];

      v = input.getVector(path, "dnDecayConst", matNum);
      mat1->numDelayedGroups = v.size();
      for (int i=0; i<v.size(); i++) {
        mat1->lambda[i] = v[i];
      }
    
      v = input.getVector(path, "dnFraction", matNum);
      for (int i=0; i<v.size(); i++) {
        mat1->beta[i] = v[i];
      }
    }
    else {
      mat1->betaEff = 0.0;
    }
  }

  v = input.getVector(path, "speed", matNum);
  if (v.size() > 0) {
    mat1->putNeutronSpeed( v );
  }

}


Material*
MaterialFactory::createNewMaterial(std::string& materialName, unsigned int block, unsigned int G, unsigned int N)
{
  Material* mat = getMaterial(materialName);
  if ( mat ) {
    std::cerr << "****** Material name collision" << std::endl;
    return mat;
  }

  mat = new Material(G, N);
  _materialMap[materialName] = mat;
  _materialID[block] = mat;
  return mat;
}


Material*
MaterialFactory::getMaterial(std::string& materialName)
{
  std::map<std::string, Material*>::iterator _matIt;
  _matIt = _materialMap.find( materialName );
  if ( _matIt == _materialMap.end() )
    {
      return NULL;
    }
  else
    {
      return _matIt->second;
    }

}

Material*
MaterialFactory::getMaterial(unsigned int& materialID)
{
  std::map<unsigned int, Material*>::iterator _matIt;
  _matIt = _materialID.find( materialID );
  if ( _matIt == _materialID.end() )
    {
      return NULL;
    }
  else
    {
      return _matIt->second;
    }

}


void
MaterialFactory::list()
{
  std::cout << std::endl << "Materials generated" << std::endl;
  for (std::map<std::string, Material*>::iterator it=_materialMap.begin(); it!=_materialMap.end(); ++it) {
    std::cout << "  " << it->first << std::endl;
  }
}
