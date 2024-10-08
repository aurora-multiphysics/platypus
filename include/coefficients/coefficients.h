#pragma once
#include "mesh_extras.hpp"
#include "MFEMContainers.h"
#include <fstream>
#include <iostream>
#include <memory>
#include <unordered_set>

namespace platypus
{
double prodFunc(double a, double b);
double fracFunc(double a, double b);

class Subdomain
{
public:
  Subdomain(std::string name_, int id_);

  std::string _name;
  int _id;
  platypus::NamedFieldsMap<mfem::Coefficient> _scalar_coefficients;
  platypus::NamedFieldsMap<mfem::VectorCoefficient> _vector_coefficients;
};

// platypus::Coefficients - stores all scalar and vector coefficients
//--SetTime
//--scalars
//--vectors

// Stores all coefficients defined over
class Coefficients
{
  double _t; // Time at which time-dependent coefficients are evaluated
public:
  Coefficients();
  ~Coefficients() = default;

  Coefficients(std::vector<Subdomain> subdomains_);

  void SetTime(double t);
  void AddGlobalCoefficientsFromSubdomains();
  void RegisterDefaultCoefficients();

  platypus::NamedFieldsMap<mfem::Coefficient> _scalars;
  platypus::NamedFieldsMap<mfem::VectorCoefficient> _vectors;
  std::vector<Subdomain> _subdomains;
};
} // namespace platypus