#pragma once

#include "MooseApp.h"

class PlatypusTestApp : public MooseApp
{
public:
  static InputParameters validParams();

  PlatypusTestApp(InputParameters parameters);
  virtual ~PlatypusTestApp();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs = false);
};
