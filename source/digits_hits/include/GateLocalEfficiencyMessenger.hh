/*----------------------
   Copyright (C): OpenGATE Collaboration

This software is distributed under the terms
of the GNU Lesser General  Public Licence (LGPL)
See LICENSE.md for further details
----------------------*/
//GND:ClassToRemove
#ifndef GateLocalEfficiencyMessenger_h
#define GateLocalEfficiencyMessenger_h 1

#include "GatePulseProcessorMessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

class GateLocalEfficiency;

class GateLocalEfficiencyMessenger: public GatePulseProcessorMessenger
{
  public:
    GateLocalEfficiencyMessenger(GateLocalEfficiency* itsPulseProcessor);
    virtual ~GateLocalEfficiencyMessenger();

    inline void SetNewValue(G4UIcommand* aCommand, G4String aString);

  private:
    G4UIcmdWithAString   *crystalEfficiencyCmd;
    G4UIcmdWithAnInteger *enableCommand;
    G4UIcmdWithAnInteger *disableCommand;
};

#endif
