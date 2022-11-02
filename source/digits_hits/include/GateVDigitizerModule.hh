/*----------------------
   Copyright (C): OpenGATE Collaboration

This software is distributed under the terms
of the GNU Lesser General  Public Licence (LGPL)
See LICENSE.md for further details
----------------------*/

// OK GND 2022

#ifndef GateVDigitizerModule_h
#define GateVDigitizerModule_h 1

#include "G4VDigitizerModule.hh"
#include "GateDigi.hh"
#include "GateDigitizer.hh"
#include "GateClockDependent.hh"
#include "GateCrystalSD.hh"

#include "globals.hh"

class GateDigitizer;

class GateVDigitizerModule : public G4VDigitizerModule, public GateClockDependent
{
public:
  
  GateVDigitizerModule(G4String DMname, G4String path, GateDigitizer *digitizer, GateCrystalSD* SD);
  GateVDigitizerModule(G4String DMname);
  virtual ~GateVDigitizerModule();
  
  //void SetInputDigiCollection();

  virtual void Digitize()=0;
  void InputCollectionID();




  //! Method overloading GateClockDependent::Describe()
  //! Print-out a description of the component
  //! Calls the pure virtual method DecribeMyself()
  virtual void Describe(size_t indent=0);

  //! Pure virtual method DecribeMyself()
  virtual void DescribeMyself(size_t indent=0) =0 ;

  inline GateDigitizer* GetDigitizer()
    { return m_digitizer; }

private:
 /* GateDigi* m_outputDigi;

  GateVDigitizerModuleMessenger *fMessenger;

  GateDigiCollection*  OutputDigiCollection;
*/
  GateDigitizer *m_digitizer;


protected:
  GateCrystalSD *m_SD;
  G4int m_outputDCID;
  G4int	m_InitDMID;

  GateDigiCollection* m_IDC;
  G4int m_DCID;
};

#endif








