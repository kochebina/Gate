/*----------------------
   Copyright (C): OpenGATE Collaboration

This software is distributed under the terms
of the GNU Lesser General  Public Licence (LGPL)
See LICENSE.md for further details
----------------------*/


#include "GateVCoincidencePulseProcessor.hh"

#include "GateTools.hh"
#include "GateCoincidencePulseProcessorChain.hh"
#include "GateCoincidenceDigiMaker.hh"
#include "GateDigitizerOld.hh"
//------------------------------------------------------------------------------------------------------
// Constructs a new pulse-processor attached to a GateDigitizerOld
GateVCoincidencePulseProcessor::GateVCoincidencePulseProcessor(GateCoincidencePulseProcessorChain* itsChain,
      	      	      	   const G4String& itsName) 
    : GateClockDependent(itsName),
/*      m_isTriCoincProc(0), //mhadi_add*/
      m_chain(itsChain)
{
  GateDigitizerOld* digitizer = GateDigitizerOld::GetInstance();

  G4cout << " in GateVCoincidencePulseProcessor call new GateCoincidenceDigiMaker\n";
  digitizer->InsertDigiMakerModule( new GateCoincidenceDigiMaker(digitizer, itsName,false) );
}
//------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------
// Method overloading GateClockDependent::Describe()
// Print-out a description of the component
// Calls the pure virtual method DecribeMyself()
void GateVCoincidencePulseProcessor::Describe(size_t indent) 
{
  GateClockDependent::Describe(indent);
  G4cout << GateTools::Indent(indent) << "Attached to:        '" << GetChain()->GetObjectName() << "'\n";
  G4cout << GateTools::Indent(indent) << "Output:             '" << GetObjectName() << "'\n";
  DescribeMyself(indent);
}
//------------------------------------------------------------------------------------------------------    
