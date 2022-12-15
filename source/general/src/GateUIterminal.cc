/*----------------------
   Copyright (C): OpenGATE Collaboration

This software is distributed under the terms
of the GNU Lesser General  Public Licence (LGPL)
See LICENSE.md for further details
----------------------*/

//LF
//#include <strstream>
#include <sstream>
//LF
#include "GateUIterminal.hh"


G4int GateUIterminal::ReceiveG4cout(const G4String& coutString)
{
  if ( coutString == "G4PhysicalVolumeModel::Validate() called.\n" )
    return 0;
  if ( coutString == "  Volume of the same name and copy number (\"world_P\", copy 0) still exists and is being used.\n" )
    return 0;
  if ( coutString == "  Be warned that this does not necessarily guarantee it's the same\n" )
    return 0;
  if ( coutString == "  volume you originally specified in /vis/scene/add/.\n" )
    return 0;

  std::cout << coutString << std::flush;
  return 0;
}


G4int GateUIterminal::ReceiveG4cerr(const G4String& cerrString)
{
	std::cerr << "[G4-cerr] " <<  cerrString << std::flush;
	// Check if this error is 'command not found' (or related) to stop Gate
    bool isMacroError = false;
    std::string::size_type i = cerrString.find("***** COMMAND NOT FOUND <", 0);
    isMacroError = isMacroError || (i != std::string::npos);
    i = cerrString.find("***** Illegal application state <", 0);
    isMacroError = isMacroError || (i != std::string::npos);
    i = cerrString.find("***** Illegal parameter (", 0);
    isMacroError = isMacroError || (i != std::string::npos);
    i = cerrString.find("***** Can not open a macro file <", 0);
    isMacroError = isMacroError || (i != std::string::npos);
    i = cerrString.find("ERROR: Can not open a macro file <", 0);
    isMacroError = isMacroError || (i != std::string::npos);

    if (isMacroError) {
      std::cerr << "[Gate] Sorry, error in a macro command : abort.\n";
      i = cerrString.find("digitzer", 0);
      bool isDigitizerMacroError = isMacroError || (i != std::string::npos);
      if (isDigitizerMacroError)
      	  {
    	  std::cerr << "------------------------ [GateSinglesDigitizer Messenger Problem]------------------------ \n";
    	  std::cerr << "Probably you try to use a command for the old version of the digitizer.\n";
    	  std::cerr << "Try </gate/digitizerMgr> commands instead, you can find the documentation here: XXXX \n"; //TODO insert the link to the documentation page
    	  std::cerr << "---------------------------------------------------------------------------------- \n";
      	  }
      exit(-1);
    }

  return 0;
}

