/*----------------------
   Copyright (C): OpenGATE Collaboration

This software is distributed under the terms
of the GNU Lesser General  Public Licence (LGPL)
See LICENSE.md for further details
----------------------*/


#include "GateCoincidenceSorterMessenger.hh"

#include "GateCoincidenceSorter.hh"
//#include "GateSystemListManager.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"

GateCoincidenceSorterMessenger::GateCoincidenceSorterMessenger(GateCoincidenceSorter* itsCoincidenceSorter)
    : GateClockDependentMessenger(itsCoincidenceSorter), m_CoincidenceSorter(itsCoincidenceSorter)
{
	 G4String DirectoryName = "/gate/digitizerMgr/CoincidenceSorter/"+m_CoincidenceSorter->GetOutputName()+"/";

	// G4cout<<"GateCoincidenceSorterMessenger constr " << itsCoincidenceSorter<<" "<< this<<G4endl;
  G4String guidance;
  G4String cmdName;


  cmdName = DirectoryName + "setWindow";
  windowCmd = new G4UIcmdWithADoubleAndUnit(cmdName,this);
  windowCmd->SetGuidance("Set time-window for coincidence");
  windowCmd->SetUnitCategory("Time");

  cmdName = DirectoryName + "setWindowJitter";
  windowJitterCmd = new G4UIcmdWithADoubleAndUnit(cmdName,this);
  windowJitterCmd->SetGuidance("Set standard deviation of window jitter");
  windowJitterCmd->SetUnitCategory("Time");

  cmdName = DirectoryName + "setOffset";
  offsetCmd = new G4UIcmdWithADoubleAndUnit(cmdName,this);
  offsetCmd->SetGuidance("Set time offset for delay coincidences");
  offsetCmd->SetUnitCategory("Time");

  cmdName = DirectoryName + "setOffsetJitter";
  offsetJitterCmd = new G4UIcmdWithADoubleAndUnit(cmdName,this);
  offsetJitterCmd->SetGuidance("Set standard deviation of offset jitter");
  offsetJitterCmd->SetUnitCategory("Time");

  cmdName = DirectoryName+"minSectorDifference";
  minSectorDiffCmd = new G4UIcmdWithAnInteger(cmdName.c_str(),this);
  minSectorDiffCmd->SetGuidance("Set the minimum sector difference for valid coincidences.");
  minSectorDiffCmd->SetParameterName("diff",false);
  minSectorDiffCmd->SetRange("diff>=1");

  cmdName = DirectoryName+"setDepth";
  setDepthCmd = new G4UIcmdWithAnInteger(cmdName.c_str(),this);
  setDepthCmd->SetGuidance("Set the depth of system-level for coincidences.");
  setDepthCmd->SetParameterName("depth",false);
  setDepthCmd->SetRange("depth>=1");

  cmdName = DirectoryName+"setPresortBufferSize";
  setPresortBufferSizeCmd = new G4UIcmdWithAnInteger(cmdName.c_str(),this);
  setPresortBufferSizeCmd->SetGuidance("Set the size of the presort buffer.");
  setPresortBufferSizeCmd->SetParameterName("size",false);
  setPresortBufferSizeCmd->SetRange("size>=32");

  cmdName = DirectoryName+"setInputCollection";
  SetInputNameCmd = new G4UIcmdWithAString(cmdName,this);
  SetInputNameCmd->SetGuidance("Set the name of the input digi channel");
  SetInputNameCmd->SetParameterName("Name",false);

  cmdName = DirectoryName+"MultiplesPolicy";
  MultiplePolicyCmd = new G4UIcmdWithAString(cmdName,this);
  MultiplePolicyCmd->SetGuidance("How to treat multiples coincidences");
  MultiplePolicyCmd->SetCandidates("killAll takeAllGoods killAllIfMultipleGoods takeWinnerOfGoods takeWinnerIfIsGood takeWinnerIfAllAreGoods keepAll keepIfAnyIsGood keepIfOnlyOneGood keepIfAllAreGoods");

  cmdName = DirectoryName+"allDigiOpenCoincGate";
  AllDigiOpenCoincGateCmd = new G4UIcmdWithABool(cmdName,this);
  AllDigiOpenCoincGateCmd->SetGuidance("Specify if a given digi can be part of two coincs");

  //For CC module
  cmdName = DirectoryName+"setTriggerOnlyByAbsorber";
  SetTriggerOnlyByAbsorberCmd = new G4UIcmdWithABool(cmdName,this);
  SetTriggerOnlyByAbsorberCmd->SetGuidance("Specify if only the digis in the absorber can open a coincidencee window");

  cmdName = DirectoryName+"setAcceptancePolicy4CC";
  SetAcceptancePolicy4CCCmd = new G4UIcmdWithAString(cmdName,this);
  SetAcceptancePolicy4CCCmd ->SetGuidance("Coincidence acceptance policy in CC");
  SetAcceptancePolicy4CCCmd ->SetCandidates("keepIfMultipleVolumeIDsInvolved keepIfMultipleVolumeNamesInvolved keepAll");

  cmdName = DirectoryName+"setEventIDCoinc";
  SetEventIDCoincCmd = new G4UIcmdWithABool(cmdName,this);
  SetEventIDCoincCmd->SetGuidance("Set to one for event identification coincidencences");


}


GateCoincidenceSorterMessenger::~GateCoincidenceSorterMessenger()
{
    delete windowCmd;
    delete offsetCmd;
    delete windowJitterCmd;
    delete offsetJitterCmd;
    delete minSectorDiffCmd;
    delete SetInputNameCmd;
    delete MultiplePolicyCmd;
    delete setPresortBufferSizeCmd;
    delete AllDigiOpenCoincGateCmd;
    delete SetTriggerOnlyByAbsorberCmd;
    delete SetAcceptancePolicy4CCCmd;
    delete SetEventIDCoincCmd;

}


void GateCoincidenceSorterMessenger::SetNewValue(G4UIcommand* aCommand, G4String newValue)
{
  if ( aCommand==windowCmd )
    { m_CoincidenceSorter->SetWindow(windowCmd->GetNewDoubleValue(newValue)); }
  else if( aCommand == windowJitterCmd )
    { m_CoincidenceSorter->SetWindowJitter(windowJitterCmd->GetNewDoubleValue(newValue)); }
  else if( aCommand == offsetCmd )
    { m_CoincidenceSorter->SetOffset(offsetCmd->GetNewDoubleValue(newValue)); }
  else if( aCommand == offsetJitterCmd )
    { m_CoincidenceSorter->SetOffsetJitter(offsetJitterCmd->GetNewDoubleValue(newValue)); }
  else if( aCommand == minSectorDiffCmd )
    { m_CoincidenceSorter->SetMinSectorDifference(minSectorDiffCmd->GetNewIntValue(newValue)); }
  else if( aCommand == setDepthCmd )
    { m_CoincidenceSorter->SetDepth(setDepthCmd->GetNewIntValue(newValue)); }
  else if( aCommand == setPresortBufferSizeCmd )
    { m_CoincidenceSorter->SetPresortBufferSize(setPresortBufferSizeCmd->GetNewIntValue(newValue)); }
  else if (aCommand == SetInputNameCmd)
    {
     m_CoincidenceSorter->SetInputName(newValue);
     m_CoincidenceSorter->SetSystem(newValue); //! Attach to the suitable system from the digitizer m_systemList (multi-system approach)
    }
  else if (aCommand == MultiplePolicyCmd)
    { m_CoincidenceSorter->SetMultiplesPolicy(newValue); }
  else if (aCommand == SetAcceptancePolicy4CCCmd)
    { m_CoincidenceSorter->SetAcceptancePolicy4CC(newValue); }
  else if (aCommand == AllDigiOpenCoincGateCmd)
    { m_CoincidenceSorter->SetAllDigiOpenCoincGate(AllDigiOpenCoincGateCmd->GetNewBoolValue(newValue)); }
  else if (aCommand == SetTriggerOnlyByAbsorberCmd)
    { m_CoincidenceSorter->SetIfTriggerOnlyByAbsorber(SetTriggerOnlyByAbsorberCmd->GetNewBoolValue(newValue));}
  else if (aCommand == SetEventIDCoincCmd)
    { m_CoincidenceSorter->SetIfEventIDCoinc(SetEventIDCoincCmd->GetNewBoolValue(newValue));}
  else
    GateClockDependentMessenger::SetNewValue(aCommand,newValue);

}
