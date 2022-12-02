/*----------------------
   Copyright (C): OpenGATE Collaboration

This software is distributed under the terms
of the GNU Lesser General  Public Licence (LGPL)
See LICENSE.md for further details
----------------------*/
//GND:ClassToRemove

#include "GateLocalEfficiency.hh"

#include "GateLocalEfficiencyMessenger.hh"
#include "GateTools.hh"
#include "GateVSystem.hh"
#include "GateVDistribution.hh"
#include "GateSystemListManager.hh"

#include "Randomize.hh"

#include "G4UnitsTable.hh"
#include <fstream>



GateLocalEfficiency::GateLocalEfficiency(GatePulseProcessorChain* itsChain,
      	      	      	      	 const G4String& itsName)
  : GateVPulseProcessor(itsChain,itsName),
    m_enabled(),
    m_efficiency(0)
{
  m_messenger = new GateLocalEfficiencyMessenger(this);
}




GateLocalEfficiency::~GateLocalEfficiency()
{
  delete m_messenger;
}

void GateLocalEfficiency::ComputeSizes()
{
   static size_t firstPass=true;
   static size_t depth=0;
   static size_t totSize=0;
   static GateVSystem* system;

   if (firstPass){
      firstPass=false;
      system= GateSystemListManager::GetInstance()->GetSystem(0);
      if (!system){
      	 G4cerr<<"[GateLocalEfficiency::ComputeSizes] Problem : no system defined\n";
	 return;
      }
      depth = system->GetTreeDepth();

      G4cout<<m_enabled.size()<<" "<< depth<<G4endl;
      if (m_enabled.size() != depth) {
	 G4cerr<<"[GateLocalEfficiency::ProcessOnePulse]Warning : enabling vector size modified (from "<<m_enabled.size()
	       <<" to "<<depth<<") and set all entries to 0\n";
	 m_enabled.resize(depth);
	 for (std::vector<G4bool>::iterator itr=m_enabled.begin();itr!=m_enabled.end();itr++) *itr=false;
      }
      totSize = system->ComputeNofSubCrystalsAtLevel(0,m_enabled);
   }
   if (m_enabled.size() != depth) {
      G4cerr<<"[GateLocalEfficiency::ProcessOnePulse]Warning : enabling vector size modified (from "<<m_enabled.size()
	    <<" to "<<depth<<") and set all entries to 0\n";
      m_enabled.resize(depth);
      for (std::vector<G4bool>::iterator itr=m_enabled.begin();itr!=m_enabled.end();itr++) *itr=false;
   }
   if (m_efficiency->MaxX() < totSize-1){
      G4cerr<<"[GateLocalEfficiency::ProcessOnePulse]Warning : efficiency table size's wrong ("<<m_efficiency->MaxX()
	    <<" instead of "<<totSize<<") disabling efficiency (all set to 1)\n";
      m_efficiency=0;
   }
}
void GateLocalEfficiency::SetMode(size_t i,G4bool val)
{
   GateVSystem* system;
   size_t depth=0;
   system= GateSystemListManager::GetInstance()->GetSystem(0);
   if (!system){
      G4cerr<<"[GateLocalEfficiency::ComputeSizes] Problem : no system defined\n";
      return;
   }
   depth = system->GetTreeDepth();
   if (m_enabled.size() != depth) {
      G4cerr<<"[GateLocalEfficiency::ProcessOnePulse]Warning : enabling vector size modified (from "<<m_enabled.size()
	    <<" to "<<depth<<") and set all entries to 0\n";
      m_enabled.resize(depth);
      for (std::vector<G4bool>::iterator itr=m_enabled.begin();itr!=m_enabled.end();itr++) *itr=false;
   }

   if (i<m_enabled.size()){
      m_enabled[i]=val;
   } else {
      G4cerr<<"[GateLocalEfficiency::SetMode] WARNING : index outside limits ("
            <<i<<">"<<m_enabled.size()<<")\n";
   }
}
void GateLocalEfficiency::ProcessOnePulse(const GatePulse* inputPulse,GatePulseList& outputPulseList)
{
   if (!m_efficiency){
      outputPulseList.push_back(new GatePulse(*inputPulse));
      return;
   }
   GateVSystem* system = GateSystemListManager::GetInstance()->GetSystem(0);
   if (!system){
      G4cerr<<"[GateLocalEfficiency::ProcessOnePulse] Problem : no system defined\n";
      return ;
   }

   static bool firstPass=true;
      if (firstPass && m_enabled.empty())  {
         ComputeSizes();
      }


   size_t ligne = system->ComputeIdFromVolID(inputPulse->GetOutputVolumeID(),m_enabled);
   G4double eff = m_efficiency->Value(ligne);

   if (CLHEP::RandFlat::shoot(0.,1.) < eff)
      outputPulseList.push_back(new GatePulse(*inputPulse));
}

void GateLocalEfficiency::DescribeMyself(size_t indent)
{
  G4cout << GateTools::Indent(indent) << "Tabular Efficiency "<< Gateendl;
}
