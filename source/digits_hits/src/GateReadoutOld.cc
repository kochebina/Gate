/*----------------------
   Copyright (C): OpenGATE Collaboration

This software is distributed under the terms
of the GNU Lesser General  Public Licence (LGPL)
See LICENSE.md for further details
----------------------*/


#include "GateReadoutOld.hh"

#include "G4UnitsTable.hh"

#include "GateOutputVolumeID.hh"
#include "GateReadoutOldMessenger.hh"
#include "GateTools.hh"
#include "GateDigitizerOld.hh"
#include "GateArrayComponent.hh"
#include "GateVSystem.hh"

//class GateVSystem;
/*
  S. Stute - June 2014: complete redesign of the ReadoutOld module and add a new policy to emulate PMT.
    - Fix bug in choosing the maximum energy pulse.We now have some temporary lists of variables to deal
      with the output pulses. These output pulses are only created at the end of the method. In previous
      versions, the output pulse was accumulating energy as long as input pulses were merged together, but
      the problem is that the comparison of energy was done between the input pulse and this output pulse
      with increasing energy. So with more than 2 pulses to be merged together, the behaviour was undefined.
    - Move all the processing into the upper method ProcessPulseList instead of using the mother version
      working into the ProcessOnePulse method. Thus the ProcessOnePulse in this class is not used anymore.
    - Create policy choice: now we can choose via the messenger between EnergyWinner and EnergyCentroid.
    - For the EnergyCentroid policy, the centroid position is computed using the crystal indices in each
      direction, doing the computation with floating point numbers, and then casting the result into
      integer indices. Using that method, we ensure the centroid position to be in a crystal (if we work
      with global position, we can fall between two crystals in the presence of gaps).
      The depth is ignored with this strategy; it is forced to be at one level above the 'crystal' level.
      If there a 'layer' level below the 'crystal' level, an energy winner strategy is adopted.

  O. Kochebina - April 2022: new messenger options are added and some minor bugs corrected
*/

GateReadoutOld::GateReadoutOld(GatePulseProcessorChain* itsChain,
      	      	      	 const G4String& itsName)
  : GateVPulseProcessor(itsChain,itsName),
    m_depth(0),
	m_policy("TakeEnergyWinner"),
	m_IsFirstEntrance(1)
{

	//G4cout<<"GateReadoutOld::GateReadoutOld"<<Gateendl;

	m_messenger = new GateReadoutOldMessenger(this);
  // S. Stute: These variables are used for the energy centroid strategy
  m_nbCrystalsX  = -1;
  m_nbCrystalsY  = -1;
  m_nbCrystalsZ  = -1;
  m_nbCrystalsXY = -1;
  m_crystalDepth = -1;
  m_systemDepth  = -1;
  m_system = NULL;
  m_crystalComponent = NULL;

  //G4cout<<"GateReadoutOld::GateReadoutOld "<< m_policy<<Gateendl;

}

GateReadoutOld::~GateReadoutOld()
{
  delete m_messenger;
}

void GateReadoutOld::SetReadoutOldParameters()
{


	//checking the if depth or ReadoutOldVolumeName are defined and that only one is set.
	if(!m_volumeName.empty() && m_depth!=0)
		GateError("***ERROR*** You can choose ReadoutOld parameter either with /setDepth OR /setReadoutOldVolume!");

	 //////////////DEPTH SETTING/////////
	 //set the previously default value for compatibility of users macros
	if(m_volumeName.empty()  && m_depth==0)
		m_depth=1; //previously default value

	 //set m_depth according user defined volume name
	 if(!m_volumeName.empty()) //only for EnergyWinner
	 	 {
		 if(m_policy =="TakeEnergyCentroid"&& !m_IsForcedDepthCentroid)
			 GateError("***ERROR*** Please, remove /setDepth or /setReadoutOldVolume for TakeEnergyCentroid policy as this parameter is set automatically. "
					 "Use /forceReadoutOldVolumeForEnergyCentroid flag if you still want to set your depth/volume for ReadoutOld.\n");



		 	 GateVSystem* m_system = this->GetChain()->GetSystem();
		 	 if (m_system==NULL) G4Exception( "GateReadoutOld::SetReadoutOldParameters", "SetReadoutOldParameters", FatalException,
	  	  	                                   "Failed to get the system corresponding to that processor chain. Abort.\n");

		 	 m_systemDepth = m_system->GetTreeDepth();

		 	 GateObjectStore* anInserterStore = GateObjectStore::GetInstance();
		 	 for (G4int i=1;i<m_systemDepth;i++)
		 	 {
		 		 GateSystemComponent* comp0= (m_system->MakeComponentListAtLevel(i))[0][0];
		 		 GateVVolume *creator = comp0->GetCreator();
		 		 GateVVolume* anInserter = anInserterStore->FindCreator(m_volumeName);

		 		 if(creator==anInserter)
		 			 m_depth=i;

	  	   }
	  	}


	 //////////////Resulting positioning SETTING/////////
	 //previously default conditions for compatibility of users macros
	/* if(m_resultingXY.empty() && m_resultingZ.empty() && m_policy =="TakeEnergyCentroid")
	 {
		 m_resultingXY="crystalCenter";
		 m_resultingZ="crystalCenter";
	 }
	 if (m_resultingXY.empty() && m_resultingZ.empty() && m_policy =="TakeEnergyWinner")
	 {
		 m_resultingXY="exactPostion";
		 m_resultingZ="exactPostion";
	 }
	 */


	if (m_policy=="TakeEnergyCentroid" && (!m_volumeName.empty()||m_depth) &&  !m_IsForcedDepthCentroid)
	 {
		GateWarning("WARNING! Commands /setDepth and /setReadoutOldVolume are ignored as Energy Centroid policy is used: "
				"the depth is forced to be at the level just above the crystal level, whatever the system used."
				"To force the depth, please, set the flag /forceReadoutOldVolumeForEnergyCentroid to true");
	 }

	if (m_policy=="TakeEnergyWinner" && m_IsForcedDepthCentroid)
		 {
		GateError("***ERROR*** Command /forceReadoutOldVolumeForEnergyCentroid can not be used for Winner policy. Abort.\n");
		 }


	if (m_policy=="TakeEnergyCentroid" )
		{

			// Find useful stuff for centroid based computation
			//m_policy = "TakeEnergyCentroid";
			// Get the system
			GateVSystem* m_system = this->GetChain()->GetSystem();
			if (m_system==NULL) G4Exception( "GateReadoutOld::ProcessPulseList", "ProcessPulseList", FatalException,
					"Failed to get the system corresponding to that processor chain. Abort.\n");
			// Get the array component corresponding to the crystal level using the name 'crystal'
			GateArrayComponent* m_crystalComponent = m_system->FindArrayComponent("crystal");
			if (m_crystalComponent==NULL) G4Exception( "GateReadoutOld::ProcessPulseList", "ProcessPulseList", FatalException,
												  "Failed to get the array component corresponding to the crystal. Abort.\n");
			// Get the number of crystals in each direction
			m_nbCrystalsZ  = m_crystalComponent->GetRepeatNumber(2);
			m_nbCrystalsY  = m_crystalComponent->GetRepeatNumber(1);
			m_nbCrystalsX  = m_crystalComponent->GetRepeatNumber(0);
			m_nbCrystalsXY = m_nbCrystalsX * m_nbCrystalsY;
			if (m_nbCrystalsX<1 || m_nbCrystalsY<1 || m_nbCrystalsZ<1)
				G4Exception( "GateReadoutOld::ProcessPulseList", "ProcessPulseList", FatalException,
						"Crystal repeater numbers are wrong !\n");
			//G4cout << "[" << GetObjectName() << "] -> Found crystal array with associated repeater: ["
			//       << m_nbCrystalsX << ";" << m_nbCrystalsY << ";" << m_nbCrystalsZ << "]\n";
			// Get tree depth of the system
			m_systemDepth = m_system->GetTreeDepth();
			//G4cout << "  Depth of the system: " << m_systemDepth << Gateendl;
			// Find the crystal depth in the system
			GateSystemComponent* this_component = m_system->GetBaseComponent();
			m_crystalDepth = 0;
			while (this_component!=m_crystalComponent && m_crystalDepth+1<m_systemDepth)
			{
				this_component = this_component->GetChildComponent(0);
				m_crystalDepth++;
			}
			if (this_component!=m_crystalComponent) G4Exception( "GateReadoutOld::ProcessPulseList", "ProcessPulseList", FatalException,
																	"Failed to get the system depth corresponding to the crystal. Abort.\n");
			// Now force m_depth to be right above the crystal depth
			//m_depth = m_crystalDepth - 1;
		if (!m_IsForcedDepthCentroid)
			{
			m_depth = m_crystalDepth - 1;
			}

		}



	if (m_policy!="TakeEnergyCentroid" && m_policy!="TakeEnergyWinner")
		G4Exception( "GateReadoutOld::SetPolicy", "SetPolicy", FatalException, "Unknown provided policy, please see the guidance. Abort.\n");

	//G4cout<<"Policy = "<< m_policy<< Gateendl;
	//G4cout<<"Depth =  "<< m_depth<<Gateendl;
	//G4cout<<"resultingXY = "<< m_resultingXY<<Gateendl;
	//G4cout<<"reulstingZ = "<< m_resultingZ<<Gateendl;

}

// S. Stute: This function is virtual (but not pure) in the mother GateVPulseProcessor.
//           We now overload it in order to be able to process all the pulses at a same time (not using the ProcessOnePulse method).
//           This is only to be able to program the option EnergyCentroid based on crystal indices.
//           So the ProcessOnePulse method right after this one is not used anymore.
//           Note: this new strategy also allowed us to correct the bug for the energyWinner policy.
GatePulseList* GateReadoutOld::ProcessPulseList(const GatePulseList* inputPulseList)
{
	if(m_IsFirstEntrance) //set parameters at the first iteration
	{
		SetReadoutOldParameters();
		m_IsFirstEntrance=0;
	}

  if (!inputPulseList)
    return 0;

  size_t n_pulses = inputPulseList->size();

  if (nVerboseLevel>1)
    G4cout << "[" << GetObjectName() << "::ProcessPulseList]: processing input list with " << n_pulses << " entries\n";
  if (!n_pulses)
    return 0;

  GatePulseList* outputPulseList = new GatePulseList(GetObjectName());
  //G4cout<<"Policy "<< m_policy<< " "<< m_depth<<" "<< m_resultingXY<<" "<< m_resultingZ<<Gateendl;

  // S. Stute: these variables are used for the energy centroid policy
  G4double* final_time = NULL;
  G4double* final_crystal_posX = NULL;
  G4double* final_crystal_posY = NULL;
  G4double* final_crystal_posZ = NULL;
  G4double* final_energy = NULL;
  G4int* final_nb_pulses = NULL;
  GatePulse** final_pulses = NULL;
  if (m_policy=="TakeEnergyCentroid")
  {
    final_time         = (G4double*)calloc(n_pulses,sizeof(G4double));
    final_crystal_posX = (G4double*)calloc(n_pulses,sizeof(G4double));
    final_crystal_posY = (G4double*)calloc(n_pulses,sizeof(G4double));
    final_crystal_posZ = (G4double*)calloc(n_pulses,sizeof(G4double));
    final_nb_pulses    = (G4int*)calloc(n_pulses,sizeof(G4int));
  }
  // S. Stute: we need energy to sum up correctly for all output pulses and affect only at the end.
  //           In previous versions, even for take Winner, the energy was affected online so the
  //           final pulse was not the winner in all cases.
  final_energy = (G4double*)calloc(n_pulses,sizeof(G4double));
  final_pulses = (GatePulse**)calloc(n_pulses,sizeof(GatePulse*));
  G4int final_nb_out_pulses = 0;

  // Start loop on input pulses
  GatePulseConstIterator iterIn;
  for (iterIn = inputPulseList->begin() ; iterIn != inputPulseList->end() ; ++iterIn)
  {
    GatePulse* inputPulse = *iterIn;

    const GateOutputVolumeID& blockID = inputPulse->GetOutputVolumeID().Top(m_depth);

    if (blockID.IsInvalid())
    {
     if (nVerboseLevel>1)
        G4cout << "[GateReadoutOld::ProcessOnePulse]: out-of-block hit for \n"
               <<  *inputPulse << Gateendl
               << " -> pulse ignored\n\n";
      continue;
    }

    // Loop inside the temporary output list to see if we have one pulse with same blockID as input
    int this_output_pulse = 0;
    for (this_output_pulse=0; this_output_pulse<final_nb_out_pulses; this_output_pulse++)
      if (final_pulses[this_output_pulse]->GetOutputVolumeID().Top(m_depth) == blockID) break;

    // Case: we found an output pulse with same blockID
    if ( this_output_pulse!=final_nb_out_pulses )
    {
      // --------------------------------------------------------------------------------
      // WinnerTakeAllPolicy (APD like)
      // --------------------------------------------------------------------------------
      if (m_policy=="TakeEnergyWinner")
      {
        // If energy is higher then replace the pulse by the new one
        if ( inputPulse->GetEnergy() > final_pulses[this_output_pulse]->GetEnergy() ) final_pulses[this_output_pulse] = inputPulse;
        final_energy[this_output_pulse] += inputPulse->GetEnergy();
      }
      // --------------------------------------------------------------------------------
      // EnergyCentroidPolicy1 (like block PMT)
      // --------------------------------------------------------------------------------
      else if (m_policy=="TakeEnergyCentroid") // Crystal element is considered to be the deepest element
      {
        // First, if the energy of this pulse is higher than the previous one, take it as the reference
        // in order to have an EnergyWinner policy at levels below the crystal, if any.
        // The final energy and crystal position will be modified at the end anyway.
        if ( inputPulse->GetEnergy() > final_pulses[this_output_pulse]->GetEnergy() ) final_pulses[this_output_pulse] = inputPulse;
        // Add the energy to get the total
        G4double energy = inputPulse->GetEnergy();
        final_energy[this_output_pulse] += energy;
        // Add the time in order to compute the mean time at the end
        final_time[this_output_pulse] += inputPulse->GetTime();
        // Get the crystal ID
        int crystal_id = inputPulse->GetComponentID(m_crystalDepth);
        // Decompose the crystal_id into X, Y and Z
        int tmp_crysXY = crystal_id % m_nbCrystalsXY;
        final_crystal_posZ[this_output_pulse] += energy * (((G4double)( crystal_id / m_nbCrystalsXY ))+0.5);
        final_crystal_posY[this_output_pulse] += energy * (((G4double)( tmp_crysXY / m_nbCrystalsX  ))+0.5);
        final_crystal_posX[this_output_pulse] += energy * (((G4double)( tmp_crysXY % m_nbCrystalsX  ))+0.5);
        // Increment the number of pulses contributing to this output pulse
        final_nb_pulses[this_output_pulse]++;
      }
      else
      {
        G4Exception( "GateReadoutOld::ProcessOnePulse", "ProcessOnePulse", FatalException, "Unknown ReadoutOld policy, this is an internal error. Abort.\n");
      }
    }
    // Case: there is no output pulse with same blockID
    else
    {
      G4double energy = inputPulse->GetEnergy();
      if (m_policy=="TakeEnergyCentroid")
      {
        // Time will be averaged then
        final_time[final_nb_out_pulses] = inputPulse->GetTime();
        // Currently there is one pulse contributing to this new pulse
        final_nb_pulses[final_nb_out_pulses] = 1;
        // Get the crystal ID
        int crystal_id = inputPulse->GetComponentID(m_crystalDepth);
        // Decompose the crystal_id into X, Y and Z
        int tmp_crysXY = crystal_id % m_nbCrystalsXY;
        final_crystal_posZ[final_nb_out_pulses] = energy * (((G4double)( crystal_id / m_nbCrystalsXY ))+0.5);
        final_crystal_posY[final_nb_out_pulses] = energy * (((G4double)( tmp_crysXY / m_nbCrystalsX  ))+0.5);
        final_crystal_posX[final_nb_out_pulses] = energy * (((G4double)( tmp_crysXY % m_nbCrystalsX  ))+0.5);
      }
      // Set the current energy
      final_energy[final_nb_out_pulses] += energy;
      // Store this pulse in the list
      final_pulses[final_nb_out_pulses] = inputPulse;
      // Increment the total number of output pulses
      final_nb_out_pulses++;
    }
  } // End for input pulse

  // S. Stute: create now the output pulse list
  for (int p=0; p<final_nb_out_pulses; p++)
  {
    // Create the pulse
    GatePulse* outputPulse = new GatePulse( final_pulses[p] );
    // Affect energy
    outputPulse->SetEnergy( final_energy[p] );
    // Special affectations for centroid policy
    if (m_policy=="TakeEnergyCentroid")
    {
      // Affect time being the mean
      outputPulse->SetTime( final_time[p] / ((G4double)final_nb_pulses[p]) );
      // Compute integer crystal indices weighted by total energy
      G4int crys_posX = ((G4int)(final_crystal_posX[p]/final_energy[p]));
      G4int crys_posY = ((G4int)(final_crystal_posY[p]/final_energy[p]));
      G4int crys_posZ = ((G4int)(final_crystal_posZ[p]/final_energy[p]));
      // Compute final crystal id
      G4int crystal_id = crys_posZ * m_nbCrystalsXY + crys_posY * m_nbCrystalsX + crys_posX;
      if((outputPulse->GetVolumeID()).GetDaughterID(m_crystalDepth) == (-1)) //check to avoid Seg. fault
      {
	  GateError("Error: not all required geometry levels and sublevels for this system are defined. "
			  "(Ex.: for cylindricalPET, the required levels are: rsector, module, submodule, crystal). Please, add them to your geometry macro");
      }

      
      // We change the level of the volumeID and the outputVolumeID corresponding to the crystal with the new crystal ID
      outputPulse->ChangeVolumeIDAndOutputVolumeIDValue(m_crystalDepth,crystal_id);
      // Change coordinates (we choose here to place the coordinates at the center of the chosen crystal)
      //outputPulse->SetGlobalPos(m_system->ComputeObjectCenter(volID));
      outputPulse->ResetGlobalPos(m_system);
      outputPulse->ResetLocalPos();
    }
    if (nVerboseLevel>1)
        G4cout << "Created new pulse for block " << outputPulse->GetOutputVolumeID().Top(m_depth) << ".\n"
               << "Resulting pulse is: \n"
               << *outputPulse << Gateendl << Gateendl ;
    outputPulseList->push_back(outputPulse);
  }

  // Free temporary variables used by the centroid policy
  if (m_policy=="TakeEnergyCentroid")
  {
    free(final_time);
    free(final_crystal_posX);
    free(final_crystal_posY);
    free(final_crystal_posZ);
    free(final_nb_pulses);
  }
  // Free other variables
  free(final_energy);
  free(final_pulses);

  if (nVerboseLevel==1)
  {
    G4cout << "[" << GetObjectName() << "::ProcessPulseList]: returning output pulse-list with " << outputPulseList->size() << " entries\n";
    GatePulseIterator iterOut;
    for (iterOut = outputPulseList->begin() ; iterOut != outputPulseList->end() ; ++iterOut)
      G4cout << **iterOut << Gateendl;
    G4cout << Gateendl;
  }

  return outputPulseList;
}

// S. Stute: obsolete method which is not used anymore. The content was moved in upper method ProcessPulseList overloaded right above.
void GateReadoutOld::ProcessOnePulse(const GatePulse* ,GatePulseList& )
{
}




void GateReadoutOld::DescribeMyself(size_t indent)
{
  G4cout << GateTools::Indent(indent) << "ReadoutOld at depth:      " << m_depth << Gateendl;
  G4cout << GateTools::Indent(indent) << "  --> policy: ";
  if (m_policy=="TakeEnergyWinner") G4cout << "TakeEnergyWinner\n";
  else if (m_policy=="TakeEnergyCentroid") G4cout << "TakeEnergyCentroid\n";
  else G4cout << "Unknown policy !\n";
}

