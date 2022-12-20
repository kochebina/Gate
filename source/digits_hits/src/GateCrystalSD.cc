/*----------------------
  Copyright (C): OpenGATE Collaboration

  This software is distributed under the terms
  of the GNU Lesser General  Public Licence (LGPL)
  See LICENSE.md for further details
  ----------------------*/

#include "GateCrystalSD.hh"

#include "../include/GateSinglesDigitizer.hh"
#include "GateDetectorConstruction.hh"
#include "GateHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4PrimaryParticle.hh"

#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4TransportationManager.hh"

#include "GateVSystem.hh"
#include "GateRotationMove.hh"
#include "GateOrbitingMove.hh"
#include "GateEccentRotMove.hh"
#include "GateSystemListManager.hh"
#include "GateVVolume.hh"

//OK GND 2022
#include "GateDigitizerMgr.hh"
//#include "GateDigitizerInitializationModule.hh"

#include "GateObjectStore.hh"
#include "GateEmittedGammaInformation.hh"
#include "GateOutputMgr.hh"
// Name of the hit collection
//const G4String GateCrystalSD::m_CrystalCollectionName = "crystalCollection"; //CrystalSD"; //has to be hardcoded somewhere. Add modifs for multiple HCs


//------------------------------------------------------------------------------
// Constructor
GateCrystalSD::GateCrystalSD(const G4String& name)
:G4VSensitiveDetector(name),
 m_system(0)
{
	// G4cout<<"GateCrystalSD::GateCrystalSD "<< name<<G4endl;
	//collectionName.insert(m_CrystalCollectionName);
	//OK GND 2022
	G4String collName=name+"Collection";
	collectionName.insert(collName);
//	collectionID = -1;

	GateDigitizerMgr* digitizerMgr=GateDigitizerMgr::GetInstance();
	//digitizerMgr->AddNewSD(this);

	GateSinglesDigitizer* digitizer = new GateSinglesDigitizer(digitizerMgr,"Singles",this);
	digitizerMgr->AddNewSinglesDigitizer(digitizer);
	//digitizerMgr->AddNewSinglesDigitizer( new GateSinglesDigitizer(digitizerMgr,"Singles",this));

	digitizerMgr->AddNewSD(this);
	//In order to get unique HCID each time when we create a new SD
	// one can take GetCollectionCapacity of G4SDManager
	// by default there is always phantomSD attached
	// thus: 1st HCID = 1
	// This is done in order to replace a block from Intialize()
	/*
	 static G4int HCID=-1;
	 if(HCID<0 ) // call only in the first event
	 			{
				HCID = G4SDManager::GetSDMpointer()->GetCollectionID(GetName()+"Collection");

			}
	 */

	HCID = G4SDManager::GetSDMpointer()->GetCollectionCapacity() ;

	//OK GND 2022 For GateToTree class adaptaion
  	GateOutputMgr::GetInstance()->RegisterNewHitsCollection(collName,false); // ??? save by default or not

}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Destructor
GateCrystalSD::~GateCrystalSD()
{
	//G4cout<<"GateCrystalSD::~GateCrystalSD()"<<G4endl;
   if(m_systemList) delete m_systemList;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
GateCrystalSD *GateCrystalSD::Clone() const {
    auto clone = new GateCrystalSD(SensitiveDetectorName);
    clone->thePathName = thePathName;
    clone->fullPathName = fullPathName;
    clone->verboseLevel = verboseLevel;
    clone->active = active;
    clone->ROgeometry = ROgeometry;
    clone->filter = filter;

    clone->m_system = m_system;
    clone->m_systemList = m_systemList;
    clone->crystalHitsCollection = crystalHitsCollection;

    return clone;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Method overloading the virtual method Initialize() of G4VSensitiveDetector
void GateCrystalSD::Initialize(G4HCofThisEvent*HCE)
{
	//G4cout<<"GateCrystalSD::Initialize("<<G4endl;
	crystalHitsCollection=new GateHitsCollection(GetName(),collectionName[0]);

	HCE->AddHitsCollection(HCID, crystalHitsCollection);

}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Implementation of the pure virtual method ProcessHits().
// This methods generates a GateHitGND and stores it into the SD's hit collection
//G4bool GateCrystalSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
G4bool GateCrystalSD::ProcessHits(G4Step*aStep, G4TouchableHistory*)
{
	//G4cout<<"GateCrystalSD::ProcessHits"<<G4endl;

  // Get the track information
  G4Track* aTrack       = aStep->GetTrack();
  G4int    trackID      = aTrack->GetTrackID();
  G4int    parentID     = aTrack->GetParentID();
  // Seb Modif 5/4/2016
  G4double trackLength  = aTrack->GetTrackLength();
  G4double trackLocalTime = aTrack->GetLocalTime();

  G4String partName     = aTrack->GetDefinition()->GetParticleName();
  G4int    PDGEncoding  = aTrack->GetDefinition()->GetPDGEncoding();

  //Get information about gamma ( generated by ExtendedVSource )
  G4int source_type = static_cast<G4int>(GateEmittedGammaInformation::SourceKind::NotDefined);
  G4int decay_type = static_cast<G4int>(GateEmittedGammaInformation::DecayModel::None);
  G4int gamma_type = static_cast<G4int>(GateEmittedGammaInformation::GammaKind::Unknown);
  G4PrimaryParticle* primary_particle = aTrack->GetDynamicParticle()->GetPrimaryParticle();
  if( primary_particle != nullptr )
  {
   GateEmittedGammaInformation* info = dynamic_cast<GateEmittedGammaInformation*>( primary_particle->GetUserInformation() );
   if ( info != nullptr )
   {
    source_type = static_cast<G4int>( info->GetSourceKind() );
    decay_type = static_cast<G4int>( info->GetDecayModel() );
    gamma_type = static_cast<G4int>( info->GetGammaKind() );   
   }
  }

  // Get the step-points
  G4StepPoint  *oldStepPoint = aStep->GetPreStepPoint(),
      	       *newStepPoint = aStep->GetPostStepPoint();


  //  Get the process name
  const G4VProcess* process = newStepPoint->GetProcessDefinedStep();
  G4String processName = ( (process != NULL) ? process->GetProcessName() : G4String() ) ;

  //  For all processes except transportation, we select the PostStepPoint volume
  //  For the transportation, we select the PreStepPoint volume
  const G4TouchableHistory* touchable;
  if ( processName == "Transportation" )
      touchable = (const G4TouchableHistory*)(oldStepPoint->GetTouchable() );
  else
      touchable = (const G4TouchableHistory*)(newStepPoint->GetTouchable() );


  GateVolumeID volumeID(touchable);


  if (volumeID.IsInvalid())
    G4Exception( "GateCrystalSD::ProcessHits", "ProcessHits", FatalException, "could not get the volume ID! Aborting!\n");

  // Get the hit global position
  //Modifs Seb 22-06-2011
  //G4ThreeVector position = oldStepPoint->GetPosition();
  G4ThreeVector position = newStepPoint->GetPosition();
  
  // Get the hit momentumDirecton
  G4ThreeVector momentumDirection = newStepPoint->GetMomentumDirection();

  // Compute the hit local position
  // (It will be in the reference frame of the PreStepPoint volume for a transportation hit)
  G4ThreeVector localPosition = volumeID.MoveToBottomVolumeFrame(position);


  // Get the scanner position and rotation angle
/*  GateSystemComponent* baseComponent = GetSystem()->GetBaseComponent();*/
  GateVSystem* system = FindSystem(volumeID);
  GateSystemComponent* baseComponent = system->GetBaseComponent();
  G4ThreeVector scannerPos = baseComponent->GetCurrentTranslation();
  G4double scannerRotAngle = 0;


  if ( baseComponent->FindRotationMove() )
    scannerRotAngle = baseComponent->FindRotationMove()->GetCurrentAngle();
  else if ( baseComponent->FindOrbitingMove() )
    scannerRotAngle = baseComponent->FindOrbitingMove()->GetCurrentAngle();
  else if ( baseComponent->FindEccentRotMove() )
    scannerRotAngle = baseComponent->FindEccentRotMove()->GetCurrentAngle();


  // deposit energy in the current step
  G4double edep = aStep->GetTotalEnergyDeposit();

 //  G4double E;

  // stepLength of the current step
  G4double stepLength = aStep->GetStepLength();
  // time of the current step
  G4double aTime = newStepPoint->GetGlobalTime();
  // Create a new crystal hit
  GateHit* aHit = new GateHit();
  //G4int id = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
 // aHit->SetEventID(id);
  // Store the data already obtained into the hit
  aHit->SetPDGEncoding( PDGEncoding );
  aHit->SetEdep( edep );
  aHit->SetStepLength( stepLength );
  aHit->SetTime( aTime );
  aHit->SetGlobalPos( position );
  aHit->SetLocalPos( localPosition );
  aHit->SetProcess( processName );
  aHit->SetTrackID( trackID );
 // Seb Modif 5/4/2016 
  aHit->SetTrackLength( trackLength );
  aHit->SetTrackLocalTime( trackLocalTime );
  aHit->SetMomentumDir( momentumDirection );
  aHit->SetParentID( parentID );
  aHit->SetVolumeID( volumeID );
  aHit->SetScannerPos( scannerPos );
  aHit->SetScannerRotAngle( scannerRotAngle );
  aHit->SetSystemID(system->GetItsNumber());
  aHit->SetSourceType( source_type );
  aHit->SetDecayType( decay_type );
  aHit->SetGammaType( gamma_type );

  // Ask the system to compute the output volume ID and store it into the hit

//Seb Modif 24/02/2009
/*  GateOutputVolumeID outputVolumeID = GetSystem()->ComputeOutputVolumeID(aHit->GetVolumeID());*/
  GateOutputVolumeID outputVolumeID = system->ComputeOutputVolumeID(aHit->GetVolumeID());
  aHit->SetOutputVolumeID(outputVolumeID);

// G4cout<<"This hit is in "<< this->GetName()<< " added to "<< crystalHitsCollection-> GetName () << " "<<G4SDManager::GetSDMpointer()->GetCollectionID(crystalHitsCollection) <<G4endl;
  // Insert the new hit into the hit collection
	//collectionID =
			crystalHitsCollection->insert( aHit );
  return true;
}
//------------------------------------------------------------------------------

//-------------------------------------------------------------------------------
//OK GND 2022
void GateCrystalSD::EndOfEvent(G4HCofThisEvent* HCE)
{
	//G4cout<<"GateCrystalSD::EndOfEvent "<<this->GetName() <<G4endl;

/*	static G4int HCID = -1;
		  if(HCID<0)
		  {
			  G4int  HCID = G4SDManager::GetSDMpointer()->GetCollectionID(GetName()+"Collection");//collectionName[0]);
			  //G4cout<<GetName()+"Collection "<< HCID<<G4endl;
			  //G4cout<< HCID <<" "<<HCE->GetHC(HCID)->GetName ()<<" "<< HCE ->GetHC(HCID)-> GetSize ()   <<G4endl;
		  }
*/

}

//-------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//! Next method underwent an important modification to be compatible with the multi-system approach
G4int GateCrystalSD::PrepareCreatorAttachment(GateVVolume* aCreator)
{

   GateVSystem* creatorSystem = GateSystemListManager::GetInstance()->FindSystemOfCreator(aCreator->GetCreator());
   if (!creatorSystem) {
      G4cout  << Gateendl << Gateendl << "[GateCrystalSD::PrepareCreatorAttachment]:\n"
            << "Volume '" << aCreator->GetObjectName() << "' does not belong to any system.\n"
            << "Your volume must belong to a system to be used with a crystalSD.\n"
            << "Attachment request ignored --> you won't have any hit output from this volume!!!\n" << Gateendl;
      return -1;
   }

   if (!m_systemList)
   {
      m_systemList = new GateSystemList;
      AddSystem(creatorSystem);
   }
   else
   {
      if((m_systemList->back()) != creatorSystem)
         AddSystem(creatorSystem);
   }

   return 0;
}
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// The next three methods were added for the multi-system approach

// Set the system to which the SD is attached
void GateCrystalSD::AddSystem(GateVSystem* aSystem)
{
   m_systemList->push_back(aSystem);
   GateDigitizerMgr::GetInstance()->AddSystem(aSystem);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
GateVSystem* GateCrystalSD::FindSystem(GateVolumeID volumeID)
{
   // MP Garcia (24/03/2014) Modif to handle imbricated boxes between the SPECThead volume and the world
    //size_t m = volumeID.size();
    //G4String hitSystemName = volumeID.GetVolume(m - 2)->GetName();
    G4String hitSystemName = volumeID.GetVolume(1)->GetName();

   size_t n = hitSystemName.size();
   hitSystemName.erase(n-5,5);

   GateVSystem* system = FindSystem(hitSystemName);

   return system;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
GateVSystem* GateCrystalSD::FindSystem(G4String& systemName)
{
   for(GateSystemIterator itr=m_systemList->begin(); itr!=m_systemList->end(); itr++)
   {
      if(systemName.compare((*itr)->GetOwnName()) == 0)
      {
         return *itr;
      }
   }

   return m_systemList->at(-1);
}
