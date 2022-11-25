
/*----------------------
  Copyright (C): OpenGATE Collaboration

  This software is distributed under the terms
  of the GNU Lesser General  Public Licence (LGPL)
  See LICENSE.md for further details
  ----------------------*/

/*!
  \class  GateEnergyResolution

  This class is not used by GATE !
  The purpose of this class is to help to create new users digitizer module(DM).

   - Create your DM by coping this class and GateDummyDigitizerMessenger class for your DM messenger
   - Places to change are marked with // ****** comment and called with "dummy" names
   - Include your module to GateDigitizerMessenger in the method DoInsertion(..)

	If you adapting some already exiting class from Old Gate Digitizer here is some of the tips
	- Digitize () is a fusion of GateVPulseProcessor::ProcessPulseList and GateXXX::ProcessOnePulse
	- pulse --> Digi
	- outputPulseList --> OutputDigiCollectionVector
	- inputPulse-->inputDigi
	- outputPulse --> m_outputDigi
	- how to adapt iterators check GateAdder class


  !!!! DO NOT FORGET TO WRITE A SHORT EXPLANATION ON WHAT DOES YOUR DM !!!!
	The example is also given in .hh
*/

#include "GateEnergyResolution.hh"
#include "GateEnergyResolutionMessenger.hh"
#include "GateDigi.hh"

#include "GateDigitizerMng.hh"
#include "GateConstants.hh"



#include "G4SystemOfUnits.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"



GateEnergyResolution::GateEnergyResolution(GateDigitizer *digitizer)
  :GateVDigitizerModule("energyResolution","digitizerMng/"+digitizer->GetSD()->GetName()+"/SinglesDigitizer/"+digitizer->m_digitizerName+"/energyResolution",digitizer,digitizer->GetSD()),
   m_reso(0),
   m_resoMin(0),
   m_resoMax(0),
   m_eref(0),
   m_slope(0),
   m_outputDigi(0),
   m_OutputDigiCollection(0),
   m_digitizer(digitizer)
 {
	G4String colName = digitizer->GetOutputName() ;
	collectionName.push_back(colName);

	new G4UnitDefinition ( "1/electronvolt", "1/eV", "Energy Slope", 1/electronvolt );
	new G4UnitDefinition ( "1/kiloelectronvolt", "1/keV", "Energy Slope", 1/kiloelectronvolt );
	new G4UnitDefinition ( "1/megaelectronvolt", "1/MeV", "Energy Slope", 1/megaelectronvolt );
	new G4UnitDefinition ( "1/gigaelectronvolt", "1/GeV", "Energy Slope", 1/gigaelectronvolt );
	new G4UnitDefinition ( "1/joule", "1/J", "Energy Slope", 1/joule );

	m_Messenger = new GateEnergyResolutionMessenger(this);
 }


GateEnergyResolution::~GateEnergyResolution()
{
  delete m_Messenger;

}




void GateEnergyResolution::Digitize()
{
	if( m_resoMin!=0 && m_resoMax!=0 && m_reso!=0)
	{
		GateError("***ERROR*** Resolution is ambiguous: you can set /resolution OR range for resolutions with /resolutionMin and /resolutionMax!");
	}

	G4String digitizerName = m_digitizer->m_digitizerName;
	G4String outputCollName = m_digitizer-> GetOutputName();

	m_OutputDigiCollection = new GateDigiCollection(GetName(),outputCollName); // to create the Digi Collection

	G4DigiManager* DigiMan = G4DigiManager::GetDMpointer();



	GateDigiCollection* IDC = 0;
	IDC = (GateDigiCollection*) (DigiMan->GetDigiCollection(m_DCID));

	GateDigi* inputDigi = new GateDigi();

  if (IDC)
     {
	  G4int n_digi = IDC->entries();

	  //loop over input digits
	  for (G4int i=0;i<n_digi;i++)
	  {
		  inputDigi=(*IDC)[i];

		  if( m_resoMin!=0 && m_resoMax!=0)
			  m_reso = G4RandFlat::shoot(m_resoMin, m_resoMax);

		  G4double energy= inputDigi->GetEnergy();
		  G4double sigma;
		  G4double resolution;

		  if (m_slope == 0 )
			  //Apply InverseSquareBlurringLaw
		  {
			  G4cout<<"InverseSquareBlurringLaw"<<G4endl;
			  resolution = m_reso * sqrt(m_eref)/ sqrt(energy);

		  }
		  else
			  //Apply LinearBlurringLaw
			  resolution = m_slope * (energy - m_eref) + m_reso;

	      sigma =(resolution*energy)/GateConstants::fwhm_to_sigma;


	      //return ((m_resolution * sqrt(m_eref)) / sqrt(energy));

		  G4double outEnergy=G4RandGauss::shoot(energy,sigma);

		  m_outputDigi = new GateDigi(*inputDigi);
		  m_outputDigi->SetEnergy(outEnergy);

		  if (nVerboseLevel>1)
		 	  G4cout << "[GateEnergyResolution::Digitize]: Created new digi from one with energy " << inputDigi->GetEnergy() << ".\n"
		 		 << "Resulting digi has energy: "<< m_outputDigi->GetEnergy() << Gateendl << Gateendl ;

		  m_OutputDigiCollection->insert(m_outputDigi);

	  }
    }
  else
    {
  	  if (nVerboseLevel>1)
  	  	G4cout << "[GateEnergyResolution::Digitize]: input digi collection is null -> nothing to do\n\n";
  	    return;
    }
  StoreDigiCollection(m_OutputDigiCollection);


}



void GateEnergyResolution::DescribeMyself(size_t indent )
{
	  G4cout << GateTools::Indent(indent) << "Resolution of " << m_reso  << " for " <<  G4BestUnit(m_eref,"Energy") << Gateendl;
;
}
