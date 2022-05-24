/*----------------------
   Copyright (C): OpenGATE Collaboration

This software is distributed under the terms
of the GNU Lesser General  Public Licence (LGPL)
See LICENSE.md for further details
----------------------*/


#ifndef GateReadout_h
#define GateReadout_h 1

#include "globals.hh"
#include <iostream>
#include <vector>
#include "G4ThreeVector.hh"

#include "GateVPulseProcessor.hh"
#include "GateVSystem.hh"
#include "GateArrayComponent.hh"
#include "GateObjectStore.hh"

//#define READOUT_POLICY_WINNER 0
//#define READOUT_POLICY_CENTROID 1

class GateReadoutMessenger;
class GateOutputVolumeID;

/*! \class  GateReadout
    \brief  Pulse-processor modelling a simple PMT readout (maximum energy wins) of a crystal-block

    - GateReadout - by Daniel.Strul@iphe.unil.ch

    - The readout is parameterised by its 'depth': pulses will be summed up if their volume IDs
      are identical up to this depth. For instance, the default depth is 1: this means that
      pulses will be considered as taking place in a same block if the first two figures
      of their volume IDs are identical

      \sa GateVPulseProcessor
*/
class GateReadout : public GateVPulseProcessor
{
  public:

    //! Constructs a new readout attached to a GateDigitizerOld
    GateReadout(GatePulseProcessorChain* itsChain,const G4String& itsName) ;

    //! Destructor
    virtual ~GateReadout() ;

    //! Implementation of the pure virtual method declared by the base class GateClockDependent
    //! print-out the attributes specific of the readout
    virtual void DescribeMyself(size_t indent);

    //! Returns the depth of the readout
    inline G4int GetDepth() const  	      	{ return m_depth; }

    //! Set the depth of the readout
    inline void  SetDepth(G4int aDepth)         { m_depth = aDepth; }

    //! Set the policy of the readout
    inline void SetPolicy(const G4String& aPolicy)  { m_policy = aPolicy; };
    inline G4String GetPolicy() const  	      	{ return m_policy; }

    //! Set the volume for the readout
    inline void SetVolumeName(const G4String& aName) { m_volumeName = aName; };
    inline G4String GetVolumeName() const  	      	{ return m_volumeName; }

    //! Set the volume for the readout even for centroid policy
     inline void ForceDepthCentroid(const G4bool& value) { m_IsForcedDepthCentroid = value; };
     inline G4bool IsDepthForcedCentroid() const  	      	{ return m_IsForcedDepthCentroid; }

    //! Set how the resulting positions should be defined
    inline void SetResultingXY(const G4String& aString) { m_resultingXY= aString;};
    inline G4String GetResultingXY() const  	      	{ return m_resultingXY; };

    inline void SetResultingZ(const G4String& aString){m_resultingZ= aString;};
    inline G4String GetResultingZ() const  	      	{ return m_resultingZ; };


    void SetReadoutParameters();


  protected:
    //! Implementation of the pure virtual method declared by the base class GateVPulseProcessor
    //! This methods processes one input-pulse
    //! It is is called by ProcessPulseList() for each of the input pulses
    //! The result of the pulse-processing is incorporated into the output pulse-list
    void ProcessOnePulse(const GatePulse* inputPulse,GatePulseList& outputPulseList);
    //! Overload the virtual (not pure) method of GateVPulseProcessor
    GatePulseList* ProcessPulseList(const GatePulseList* inputPulseList);



  private:
    //! The default is the one parameter that defines how a readout works:
    //! pulses will be summed up if their volume IDs are identical up to this depth.
    //! For instance, the default depth is 1: this means that pulses will be considered as
    //! taking place in a same block if the first two figures of their volume IDs are identical
    G4int m_depth;

    //! S. Stute: add an option to choose the policy of the readout (using two define integers; see the beginning of this file)
    //G4int m_policy;
    G4String m_policy;
    GateVSystem* m_system;
    G4int m_nbCrystalsX;
    G4int m_nbCrystalsY;
    G4int m_nbCrystalsZ;
    G4int m_nbCrystalsXY;
    G4int m_systemDepth;
    G4int m_crystalDepth;
    GateArrayComponent* m_crystalComponent;

    G4String m_volumeName;
    G4bool m_IsForcedDepthCentroid;

    G4String m_resultingXY;
    G4String m_resultingZ;
    G4bool   m_IsFirstEntrance;//Entrance

    std::vector<int> numberOfComponentForLevel; //!< Table of number of element for each geometric level
    G4int numberOfHigherLevels ;  //!< number of geometric level higher than the one chosen by the user
    G4int numberOfLowerLevels ;  //!< number of geometric level higher than the one chosen by the user
    GateReadoutMessenger *m_messenger;	  //!< Messenger for this readout
};


#endif
