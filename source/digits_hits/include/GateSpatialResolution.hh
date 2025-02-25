/*----------------------
   Copyright (C): OpenGATE Collaboration

This software is distributed under the terms
of the GNU Lesser General  Public Licence (LGPL)
See LICENSE.md for further details
----------------------*/

// OK GND 2022

/*! \class  GateSpatialResolution
    \brief  GateSpatialResolution does some dummy things with input digi
    to create output digi

    \sa GateSpatialResolution, GateSpatialResolutionMessenger
*/

#ifndef GateSpatialResolution_h
#define GateSpatialResolution_h 1

#include "GateVDigitizerModule.hh"
#include "GateDigi.hh"
#include "GateClockDependent.hh"
#include "GateCrystalSD.hh"

#include "globals.hh"

#include "GateSpatialResolutionMessenger.hh"

#include "G4VoxelLimits.hh"
#include "G4TouchableHistoryHandle.hh"
#include "GateSinglesDigitizer.hh"

class GateVDistribution;

class GateSpatialResolution : public GateVDigitizerModule

{
public:
  
  GateSpatialResolution(GateSinglesDigitizer *digitizer, G4String name);
  ~GateSpatialResolution();
  
  void Digitize() override;



  //! These functions return the resolution in use.
    G4double GetFWHM()				{ return m_fwhm; }
    GateVDistribution* GetFWHMxdistrib()	{ return m_fwhmXdistrib; }
    GateVDistribution* GetFWHMydistrib()    	{ return m_fwhmYdistrib; }
    GateVDistribution* GetFWHMzdistrib()    	{ return m_fwhmZdistrib; }
    G4String GetNameAxis()				   { return m_nameAxis;     }
    GateVDistribution* GetFWHMDistrib2D()	{ return m_fwhmDistrib2D; }

    G4double GetFWHMx()         { return m_fwhmX; }
    G4double GetFWHMy()			{ return m_fwhmY; }
    G4double GetFWHMz()			{ return m_fwhmZ; }

    //! These functions set the spresolution of a gaussian spblurring.
    /*!
      If you want a resolution of 10%, SetSpresolution(0.1)
    */
    void SetFWHM(G4double val)   { m_fwhm = val;  }
    void SetFWHMxdistrib(GateVDistribution* dist)  { m_fwhmXdistrib= dist; }
    void SetFWHMydistrib(GateVDistribution* dist)  { m_fwhmYdistrib = dist; }
    void SetFWHMzdistrib(GateVDistribution* dist)  { m_fwhmZdistrib = dist; }


    void SetNameAxis(const G4String& name) {m_nameAxis=name;}
    void SetFWHMDistrib2D(GateVDistribution* dist)  { m_fwhmDistrib2D= dist;}

    void SetFWHMx(G4double val)   { m_fwhmX = val;  }
    void SetFWHMy(G4double val)   { m_fwhmY = val;  }
    void SetFWHMz(G4double val)   { m_fwhmZ = val;  }
    void SetSpatialResolutionParameters();
    inline void ConfineInsideOfSmallestElement(const G4bool& value) { m_IsConfined = value; };
    inline G4bool IsConfinedInsideOfSmallestElement() const  	      	{ return m_IsConfined; }
    inline void SetUseTruncatedGaussian(const G4bool& value) { m_UseTruncatedGaussian = value; };
    inline G4bool GetUseTruncatedGaussian(const G4bool& value) 	      	{ return m_UseTruncatedGaussian; }

    void UpdatePos(G4double ,G4double ,G4double );
    void LocateOutputDigi(GateDigi* inputDigi, G4double PxNew,G4double PyNew,G4double PzNew);

    void UpdateVolumeID();


    //! Implementation of the pure virtual method declared by the base class GateClockDependent
    //! print-out the attributes specific of the blurring
    void DescribeMyself(size_t );

protected:
    G4double m_fwhm;


    G4double m_fwhmX;
    G4double m_fwhmY;
    G4double m_fwhmZ;

    GateVDistribution* m_fwhmXdistrib;
    GateVDistribution* m_fwhmYdistrib;
    GateVDistribution* m_fwhmZdistrib;

    GateVDistribution* m_fwhmDistrib2D;

    G4String m_nameAxis;


    G4bool m_IsConfined;
    G4bool m_UseTruncatedGaussian;
    G4Navigator* m_Navigator;
    G4TouchableHistoryHandle m_Touchable;



private:

   G4int m_systemDepth;


  GateDigi* m_outputDigi;;
  GateSpatialResolutionMessenger *m_Messenger;

  GateDigiCollection*  m_OutputDigiCollection;

  GateSinglesDigitizer *m_digitizer;
  G4bool   m_IsFirstEntrance;
  G4VoxelLimits limits;
  G4double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;
  G4AffineTransform at;
};

#endif








