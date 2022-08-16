/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Version: 1.0 $
// $Date: 2022-03-25 14:32:14
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SPSW03.h,v $

// Written: Chuandong Xie
//			Xi'an University of Architecture and Techonology 
// Create: 2022-03-25 14:32:14
// Mail: chuandongxie@xauat.edu.cn
// 
// Description: This file contains the class definition for SPSW03.
//				SPSW03 is developed based on SPSW02 but considers 
//				imperfection of infill plates
//-----------------------------------------------------------------------
//                  SPSW03 by Chuandong Xie (2022)
//-----------------------------------------------------------------------


#ifndef SPSW03_h
#define SPSW03_h
#include <UniaxialMaterial.h>
class SPSW03:public UniaxialMaterial
{
public:
	SPSW03(int tag, double fpy, double nFac, double E0, double b, double t, double hs,
		double l, double R, double epsPCFac, double pstCapEFac, double gama, double c, double resFac,
		int mode, double impf);

	SPSW03(int tag, double E0, double b, double _FTS, double _FCS, double _cmpUnldngEFac,
		double _sigTEfac, double _sigTFfac, double _epsTFfac, double R, double epsPCFac,
		double pstCapEFac, double gama, double c, double resFac);
	
	SPSW03(void);
	virtual ~SPSW03();

	const char* getClassType(void) const {return "SPSW03";};

	double getInitialTangent();
	UniaxialMaterial* getCopy();

	virtual int setTrialStrain(double strain, double strainRate = 0.0);
	virtual double getStrain();
	virtual double getStress();
	virtual double getTangent();
	virtual int commitState();
	virtual int revertToLastCommit();
	virtual int revertToStart();
	virtual int sendSelf(int commitTag, Channel& theChannel);
	virtual int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker & theBroker);
	virtual void Print(OPS_Stream & s, int flag = 0);

    // virtual Response *setResponse (const char **argv, int argc, 
	// OPS_Stream &theOutputStream);
    // virtual int getResponse (int responseID, Information &matInformation);
	virtual double getEnergy();
	virtual double getInitYieldStrain() {return FTS/E0;}

protected:
	
private:
	// FUNCTIONs --------------------------------------------------------------------------
	void updateDamage();
	void Calc_sigcr ();
	void MenegottoPinto (double epsc, double bT, double R, double& sigc, double& ec);
	//void MenegottoPinto (double epsc, double bT, double R, double& sigc, double& ec, double er, double sr, double e0, double s0);
	//matpar : PLATE PROPERTIES

	//  MATERIAL INPUTS -------------------------------------------------------------------
	double E0;			//initial stiffness
	double b;			//hardening ratio

	// '-geom'
	double fpy;			//yield stress
	double nFac;		// negative strength
	double t;			//plate thickness
	double hs;			//plate height
	double l;			//plate length

	// '-params'
	double FTS, FCS;		//yield stresses before damage
	double Fts	;			//trial yield stress of a tension strip (damage not commited)
	double Fcs	;			//yield stress of a compression strip
	double cmpUnldngEFac;
	double sigTEfac, sigTFfac, epsTFfac;

	// '-R'
	double R;

	// '-damage'
	double epsPCFac;		//ratio between post cap strain and yield strain
	double pstcpEFac;		//ratio between post cap stifness and initial stiffness
	double resFac;
	double gama, FailEnerg;			//damage parameters
	double c;				// exponent defining the rate of deterioration. Resonable range: 1.0~2.0.
							// A value of 2.0 will slow sown early deteriorarion and accelerate deterioration 
							// in later cycles, whereas a value of 1.0 implies an almost constant rate of deterioration.

	// '-IPF'
	int mode;				// mode 1 for half sine or mode 2 for full sine
	double impf;			// imperfection ratio, e.g. 1000 for 1/1000, default: 1000

	// INTERNAL VARIABLES -----------------------------------------------------------------
	double epsmax;		//
	double sigmax;		//
	double epss0;		//
	double sigs0;		//
	double epsr;		//
	double sigr;		//
	double epsTF;		//
	double plstr;		//
	double sig;		//
	double e;		//
	double eps;		//
	//bool capStrFlg;		//flag to indicate whether capping strain has been passed or not
	bool givenParams;	// params or geom
	double lStrip;		// strip length
	double outDeform;	// out-of-plane deformation
	double lImpf;		// strip length considering imperfection
	double epsImpf;		// extra strain brought by imperfection
	double E0Impf;		// actual stiffness after considering imperfection


	// STORED VALUES ----------------------------------------------------------------------
	int kon;		// loading state
					/*	
						kon == 0 (value at start): elastic with (fabs(eps) < epsyn || eps < -epsyn)
						kon == 11: comp. branch
						kon == 12: reversal from comp. when sarts at eps >  plstr - 2.*Fcs/(cmpUnldngEFac*E0) (pinching branch excluded)
						kon == 13: reversal from comp. when sarts at eps <  plstr - 2.*Fcs/(cmpUnldngEFac*E0)  (pinching branch included)
						kon == 21: Menegotto-Pinto tension loading
					*/

	// VARIABLES FOR LOAD HISTORY MEMORY --------------------------------------------------
	double epsmaxP;		//max eps in tension
	double sigmaxP;		//max sig in tension
	double epss0P;		//eps at asymptotes intersection
	double sigs0P;		//sig at ...
	double epssrP;		//eps at last inversion point
	double sigsrP;		//sig at ...
	double epsTFP;		//eps at tension field redevelopment point
	double plstrP;		//plastic strain
	int konP;			//index for loading-unloading
	double epsP;			//eps at previous converged step
	double sigP;			//stress at ...
	double eP;			//stiffness modulus at ...
	double excurEnergP, totalEnergP, betaP;
	double E0ImpfP;

	// TEMPORARY VARIABLES STORED HERE FOR EFFICIENCY -------------------------------------
	double excurEnerg, totalEnerg, beta;
};
#endif
