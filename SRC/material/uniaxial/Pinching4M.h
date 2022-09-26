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

// $vision: 1.0 $
// $Date: 2022-08-25 11:13:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Pinching4M.cpp,v $

// Written: Chuandong Xie (chuandongxie@xauat.edu.cn)
// Created: August 2022
// Update:  August 2022
//
// Description: This file contains the class implementation for 
// Pinching4M material based on the Pinching4 material 
// (vision 1.3, Created by NM (nmitra@u.washington.edu)). The original material
// is defined by 4 points on the positive and negative envelopes and a bunch of 
// damage parameters and accounts for 3 types of damage rules: Strength degradation,
// Stiffness degradation and unloading stiffness degradation. This version is trying
// to redefine the points at backbone curves, as well as unloading and reloading points.

#ifndef Pinching4M_h
#define Pinching4M_h

#include <UniaxialMaterial.h>
#include <FileStream.h>
#include <OPS_Stream.h>
#include <Vector.h>

class Pinching4M : public UniaxialMaterial
{
public:
	Pinching4M(int tag,
		double stress1p, double strain1p, double stress2p, double strain2p,
		double stress3p, double strain3p, double stress4p, double strain4p,
		double stress1n, double strain1n, double stress2n, double strain2n,
		double stress3n, double strain3n, double stress4n, double strain4n,
		double rDispP, double rForceP, double uForceP,
		double rDispN, double rForceN, double uForceN,
		double gammaK1, double gammaK2, double gammaK3,
		double gammaK4, double gammaKLimit,
		double gammaD1, double gammaD2, double gammaD3,
		double gammaD4, double gammaDLimit,
		double gammaF1, double gammaF2, double gammaF3,
		double gammaF4, double gammaFLimit, double gammaE, int DmgCyc);

	Pinching4M(int tag,
		double stress1p, double strain1p, double stress2p, double strain2p,
		double stress3p, double strain3p, double stress4p, double strain4p,
		double rDispP, double rForceP, double uForceP,
		double gammaK1, double gammaK2, double gammaK3,
		double gammaK4, double gammaKLimit,
		double gammaD1, double gammaD2, double gammaD3,
		double gammaD4, double gammaDLimit,
		double gammaF1, double gammaF2, double gammaF3,
		double gammaF4, double gammaFLimit, double gammaE, int DmgCyc);

	Pinching4M();
	~Pinching4M();

	int setTrialStrain(double strain, double strainRate = 0.0);
	double getStrain(void);
	double getStress(void);
	double getTangent(void);
	double getInitialTangent(void);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	UniaxialMaterial* getCopy(void);

	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel,
		FEM_ObjectBroker& theBroker);

	void Print(OPS_Stream& s, int flag = 0);

	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int parameterID, Information& info);

protected:

private:
	// Backbone parameters
	double stress1p; double strain1p; double stress2p; double strain2p;
	double stress3p; double strain3p; double stress4p; double strain4p;
	double stress1n; double strain1n; double stress2n; double strain2n;
	double stress3n; double strain3n; double stress4n; double strain4n;
	double stressyp; double stressyn; // Yielding stress
	Vector envlpPosStress; Vector envlpPosStrain;
	Vector envlpNegStress; Vector envlpNegStrain;

	int tagMat;  // material tag

	// Damage parameters
	double gammaK1; double gammaK2; double gammaK3; double gammaK4; double gammaKLimit;
	double gammaD1; double gammaD2; double gammaD3; double gammaD4; double gammaDLimit;
	double gammaF1; double gammaF2; double gammaF3; double gammaF4; double gammaFLimit;
	double gammaE;
	double TnCycle, CnCycle;	// number of cycles contributing to damage calculation
	int DmgCyc;					// flag for indicating whether no. of cycles are to be used for damage calculation

	// Unloading-reloading parameters
	double rDispP; double rForceP; double uForceP;
	double rDispN; double rForceN; double uForceN;

	Vector state3Stress; Vector state3Strain; Vector state4Stress; Vector state4Strain;
	// 0: Low state;
	// 1: Reloading point;
	// 2:
	// 3: high state;
	Vector envlpPosDamgdStress; Vector envlpNegDamgdStress;

	// Trial state variables
	double Tstress;
	double Tstrain;
	double Ttangent;

	// Converged material history parameters
	int Cstate;
	double Cstrain;
	double Cstress;
	double CstrainRate;
	double lowCstateStrain;
	double lowCstateStress;
	double hghCstateStrain;
	double hghCstateStress;
	double CminStrainDmnd; // Minimum historic deformation demands
	double CmaxStrainDmnd; // Maximum historic deformation demands
	double Cenergy;			//
	double CgammaK;
	double CgammaD;
	double CgammaF;
	double gammaKUsed;
	double gammaFUsed;

	// Trial material history parameters
	int Tstate;
	double dstrain;
	double TstrainRate;
	double lowTstateStrain;
	double lowTstateStress;
	double hghTstateStrain; // Largest strain during the loading history
	double hghTstateStress; // Largest stress during the loading history
	double TminStrainDmnd;
	double TmaxStrainDmnd;
	double Tenergy;			//
	double TgammaK;
	double TgammaD;
	double TgammaF;

	// Strength and stiffness parameters;
	double kElasticPos;
	double kElasticNeg;
	double kElasticPosDamgd;
	double kElasticNegDamgd;
	double uMaxDamgd;
	double uMinDamgd;
	double kPostYieldPosDamgd; // Post-yielding stiffness of positive branch
	double kPostYieldNegDamgd; // Post-yielding stiffness of negative branch

	// Energy parameters
	double energyCapacity;
	double kunload;
	double elasticStrainEnergy;

	void SetEnvelope(void);							// Set the initial backbone envelope for the material based upon the input by the user
	void getstate(double, double);					// Determine the state of the material based upon the material history and current stress demand
	double posEnvlpStress(double);					// Return positive damage stress of the material
	double posEnvlpTangent(double);					// Return positive tangent of the material
	double negEnvlpStress(double);					// Return negative damaged stress of the material
	double negEnvlpTangent(double);					// Return negative tangent of the material
	void getState3(Vector&, Vector&, double);		// Form the backbone envelope of state 3
	void getState4(Vector&, Vector&, double);		// Form the backbone envelope of state 4
	double Envlp3Tangent(Vector, Vector, double);	// Determine the tangent of the envelope at state 3
	double Envlp3Stress(Vector, Vector, double);	// Determine the stress of the envelope at state 3
	double Envlp4Tangent(Vector, Vector, double);	// Determine the tangent of the envelope at state 4
	double Envlp4Stress(Vector, Vector, double);	// Determine the stress of the envelope at state 4
	void updateDmg(double, double);					// Determine the damages at a particular state of the material
	
#ifdef _G3DEBUG
FileStream* fg;
#endif

};
#endif
