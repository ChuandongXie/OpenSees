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

#include <Pinching4M.h>
#include <OPS_Globals.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <OPS_Stream.h>
#include <string.h>
#include <elementAPI.h>
#include <Information.h>
#include <Parameter.h>

void* OPS_Pinching4M(void)
{
	int numdata = OPS_GetNumRemainingInputArgs();
	if (numdata != 40 && numdata != 29) {
		opserr << "WARNING: Insufficient arguments\n";
		return 0;
	}

	int tag;
	numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
		return 0;
	}

	UniaxialMaterial* theMaterial = 0;
	int tDmg = -1;
	if (OPS_GetNumRemainingInputArgs() == 39) {
		double data[38];
		numdata = 38;
		if (OPS_GetDoubleInput(&numdata, data)) {
			return 0;
		}
		const char* type = OPS_GetString();
		if (strcmp(type, "cycle") == 0 || strcmp(type, "Cycle") == 0 || strcmp(type, "damageCycle") == 0 || strcmp(type, "DamageCycle") == 0) {
			tDmg = 1;
		}
		else if (strcmp(type, "energy") == 0 || strcmp(type, "Energy") == 0 || strcmp(type, "damageEnergy") == 0 || strcmp(type, "DamageEnergy") == 0) {
			tDmg = 0;
		}
		else {
			opserr << "WARNING invalid type of damage calculation specified\n";
			opserr << "Pinching4M material: " << tag << endln;
			return 0;
		}

		theMaterial = new Pinching4M(tag,
			data[0], data[1], data[2], data[3],
			data[4], data[5], data[6], data[7],
			data[8], data[9], data[10], data[11],
			data[12], data[13], data[14], data[15],
			data[16], data[17], data[18],
			data[19], data[20], data[21],
			data[22], data[23], data[24],
			data[25], data[26],
			data[27], data[28], data[29],
			data[30], data[31],
			data[32], data[33], data[34],
			data[35], data[36], data[37], tDmg);
	}
	else if (OPS_GetNumRemainingInputArgs() == 28) {
		double data[27];
		numdata = 27;
		if (OPS_GetDoubleInput(&numdata, data)) {
			return 0;
		}
		const char* type = OPS_GetString();
		if (strcmp(type, "cycle") == 0 || strcmp(type, "Cycle") == 0 || strcmp(type, "damageCycle") == 0 || strcmp(type, "DamageCycle") == 0) {
			tDmg = 1;
		}
		else if (strcmp(type, "energy") == 0 || strcmp(type, "Energy") == 0 || strcmp(type, "damageEnergy") == 0 || strcmp(type, "DamageEnergy") == 0) {
			tDmg = 0;
		}
		else {
			opserr << "WARNING invalid type of damage calculation specified\n";
			opserr << "Pinching4M material: " << tag << endln;
			return 0;
		}

		theMaterial = new Pinching4M(tag,
			data[0], data[1], data[2], data[3],
			data[4], data[5], data[6], data[7],
			data[8], data[9], data[10],
			data[11], data[12], data[13],
			data[14], data[15],
			data[16], data[17], data[18],
			data[19], data[20],
			data[21], data[22], data[23],
			data[24], data[25], data[26], tDmg);
	}

	if (theMaterial == 0) {
		opserr << "WARNING: failed to create Pinching4M material\n";
		return 0;
	}
	return theMaterial;
}

Pinching4M::Pinching4M(int tag,
	double f1p, double d1p, double f2p, double d2p,
	double f3p, double d3p, double f4p, double d4p,
	double f1n, double d1n, double f2n, double d2n,
	double f3n, double d3n, double f4n, double d4n,
	double mdp, double mfp, double msp,
	double mdn, double mfn, double msn,
	double gk1, double gk2, double gk3, double gk4, double gklim,
	double gd1, double gd2, double gd3, double gd4, double gdlim,
	double gf1, double gf2, double gf3, double gf4, double gflim, double ge, int dc) :
	UniaxialMaterial(tag, MAT_TAG_Pinching4M),
	stress1p(f1p), strain1p(d1p), stress2p(f2p), strain2p(d2p),
	stress3p(f3p), strain3p(d3p), stress4p(f4p), strain4p(d4p),
	stress1n(f1n), strain1n(d1n), stress2n(f2n), strain2n(d2n),
	stress3n(f3n), strain3n(d3n), stress4n(f4n), strain4n(d4n),
	envlpPosStress(6), envlpPosStrain(6), envlpNegStress(6), envlpNegStrain(6), tagMat(tag),
	gammaK1(gk1), gammaK2(gk2), gammaK3(gk3), gammaK4(gk4), gammaKLimit(gklim),
	gammaD1(gd1), gammaD2(gd2), gammaD3(gd3), gammaD4(gd4), gammaDLimit(gdlim),
	gammaF1(gf1), gammaF2(gf2), gammaF3(gf3), gammaF4(gf4), gammaFLimit(gflim),
	gammaE(ge), TnCycle(0.0), CnCycle(0.0), DmgCyc(dc),
	rDispP(mdp), rForceP(mfp), uForceP(msp), rDispN(mdn), rForceN(mfn), uForceN(msn),
	state3Stress(4), state3Strain(4), state4Stress(4), state4Strain(4),
	envlpPosDamgdStress(6), envlpNegDamgdStress(6)
{
	bool error = false;

	// Positive backbone parameters
	if (strain1p <= 0.0)
		error = true;
	if (strain2p <= 0.0)
		error = true;
	if (strain3p <= 0.0)
		error = true;
	if (strain4p <= 0.0)
		error = true;

	// Negative backbone parameters
	if (strain1n >= 0.0)
		error = true;
	if (strain2n >= 0.0)
		error = true;
	if (strain3n >= 0.0)
		error = true;
	if (strain4n >= 0.0)
		error = true;

	if (error) {
		opserr << "ERROR: -- input backbone is not unique(one-to-one), Pinching4M::Pinching4M" << "\a";
	}

	// Yielding stress
	stressyp = strain1p;
	stressyn = strain1n;

	// Initialize envelope variables
	envlpPosStress.Zero();
	envlpPosStrain.Zero();
	envlpNegStress.Zero();
	envlpNegStrain.Zero();

	// Initialize state3/state4 variables
	state3Stress.Zero();
	state3Strain.Zero();
	state4Stress.Zero();
	state4Strain.Zero();

	energyCapacity = 0.0;
	kunload = 0.0;
	elasticStrainEnergy = 0.0;

	// Post yielding stiffness
	kPostYieldPosDamgd = 0.0;
	kPostYieldNegDamgd = 0.0;

	// Set envelope slopes
	this->SetEnvelope();

	envlpPosDamgdStress = envlpPosStress;
	envlpNegDamgdStress = envlpNegStress;

	// Initialize history variables
	this->revertToStart();
	this->revertToLastCommit();
}

Pinching4M::Pinching4M(int tag,
	double f1p, double d1p, double f2p, double d2p,
	double f3p, double d3p, double f4p, double d4p,
	double mdp, double mfp, double msp,
	double gk1, double gk2, double gk3, double gk4, double gklim,
	double gd1, double gd2, double gd3, double gd4, double gdlim,
	double gf1, double gf2, double gf3, double gf4, double gflim, double ge, int dc) :
	UniaxialMaterial(tag, MAT_TAG_Pinching4M),
	stress1p(f1p), strain1p(d1p), stress2p(f2p), strain2p(d2p),
	stress3p(f3p), strain3p(d3p), stress4p(f4p), strain4p(d4p),
	envlpPosStress(6), envlpPosStrain(6), envlpNegStress(6), envlpNegStrain(6), tagMat(tag),
	gammaK1(gk1), gammaK2(gk2), gammaK3(gk3), gammaK4(gk4), gammaKLimit(gklim),
	gammaD1(gd1), gammaD2(gd2), gammaD3(gd3), gammaD4(gd4), gammaDLimit(gdlim),
	gammaF1(gf1), gammaF2(gf2), gammaF3(gf3), gammaF4(gf4), gammaFLimit(gflim),
	gammaE(ge), TnCycle(0.0), CnCycle(0.0), DmgCyc(dc),
	rDispP(mdp), rForceP(mfp), uForceP(msp),
	state3Stress(4), state3Strain(4), state4Stress(4), state4Strain(4),
	envlpPosDamgdStress(6), envlpNegDamgdStress(6)
{
	bool error = false;

	// Positive backbone parameters
	if (strain1p <= 0.0)
		error = true;
	if (strain2p <= 0.0)
		error = true;
	if (strain3p <= 0.0)
		error = true;
	if (strain4p <= 0.0)
		error = true;

	if (error) {
		opserr << "ERROR: -- input backbone is not unique(one-to-one), Pinching4M::Pinching4M" << "\a";
	}

	// Set negative variables
	strain1n = -strain1p;
	stress1n = -stress1p;
	strain2n = -strain2p;
	stress2n = -stress2p;
	strain3n = -strain3p;
	stress3n = -stress3p;
	strain4n = -strain4p;
	stress4n = -stress4p;
	rDispN = rDispP;
	rForceN = rForceP;
	uForceN = uForceP;

	// Yielding stress
	stressyp = strain1p;
	stressyn = strain1n;

	// Initialize envelope variables
	envlpPosStress.Zero();
	envlpPosStrain.Zero();
	envlpNegStress.Zero();
	envlpNegStrain.Zero();

	// Initialize state3/state4 variables
	state3Stress.Zero(); // Four values inside
	state3Strain.Zero();
	state4Stress.Zero();
	state4Strain.Zero();

	energyCapacity = 0.0;
	kunload = 0.0;
	elasticStrainEnergy = 0.0;

	// Post yielding stiffness
	kPostYieldPosDamgd = 0.0;
	kPostYieldNegDamgd = 0.0;

	// Set envelope slopes
	this->SetEnvelope();

	envlpPosDamgdStress = envlpPosStress;
	envlpNegDamgdStress = envlpNegStress;

	// Initialize history variables
	this->revertToStart();
	this->revertToLastCommit();
}


Pinching4M::Pinching4M() :
	UniaxialMaterial(0, MAT_TAG_Pinching4M),
	stress1p(0.0), strain1p(0.0), stress2p(0.0), strain2p(0.0),
	stress3p(0.0), strain3p(0.0), stress4p(0.0), strain4p(0.0),
	stress1n(0.0), strain1n(0.0), stress2n(0.0), strain2n(0.0),
	stress3n(0.0), strain3n(0.0), stress4n(0.0), strain4n(0.0),
	gammaK1(0.0), gammaK2(0.0), gammaK3(0.0), gammaKLimit(0.0),
	gammaD1(0.0), gammaD2(0.0), gammaD3(0.0), gammaDLimit(0.0),
	gammaF1(0.0), gammaF2(0.0), gammaF3(0.0), gammaFLimit(0.0), gammaE(0.0),
	rDispP(0.0), rForceP(0.0), uForceP(0.0), rDispN(0.0), rForceN(0.0), uForceN(0.0)
{

}

Pinching4M::~Pinching4M()
{
#ifdef _G3DEBUG
	fg->close();
#endif
}

void Pinching4M::SetEnvelope(void)
{
	double kPos = stress1p / strain1p;			// stiffness
	double kNeg = stress1n / strain1n;			// stiffness
	double k = (kPos > kNeg) ? kPos : kNeg;		// Get kPos as k if kPos > kNeg
	double u = (strain1p > -strain1n) ? 1e-4 * strain1p : -1e-4 * strain1n;	// Initialize envelope start strain

	envlpPosStrain(0) = u;
	envlpPosStress(0) = u * k;
	envlpNegStrain(0) = -u;
	envlpNegStress(0) = -u * k;

	envlpPosStrain(1) = strain1p;
	envlpPosStrain(2) = strain2p;
	envlpPosStrain(3) = strain3p;
	envlpPosStrain(4) = strain4p;

	envlpNegStrain(1) = strain1n;
	envlpNegStrain(2) = strain2n;
	envlpNegStrain(3) = strain3n;
	envlpNegStrain(4) = strain4n;

	envlpPosStress(1) = stress1p;
	envlpPosStress(2) = stress2p;
	envlpPosStress(3) = stress3p;
	envlpPosStress(4) = stress4p;

	envlpNegStress(1) = stress1n;
	envlpNegStress(2) = stress2n;
	envlpNegStress(3) = stress3n;
	envlpNegStress(4) = stress4n;

	double k1 = (stress4p - stress3p) / (strain4p - strain3p);	// Degradation stiffness at positive branch
	double k2 = (stress4n - stress3n) / (strain4n - strain3n);	// Degradation stiffness at negative branch

	envlpPosStrain(5) = 1e+6 * strain4p;
	envlpPosStress(5) = (k1 > 0.0) ? stress4p + k1 * (envlpPosStrain(5) - strain4p) : stress4p * 1.01; // Change 1.0 to 1.01
	envlpNegStrain(5) = 1e+6 * strain4n;
	envlpNegStress(5) = (k2 > 0.0) ? stress4n + k2 * (envlpNegStrain(5) - strain4n) : stress4n * 1.01; // Change 1.0 to 1.01

	// No pinching branches


	// Define critical material properties
	kElasticPos = envlpPosStress(1) / envlpPosStrain(1); // Elastic stiffness of positive branch
	kElasticNeg = envlpNegStress(1) / envlpNegStrain(1); // Elastic stiffness of negative branch

	double energypos = 0.5 * envlpPosStrain(0) * envlpPosStress(0);

	for (int jt = 0; jt < 4; jt++) {
		energypos += 0.5 * (envlpPosStress(jt) + envlpPosStress(jt + 1)) * (envlpPosStrain(jt + 1) - envlpPosStrain(jt)); // Area
	}

	double energyneg = 0.5 * envlpNegStrain(0) * envlpNegStress(0);

	for (int jy = 0; jy < 4; jy++) {
		energyneg += 0.5 * (envlpNegStress(jy) + envlpNegStress(jy + 1)) * (envlpNegStrain(jy + 1) - envlpNegStrain(jy)); // Area
	}

	double max_energy = (energypos > energyneg) ? energypos : energyneg;
	energyCapacity = gammaE * max_energy; // Total energy dissipation capacity
}


int Pinching4M::revertToStart(void)
{
	// Set converge material parameters
	Cstate = 0;
	Cstrain = 0.0;
	Cstress = 0.0;
	CstrainRate = 0.0;
	lowCstateStrain = envlpNegStrain(0);
	lowCstateStress = envlpNegStress(0);
	hghCstateStrain = envlpPosStrain(0);
	hghCstateStress = envlpPosStress(0);
	CminStrainDmnd = envlpNegStrain(1); // After yielding point damage could happen
	CmaxStrainDmnd = envlpPosStrain(1);
	Cenergy = 0.0;
	CgammaK = 0.0;
	CgammaD = 0.0;
	CgammaF = 0.0;
	CnCycle = 0.0;
	Ttangent = envlpPosStress(0) / envlpPosStrain(0);
	dstrain = 0.0;
	gammaKUsed = 0.0;
	gammaFUsed = 0.0;

	kElasticPosDamgd = kElasticPos;
	kElasticNegDamgd = kElasticNeg;
	uMaxDamgd = CmaxStrainDmnd;
	uMinDamgd = CminStrainDmnd;

	return 0;
}

int Pinching4M::revertToLastCommit(void)
{
	Tstate = Cstate;
	TstrainRate = CstrainRate;
	lowTstateStrain = lowCstateStrain;
	lowTstateStress = lowCstateStress;
	hghTstateStrain = hghCstateStrain;
	hghTstateStress = hghCstateStress;
	TminStrainDmnd = CminStrainDmnd;
	TmaxStrainDmnd = CmaxStrainDmnd;
	Tenergy = Cenergy;
	Tstrain = Cstrain;
	Tstress = Cstress;
	TgammaD = CgammaD;
	TgammaK = CgammaK;
	TgammaF = CgammaF;
	TnCycle = CnCycle;
	return 0;
}

int Pinching4M::setTrialStrain(double strain, double CstrainRate)
{
	// Set trial parameter
	Tstrain = strain; // Get strain
	Tstate = Cstate;
	Tenergy = Cenergy;
	lowTstateStrain = lowCstateStrain;
	hghTstateStrain = hghCstateStrain;
	lowTstateStress = lowCstateStress;
	hghTstateStress = hghCstateStress;
	TminStrainDmnd = CminStrainDmnd;
	TmaxStrainDmnd = CmaxStrainDmnd;
	TgammaF = CgammaF;
	TgammaK = CgammaK;
	TgammaD = CgammaD;

	dstrain = Tstrain - Cstrain; // Get delta strain
	if (dstrain<1e-12 && dstrain>-1e-12) {
		dstrain = 0.0;
	}

	// Determine new state if there is a change in state
	getstate(Tstrain, dstrain);

	switch (Tstate)
	{
	case 0:
		Ttangent = envlpPosStress(0) / envlpPosStrain(0); // Use positive branch elastic stiffness for trial
		Tstress = Ttangent * Tstrain;
		break;
	case 1:
		Tstress = posEnvlpStress(strain); // Use positive branch data
		Ttangent = posEnvlpTangent(strain);
		break;
	case 2:
		Ttangent = negEnvlpTangent(strain); // Use negative branch data
		Tstress = negEnvlpStress(strain);
		break;
	case 3:
		kunload = (hghTstateStrain < 0.0) ? kElasticNegDamgd : kElasticPosDamgd;
		state3Strain(0) = lowTstateStrain; // Update state 3 threshold
		state3Strain(3) = hghTstateStrain;
		state3Stress(0) = lowTstateStress;
		state3Stress(3) = hghTstateStress;

		getState3(state3Strain, state3Stress, kunload); // Update state3 branch
		Ttangent = Envlp3Tangent(state3Strain, state3Stress, strain); // Use state3 branch data
		Tstress = Envlp3Stress(state3Strain, state3Stress, strain);

		break;
	case 4:
		kunload = (lowTstateStrain < 0.0) ? kElasticNegDamgd : kElasticPosDamgd;
		state4Strain(0) = lowTstateStrain;
		state4Strain(3) = hghTstateStrain;
		state4Stress(0) = lowTstateStress;
		state4Stress(3) = hghTstateStress;

		getState4(state4Strain, state4Stress, kunload); // Update state4 branch
		Ttangent = Envlp4Tangent(state4Strain, state4Stress, strain); // Use state4 branch data
		Tstress = Envlp4Stress(state4Strain, state4Stress, strain);
		break;
	}

	double denergy = 0.5 * (Tstress + Cstress) * dstrain;
	elasticStrainEnergy = (Tstrain > 0.0) ? 0.5 * Tstress / kElasticPosDamgd * Tstress : 0.5 * Tstress / kElasticNegDamgd * Tstress;

	Tenergy = Cenergy + denergy; // Update damage

	updateDmg(Tstrain, dstrain);
	return 0;
}

void Pinching4M::getstate(double u, double du) // Trail strain and delta strain
{
	int cid = 0;
	int cis = 0;
	int newState = 0;
	if (du * CstrainRate <= 0.0) { // strainRate = 0.0
		cid = 1;
	}
	if (u < lowTstateStrain || u > hghTstateStrain || cid) { // Strain lager than minimum value
		if (Tstate == 0) { // First running
			if (u > hghTstateStrain) { // Positive branch
				cis = 1;
				newState = 1;
				lowTstateStrain = envlpPosStrain(0); // Use postive branch data
				lowTstateStress = envlpPosStress(0);
				hghTstateStrain = envlpPosStrain(5);
				hghTstateStress = envlpPosStress(5);
			}
			else if (u < lowTstateStrain) { // Negative branch
				cis = 1;
				newState = 2;
				lowTstateStrain = envlpNegStrain(5); // Use negative branch data
				lowTstateStress = envlpNegStress(5);
				hghTstateStrain = envlpNegStrain(0);
				hghTstateStress = envlpNegStress(0);
			}
		}
		else if (Tstate == 1 && du < 0.0) { // Unload from positive branch
			cis = 1;
			if (Cstrain > TmaxStrainDmnd) { // Update TmaxStrainDmnd
				TmaxStrainDmnd = u - du; // Is this correct? du is negative value
			}
			if (TmaxStrainDmnd < uMaxDamgd) { // maxStrainDmnd should less than uMaxDamgd
				TmaxStrainDmnd = uMaxDamgd;
			}
			if (u < uMinDamgd) {	// u less than damage strain
				newState = 2;		// Assign negative branch
				gammaFUsed = CgammaF;
				for (int i = 0; i <= 5; i++) {
					envlpNegDamgdStress(i) = envlpNegStress(i) * (1.0 - gammaFUsed); // Update damaged branch every time
				}
				kPostYieldNegDamgd = (envlpNegDamgdStress(2) - envlpNegDamgdStress(1)) / (envlpNegStrain(2) - envlpNegStrain(1)); // Post-yielding damaged stiffness at negative branch
				lowTstateStrain = envlpNegStrain(5);
				lowTstateStress = envlpNegStress(5);
				hghTstateStrain = envlpNegStrain(0);
				hghTstateStress = envlpNegStress(0);
			}
			else {
				newState = 3; // Negative pinching branch
				lowTstateStrain = uMinDamgd; // Update lowTstateStrain
				gammaFUsed = CgammaF;
				for (int i = 0; i <= 5; i++) {
					envlpNegDamgdStress(i) = envlpNegStress(i) * (1.0 - gammaFUsed); // Update damaged branch every time
				}
				kPostYieldNegDamgd = (envlpNegDamgdStress(2) - envlpNegDamgdStress(1)) / (envlpNegStrain(2) - envlpNegStrain(1));
				lowTstateStress = negEnvlpStress(uMinDamgd); // Update lowTstateStress
				hghTstateStrain = Cstrain; // Set hghTstateStrain as history strain
				hghTstateStress = Cstress; // Set hghTstateStress as history strain
			}
			gammaKUsed = CgammaK;
			kElasticPosDamgd = kElasticPos * (1.0 - gammaKUsed); // Update damaged stiffness
		}
		else if (Tstate == 2 && du > 0.0) { // Reload from negative branch
			cis = 1;
			if (Cstrain < TminStrainDmnd) {
				TminStrainDmnd = Cstrain; // Set history strain as TminStrainDmnd
			}
			if (TminStrainDmnd > uMinDamgd) {
				TminStrainDmnd = uMinDamgd; // TminStrainDmnd should be less than uMinDamgd
			}
			if (u > uMaxDamgd) { // Larger than pinching threshold
				newState = 1; // Positive branch
				gammaFUsed = CgammaF;
				for (int i = 0; i <= 5; i++) {
					envlpPosDamgdStress(i) = envlpPosStress(i) * (1.0 - gammaFUsed); // Update damaged branch every time
				}
				kPostYieldPosDamgd = (envlpPosDamgdStress(2) - envlpPosDamgdStress(1)) / (envlpPosStrain(2) - envlpPosStrain(1));
				lowTstateStrain = envlpPosStrain(0);
				lowTstateStress = envlpPosStress(0);
				hghTstateStrain = envlpPosStrain(5);
				hghTstateStress = envlpPosStress(5);
			}
			else {
				newState = 4; // Positive pinching branch
				lowTstateStrain = Cstrain;
				lowTstateStress = Cstress;
				hghTstateStrain = uMaxDamgd;
				gammaFUsed = CgammaF; // Update gammaFUsed
				for (int i = 0; i <= 5; i++) {
					envlpPosDamgdStress(i) = envlpPosStress(i) * (1.0 - gammaFUsed); // Update damaged envelope
				}
				kPostYieldPosDamgd = (envlpPosDamgdStress(2) - envlpPosDamgdStress(1)) / (envlpPosStrain(2) - envlpPosStrain(1));
				hghTstateStress = posEnvlpStress(uMaxDamgd); // Update hghTstateStress
			}
			gammaKUsed = CgammaK;
			kElasticNegDamgd = kElasticNeg * (1.0 - gammaKUsed); // Update damaged stiffness
		}
		else if (Tstate == 3) {
			if (u < lowTstateStrain) { // u not in the range of pinching strain
				cis = 1;
				newState = 2; // Negative branch
				lowTstateStrain = envlpNegStrain(5); // Update state strain
				hghTstateStrain = envlpNegStrain(0);
				lowTstateStress = envlpNegDamgdStress(5);
				hghTstateStress = envlpNegDamgdStress(0);
			}
			else if (u > uMaxDamgd && du > 0.0) { // u not in the range of pinching strain
				cis = 1;
				newState = 1; // Positive branch
				lowTstateStrain = envlpPosStrain(0); // Update state strain
				lowTstateStress = envlpPosStress(0);
				hghTstateStrain = envlpPosStrain(5);
				hghTstateStress = envlpPosStress(5);
			}
			else if (du > 0.0) { // Positive pinching branch from negative pinching branch
				cis = 1;
				newState = 4;
				lowTstateStrain = Cstrain;
				lowTstateStress = Cstress;
				hghTstateStrain = uMaxDamgd;
				gammaFUsed = CgammaF;
				for (int i = 0; i <= 5; i++) {
					envlpPosDamgdStress(i) = envlpPosStress(i) * (1.0 - gammaFUsed);
				}
				kPostYieldPosDamgd = (envlpPosDamgdStress(2) - envlpPosDamgdStress(1)) / (envlpPosStrain(2) - envlpPosStrain(1));
				hghTstateStress = posEnvlpStress(uMaxDamgd);
				gammaKUsed = CgammaK;
				kElasticNegDamgd = kElasticNeg * (1.0 - gammaKUsed);
			} // Continue on state 3
		}
		else if (Tstate == 4) {
			if (u > hghTstateStrain) { // Not in the range of pinching strain
				cis = 1;
				newState = 1; // Positive branch
				lowTstateStrain = envlpPosStrain(0); // Update state variables
				lowTstateStress = envlpPosDamgdStress(0);
				hghTstateStrain = envlpPosStrain(5);
				hghTstateStress = envlpPosDamgdStress(5);
			}
			else if (u < uMinDamgd && du < 0.0) { // Not in the range of pinching strain
				cis = 1;
				newState = 2; // Negative branch from positive pinching branch
				lowTstateStrain = envlpNegStrain(5); // Update state variables
				lowTstateStress = envlpNegDamgdStress(5);
				hghTstateStrain = envlpNegStrain(0);
				hghTstateStress = envlpNegDamgdStress(0);
			}
			else if (du < 0.0) { // Negative pinching branch from positive pinching branch
				cis = 1;
				newState = 3;
				lowTstateStrain = uMinDamgd;
				gammaFUsed = CgammaF;
				for (int i = 0; i <= 5; i++) {
					envlpNegDamgdStress(i) = envlpNegStress(i) * (1.0 - gammaFUsed);
				}
				kPostYieldNegDamgd = (envlpNegDamgdStress(2) - envlpNegDamgdStress(1)) / (envlpNegStrain(2) - envlpNegStrain(1));
				lowTstateStress = negEnvlpStress(uMinDamgd);
				hghTstateStrain = Cstrain;
				hghTstateStress = Cstress;
				gammaKUsed = CgammaK;
				kElasticPosDamgd = kElasticPos * (1.0 - gammaKUsed);
			} // Continue on state 4
		}
	}
	if (cis) {
		Tstate = newState; // update Tstate
	}
}

double Pinching4M::posEnvlpStress(double u)
{
	double k = 0.0;
	bool isGet = false;
	int i = 0;
	double f = 0.0;
	while (k == 0.0 && i <= 4 && isGet == false) {
		if (u <= envlpPosStrain(i + 1)) { // Look for current position
			k = (envlpPosDamgdStress(i + 1) - envlpPosDamgdStress(i)) / (envlpPosStrain(i + 1) - envlpPosStrain(i));
			f = envlpPosDamgdStress(i) + (u - envlpPosStrain(i)) * k;
			isGet = true;
		}
		i++;
	}
	//if (k == 0.0) {
	//	k = (envlpPosDamgdStress(5) - envlpPosDamgdStress(4)) / (envlpPosStrain(5) - envlpPosStrain(4));
	//	f = envlpPosDamgdStress(5) + k * (u - envlpPosStrain(5));
	//}
	return f;
}

double Pinching4M::posEnvlpTangent(double u)
{
	double k = 0.0;
	bool isGet = false;
	int i = 0;
	while (k == 0.0 && i <= 4 && isGet == false) {
		if (u <= envlpPosStrain(i + 1)) {
			k = (envlpPosDamgdStress(i + 1) - envlpPosDamgdStress(i)) / (envlpPosStrain(i + 1) - envlpPosStrain(i));
			isGet = true;
		}
		i++;
	}
	//if (k == 0.0) {
	//	k = (envlpPosDamgdStress(5) - envlpPosDamgdStress(4)) / (envlpPosStrain(5) - envlpPosStrain(4));
	//}
	return k;
}

double Pinching4M::negEnvlpStress(double u)
{
	double k = 0.0;
	bool isGet = false;
	int i = 0;
	double f = 0.0;
	while (k == 0.0 && i <= 4 && isGet == false) {
		if (u >= envlpNegStrain(i + 1)) { // Look for current position
			k = (envlpNegDamgdStress(i) - envlpNegDamgdStress(i + 1)) / (envlpNegStrain(i) - envlpNegStrain(i + 1));
			f = envlpNegDamgdStress(i + 1) + (u - envlpNegStrain(i + 1)) * k;
			isGet = true;
		}
		i++;
	}
	//if (k == 0.0) {
	//	k = (envlpNegDamgdStress(4) - envlpNegDamgdStress(5)) / (envlpNegStrain(4) - envlpNegStrain(5));
	//	f = envlpNegDamgdStress(5) + k * (u - envlpNegStrain(5));
	//}
	return f;
}

double Pinching4M::negEnvlpTangent(double u)
{
	double k = 0.0;
	bool isGet = false;
	int i = 0;
	while (k == 0.0 && i <= 4 && isGet == false) {
		if (u >= envlpNegStrain(i + 1)) {
			k = (envlpNegDamgdStress(i) - envlpNegDamgdStress(i + 1)) / (envlpNegStrain(i) - envlpNegStrain(i + 1));
			isGet = true;
		}
		i++;
	}
	/*if (k == 0.0) {
		k = (envlpNegDamgdStress(4) - envlpNegDamgdStress(5)) / (envlpNegStrain(4) - envlpNegStrain(5));
	}*/
	return k;
}

void Pinching4M::getState3(Vector& state3Strain, Vector& state3Stress, double kunload)
{
	// Calculation of the unloading point
	double strainNegUnlStart = hghTstateStrain;         // Strain at unloading point
	double stressNegUnlStart = hghTstateStress;         // Stress at unloading point
	double strainNegY = envlpNegStrain(1);          // Strain at yielding point
	double stressNegY = envlpNegDamgdStress(1); // Stress at yielding point
	double unloadNegStressFul = kunload / (kunload - kPostYieldNegDamgd) * (kPostYieldNegDamgd * (strainNegUnlStart - strainNegY - stressNegUnlStart / kunload) + stressNegY);
	double unloadNegStress = unloadNegStressFul * uForceN;
	// Check unloadNegStress
	if (unloadNegStress > 0) {
		unloadNegStress = 0.0;
	}
	state3Stress(2) = unloadNegStress;
	state3Strain(2) = hghTstateStrain - (hghTstateStress - state3Stress(2)) / kunload;

	// Calculation of the reloading point
	double strainS3Max = lowTstateStrain; // Maximum point at negative pinching branch
	double stressS3Max = lowTstateStress;
	double stressS3RelMin = unloadNegStress; // Minimum point for reloading
	double strainS3RelMin = stressS3RelMin / (stressNegY / strainNegY);
	double stressS3RelMax = unloadNegStress; // Maximum point for reloading
	double strainS3RelMax = (stressS3RelMax - stressS3Max) / kunload + strainS3Max;
	double reloadNegStress = rForceN * (stressS3Max - unloadNegStress) + unloadNegStress; // Stress at the reloading point
	double reloadNegStrainMin = (reloadNegStress - stressS3RelMin) / (stressNegY / strainNegY) + strainS3RelMin; // Strain for the start point of reloading point
	double reloadNegStrainMax = (reloadNegStress - stressS3RelMax) / kunload + strainS3RelMax; // Strain for the end point of reloading
	/*if (reloadNegStress < stressNegY) {
		reloadNegStrainMin = (reloadNegStress - stressNegY) / kPostYieldNegDamgd + strainNegY;
	}*/
	double reloadNegStrain = rDispN * (reloadNegStrainMax - reloadNegStrainMin) + reloadNegStrainMin;
	state3Strain(1) = reloadNegStrain;
	state3Stress(1) = reloadNegStress;

	// Check the reloading point
	if ((abs(state3Strain(0) - state3Strain(1)) < 1e-8) && (abs(state3Stress(0) == state3Stress(1)) < 1e-8)) { // Same point
		double du = state3Strain(0) - state3Strain(2);
		double df = state3Stress(0) - state3Stress(2);
		state4Strain(1) = state4Strain(2) + 0.5 * du;
		state4Stress(1) = state4Stress(2) + 0.5 * df;
	}


	/*double kmax = (kunload > kElasticNegDamgd) ? kunload : kElasticNegDamgd;
	if (state3Strain(0) * state3Strain(3) < 0.0) { // Trilinear unload reload path expected, first define point for reloading
		state3Strain(1) = lowTstateStrain * rDispN;

		if (rForceN - uForceN > 1e-8) { // rForceN > uForceN
			state3Stress(1) = lowTstateStress * rForceN;
		}
		else { // rForceN < uForceN
			if (TminStrainDmnd < envlpNegStrain(3)) {
				double st1 = lowTstateStress * uForceN * (1.0 + 1e-6);
				double st2 = envlpNegDamgdStress(4) * (1.0 + 1e-6);
				state3Stress(1) = (st1 < st2) ? st1 : st2;
			}
			else {
				double st1 = envlpNegDamgdStress(3) * uForceN * (1.0 + 1e-6);
				double st2 = envlpNegDamgdStress(4) * (1.0 + 1e-6);
				state3Stress(1) = (st1 < st2) ? st1 : st2;
			}
		}

		// If reload stiffness exceeds unload stiffness, reduce reload stiffness to make it equal to unload stiffness
		if ((state3Stress(1) - state3Stress(0)) / (state3Strain(1) - state3Strain(0)) > kElasticNegDamgd) {
			state3Strain(1) = lowTstateStrain + (state3Stress(1) - state3Stress(0)) / kElasticNegDamgd;
		}

		// Check that reloading point is not behind point 4
		if (state3Strain(1) > state3Strain(3)) {
			// Path taken to be a straight line between points 1 and 4
			double du = state3Strain(3) - state3Strain(0);
			double df = state3Stress(3) - state3Stress(0);
			state3Strain(1) = state3Strain(0) + 0.33 * du;
			state3Strain(2) = state3Strain(0) + 0.67 * du;
			state3Stress(1) = state3Stress(0) + 0.33 * df;
			state3Stress(2) = state3Stress(0) + 0.67 * df;
		}
		else {
			if (TminStrainDmnd < envlpNegStrain(3)) {
				state3Stress(2) = uForceN * envlpNegDamgdStress(4);
				//state3Stress(2) = uForceN * envlpDamgdStressUnlodFroNegFul(4);
			}
			else {
				state3Stress(2) = uForceN * envlpNegDamgdStress(3);
				//state3Stress(2) = uForceN * envlpDamgdStressUnlodFroNegFul(3);
			}
			state3Strain(2) = hghTstateStrain - (hghTstateStress - state3Stress(2)) / kunload;

			if (state3Strain(2) > state3Strain(3)) {
				// Point3 should be along a line between 2 and 4
				double du = state3Strain(3) - state3Strain(1);
				double df = state3Stress(3) - state3Stress(1);
				state3Strain(2) = state3Strain(1) + 0.5 * du;
				state3Stress(2) = state3Stress(1) + 0.5 * df;
			}
			else if ((state3Stress(2) - state3Stress(1)) / (state3Strain(2) - state3Strain(1)) > kmax) {
				// Linear unload-reload path expected
				double du = state3Strain(3) - state3Strain(0);
				double df = state3Stress(3) - state3Stress(0);
				state3Strain(1) = state3Strain(0) + 0.33 * du;
				state3Strain(2) = state3Strain(0) + 0.67 * du;
				state3Stress(1) = state3Stress(0) + 0.33 * df;
				state3Stress(2) = state3Stress(0) + 0.67 * df;
			}
			else if ((state3Strain(2) < state3Strain(1)) || ((state3Stress(2) - state3Stress(1)) / (state3Strain(2) - state3Strain(1)) < 0)) {
				if (state3Strain(2) < 0.0) {
					// Point 3 should be along a line between 2 and 4
					double du = state3Strain(3) - state3Strain(1);
					double df = state3Stress(3) - state3Stress(1);
					state3Strain(2) = state3Strain(1) + 0.5 * du;
					state3Stress(2) = state3Stress(1) + 0.5 * df;
				}
				else if (state3Strain(1) > 0.0) {
					// Point 2 should be along a line between 1 and 3
					double du = state3Strain(2) - state3Strain(0);
					double df = state3Stress(2) - state3Stress(0);
					state3Strain(1) = state3Strain(0) + 0.5 * du;
					state3Stress(1) = state3Stress(0) + 0.5 * df;
				}
				else {
					double avgforce = 0.5 * (state3Stress(2) + state3Stress(1));
					double dfr = 0.0;
					if (avgforce < 0.0) {
						dfr = -avgforce / 100;
					}
					else {
						dfr = avgforce / 100;
					}
					double slope12 = (state3Stress(1) - state3Stress(0)) / (state3Strain(1) - state3Strain(0));
					double slope34 = (state3Stress(3) - state3Stress(2)) / (state3Strain(3) - state3Strain(2));
					state3Stress(1) = avgforce - dfr;
					state3Stress(2) = avgforce + dfr;
					state3Strain(1) = state3Strain(0) + (state3Stress(1) - state3Stress(0)) / slope12;
					state3Strain(2) = state3Strain(3) - (state3Stress(3) - state3Stress(2)) / slope34;
				}
			}
		}

	}
	else { // Linear unload reload path is expected
		double du = state3Strain(3) - state3Strain(0);
		double df = state3Stress(3) - state3Stress(0);
		state3Strain(1) = state3Strain(0) + 0.33 * du;
		state3Strain(2) = state3Strain(0) + 0.67 * du;
		state3Stress(1) = state3Stress(0) + 0.33 * df;
		state3Stress(2) = state3Stress(0) + 0.67 * df;
	}
	double checkSlope = state3Stress(0) / state3Strain(0);
	double slope = 0.0;
	// Final check
	int i = 0;
	while (i < 3) {
		double du = state3Strain(i + 1) - state3Strain(i);
		double df = state3Stress(i + 1) - state3Stress(i);
		if (du < 0.0 || df < 0.0) {
			double du = state3Strain(3) - state3Strain(0);
			double df = state3Stress(3) - state3Stress(0);
			state3Strain(1) = state3Strain(0) + 0.33 * du;
			state3Strain(2) = state3Strain(0) + 0.67 * du;
			state3Stress(1) = state3Stress(0) + 0.33 * df;
			state3Stress(2) = state3Stress(0) + 0.67 * df;
			slope = df / du;
			i = 3;
		}
		if (slope > 1e-8 && slope < checkSlope) {
			state3Strain(1) = 0.0;
			state3Stress(1) = 0.0;
			state3Strain(2) = state3Strain(3) / 2;
			state3Stress(2) = state3Stress(3) / 2;
		}
		i++;
	}*/
}

void Pinching4M::getState4(Vector& state4Strain, Vector& state4Stress, double kunload)
{
	// Calculation of the unloading point
	double strainPosUnlStart = lowTstateStrain;         // Strain at unloading point
	double stressPosUnlStart = lowTstateStress;         // Stress at unloading point
	double strainPosY = envlpPosStrain(1);          // Strain at yielding point
	double StressPosY = envlpPosDamgdStress(1); // Stress at yielding point
	double unloadPosStressFul = kunload / (kunload - kPostYieldPosDamgd) * (kPostYieldPosDamgd * (strainPosUnlStart - strainPosY - stressPosUnlStart / kunload) + StressPosY);
	double unloadPosStress = unloadPosStressFul * uForceP;
	// Check unloadNegStress
	if (unloadPosStress < 0) {
		unloadPosStress = 0.0;
	}
	state4Stress(1) = unloadPosStress;
	state4Strain(1) = lowTstateStrain + (-lowTstateStress + state4Stress(1)) / kunload;

	// Calculation of the reloading point
	double strainS4Max = hghTstateStrain; // Maximum point at positive pinching branch
	double stressS4Max = hghTstateStress;
	double stressS4RelMin = unloadPosStress; // Minimum point for reloading
	double strainS4RelMin = stressS4RelMin / (StressPosY / strainPosY);
	double stressS4RelMax = unloadPosStress; // Maximum point for reloading
	double strainS4RelMax = (stressS4RelMax - stressS4Max) / kunload + strainS4Max;
	double reloadPosStress = rForceN * (stressS4Max - unloadPosStress) + unloadPosStress; // Stress at the reloading point
	double reloadPosStrainMin = (reloadPosStress - stressS4RelMin) / (StressPosY / strainPosY) + strainS4RelMin;
	double reloadPosStrainMax = (reloadPosStress - stressS4RelMax) / kunload + strainS4RelMax;
	/*if (reloadPosStress > StressPosY) {
		reloadPosStrainMin = (reloadPosStress - StressPosY) / kPostYieldPosDamgd + strainPosY;
	}*/
	double reloadPosStrain = rDispN * (reloadPosStrainMax - reloadPosStrainMin) + reloadPosStrainMin;
	state4Strain(2) = reloadPosStrain;
	state4Stress(2) = reloadPosStress;
	
	// Check the reloading point
	if ((abs(state4Strain(2) - state4Strain(3)) < 1e-8) && (abs(state4Stress(2) - state4Stress(3)) < 1e-8)) { // Same point
		double du = state4Strain(3) - state4Strain(1);
		double df = state4Stress(3) - state4Stress(1);
		state4Strain(2) = state4Strain(1) + 0.5 * du;
		state4Stress(2) = state4Stress(1) + 0.5 * df;
	}

	/*
	double kmax = (kunload > kElasticPosDamgd) ? kunload : kElasticPosDamgd;
	if (state4Strain(0) * state4Strain(3) < 0.0) {
		// Trilinear unload reload path epected
		state4Strain(2) = hghTstateStrain * rDispP;
		if (uForceP == 0.0) {
			state4Stress(2) = hghTstateStress * rForceP;
		}
		else if (rForceP - uForceP > 1e-8) {
			state4Stress(2) = hghTstateStress * rForceP;
		}
		else {
			if (TmaxStrainDmnd > envlpPosStrain(3)) {
				double st1 = hghTstateStress * uForceP * (1.0 + 1e-6);
				double st2 = envlpPosDamgdStress(4) * (1.0 + 1e-6);
				state4Stress(2) = (st1 > st2) ? st1 : st2;
			}
			else {
				double st1 = envlpPosDamgdStress(3) * uForceP * (1.0 + 1e-6);
				double st2 = envlpPosDamgdStress(4) * (1.0 + 1e-6);
				state4Stress(2) = (st1 > st2) ? st1 : st2;
			}
		}
		// If reload stiffness exceeds unload stiffness, reduce reload stiffness to make it equal to unload stiffness
		if ((state4Stress(3) - state4Stress(2)) / (state4Strain(3) - state4Strain(2)) > kElasticPosDamgd) {
			state4Strain(2) = hghTstateStrain - (state4Stress(3) - state4Stress(2)) / kElasticPosDamgd;
		}
		// Check that reloading point is not behind point 1
		if (state4Strain(2) < state4Strain(0)) {
			// Path taken to be a straight line between points 1 and 4
			double du = state4Strain(3) - state4Strain(0);
			double df = state4Stress(3) - state4Stress(0);
			state4Strain(1) = state4Strain(0) + 0.33 * du;
			state4Strain(2) = state4Strain(0) + 0.67 * du;
			state4Stress(1) = state4Stress(0) + 0.33 * df;
			state4Stress(2) = state4Stress(0) + 0.67 * df;
		}
		else {
			if (TmaxStrainDmnd > envlpPosStrain(3)) {
				state4Stress(1) = uForceP * envlpPosDamgdStress(4);
			}
			else {
				state4Stress(1) = uForceP * envlpPosDamgdStress(3);
			}
			state4Strain(1) = lowTstateStrain + (-lowTstateStress + state4Stress(1)) / kunload;
			if (state4Strain(1) < state4Strain(0)) {
				// Point 2 should be along a line between 1 and 3
				double du = state4Strain(2) - state4Strain(0);
				double df = state4Stress(2) - state4Stress(0);
				state4Strain(1) = state4Strain(0) + 0.5 * du;
				state4Stress(1) = state4Stress(0) + 0.5 * df;
			}
			else if ((state4Stress(2) - state4Stress(1)) / (state4Strain(2) - state4Strain(1)) > kmax) {
				// linear unload-reload path expected
				double du = state4Strain(3) - state4Strain(0);
				double df = state4Stress(3) - state4Stress(0);
				state4Strain(1) = state4Strain(0) + 0.33 * du;
				state4Strain(2) = state4Strain(0) + 0.67 * du;
				state4Stress(1) = state4Stress(0) + 0.33 * df;
				state4Stress(2) = state4Stress(0) + 0.67 * df;
			}
			else if ((state4Strain(2) < state4Strain(1)) || ((state4Stress(2) - state4Stress(1)) / (state4Strain(2) - state4Strain(1)) < 0)) {
				if (state4Strain(1) > 0.0) {
					// Point 2 should be along a line between 1 and 3
					double du = state4Strain(2) - state4Strain(0);
					double df = state4Stress(2) - state4Stress(0);
					state4Strain(1) = state4Strain(0) + 0.5 * du;
					state4Stress(1) = state4Stress(0) + 0.5 * df;
				}
				else if (state4Strain(2) < 0.0) {
					// Point 2 should be along a line between 2 and 4
					double du = state4Strain(3) - state4Strain(1);
					double df = state4Stress(3) - state4Stress(1);
					state4Strain(2) = state4Strain(1) + 0.5 * du;
					state4Stress(2) = state4Stress(1) + 0.5 * df;
				}
				else {
					double avgforce = 0.5 * (state4Stress(2) + state4Stress(1));
					double dfr = 0.0;
					if (avgforce < 0.0) {
						dfr = -avgforce / 100;
					}
					else {
						dfr = avgforce / 100;
					}
					double slope12 = (state4Stress(1) - state4Stress(0)) / (state4Strain(1) - state4Strain(0));
					double slope34 = (state4Stress(3) - state4Stress(2)) / (state4Strain(3) - state4Strain(2));
					state4Stress(1) = avgforce - dfr;
					state4Stress(2) = avgforce + dfr;
					state4Strain(1) = state4Strain(0) + (state4Stress(1) - state4Stress(0)) / slope12;
					state4Strain(2) = state4Strain(3) - (state4Stress(3) - state4Stress(2)) / slope34;
				}
			}
		}
	}
	else {
		// Linear unload reload path is expected
		double du = state4Strain(3) - state4Strain(0);
		double df = state4Stress(3) - state4Stress(0);
		state4Strain(1) = state4Strain(0) + 0.33 * du;
		state4Strain(2) = state4Strain(0) + 0.67 * du;
		state4Stress(1) = state4Stress(0) + 0.33 * df;
		state4Stress(2) = state4Stress(0) + 0.67 * df;
	}

	double checkSlope = state4Stress(0) / state4Strain(0);
	double slope = 0.0;

	// Final check
	int i = 0;
	while (i < 3) {
		double du = state4Strain(i + 1) - state4Strain(i);
		double df = state4Stress(i + 1) - state4Stress(i);
		if (du < 0.0 || df < 0.0) {
			double du = state4Strain(3) - state4Strain(0);
			double df = state4Stress(3) - state4Stress(0);
			state4Strain(1) = state4Strain(0) + 0.33 * du;
			state4Strain(2) = state4Strain(0) + 0.67 * du;
			state4Stress(1) = state4Stress(0) + 0.33 * df;
			state4Stress(2) = state4Stress(0) + 0.67 * df;
			slope = df / du;
			i = 3;
		}
		if (slope > 1e-8 && slope < checkSlope) {
			state4Strain(1) = 0.0;
			state4Stress(1) = 0.0;
			state4Strain(2) = state4Strain(3) / 2;
			state4Stress(2) = state4Stress(3) / 2;
		}
		i++;
	}*/
}


double Pinching4M::Envlp3Tangent(Vector s3Strain, Vector s3Stress, double u)
{
	double k = 0.0;
	int i = 0;
	while ((k == 0.0 || i <= 2) && (i <= 2))
	{
		if (u >= s3Strain(i)) {
			k = (s3Stress(i + 1) - s3Stress(i)) / (s3Strain(i + 1) - s3Strain(i));
		}
		i++;
	}
	//if (k == 0.0) {
	//	if (u < s3Strain(0)) {
	//		i = 0;
	//	}
	//	else {
	//		i = 2;
	//	}
	//	k = (s3Stress(i + 1) - s3Stress(i)) / (s3Strain(i + 1) - s3Strain(i));
	//}
	return k;
}

double Pinching4M::Envlp3Stress(Vector s3Strain, Vector s3Stress, double u)
{
	double k = 0.0;
	int i = 0;
	double f = 0.0;
	while ((k == 0.0 || i <= 2) && (i <= 2))
	{
		if (u >= s3Strain(i)) {
			k = (s3Stress(i + 1) - s3Stress(i)) / (s3Strain(i + 1) - s3Strain(i));
			f = s3Stress(i) + (u - s3Strain(i)) * k;
		}
		i++;
	}
	//if (k == 0.0) {
	//	if (u < s3Strain(0)) {
	//		i = 0;
	//	}
	//	else {
	//		i = 2;
	//	}
	//	k = (s3Stress(i + 1) - s3Stress(i)) / (s3Strain(i + 1) - s3Strain(i));
	//	f = s3Stress(i) + (u - s3Strain(i)) * k;
	//}
	return f;
}

double Pinching4M::Envlp4Tangent(Vector s4Strain, Vector s4Stress, double u)
{
	double k = 0.0;
	int i = 0;
	while ((k == 0.0 || i <= 2) && (i <= 2))
	{
		if (u >= s4Strain(i)) {
			k = (s4Stress(i + 1) - s4Stress(i)) / (s4Strain(i + 1) - s4Strain(i));
		}
		i++;
	}
	//if (k == 0.0) {
	//	if (u < s4Strain(0)) {
	//		i = 0;
	//	}
	//	else {
	//		i = 2;
	//	}
	//	k = (s4Stress(i + 1) - s4Stress(i)) / (s4Strain(i + 1) - s4Strain(i));
	//}
	return k;
}

double Pinching4M::Envlp4Stress(Vector s4Strain, Vector s4Stress, double u)
{
	double k = 0.0;
	int i = 0;
	double f = 0.0;
	while ((k == 0.0 || i <= 2) && (i <= 2))
	{
		if (u >= s4Strain(i)) {
			k = (s4Stress(i + 1) - s4Stress(i)) / (s4Strain(i + 1) - s4Strain(i));
			f = s4Stress(i) + (u - s4Strain(i)) * k;
		}
		i++;
	}
	//if (k == 0.0) {
	//	if (u < s4Strain(0)) {
	//		i = 0;
	//	}
	//	else {
	//		i = 2;
	//	}
	//	k = (s4Stress(i + 1) - s4Stress(i)) / (s4Strain(i + 1) - s4Strain(i));
	//	f = s4Stress(i) + (u - s4Strain(i)) * k;
	//}
	return f;
}



void Pinching4M::updateDmg(double strain, double dstrain)
{
	double tes = 0.0;
	double umaxAbs = (TmaxStrainDmnd > -TminStrainDmnd) ? TmaxStrainDmnd : -TminStrainDmnd; // Maximum historic deformation demand
	double uultAbs = (envlpPosStrain(4) > -envlpNegStrain(4)) ? envlpPosStrain(4) : -envlpNegStrain(4); // Positive / negative deformations that define failure
	TnCycle = CnCycle + fabs(dstrain) / (4 * umaxAbs);
	if ((strain < uultAbs && strain >-uultAbs) && Tenergy < energyCapacity) // Under the range of damage
	{
		TgammaK = gammaK1 * pow((umaxAbs / uultAbs), gammaK3);
		TgammaD = gammaD1 * pow((umaxAbs / uultAbs), gammaD3);
		TgammaF = gammaF1 * pow((umaxAbs / uultAbs), gammaF3);

		if (Tenergy > elasticStrainEnergy && DmgCyc == 0) {
			tes = ((Tenergy - elasticStrainEnergy) / energyCapacity);
			TgammaK = TgammaK + gammaK2 * pow(tes, gammaK4);
			TgammaD = TgammaD + gammaD2 * pow(tes, gammaD4);
			TgammaF = TgammaF + gammaF2 * pow(tes, gammaF4);
		}
		else if (DmgCyc == 1) {
			TgammaK = TgammaK + gammaK2 * pow(TnCycle, gammaK4);
			TgammaD = TgammaD + gammaD2 * pow(TnCycle, gammaD4);
			TgammaF = TgammaF + gammaF2 * pow(TnCycle, gammaF4);
		}
		double kminP = (posEnvlpStress(TmaxStrainDmnd) / TmaxStrainDmnd);
		double kminN = (negEnvlpStress(TminStrainDmnd) / TminStrainDmnd);
		double kmin = ((kminP / kElasticPos) > (kminN / kElasticNeg)) ? (kminP / kElasticPos) : (kminN / kElasticNeg);
		double gammaKLimEnv = (0.0 > (1.0 - kmin)) ? 0.0 : (1.0 - kmin);

		double k1 = (TgammaK < gammaKLimit) ? TgammaK : gammaKLimit;
		TgammaK = (k1 < gammaKLimEnv) ? k1 : gammaKLimEnv;
		TgammaD = (TgammaD < gammaDLimit) ? TgammaD : gammaDLimit;
		TgammaF = (TgammaF < gammaFLimit) ? TgammaF : gammaFLimit;
	}
	else if (strain < uultAbs && strain > -uultAbs) { // Failure
		double kminP = (posEnvlpStress(TmaxStrainDmnd) / TmaxStrainDmnd);
		double kminN = (negEnvlpStress(TminStrainDmnd) / TminStrainDmnd);
		double kmin = ((kminP / kElasticPos) >= (kminN / kElasticNeg)) ? (kminP / kElasticPos) : (kminN / kElasticNeg);
		double gammaKLimEnv = (0.0 > (1.0 - kmin)) ? 0.0 : (1.0 - kmin);

		TgammaK = (gammaKLimit < gammaKLimEnv) ? gammaKLimit : gammaKLimEnv;
		TgammaD = gammaDLimit;
		TgammaF = gammaFLimit;
	}
}

double Pinching4M::getStrain(void)
{
	return Tstrain;
}

double Pinching4M::getStress(void)
{
	return Tstress;
}

double Pinching4M::getTangent(void)
{
	return Ttangent;
}

double Pinching4M::getInitialTangent(void)
{
	return envlpPosStress(0) / envlpPosStrain(0);
}

int Pinching4M::commitState(void) {
	Cstate = Tstate;

	if (dstrain > 1e-12 || dstrain < -(1e-12)) {
		CstrainRate = dstrain;
	}
	else {
		CstrainRate = TstrainRate;
	}

	lowCstateStrain = lowTstateStrain;
	lowCstateStress = lowTstateStress;
	hghCstateStrain = hghTstateStrain;
	hghCstateStress = hghTstateStress;
	CminStrainDmnd = TminStrainDmnd;
	CmaxStrainDmnd = TmaxStrainDmnd;
	Cenergy = Tenergy;

	Cstress = Tstress;
	Cstrain = Tstrain;

	CgammaK = TgammaK;
	CgammaD = TgammaD;
	CgammaF = TgammaF;

	// Define adjusted strength and stiffness parameters
	kElasticPosDamgd = kElasticPos * (1 - gammaKUsed);
	kElasticNegDamgd = kElasticNeg * (1 - gammaKUsed);

	uMaxDamgd = TmaxStrainDmnd * (1 + CgammaD);
	uMinDamgd = TminStrainDmnd * (1 + CgammaD);

	envlpPosDamgdStress = envlpPosStress * (1 - gammaFUsed);
	envlpNegDamgdStress = envlpNegStress * (1 - gammaFUsed);

	CnCycle = TnCycle; // Number of cycles of loading

#ifdef _G3DEBUG
	(*fg) << tagMat << " " << CgammaF << " " << CgammaK << " " << CgammaD << endln;
#endif
	return 0;
}

UniaxialMaterial* Pinching4M::getCopy(void)
{
	Pinching4M* theCopy = new Pinching4M(this->getTag(),
		stress1p, strain1p, stress2p, strain2p, stress3p, strain3p, stress4p, strain4p,
		stress1n, strain1n, stress2n, strain2n, stress3n, strain3n, stress4n, strain4n,
		rDispP, rForceP, uForceP, rDispN, rForceN, uForceN,
		gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit,
		gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit,
		gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE, DmgCyc);

	theCopy->rDispN = rDispN;
	theCopy->rDispP = rDispP;
	theCopy->rForceN = rForceN;
	theCopy->rForceP = rForceP;
	theCopy->uForceN = uForceN;
	theCopy->uForceP = uForceP;

	// Trial state variables
	theCopy->Tstress = Tstress;
	theCopy->Tstrain = Tstrain;
	theCopy->Ttangent = Ttangent;

	// Coverged material history parameters
	theCopy->Cstate = Cstate;
	theCopy->Cstrain = Cstrain;
	theCopy->Cstress = Cstress;
	theCopy->CstrainRate = CstrainRate;

	theCopy->lowCstateStrain = lowCstateStrain;
	theCopy->lowCstateStress = lowCstateStress;
	theCopy->hghCstateStrain = hghCstateStrain;
	theCopy->hghCstateStress = hghCstateStress;
	theCopy->CminStrainDmnd = CminStrainDmnd;
	theCopy->CmaxStrainDmnd = CmaxStrainDmnd;
	theCopy->Cenergy = Cenergy;
	theCopy->CgammaK = CgammaK;
	theCopy->CgammaD = CgammaD;
	theCopy->CgammaF = CgammaF;
	theCopy->CnCycle = CnCycle;
	theCopy->gammaKUsed = gammaKUsed;
	theCopy->gammaFUsed = gammaFUsed;
	theCopy->DmgCyc = DmgCyc;

	// Trial material history parameters
	theCopy->Tstate = Tstate;
	theCopy->dstrain = dstrain;
	theCopy->lowTstateStrain = lowTstateStrain;
	theCopy->lowTstateStress = lowTstateStress;
	theCopy->hghTstateStrain = hghTstateStrain;
	theCopy->hghTstateStress = hghTstateStress;
	theCopy->TminStrainDmnd = TminStrainDmnd;
	theCopy->TmaxStrainDmnd = TmaxStrainDmnd;
	theCopy->Tenergy = Tenergy;
	theCopy->TgammaK = TgammaK;
	theCopy->TgammaD = TgammaD;
	theCopy->TgammaF = TgammaF;
	theCopy->TnCycle = TnCycle;

	// Strength and stiffness parameters
	theCopy->kElasticPos = kElasticPos;
	theCopy->kElasticNeg = kElasticNeg;
	theCopy->kElasticPosDamgd = kElasticPosDamgd;
	theCopy->kElasticNegDamgd = kElasticNegDamgd;
	theCopy->uMaxDamgd = uMaxDamgd;
	theCopy->uMinDamgd = uMinDamgd;

	// Yielding stress
	theCopy->stressyp = strain1p;
	theCopy->stressyn = strain1n;

	for (int i = 0; i < 6; i++)
	{
		theCopy->envlpPosStrain(i) = envlpPosStrain(i);
		theCopy->envlpPosStress(i) = envlpPosStress(i);
		theCopy->envlpNegStrain(i) = envlpNegStrain(i);
		theCopy->envlpNegStress(i) = envlpNegStress(i);
		theCopy->envlpNegDamgdStress(i) = envlpNegDamgdStress(i);
		theCopy->envlpPosDamgdStress(i) = envlpPosDamgdStress(i);
	}

	for (int j = 0; j < 4; j++)
	{
		theCopy->state3Strain(j) = state3Strain(j);
		theCopy->state3Stress(j) = state3Stress(j);
		theCopy->state4Strain(j) = state4Strain(j);
		theCopy->state4Stress(j) = state4Stress(j);
	}

	theCopy->energyCapacity = energyCapacity;
	theCopy->kunload = kunload;
	theCopy->elasticStrainEnergy = elasticStrainEnergy;

	return theCopy;
}

int Pinching4M::sendSelf(int commitTag, Channel& theChannel)
{
	return -1;
}

int Pinching4M::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	return -1;
}

void Pinching4M::Print(OPS_Stream& s, int flag)
{
	s << "Pinching4M, tag: " << this->getTag() << endln;
	s << "strain: " << Tstrain << endln;
	s << "stress: " << Tstress << endln;
	s << "state: " << Tstate << endln;
}

int Pinching4M::setParameter(const char** argv, int argc, Parameter& param) {
	// Parameters for backbone control points
	if (strcmp(argv[0], "f1p") == 0 || strcmp(argv[0], "stress1p") == 0) {
		param.setValue(stress1p);
		return param.addObject(1, this);
	}
	if (strcmp(argv[0], "d1p") == 0 || strcmp(argv[0], "strain1p") == 0) {
		param.setValue(strain1p);
		return param.addObject(2, this);
	}
	if (strcmp(argv[0], "f2p") == 0 || strcmp(argv[0], "stress2p") == 0) {
		param.setValue(stress2p);
		return param.addObject(3, this);
	}
	if (strcmp(argv[0], "d2p") == 0 || strcmp(argv[0], "strain2p") == 0) {
		param.setValue(strain2p);
		return param.addObject(4, this);
	}
	if (strcmp(argv[0], "f3p") == 0 || strcmp(argv[0], "stress3p") == 0) {
		param.setValue(stress3p);
		return param.addObject(5, this);
	}
	if (strcmp(argv[0], "d3p") == 0 || strcmp(argv[0], "strain3p") == 0) {
		param.setValue(strain3p);
		return param.addObject(6, this);
	}
	if (strcmp(argv[0], "f4p") == 0 || strcmp(argv[0], "stress4p") == 0) {
		param.setValue(stress4p);
		return param.addObject(7, this);
	}
	if (strcmp(argv[0], "d4p") == 0 || strcmp(argv[0], "strain4p") == 0) {
		param.setValue(strain4p);
		return param.addObject(8, this);
	}
	if (strcmp(argv[0], "f1n") == 0 || strcmp(argv[0], "stress1n") == 0) {
		param.setValue(stress1n);
		return param.addObject(9, this);
	}
	if (strcmp(argv[0], "d1n") == 0 || strcmp(argv[0], "strain1n") == 0) {
		param.setValue(strain1n);
		return param.addObject(10, this);
	}
	if (strcmp(argv[0], "f2n") == 0 || strcmp(argv[0], "stress2n") == 0) {
		param.setValue(stress2n);
		return param.addObject(11, this);
	}
	if (strcmp(argv[0], "d2n") == 0 || strcmp(argv[0], "strain2n") == 0) {
		param.setValue(strain2n);
		return param.addObject(12, this);
	}
	if (strcmp(argv[0], "f3n") == 0 || strcmp(argv[0], "stress3n") == 0) {
		param.setValue(stress3n);
		return param.addObject(13, this);
	}
	if (strcmp(argv[0], "d3n") == 0 || strcmp(argv[0], "strain3n") == 0) {
		param.setValue(strain3n);
		return param.addObject(14, this);
	}
	if (strcmp(argv[0], "f4n") == 0 || strcmp(argv[0], "stress4n") == 0) {
		param.setValue(stress4n);
		return param.addObject(15, this);
	}
	if (strcmp(argv[0], "d4n") == 0 || strcmp(argv[0], "strain4n") == 0) {
		param.setValue(strain4n);
		return param.addObject(16, this);
	}

	// Parameters for hysteretic rules
	if (strcmp(argv[0], "rDispP") == 0) {
		param.setValue(rDispP);
		return param.addObject(17, this);
	}
	if (strcmp(argv[0], "rForceP") == 0) {
		param.setValue(rForceP);
		return param.addObject(18, this);
	}
	if (strcmp(argv[0], "uForceP") == 0) {
		param.setValue(uForceP);
		return param.addObject(19, this);
	}
	if (strcmp(argv[0], "rDispN") == 0) {
		param.setValue(rDispN);
		return param.addObject(20, this);
	}
	if (strcmp(argv[0], "rForceN") == 0) {
		param.setValue(rForceN);
		return param.addObject(21, this);
	}
	if (strcmp(argv[0], "uForceN") == 0) {
		param.setValue(uForceN);
		return param.addObject(22, this);
	}

	return -1;
}

int Pinching4M::updateParameter(int parameterID, Information& info)
{
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		this->stress1p = info.theDouble;
		break;
	case 2:
		this->strain1p = info.theDouble;
		break;
	case 3:
		this->stress2p = info.theDouble;
		break;
	case 4:
		this->strain2p = info.theDouble;
		break;
	case 5:
		this->stress3p = info.theDouble;
		break;
	case 6:
		this->strain3p = info.theDouble;
		break;
	case 7:
		this->stress4p = info.theDouble;
		break;
	case 8:
		this->strain4p = info.theDouble;
		break;
	case 9:
		this->stress1n = info.theDouble;
		break;
	case 10:
		this->strain1n = info.theDouble;
		break;
	case 11:
		this->stress2n = info.theDouble;
		break;
	case 12:
		this->strain2n = info.theDouble;
		break;
	case 13:
		this->stress3n = info.theDouble;
		break;
	case 14:
		this->strain3n = info.theDouble;
		break;
	case 15:
		this->stress4n = info.theDouble;
		break;
	case 16:
		this->strain4n = info.theDouble;
		break;
	case 17:
		this->rDispP = info.theDouble;
		break;
	case 18:
		this->rForceP = info.theDouble;
		break;
	case 19:
		this->uForceP = info.theDouble;
		break;
	case 20:
		this->rDispN = info.theDouble;
		break;
	case 21:
		this->rForceN = info.theDouble;
		break;
	case 22:
		this->uForceN = info.theDouble;
		break;
	default:
		return -1;
	}

	// Update the envelope
	this->SetEnvelope();

	return 0;
}