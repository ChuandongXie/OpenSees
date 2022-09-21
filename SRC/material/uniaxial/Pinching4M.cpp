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

void* OPS_Pinching4M()
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

	UniaxialMaterial* mat = 0;
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

		mat = new Pinching4M(tag,
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

		mat = new Pinching4M(tag,
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

	if (mat == 0) {
		opserr << "WARNING: failed to create Pinching4M material\n";
		return 0;
	}
	return mat;
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
	double gf1, double gf2, double gf3, double gf4, double gflim, double ge, int dc):
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

	envlpPosStress.Zero(); envlpPosStrain.Zero(); envlpNegStress.Zero(); envlpNegStrain.Zero();
	energyCapacity = 0.0; kunload = 0.0; elasticStrainEnergy = 0.0;
}