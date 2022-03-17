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

// $Revision:
// $Date:
// $Source:

// Written: Chuandong Xie at Xi'an University of Architecture & Technology
// Version: 0.1 beta
// Created: 2020/12/01 20:58:02
// 
//
// Description: This file contains the class implementation of SteelBuckling.
// SteelBuckling provides the abstraction of an steel squares under cyclic loading
// considering buckling.
// Notes: 1. The material is for zeroLength element only.
//        2. This is a beta version. Further refine work is needed. 
//-----------------------------------------------------------------------
//              Steel bar considering buckling
//      Developed  by Chuandong Xie, (chuandongxie@xauat.edu.com)   (2020)
//-----------------------------------------------------------------------

#include <math.h>
#include <stdlib.h>
#include <SteelBuckling.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <elementAPI.h>
#include <OPS_Globals.h>

#ifdef _USRDLL
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MAXOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

static int numSteelBuckling = 0;

void*
OPS_SteelBuckling()
{
	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial* theMaterial = 0;

	int iData[1];
	double dData[14];

	int argc = 1,
		numData = 14;

	if (numSteelBuckling == 0)
	{
		opserr << "SteelBuckling uniaxial material - Developed by Chuandong Xie\n";
		numSteelBuckling = 1;
	}

	if (OPS_GetIntInput(&argc, iData) != 0) {
		opserr << "WARNING invalid uniaxialMaterial SteelBuckling tag\n" << endln;
		return 0;
	}

	// Read the basic parameters
	argc = OPS_GetNumRemainingInputArgs();

	if (argc != numData) {
		opserr << "Invalid $argc, want: uniaxialMaterial SteelBuckling " << iData[0] << " bar_t? bar_h? bar_l? bar_n? f_y? E_0? b? f_u? R_0? cR_1? cR_2? R_u? alpha? R_ur?" << endln;
		return 0;
	}

	if (OPS_GetDoubleInput(&argc, dData) != 0) {
		opserr << "Invalid #args, want: uniaxialMaterial SteelBuckling " << iData[0] << " bar_t? bar_h? bar_l? bar_n? f_y? E_0? b? f_u? R_0? cR_1? cR_2? R_u? alpha? R_ur?" << endln;
		return 0;
	}

	// Parsing was successful, allocate the material
	theMaterial = new SteelBuckling(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4],
		dData[5], dData[6], dData[7], dData[8], dData[9], dData[10], dData[11],
		dData[12], dData[13]);

	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type SteelBuckling Material\n";
		return 0;
	}

return theMaterial;
}


SteelBuckling::SteelBuckling(int tag,
	double _bar_t, double _bar_h, double _bar_l, double _bar_n,
	double _f_y, double _f_u, double _E_0, double _b,
	double _R_0, double _cR_1, double _cR_2, double _R_u,
	double _alpha, double _R_ur) :
	UniaxialMaterial(tag, MAT_TAG_SteelBuckling),
	bar_t(_bar_t), bar_h(_bar_h), bar_l(_bar_l), bar_n(_bar_n),
	f_y(_f_y), f_u(_f_u), E_0(_E_0), b(_b),
	R_0(_R_0), cR_1(_cR_1), cR_2(_cR_2), R_u(_R_u),
	alpha(_alpha), R_ur(_R_ur)
{
	// initial parameters
	state_P = 0.0;
	eps = 0.0;
	sig = 0.0;

	eps_y0 = f_y / E_0;
	eps_l = (f_u - f_y) / (b * E_0) + eps_y0;
		
	eps_tr = 9.0 * eps_y0;
	eps_pl = 0.0;
	eps_plP = 0.0;
	eps_maxP = eps_y0;
	eps_inc = 0.00000000001;

	E_uP = E_0;
	// lambda_p = sqrt((f_y * bar_l) / (100.0 * bar_h)); // for zerolength element
	lambda_p = sqrt((f_y / bar_t / bar_h / bar_n) / (100.0)) * bar_l / bar_h;

	E_t = 0.0;
	E_sec = 0.0;
	E_pb = 0.0;
	E_ur = 0.0;
	E_ur1 = 0.0;
	E_ur2 = 0.0;
	E_ru = 0.0;

	//
	//this->revertToStart();
	//this->revertToLastCommit();
	bool sc_simplifier = false; // for simplify when the steel square is very short
}


SteelBuckling::SteelBuckling(void):
	UniaxialMaterial(0, MAT_TAG_SteelBuckling)
{
	state = 0;
}

SteelBuckling::~SteelBuckling(void)
{
	// Does nothing
}

UniaxialMaterial*
SteelBuckling::getCopy(void)
{
	SteelBuckling* theCopy = new SteelBuckling(this->getTag(), 
		bar_t, bar_h, bar_l, bar_n,
		f_y, f_u, E_0, b, R_0, cR_1, cR_2, R_u, alpha, R_ur);

	// current round;
	theCopy->state = state;
	theCopy->eps = eps;
	theCopy->sig = sig;
	theCopy->E = E;
	theCopy->eps_pl = eps_pl;
	theCopy->eps_max = eps_max;

	theCopy->eps_t0 = eps_t0;
	theCopy->sig_t0 = sig_t0;
	theCopy->eps_ty = eps_ty;
	theCopy->sig_ty = sig_ty;

	theCopy->eps_u0 = eps_u0;
	theCopy->sig_u0 = sig_u0;
	theCopy->eps_ul = eps_ul;
	theCopy->E_u = E_u;

	theCopy->eps_c0 = eps_c0;
	theCopy->sig_c0 = sig_c0;
	theCopy->eps_cy = eps_cy;
	theCopy->sig_cy = sig_cy;
	theCopy->E_sec = E_sec;

	theCopy->eps_ur0 = eps_ur0;
	theCopy->sig_ur0 = sig_ur0;
	theCopy->eps_ur1 = eps_ur1;
	theCopy->sig_ur1 = sig_ur1;

	// previous round
	theCopy->state_P = state_P;
	theCopy->eps_P = eps_P;
	theCopy->sig_P = sig_P;
	theCopy->E_P = E_P;
	theCopy->eps_plP = eps_plP;
	theCopy->eps_maxP = eps_maxP;

	theCopy->eps_t0P = eps_t0P;
	theCopy->sig_t0P = sig_t0P;
	theCopy->eps_tyP = eps_tyP;
	theCopy->sig_tyP = sig_tyP;

	theCopy->eps_u0P = eps_u0P;
	theCopy->sig_u0P = sig_u0P;
	theCopy->eps_ulP = eps_ulP;
	theCopy->E_uP = E_uP;

	theCopy->eps_c0P = eps_c0P;
	theCopy->sig_c0P = sig_c0P;
	theCopy->eps_cyP = eps_cyP;
	theCopy->sig_cyP = sig_cyP;
	theCopy->E_secP = E_secP;

	theCopy->eps_ur0P = eps_ur0P;
	theCopy->sig_ur0P = sig_ur0P;
	theCopy->eps_ur1P = eps_ur1P;
	theCopy->sig_ur1P = sig_ur1P;

	return theCopy;
}

double
SteelBuckling::getInitialTangent(void)
{
	return E_0;
}

int
SteelBuckling::setTrialStrain(double trialStrain, double strainRate)
{
	eps = trialStrain;
	delta_eps = eps - eps_P;		// delta epsilon
	eps_max = eps_maxP;
	state = state_P;

	if (delta_eps > 0.0) {			// tensile loading
		if (state == 0 || state == 10) {				// initial tensile loading
			if (fabs(delta_eps) < 10.0 * DBL_EPSILON) {
				E = E_0;
				sig = 0.0;
				state = 10;
				return 0;
			}
			else {
				calTEN();
			}
		}
		else if (state == 1) {		// continued tensile loading
			calTEN();
		}
		else if (state == 2) {		// reloading from unloading
			calTEN();
		}
		else if (state == 3) {		// reloading from compression
			calTEN();
		}
		else if (state == 4) {		// reloading from post-buckling
			if (fabs(eps - eps_max) < eps_tr) {	// tension from post-buckling
				calTEN();
			}
			else {					// unloading-reloading from post-buckling
				calUR();
			}
		}
		else if (state == 5) {		// continue unloading and reloading
			calUR();
		}
		else if (state == 6) {		// loading from unloading
			calRU();
		}
		else {
			opserr << "Invalid State for tension in setTrialStrain " << endln;
			return -1;
		}
	}
	else {
		if (state == 0 || state == 10) {		// initial compression
			if (fabs(delta_eps) < 10.0 * DBL_EPSILON) {
				E = E_0;
				sig = 0.0;
				state = 10;
				return 0;
			}
			else {
				calCOM();
			}
		}
		else if (state == 1) {		// unloading from tensile loading
			calUN();
		}
		else if (state == 2) {		// unloading
			calUN();
		}
		else if (state == 3) {		// continue compression
			calCOM();
		}
		else if (state == 4) {		// continue post-buckling compression
			calPB();
		}
		else if (state == 5) {		// unloading from reloading
			calRU();
		}
		else if (state == 6) {		// continue unloading from reloading
			calRU();
		}
		else {
			opserr << "Invalid State for tension in setTrialStrain " << endln;
			return -1;
		}
	}
	//opserr << "(eps,sig) = (" << eps << ", " << sig << ")" << endln;

	return 0;
}

int
SteelBuckling::calTEN(void)
{
	eps_pl = eps_plP;			// update plastic strain
	if (state == 0) {				// initial tensile loading
		eps_t0 = 0.0;
		sig_t0 = 0.0;
		eps_ty = eps_y0;
		sig_ty = f_y;
	}
	else if (state == 1) {			// continued tensile loading
		eps_t0 = eps_t0P;
		sig_t0 = sig_t0P;
		eps_ty = eps_tyP;
		sig_ty = sig_tyP;
	}
	else if (state == 2 || state == 3) {
		eps_t0 = eps_P;
		sig_t0 = sig_P;
		eps_ty = (f_y * (1 - b) - sig_t0 + E_0 * eps_t0) / (E_0 * (1 - b));
		sig_ty = f_y * (1 - b) + eps_ty * E_0 * b;
	}
	else if (state == 4) {
		eps_t0 = eps_P;
		sig_t0 = sig_P;
		eps_ty = (f_y * (1 - b) - sig_t0 + E_0 * eps_t0) / (E_0 * (1 - b));
		sig_ty = f_y * (1 - b) + eps_ty * E_0 * b;
	}
	else if (state == 5) {			// tensile loading from unloading and reloading from post-buckling
		eps_t0 = eps_P;				// set eps_t0, sig_t0 in last unloading-reloading state
		sig_t0 = sig_P;
		eps_ty = eps_tyP;
		sig_ty = sig_tyP;
	}
	else {
		opserr << "Invalid State for calTEN " << endln;
		return -1;
	}
	// Menegotto-Pinto steel model
	eps_t = eps;
	sig_t = calTENcore(eps_t, eps_pl, eps_t0, sig_t0, eps_ty, sig_ty);
	sig_t_inc = calTENcore(eps_t + eps_inc, eps_pl, eps_t0, sig_t0, eps_ty, sig_ty) - sig_t;
	E_t = sig_t_inc / eps_inc;

	// update values
	eps_max = std::max(eps_t, eps_max);
	sig = sig_t;
	E = E_t;
	state = 1;

	return 0;
}

int
SteelBuckling::calUN(void)
{
	if (state == 1) {					// unloading from tension loading
		eps_pl = std::max(eps_max - eps_y0, 0.0);
		E_u = E_0 * (0.82 + 1.0 / (5.55 + 1000 * eps_pl));
		eps_u0 = eps_P;			
		sig_u0 = sig_P;
		eps_ul = (sig_ul - sig_u0) / E_u + eps_u0;
	}
	else if (state == 2) {				// continue unloading
		eps_pl = eps_plP;
		E_u = E_uP;
		eps_u0 = eps_u0P;
		sig_u0 = sig_u0P;
		eps_ul = eps_ulP;
	}
	else {
		opserr << "Invalid State for calUN " << endln;
		return -1;
	}
	eps_u = eps;
	sig_u = (eps_u - eps_u0) * E_u + sig_u0;

	// check if sig_u less than sig_ul
	if (sig_u < sig_ul) {
		eps_plP = eps_pl;
		calCOM();
	}
	else {
		// update values
		sig = sig_u;
		E = E_u;
		state = 2;
	}

	return 0;
}

int
SteelBuckling::calCOM(void)
{
	eps_pl = eps_plP;
	if (state == 0) {				// initial compression
		eps_c0 = 0.0;
		sig_c0 = 0.0;
		E_sec = E_0 / (1.0 + pow(eps_pl / 0.005, 0.5));
		E_sec = std::max(E_sec, 0.2 * E_0);
		sig_cy = -f_y / (1 + 0.7 * (eps_pl / bar_l) * lambda_p);
		eps_cy = eps_c0 - (sig_c0 - sig_cy) / E_sec;
	}
	else if (state == 1 || state == 2) {			// compression from tension or unloading
		eps_c0 = eps_ul;
		sig_c0 = sig_ul;
		E_sec = E_0 / (1.0 + pow(eps_pl / 0.005, 0.5));
		E_sec = std::max(E_sec, 0.2 * E_0);
		sig_cy = -f_y / (1 + 0.7 * (eps_pl / bar_l) * lambda_p);
		eps_cy = eps_c0 - (sig_c0 - sig_cy) / E_sec;
	}
	else if (state == 3) {			// compression continued
		eps_c0 = eps_c0P;
		sig_c0 = sig_c0P;
		E_sec = E_secP;
		sig_cy = sig_cyP;
		eps_cy = eps_cyP;
	}
	else if (state == 5) {			// compression (unloading) from unloading-reloading state
		eps_c0 = eps_P;
		sig_c0 = sig_P;
		E_sec = E_secP;
		eps_cy = eps_cyP;
		sig_cy = sig_cyP;
	}
	else {
		opserr << "Invalid State for calUN " << endln;
		return -1;
	}
	
	// calculate stress
	eps_c = eps;
	sig_c = (eps_c - eps_c0) * E_sec + sig_c0;

	if (sig_c < sig_cy) {
		eps_cyP = eps_cy;
		sig_cyP = sig_cy;
		calPB();
	}
	else {
		sig = sig_c;
		E = E_sec;
		state = 3;
	}

	return 0;
}

int
SteelBuckling::calPB(void)
{
	eps_cy = eps_cyP;
	sig_cy = sig_cyP;
	eps_pb = eps;
	sig_pb = calPBcore(eps_pb, eps_cy, sig_cy);
	sig_pb_inc = sig_pb - calPBcore(eps_pb - eps_inc, eps_cy, sig_cy);
	E_pb = sig_pb_inc / eps_inc;

	// update values
	sig = sig_pb;
	E = E_pb;
	state = 4;

	return 0;
}

int
SteelBuckling::calUR(void)
{
	if (state == 4) {
		eps_ur0 = eps_P; // reversal point
		sig_ur0 = sig_P;

		if (eps_u0P >= eps_y0) { // yielding point
			eps_ty = eps_u0P;
			sig_ty = sig_u0P;
		}
		else {
			eps_ty = eps_y0;
			sig_ty = f_y;
		}

		E_u = E_uP;
		E_ur0 = alpha * E_u;
		
		sig_ur1 = 0.15 * sig_ty; // get ur1 point
		eps_ur1 = (sig_ur1 - sig_ur0) / E_ur0 + eps_ur0;
		E_ur1 = 20.0 / lambda_p * (sig_ty - sig_ur1) / (eps_ty - eps_ur1);
		E_ur2 = lambda_p / 20.0 * (sig_ty - sig_ur1) / (eps_ty - eps_ur1); // tangent of piecewise part
		
		eps_ury = (-E_ur1 * eps_ur1 + E_ur2 * eps_ty + sig_ur1 - sig_ty) / (E_ur2 - E_ur1);
		sig_ury = (-E_ur1 * E_ur2 * eps_ur1 + E_ur1 * E_ur2 * eps_ty + E_ur2 * sig_ur1 - E_ur1 * sig_ty) / (E_ur2 - E_ur1);
		b_ur = E_ur2 / E_ur1;
	}
	else if (state == 5 || state == 6) {
		eps_ur0 = eps_ur0P;
		sig_ur0 = sig_ur0P;
		eps_ur1 = eps_ur1P;
		sig_ur1 = sig_ur1P;
		eps_ty = eps_tyP;
		sig_ty = sig_tyP;
		E_ur0 = E_ur0P;
		E_ur1 = E_ur1P;
		E_ur2 = E_ur2P;
		eps_ury = eps_uryP;
		sig_ury = sig_uryP;
		b_ur = b_urP;
	}
	else {
		opserr << "Invalid state in calUR()\n";
		return -1;
	}

	eps_ur = eps; // obtain current strain

	if (eps_ur < eps_ur1) { // at ur1 state, a linear relationship will be used
		sig_ur = E_ur0 * (eps_ur - eps_ur0) + sig_ur0;
		// update
		sig = sig_ur;
		E = E_ur0;
		state = 5;
		
	}
	else {
		if (sc_simplifier == false) { // for simplify when the steel square is very short
			sig_ur = calURcore(eps_ur, b_ur, R_ur, eps_ur1, eps_ury, sig_ur1, sig_ury);
			sig_ur_inc = calURcore(eps_ur + eps_inc, b_ur, R_ur, eps_ur1, eps_ury, sig_ur1, sig_ury);

			if (sig_ur < 0.85 * sig_ty) {
				sig = sig_ur;
				E = sig_ur_inc / eps_inc;
				state = 5;
			}
			else {
				eps_tyP = eps_ty;
				sig_tyP = sig_ty;
				state = 5;
				calTEN();
			}
		}
		else {
			eps_tyP = eps_ty;
			sig_tyP = sig_ty;
			state = 5;
			calTEN();
		}
	}
	
	return 0;
}

int
SteelBuckling::calRU(void)
{
	if (state == 5) {
		eps_ru0 = eps_P;
		sig_ru0 = sig_P;
		E_ru = E_secP;

		eps_cy = eps_cyP;
		sig_cy = sig_cyP;
		eps_ru_trial = eps;

		do {		// obtain eps_rul, sig_rul
			sig_ru_trial = (eps_ru_trial - eps_ru0) * E_ru + sig_ru0;
			sig_rul_trial = calPBcore(eps_ru_trial, eps_cy, sig_cy);
			eps_ru_trial -= -0.00001;
		} while (sig_ru_trial < sig_rul_trial);

		eps_rul = eps_ru_trial;
		sig_rul = sig_rul_trial;
	}
	else {
		eps_ru0 = eps_ru0P;
		sig_ru0 = sig_ru0P;
		eps_rul = eps_rulP;
		sig_rul = sig_rulP;
		E_ru = E_ruP;
	}

	eps_ru = eps;
	sig_ru = (eps_ru - eps_ru0) * E_ru + sig_ru0;

	if (sig_ru < sig_rul) {
		calPB();
	}
	else if (sig_ru > sig_ru0) {
		calUR();
	}
	else {
		// update values
		sig = sig_ru;
		E = E_ru;
		state = 6;
	}

	return 0;
}

double 
SteelBuckling::calTENcore(double eps_t, double eps_pl, double eps_t0, double sig_t0, double eps_ty, double sig_ty)
{
	xi = eps_pl / f_y;
	R_y = R_0 * (1 - (cR_1 * xi) / (cR_2 + xi));
	eps_rat = (eps_t - eps_t0) / (eps_ty - eps_t0);
	eps_bar = (eps_t - eps_t0) / (eps_l - eps_t0);
	// calculate the stress corresponding to eps_t
	sig_D = eps_rat * (b / (pow(2.0 + pow(eps_bar, R_u), 1 / R_u)) + (1 - b) / (pow((1.0 + pow(eps_rat, R_y)), 1 / R_y))) * (sig_ty - sig_t0) + sig_t0;

	return sig_D;
}

double 
SteelBuckling::calPBcore(double eps_pb, double eps_cy, double sig_cy)
{
	if (sc_simplifier == false) {
		sig_D = (bar_h * sig_cy) / (pow((2 * (eps_cy - eps_pb) * bar_l * bar_l), 1.0 / 2.0) + bar_h);
		//sig_D = (1.0 / (pow(2.0, .5))) * (bar_h / bar_l) * (1.0 / pow(eps_pb, .5)) * sig_cy;
		// TODO
	}
	else {
		sig_D = sig_cy;
	}

	return sig_D;
}

double
SteelBuckling::calURcore(double eps_ur, double b_ur, double R_ur, double eps_ur1, double eps_ury, double sig_ur1, double sig_ury) {
	eps_ur_rat = (eps_ur - eps_ur1) / (eps_ury - eps_ur1);
	sig_D = (b_ur * eps_ur_rat + ((1.0 - b_ur) * eps_ur_rat) / (pow(1.0 + pow(eps_ur_rat, R_ur), 1.0 / R_ur))) * (sig_ury - sig_ur1) + sig_ur1;
	return sig_D;
}

double
SteelBuckling::getStrain(void)
{
	return eps;
}

double
SteelBuckling::getStress(void)
{
	return sig;
}

double
SteelBuckling::getTangent(void)
{
	return E;
}

int
SteelBuckling::commitState(void)
{
	state_P = state;
	eps_P = eps;
	sig_P = sig;
	E_P = E;

	eps_plP = eps_pl;
	eps_maxP = eps_max;

	eps_t0P = eps_t0;
	sig_t0P = sig_t0;
	eps_tyP = eps_ty;
	sig_tyP = sig_ty;

	eps_u0P = eps_u0;
	sig_u0P = sig_u0;
	eps_ulP = eps_ul;
	E_uP = E_u;

	eps_c0P = eps_c0;
	sig_c0P = sig_c0;
	eps_cyP = eps_cy;
	sig_cyP = sig_cy;
	E_secP = E_sec;

	eps_ur0P = eps_ur0;
	sig_ur0P = sig_ur0;
	eps_ur1P = eps_ur1;
	sig_ur1P = sig_ur1;

	eps_uryP = eps_ury;
	sig_uryP = sig_ury;
	b_urP = b_ur;

	eps_ru0P = eps_ru0;
	sig_ru0P = sig_ru0;
	eps_rulP = eps_rul;
	sig_rulP = sig_rul;
	E_ruP = E_ru;

	E_ur1P = E_ur1; // add 2021-05-25
	E_ur2P = E_ur2;
	E_ur0P = E_ur0;

	return 0;
}

int
SteelBuckling::revertToLastCommit(void)
{
	state = state_P;
	eps = eps_P;
	sig = sig_P;
	E = E_P;

	eps_pl = eps_plP;
	eps_max = eps_maxP;
	
	eps_t0 = eps_t0P;
	sig_t0 = sig_t0P;
	eps_ty = eps_tyP;
	sig_ty = sig_tyP;

	eps_u0 = eps_u0P;
	sig_u0 = sig_u0P;
	eps_ul = eps_ulP;
	E_u = E_uP;

	eps_c0 = eps_c0P;
	sig_c0 = sig_c0P;
	eps_cy = eps_cyP;
	sig_cy = sig_cyP;
	E_sec = E_secP;

	eps_ur0 = eps_ur0P;
	sig_ur0 = sig_ur0P;
	eps_ur1 = eps_ur1P;
	sig_ur1 = sig_ur1P;

	eps_ury = eps_uryP;
	sig_ury = sig_uryP;
	b_ur = b_urP;

	eps_ru0 = eps_ru0P;
	sig_ru0 = sig_ru0P;
	eps_rul = eps_rulP;
	sig_rul = sig_rulP;
	E_ru = E_ruP;

	E_ur1 = E_ur1P; // add 2021-05-25
	E_ur2 = E_ur2P;
	E_ur0 = E_ur0P;

	return 0;
}

int
SteelBuckling::revertToStart(void)
{
	state_P = 0;
	eps_P = 0.0;
	sig_P = 0.0;
	eps = 0.0;
	sig = 0.0;
	E_P = E_0;

	eps_plP = 0.0;
	eps_maxP = f_y / E_0;

	eps_t0P = 0.0;
	sig_t0P = 0.0;
	eps_tyP = f_y / E_0;
	sig_tyP = f_y;

	eps_u0P = 0.0;
	sig_u0P = 0.0;
	eps_ul = 0.0;
	E_uP = E_0;

	eps_c0P = eps_c0;
	sig_c0P = sig_c0;
	eps_cyP = -f_y / E_0;
	sig_cyP = -f_y;
	E_secP = E_0;

	eps_ur0P = 0.0;
	sig_ur0P = 0.0;
	eps_ur1P = 0.0;
	sig_ur1P = 0.0;
	eps_uryP = 0.0;
	sig_uryP = 0.0;
	b_urP = 0.0;

	eps_ru0P = 0.0;
	sig_ru0P = 0.0;
	eps_rulP = 0.0;
	sig_rulP = 0.0;
	E_ruP = 0.0;

	E_ur1P = 0.0;  // add 2021-05-25
	E_ur2P = 0.0;
	E_ur0P = 0.0;

	return 0;
}

int
SteelBuckling::sendSelf(int commitTag, Channel& theChannel)
{
	static Vector data(49);
	data(0) = this->getTag();

	// material properties
	data(1) = bar_t;
	data(2) = bar_h;
	data(3) = bar_l;
	data(4) = bar_n;
	data(5) = f_y;
	data(6) = f_u;
	data(7) = E_0;
	data(8) = b;
	data(9) = R_0;
	data(10) = cR_1;
	data(11) = cR_2;
	data(12) = R_u;
	data(13) = alpha;
	data(14) = R_ur;

	// history variables from last converged state
	data(15) = state_P;
	data(16) = eps_P;
	data(17) = sig_P;
	data(18) = E_P;
	
	data(19) = eps_plP;
	data(20) = eps_maxP;

	data(21) = eps_t0P;
	data(22) = sig_t0P;
	data(23) = eps_tyP;
	data(24) = sig_tyP;

	data(25) = eps_u0P;
	data(26) = sig_u0P;
	data(27) = eps_ulP;
	data(28) = E_uP;

	data(29) = eps_c0P;
	data(30) = sig_c0P;
	data(31) = eps_cyP;
	data(32) = sig_cyP;
	data(33) = E_secP;

	data(34) = eps_ur0P;
	data(35) = sig_ur0P;
	data(36) = eps_ur1P;
	data(37) = sig_ur1P;
	data(38) = eps_uryP;
	data(39) = sig_uryP;
	data(40) = b_urP;
	
	data(41) = eps_ru0P;
	data(42) = sig_ru0P;
	data(43) = eps_rulP;
	data(44) = sig_rulP;
	data(45) = E_ruP;

	data(46) = E_ur1P; // add 2021-05-25
	data(47) = E_ur2P;
	data(48) = E_ur0P;

	if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "SteelBuckling::sendSelf() - failed to sendSelf\n";
		return -1;
	}

	return 0;
}

int
SteelBuckling::recvSelf(int commitTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	static Vector data(49);

	if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "SteelBuckling::recvSelf() - failed to recvSelf\n";
		return -1;
	}

	this->setTag(int(data(0)));

	// material properties
	bar_t = data(1);
	bar_h = data(2);
	bar_l = data(3);
	bar_n = data(4);
	f_y = data(5);
	f_u = data(6);
	E_0 = data(7);
	b = data(8);
	R_0 = data(9);
	cR_1 = data(10);
	cR_2 = data(11);
	R_u = data(12);
	alpha = data(13);
	R_ur = data(14);

	// history variables from last converged state
	state_P = data(15);
	eps_P = data(16);
	sig_P = data(17);
	E_P = data(18);

	eps_plP = data(19);
	eps_maxP = data(20);

	eps_t0P = data(21);
	sig_t0P = data(22);
	eps_tyP = data(23);
	sig_tyP = data(24);

	eps_u0P = data(25);
	sig_u0P = data(26);
	eps_ulP = data(27);
	E_uP = data(28);

	eps_c0P = data(29);
	sig_c0P = data(30);
	eps_cyP = data(31);
	sig_cyP = data(32);
	E_secP = data(33);

	eps_ur0P = data(34);
	sig_ur0P = data(35);
	eps_ur1P = data(36);
	sig_ur1P = data(37);
	eps_uryP = data(38);
	sig_uryP = data(39);
	b_urP = data(40);
	

	eps_ru0P = data(41);
	sig_ru0P = data(42);
	eps_rulP = data(43);
	sig_rulP = data(44);
	E_ruP = data(45);

	E_ur1P = data(46); // add 2021 - 05 - 25
	E_ur2P = data(47);
	E_ur0P = data(48);

	E = E_P;
	eps = eps_P;
	sig = sig_P;

	return 0;
}


void
SteelBuckling::Print(OPS_Stream& s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "SteelBuckling tag: " << this->getTag() << endln;
		s << " bar_t : " << bar_t << ", ";
		s << " bar_h : " << bar_h << ", ";
		s << " bar_l : " << bar_l << ", ";
		s << " bar_n : " << bar_n << ", ";
		s << " f_y : " << f_y << ", ";
		s << " f_u : " << f_u << ", ";
		s << " E_0 : " << E_0 << ", ";
		s << " b : " << b << ", ";
		s << " R_0 : " << ", ";
		s << " cR_1 : " << cR_1 << ", ";
		s << " cR_2 : " << cR_2 << ", ";
		s << " R_u : " << R_u << ", ";
		s << " alpha : " << alpha << ", ";
		s << " R_ur : " << R_ur << ", ";
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"SteelBuckling\", ";
		s << "\"bar_t\": " << bar_t << ", ";
		s << "\"bar_h\": " << bar_h << ", ";
		s << "\"bar_l\": " << bar_l << ", ";
		s << "\"bar_n\": " << bar_n << ", ";
		s << "\"f_y\": " << f_y << ", ";
		s << "\"f_u\": " << f_u << ", ";
		s << "\"E_0\": " << E_0 << ", ";
		s << "\"b\": " << b << ", ";
		s << "\"R_0\": " << R_0 << ", ";
		s << "\"cR_1\": " << cR_1 << ", ";
		s << "\"cR_2\": " << cR_2 << ", ";
		s << "\"R_u\": " << R_u << ", ";
		s << "\"alpha\": " << alpha << ", ";
		s << "\"R_ur\": " << R_ur << ", ";
	}
}

// not include sensitivity