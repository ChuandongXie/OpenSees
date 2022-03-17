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
//      Developed  by Chuandong Xie, (chuandongxie@163.com)   (2020)
//-----------------------------------------------------------------------

#ifndef SteelBuckling_h
#define SteelBuckling_h

#include <UniaxialMaterial.h>

class SteelBuckling : public UniaxialMaterial
{
public:
	SteelBuckling(int tag,
		double bar_t, double bar_h, double bar_l, double bar_n,
		double f_y, double f_u, double E_0, double b,
		double R_0, double cR_1, double cR_2, double R_u,
		double alpha, double R_ur);

	SteelBuckling(void);
	virtual ~SteelBuckling();

	const char* getClassType(void) const { return "SteelBuckling"; };

	double getInitialTangent(void);
	UniaxialMaterial* getCopy(void);

	int calTEN(void);
	int calUN(void);
	int calCOM(void);
	int calPB(void);
	int calUR(void);
	int calRU(void);

	double calTENcore(double eps_t, double eps_pl, double eps_t0, double sig_t0, double eps_ty, double sig_ty);
	double calPBcore(double eps_pb, double eps_cy, double sig_cy);
	double calURcore(double eps_ur, double b_ur, double R_ur, double eps_ur1, double eps_ury, double sig_ur1, double sig_ury);

	int setTrialStrain(double strain, double strainRate = 0.0);
	double getStrain(void);
	double getStress(void);
	double getTangent(void);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel &theChannel,
		FEM_ObjectBroker& theBroker);

	void Print(OPS_Stream& s, int flag = 0);

	//int setParameter(const char** argv, int argc, Parameter& param);
	//int updateParameter(int parameterID, Information& info);

protected:

private:
	//  MATERIAL INPUTS -------------------------------------------------------------------

	// steel bar parameter
	double bar_t;		// thickness of bar
	double bar_h;		// width of bar
	double bar_l;		// length of bar
	double bar_n;		// number of bars

	// material parameter
	double f_y;			// yield strength
	double f_u;			// ultimate strength
	double E_0;			// initial stiffness
	double b;			// hardening ratio
	double R_0;			// exp transition elastic-plastic  -- same as Steel02
	double cR_1;		// coefficient for changing R_0 to R  -- same as Steel02
	double cR_2;		// coefficient for changing R_0 to R  -- same as Steel02
	double R_u;			// exp transition kinematic hardening to perfectly plastic -- same as Steel4
	double alpha;		// coefficient for reloading from post-buckling // unloading tangent
	double R_ur;		// exp transition at re-tension
	bool sc_simplifier; // for simplify when the steel square is very short

	// FIXED PROPERTIES -------------------------------------------------------------------
	double eps_y0;		// initial yield strain
	double eps_inc;		// eps increment used for tangent evaluation = 1e-7
	double eps_tr;		// eps determining whether state 5 or state 1 is used
	double eps_l;		// eps at the intersection of the hardening and perfectly plastic asymptotes

	// INTERNAL VARIABLES -----------------------------------------------------------------
	double eps;			// eps at current step
	double sig;			// sig at current step
	double E;			// stiffness (tangent slope) at the current eps-sig point
	double eps_max;		// maximum strain during the total load history
	double eps_pl;		// accumulated plastic strain in the current half-cycle
	double R_y;			// exp transition elastic to kinematic hardening

	double eps_t0;		// eps at previous load reversal point
	double sig_t0;		// sig at previous load reversal point
	double eps_ty;		// eps at intersection of hardening and linear elastic asymptotes
	double sig_ty;		// sig at intersection of hardening and linear elastic asymptotes
	double eps_t;		// eps at tension evelope
	double sig_t;		// sig at tension evelope
	double sig_t_inc;	// sig increment at tension evelope
	double E_t;			// stiffness at tension evelop

	double eps_u0;		// eps at tension-unloading reversal point
	double sig_u0;		// sig at tension-unloading reversal point
	double eps_ul;		// eps at unloading-compression reversal point
	const double sig_ul = 0.0;		// sig at unloading-compression reversal point
	double eps_u;		// eps at unloading evelope
	double sig_u;		// sig at unloading evelope
	double E_u;		// stiffness (tangent slope) at unloading from tensile loading

	double eps_c0;		// eps at unloading-compression reversal point
	double sig_c0;		// sig at unloading-compression reversal point
	double eps_cy;		// eps at compression-buckling reversal point
	double sig_cy;		// sig at compression-buckling reversal point
	double eps_c;		// eps at compression evelope
	double sig_c;		// sig at compression evelope
	double E_sec;		// stiffness (tangent slope) at compressive loading from initial state or unloading

	double eps_pb;		// eps at post-buckling evelope
	double sig_pb;		// sig at post-buckling evelope
	double sig_pb_inc;	// sig increment at post-buckling evelope
	double E_pb;		// stiffness at post-buckling evelope

	double eps_ur;		// eps at unloading-reloading evelope
	double sig_ur;		// sig at unloading-reloading evelope
	double eps_ur0;		// eps at first unloading-reloading evelope
	double sig_ur0;		// sig at first unloading-reloading evelope
	double eps_ur1;		// eps at second unloading-reloading evelope
	double sig_ur1;		// sig at second unloading-reloading evelope
	double eps_ury;
	double sig_ury;
	double b_ur;
	double eps_ur_rat;
	double eps_r1;		// eps at first reverse point 
	double sig_r1;		// sig at first reverse point
	double eps_r2;		// eps at second reverse point
	double sig_r2;		// sig at second reverse point
	double sig_ur_inc;	// sig increment at unloading-reloading evelope
	double E_ur1;		// initial stiffness for post-buckling unloading
	double E_ur2;		// stiffness after first reverse point in post-buckling unloading
	double E_ur0;
	double E_ur;		// stiffness at unloading-reloading evelope

	double eps_ru0;		// eps at unloading-reloading reversal point
	double sig_ru0;		// sig at unloading-reloading reversal point
	double eps_rul;		// eps at unloading-post-buckling reversal point
	double sig_rul;		// sig at unloading-post-buckling reversal point
	double eps_ru;		// eps at unloading from reloading evelope
	double sig_ru;		// sig at unloading from reloading evelope
	double E_ru;		// stiffness at unloading from reloading evelope
	
	// STORED VALUES ----------------------------------------------------------------------
	int state;			// loading state
						/*
						state 0 -- initial state
						state 1 -- tensile loading
						state 2 -- unloading from tensile loading
						state 3 -- compressive loading
						state 4 -- post-buckling loading
						state 5 -- reloading from post-buckling loading
						state 6 -- unloading from reloading
						state 10
						*/

	// VARIABLES FOR LOAD HISTORY MEMORY --------------------------------------------------
	int state_P;
	double eps_P;
	double sig_P;
	double E_P;

	double eps_plP;
	double eps_maxP;

	double eps_t0P;
	double sig_t0P;
	double eps_tyP;
	double sig_tyP;

	double eps_u0P;
	double sig_u0P;
	double eps_ulP;
	double E_uP;

	double eps_c0P;
	double sig_c0P;
	double eps_cyP;
	double sig_cyP;
	double E_secP;

	double eps_ur0P;
	double sig_ur0P;
	double eps_ur1P;
	double sig_ur1P;
	double E_ur1P;
	double E_ur2P;
	double E_ur0P;
	double eps_uryP;
	double sig_uryP;
	double b_urP;

	double eps_ru0P;
	double sig_ru0P;
	double eps_rulP;
	double sig_rulP;
	double E_ruP;

	// TEMPORARY VARIABLES STORED HERE FOR EFFICIENCY -------------------------------------
	double sig_D;		// temp sig
	double delta_eps;	// delta epsilon during the last load step
	double xi;			// xi value for R_y calculation
	double eps_rat;		// eps*
	double eps_bar;		// eps_bar
	double sig_inc;		// sig increment used for tangent evaluation
	double eps_plD;		// temp variable for plastic strain calculation
	double lambda_p;	// slenderness of steel bar
	double temp1, temp2, temp3, temp4;		// temp polynormal coefficients
	double eps_ru_trial;	// trial value
	double sig_ru_trial;	// trial value
	double sig_rul_trial;	// trial value

	
};
#endif
