/*
 * TestHippoWntODE.cpp
 *
 *  Created on: 1 Nov 2017
 *      Author: daniel
 */
/*
 * HippoWntContactInhibitionODESystem.cpp
 *
 *  Created on: 1 Nov 2017
 *      Author: daniel
 */

#include "TestHippoWntODE.hpp"
#include "CellwiseOdeSystemInformation.hpp"
#include "IsNan.hpp"
#include "Debug.hpp"

TestHippoWntODE::TestHippoWntODE(std::vector<double> stateVariables)
		: AbstractOdeSystem(16)
		{

		mpSystemInfo.reset(new CellwiseOdeSystemInformation<TestHippoWntODE>);

		Init();
}

TestHippoWntODE::~TestHippoWntODE()
{
    // Do nothing
}

void TestHippoWntODE::Init()
{
	// Swat (2004) parameters
	double k1 = 1.0;
	double k2 = 1.6;
	double k3 = 0.05;
	double k16 = 0.4;
	double k34 = 0.04;
	double k43 = 0.01;
	double k61 = 0.3;
	double k23 = 0.3;
	double a = 0.04;
	double J11 = 0.5;
	double J12 = 5.0;
	double J61 = 5.0;
	double J62 = 8.0;
	double J13 = 0.002;
	double J63 = 2.0;
	double Km1 = 0.5;
	double Km2 = 4.0;
	double Km4 = 0.3;
	double kp = 0.05;
	double phi_pRb = 0.005;
	double phi_E2F1 = 0.1;
	double phi_CycDi = 0.023;
	double phi_CycDa = 0.03;
	double phi_pRbp = 0.06;

	// Value of the mitogenic factor to make the van Leeuwen model influence cell cycle just the same
	double mitogenic_factorF = 1.0/25.0;

	// Non-dimensionalise parameters
	mk2d = k2/(Km2*phi_E2F1);
	mk3d = k3*mitogenic_factorF/(Km4*phi_E2F1);
	mk34d = k34/phi_E2F1;
	mk43d = k43/phi_E2F1;
	mk23d = k23*Km2/(Km4*phi_E2F1);
	mad = a/Km2;
	mJ11d = J11*phi_E2F1/k1;
	mJ12d = J12*phi_E2F1/k1;
	mJ13d = J13*phi_E2F1/k1;
	mJ61d = J61*phi_E2F1/k1;
	mJ62d = J62*phi_E2F1/k1;
	mJ63d = J63*phi_E2F1/k1;
	mKm1d = Km1/Km2;
	mkpd = kp/(Km2*phi_E2F1);
	mphi_r = phi_pRb/phi_E2F1;
	mphi_i = phi_CycDi/phi_E2F1;
	mphi_j = phi_CycDa/phi_E2F1;
	mphi_p = phi_pRbp/phi_E2F1;
	mk16d = k16*Km4/phi_E2F1;
	mk61d = k61/phi_E2F1;
	mPhiE2F1 = phi_E2F1;

	// Initialize van Leeuwen model parameters
	mSa = 20;   //  nM/h
	mSca = 250; //  (nMh)^-1
	mSc = 25;   //  nM/h
	mSct = 30;  //  (nMh)^-1
	mSd = 100;  //  h^-1
	mSt = 10;   //  nM/h
	mSx = 10;   //  nM/h
	mSy = 10;   //  h^-1
	mDa = 2;    //  h^-1
	mDca = 350; //  h^-1
	mDc = 1;    //  h^-1
	mDct = 750; //  h^-1
	mDd = 5;    //  h^-1
	mDdx = 5;   //  h^-1
	mDt = 0.4;  //  h^-1
	mDu = 50;   //  h^-1
	mDx = 100;  //  h^-1
	mDy = 1;    //  h^-1
	mKd = 5;    //  nM
	mKt = 50;   //  nM
	mPu = 100;  //  h^-1
	mXiD = 5;   //  h^-1
	mXiDx = 5;  //  h^-1
	mXiX = 200; //  h^-1
	mKh = 100;

}

void TestHippoWntODE::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{

	Init();
	double r = rY[0];
	double e = rY[1];
	double i = rY[2];
	double j = rY[3];
	double p = rY[4];

	double dx1 = 0.0;
	double dx2 = 0.0;
	double dx3 = 0.0;
	double dx4 = 0.0;
	double dx5 = 0.0;

	// Bit back-to-front, but work out the Wnt section first...

	// Variables
	double D = rY[5];
	double X = rY[6];
	double Cu = rY[7];
	double Cc = rY[8];
	double A = rY[9];
	double Ca = rY[10];
	double T = rY[11];
	double Ct = rY[12];
	double Y = rY[13];
	double stimulus_wnt = rY[14];
	double Ch = rY[15];
	//PRINT_VARIABLE(stimulus_wnt);

	// Totals
	double Cf = Cc+Ch;

	double d_d_hat = mDd + mXiD*stimulus_wnt;
	double d_d_x_hat = mDdx + mXiDx*stimulus_wnt;
	double d_x_hat = mDx + mXiX*stimulus_wnt;

	double sigma_D = 0.0;   // for healthy cells =0.5 APC1, =1.0 APC2

	// Now the cell cycle stuff...

	// dr
	dx1 = e/(mKm1d+e)*mJ11d/(mJ11d+r)*mJ61d/(mJ61d+p) - mk16d*r*j+mk61d*p-mphi_r*r;
	// de
	dx2 = mkpd+mk2d*(mad*mad+e*e)/(1+e*e)*mJ12d/(mJ12d+r)*mJ62d/(mJ62d+p) - e;
	// di - changed to include Ct+Mt - transcriptional beta-catenin
	dx3 = mk3d*(Ct) + mk23d*e*mJ13d/(mJ13d+r)*mJ63d/(mJ63d+p) + mk43d*j - mk34d*i*j/(1+j) - mphi_i*i;
	// dj
	dx4 = mk34d*i*j/(1+j) - (mk43d+mphi_j)*j;
	// dp
	dx5 = mk16d*r*j - mk61d*p - mphi_p*p;

	double factor = mPhiE2F1*60.0;  // Convert non-dimensional d/dt s to d/dt in hours.

    rDY[0] = dx1*factor;
    rDY[1] = dx2*factor;
    rDY[2] = dx3*factor;
    rDY[3] = dx4*factor;
    rDY[4] = dx5*factor;


    double mean_vol = 0.68;

    mHealthyCellCIThreshold = 1;
    mCellVolume = 1;

    double mPh_scale  = 0.0 + ((80)/(1+exp(50*(mCellVolume-(mHealthyCellCIThreshold*mean_vol)))));
    mPh = mPh_scale;
    //PRINT_VARIABLE(mPh);

       // The adapted van Leeuwen and Ward ODE system
	  rDY[5] = (1.0-sigma_D)*mSd*X - (d_d_hat + d_d_x_hat)*D; // D

	  rDY[6] = mSx - (1.0-sigma_D)*mSd*X - d_x_hat*X + d_d_x_hat*D; // X

	  rDY[7] = (mPu*D*Cf)/(Cf+mKd) - mDu*Cu; // Cu

	  rDY[8] = mSc + mDca*Ca + mDct*Ct - (mSca*A + mSct*T + mDc)*Cc
			 - (mPu*D*Cc)/(Cf+mKd) - ((mPh*Cc)/(Cf+mKh)); // Cc

	  rDY[9] = mSa + mDca*(Ca) - ((mSca*Ch) + mDa + (mSca*Cc))*A; // A

	  rDY[10] = (mSca*Cc*A) - (mDca*Ca) + (mSca*Ch*A); // Ca

	  rDY[11] = mSt + (mDct*(Ct)) - (mSct*(Cf)*T) - (mDt*T); // T

	  rDY[12] = (mSct*Cc*T) - (mDct*Ct); // Ct

	  rDY[13] = ((mSy*(Ct))/(Ct + mKt)) - (mDy*Y); // Y

	  rDY[14] = 0.0;  // don't interfere with Wnt stimulus

	  rDY[15] = ((mPh*Cc)/(Cf+mKh)) - (mSca*Ch*A) - ((mPu*D*Ch)/(Cf+mKd)); // Cy

}


template<>
void CellwiseOdeSystemInformation<TestHippoWntODE>::Initialise()
{
    this->mVariableNames.push_back("pRb");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(7.357000000000000e-01);

    this->mVariableNames.push_back("E2F1");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(1.713000000000000e-01);

    this->mVariableNames.push_back("CycD_i");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(6.900000000000001e-02);

    this->mVariableNames.push_back("CycD_a");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(3.333333333333334e-03);

    this->mVariableNames.push_back("pRb_p");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(1.000000000000000e-04);

    this->mVariableNames.push_back("D");  // Destruction complex (APC/Axin/GSK3B)
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("X");  // Axin
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Cu"); // beta-catenin to be ubiquitinated
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Cc"); // Open form beta-catenin
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("A");  // 'Free' adhesion molecules
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Ca"); // Co-A Adhesion complex
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled indouble  later

    this->mVariableNames.push_back("T"); // `Free' transcription molecules (TCF)
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Ct"); // Co-T open form beta-catenin/TCF
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Y"); // Wnt target protein
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Sw"); // Wnt stimulus
    this->mVariableUnits.push_back("nM");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mVariableNames.push_back("Cy"); // Wnt stimulus
	this->mVariableUnits.push_back("nM");
	this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(TestHippoWntODE)



