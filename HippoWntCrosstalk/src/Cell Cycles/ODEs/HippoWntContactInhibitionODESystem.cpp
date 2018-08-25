/*
 * HippoWntContactInhibitionODESystem.cpp
 *
 *  Created on: 1 Nov 2017
 *      Author: daniel
 */

#include "HippoWntContactInhibitionODESystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"
#include "IsNan.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "Debug.hpp"

HippoWntContactInhibitionODESystem::HippoWntContactInhibitionODESystem(double wntLevel,
																		boost::shared_ptr<AbstractCellMutationState> pMutationState,
																		std::vector<double> stateVariables)
		: AbstractOdeSystem(16),
		mpMutationState(pMutationState),
		mWntLevel(wntLevel)
		{

		mpSystemInfo.reset(new CellwiseOdeSystemInformation<HippoWntContactInhibitionODESystem>);

		/**
		* State variables are
		*
		*  0. r = pRb
		*  1. e = E2F1 (This is the S-phase indicator)
		*  2. i = CycD (inactive)
		*  3. j = CycD (active)
		*  4. p = pRb-p
		*  5. D = APC destruction complex
		*  6. X = Axin
		*  7. Cu = Beta Cat marked for ubiquitination
		*  8. Cc = Open form Beta Cat
		*  9. A = Free Adhesion molecules
		*  10. Ca = BetaCat/Adhesion
		*  11. T = free TCF
		*  12. Ct = Open BetaCat/TCF
		*  13. Y = Wnt Target protein
		*  14. Wnt level
		*  15. Cy = Beta cat/pYAP
		*/
		Init(); // set up parameter values

		double d_d_hat = mDd + mXiD*wntLevel;
		double d_d_x_hat = mDdx + mXiDx*wntLevel;
		double d_x_hat = mDx + mXiX*wntLevel;
		double p_c_hat = mPc + mXiC*wntLevel;

		double sigma_D = 0.0; // for healthy cells

		if (!mpMutationState)
		{
		// No mutations specified
		}
		else if (mpMutationState->IsType<ApcOneHitCellMutationState>())
		{
			// 0.5
		sigma_D = 0.5;
		}
		else if (mpMutationState->IsType<ApcTwoHitCellMutationState>())
		{
			//1
		sigma_D = 1.0;
		}
		else if (mpMutationState->IsType<BetaCateninOneHitCellMutationState>())
		{
		//sigma_B = 0.5; // never used!
		}
		// Other mutations have no effect.

		// Cell-specific initial conditions
		double steady_D = ((1.0-sigma_D)*mSd*mSx)/((1.0-sigma_D)*mSd*d_d_hat + d_x_hat*(d_d_hat + d_d_x_hat));
		SetDefaultInitialCondition(5, steady_D); // Destruction complex (APC/Axin/GSK3B)

		double temp = (mSx*(d_d_hat+d_d_x_hat))/((1.0-sigma_D)*mSd*d_d_hat+d_x_hat*(d_d_hat+d_d_x_hat));
		SetDefaultInitialCondition(6, temp);  // Axin

		double steady_Cf = ((mSc-mDc*mKd - mPu*steady_D)+sqrt(SmallPow((mSc-mDc*mKd - mPu*steady_D),2) + (4.0*mSc*mDc*mKd)))/(2.0*mDc);
		temp = (mPu*steady_D*steady_Cf)/(mDu*(steady_Cf+mKd));
		SetDefaultInitialCondition(7, temp); // beta-catenin to be ubiquitinated

		double theta = mDc + (mPu*steady_D)/(steady_Cf + mKd);

		double steady_Co = ( mSc - theta*mKc + sqrt(4.0*mSc*theta*mKc + SmallPow((mSc - theta*mKc),2)) )/(2.0*theta);
		SetDefaultInitialCondition(8, steady_Co); // Open form beta-catenin

		SetDefaultInitialCondition(9, mSa/mDa);  // 'Free' adhesion molecules

		SetDefaultInitialCondition(10, mSa*mSca*steady_Co/(mDa*mDca)); // Co-A Adhesion complex

		SetDefaultInitialCondition(11, mSt/mDt); // `Free' transcription molecules (TCF)

		SetDefaultInitialCondition(12, mSct*mSt*steady_Co/(mDt*mDct)); // Co-T open form beta-catenin/TCF

		temp = (mSct*mSt*mSy*steady_Cf)/(mDy*(mSct*mSt*steady_Cf + mDct*mDt*mKt));
		SetDefaultInitialCondition(13, temp); // Wnt target protein

		SetDefaultInitialCondition(14, wntLevel); // Wnt stimulus

		SetDefaultInitialCondition(15, 0.0); // YAP-P/Beta-catenin

		if (stateVariables != std::vector<double>())
		{
		SetStateVariables(stateVariables);
		}
}

void HippoWntContactInhibitionODESystem::SetMutationState(boost::shared_ptr<AbstractCellMutationState> pMutationState)
{
    mpMutationState = pMutationState;
}

HippoWntContactInhibitionODESystem::~HippoWntContactInhibitionODESystem()
{
    // Do nothing
}

void HippoWntContactInhibitionODESystem::Init()
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
	mKc = 200;  //  nM
	mKd = 5;    //  nM
	mKt = 50;   //  nM
	mPc = 0.0;  //  h^-1
	mPu = 100;  //  h^-1
	mXiD = 5;   //  h^-1
	mXiDx = 5;  //  h^-1
	mXiX = 200; //  h^-1
	mKh = 100;
	mPh = 0;

}

void HippoWntContactInhibitionODESystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
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

	// Totals
	double Cf = Cc+Ch;

	double d_d_hat = mDd + mXiD*stimulus_wnt;
	double d_d_x_hat = mDdx + mXiDx*stimulus_wnt;
	double d_x_hat = mDx + mXiX*stimulus_wnt;

	double sigma_D = 0.0;   // for healthy cells
	//double sigma_B = 0.0;   // for healthy cells

	if (!mpMutationState)
	{
		// No mutations specified
	}
	else
	{
		sigma_D = 1.0;
	}

	// Other mutations have no effect.

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


   if (mpMutationState->IsType<ApcTwoHitCellMutationState>()){
	   double mPh_scale = 0;
	   mPh = mPh_scale;
   }
   else {
	   if (mCellVolume<=0.2){
	   double mPh_scale  = 0.0 + ((75)/(1+exp(50*(mCellVolume-((mHealthyCellCIThreshold-0.1)*mEquilibriumVolume/10)))));
	   mPh = mPh_scale;
	   }
	   else{
	   double mPh_scale  = 0.0 + ((75)/(1+exp(50*(mCellVolume-((mHealthyCellCIThreshold-0.1)*mEquilibriumVolume)))));
	   mPh = mPh_scale;
	   }
   }

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

  rDY[15] = ((mPh*Cc)/(Cf+mKh)) - (mSca*Ch*A); // Cy

}

const boost::shared_ptr<AbstractCellMutationState> HippoWntContactInhibitionODESystem::GetMutationState() const
{
    return mpMutationState;
}

bool HippoWntContactInhibitionODESystem::CalculateStoppingEvent(double time, const std::vector<double>& rY)
{
    std::vector<double> dy(rY.size());
    EvaluateYDerivatives(time, rY, dy);

    return (rY[1] > 1.0 && dy[1] > 0.0);
}

double HippoWntContactInhibitionODESystem::CalculateRootFunction(double time, const std::vector<double>& rY)
{
    return rY[1] - 1.0;
}

template<>
void CellwiseOdeSystemInformation<HippoWntContactInhibitionODESystem>::Initialise()
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

double HippoWntContactInhibitionODESystem::GetWntLevel() const
{
    return mWntLevel;
}

void HippoWntContactInhibitionODESystem::SetCellParameters(double cell_volume, double healthy_threshold, double mutant_threshold, double equilibrium_volume)
{
	mCellVolume = cell_volume;
	mHealthyCellCIThreshold = healthy_threshold;
	mMutantCellCIThreshold = mutant_threshold;
	mEquilibriumVolume = equilibrium_volume;
}

double HippoWntContactInhibitionODESystem::GetCellVolume() const
{
	return mCellVolume;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(HippoWntContactInhibitionODESystem)

