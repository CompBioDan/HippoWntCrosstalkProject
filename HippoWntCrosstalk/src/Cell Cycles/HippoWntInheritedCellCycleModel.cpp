/*
 * HippoWntInheritedCellCycleModel.cpp
 *
 *  Created on: 22 Feb 2018
 *      Author: chastebox
 */
/*
 * HippoWntInheritedCellCycleModel.cpp
 *
 *  Created on: 1 Nov 2017
 *      Author: daniel
 */

#include "UblasIncludes.hpp"
#include "HippoWntInheritedCellCycleModel.hpp"
#include "Debug.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "RandomNumberGenerator.hpp"
#include <cxxtest/TestSuite.h>

HippoWntInheritedCellCycleModel::HippoWntInheritedCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
   : AbstractWntOdeBasedCellCycleModel(pOdeSolver),
	 mMutantCellCIThreshold(DOUBLE_UNSET),
	 mHealthyCellCIThreshold(DOUBLE_UNSET),
	 mEquilibriumVolume(DOUBLE_UNSET),
	 mWntThreshold(DOUBLE_UNSET),
	 mWntDivisionValue(DOUBLE_UNSET),
	 mWntNoise(DOUBLE_UNSET),
	 mMotherWnt(DOUBLE_UNSET)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<HippoWntInheritedCellCycleModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        // Chaste solvers always check for stopping events, CVODE needs to be instructed to do so
        mpOdeSolver->CheckForStoppingEvents();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<HippoWntInheritedCellCycleModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.00005);
#endif //CHASTE_CVODE
    }
}


void HippoWntInheritedCellCycleModel::InitialiseOdeSystem(double wntConcentration, boost::shared_ptr<AbstractCellMutationState> pMutationState)
{
    mpOdeSystem = new HippoWntContactInhibitionODESystem(wntConcentration, pMutationState);

}

// Initialise the cell cycle model at the start of a simulation. It is called precisely one per cell setup in the initial population.
void HippoWntInheritedCellCycleModel::Initialise()
{
    assert(mpOdeSystem == NULL);
    assert(mpCell != NULL);

    // Pass this time step's Wnt stimulus into the solver as a constant over this timestep.
    double base_wnt = GetWntLevel();
	double threshold_new = 1-mWntThreshold;
	double wnt_level;

	if (base_wnt < threshold_new){
		wnt_level = 0;
	}
	else{
		wnt_level = 1;
	}

	if (mpCell->GetMutationState()->IsType<ApcTwoHitCellMutationState>()){
		wnt_level = 1;
	}

    InitialiseOdeSystem(wnt_level, mpCell->GetMutationState());

    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());

    ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();

}

// Used to create new copies of this cell cycle model
AbstractCellCycleModel* HippoWntInheritedCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell-cycle model
	HippoWntInheritedCellCycleModel* p_model = new HippoWntInheritedCellCycleModel(mpOdeSolver);

    p_model->SetBirthTime(mBirthTime);
    p_model->SetDimension(mDimension);
    p_model->SetMinimumGapDuration(mMinimumGapDuration);
    p_model->SetStemCellG1Duration(mStemCellG1Duration);
    p_model->SetTransitCellG1Duration(mTransitCellG1Duration);
    p_model->SetSDuration(mSDuration);
    p_model->SetG2Duration(mG2Duration);
    p_model->SetMDuration(mMDuration);
    p_model->SetDivideTime(mDivideTime);
    p_model->SetFinishedRunningOdes(mFinishedRunningOdes);
    p_model->SetG2PhaseStartTime(mG2PhaseStartTime);
    p_model->SetLastTime(mLastTime);
    p_model->SetHealthyCIThreshold(mHealthyCellCIThreshold);
    p_model->SetMutantCIThreshold(mMutantCellCIThreshold);
    p_model->SetEquilibriumVolume(mEquilibriumVolume);
    p_model->SetWntThreshold(mWntThreshold);
    p_model->SetWntNoise(mWntNoise);

    /*
     * Create the new cell-cycle model's ODE system and use the current values
     * of the state variables in mpOdeSystem as an initial condition.
     */

    assert(mpOdeSystem);
    double base_wnt = GetWntLevel();
	double threshold_new = 1-mWntThreshold;
	double wnt_level;

	if (base_wnt < threshold_new){
		wnt_level = 0;
	}
	else{
		wnt_level = 1;
	}

	if (mpCell->GetMutationState()->IsType<ApcTwoHitCellMutationState>()){
		mpOdeSystem->rGetStateVariables()[14] = 1;
	}

    p_model->SetOdeSystem(new HippoWntContactInhibitionODESystem(wnt_level, mpCell->GetMutationState()));
    p_model->SetStateVariables(mpOdeSystem->rGetStateVariables());


    mpOdeSystem->rGetStateVariables()[14] = mMotherWnt*(0.5-mWntDivisionValue);
    RandomNumberGenerator::Instance()->Reseed(time(NULL));

    // Just trying to put the thresholds in again
	double cell_volume = mpCell->GetCellData()->GetItem("volume"); // CURRENT CELL VOUME
	double h_thresh = mHealthyCellCIThreshold; // HEALTHY CELL CI THRESHOLD
	double m_thresh = mMutantCellCIThreshold; // MUTANT CELL CI THRESHOLD

	static_cast<HippoWntContactInhibitionODESystem*>(mpOdeSystem)->SetCellParameters(cell_volume, h_thresh, m_thresh, mEquilibriumVolume); // Set the parameters for the CC

    return p_model;
}


// Reset the cell cycle model and change the wnt concentration
void HippoWntInheritedCellCycleModel::ResetForDivision(){

	AbstractWntOdeBasedCellCycleModel::ResetForDivision();

	// Pass this time step's Wnt stimulus into the solver as a constant over this timestep.
	double base_wnt = GetWntLevel();
	double threshold_new = 1-mWntThreshold;

	RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
	mWntDivisionValue = p_gen->NormalRandomDeviate(0,mWntNoise); // Normally distributed noise term for wnt division bias
	mMotherWnt = mpOdeSystem->rGetStateVariables()[14]; // Keep track of mother wnt level for distirbution
	if (base_wnt < threshold_new){
		mpOdeSystem->rGetStateVariables()[14] = (0.5+mWntDivisionValue)*mpOdeSystem->rGetStateVariables()[14];
	}
	else{
		mpOdeSystem->rGetStateVariables()[14] = 1;
	}

	if (mpCell->GetMutationState()->IsType<ApcTwoHitCellMutationState>()){
		mpOdeSystem->rGetStateVariables()[14] = 1;
	}
	//mpOdeSystem->rGetStateVariables()[14] = wnt_level;

	ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();

}


void HippoWntInheritedCellCycleModel::AdjustOdeParameters(double currentTime){

	if ((mMutantCellCIThreshold == DOUBLE_UNSET) || (mHealthyCellCIThreshold == DOUBLE_UNSET))
	{
			 EXCEPTION("The member variables mQuiescentVolumeFraction and mEquilibriumVolume have not yet been set.");
	}

	HippoWntInheritedCellCycleModel* p_model = static_cast<HippoWntInheritedCellCycleModel*>(mpCell->GetCellCycleModel()); // Static cast to get cellcycle info

	double cell_volume = mpCell->GetCellData()->GetItem("volume"); // CURRENT CELL VOUME

	double h_thresh = mHealthyCellCIThreshold; // HEALTHY CELL CI THRESHOLD
	double m_thresh = mMutantCellCIThreshold; // MUTANT CELL CI THRESHOLD

	static_cast<HippoWntContactInhibitionODESystem*>(mpOdeSystem)->SetCellParameters(cell_volume, h_thresh, m_thresh,mEquilibriumVolume);

	// Use the cell's current mutation status as another input
	static_cast<HippoWntContactInhibitionODESystem*>(mpOdeSystem)->SetMutationState(mpCell->GetMutationState());

	double base_wnt = GetWntLevel();
	double threshold_new = 1-mWntThreshold;
	if (base_wnt >= threshold_new){
		mpOdeSystem->rGetStateVariables()[14] = 1;
	}

	if (mpCell->GetMutationState()->IsType<ApcTwoHitCellMutationState>()){
		mpOdeSystem->rGetStateVariables()[14] = 1;
	}


	ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();
}


void HippoWntInheritedCellCycleModel::SetCIThresholds(double healthy_cells, double mutant_cells){

	mMutantCellCIThreshold = mutant_cells;
	mHealthyCellCIThreshold = healthy_cells;

}

void HippoWntInheritedCellCycleModel::SetHealthyCIThreshold(double healthy_cells){

	mHealthyCellCIThreshold = healthy_cells;

}

void HippoWntInheritedCellCycleModel::SetMutantCIThreshold(double mutant_cells){

	mMutantCellCIThreshold = mutant_cells;

}

double HippoWntInheritedCellCycleModel::GetHealthyCIThreshold(){

	double threshold = mHealthyCellCIThreshold;
	return threshold;

}
double HippoWntInheritedCellCycleModel::GetMutantCIThreshold(){

	double threshold = mMutantCellCIThreshold;
	return threshold;

}
double HippoWntInheritedCellCycleModel::GetWntConcentration(){

	assert(mpOdeSystem != NULL);
	double wnt = mpOdeSystem->rGetStateVariables()[14];
	return wnt;

}

void HippoWntInheritedCellCycleModel::SetEquilibriumVolume(double equilibrium_volume){

	mEquilibriumVolume = equilibrium_volume;

}

void HippoWntInheritedCellCycleModel::SetWntNoise(double wnt_noise){

	mWntNoise = wnt_noise;

}

void HippoWntInheritedCellCycleModel::ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel()
{
    assert(mpOdeSystem != NULL);
    assert(mpCell != NULL);

    double beta_catenin_level =   mpOdeSystem->rGetStateVariables()[12];

    /*
	* This method is usually called within a CellBasedSimulation, after the CellPopulation
	* has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
	* CellPropertyRegistry::Instance() here when setting the CellProliferativeType, we
	* would be creating a new CellPropertyRegistry. In this case the cell proliferative
	* type counts, as returned by AbstractCellPopulation::GetCellProliferativeTypeCount(),
	* would be incorrect. We must therefore access the CellProliferativeType via the cell's
	* CellPropertyCollection.
	*/
    // For mitogenic stimulus of 1/25.0 in Wnt equations

    if (beta_catenin_level < 9.2)
    {
        boost::shared_ptr<AbstractCellProperty> p_diff_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_diff_type);
    }
    else
    {
        boost::shared_ptr<AbstractCellProperty> p_transit_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_transit_type);
    }
}


void HippoWntInheritedCellCycleModel::SetWntThreshold(double wnt_threshold){

	mWntThreshold = wnt_threshold;

}

double HippoWntInheritedCellCycleModel::GetWntThreshold(){

	return mWntThreshold;

}


double HippoWntInheritedCellCycleModel::GetMembraneBoundBetaCateninLevel()
{
    return mpOdeSystem->rGetStateVariables()[10];
}

double HippoWntInheritedCellCycleModel::GetCytoplasmicBetaCateninLevel()
{
    return  mpOdeSystem->rGetStateVariables()[7] + mpOdeSystem->rGetStateVariables()[8]
          + mpOdeSystem->rGetStateVariables()[15];
}

double HippoWntInheritedCellCycleModel::GetNuclearBetaCateninLevel()
{
    return mpOdeSystem->rGetStateVariables()[12];
}

void HippoWntInheritedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output

    // Call method on direct parent class
    AbstractWntOdeBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(HippoWntInheritedCellCycleModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(HippoWntInheritedCellCycleModel)





