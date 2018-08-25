/*
 * HippoDynamicWntCellCycleModel.cpp
 *
 *  Created on: 1 Nov 2017
 *      Author: daniel
 */

#include "UblasIncludes.hpp"
#include "HippoDynamicWntCellCycleModel.hpp"
#include "Debug.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include <cxxtest/TestSuite.h>

HippoDynamicWntCellCycleModel::HippoDynamicWntCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
   : AbstractWntOdeBasedCellCycleModel(pOdeSolver),
	 mMutantCellCIThreshold(DOUBLE_UNSET),
	 mHealthyCellCIThreshold(DOUBLE_UNSET),
	 mEquilibriumVolume(DOUBLE_UNSET),
	 mWntThreshold(DOUBLE_UNSET)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<HippoDynamicWntCellCycleModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        // Chaste solvers always check for stopping events, CVODE needs to be instructed to do so
        mpOdeSolver->CheckForStoppingEvents();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<HippoDynamicWntCellCycleModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.00005);
#endif //CHASTE_CVODE
    }
}


void HippoDynamicWntCellCycleModel::InitialiseOdeSystem(double wntConcentration, boost::shared_ptr<AbstractCellMutationState> pMutationState)
{
    mpOdeSystem = new HippoWntContactInhibitionODESystem(wntConcentration, pMutationState);
}

AbstractCellCycleModel* HippoDynamicWntCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell-cycle model
	HippoDynamicWntCellCycleModel* p_model = new HippoDynamicWntCellCycleModel(mpOdeSolver);

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
		wnt_level = (base_wnt-threshold_new)/mWntThreshold;
	}

    p_model->SetOdeSystem(new HippoWntContactInhibitionODESystem(wnt_level, mpCell->GetMutationState()));
    p_model->SetStateVariables(mpOdeSystem->rGetStateVariables());

	double cell_volume = mpCell->GetCellData()->GetItem("volume"); // CURRENT CELL VOUME
	double h_thresh = mHealthyCellCIThreshold; // HEALTHY CELL CI THRESHOLD
	double m_thresh = mMutantCellCIThreshold; // MUTANT CELL CI THRESHOLD

	static_cast<HippoWntContactInhibitionODESystem*>(mpOdeSystem)->SetCellParameters(cell_volume, h_thresh, m_thresh, mEquilibriumVolume); // Set the parameters for the CC

    return p_model;
}


void HippoDynamicWntCellCycleModel::AdjustOdeParameters(double currentTime){

	if ((mMutantCellCIThreshold == DOUBLE_UNSET) || (mHealthyCellCIThreshold == DOUBLE_UNSET))
	{
	         EXCEPTION("The member variables mQuiescentVolumeFraction and mEquilibriumVolume have not yet been set.");
	}

	//TS_ASSERT_LESS_THAN(mEquilibriumVolume,1);
	//TS_ASSERT_LESS_THAN(mHealthyCellCIThreshold,1);
	//TS_ASSERT_LESS_THAN(mMutantCellCIThreshold,1);
	//TS_ASSERT_LESS_THAN(mWntThreshold,1);

	HippoDynamicWntCellCycleModel* p_model = static_cast<HippoDynamicWntCellCycleModel*>(mpCell->GetCellCycleModel()); // Static cast to get cellcycle info

	double cell_volume = mpCell->GetCellData()->GetItem("volume"); // CURRENT CELL VOUME

	double h_thresh = mHealthyCellCIThreshold; // HEALTHY CELL CI THRESHOLD
	double m_thresh = mMutantCellCIThreshold; // MUTANT CELL CI THRESHOLD

	static_cast<HippoWntContactInhibitionODESystem*>(mpOdeSystem)->SetCellParameters(cell_volume, h_thresh, m_thresh,mEquilibriumVolume);

	// Pass this time step's Wnt stimulus into the solver as a constant over this timestep.
	double base_wnt = GetWntLevel();
	double threshold_new = 1-mWntThreshold;
	double wnt_level = 0;

	if (base_wnt < threshold_new){
		wnt_level = 0;
	}
	else{
		wnt_level = (base_wnt-threshold_new)/mWntThreshold;
	}
	mpOdeSystem->rGetStateVariables()[14] = wnt_level;

	// Use the cell's current mutation status as another input
	static_cast<HippoWntContactInhibitionODESystem*>(mpOdeSystem)->SetMutationState(mpCell->GetMutationState());

	ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();
}

void HippoDynamicWntCellCycleModel::SetHealthyCIThreshold(double healthy_cells){

	mHealthyCellCIThreshold = healthy_cells;

}

void HippoDynamicWntCellCycleModel::SetMutantCIThreshold(double mutant_cells){

	mMutantCellCIThreshold = mutant_cells;

}

double HippoDynamicWntCellCycleModel::GetHealthyCIThreshold(){

	double threshold = mHealthyCellCIThreshold;
	return threshold;

}
double HippoDynamicWntCellCycleModel::GetMutantCIThreshold(){

	double threshold = mMutantCellCIThreshold;
	return threshold;

}
double HippoDynamicWntCellCycleModel::GetE2F1Level(){

	assert(mpOdeSystem != NULL);
	double e2f1 = mpOdeSystem->rGetStateVariables()[2];
	return e2f1;

}

void HippoDynamicWntCellCycleModel::SetEquilibriumVolume(double equilibrium_volume){

	mEquilibriumVolume = equilibrium_volume;

}

void HippoDynamicWntCellCycleModel::ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel()
{
    assert(mpOdeSystem != NULL);
    assert(mpCell != NULL);

    double beta_catenin_level =   mpOdeSystem->rGetStateVariables()[12];

    // For mitogenic stimulus of 1/25.0 in Wnt equations
    if (beta_catenin_level < 10.188)
    {
        /*
         * This method is usually called within a CellBasedSimulation, after the CellPopulation
         * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
         * CellPropertyRegistry::Instance() here when setting the CellProliferativeType, we
         * would be creating a new CellPropertyRegistry. In this case the cell proliferative
         * type counts, as returned by AbstractCellPopulation::GetCellProliferativeTypeCount(),
         * would be incorrect. We must therefore access the CellProliferativeType via the cell's
         * CellPropertyCollection.
         */

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


void HippoDynamicWntCellCycleModel::Initialise()
{
    assert(mpOdeSystem == NULL);
    assert(mpCell != NULL);

    // Pass this time step's Wnt stimulus into the solver as a constant over this timestep.
	double base_wnt = GetWntLevel();
	double threshold_new = 1-mWntThreshold;
	double wnt_level = 0;

	if (base_wnt < threshold_new){
		wnt_level = 0;
	}
	else{
		wnt_level = (base_wnt-threshold_new)/mWntThreshold;
	}

    InitialiseOdeSystem(wnt_level, mpCell->GetMutationState());

    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());

    ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();

}

void HippoDynamicWntCellCycleModel::SetWntThreshold(double wnt_threshold){

	mWntThreshold = wnt_threshold;

}

double HippoDynamicWntCellCycleModel::GetWntThreshold(){

	return mWntThreshold;

}

double HippoDynamicWntCellCycleModel::GetMembraneBoundBetaCateninLevel()
{
    return mpOdeSystem->rGetStateVariables()[10];
}

double HippoDynamicWntCellCycleModel::GetCytoplasmicBetaCateninLevel()
{
    return  mpOdeSystem->rGetStateVariables()[7] + mpOdeSystem->rGetStateVariables()[8]
          + mpOdeSystem->rGetStateVariables()[15];
}

double HippoDynamicWntCellCycleModel::GetNuclearBetaCateninLevel()
{
    return mpOdeSystem->rGetStateVariables()[12];
}

void HippoDynamicWntCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output

    // Call method on direct parent class
    AbstractWntOdeBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(HippoDynamicWntCellCycleModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(HippoDynamicWntCellCycleModel)

