/*
 * HippoStaticWntCellCycleModel.hpp
 *
 *  Created on: 1 Nov 2017
 *      Author: daniel
 */

#ifndef PROJECTS_DANIELW_SRC_CELL_CYCLES_HippoStaticWntCellCycleModel_HPP_
#define PROJECTS_DANIELW_SRC_CELL_CYCLES_HippoStaticWntCellCycleModel_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <cfloat>

#include "AbstractOdeSystem.hpp"
#include "AbstractWntOdeBasedCellCycleModel.hpp"
#include "HippoWntContactInhibitionODESystem.hpp"
#include "AbstractCellMutationState.hpp"

/**
 * An ODE dependent cell cycle model that depends on the Wnt concentration as well as a Hippo signal dependent on cell volume
 */
class HippoStaticWntCellCycleModel : public AbstractWntOdeBasedCellCycleModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractWntOdeBasedCellCycleModel>(*this);
		archive & mHealthyCellCIThreshold;
		archive & mMutantCellCIThreshold;
		archive & mEquilibriumVolume;
		archive & mWntThreshold;
    }

    double mHealthyCellCIThreshold; // Contact Inhibition volume threshold for healthy (wild type) cells
    double mMutantCellCIThreshold; // Contact inhibition volume threshold for our mutant cells
    double mEquilibriumVolume;
    double mWntThreshold;
    /**
     * Called by Initialise() and UpdateCellProliferativeType() only.
     * Updates mCellProliferativeType to match mpOdeSystem's beta-catenin levels
     *
     * This carries out the work for UpdateCellProliferativeType();
     * But does not check the current time so it can be used by the initialise method.
     */
    void ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();

    /**
     * Adjust any ODE parameters needed before solving until currentTime.
     *
     * @param currentTime  the time up to which the system will be solved.
     */
    void AdjustOdeParameters(double currentTime);

public:

    /**
     * Default constructor.
     *
     * @param pOdeSolver An optional pointer to a cell-cycle model ODE solver object (allows the use of different ODE solvers)
     */
    HippoStaticWntCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * See AbstractCellCycleModel::Initialise()
     *
     * In this case we set up a new ODE system for a daughter cell.
     */
    void Initialise();

    /**
     * @return the level of membrane bound beta-catenin. To be used in cell-cell adhesion calculations.
     */
    double GetMembraneBoundBetaCateninLevel();

    /**
     * @return the level of cytoplasmic beta-catenin (including ubiquitinated - awaiting degradation)
     */
    double GetCytoplasmicBetaCateninLevel();

    /**
     * @return the level of nuclear beta-catenin. To be used in transcription
     */
    double GetNuclearBetaCateninLevel();

    /**
     * Pure virtual method to be implemented in concrete classes, which
     * should should allocate the mOdeSystem variable
     *
     * @param wntConcentration Wnt concentration
     * @param pMutationState Mutation state
     */
    void InitialiseOdeSystem(double wntConcentration, boost::shared_ptr<AbstractCellMutationState> pMutationState);

    /**
	 * Overridden builder method to create new copies of
	 * this cell-cycle model.
	 * @return the new cell-cycle model
	 */
	AbstractCellCycleModel* CreateCellCycleModel();

	void SetCIThresholds(double healthy_cell_threshold, double mutant_cell_threshold);

	void SetWntThreshold(double wnt_threshold);
	double GetWntThreshold();

	double GetHealthyCIThreshold();
	double GetMutantCIThreshold();

	void SetHealthyCIThreshold(double healthy_cell_threshold);
	void SetMutantCIThreshold(double mutant_cell_threshold);

	double GetE2F1Level();

	void SetEquilibriumVolume(double equilibrium_volume);
	void ResetForDivision();

    /**
     * Outputs cell-cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(HippoStaticWntCellCycleModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(HippoStaticWntCellCycleModel)

#endif /* PROJECTS_DANIELW_SRC_CELL_CYCLES_HippoStaticWntCellCycleModel_HPP_ */
