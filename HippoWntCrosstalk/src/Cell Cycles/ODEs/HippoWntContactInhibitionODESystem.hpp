/*
 * HippoWntContactInhibitionODESystem.hpp
 *
 *  Created on: 1 Nov 2017
 *      Author: daniel
 */

#ifndef PROJECTS_DANIELW_SRC_CELL_CYCLES_ODES_HIPPOWNTCONTACTINHIBITIONODESYSTEM_HPP_
#define PROJECTS_DANIELW_SRC_CELL_CYCLES_ODES_HIPPOWNTCONTACTINHIBITIONODESYSTEM_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <cmath>
#include <iostream>

#include "AbstractOdeSystem.hpp"
#include "AbstractCellMutationState.hpp"
#include "MathsCustomFunctions.hpp"

/**
 * Represents a modified version of the van Leeuwen et al. (2007) system of ODEs
 * [doi:10.1016/j.jtbi.2007.01.019]
 * coupled to the Swat et al. cell-cycle model equations.
 * [doi:10.1093/bioinformatics/bth110] which includes a form of beta-catenin
 * that has coupled with phosphorylated YAP preventing nuclear sequestration
 *
 * The variables are
 *
 *   0. r = pRb
 *   1. e = E2F1 (This is the S-phase indicator)
 *   2. i = CycD (inactive)
 *   3. j = CycD (active)
 *   4. p = pRb-p
 *   5. D = APC destruction complex
 *   6. X = Axin
 *   7. Cu = Beta Cat marked for ubiquitination
 *   8. Co = Open form Beta Cat
 *   9. A = Free Adhesion molecules
 *   10. Ca = BetaCat/Adhesion
 *   11. T = free TCF
 *   12. Cot = Open BetaCat/TCF
 *   13. Y = Wnt Target protein
 *   14. Wnt level
 *   15. Cy = Beta-cat/pYAP
 */
class HippoWntContactInhibitionODESystem : public AbstractOdeSystem
{
private:

    /**
     * Parameters for the Swat et al. (2004) model
     */

    /** Dimensional parameter k_2. */
    double mk2d;
    /** Dimensional parameter k_3. */
    double mk3d;
    /** Dimensional parameter k_34. */
    double mk34d;
    /** Dimensional parameter k_2. */
    double mk43d;
    /** Dimensional parameter k_23. */
    double mk23d;
    /** Dimensional parameter a. */
    double mad;
    /** Dimensional parameter J_11. */
    double mJ11d;
    /** Dimensional parameter J_12. */
    double mJ12d;
    /** Dimensional parameter J_13. */
    double mJ13d;
    /** Dimensional parameter J_13. */
    double mJ61d;
    /** Dimensional parameter J_62. */
    double mJ62d;
    /** Dimensional parameter J_63. */
    double mJ63d;
    /** Dimensional parameter K_m1. */
    double mKm1d;
    /** Dimensional parameter k_p. */
    double mkpd;
    /** Dimensionless parameter phi_r. */
    double mphi_r;
    /** Dimensionless parameter phi_i. */
    double mphi_i;
    /** Dimensionless parameter phi_j. */
    double mphi_j;
    /** Dimensionless parameter phi_p. */
    double mphi_p;
    /** Dimensional parameter k_16. */
    double mk16d;
    /** Dimensional parameter k_61. */
    double mk61d;
    /** Dimensionless parameter phi_E2F1. */
    double mPhiE2F1;

    /**
     * Parameters for the Van Leeuwen et al. (2007) model and Ward et al. (2017) model
     */

    /** Dimensionless parameter s_A. */
    double mSa;
    /** Dimensionless parameter s_CA. */
    double mSca;
    /** Dimensionless parameter s_C. */
    double mSc;
    /** Dimensionless parameter s_CT. */
    double mSct;
    /** Dimensionless parameter s_D. */
    double mSd;
    /** Dimensionless parameter s_T. */
    double mSt;
    /** Dimensionless parameter s_X. */
    double mSx;
    /** Dimensionless parameter s_Y. */
    double mSy;
    /** Dimensionless parameter d_A. */
    double mDa;
    /** Dimensionless parameter d_CA. */
    double mDca;
    /** Dimensionless parameter d_C. */
    double mDc;
    /** Dimensionless parameter d_CT. */
    double mDct;
    /** Dimensionless parameter d_D. */
    double mDd;
    /** Dimensionless parameter d_Dx. */
    double mDdx;
    /** Dimensionless parameter d_T. */
    double mDt;
    /** Dimensionless parameter d_U. */
    double mDu;
    /** Dimensionless parameter d_X. */
    double mDx;
    /** Dimensionless parameter d_Y. */
    double mDy;
    /** Dimensionless parameter K_c. */
    double mKc;
    /** Dimensionless parameter K_D. */
    double mKd;
    /** Dimensionless parameter K_T. */
    double mKt;
    /** Dimensionless parameter p_c. */
    double mPc;
    /** Dimensionless parameter p_u. */
    double mPu;
    /** Dimensionless parameter xi_D. */
    double mXiD;
    /** Dimensionless parameter xi_Dx. */
    double mXiDx;
    /** Dimensionless parameter xi_X. */
    double mXiX;
    /** Dimensionless parameter xi_C. */
    double mXiC;

    double mScaScaled;
    double mSctScaled;
    double mPh;
    double mKh;

    /** The mutation state of the cell. */
    boost::shared_ptr<AbstractCellMutationState> mpMutationState;

    /** The Wnt level (this affects the ODE system). */
    double mWntLevel;

    /** The volume of the cell */
    double mCellVolume;
    double mCellAge;
    double mHealthyCellCIThreshold;
    double mMutantCellCIThreshold;
    double mProliferativeType;
    double mEquilibriumVolume;

    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

public:

    /**
     * Constructor.
     * @param wntLevel is a non-dimensional Wnt value between 0 and 1. This sets up the Wnt pathway in its steady state.
     * @param pMutationState cell mutation; some affect the ODE system
     * @param stateVariables optional initial conditions for state variables (only used in archiving)
     */
    HippoWntContactInhibitionODESystem(double wntLevel = 0.0,
                                            boost::shared_ptr<AbstractCellMutationState> pMutationState=boost::shared_ptr<AbstractCellMutationState>(),
                                            std::vector<double> stateVariables=std::vector<double>());

    /**
     * Destructor.
     */
    ~HippoWntContactInhibitionODESystem();

    /**
     * Initialise parameter values.
     */
    void Init();

    /**
     * Set the mutation state of the cell.
     *
     * This should be called by the relevant cell-cycle model before any solving
     * of the ODE system (as it is used to evaluate the Y derivatives).
     *
     * @param pMutationState the mutation state.
     */
    void SetMutationState(boost::shared_ptr<AbstractCellMutationState> pMutationState);

    /**
     * Called by the archive function on the Wnt cell-cycle model.
     *
     * @return #mpMutationState the mutation state of the cell.
     */
    const boost::shared_ptr<AbstractCellMutationState> GetMutationState() const;

    /**
     * Compute the RHS of the system of ODEs.
     *
     * Returns a vector representing the RHS of the ODEs at each time step, y' = [y1' ... yn'].
     * An ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time used to evaluate the RHS.
     * @param rY value of the solution vector used to evaluate the RHS.
     * @param rDY filled in with the resulting derivatives (using van Leeuwen et al. (2007) system of equations)
     */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

    /**
     * Calculate whether the conditions for the cell cycle to finish have been met.
     *
     * @param time at which to calculate whether the stopping event has occurred
     * @param rY value of the solution vector used to evaluate the RHS
     *
     * @return whether or not stopping conditions have been met
     */
    bool CalculateStoppingEvent(double time, const std::vector<double>& rY);

    /**
     * When using CVODE this function is called instead of CalculateStoppingEvent.
     * It allows the point at which rY[1] reaches 1 to be found to greater precision.
     *
     * @param time at which to calculate whether the stopping event has occurred
     * @param rY value of the solution vector used to evaluate the RHS
     *
     * @return How close we are to the root of the stopping condition
     */
    double CalculateRootFunction(double time, const std::vector<double>& rY);

    /**
     * @return #mWntLevel.
     */
    double GetWntLevel() const;

    /**
     * @return #mCellVolme.
     */
    void SetCellParameters(double cell_volume, double healthy_threshold, double mutant_threshold, double equilibrium_volume);

    double GetCellVolume() const;

    /**
     * @return #mHypothesis.
     */
    unsigned GetHypothesis() const;
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(HippoWntContactInhibitionODESystem)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a AdaptedLeeuwenCICellCycleOdeSystem.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const HippoWntContactInhibitionODESystem * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const double wnt_level = t->GetWntLevel();
    ar & wnt_level;

    const boost::shared_ptr<AbstractCellMutationState> p_mutation_state = t->GetMutationState();
    ar & p_mutation_state;

    const std::vector<double> state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}

/**
 * De-serialize constructor parameters and initialise a AdaptedLeeuwenCICellCycleOdeSystem.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, HippoWntContactInhibitionODESystem * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    double wnt_level;
    ar & wnt_level;

    boost::shared_ptr<AbstractCellMutationState> p_mutation_state;
    ar & p_mutation_state;

    std::vector<double> state_variables;
    ar & state_variables;

    // Invoke inplace constructor to initialise instance
    ::new(t)HippoWntContactInhibitionODESystem(wnt_level, p_mutation_state, state_variables);
}
}
} // namespace ...

#endif /* PROJECTS_DANIELW_SRC_CELL_CYCLES_ODES_HIPPOWNTCONTACTINHIBITIONODESYSTEM_HPP_ */
