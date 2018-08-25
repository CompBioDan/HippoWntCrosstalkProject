/*
 * NanogTrackingModifier.cpp
 *
 *  Created on: 14 Dec 2016
 *      Author: daniel
 */

#include "WntTrackingModifier.hpp"

#include "HippoWntInheritedCellCycleModel.hpp"

template<unsigned DIM>
WntTrackingModifier<DIM>::WntTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{}

template<unsigned DIM>
WntTrackingModifier<DIM>::~WntTrackingModifier()
{}

template<unsigned DIM>
void WntTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void WntTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void WntTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    // Loop over all of the cells within the population and fetch the relevant Nanog concentrations from the NanogInheritedCellCycleModel
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
    	// Create a static cast to the cell cycle model in order to fetch the nanog concentration to be allocated to celldata
        HippoWntInheritedCellCycleModel* p_model = static_cast<HippoWntInheritedCellCycleModel*>(cell_iter->GetCellCycleModel());
        double this_wnt = p_model->GetWntConcentration();
        //double this_e2f1 = p_model->GetWntLevel();
        // Note that the state variables must be in the same order as listed in HippoWntOdeSystem
        cell_iter->GetCellData()->SetItem("wnt", this_wnt);
    }
}

template<unsigned DIM>
void WntTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class WntTrackingModifier<1>;
template class WntTrackingModifier<2>;
template class WntTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(WntTrackingModifier)



