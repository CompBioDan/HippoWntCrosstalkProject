/**
 * CellNanogWriter.cpp
 *
 * Cell writer class used to track the Nanog concentration levels for all the cells in "TestAgentBasedNanogExpression.hpp".
 *
 * The output of this class is a data file (cell_nanog.dat) that contains the cell type, cell ID, cell position and nanog
 * concentration levels within each individual cell.
 *
 * Format of data output: (1) Time (2) Cell Type (3) Cell ID (4) Cell X Position (5) Cell Y Position (6) Nanog Concentration
 *
 * Written by: Daniel Ward - dw0293@bristol.ac.uk
 *
 **/

#include "CellODEWriter.hpp"

#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellODEWriter<ELEMENT_DIM, SPACE_DIM>::CellODEWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("cell_ODE.dat")
{
    this->mVtkCellDataName = "ODE Level";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellODEWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double e2f1 = pCell->GetCellData()->GetItem("e2f1");
    return e2f1;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellODEWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{

	// Output a marker depicting the cell type (0 = differentiated cell, 1 = proliferating cell, 2 = other cell)
	if (pCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()){
		*this->mpOutStream << 0 << " ";
	}else if (pCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>()){
		*this->mpOutStream << 1 << " ";
	}else {
		*this->mpOutStream << 2 << " ";
	}

    // Output this cell's ID
    unsigned cell_id = pCell->GetCellId();
    *this->mpOutStream << cell_id << " ";

    // Output the position of this cell's centre
    c_vector<double, SPACE_DIM> centre_location = pCellPopulation->GetLocationOfCellCentre(pCell);
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << centre_location[i] << " ";
    }

    // Output this cell's level of nanog
    double e2f1 = pCell->GetCellData()->GetItem("e2f1");
    *this->mpOutStream << e2f1 << " ";

    // Output this cell's level of nanog
	double vol = pCell->GetCellData()->GetItem("volume");
	*this->mpOutStream << vol << " ";


}

// Explicit instantiation
template class CellODEWriter<1,1>;
template class CellODEWriter<1,2>;
template class CellODEWriter<2,2>;
template class CellODEWriter<1,3>;
template class CellODEWriter<2,3>;
template class CellODEWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellODEWriter)
