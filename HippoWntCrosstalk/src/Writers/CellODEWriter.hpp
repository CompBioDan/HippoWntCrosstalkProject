/*
 * CellODEWriter.cpp
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
 */

#ifndef PROJECTS_NANOGSTEMCELLS2016_SRC_CELLNANOGWRITER_HPP_
#define PROJECTS_NANOGSTEMCELLS2016_SRC_CELLNANOGWRITER_HPP_

#include "AbstractCellWriter.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A class written using the visitor pattern for writing to file, for each cell,
 * the level of Nanog protein within the cell.
 *
 * The output file is called cell_nanog.dat by default. If VTK is switched on,
 * then the writer also specifies the VTK output for each cell, which is stored in
 * the VTK cell data "Nanog Level" by default.
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellODEWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
    /** Needed for serialization. */
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
        archive & boost::serialization::base_object<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    CellODEWriter();

    /**
     * Overridden GetCellDataForVtkOutput() method.
     *
     * Get a double associated with a cell. This method reduces duplication
     * of code between the methods VisitCell() and AddVtkData().
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell
     *
     * @return data associated with the cell
     */
    double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Overridden VisitCell() method.
     *
     * Visit a cell and write its delta and notch data.
     *
     * Outputs a line of space-separated values of the form:
     * ...[cell type] [cell id] [x-pos] [y-pos] [z-pos] [nanog]...
     * with [y-pos] and [z-pos] included for 2 and 3 dimensional simulations, respectively.
     *
     * This is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell
     */
    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellODEWriter)

#endif /* CellODEWriter_HPP_ */
