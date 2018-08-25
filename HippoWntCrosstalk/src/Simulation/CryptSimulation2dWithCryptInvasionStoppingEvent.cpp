/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "CryptSimulation2dWithCryptInvasionStoppingEvent.hpp"

CryptSimulation2dWithCryptInvasionStoppingEvent::CryptSimulation2dWithCryptInvasionStoppingEvent(AbstractCellPopulation<2>& rTissue,
                                                double cryptHeight,
                                                bool deleteTissueAndForceCollection,
                                                bool initialiseCells)
  : CryptSimulation2d(rTissue, deleteTissueAndForceCollection, initialiseCells),
    mCheckForStoppingEvent(false),
    mCryptHeight(cryptHeight)
{
}

bool CryptSimulation2dWithCryptInvasionStoppingEvent::StoppingEventHasOccurred()
{
    bool has_stopping_event_occurred = false;

    // Only check for stopping events if told to...
    if (mCheckForStoppingEvent)
    {
        // ...and to speed things up, only check every mSamplingTimestepMultiple timesteps
        if (SimulationTime::Instance()->GetTimeStepsElapsed()%mSamplingTimestepMultiple==0)
        {
            boost::shared_ptr<CellPropertyRegistry> p_registry = rGetCellPopulation().GetCellPropertyRegistry();

            // Get the number of cells
            unsigned num_cells = rGetCellPopulation().GetNumRealCells();
            unsigned num_wild_type_cells = p_registry->Get<WildTypeCellMutationState>()->GetCellCount();
            unsigned num_labelled_cells = p_registry->Get<CellLabel>()->GetCellCount();
            
            // Stop the simulation if there are no healthy unlabelled cells left, or no mutant/healthy labelled cells left
            if (    (num_wild_type_cells == 0 )
                 || (num_wild_type_cells == num_cells && num_labelled_cells==0        )
                 || (num_wild_type_cells == num_cells && num_labelled_cells==num_cells) )
            {
                has_stopping_event_occurred = true;
            }
        }
    }
    return has_stopping_event_occurred;
}

void CryptSimulation2dWithCryptInvasionStoppingEvent::SetCheckForStoppingEvent(bool checkForStoppingEvent)
{
    mCheckForStoppingEvent = checkForStoppingEvent;
}

bool CryptSimulation2dWithCryptInvasionStoppingEvent::GetCheckForStoppingEvent()
{
    return mCheckForStoppingEvent;
}

double CryptSimulation2dWithCryptInvasionStoppingEvent::GetCryptHeight()
{
    return mCryptHeight;
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CryptSimulation2dWithCryptInvasionStoppingEvent)
