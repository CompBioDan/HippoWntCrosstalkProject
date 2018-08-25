/*

Copyright (c) 20017-2018, University of Bristol.
All rights reserved.

University of Bristol means the Chancellor, Masters and Scholars of the
University of Bristol.

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

#ifndef PROJECTS_DANIELW_TEST_TESTCRYPTINVASIONDIVISIONWNT_HPP_
#define PROJECTS_DANIELW_TEST_TESTCRYPTINVASIONDIVISIONWNT_HPP_


#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include <cxxtest/TestSuite.h>
#include "CryptSimulation2dWithCryptInvasionStoppingEvent.hpp"

#include "LinearSpringWithVariableSpringConstantsForce.hpp"

#include "CryptCellsGenerator.hpp"

#include "HoneycombMeshGenerator.hpp"
#include "CellLabel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "Version.hpp"
#include "FileFinder.hpp"

#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"

#include "VolumeTrackingModifier.hpp"

#include "SloughingCellKiller.hpp"

#include "CellAgesWriter.hpp"
#include "CellODEWriter.hpp"
#include "WntTrackingModifier.hpp"

#include "HippoWntInheritedCellCycleModel.hpp"

/*
 * A single crypt invasion experiment.
 */
class TestCryptInvasionDivisionWnt : public AbstractCellBasedTestSuite
{
public:

    /*
     * Here we take a crypt in dynamic equilibrium and bestow a mutation on a single cell.
     *
     * We must choose:
     *
     *  - at what height up the crypt (as a fraction of the crypt length) the mutation occurs
     *  - the type of mutation (this affects the cell's proliferative and/or adhesive properties)
     *  - whether we assume that each cell's damping constant is proportional to its area
     *
     * We then simulate the crypt, keeping track of the mutant cell and its progeny.
     *
     * The simulation is stopped as soon as a stopping event occurs; this is implemented in the
     * simulation method StoppingEventHasOccurred(), which is overridden in the new class
     * CryptSimulation2dWithCryptInvasionStoppingEvent.
     *
     * We then record the initial conditions of mutation, which of the stopping events occurred
     * and how long the simulation ran for.
     */
    void TestForCryptInvasion() throw (Exception)
    {
        // Set parameters for this compile
        unsigned max_num_experiments_per_success = 1000; // Maximum number of successive experiments to simulate per success
        unsigned max_experiments_possible = 15000; // Never do more than this number of experiments.
        const unsigned num_successes_aim = 500; // Stop simulation after this number of experiments.
        unsigned num_successes = 0;             // No successes yet...
        double load_time = 300;                 // time at which each simulation must be first loaded
        double maximum_duration = 1000;        // time we are prepared to wait for each simulation to end

        // Compilation information
        std::cout << "Compiled from Chaste revision number: " << ChasteBuildInfo::GetVersionString() << "\n\n";

        CommandLineArguments* p_args = CommandLineArguments::Instance();
        unsigned argc = *(p_args->p_argc); // has the number of arguments, and
        char **argv = *(p_args->p_argv); // is a char** of them.
        std::cout << "#" << argc-1 << " arguments supplied.\n" << std::flush;

        if (argc != 4)
        {
            std::cerr << "TestCryptInvasionDivisionWnt::Please input three arguments\n"
                         "* Healthy cell contact inhibition volume threshold\n"
                         "* Wnt reservoir size\n"
            			 "* Wnt division bias (bias towards mother cell when dividing Wnt following mitosis\n" << std::flush;
            return;
        }

        double bucket_size = 0.05;

        double eq_volume = 0.8; // cell equilibrium volume

		double healthy_cell_ci_threshold = atof(argv[1]);
		double mutant_cell_ci_threshold = 0.0;
		double wnt_reservoir = atof(argv[2]);
		double wnt_division_bias = atof(argv[3]);

		std::cout << "Healthy contact inhibition threshold = " << healthy_cell_ci_threshold << "\n" << std::flush;
		std::cout << "Wnt threshold level = " << wnt_reservoir << "\n" << std::flush;
		std::cout << "Wnt division bias = " << wnt_division_bias << "\n" << std::flush;

        // Change directory_to_copy_from according to which archive
        // (i.e. choice of geometry and cell cycle model) is used
        FileFinder directory_to_copy_from("projects/DanielW/test/data/SteadyStateHippoInherited/sunter3_archive", RelativeTo::ChasteSourceRoot);
        TS_ASSERT(directory_to_copy_from.IsDir());


        // Set output directory
		std::string arg1_string(argv[1]);
		std::string arg2_string(argv[2]);
		std::string arg3_string(argv[3]);
		std::string simParams = arg1_string + "_" + arg2_string + "_" + arg3_string;
		std::string output_directory = "HippoCryptInvasionDivisionWnt_" + simParams;
		std::string mainCISolutions = "/home/daniel/Documents/chasteSolutions/HippoCryptInvasionDivisionWnt_";

        // Create results file handler
        OutputFileHandler results_handler(output_directory, false);

        // Where the archive will be copied to
        std::string directory_to_copy_to = results_handler.GetChasteTestOutputDirectory() + output_directory;

        // Make a folder to put archives into
        std::string command = "mkdir -p " + directory_to_copy_to + "/archive";
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value, 0);

        // Remove existing archives to avoid confusion!
        command = "rm -rf " + directory_to_copy_to + "/archive/*";

        // Copy archive to output directory and rename so that it is recognised by CellBasedSimulationArchiver
        command = "cp -rf --remove-destination " + directory_to_copy_from.GetAbsolutePath() + "/* " + directory_to_copy_to + "/archive";
        return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value, 0);

        // Create overall results file
        std::string overall_results_filename = "overall_results" + arg1_string + "_" + arg2_string + "_" + arg3_string + ".dat";
        out_stream overall_results_file = results_handler.OpenOutputFile(overall_results_filename);
        // Create data file for analysis afterwards
        std::string data_results_filename = "just_data" + arg1_string + "_" + arg2_string + "_" + arg3_string + ".dat";
        out_stream data_results_file = results_handler.OpenOutputFile(data_results_filename);

        // Create matlab results data file
		std::string newDataFile1 = "/home/daniel/Documents/chasteSolutions/matlabDataFileHippoCryptInvasionDivisionWnt";
		std::string newDataFile2 = "/home/daniel/Documents/chasteSolutions/matlabDataFileHippoCryptInvasionDivisionWnt/Washed_" + simParams;
		std::string newDataFile3 = "/home/daniel/Documents/chasteSolutions/matlabDataFileHippoCryptInvasionDivisionWnt/Taken_" + simParams;
		std::string newDataFile4 = "/home/daniel/Documents/chasteSolutions/matlabDataFileHippoCryptInvasionDivisionWnt/No_" + simParams;
		boost::filesystem::path dir1(newDataFile1);
		boost::filesystem::path dir2(newDataFile2);
		boost::filesystem::path dir3(newDataFile3);
		boost::filesystem::path dir4(newDataFile4);
		boost::filesystem::create_directory(dir1);
		boost::filesystem::create_directory(dir2);
		boost::filesystem::create_directory(dir3);
		boost::filesystem::create_directory(dir4);
		unsigned wCount = 0; unsigned tCount = 0; unsigned nCount = 0;

        // Run experiments
        bool can_reuse_previous_simulation = false;
        double previous_simulation_end_time = load_time;
        unsigned max_num_experiments = max_num_experiments_per_success; // More experiments allowed after a success below.
        for (unsigned i=0; i<max_num_experiments; i++)
        {
            std::cout << "EXPERIMENT = " << i << "\n" << std::flush;
            *overall_results_file << "EXPERIMENT = " << i << "\n" << std::flush;

            // Load new initial condition for simulation
            CryptSimulation2dWithCryptInvasionStoppingEvent* p_simulator;

            // If the previous simulation ended with the mutant population being washed
            // out of the crypt, then we can re-use this as our next initial condition
            if (!can_reuse_previous_simulation)
            {
                // Otherwise, we need to load the original archive and run for a random length
                // of time in order to get a new initial condition
                p_simulator = CellBasedSimulationArchiver<2, CryptSimulation2dWithCryptInvasionStoppingEvent>::Load(output_directory, load_time);

                // We must reseed the random number generator, which was archived along with the simulation
                RandomNumberGenerator::Instance()->Reseed(time(NULL));

                previous_simulation_end_time = load_time + 50.0*RandomNumberGenerator::Instance()->ranf();

                std::vector<boost::shared_ptr<Cell> > cells;

                for (AbstractCellPopulation<2>::Iterator cell_iter = p_simulator->rGetCellPopulation().Begin();
					cell_iter != p_simulator->rGetCellPopulation().End();
					++cell_iter)
				{
                		HippoWntInheritedCellCycleModel* p_Model = static_cast<HippoWntInheritedCellCycleModel*> (cell_iter->GetCellCycleModel());
						p_Model ->SetHealthyCIThreshold(healthy_cell_ci_threshold);
						p_Model ->SetMutantCIThreshold(mutant_cell_ci_threshold);
						p_Model ->SetEquilibriumVolume(eq_volume);
						p_Model ->SetWntThreshold(wnt_reservoir);
						p_Model ->SetWntNoise(wnt_division_bias);

				}

				p_simulator->SetEndTime(previous_simulation_end_time);
				p_simulator->SetOutputDirectory(output_directory);
				p_simulator->Solve();

				CellBasedSimulationArchiver<2, CryptSimulation2dWithCryptInvasionStoppingEvent>::Save(p_simulator);
            }

            p_simulator = CellBasedSimulationArchiver<2, CryptSimulation2dWithCryptInvasionStoppingEvent>::Load(output_directory, previous_simulation_end_time);
            p_simulator->SetOutputDirectory(output_directory);

            // Set effect of mutation on adhesion
            MeshBasedCellPopulation<2>* p_population = static_cast<MeshBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation()));
            double mutant_damping_constant = p_population->GetDampingConstantNormal();
            p_population->SetDampingConstantMutant(mutant_damping_constant);

            for (AbstractCellPopulation<2>::Iterator cell_iter = p_simulator->rGetCellPopulation().Begin();
					cell_iter != p_simulator->rGetCellPopulation().End();
					++cell_iter)
			{

					HippoWntInheritedCellCycleModel* p_Model = static_cast<HippoWntInheritedCellCycleModel*> (cell_iter->GetCellCycleModel());
					p_Model ->SetHealthyCIThreshold(healthy_cell_ci_threshold);
					p_Model ->SetMutantCIThreshold(mutant_cell_ci_threshold);
					p_Model ->SetEquilibriumVolume(eq_volume);
					p_Model ->SetWntThreshold(wnt_reservoir);
					p_Model ->SetWntNoise(wnt_division_bias);

			}

            // Configure simulation
            p_simulator->SetEndTime(previous_simulation_end_time + maximum_duration);
            p_simulator->SetSamplingTimestepMultiple(960);
            p_simulator->SetCheckForStoppingEvent(true);

            // Set dependence of damping constants on cell areas
            p_population->SetAreaBasedDampingConstant(false);

            double crypt_length = p_simulator->GetCryptHeight();

            // Compute minimum and maximum height of initial mutation
            double minimum_mutation_height = crypt_length;
            double maximum_mutation_height = (bucket_size)*crypt_length;

            // Iterate over the tissue and generate a set of real cells located in the specified bucket
            std::vector<boost::shared_ptr<Cell> > cells_in_bucket;
            for (AbstractCellPopulation<2>::Iterator cell_iter = p_simulator->rGetCellPopulation().Begin();
                 cell_iter != p_simulator->rGetCellPopulation().End();
                 ++cell_iter)
            {
                double cell_height = p_simulator->rGetCellPopulation().GetLocationOfCellCentre(*cell_iter)[1];
                if ( cell_height >= minimum_mutation_height && cell_height < maximum_mutation_height )
                {
                    cells_in_bucket.push_back(*cell_iter);
                }
            }

            // Now choose a random cell from this set and give it a mutation
            unsigned random_index = RandomNumberGenerator::Instance()->randMod(cells_in_bucket.size());
            boost::shared_ptr<Cell> p_mutant_cell = cells_in_bucket[random_index];

            // Sanity check
            if (   p_simulator->rGetCellPopulation().GetLocationOfCellCentre(p_mutant_cell)[1] < minimum_mutation_height
                || p_simulator->rGetCellPopulation().GetLocationOfCellCentre(p_mutant_cell)[1] > maximum_mutation_height)
            {
                EXCEPTION("The height of the initial mutation is not within the prescribed limits");
            }

            // Set state of initial mutation
            boost::shared_ptr<AbstractCellProperty> p_mutation_state;

			p_mutation_state = p_simulator->rGetCellPopulation().GetCellPropertyRegistry()->Get<ApcTwoHitCellMutationState>();
			p_mutant_cell->SetMutationState(p_mutation_state);

            // Store the initial height and corresponding location index of the mutant cell
            unsigned ancestor_cell = p_simulator->rGetCellPopulation().GetLocationIndexUsingCell(p_mutant_cell);
            double mutant_cell_height = p_simulator->rGetCellPopulation().GetLocationOfCellCentre(p_mutant_cell)[1];

            // Record details of initial mutation
            *overall_results_file << "Initial mutation occurred at node " << ancestor_cell << "\n" << std::flush;

            *overall_results_file << "Experiment start time = " << SimulationTime::Instance()->GetTime() << "\n" << std::flush;

            // Run simulation
            p_simulator->SetCheckForStoppingEvent(false);
            p_simulator->Solve();

            // Record actual duration of experiment (a stopping event may have occurred)
            double actual_experiment_duration = SimulationTime::Instance()->GetTime() - previous_simulation_end_time;
            *overall_results_file << "Experiment duration = " << actual_experiment_duration << "\n" << std::flush;
            *data_results_file << actual_experiment_duration << "\n" << std::flush;

            // Record how many cells are in the crypt at the end of the experiment
			unsigned num_cells_at_end = p_simulator->rGetCellPopulation().GetNumRealCells();
			*overall_results_file << "Number of cells in crypt at end = " << num_cells_at_end << "\n" << std::flush;
			*data_results_file << num_cells_at_end << "\n" << std::flush;

            // Record what happened to the mutant population
            // (note that this is only updated every mSamplingTimestepMultiple in the Solve() method)
            boost::shared_ptr<CellPropertyRegistry> p_registry = p_simulator->rGetCellPopulation().GetCellPropertyRegistry();
            unsigned total_num_cells = p_simulator->rGetCellPopulation().GetNumRealCells();
            unsigned num_wild_type_cells = p_registry->Get<WildTypeCellMutationState>()->GetCellCount();
            unsigned num_labelled_cells = p_registry->Get<CellLabel>()->GetCellCount();

            // Save simulation and tidy up
            CellBasedSimulationArchiver<2, CryptSimulation2dWithCryptInvasionStoppingEvent>::Save(p_simulator);
            delete p_simulator;

            if (num_wild_type_cells == 0 || num_labelled_cells==total_num_cells)
			{
				// If there are no healthy cells left, then the mutant progeny must have taken
				// over the entire crypt, so we must go back to the original initial condition
				// for the next simulation
				*overall_results_file << "Mutant population taken over crypt\n" << std::flush;
				*data_results_file << 2 << "\n" << std::flush;
				num_successes++;  // Keep track of how many successes we have had.
				max_num_experiments += max_num_experiments_per_success;   // Allow more experiments after this success...

				std::ostringstream load_t;
				load_t << previous_simulation_end_time;
				std::string itrack = boost::lexical_cast<std::string>(tCount);

				// Rename the ages file so we can read into matlab later
				std::string oldAgeFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/cellages.dat";
				std::string newAgeFile = newDataFile3 + "/cellages_" + itrack + ".dat";
				std::rename(oldAgeFile.c_str(),newAgeFile.c_str());
				// Rename the volumes file so we can read into matlab later
				std::string oldVolFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/cellareas.dat";
				std::string newVolFile = newDataFile3 + "/cellareas_" + itrack + ".dat";
				std::rename(oldVolFile.c_str(),newVolFile.c_str());
				// Rename the cell types file so we can read it into the matlab script later
				std::string oldTypeFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/results.vizcelltypes";
				std::string newTypeFile = newDataFile3 + "/results_" + itrack + ".vizcelltypes";
				std::rename(oldTypeFile.c_str(),newTypeFile.c_str());
				// Rename the cell types file so we can read it into the matlab script later
				std::string oldPhaseFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/results.vizcellphases";
				std::string newPhaseFile = newDataFile3 + "/results_" + itrack + ".vizcellphases";
				std::rename(oldPhaseFile.c_str(),newPhaseFile.c_str());

				std::string oldWntFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/cellwnt.dat";
				std::string newWntFile = newDataFile3 + "/cellwnt_" + itrack + ".dat";
				std::rename(oldWntFile.c_str(),newWntFile.c_str());

				tCount++;

				can_reuse_previous_simulation = false;
				previous_simulation_end_time = load_time;

				if (max_num_experiments > max_experiments_possible)
				{
					max_num_experiments = max_experiments_possible;
				}
				if (num_successes >= num_successes_aim)
				{
					break; // End experiments
				}
			}
			else if (num_wild_type_cells == total_num_cells && num_labelled_cells == 0)
			{
				// If there are no mutant cells left, then the mutant progeny must have been
				// washed out of the crypt, so we can use the end of this simulation as the
				// initial condition for the next simulation
				*overall_results_file << "Mutant population washed out of crypt\n" << std::flush;
				*data_results_file << 1 << "\n" << std::flush;

				std::ostringstream load_t;
				load_t << previous_simulation_end_time;
				std::string itrack = boost::lexical_cast<std::string>(wCount);

				// Rename the ages file so we can read into matlab later
				std::string oldAgeFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/cellages.dat";
				std::string newAgeFile = newDataFile2 + "/cellages_" + itrack + ".dat";
				std::rename(oldAgeFile.c_str(),newAgeFile.c_str());
				// Rename the volumes file so we can read into matlab later
				std::string oldVolFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/cellareas.dat";
				std::string newVolFile = newDataFile2 + "/cellareas_" + itrack + ".dat";
				std::rename(oldVolFile.c_str(),newVolFile.c_str());
				 // Rename the cell types file so we can read it into the matlab script later
				std::string oldTypeFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/results.vizcelltypes";
				std::string newTypeFile = newDataFile2 + "/results_" + itrack + ".vizcelltypes";
				std::rename(oldTypeFile.c_str(),newTypeFile.c_str());
				// Rename the cell types file so we can read it into the matlab script later
				std::string oldPhaseFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/results.vizcellphases";
				std::string newPhaseFile = newDataFile2 + "/results_" + itrack + ".vizcellphases";
				std::rename(oldPhaseFile.c_str(),newPhaseFile.c_str());

				std::string oldWntFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/cellwnt.dat";
				std::string newWntFile = newDataFile2 + "/cellwnt_" + itrack + ".dat";
				std::rename(oldWntFile.c_str(),newWntFile.c_str());
				wCount++;

				can_reuse_previous_simulation = true;
				previous_simulation_end_time = SimulationTime::Instance()->GetTime();
				if (num_successes >= num_successes_aim)
				{
					break; // End experiments
				}
			}
			else
			{
				// If both healthy and mutant cells are still present in the crypt at the end
				// of the maximum allowed duration, then we must go back to the original initial
				// condition for the next simulation
				*overall_results_file << "No result\n" << std::flush;

				std::ostringstream load_t;
				load_t << previous_simulation_end_time;
				std::string itrack = boost::lexical_cast<std::string>(nCount);

				 // Rename the ages file so we can read into matlab later
				std::string oldAgeFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/cellages.dat";
				std::string newAgeFile = newDataFile4 + "/cellages_" + itrack + ".dat";
				std::rename(oldAgeFile.c_str(),newAgeFile.c_str());
				// Rename the volumes file so we can read into matlab later
				std::string oldVolFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/cellareas.dat";
				std::string newVolFile = newDataFile4 + "/cellareas_" + itrack + ".dat";
				std::rename(oldVolFile.c_str(),newVolFile.c_str());
				// Rename the cell types file so we can read it into the matlab script later
				std::string oldTypeFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/results.vizcelltypes";
				std::string newTypeFile = newDataFile4 + "/results_" + itrack + ".vizcelltypes";
				std::rename(oldTypeFile.c_str(),newTypeFile.c_str());
				// Rename the cell types file so we can read it into the matlab script later
				std::string oldPhaseFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/results.vizcellphases";
				std::string newPhaseFile = newDataFile4 + "/results_" + itrack + ".vizcellphases";
				std::rename(oldPhaseFile.c_str(),newPhaseFile.c_str());

				std::string oldWntFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/cellwnt.dat";
				std::string newWntFile = newDataFile4 + "/cellwnt_" + itrack + ".dat";
				std::rename(oldWntFile.c_str(),newWntFile.c_str());
				nCount++;

				can_reuse_previous_simulation = false;
				previous_simulation_end_time = load_time;
				if (num_successes >= num_successes_aim)
				{
					break; // End experiments
				}
			}
		}

        // Close results files and tidy up
        *overall_results_file << "EXPERIMENTS COMPLETE\n" << std::flush;
        overall_results_file->close();
        WntConcentration<2>::Destroy();
    }

};


#endif /* PROJECTS_DANIELW_TEST_TESTCRYPTINVASIONDIVISIONWNT_HPP_ */
