/*
 * TestDivisionBasedWntGradient.hpp
 *
 *  Created on: 24 Feb 2018
 *      Author: chastebox
 */

#ifndef PROJECTS_DANIELW_TEST_TESTDIVISIONBASEDWNTGRADIENT_HPP_
#define PROJECTS_DANIELW_TEST_TESTDIVISIONBASEDWNTGRADIENT_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"


#include "LinearSpringWithVariableSpringConstantsForce.hpp"
#include "CryptSimulation2dWithCryptInvasionStoppingEvent.hpp"

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
#include "WntTrackingModifier.hpp"

#include "CellAgesWriter.hpp"
#include "CellODEWriter.hpp"

#include "HippoWntInheritedCellCycleModel.hpp"
#include "SloughingCellKiller.hpp"

/*
 * Wnt gradient generation test, using a division based spreading of Wnt.
 */
class TestDivisionBasedWntGradient : public AbstractCellBasedTestSuite
{
public:

    /*
     * Here we test the formation of a Wnt gradient through the division based spreading of Wnt from a
     * Wnt high region at the base of the crypt spreading towards the top of the crypt.
     *
     * We must choose:
     *
     *  - The size of the reservoir of Wnt within the crypt (as a % of the crypt axis)
     *  - The bias of Wnt division between the daughter cells following mitosis
     *  - The threshold volume at which cells undergo contact inhibition
     *
     * We then simulate the crypt.
     *
     */
    void TestForWntGradient() throw (Exception)
    {
        // Set parameters for this compile
        unsigned max_num_experiments = 150; // Never do more than this number of experiments.
        double load_time = 300;   // time at which each simulation must be first loaded
        double maximum_duration = 1000;  // time to run the simulation for

        // Compilation information
        std::cout << "Compiled from Chaste revision number: " << ChasteBuildInfo::GetVersionString() << "\n\n";

        CommandLineArguments* p_args = CommandLineArguments::Instance();
        unsigned argc = *(p_args->p_argc); // has the number of arguments, and
        char **argv = *(p_args->p_argv); // is a char** of them.
        std::cout << "#" << argc-1 << " arguments supplied.\n" << std::flush;

        if (argc != 4)
        {
            std::cerr << "TestCryptInvasion::Please input three arguments\n"
            			 "* Wnt region size (% of crypt)\n"
                         "* Bias to wnt allocation on division\n"
                         "* Cell contact inhibition volume threshold\n" << std::flush;
            return;
        }

        // Set the input variables
        double wnt_threshold_value = atof(argv[1]);
        double wnt_division_bias = atof(argv[2]);
        double healthy_cell_ci_threshold = atof(argv[3]);
        double mutant_cell_ci_threshold = atof(argv[3]);

        std::cout << "Wnt threshold = " << atof(argv[1]) << "\n" << std::flush;
        std::cout << "Wnt bias on division = " << atof(argv[2]) << "\n" << std::flush;
        std::cout << "Contact inhibition threshold = " << healthy_cell_ci_threshold << "\n" << std::flush;
        std::string archive_directory_to_copy;

        // Change directory_to_copy_from according to which archive
        // (i.e. choice of geometry and cell cycle model) is used
        archive_directory_to_copy = "projects/DanielW/test/data/SteadyStateHippoInherited/sunter3_archive";

        FileFinder directory_to_copy_from(archive_directory_to_copy, RelativeTo::ChasteSourceRoot);
        TS_ASSERT(directory_to_copy_from.IsDir());

        // Set output directory
        std::string arg1_string(argv[1]);
        std::string arg2_string(argv[2]);
        std::string arg3_string(argv[3]);
        std::string simParams = arg1_string + "_" + arg2_string + "_" + arg3_string;
        std::string output_directory = "HippoWntGradientTest_" + simParams; // output directory name includes the input variables
        std::string mainSolutions = "/home/daniel/Documents/chasteSolutions/HippoWntGradientTest_"; // set the output directory


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

        // Create matlab results data file
		std::string newDataFile1 = "/home/daniel/Documents/chasteSolutions/matlabDataFileHippoWntGradient";
		std::string all_save_file = newDataFile1 + "/" + simParams;
		boost::filesystem::path dir1(newDataFile1);
		boost::filesystem::path dir2(all_save_file);
		boost::filesystem::create_directory(dir1);
		boost::filesystem::create_directory(dir2);


        // Run experiments
        bool can_reuse_previous_simulation = false;
        double previous_simulation_end_time = load_time;
        for (unsigned i=0; i<max_num_experiments; i++)
        {
            std::cout << "EXPERIMENT = " << i << "\n" << std::flush;

            // Load new initial condition for simulation
            CryptSimulation2dWithCryptInvasionStoppingEvent* p_simulator;

            p_simulator = CellBasedSimulationArchiver<2, CryptSimulation2dWithCryptInvasionStoppingEvent>::Load(output_directory, previous_simulation_end_time);
            p_simulator->SetOutputDirectory(output_directory);

            // Configure simulation
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
				p_Model ->SetEquilibriumVolume(0.8);
				p_Model ->SetWntThreshold(wnt_threshold_value);
				p_Model ->SetWntNoise(wnt_division_bias);

            }

            // Configure simulation
            p_simulator->SetEndTime(previous_simulation_end_time + maximum_duration);
            p_simulator->SetSamplingTimestepMultiple(960);
            p_simulator->SetCheckForStoppingEvent(false);

            // Set dependence of damping constants on cell areas
            p_population->SetAreaBasedDampingConstant(false);

            double crypt_length = p_simulator->GetCryptHeight();

            // Run simulation
            p_simulator->SetCheckForStoppingEvent(false);
            p_simulator->Solve();

            // Save simulation and tidy up
            CellBasedSimulationArchiver<2, CryptSimulation2dWithCryptInvasionStoppingEvent>::Save(p_simulator);
            delete p_simulator;

			std::ostringstream load_t;
			load_t << previous_simulation_end_time;
			std::string itrack = boost::lexical_cast<std::string>(i);

			// Rename the ages file so we can read into matlab later
			std::string oldAgeFile = mainSolutions + simParams + "/results_from_time_" + load_t.str() + "/cellages.dat";
			std::string newAgeFile = all_save_file + "/cellages_" + itrack + ".dat";
			std::rename(oldAgeFile.c_str(),newAgeFile.c_str());
			// Rename the volumes file so we can read into matlab later
			std::string oldVolFile = mainSolutions + simParams + "/results_from_time_" + load_t.str() + "/cellareas.dat";
			std::string newVolFile = all_save_file + "/cellareas_" + itrack + ".dat";
			std::rename(oldVolFile.c_str(),newVolFile.c_str());
			 // Rename the cell types file so we can read it into the matlab script later
			std::string oldTypeFile = mainSolutions + simParams + "/results_from_time_" + load_t.str() + "/results.vizcelltypes";
			std::string newTypeFile = all_save_file + "/results_" + itrack + ".vizcelltypes";
			std::rename(oldTypeFile.c_str(),newTypeFile.c_str());
			// Rename the cell types file so we can read it into the matlab script later
			std::string oldPhaseFile = mainSolutions + simParams + "/results_from_time_" + load_t.str() + "/results.vizcellphases";
			std::string newPhaseFile = all_save_file + "/results_" + itrack + ".vizcellphases";
			std::rename(oldPhaseFile.c_str(),newPhaseFile.c_str());

			std::string oldWntFile = mainSolutions + simParams + "/results_from_time_" + load_t.str() + "/cellwnt.dat";
			std::string newWntFile = all_save_file + "/cellwnt_" + itrack + ".dat";
			std::rename(oldWntFile.c_str(),newWntFile.c_str());

			can_reuse_previous_simulation = true;
			previous_simulation_end_time = SimulationTime::Instance()->GetTime();

        }

    }

};

#endif /* PROJECTS_DANIELW_TEST_TESTDIVISIONBASEDWNTGRADIENT_HPP_ */
