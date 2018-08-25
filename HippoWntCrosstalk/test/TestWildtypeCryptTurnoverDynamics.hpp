/*
 * TestCryptTurnoverDynamics.hpp
 *
 *  Created on: 12 Jan 2018
 *      Author: daniel
 */

#ifndef PROJECTS_DANIELW_TEST_TESTCRYPTTURNOVERDYNAMICS_HPP_
#define PROJECTS_DANIELW_TEST_TESTCRYPTTURNOVERDYNAMICS_HPP_

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

#include "SimpleWntCellCycleModel.hpp"
#include "HippoDynamicWntCellCycleModel.hpp"
#include "HippoStaticWntCellCycleModel.hpp"

/*
 * A multiple crypt turnover test
 */
class TestCryptTurnoverDynamics: public AbstractCellBasedTestSuite
{
public:

    /*
     * Testing the wild-type turnover dynamics of the crypt with the different models of Wnt (static/dynamic) with Hippo CI
     *
     */
    void TestWildtypeTurnover() throw (Exception)
    {
        // Set parameters for this compile
        const unsigned max_num_experiments = 100; // Stop simulation after this number of experiments.
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
            std::cerr << "TestCryptInvasion::Please input three arguments\n"
                         " Cell-cycle type: 1=StaticWnt 2=DynamicWnt\n"
            			 "* Wnt threshold value\n"
                         "* Healthy cell contact inhibition volume threshold\n"  << std::flush;
            return;
        }

        // Used to be inputs now fixed
        double cc_class = atof(argv[1]);
        double wnt_threshold_value = atof(argv[2]);
        double healthy_cell_ci_threshold = atof(argv[3]);
        double mutant_cell_ci_threshold = atof(argv[3]);
		double eq_vol = 0.8;

        std::cout << "Cell-cycle class (1=STATIC,2=DYNAMIC) = " << atoi(argv[1]) << "\n" << std::flush;
        std::cout << "Wnt threshold = " << atof(argv[2]) << "\n" << std::flush;
        std::cout << "Healthy contact inhibition threshold = " << healthy_cell_ci_threshold << "\n" << std::flush;
        std::string archive_directory_to_copy;


        // Change directory_to_copy_from according to which archive
        // (i.e. choice of geometry and cell cycle model) is used
        if (cc_class==1){
        	archive_directory_to_copy = "projects/DanielW/test/data/SteadyStateHippoStatic/sunter3_archive";
        }else if (cc_class==2){
        	archive_directory_to_copy = "projects/DanielW/test/data/SteadyStateHippoDynamic/sunter3_archive";
        }

        //std::string archive_directory_to_copy = "projects/DanielW/test/data/SteadyStateHippoStatic/sunter3_archive";
        FileFinder directory_to_copy_from(archive_directory_to_copy, RelativeTo::ChasteSourceRoot);
        TS_ASSERT(directory_to_copy_from.IsDir());

        // Set output directory
        std::string arg1_string(argv[1]);
        std::string arg2_string(argv[2]);
        std::string arg3_string(argv[3]);
        std::string simParams = arg1_string + "_" + arg2_string + "_" + arg3_string;
        std::string output_directory = "HippoCryptTurnover_" + simParams;
        std::string mainCISolutions = "/home/daniel/Documents/chasteSolutions/HippoCryptTurnover_"; // main output folder


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
		std::string newDataFile1 = "/home/daniel/Documents/chasteSolutions/matlabDataFileHippoTurnover"; // output saving for multiple runs
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
            	if (cc_class==1){
					HippoStaticWntCellCycleModel* p_Model = static_cast<HippoStaticWntCellCycleModel*> (cell_iter->GetCellCycleModel());
					p_Model ->SetHealthyCIThreshold(healthy_cell_ci_threshold);
					p_Model ->SetMutantCIThreshold(mutant_cell_ci_threshold);
					p_Model ->SetEquilibriumVolume(eq_vol);
					p_Model ->SetWntThreshold(wnt_threshold_value);

				} else if (cc_class==2){
					HippoDynamicWntCellCycleModel* p_Model = static_cast<HippoDynamicWntCellCycleModel*> (cell_iter->GetCellCycleModel());
					p_Model ->SetHealthyCIThreshold(healthy_cell_ci_threshold);
					p_Model ->SetMutantCIThreshold(mutant_cell_ci_threshold);
					p_Model ->SetEquilibriumVolume(eq_vol);
					p_Model ->SetWntThreshold(wnt_threshold_value);
				}
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
            p_simulator->SetCheckForStoppingEvent(false);

            // Record what happened to the mutant population
            // (note that this is only updated every mSamplingTimestepMultiple in the Solve() method)
            boost::shared_ptr<CellPropertyRegistry> p_registry = p_simulator->rGetCellPopulation().GetCellPropertyRegistry();
            unsigned total_num_cells = p_simulator->rGetCellPopulation().GetNumRealCells();
            unsigned num_wild_type_cells = p_registry->Get<WildTypeCellMutationState>()->GetCellCount();
            unsigned num_labelled_cells = p_registry->Get<CellLabel>()->GetCellCount();

            // Save simulation and tidy up
            CellBasedSimulationArchiver<2, CryptSimulation2dWithCryptInvasionStoppingEvent>::Save(p_simulator);
            delete p_simulator;

			std::ostringstream load_t;
			load_t << previous_simulation_end_time;
			std::string itrack = boost::lexical_cast<std::string>(i);

			// Rename the ages file so we can read into matlab later
			std::string oldAgeFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/cellages.dat";
			std::string newAgeFile = all_save_file + "/cellages_" + itrack + ".dat";
			std::rename(oldAgeFile.c_str(),newAgeFile.c_str());
			// Rename the volumes file so we can read into matlab later
			std::string oldVolFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/cellareas.dat";
			std::string newVolFile = all_save_file + "/cellareas_" + itrack + ".dat";
			std::rename(oldVolFile.c_str(),newVolFile.c_str());
			 // Rename the cell types file so we can read it into the matlab script later
			std::string oldTypeFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/results.vizcelltypes";
			std::string newTypeFile = all_save_file + "/results_" + itrack + ".vizcelltypes";
			std::rename(oldTypeFile.c_str(),newTypeFile.c_str());
			// Rename the cell types file so we can read it into the matlab script later
			std::string oldPhaseFile = mainCISolutions + simParams + "/results_from_time_" + load_t.str() + "/results.vizcellphases";
			std::string newPhaseFile = all_save_file + "/results_" + itrack + ".vizcellphases";
			std::rename(oldPhaseFile.c_str(),newPhaseFile.c_str());

			can_reuse_previous_simulation = true;

        }

        WntConcentration<2>::Destroy();

    }

};

#endif /* PROJECTS_DANIELW_TEST_TESTCRYPTTURNOVERDYNAMICS_HPP_ */
