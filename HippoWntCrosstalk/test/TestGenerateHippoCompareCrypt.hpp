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


#ifndef PROJECTS_HIPPOWNTCROSSTALK_TEST_TESTGENERATEHIPPOCOMPARECRYPT_HPP_
#define PROJECTS_HIPPOWNTCROSSTALK_TEST_TESTGENERATEHIPPOCOMPARECRYPT_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "SunterSetup.hpp"
#include "SmartPointers.hpp"
#include "CryptCellsGenerator.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CryptSimulation2dWithCryptInvasionStoppingEvent.hpp"

#include "GeneralisedLinearSpringForce.hpp"

#include "VolumeTrackingModifier.hpp"
#include "WntTrackingModifier.hpp"

#include "SimpleWntCellCycleModel.hpp"
#include "HippoStaticWntCellCycleModel.hpp"
#include "HippoDynamicWntCellCycleModel.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisOne.hpp"
#include "HippoWntInheritedCellCycleModel.hpp"

#include "CellAgesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"

#include "SloughingCellKiller.hpp"

class TestGenerateHippoCompareCrypt : public CxxTest::TestSuite
{
public:

    /**
     * This test generates steady state crypts for use in crypt invasion simulations.
     *
     * The setup is done with the sunter_index (which gives different geometries;
     * 1 and 3 referring to different sites in the mouse colon)
     *
     * The archives which are created should be saved and used as an input into other simulations.
     *
     * This can be used to create steady states for different models
     * (1) Wnt network only dependent cell-cycle model (based on Van Leeuwen et al (2009))
     * (2) Hippo dependent cell-cycle model with dynamically updating Wnt signal (implementing contact-inhibition)
     * (3) Hippo dependent cell-cycle model with a static (birth defined) Wnt signal
     * (4) Division based spreading of Wnt
     */
	 void TestGenerateSteadyStateResultsForContactInhibitionExperiments() throw(Exception)
	{
		std::cout << "\nGenerating steady state archives...\n" << std::flush;

		// THIS IS A HACK WHICH HIJACKS PETSC TO GET ARGUMENTS INTO A TEST!
		CommandLineArguments* p_args = CommandLineArguments::Instance();
		unsigned argc = *(p_args->p_argc); // has the number of arguments, and
		char **argv = *(p_args->p_argv); // is a char** of them.
		std::cout << "#" << argc-1 << " arguments supplied.\n" << std::flush;

			   if (argc !=  3)
			   {
				   std::cerr << "TestGenerateCryptSteadyState::Please input two arguments\n"
								"* Required steady state (1=Wnt, 2=WntHippoDynamic, 3=WntHippoStatic, 4=DivisionWnt) \n"
								"* Wnt threshold level (between 0 & 1)\n" << std::flush;
				   return;
			   }

		double steady_state_class = atof(argv[1]);
		double wnt_threshold_value = atof(argv[2]);

		std::cout << "Chosen steady state type = " << steady_state_class << "\n" << std::flush;
		std::cout << "Wnt threshold = " << wnt_threshold_value << "\n" << std::flush;

		unsigned sunter_index = 3u; // the geometry to use (see above)

		double end_of_simulation = 300.0; // hours

		SunterSetup sunter_setup(sunter_index); // SunterGeometry and parameters set here
		CylindricalHoneycombMeshGenerator generator = sunter_setup.SetSunterParametersAndGetMeshGenerator();

		// Store the height of the crypt, for use by the SloughingCellKiller
		double crypt_height = generator.GetDomainDepth();

		// Create cylindrical mesh
		Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

		// Get location indices corresponding to real cells in mesh
		std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

		// Set up instance of SimulationTime singleton
		SimulationTime* p_simulation_time = SimulationTime::Instance();
		p_simulation_time->SetStartTime(0.0);

		// Set up each cell with a simple Wnt-based cell cycle model
		std::vector<boost::shared_ptr<Cell> > cells;

		if (steady_state_class == 1){
			CryptCellsGenerator<VanLeeuwen2009WntSwatCellCycleModelHypothesisOne> cell_generator;
			cell_generator.Generate(cells, p_mesh, location_indices, true);
		}
		else if (steady_state_class == 2){
			CryptCellsGenerator<HippoDynamicWntCellCycleModel> cell_generator;
			cell_generator.Generate(cells, p_mesh, location_indices, true);
		}
		else if (steady_state_class == 3){
			CryptCellsGenerator<HippoStaticWntCellCycleModel> cell_generator;
			cell_generator.Generate(cells, p_mesh, location_indices, true);
		}
		else if (steady_state_class == 4){
			CryptCellsGenerator<HippoWntInheritedCellCycleModel> cell_generator;
			cell_generator.Generate(cells, p_mesh, location_indices, true);
		}

		// cell_generator.Generate(cells, p_mesh, location_indices, true);
		sunter_setup.SetUpCellCycleModelParameters(cells);

		// Create crypt facade
		MeshBasedCellPopulation<2> crypt(*p_mesh, cells);


		// Set up instance of WntConcentration singleton and associate it with crypt
		WntConcentration<2>::Instance()->SetType(LINEAR);
		WntConcentration<2>::Instance()->SetCellPopulation(crypt);
		WntConcentration<2>::Instance()->SetCryptLength(crypt_height);


		// Create crypt simulation
		CryptSimulation2dWithCryptInvasionStoppingEvent simulator(crypt, crypt_height);

		if (steady_state_class == 2){
			for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
				cell_iter != simulator.rGetCellPopulation().End();
				++cell_iter)
			{
				HippoDynamicWntCellCycleModel* p_Model = static_cast<HippoDynamicWntCellCycleModel*> (cell_iter->GetCellCycleModel());
				p_Model ->SetHealthyCIThreshold(0.0); // set cell healthy ci volume threshold
				p_Model ->SetMutantCIThreshold(0.0); // set cell mutant ci volume threshold
				p_Model ->SetEquilibriumVolume(0.8); // set cell equilibrium volume
				p_Model ->SetWntThreshold(wnt_threshold_value); // set size of Wnt gradient in crypt
			}
		}
		if (steady_state_class == 3){
			for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
				cell_iter != simulator.rGetCellPopulation().End();
				++cell_iter)
			{
				HippoStaticWntCellCycleModel* p_Model = static_cast<HippoStaticWntCellCycleModel*> (cell_iter->GetCellCycleModel());
				p_Model ->SetHealthyCIThreshold(0.0); // set cell healthy ci volume threshold
				p_Model ->SetMutantCIThreshold(0.0); // set cell mutant ci volume threshold
				p_Model ->SetEquilibriumVolume(0.8); // set equilibrium volume of the cell
				p_Model ->SetWntThreshold(wnt_threshold_value); // set the size of the Wnt gradietn in the crypt

			}
	   }
		if (steady_state_class == 4){
			for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
				cell_iter != simulator.rGetCellPopulation().End();
				++cell_iter)
			{
				HippoWntInheritedCellCycleModel* p_Model = static_cast<HippoWntInheritedCellCycleModel*> (cell_iter->GetCellCycleModel());
				p_Model ->SetHealthyCIThreshold(0.0); // set healthy cell ci threshold volume
				p_Model ->SetMutantCIThreshold(0.0); // set mutant cell ci threshold volume
				p_Model ->SetEquilibriumVolume(0.8); // set cell equilibrium volume
				p_Model ->SetWntThreshold(wnt_threshold_value); // set threshold for Wnt reservoir in crypt
				p_Model ->SetWntNoise(0.0); // set bias in the division of Wnt following mitotis

			}
	   }

		// Set where to output simulation results
		// Set output directory depending on the steady state type we chose
		if (steady_state_class == 1){
		std::string output_directory = "WntCryptSetup";
		simulator.SetOutputDirectory(output_directory);
		}
		else if (steady_state_class == 2){
		std::string output_directory = "HippoDynamicCryptSetup";
		simulator.SetOutputDirectory(output_directory);
		}
		else if (steady_state_class == 3){
		std::string output_directory = "HippoStaticCryptSetup";
		simulator.SetOutputDirectory(output_directory);
		}
		else if (steady_state_class == 4){
		std::string output_directory = "HippoInheritedWntCryptSetup";
		simulator.SetOutputDirectory(output_directory);
		}

		// Set up force law
		MAKE_PTR(GeneralisedLinearSpringForce<2>, p_meineke_force);
		p_meineke_force->SetMeinekeSpringStiffness(30.0);
		simulator.AddForce(p_meineke_force);

		// Volume tracking modifier
		MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
		simulator.AddSimulationModifier(p_modifier);

		MAKE_PTR(WntTrackingModifier<2>,p_wnt_modifier);
		simulator.AddSimulationModifier(p_wnt_modifier);

		// Set simulation to output cell types
		crypt.AddCellWriter<CellProliferativeTypesWriter>();
		crypt.AddCellWriter<CellMutationStatesWriter>();
		crypt.AddCellWriter<CellProliferativePhasesWriter>();
		crypt.AddCellWriter<CellAgesWriter>(); // outputs the cell ages to a cellages.dat
		crypt.AddCellWriter<CellVolumesWriter>(); // outputs the cell volumes (areas in 2D) to cellareas.dat

		// Set length of simulation
		simulator.SetEndTime(end_of_simulation);

		// Only write results to file every hour
		simulator.SetSamplingTimestepMultiple(120);

		// Set up sloughing cell killer and pass in to simulation
		MAKE_PTR_ARGS(SloughingCellKiller<2>, p_cell_killer, (&simulator.rGetCellPopulation(), crypt_height));
		simulator.AddCellKiller(p_cell_killer);

		// Do not give mutant cells any different movement properties to normal ones
		crypt.SetDampingConstantMutant(crypt.GetDampingConstantNormal());

		// A small random upward force is applied to cells at the base to
		// prevent an unstable equilibrium with many cells with compressed springs
		// building up on y=0.
		simulator.UseJiggledBottomCells();

		// Don't check for stopping events whilst getting to
		// a steady state.
		simulator.SetCheckForStoppingEvent(false);

		// Run simulation
		simulator.Solve();

		// Save simulation
		CellBasedSimulationArchiver<2, CryptSimulation2dWithCryptInvasionStoppingEvent>::Save(&simulator);


	}
};

#endif /* PROJECTS_HIPPOWNTCROSSTALK_TEST_TESTGENERATEHIPPOCOMPARECRYPT_HPP_ */
