/*

Copyright (c) 2005-2019, University of Oxford.
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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */

#ifndef TESTDELTANOTCHREPORTERPROTRUSIONLIMI_HPP_
#define TESTDELTANOTCHREPORTERPROTRUSIONLIMI_HPP_

/*
 * = An example showing how to run Delta/Notch simulations with cellular protrusions =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to simulate a growing cell monolayer culture
 * into which a model of Delta/Notch/Reporter signalling with mutual inactivation is incorporated. 
 * This model was developed by Sprinzak et al. (Sprinzak D, Lakhanpal A, LeBon L, Garcia-Ojalvo J, Elowitz MB (2011) 
 * "Mutual Inactivation of Notch Receptors and Ligands Facilitates Developmental Patterning". 
 * PLoS Comput Biol 7(6): e1002069. https://doi.org/10.1371/journal.pcbi.1002069) 
 * This comprises three ODEs to describe the evolution in concentrations of Delta, Notch, and a 
 * Reporter of Notch Intracellular Domain in each cell. 
 * The possibility of signalling through cellular protrusions is included, as based on models by
 * Baum et al. (2016) and Bajpai et al. (2020). The ODEs
 * for Notch and Delta include a reaction term that depends on the mean Delta and mean Notch
 * concentration among cells in junctional contact and in protrusional. 
 * Thus in this simulation each cell needs to be able to access information about its neighbours 
 * and its protrusional contacts. We use the {{{CellData}}} class to facilitate this, 
 * and introduce a new CellProperty for protruional contacts. 
 * A tracking modifier "DeltaNotchReporterProtrusionTrackingModifier" is introduced to update 
 * the relevant parameters at each time step.
 *
 * EMPTYLINE
 *
 * == The test ==
 *
 * EMPTYLINE
 *
 * As in previous tutorials, we begin by including the necessary header files. We have
 * encountered these files already. Recall that often, either {{{CheckpointArchiveTypes.hpp}}}
 * or {{{CellBasedSimulationArchiver.hpp}}} must be included the first Chaste header.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellAgesWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"

#include "DeltaNotchReporterProtrusionSrnModelLimi.hpp"
#include "DeltaNotchReporterProtrusionTrackingModifier.hpp"

class TestDeltaNotchReporterProtrusionLimi : public AbstractCellBasedTestSuite
{
public:

    void TestVertexBasedMonolayerWithDeltaNotch()
    {
        /* We include the next line because Vertex simulations cannot be run in parallel */
        EXIT_IF_PARALLEL;

        /* First we create a vertex mesh to run our simulations in. A relatively large mesh is required
        for the patterning impact to appear. 
        */
        HoneycombVertexMeshGenerator generator(21, 21);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        
        // // option to create a cylindrical vertex mesh for periodicity in the x-direction. 
        // CylindricalHoneycombVertexMeshGenerator generator(21, 21);
        // Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        /* We then create some cells, each with a cell-cycle model, {{{UniformG1GenerationalCellCycleModel}}} and a subcellular reaction network model
         * {{{DeltaNotchReporterProtrusionSrnModelLimi}}}, which incorporates a Delta/Notch/Reporter ODE system 
         * with protrusional contacts and mutual inactivation,
         * as developed by Sprinzak et al (2011) and Baum et al (2016)
         * We choose to initialise the concentrations of Notch and Delta to random levels in [0, 1] 
         * in each cell. Similarly the concentration of Reporter in each cell is randomly initialized, but
         * to a very low value.
         * In this example we choose to make each cell differentiated,
         * so that no cell division occurs. */
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();
            p_cc_model->SetDimension(2);

            /* We choose to initialise the concentrations of Notch and Delta to random levels 
            in each cell. Similarly the concentration of Reporter in each cell is randomly initialized, but
            to a very low value.*/
            std::vector<double> initial_conditions;
            initial_conditions.push_back(0.1*RandomNumberGenerator::Instance()->ranf());
            initial_conditions.push_back(0.1*RandomNumberGenerator::Instance()->ranf());
            initial_conditions.push_back(0.01*RandomNumberGenerator::Instance()->ranf());
            DeltaNotchReporterProtrusionSrnModelLimi* p_srn_model = new DeltaNotchReporterProtrusionSrnModelLimi();
            p_srn_model->SetInitialConditions(initial_conditions);

            /* We additionally add the Cell Property protrusion contacts which contains the indices of all
            cells that this cell is in protrusion-mediated contact with. */
            CellPropertyCollection collection;
            MAKE_PTR(CellProtrusionContacts, p_protrusion_contacts);
            collection.AddProperty(p_protrusion_contacts);
            
            CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model, false, collection));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        /* Using the vertex mesh and cells, we create a cell-based population object, and specify which results to
         * output to file. */
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        /* We are now in a position to create and configure the cell-based simulation object, pass a force law to it,
         * and run the simulation. We can make the simulation run for longer to see more patterning by increasing the end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("AllDirL1Tl1.0inN0.9P0.1outN0.01P0.1end10");
        simulator.SetSamplingTimestepMultiple(20);
        simulator.SetEndTime(10.0);

        /* Then, we define the modifier class, which automatically updates the properties of the cells and passes it to the simulation.*/
        MAKE_PTR(DeltaNotchReporterProtrusionTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        /* This modifier assigns target areas to each cell, which are required by the {{{NagaiHondaForce}}}.
         */
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);
        simulator.Solve();
    }

};

#endif /*TESTDELTANOTCHREPORTERPROTRUSIONLIMI_HPP_*/
