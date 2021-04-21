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

#ifndef TESTDELTANOTCHREPORTERPROTRUSIONLIMIDIVISIONCONTROL_HPP_
#define TESTDELTANOTCHREPORTERPROTRUSIONLIMIDIVISIONCONTROL_HPP_

/*
 * = Delta/Notch simulations with cellular protrusions and Notch-dependent cell division =
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
 * 
 * EMPTYLINE
 * 
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
 * This test is used as a control for Notch-dependent cell division, as suggested by Hunter et al. (2016)
 * Ginger L. Hunter, Zena Hadjivasiliou, Hope Bonin, Li He, Norbert Perrimon, Guillaume Charras, Buzz Baum; 
 * "Coordinated control of Notch/Delta signalling and cell cycle progression drives lateral inhibition-mediated tissue patterning". 
 * Development 1 July 2016; 143 (13): 2305–2310. doi: https://doi.org/10.1242/dev.134213
 * For comparison of results obtained with Notch-dependent cell division, we assume that cells have a fixed
 * probability of division at each time step. This is achieved through Chaste's base {{{FixedProbabilityCellCycleModel}}}
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

#include "FixedProbabilityCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DeltaNotchReporterProtrusionSrnModelLimi.hpp"
#include "DeltaNotchReporterProtrusionDifferentiationTrackingModifier.hpp"

class TestDeltaNotchReporterProtrusionLimiDivisionControl : public AbstractCellBasedTestSuite
{
public:

    void TestVertexBasedMonolayerWithDeltaNotch()
    {
        /* We include the next line because Vertex simulations cannot be run in parallel */
        EXIT_IF_PARALLEL;

        /* First we create a vertex mesh to run our simulations in. A large mesh is required
        for the patterning impact to appear.  
        */
        HoneycombVertexMeshGenerator generator(12, 12);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // // option to create a cylindrical vertex mesh for periodicity in the x-direction. 
        // CylindricalHoneycombVertexMeshGenerator generator(12, 12);
        // Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        /* We then create some cells, each with a subcellular reaction network model
         * {{{DeltaNotchReporterProtrusionSrnModelLimi}}}, which incorporates a Delta/Notch/Reporter ODE system 
         * with protrusional contacts and mutual inactivation,
         * as developed by Sprinzak et al (2011) and Baum et al (2016)
         * We choose to initialise the concentrations of Notch and Delta to random levels in [0, 1] 
         * in each cell. Similarly the concentration of Reporter in each cell is randomly initialized, but
         * to a very low value.
         * */
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            /* We begin by assigning the cell a CellCycle model which passes a fixed probability of division
            * at each time step to the cells. No cell phases are modelled. */
            FixedProbabilityCellCycleModel* p_transit_model = new FixedProbabilityCellCycleModel();
            
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
            
            CellPtr p_cell(new Cell(p_state, p_transit_model, p_srn_model, false, collection));
            TS_ASSERT_THROWS_NOTHING(p_cell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>());

            p_cell->SetCellProliferativeType(p_transit_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        /* Using the vertex mesh and cells, we create a cell-based population object, and specify which results to
         * output to file. */
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        /* We are now in a position to create and configure the cell-based simulation object, pass a force law to it,
         * and run the simulation. We can make the simulation run for longer to see more patterning by increasing the end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Attempt6AllDirL1Tl0.25");
        simulator.SetSamplingTimestepMultiple(20);
        simulator.SetEndTime(10.0);

        /* Then, we define the modifier class, which automatically updates the properties of the cells and passes it to the simulation.*/
        MAKE_PTR(DeltaNotchReporterProtrusionDifferentiationTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        /* This modifier assigns target areas to each cell, which are required by the {{{NagaiHondaForce}}}.
         */
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        p_growth_modifier->SetGrowthDuration(1.0);
        simulator.AddSimulationModifier(p_growth_modifier);
        simulator.Solve();
    }

};

#endif /*TESTDELTANOTCHREPORTERPROTRUSIONLIMIDIVISIONCONTROL_HPP_*/
