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

#ifndef TESTDELTANOTCHREPORTER12x12PERIODICLIMIPREPATTERN_HPP_
#define TESTDELTANOTCHREPORTER12x12PERIODICLIMIPREPATTERN_HPP_

/*
 * = An example showing how to run Delta/Notch simulations =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to simulate a growing cell monolayer culture
 * into which a simple model of Delta/Notch signalling is incorporated. This model was developed
 * by Collier et al. ("Pattern formation by lateral inhibition with feedback: a mathematical
 * model of delta-notch intercellular signalling", J. Theor. Biol. 183:429-446) and comprises
 * two ODEs to describe the evolution in concentrations of Delta and Notch in each cell. The ODE
 * for Notch includes a reaction term that depends on the mean Delta concentration among neighbouring
 * cells. Thus in this simulation each cell needs to be able to access information about its
 * neighbours. We use the {{{CellData}}} class to facilitate this, and introduce a subclass
 * of {{{OffLatticeSimulation}}} called {{{DeltaNotchOffLatticeSimulation}}} to handle the updating
 * of {{{CellData}}} at each time step as cell neighbours change.
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
/*
 * The next header file defines a simple subcellular reaction network model that includes the functionality
 * for solving each cell's Delta/Notch signalling ODE system at each time step, using information about neighbouring
 * cells through the {{{CellData}}} class.
 */
#include "DeltaNotchReporterSrnModelLimi.hpp"
/*
 * The next header defines the simulation class modifier corresponding to the Delta-Notch SRN model.
 * This modifier leads to the {{{CellData}}} cell property being updated at each timestep to deal with Delta-Notch signalling.
 */
#include "DeltaNotchReporterTrackingModifier.hpp"

/* Having included all the necessary header files, we proceed by defining the test class.
 */
class TestDeltaNotchReporter12x12PeriodicLimiPrepattern : public AbstractCellBasedTestSuite
{
public:

    /*
     * EMPTYLINE
     *
     * == Test 1: a periodic vertex-based monolayer with Delta/Notch signalling ==
     *
     * EMPTYLINE
     *
     * In the first test, we demonstrate how to simulate a monolayer that incorporates
     * Delta/Notch signalling, using a vertex-based approach.
     */
    void TestVertexBasedMonolayerWithDeltaNotch()
    {
        /* We include the next line because Vertex simulations cannot be run in parallel */
        EXIT_IF_PARALLEL;

        /* First we create a 7x7 cylindrical vertex mesh for periodicity in the x-direction. */
        CylindricalHoneycombVertexMeshGenerator generator(12, 12);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        /* We then create some cells, each with a cell-cycle model, {{{UniformG1GenerationalCellCycleModel}}} and a subcellular reaction network model
         * {{{DeltaNotchSrnModel}}}, which
         * incorporates a Delta/Notch ODE system, here we use the hard coded initial conditions of 1.0 and 1.0.
         * In this example we choose to make each cell differentiated,
         * so that no cell division occurs. */
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();
            p_cc_model->SetDimension(2);

            /* For prepatterns, only stripes appear to affect the final pattern. Due to the 12x12
            form of the grid, strips can be made by taking any value modulo 3, 4 or 6 (divisors of 12).
            Thickness of the prepattern stripes is determined by the value of the remainder.*/
            std::vector<double> initial_conditions;
            // mod 2 striped prepattern           
            // if ((elem_index%2)<1){
            //     initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            //   }
            //   else{
            //     initial_conditions.push_back(9.0 + 1.0*RandomNumberGenerator::Instance()->ranf());
            //   }
            // mod 3 striped prepattern
              //  if ((elem_index%3)<2){
              //     initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
              //   }
              //   else{
              //     initial_conditions.push_back(9.0 + 1.0*RandomNumberGenerator::Instance()->ranf());
              //   }
            // // mod 4 striped prepattern
            //    if ((elem_index%4)<3){
            //       initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            //     }
            //     else{
            //       initial_conditions.push_back(9.0 + 1.0*RandomNumberGenerator::Instance()->ranf());
            //     }
            // // mod 6 striped prepattern
            //   if ((elem_index%6)<4){
            //     initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            //   }
            //   else{
            //     initial_conditions.push_back(9.0 + 1.0*RandomNumberGenerator::Instance()->ranf());
            //   }
            // // offset 2x2 patches prepattern, lines thickness 1, 2
            // if ( ((elem_index%4)<2)&&((int)(elem_index/12)<2) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else if ( ((elem_index%4)<2)&&((int)(elem_index/12)<8) &&((int)(elem_index/12)>5) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else if ( ((elem_index%4)>1)&&((int)(elem_index/12)<5) &&((int)(elem_index/12)>2) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else if ( ((elem_index%4)>1)&&((int)(elem_index/12)<11) &&((int)(elem_index/12)>8) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else{
            //   initial_conditions.push_back(9.0 + 1.0*RandomNumberGenerator::Instance()->ranf());
            // }
            // // offset 2x2 patches prepattern, lines thickness 1
            // if ( (((elem_index + 1)%3)<2)&&((int)(elem_index/12)<2) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else if ( (((elem_index + 1)%3)<2)&&((int)(elem_index/12)<8) &&((int)(elem_index/12)>5) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else if ( ((elem_index%3)>0)&&((int)(elem_index/12)<5) &&((int)(elem_index/12)>2) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else if ( ((elem_index%3)>0)&&((int)(elem_index/12)<11) &&((int)(elem_index/12)>8) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else{
            //   initial_conditions.push_back(9.0 + 1.0*RandomNumberGenerator::Instance()->ranf());
            // }
            // // regular 2x2 patches prepattern, thickness 2 line between
            // if ( ((elem_index%4)<2)&&((int)(elem_index/12)<2) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else if ( ((elem_index%4)<2)&&((int)(elem_index/12)<6) &&((int)(elem_index/12)>3) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else if ( ((elem_index%4)<2)&&((int)(elem_index/12)<10) &&((int)(elem_index/12)>7) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else{
            //   initial_conditions.push_back(9.0 + 1.0*RandomNumberGenerator::Instance()->ranf());
            // }
            // // regular 2x2 patches prepattern, thickness 1 line between
            // if ( (((elem_index)%3)<2)&&((int)(elem_index/12)<2) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else if ( ((elem_index%3)<2)&&((int)(elem_index/12)<8) &&((int)(elem_index/12)>5) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else if ( ((elem_index%3)<2)&&((int)(elem_index/12)<5) &&((int)(elem_index/12)>2) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else if ( ((elem_index%3)<2)&&((int)(elem_index/12)<11) &&((int)(elem_index/12)>8) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else{
            //   initial_conditions.push_back(9.0 + 1.0*RandomNumberGenerator::Instance()->ranf());
            // }
            // // regular 3x3 patches prepattern, thickness 1 line between
            // if ( ((elem_index%4)<3)&&((int)(elem_index/12)<3) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else if ( ((elem_index%4)<3)&&((int)(elem_index/12)<7) &&((int)(elem_index/12)>3) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else if ( ((elem_index%4)<3)&&((int)(elem_index/12)<11) &&((int)(elem_index/12)>7) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else{
            //   initial_conditions.push_back(9.0 + 1.0*RandomNumberGenerator::Instance()->ranf());
            // }
            // // regular 4x4 patches prepattern, thickness 2 line between
            // if ( ((elem_index%6)<4)&&((int)(elem_index/12)<4) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else if ( ((elem_index%6)<4)&&((int)(elem_index/12)<10) &&((int)(elem_index/12)>5) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else{
            //   initial_conditions.push_back(9.0 + 1.0*RandomNumberGenerator::Instance()->ranf());
            // }
            // // offset 4x4 patches prepattern, thickness 2 line between
            // if ( ((elem_index%6)<4)&&((int)(elem_index/12)<4) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else if ( ((elem_index%6)>1)&&((int)(elem_index/12)<10) &&((int)(elem_index/12)>5) ){
            //   initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // }
            // else{
            //   initial_conditions.push_back(9.0 + 1.0*RandomNumberGenerator::Instance()->ranf());
            // }
            // regular 5x5 patches prepattern, thickness 1 line between
            if ( ((elem_index%6)<5)&&((int)(elem_index/12)<5) ){
              initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            }
            else if ( ((elem_index%6)<5)&&((int)(elem_index/12)<11) &&((int)(elem_index/12)>5) ){
              initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            }
            else{
              initial_conditions.push_back(9.0 + 1.0*RandomNumberGenerator::Instance()->ranf());
            }
            initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf()); 
            initial_conditions.push_back(0.1*RandomNumberGenerator::Instance()->ranf());
            DeltaNotchReporterSrnModelLimi* p_srn_model = new DeltaNotchReporterSrnModelLimi();
            p_srn_model->SetInitialConditions(initial_conditions);

            CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));
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
        simulator.SetOutputDirectory("12x12Periodic5x5RegularPatchesLt1Notch");
        simulator.SetSamplingTimestepMultiple(40);
        simulator.SetEndTime(20.0);

        /* Then, we define the modifier class, which automatically updates the values of Delta and Notch within the cells in {{{CellData}}} and passes it to the simulation.*/
        MAKE_PTR(DeltaNotchReporterTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        /* This modifier assigns target areas to each cell, which are required by the {{{NagaiHondaForce}}}.
         */
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);
        simulator.Solve();
    }

    /*
     * EMPTYLINE
     *
     * To visualize the results, use Paraview. See the UserTutorials/VisualizingWithParaview tutorial for more information.
     *
     * Load the file {{{/tmp/$USER/testoutput/TestVertexBasedMonolayerWithDeltaNotch/results_from_time_0/results.pvd}}}.
     *
     * EMPTYLINE
     *
     * == Test 2 - a node-based monolayer with Delta/Notch signalling ==
     *
     * EMPTYLINE
     *
     * In the next test we run a similar simulation as before, but this time with node-based
     * 'overlapping spheres' model.
     */
};

#endif /*TESTDELTANOTCHREPORTER12x12PERIODICLIMIPREPATTERN_HPP_*/
