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

#ifndef TESTREPORTERDEPENDENTCELLCYCLEMODELS_HPP_
#define TESTREPORTERDEPENDENTCELLCYCLEMODELS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <fstream>

#include "AbstractCellBasedTestSuite.hpp"
#include "ReporterDependentCellCycleModel.hpp"
#include "CellCycleTimesGenerator.hpp"
#include "CellLabel.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "Debug.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FileComparison.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"


//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestSimpleCellCycleModels : public AbstractCellBasedTestSuite
{
public:
void TestReporterDependentCellCycleModel()
    {
        TS_ASSERT_THROWS_NOTHING(ReporterDependentCellCycleModel cell_model3);

        ReporterDependentCellCycleModel* p_diff_model = new ReporterDependentCellCycleModel;
        ReporterDependentCellCycleModel* p_transit_model = new ReporterDependentCellCycleModel;

        TS_ASSERT_DELTA(p_transit_model->GetDissociationConstant(), 100.0, 1e-9);
        TS_ASSERT_DELTA(p_transit_model->GetHillCoefficient(), 5.0, 1e-9);
        TS_ASSERT_DELTA(p_transit_model->GetMinimumDivisionSimulationTime(), 2.0, 1e-9);

        // Change parameters for this model
        p_transit_model->SetDissociationConstant(200.0);
        p_transit_model->SetMinimumDivisionSimulationTime(1.0);
        TS_ASSERT_DELTA(p_transit_model->GetDissociationConstant(), 200.0, 1e-9);
        TS_ASSERT_DELTA(p_transit_model->GetMinimumDivisionSimulationTime(), 1.0, 1e-9);

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        CellPtr p_transit_cell(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell->SetCellProliferativeType(p_transit_type);
        p_transit_cell->InitialiseCellCycleModel();

        CellPtr p_diff_cell(new Cell(p_healthy_state, p_diff_model));
        p_diff_cell->SetCellProliferativeType(p_diff_type);
        p_diff_cell->InitialiseCellCycleModel();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, num_steps);

        p_transit_cell->GetCellData()->SetItem("reporter", 0.01);
        TS_ASSERT_EQUALS(p_transit_cell->ReadyToDivide(), false);

        p_diff_cell->GetCellData()->SetItem("reporter", 0.01);
        TS_ASSERT_EQUALS(p_diff_cell->ReadyToDivide(), false);

        for (unsigned i = 0; i < num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            p_transit_cell->GetCellData()->SetItem("reporter", 1000.0);
        }

        double time = p_simulation_time->GetTime();
        std::cout << "End time: " << time << std::endl;

        double dt = p_simulation_time->GetTimeStep();
        std::cout << "Time step: " << dt << std::endl; 
        // note that the time step dependence must be erased to pass the test since otherwise
        // the probability of division at each time step is 0.1, so this would need to be reevaluated        

        TS_ASSERT_EQUALS(p_transit_cell->ReadyToDivide(), true);

        TS_ASSERT_EQUALS(p_diff_cell->ReadyToDivide(), false);

        TS_ASSERT_DELTA(p_transit_model->GetAge(), p_simulation_time->GetTime(), 1e-9);
        TS_ASSERT_DELTA(p_diff_model->GetAge(), p_simulation_time->GetTime(), 1e-9);

        // Check that cell division correctly resets the cell cycle phase
        CellPtr p_transit_cell2 = p_transit_cell->Divide();
        ReporterDependentCellCycleModel* p_transit_model2 = static_cast<ReporterDependentCellCycleModel*>(p_transit_cell2->GetCellCycleModel());

        TS_ASSERT_EQUALS(p_transit_model2->ReadyToDivide(), false);
        TS_ASSERT_DELTA(p_transit_model2->GetDissociationConstant(), 200.0, 1e-9);
        TS_ASSERT_DELTA(p_transit_model2->GetMinimumDivisionSimulationTime(), 1.0, 1e-9);

        TS_ASSERT_EQUALS(p_transit_cell2->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>(), true);
        TS_ASSERT_EQUALS(p_transit_cell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>(), false);

    }
    void TestArchiveReporterDependentCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ReporterDependentCellCycleModel.arch";

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new ReporterDependentCellCycleModel;

            p_model->SetDimension(2);
            p_model->SetBirthTime(-1.0);
            static_cast<ReporterDependentCellCycleModel*>(p_model)->SetDissociationConstant(0.5);
            static_cast<ReporterDependentCellCycleModel*>(p_model)->SetHillCoefficient(0.1);
            static_cast<ReporterDependentCellCycleModel*>(p_model)->SetMinimumDivisionSimulationTime(0.1);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_model;

            delete p_model;
            SimulationTime::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractCellCycleModel* p_model2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_model2;

            // Check private data has been restored correctly
            TS_ASSERT_DELTA(p_model2->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(p_model2->GetAge(), 1.0, 1e-12);
            TS_ASSERT_EQUALS(p_model2->GetDimension(), 2u);
            TS_ASSERT_DELTA(static_cast<ReporterDependentCellCycleModel*>(p_model2)->GetDissociationConstant(), 0.5, 1e-9);
            TS_ASSERT_DELTA(static_cast<ReporterDependentCellCycleModel*>(p_model2)->GetMinimumDivisionSimulationTime(), 0.1, 1e-9);

            // Avoid memory leaks
            delete p_model2;
        }
    }
};

#endif /*TESTREPORTERDEPENDENTCELLCYCLEMODELS_HPP_*/
