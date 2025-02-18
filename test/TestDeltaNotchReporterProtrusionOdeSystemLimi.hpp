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

#ifndef TESTDELTANOTCHREPORTERPROTRUSIONODESYSTEMLIMI_HPP_
#define TESTDELTANOTCHREPORTERPROTRUSIONODESYSTEMLIMI_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <vector>
#include <iostream>

#include "OutputFileHandler.hpp"
#include "DeltaNotchReporterProtrusionOdeSystemLimi.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "RungeKuttaFehlbergIvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "CvodeAdaptor.hpp"
#include "Timer.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

// Tests system of Delta/Notch/Reporter with mutual inactivation and protrusions.
class TestDeltaNotchReporterProtrusionOdeSystemLimi : public CxxTest::TestSuite
{
public:

    void TestDeltaNotchReporterProtrusionOdeSystemSetup()
    {
#ifdef CHASTE_CVODE
        DeltaNotchReporterProtrusionOdeSystemLimi ode_system;

        double h_value = 1.0; // Was 0.0001 for some utterly bizarre reason
        CvodeAdaptor cvode_solver;
        OdeSolution solutions;

        std::vector<double> initial_conditions = ode_system.GetInitialConditions();
        unsigned number_of_parameters = ode_system.GetNumberOfParameters();

        TS_ASSERT_DELTA(initial_conditions[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(initial_conditions[1], 1.0, 1e-6);
        TS_ASSERT_DELTA(initial_conditions[2], 1.0, 1e-6);
        TS_ASSERT_DELTA(number_of_parameters, 9, 1e-6);

        Timer::Reset();
        solutions = cvode_solver.Solve(&ode_system, initial_conditions, 0.0, 100.0, h_value, h_value);
        Timer::Print("1. Cvode for 100 hours");

        // Test solutions are OK for a small time increase...
        int end = solutions.rGetSolutions().size() - 1;
        std::cout << "End time is:" << end << std::endl;
        // Tests the simulation is ending at the right time...(going into S phase at 7.8 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 100, 1e-2);

        // Decent results
        // TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 0.9615, 1e-4);
        // TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 0.0107, 1e-4);
#else
        std::cout << "CVODE is not enabled. " << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_CVODE
    }

    void TestArchiving()
    {
#ifdef CHASTE_CVODE
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "delta_notch_ode.arch";

        {
            std::vector<double> state_variables;
            state_variables.push_back(3.0);
            state_variables.push_back(4.0);
            state_variables.push_back(5.0);

            DeltaNotchReporterProtrusionOdeSystemLimi ode_system(state_variables);

            ode_system.SetDefaultInitialCondition(1, 3.25);
            ode_system.SetDefaultInitialCondition(2, 4.0);

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();

            // These are the initial conditions hard-coded in the constructor.
            TS_ASSERT_EQUALS(initial_conditions.size(), 3u);
            TS_ASSERT_DELTA(initial_conditions[0], 1.0, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[1], 3.25, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[2], 4.0, 1e-6);

            ode_system.SetParameter(0, 1.0);
            ode_system.SetParameter(1, 10.0);

            double mean_notch = ode_system.GetParameter("Mean Notch");
            double mean_delta = ode_system.GetParameter("Mean Delta");
            std::cout << "Parameter Mean Notch is " << mean_notch <<std::endl;
            std::cout << "Parameter Mean Delta is " << mean_delta <<std::endl;

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive ODE system
            AbstractOdeSystem* const p_const_ode_system = &ode_system;
            output_arch << p_const_ode_system;
        }

        {
            AbstractOdeSystem* p_ode_system;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_ode_system;

            // Check that archiving worked correctly
            std::vector<double> initial_conditions = p_ode_system->GetInitialConditions();

            TS_ASSERT_EQUALS(initial_conditions.size(), 3u);
            TS_ASSERT_DELTA(initial_conditions[0], 1.0, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[1], 1.0, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[2], 1.0, 1e-6);

            double mean_notch = p_ode_system->GetParameter("Mean Notch");
            double mean_delta = p_ode_system->GetParameter("Mean Delta");

            double protrusion_length = p_ode_system->GetParameter("Protrusion Length");
            double protrusion_angular_opening = p_ode_system->GetParameter("Protrusion Angular Opening");

            TS_ASSERT_DELTA(mean_notch, 1.0, 1e-3);
            TS_ASSERT_DELTA(mean_delta, 10.0, 1e-3);

            // check that this matches the current parameters
            TS_ASSERT_DELTA(protrusion_length, 2.0, 1e-3);
            TS_ASSERT_DELTA(protrusion_angular_opening, 0.5*M_PI, 1e-3);

            // Get state variables
            double var1 = p_ode_system->GetStateVariable(0);
            double var2 = p_ode_system->GetStateVariable(1);
            double var3 = p_ode_system->GetStateVariable(2);

            TS_ASSERT_DELTA(var1, 3.0, 1e-3);
            TS_ASSERT_DELTA(var2, 4.0, 1e-3);
            TS_ASSERT_DELTA(var3, 5.0, 1e-3);

            // Tidy up
            delete p_ode_system;
        }
#endif //CHASTE_CVODE
    }

    void TestSetStateVariables()
    {
#ifdef CHASTE_CVODE

        std::vector<double> state_vars;
        state_vars.push_back(0.0);
        state_vars.push_back(1.0);
        state_vars.push_back(1.0);
        DeltaNotchReporterProtrusionOdeSystemLimi ode_system(state_vars);

        TS_ASSERT_EQUALS(ode_system.GetStateVariable(0),0.0);
        TS_ASSERT_EQUALS(ode_system.GetStateVariable(1),1.0);
        TS_ASSERT_EQUALS(ode_system.GetStateVariable(2),1.0);

#endif //CHASTE_CVODE
   }
};

#endif /*TESTDELTANOTCHREPORTERPROTRUSIONODESYSTEMLIMI_HPP_*/
