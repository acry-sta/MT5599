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

#include "DeltaNotchReporterSrnModelLimi.hpp"

DeltaNotchReporterSrnModelLimi::DeltaNotchReporterSrnModelLimi(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeSrnModel(3, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<DeltaNotchReporterSrnModelLimi, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<DeltaNotchReporterSrnModelLimi, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

DeltaNotchReporterSrnModelLimi::DeltaNotchReporterSrnModelLimi(const DeltaNotchReporterSrnModelLimi& rModel)
    : AbstractOdeSrnModel(rModel)
{
    /*
     * Set each member variable of the new SRN model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new SRN model's member variables
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new SRN model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new SRN model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     */

    assert(rModel.GetOdeSystem());
    SetOdeSystem(new DeltaNotchReporterOdeSystemLimi(rModel.GetOdeSystem()->rGetStateVariables()));
}

AbstractSrnModel* DeltaNotchReporterSrnModelLimi::CreateSrnModel()
{
    return new DeltaNotchReporterSrnModelLimi(*this);
}

void DeltaNotchReporterSrnModelLimi::SimulateToCurrentTime()
{
    // Custom behaviour
    UpdateDeltaNotchReporter();

    // Run the ODE simulation as needed
    AbstractOdeSrnModel::SimulateToCurrentTime();
}

void DeltaNotchReporterSrnModelLimi::Initialise()
{
    AbstractOdeSrnModel::Initialise(new DeltaNotchReporterOdeSystemLimi);
}

void DeltaNotchReporterSrnModelLimi::UpdateDeltaNotchReporter()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    double mean_notch = mpCell->GetCellData()->GetItem("mean notch");
    mpOdeSystem->SetParameter("Mean Notch", mean_notch);

    double mean_delta = mpCell->GetCellData()->GetItem("mean delta");
    mpOdeSystem->SetParameter("Mean Delta", mean_delta);
}

double DeltaNotchReporterSrnModelLimi::GetNotch()
{
    assert(mpOdeSystem != nullptr);
    double notch = mpOdeSystem->rGetStateVariables()[0];
    return notch;
}

double DeltaNotchReporterSrnModelLimi::GetDelta()
{
    assert(mpOdeSystem != nullptr);
    double delta = mpOdeSystem->rGetStateVariables()[1];
    return delta;
}

double DeltaNotchReporterSrnModelLimi::GetReporter()
{
    assert(mpOdeSystem != nullptr);
    double reporter = mpOdeSystem->rGetStateVariables()[2];
    return reporter;
}

double DeltaNotchReporterSrnModelLimi::GetMeanNeighbouringNotch()
{
    assert(mpOdeSystem != nullptr);
    double mean_neighbouring_notch = mpOdeSystem->GetParameter("Mean Notch");
    return mean_neighbouring_notch;
}

double DeltaNotchReporterSrnModelLimi::GetMeanNeighbouringDelta()
{
    assert(mpOdeSystem != nullptr);
    double mean_neighbouring_delta = mpOdeSystem->GetParameter("Mean Delta");
    return mean_neighbouring_delta;
}

void DeltaNotchReporterSrnModelLimi::OutputSrnModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DeltaNotchReporterSrnModelLimi)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(DeltaNotchReporterSrnModelLimi)
