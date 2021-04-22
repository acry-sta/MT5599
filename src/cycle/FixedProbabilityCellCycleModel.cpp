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

#include "FixedProbabilityCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

FixedProbabilityCellCycleModel::FixedProbabilityCellCycleModel()
    : AbstractCellCycleModel(),
      mDivisionProbability(0.05),
      mMinimumDivisionSimulationTime(2.0)
{
}

FixedProbabilityCellCycleModel::FixedProbabilityCellCycleModel(const FixedProbabilityCellCycleModel& rModel)
   : AbstractCellCycleModel(rModel),
     mDivisionProbability(rModel.mDivisionProbability),
     mMinimumDivisionSimulationTime(rModel.mMinimumDivisionSimulationTime)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     */
}

bool FixedProbabilityCellCycleModel::ReadyToDivide()
{
    assert(mpCell != nullptr);

    if (!mReadyToDivide)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        double time = p_simulation_time->GetTime();
        if (time > mMinimumDivisionSimulationTime) // cells only allowed to divide after min sim time
        {
            double dt = p_simulation_time->GetTimeStep();
            if (!(mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()))
            {
                RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
                if (p_gen->ranf() < mDivisionProbability*dt) // generate a random number to compare with fixed probability
                {
                    mReadyToDivide = true;
                }
            }
        }
    }
    return mReadyToDivide;
}

void FixedProbabilityCellCycleModel::ResetForDivision()
{
    AbstractCellCycleModel::ResetForDivision();
    boost::shared_ptr<AbstractCellProperty> p_diff_type =
        mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
    mpCell->SetCellProliferativeType(p_diff_type);
    mpCell->SetCellProliferativeType(p_diff_type);
}

void FixedProbabilityCellCycleModel::InitialiseDaughterCell()
{
    boost::shared_ptr<AbstractCellProperty> p_diff_type =
        mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
}

AbstractCellCycleModel* FixedProbabilityCellCycleModel::CreateCellCycleModel()
{
    return new FixedProbabilityCellCycleModel(*this);
}

void FixedProbabilityCellCycleModel::SetDivisionProbability(double divisionProbability)
{
    mDivisionProbability = divisionProbability;
}

double FixedProbabilityCellCycleModel::GetDivisionProbability()
{
    return mDivisionProbability;
}

void FixedProbabilityCellCycleModel::SetMinimumDivisionSimulationTime(double minimumDivisionSimulationTime)
{
    mMinimumDivisionSimulationTime = minimumDivisionSimulationTime;
}

double FixedProbabilityCellCycleModel::GetMinimumDivisionSimulationTime()
{
    return mMinimumDivisionSimulationTime;
}

double FixedProbabilityCellCycleModel::GetAverageTransitCellCycleTime()
{
    return 1.0/mDivisionProbability; // probably not correct
}

double FixedProbabilityCellCycleModel::GetAverageStemCellCycleTime()
{
    return 1.0/mDivisionProbability; // stem cells are not included, this is redundant
}

void FixedProbabilityCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<DivisionProbability>" << mDivisionProbability << "</DivisionProbability>\n";
    *rParamsFile << "\t\t\t<MinimumDivisionSimulationTime>" << mMinimumDivisionSimulationTime << "</MinimumDivisionSimulationTime>\n";

    // Call method on direct parent class
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(FixedProbabilityCellCycleModel)
