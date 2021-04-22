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

#include "ReporterDependentCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

ReporterDependentCellCycleModel::ReporterDependentCellCycleModel()
    : AbstractCellCycleModel(),
      mDissociationConstant(100.0),
      mHillCoefficient(5.0),
      mMinimumDivisionSimulationTime(2.0)
{
}

ReporterDependentCellCycleModel::ReporterDependentCellCycleModel(const ReporterDependentCellCycleModel& rModel)
   : AbstractCellCycleModel(rModel),
     mDissociationConstant(rModel.mDissociationConstant),
     mHillCoefficient(rModel.mHillCoefficient),
     mMinimumDivisionSimulationTime(rModel.mMinimumDivisionSimulationTime)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     */
}

bool ReporterDependentCellCycleModel::ReadyToDivide()
{
    assert(mpCell != nullptr);

    if (!mReadyToDivide)
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        double time = p_simulation_time->GetTime();
        if (time > mMinimumDivisionSimulationTime)
        {
            double reporter = mpCell->GetCellData()->GetItem("reporter");
            double dt = p_simulation_time->GetTimeStep();
            if (!(mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()))
            {
                RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
                double division_probability;
                // calculate division probability based on Hunter et al (2016)
                division_probability = std::pow(reporter, mHillCoefficient)/(std::pow(mDissociationConstant, mHillCoefficient) + std::pow(reporter, mHillCoefficient));
                if (p_gen->ranf() < division_probability*dt)
                {
                    mReadyToDivide = true;
                }
            }
        }
    }
    return mReadyToDivide;
}

void ReporterDependentCellCycleModel::ResetForDivision()
{
    AbstractCellCycleModel::ResetForDivision();
    boost::shared_ptr<AbstractCellProperty> p_diff_type =
        mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
    mpCell->SetCellProliferativeType(p_diff_type);
    mpCell->SetCellProliferativeType(p_diff_type);
}

void ReporterDependentCellCycleModel::InitialiseDaughterCell()
{
    boost::shared_ptr<AbstractCellProperty> p_diff_type =
        mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
}

AbstractCellCycleModel* ReporterDependentCellCycleModel::CreateCellCycleModel()
{
    return new ReporterDependentCellCycleModel(*this);
}

void ReporterDependentCellCycleModel::SetDissociationConstant(double dissociationConstant)
{
    mDissociationConstant = dissociationConstant;
}

double ReporterDependentCellCycleModel::GetDissociationConstant()
{
    return mDissociationConstant;
}

void ReporterDependentCellCycleModel::SetHillCoefficient(double hillCoefficient)
{
    mHillCoefficient = hillCoefficient;
}

double ReporterDependentCellCycleModel::GetHillCoefficient()
{
    return mHillCoefficient;
}

void ReporterDependentCellCycleModel::SetMinimumDivisionSimulationTime(double minimumDivisionSimulationTime)
{
    mMinimumDivisionSimulationTime = minimumDivisionSimulationTime;
}

double ReporterDependentCellCycleModel::GetMinimumDivisionSimulationTime()
{
    return mMinimumDivisionSimulationTime;
}

double ReporterDependentCellCycleModel::GetAverageTransitCellCycleTime()
{
    return 0.0; // this is not determined
}

double ReporterDependentCellCycleModel::GetAverageStemCellCycleTime()
{
    return 0.0; // stem cells are not modelled so this is redundant
}

void ReporterDependentCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<DissociationConstant>" << mDissociationConstant << "</DissociationConstant>\n";
    *rParamsFile << "\t\t\t<HillCoefficient>" << mHillCoefficient << "</HillCoefficient>\n";
    *rParamsFile << "\t\t\t<MinimumDivisionSimulationTime>" << mMinimumDivisionSimulationTime << "</MinimumDivisionSimulationTime>\n";

    // Call method on direct parent class
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ReporterDependentCellCycleModel)
