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

#include "DeltaNotchReporterProtrusionOdeSystemLimi.hpp"
#include "CellwiseOdeSystemInformation.hpp"

DeltaNotchReporterProtrusionOdeSystemLimi::DeltaNotchReporterProtrusionOdeSystemLimi(std::vector<double> stateVariables)
    : AbstractOdeSystem(3)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<DeltaNotchReporterProtrusionOdeSystemLimi>);

    /**
     * The state variables are as follows:
     *
     * 0 - Notch concentration for this cell
     * 1 - Delta concentration for this cell
     *
     * We store the last state variable so that it can be written
     * to file at each time step alongside the others, and visualized.
     */

    SetDefaultInitialCondition(0, 1.0); // soon overwritten
    SetDefaultInitialCondition(1, 1.0); // soon overwritten
    SetDefaultInitialCondition(2, 1.0); // soon overwritten

    this->mParameters.push_back(0.5); // mean notch initial value
    this->mParameters.push_back(0.5); // mean delta initial value

    this->mParameters.push_back(0.5); // protrusion-mediated mean notch initial value
    this->mParameters.push_back(0.5); // protrusion-mediated mean delta initial value


    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

DeltaNotchReporterProtrusionOdeSystemLimi::~DeltaNotchReporterProtrusionOdeSystemLimi()
{
}

void DeltaNotchReporterProtrusionOdeSystemLimi::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    double notch = rY[0];
    double delta = rY[1];
    double reporter = rY[2];
    double mean_notch = this->mParameters[0]; // Shorthand for "this->mParameter("Mean Notch");"
    double mean_delta = this->mParameters[1]; // Shorthand for "this->mParameter("Mean Delta");"
    double protrusion_notch = this->mParameters[2]; // Shorthand for "this->mParameter("Mean Notch");"
    double protrusion_delta = this->mParameters[3]; // Shorthand for "this->mParameter("Mean Delta");"

    // we allow for weighting of the different methods of delta and notch arriving at the cell
    // to account for possible differences in efficiency of signalling mechanisms; we denote wn
    // the weight for neighboring cell signals, and wp the weight for protrusion-mediated signalling. 
    double weighted_delta_in = 1.0 * mean_notch + 0.1 * protrusion_notch
    double weighted_notch_in = 1.0 * mean_delta + 0.1 * protrusion_delta 
    
    // total amount of bound delta that leads 
    // weights qn and qp are bounded by wn and wp. qn/wn is the proportion of Delta molecules bound
    // in trans that leads to a Notch signal in the receiving cell
    double weighted_delta_out = 0.01 * mean_delta + 0.025 * protrusion_delta

    // The next two lines define the ODE system by Baum et al. (2016)
    rDY[0] = 100.0 - notch - notch*weighted_delta_in/2.0 - notch*delta/0.5;  // d[Notch]/dt
    rDY[1] = 500.0/(1.0 + reporter*reporter) - delta - delta*weighted_notch_in/2.0 - notch*delta/0.5;   // d[Delta]/dt
    rDY[2] = 300000.0*notch*weighted_delta_out*notch*weighted_delta_out/(1000000.0 + notch*weighted_delta_out*notch*weighted_delta_out) - reporter;   // d[Reporter]/dt
}

template<>
void CellwiseOdeSystemInformation<DeltaNotchReporterProtrusionOdeSystemLimi>::Initialise()
{
    this->mVariableNames.push_back("Notch");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("Delta");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("Reporter");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mParameterNames.push_back("Mean Notch");
    this->mParameterUnits.push_back("non-dim");
    
    this->mParameterNames.push_back("Mean Delta");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Protrusion Notch");
    this->mParameterUnits.push_back("non-dim");
    
    this->mParameterNames.push_back("Protrusion Delta");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Protrusion Notch");
    this->mParameterUnits.push_back("non-dim");
    
    this->mParameterNames.push_back("Protrusion Delta");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DeltaNotchReporterProtrusionOdeSystemLimi)
