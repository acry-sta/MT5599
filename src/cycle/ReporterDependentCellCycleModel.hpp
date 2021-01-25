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

#ifndef REPORTERDEPENDENTCELLCYCLEMODEL_HPP_
#define REPORTERDEPENDENTCELLCYCLEMODEL_HPP_

#include "AbstractCellCycleModel.hpp"

/**
 * Simple cell-cycle model where mature non-differentiated cells have a specified probability of
 * dividing per hour.
 *
 * The class includes two parameters: the first, mDivisionProbability, defines the probability
 * of dividing per hour; the second, mMinimumDivisionAge, defines a minimum age at which cells
 * may divide.
 */
class ReporterDependentCellCycleModel : public AbstractCellCycleModel
{
private:

    friend class boost::serialization::access;

    /**
     * Boost Serialization method for archiving/checkpointing
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        archive & mDissociationConstant;
        archive & mHillCoefficient;
        archive & mMinimumDivisionSimulationTime;
    }

protected:
    /**
     * Dissociation constant K in denominator of division probability.
     * This determines the range of values for which division is likely (Hunter et al.)
     * Defaults to 100.
     */
    double mDissociationConstant;

    /**
     * Hill coefficient q used in determining the probability of division.
     * The larger q, the smaller R must be to induce division. Defaults to 5.
     */
    double mHillCoefficient;

    /**
     * Minimum time simulation must run before cells may divide.
     * This is required for correct setup of Delta-Notch interactions to take place. 
     * Defaults to 2 hours.
     */
    double mMinimumDivisionSimulationTime;

    /**
     * Protected copy-constructor for use by CreateCellCycleModel.
     *
     * The only way for external code to create a copy of a cell cycle model
     * is by calling that method, to ensure that a model of the correct subclass
     * is created. This copy-constructor helps subclasses to ensure that all
     * member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a
     * daughter cell upon cell division. Note that the parent cell cycle model
     * will have had ResetForDivision() called just before CreateCellCycleModel()
     * is called, so performing an exact copy of the parent is suitable behaviour.
     * Any daughter-cell-specific initialisation can be done in InitialiseDaughterCell().
     *
     * @param rModel the cell cycle model to copy.
     */
    ReporterDependentCellCycleModel(const ReporterDependentCellCycleModel& rModel);

public:

    /**
     * Constructor.
     */
    ReporterDependentCellCycleModel();

    /**
     * Overridden ReadyToDivide() method.
     * 
     * Division is dependent on Reporter expression in the cell, and can only occur
     * once the simulation has run a certain minimum time. If the simulation time is 
     * greater than mMinimumDivisionSimulationTime, then we draw a uniform
     * random number r ~ U[0,1] and retrieve the cell's reporter expression R.
     * If r < (R^q)/(K^q + R^q)*dt, where dt is the
     * simulation time step, then the cell is ready to divide and we return true.
     * Otherwise, the cell is not yet ready to divide and we return false.
     *
     * @return whether the cell is ready to divide.
     */
    bool ReadyToDivide();

    /**
     * Overriden InitialiseDaughterCell method, this allows us to mark any cells that have divided
     * as differentiated cells that no longer divide.
     */
    void InitialiseDaughterCell();

    /**
     * Overriden ResetForDivision method, this allows us to mark any cells that have divided
     * as differentiated cells that no longer divide.
     */
    void ResetForDivision();

    /**
     * Overridden builder method to create new instances of
     * the cell-cycle model.
     *
     * @return new cell-cycle model
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Set the value of mDivisionProbability.
     *
     * @param dissociationConstant the new value of mDivisionProbability
     */
    void SetDissociationConstant(double dissociationConstant);

    /**
     * Get mDissociationConstant.
     *
     * @return mDissociationConstant
     */
    double GetDissociationConstant();

    /**
     * Set the value of mHillCoefficient.
     *
     * @param hillCoefficient the new value of mHillCoefficient
     */
    void SetHillCoefficient(double hillCoefficient);

    /**
     * Get mHillCoefficient.
     *
     * @return mHillCoefficient
     */
    double GetHillCoefficient();

    /**
     * Set the value of mMinimumDivisionSimulationTime.
     *
     * @param minimumDivisionSimulationTime the new value of mMinimumDivisionSimulationTime
     */
    void SetMinimumDivisionSimulationTime(double minimumDivisionSimulationTime);

    /**
     * Get mMinimumDivisionSimulationTime.
     *
     * @return mMinimumDivisionSimulationTime
     */
    double GetMinimumDivisionSimulationTime();

    /**
     * Overridden GetAverageTransitCellCycleTime() method.
     *
     * @return the average cell cycle time for cells of transit proliferative type
     */
    double GetAverageTransitCellCycleTime();

    /**
     * Overridden GetAverageStemCellCycleTime() method.
     *
     * @return the average cell cycle time for cells of stem proliferative type
     */
    double GetAverageStemCellCycleTime();

    /**
     * Overridden OutputCellCycleModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(ReporterDependentCellCycleModel)

#endif // REPORTERDEPENDENTCELLCYCLEMODEL_HPP_
