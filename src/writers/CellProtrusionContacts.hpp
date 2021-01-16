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

#ifndef CELLPROTRUSIONCONTACTS_HPP_
#define CELLPROTRUSIONCONTACTS_HPP_

#include <boost/shared_ptr.hpp>
#include <set>

#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include "Exception.hpp"

/**
 * CellProtrusionContacts cell property class.
 *
 * This cell property allows each cell to store the indices of the cells it is in protrusion-mediated contact with,
 * based on a given length, angle, and angle of opening of a protrusion arc, as well as an activation threshold to
 * control for minimum overlap of protrusions, to be used in Delta Notch models. Protrusion details will be given in
 * the SRN model. Other classes may interrogate or modify the values stored in this class.
 *
 * Within the Cell constructor, an empty CellProtrusionContacts object is created and passed to the Cell
 * (unless there is already a CellProtrusionContacts object present in mCellPropertyCollection).
 */
class CellProtrusionContacts : public AbstractCellProperty
{
private:

    /**
     * The set of protrusion contact indices.
     */
    std::set<unsigned> mCellProtrusionContacts; 

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellProperty>(*this);
        archive & mCellProtrusionContacts;
    }

public:
   
    /**
     * We need the empty virtual destructor in this class to ensure Boost
     * serialization works correctly with static libraries.
     */
    virtual ~CellProtrusionContacts();

    /**
     * This assigns the set of protrusion contacts. 
     */
    void SetCellProtrusionContacts(const std::set<unsigned> indices);

    /**
     * @return data.
     */
    std::set<unsigned> GetCellProtrusionContacts() const;

    /**
     * @return number of contacts overall
     */
    unsigned GetNumCellProtrusionContacts() const;
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellProtrusionContacts)

#endif /* CELLPROTRUSIONCONTACTS_HPP_ */
