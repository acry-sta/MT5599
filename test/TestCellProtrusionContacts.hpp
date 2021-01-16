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

#ifndef TESTCELLPROTRUSIONCONTACTS_HPP_
#define TESTCELLPROTRUSIONCONTACTS_HPP_

#include <cxxtest/TestSuite.h>

#include <fstream>
#include <iostream>
#include <boost/shared_ptr.hpp>

#include "Cell.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "CellPropertyRegistry.hpp"
#include "CellProtrusionContacts.hpp"
#include "SmartPointers.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"

#define NEW_PROP(type, name) boost::shared_ptr<AbstractCellProperty> name(new type)

class TestCellProtrusionContacts : public AbstractCellBasedTestSuite
{
public:

    void TestPropertyCollection()
    {
        CellPropertyCollection collection;
        TS_ASSERT_EQUALS(collection.GetSize(), 0u);

        // Test that the CellPropertyRegistry assigned to the CellPropertyCollection defaults to CellPropertyRegistry::Instance().
        TS_ASSERT_EQUALS(collection.GetCellPropertyRegistry(),CellPropertyRegistry::Instance());
        MAKE_PTR(WildTypeCellMutationState, p_wt_mutation);
        MAKE_PTR(CellProtrusionContacts, p_protrusion_contacts);

        TS_ASSERT_EQUALS(p_wt_mutation->GetIdentifier(), "WildTypeCellMutationState");

        collection.AddProperty(p_protrusion_contacts);
        collection.AddProperty(p_wt_mutation);

        // Test we can't add the same *object* twice
        TS_ASSERT_THROWS_THIS(collection.AddProperty(p_protrusion_contacts),
                              "That property object is already in the collection.");
        MAKE_PTR(WildTypeCellMutationState, p_wt_mutation_2);
        collection.AddProperty(p_wt_mutation_2);
        collection.RemoveProperty(p_wt_mutation_2);

        
        // Check the contents
        TS_ASSERT_EQUALS(collection.GetSize(), 2u);
        // ...by object
        TS_ASSERT_EQUALS(collection.HasProperty(p_wt_mutation), true);
        TS_ASSERT_EQUALS(collection.HasProperty(p_protrusion_contacts), true);
        MAKE_PTR(ApcOneHitCellMutationState, p_apc1_mutation);
        TS_ASSERT_EQUALS(collection.HasProperty(p_apc1_mutation), false);
        // ...by type
        TS_ASSERT_EQUALS(collection.HasProperty<WildTypeCellMutationState>(), true);
        TS_ASSERT_EQUALS(collection.HasProperty<CellProtrusionContacts>(), true);
        // ..by subclass
        TS_ASSERT_EQUALS(collection.HasPropertyType<AbstractCellProperty>(), true);
        TS_ASSERT_EQUALS(collection.HasPropertyType<AbstractCellMutationState>(), true);
        //TS_ASSERT_EQUALS(!collection.HasProperty<AbstractCellMutationState>(), false); <-- This won't compile (yet)
        // ..by iteration

        // Remove property
        collection.RemoveProperty<WildTypeCellMutationState>();
        TS_ASSERT_EQUALS(collection.HasProperty<WildTypeCellMutationState>(), false);
        collection.RemoveProperty(p_protrusion_contacts);
        TS_ASSERT_EQUALS(collection.HasProperty<CellProtrusionContacts>(), false);
        TS_ASSERT_THROWS_THIS(collection.RemoveProperty<WildTypeCellMutationState>(),
                              "Collection does not contain the given property type.");
        TS_ASSERT_THROWS_THIS(collection.RemoveProperty(p_protrusion_contacts),
                              "Collection does not contain the given property.");

        TS_ASSERT_EQUALS(collection.GetSize(), 0u);

        collection.AddProperty(p_wt_mutation);
        collection.AddProperty(p_protrusion_contacts);
        CellPropertyCollection properties = collection.GetPropertiesType<AbstractCellProperty>();
        TS_ASSERT_EQUALS(properties.GetSize(), 2u);

        CellPropertyCollection mutations = collection.GetPropertiesType<AbstractCellMutationState>();
        TS_ASSERT_EQUALS(mutations.GetSize(), 1u);
        CellPropertyCollection::Iterator it = mutations.Begin();
        TS_ASSERT_EQUALS(collection.HasProperty(*it), true);
        TS_ASSERT_EQUALS((*it)->IsSubType<AbstractCellMutationState>(), true);
        TS_ASSERT( (*it)->IsType<WildTypeCellMutationState>() );
    }
        
    void TestCellData()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(25, 2);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        // boost::shared_ptr<AbstractCellProperty> p_protrusion_contacts(CellPropertyRegistry::Instance()->Get<CellProtrusionContacts>());
        MAKE_PTR(CellProtrusionContacts, p_protrusion_contacts);


        FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
        CellPropertyCollection collection;
        collection.AddProperty(p_protrusion_contacts);
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
        p_cell->SetCellProliferativeType(p_type);
        p_cell->InitialiseCellCycleModel();
        
        //Before adding the CellData to the cell
        CellPropertyCollection cell_protrusion_contacts_collection = p_cell->rGetCellPropertyCollection().GetPropertiesType<CellProtrusionContacts>();
        unsigned size = cell_protrusion_contacts_collection.GetSize();
        boost::shared_ptr<CellProtrusionContacts> p_cell_protrusion_contacts = boost::static_pointer_cast<CellProtrusionContacts>(cell_protrusion_contacts_collection.GetProperty());
        TS_ASSERT_THROWS_NOTHING(p_cell_protrusion_contacts);
        std::cout << "Size after creating pointer: " << size;


        std::set<unsigned> indices_example;
        indices_example.insert(1);
        indices_example.insert(3);
        p_cell_protrusion_contacts->SetCellProtrusionContacts(indices_example);
        size = cell_protrusion_contacts_collection.GetSize();
        std::cout << "Size after adding {1, 3} to one cell: " << size;

        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);

        CellPtr p_cell2 = p_cell->Divide();

        std::set<unsigned> replacement_indices_example;
        replacement_indices_example.insert(2);
        size = cell_protrusion_contacts_collection.GetSize();
        std::cout << "Size after cell division: " << size;
        boost::shared_ptr<CellProtrusionContacts> p_parentcell_data = boost::static_pointer_cast<CellProtrusionContacts>(cell_protrusion_contacts_collection.GetProperty());
        p_parentcell_data->SetCellProtrusionContacts(replacement_indices_example);
        size = cell_protrusion_contacts_collection.GetSize();
        std::cout << "Size after replacing parent set with {2}: " << size;

        TS_ASSERT_EQUALS(p_cell_protrusion_contacts->GetCellProtrusionContacts(), replacement_indices_example);

        CellPropertyCollection cell2_protrusioncontacts_collection = p_cell2->rGetCellPropertyCollection().GetPropertiesType<CellProtrusionContacts>();
        boost::shared_ptr<CellProtrusionContacts> p_daughtercell_protrusion_contacts = boost::static_pointer_cast<CellProtrusionContacts>(cell2_protrusioncontacts_collection.GetProperty());

        TS_ASSERT_EQUALS(p_daughtercell_protrusion_contacts->GetCellProtrusionContacts(), indices_example);
    }

};

#endif /*TESTCELLPROTRUSIONCONTACTS_HPP_*/
