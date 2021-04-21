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

#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"

#include "DifferentiatedCellProliferativeType.hpp"

// check that adding protrusional contacts as a cell property similar to the
// nearest neighbours cell property functions as expected.
class TestCellProtrusionContacts : public AbstractCellBasedTestSuite
{
public:

    void TestCellProtrusionContactsProperty()
    {
        MAKE_PTR(CellProtrusionContacts, p_property);

        TS_ASSERT_EQUALS(p_property->GetCellCount(), 0u);
        p_property->IncrementCellCount();
        TS_ASSERT_EQUALS(p_property->GetCellCount(), 1u);
        p_property->DecrementCellCount();
        TS_ASSERT_EQUALS(p_property->GetCellCount(), 0u);
        TS_ASSERT_THROWS_THIS(p_property->DecrementCellCount(),
                "Cannot decrement cell count: no cells have this cell property");

        TS_ASSERT_EQUALS(p_property->IsType<WildTypeCellMutationState>(), false);
        TS_ASSERT_EQUALS(p_property->IsType<CellProtrusionContacts>(), true);

        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "property.arch";

        {
            AbstractCellProperty* const p_const_property = new CellProtrusionContacts();
            p_const_property->IncrementCellCount();

            TS_ASSERT_EQUALS(p_const_property->GetCellCount(), 1u);
            TS_ASSERT_THROWS_NOTHING(dynamic_cast<CellProtrusionContacts*>(p_const_property)->GetCellProtrusionContacts());

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_const_property;

            delete p_const_property;
        }

        {
            AbstractCellProperty* p_arch_property;

            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_arch_property;

            TS_ASSERT_EQUALS(p_arch_property->GetCellCount(), 1u);

            CellProtrusionContacts* p_real_property = dynamic_cast<CellProtrusionContacts*>(p_arch_property);
            TS_ASSERT(p_real_property != NULL);

            delete p_arch_property;
        }
    }


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

        CellPropertyCollection cell_protrusion_contacts_collection = collection.GetPropertiesType<CellProtrusionContacts>();
        unsigned size = cell_protrusion_contacts_collection.GetSize();
        std::cout << "When added directly, collection has size: " << size;


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
        boost::shared_ptr<AbstractCellProperty> p_protrusion_contacts(CellPropertyRegistry::Instance()->Get<CellProtrusionContacts>());

        FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
        CellPropertyCollection collection;
        collection.AddProperty(p_protrusion_contacts);
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_model, NULL, false, collection));
        p_cell->SetCellProliferativeType(p_type);
        p_cell->InitialiseCellCycleModel();
        
        //Before adding the CellData to the cell
        CellPropertyCollection cell_protrusion_contacts_collection = p_cell->rGetCellPropertyCollection().GetPropertiesType<CellProtrusionContacts>();
        unsigned size = cell_protrusion_contacts_collection.GetSize();

        TS_ASSERT_THROWS_NOTHING(boost::static_pointer_cast<CellProtrusionContacts>(cell_protrusion_contacts_collection.GetProperty()));
        std::set<unsigned> indices_example;
        indices_example.insert(1);
        indices_example.insert(3);
        boost::shared_ptr<CellProtrusionContacts> p_cell_protrusion_contacts = boost::static_pointer_cast<CellProtrusionContacts>(cell_protrusion_contacts_collection.GetProperty());
        p_cell_protrusion_contacts->SetCellProtrusionContacts(indices_example);
        size = cell_protrusion_contacts_collection.GetSize();
        TS_ASSERT_EQUALS(p_cell_protrusion_contacts->GetCellProtrusionContacts(), indices_example);

        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);

        CellPtr p_cell2 = p_cell->Divide();

        // Copy all cell data (note we create a new object not just copying the pointer)
        CellPropertyCollection daughter_property_collection = p_cell2->rGetCellPropertyCollection().GetPropertiesType<CellProtrusionContacts>();
        TS_ASSERT(daughter_property_collection.HasPropertyType<CellProtrusionContacts>());
 
        // Get the existing copy of the cell data and remove it from the daughter cell
        daughter_property_collection.RemoveProperty(p_cell_protrusion_contacts);

        // Create a new cell data object using the copy constructor and add this to the daughter cell
        MAKE_PTR_ARGS(CellProtrusionContacts, p_daughter_cell_data, (*p_cell_protrusion_contacts));
        daughter_property_collection.AddProperty(p_daughter_cell_data);

        std::set<unsigned> replacement_indices_example;
        replacement_indices_example.insert(2);
        CellPropertyCollection parent_cell_property_collection_now = p_cell->rGetCellPropertyCollection().GetPropertiesType<CellProtrusionContacts>();
        boost::shared_ptr<CellProtrusionContacts> p_parent_cell_protrusion_contacts_now = boost::static_pointer_cast<CellProtrusionContacts>(parent_cell_property_collection_now.GetProperty());
        p_parent_cell_protrusion_contacts_now->SetCellProtrusionContacts(replacement_indices_example);

        TS_ASSERT_EQUALS(p_cell_protrusion_contacts->GetCellProtrusionContacts(), replacement_indices_example);
        TS_ASSERT_EQUALS(p_parent_cell_protrusion_contacts_now->GetCellProtrusionContacts(), replacement_indices_example);

        boost::shared_ptr<CellProtrusionContacts> p_daughtercell_protrusion_contacts = boost::static_pointer_cast<CellProtrusionContacts>(daughter_property_collection.GetProperty());

        p_daughtercell_protrusion_contacts->SetCellProtrusionContacts(indices_example);

        TS_ASSERT_EQUALS(p_daughtercell_protrusion_contacts->GetCellProtrusionContacts(), indices_example);
        TS_ASSERT_EQUALS(p_parent_cell_protrusion_contacts_now->GetCellProtrusionContacts(), replacement_indices_example);

    }

    void TestCellDataOtherWayAround()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(25, 2);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_protrusion_contacts(CellPropertyRegistry::Instance()->Get<CellProtrusionContacts>());

        FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
        CellPropertyCollection collection;
        collection.AddProperty(p_protrusion_contacts);
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_model, NULL, false, collection));
        p_cell->SetCellProliferativeType(p_type);
        p_cell->InitialiseCellCycleModel();
        
        //Before adding the CellData to the cell
        CellPropertyCollection cell_protrusion_contacts_collection = p_cell->rGetCellPropertyCollection().GetPropertiesType<CellProtrusionContacts>();
        unsigned size = cell_protrusion_contacts_collection.GetSize();

        TS_ASSERT_THROWS_NOTHING(boost::static_pointer_cast<CellProtrusionContacts>(cell_protrusion_contacts_collection.GetProperty()));
        std::set<unsigned> indices_example;
        indices_example.insert(1);
        indices_example.insert(3);
        boost::shared_ptr<CellProtrusionContacts> p_cell_protrusion_contacts = boost::static_pointer_cast<CellProtrusionContacts>(cell_protrusion_contacts_collection.GetProperty());
        p_cell_protrusion_contacts->SetCellProtrusionContacts(indices_example);
        size = cell_protrusion_contacts_collection.GetSize();
        TS_ASSERT_EQUALS(p_cell_protrusion_contacts->GetCellProtrusionContacts(), indices_example);

        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);

        CellPtr p_cell2 = p_cell->Divide();

        // Copy all cell data (note we create a new object not just copying the pointer)
        CellPropertyCollection daughter_property_collection = p_cell2->rGetCellPropertyCollection().GetPropertiesType<CellProtrusionContacts>();
        TS_ASSERT(daughter_property_collection.HasPropertyType<CellProtrusionContacts>());
 
        // Get the existing copy of the cell data and remove it from the daughter cell
        daughter_property_collection.RemoveProperty(p_cell_protrusion_contacts);

        // Create a new cell data object using the copy constructor and add this to the daughter cell
        MAKE_PTR_ARGS(CellProtrusionContacts, p_daughter_cell_data, (*p_cell_protrusion_contacts));
        daughter_property_collection.AddProperty(p_daughter_cell_data);

        std::set<unsigned> replacement_indices_example;
        replacement_indices_example.insert(2);
        CellPropertyCollection parent_cell_property_collection_now = p_cell->rGetCellPropertyCollection().GetPropertiesType<CellProtrusionContacts>();
        boost::shared_ptr<CellProtrusionContacts> p_parent_cell_protrusion_contacts_now = boost::static_pointer_cast<CellProtrusionContacts>(parent_cell_property_collection_now.GetProperty());

        TS_ASSERT_EQUALS(p_cell_protrusion_contacts->GetCellProtrusionContacts(), indices_example);
        TS_ASSERT_EQUALS(p_parent_cell_protrusion_contacts_now->GetCellProtrusionContacts(), indices_example);

        boost::shared_ptr<CellProtrusionContacts> p_daughtercell_protrusion_contacts = boost::static_pointer_cast<CellProtrusionContacts>(daughter_property_collection.GetProperty());

        p_daughtercell_protrusion_contacts->SetCellProtrusionContacts(replacement_indices_example);

        TS_ASSERT_EQUALS(p_daughtercell_protrusion_contacts->GetCellProtrusionContacts(), replacement_indices_example);
        TS_ASSERT_EQUALS(p_parent_cell_protrusion_contacts_now->GetCellProtrusionContacts(), indices_example);

    }

    void TestCellVecData()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(25, 2);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
        p_cell->SetCellProliferativeType(p_type);
        p_cell->InitialiseCellCycleModel();

        // Add property CellProtrusionContacts
        MAKE_PTR(CellProtrusionContacts, p_protrusion_contacts);
        p_cell->AddCellProperty(p_protrusion_contacts);
        TS_ASSERT(p_cell->rGetCellPropertyCollection().HasPropertyType<CellProtrusionContacts>());

        CellPropertyCollection parent_cell_property_collection = p_cell->rGetCellPropertyCollection().GetPropertiesType<CellProtrusionContacts>();
        boost::shared_ptr<CellProtrusionContacts> p_parent_cell_protrusion_contacts = boost::static_pointer_cast<CellProtrusionContacts>(parent_cell_property_collection.GetProperty());

        std::set<unsigned> indices_example;
        indices_example.insert(1);
        indices_example.insert(3);
        p_parent_cell_protrusion_contacts->SetCellProtrusionContacts(indices_example);

        p_simulation_time->IncrementTimeOneStep();
        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);

        CellPtr p_daughter_cell = p_cell->Divide();

        std::set<unsigned> replacement_indices_example;
        replacement_indices_example.insert(2);
        CellPropertyCollection parent_cell_property_collection_now = p_cell->rGetCellPropertyCollection().GetPropertiesType<CellProtrusionContacts>();
        boost::shared_ptr<CellProtrusionContacts> p_parent_cell_protrusion_contacts_now = boost::static_pointer_cast<CellProtrusionContacts>(parent_cell_property_collection_now.GetProperty());

        p_parent_cell_protrusion_contacts_now->SetCellProtrusionContacts(replacement_indices_example);

        TS_ASSERT_EQUALS(p_parent_cell_protrusion_contacts_now->GetCellProtrusionContacts(), replacement_indices_example);

        CellPropertyCollection daughter_cell_property_collection = p_daughter_cell->rGetCellPropertyCollection().GetPropertiesType<CellProtrusionContacts>();
        boost::shared_ptr<CellProtrusionContacts> p_daughter_cell_protrusion_contacts = boost::static_pointer_cast<CellProtrusionContacts>(daughter_cell_property_collection.GetProperty());

        TS_ASSERT_EQUALS(p_daughter_cell_protrusion_contacts->GetCellProtrusionContacts(), indices_example);
    
    }

    void TestDifferentCellProtrusionContacts()
    {
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(25, 2);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        for (unsigned elem_index=0; elem_index<2; elem_index++)
        {
            FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();

            /* We additionally add the Cell Property protrusion contacts which contains the indices of all
            cells that this cell is in protrusion-mediated contact with. */
            CellPropertyCollection collection;
            MAKE_PTR(CellProtrusionContacts, p_protrusion_contacts); // okay......... really?
            collection.AddProperty(p_protrusion_contacts);
            
            CellPtr p_cell(new Cell(p_state, p_cell_model, NULL, false, collection));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        std::set<unsigned> indices_example;
        indices_example.insert(1);
        indices_example.insert(3);

        std::set<unsigned> replacement_indices_example;
        replacement_indices_example.insert(2);

        CellPtr p_cell1 = cells[0];
        CellPropertyCollection cell1_protrusion_contacts_collection = p_cell1->rGetCellPropertyCollection().GetPropertiesType<CellProtrusionContacts>();
        boost::shared_ptr<CellProtrusionContacts> p_cell1_protrusion_contacts = boost::static_pointer_cast<CellProtrusionContacts>(cell1_protrusion_contacts_collection.GetProperty());
        p_cell1_protrusion_contacts->SetCellProtrusionContacts(indices_example);
        TS_ASSERT_EQUALS(p_cell1_protrusion_contacts->GetCellProtrusionContacts(), indices_example);

        CellPtr p_cell2 = cells[1];
        CellPropertyCollection cell2_protrusion_contacts_collection = p_cell2->rGetCellPropertyCollection().GetPropertiesType<CellProtrusionContacts>();
        boost::shared_ptr<CellProtrusionContacts> p_cell2_protrusion_contacts = boost::static_pointer_cast<CellProtrusionContacts>(cell2_protrusion_contacts_collection.GetProperty());
        p_cell2_protrusion_contacts->SetCellProtrusionContacts(replacement_indices_example);
        TS_ASSERT_EQUALS(p_cell2_protrusion_contacts->GetCellProtrusionContacts(), replacement_indices_example);
        TS_ASSERT_EQUALS(p_cell1_protrusion_contacts->GetCellProtrusionContacts(), indices_example);

    }

};

#endif /*TESTCELLPROTRUSIONCONTACTS_HPP_*/
