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

#include "DeltaNotchReporterProtrusionTrackingModifier.hpp"
#include "DeltaNotchReporterProtrusionSrnModelLimi.hpp"

template<unsigned DIM>
DeltaNotchReporterProtrusionTrackingModifier<DIM>::DeltaNotchReporterProtrusionTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
DeltaNotchReporterProtrusionTrackingModifier<DIM>::~DeltaNotchReporterProtrusionTrackingModifier()
{
}

template<unsigned DIM>
void DeltaNotchReporterProtrusionTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void DeltaNotchReporterProtrusionTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */

    // We first add the property CellProtrusionContacts cells to each cell using the designated helper class
    SetProtrusionNeighbours(rCellPopulation);

    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void DeltaNotchReporterProtrusionTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    // First recover each cell's Notch, Delta, and Reporter concentrations from the ODEs and store in CellData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        DeltaNotchReporterProtrusionSrnModelLimi* p_model = static_cast<DeltaNotchReporterProtrusionSrnModelLimi*>(cell_iter->GetSrnModel());
        double this_delta = p_model->GetDelta();
        double this_notch = p_model->GetNotch();
        double this_reporter = p_model->GetReporter();

        // Note that the state variables must be in the same order as listed in ode system file
        cell_iter->GetCellData()->SetItem("notch", this_notch);
        cell_iter->GetCellData()->SetItem("delta", this_delta);
        cell_iter->GetCellData()->SetItem("reporter", this_reporter);
    }

    // Next iterate over the population to compute and store the concentration of each substance received
    // by neighbouring contact ("mean") or by protrusional contact ("protrusion") in CellData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get the set of neighbouring location indices
        std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringLocationIndices(*cell_iter);
        
        // if (*cell_iter->HasCellProperty<CellProtrusionContacts>()) 
        // p_cell->rGetCellPropertyCollection().HasProperty(p_wild_type) // isn't working? confused
        // p_cell->rGetCellPropertyCollection().GetPropertiesType<CellProtrusionContacts>()
        if (*cell_iter->rGetCellPropertyCollection().HasPropertyType<CellProtrusionContacts>())
        {
        std::set<unsigned> protrusion_contact_indices = GetCellProtrusionContactsReference(*cell_iter)->GetCellProtrusionContacts();
        }
        else
        {
           EXCEPTION("Protrusion Contacts cell property not implemented");
        }
        
        // Compute this cell's average neighbouring Notch and Delta concentration and store in CellData
        if (!neighbour_indices.empty())
        {
            double mean_notch = 0.0;
            double mean_delta = 0.0;
            for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                 iter != neighbour_indices.end();
                 ++iter)
            {
                CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                double this_notch = p_cell->GetCellData()->GetItem("notch");
                double this_delta = p_cell->GetCellData()->GetItem("delta");
                mean_notch += this_notch/neighbour_indices.size();
                mean_delta += this_delta/neighbour_indices.size();
            } 
            cell_iter->GetCellData()->SetItem("mean notch", mean_notch);
            cell_iter->GetCellData()->SetItem("mean delta", mean_delta);
        }
        else
        {
            // If this cell has no neighbours, such as an isolated cell in a CaBasedCellPopulation, store 0.0 for the cell data
            cell_iter->GetCellData()->SetItem("mean notch", 0.0);
            cell_iter->GetCellData()->SetItem("mean delta", 0.0);
        }

        // Compute this cell's average protrusion-mediated received Notch and Delta concentration and store in CellData
        if (!protrusion_contact_indices.empty())
        {
            double protrusion_notch = 0.0;
            double protrusion_delta = 0.0;
            for (std::set<unsigned>::iterator iter = protrusion_contact_indices.begin();
                 iter != protrusion_contact_indices.end();
                 ++iter)
            {
                CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                double this_notch = p_cell->GetCellData()->GetItem("notch");
                double this_delta = p_cell->GetCellData()->GetItem("delta");
                protrusion_notch += this_notch/protrusion_contact_indices.size(); // currently unweighted by protrusion overlap
                protrusion_delta += this_delta/protrusion_contact_indices.size();
            } 
            cell_iter->GetCellData()->SetItem("protrusion notch", protrusion_notch);
            cell_iter->GetCellData()->SetItem("protrusion delta", protrusion_delta);
        }
        else
        {
            // If this cell has no neighbours, such as an isolated cell in a CaBasedCellPopulation, store 0.0 for the cell data
            cell_iter->GetCellData()->SetItem("protrusion notch", 0.0);
            cell_iter->GetCellData()->SetItem("protrusion delta", 0.0);
        }
    }
}

template<unsigned DIM>
void DeltaNotchReporterProtrusionTrackingModifier<DIM>::SetProtrusionNeighbours(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // this will associate for each cell a set of the indices of all the cells that this cell is in protrusion-mediated
    // contact with, which depends on the length of the protrusion, its angle, its angle of opening,
    // and the activation threshold which we choose

    // get reference to mesh
    AbstractMesh<DIM, DIM>& r_mesh = rCellPopulation.rGetMesh();
    
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter) // dereferencing iterator gives current cell
    {
        // make empty set for indices
        std::set<unsigned> protrusion_contact_indices;
        // retrieve protrusion details 
        DeltaNotchReporterProtrusionSrnModelLimi* p_model = static_cast<DeltaNotchReporterProtrusionSrnModelLimi*>(cell_iter->GetSrnModel());
        double this_protrusion_length = p_model->GetProtrusionLength();
        double this_protrusion_tip_length = p_model->GetProtrusionTipLength();
        double this_protrusion_angle = p_model->GetProtrusionAngle();
        double this_protrusion_angular_opening = p_model->GetProtrusionAngularOpening();
        double angular_activation_threshold = p_model->GetAngularActivationThreshold();
        // get location index of cell we are looking at
        unsigned index_A = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        
        // Find the indices of the elements owned by this node
        for (typename AbstractCellPopulation<DIM>::Iterator possible_contact_cell_iter = rCellPopulation.Begin();
         possible_contact_cell_iter != rCellPopulation.End();
         ++possible_contact_cell_iter) // iterate over all cells to find all possible contacts
        {
            // get location index of possible contact
            unsigned index_B = rCellPopulation.GetLocationIndexUsingCell(*possible_contact_cell_iter);
                    
            c_vector<double, DIM> vector_between_cells = r_mesh.GetVectorFromAtoB(r_mesh.GetNode(index_A)->rGetLocation(), r_mesh.GetNode(index_B)->rGetLocation());
            c_vector<double, DIM> unit_vector_between_cells = vector_between_cells/norm_2(vector_between_cells);
            double distance_between_nodes = norm_2(vector_between_cells);
            if ( ((2*(this_protrusion_length-this_protrusion_tip_length))<distance_between_nodes) && 
            (distance_between_nodes < (2*(this_protrusion_length+this_protrusion_tip_length))) )
            {
                if (sin(this_protrusion_angular_opening)*sin(this_protrusion_angular_opening) > angular_activation_threshold)
                {
                    protrusion_contact_indices.insert(index_B);
                }
                else
                {
                    // calculate the dot product and divide by 2
                    double dot_product_1;
                    dot_product_1 = unit_vector_between_cells[0]*cos(this_protrusion_angle) + unit_vector_between_cells[1]*sin(this_protrusion_angle);
                    double dot_product_2;
                    dot_product_2 = unit_vector_between_cells[0]*cos(this_protrusion_angle) + unit_vector_between_cells[1]*sin(this_protrusion_angle);

                    if ( (0.5*(dot_product_1*dot_product_1 + dot_product_2*dot_product_2)) > angular_activation_threshold)
                    {
                        protrusion_contact_indices.insert(index_B);
                    }
                }
            }
            GetCellProtrusionContactsReference(*cell_iter)->SetCellProtrusionContacts(protrusion_contact_indices);
        }
    }
}

template<unsigned DIM>
boost::shared_ptr<CellProtrusionContacts> DeltaNotchReporterProtrusionTrackingModifier<DIM>::GetCellProtrusionContactsReference(CellPtr p_cell)
{
    CellPropertyCollection cell_protrusion_contacts_collection = p_cell->rGetCellPropertyCollection().GetPropertiesType<CellProtrusionContacts>();
    return boost::static_pointer_cast<CellProtrusionContacts>(cell_protrusion_contacts_collection.GetProperty());
}


template<unsigned DIM>
void DeltaNotchReporterProtrusionTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class DeltaNotchReporterProtrusionTrackingModifier<1>;
template class DeltaNotchReporterProtrusionTrackingModifier<2>;
template class DeltaNotchReporterProtrusionTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DeltaNotchReporterProtrusionTrackingModifier)
