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

#include "DeltaNotchReporterProtrusionDifferentiationTrackingModifier.hpp"
#include "DeltaNotchReporterProtrusionSrnModelLimi.hpp"
#include "TransitCellProliferativeType.hpp"

template<unsigned DIM>
DeltaNotchReporterProtrusionDifferentiationTrackingModifier<DIM>::DeltaNotchReporterProtrusionDifferentiationTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
DeltaNotchReporterProtrusionDifferentiationTrackingModifier<DIM>::~DeltaNotchReporterProtrusionDifferentiationTrackingModifier()
{
}

template<unsigned DIM>
void DeltaNotchReporterProtrusionDifferentiationTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void DeltaNotchReporterProtrusionDifferentiationTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
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
void DeltaNotchReporterProtrusionDifferentiationTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        if(((*cell_iter)->GetCellProliferativeType())->IsSame(p_transit_type))
            {
            // Get the set of neighbouring location indices
            std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringLocationIndices(*cell_iter);
            
            // Get the set of protrusion contacts location indices
            std::set<unsigned> protrusion_contact_indices = GetCellProtrusionContactsReference(*cell_iter)->GetCellProtrusionContacts();
            
            // Compute this cell's average neighbouring Notch and Delta concentration and store in CellData
            if (!neighbour_indices.empty())
            {
                double mean_notch = 0.0;
                double mean_delta = 0.0;
                unsigned count_transit_type_neighbours = 0;
                for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                    iter != neighbour_indices.end();
                    ++iter)
                {
                    CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                    if(((*cell_iter)->GetCellProliferativeType())->IsSame(p_transit_type))
                    {
                        double this_notch = p_cell->GetCellData()->GetItem("notch");
                        double this_delta = p_cell->GetCellData()->GetItem("delta");
                        ++count_transit_type_neighbours;
                        mean_notch += this_notch;
                        mean_delta += this_delta;
                    }
                } 
                if(count_transit_type_neighbours!=0)
                {
                    mean_notch = mean_notch / (double)count_transit_type_neighbours;
                    mean_delta = mean_delta / (double)count_transit_type_neighbours;
                    cell_iter->GetCellData()->SetItem("mean notch", mean_notch);
                    cell_iter->GetCellData()->SetItem("mean delta", mean_delta);
                }
                else
                {
                    // If this cell has no transit type neighbours as all neighbours have divided, store 0.0 for the cell data
                    cell_iter->GetCellData()->SetItem("mean notch", 0.0);
                    cell_iter->GetCellData()->SetItem("mean delta", 0.0);
                }
                
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
                unsigned count_transit_type_protrusion_neighbours = 0;
                for (std::set<unsigned>::iterator iter = protrusion_contact_indices.begin();
                    iter != protrusion_contact_indices.end();
                    ++iter)
                {
                    CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                    if(((*cell_iter)->GetCellProliferativeType())->IsSame(p_transit_type))
                    {
                        double this_notch = p_cell->GetCellData()->GetItem("notch");
                        double this_delta = p_cell->GetCellData()->GetItem("delta");
                        ++count_transit_type_protrusion_neighbours;
                        protrusion_notch += this_notch;
                        protrusion_delta += this_delta;
                    }
                } 
                if(count_transit_type_protrusion_neighbours!=0)
                {
                    protrusion_notch = protrusion_notch / (double)count_transit_type_protrusion_neighbours;
                    protrusion_delta = protrusion_delta / (double)count_transit_type_protrusion_neighbours;
                    cell_iter->GetCellData()->SetItem("protrusion notch", protrusion_notch);
                    cell_iter->GetCellData()->SetItem("protrusion delta", protrusion_delta);
                }
                else
                {
                    // If this cell has no transit type protrusion neighbours as all neighbours have divided, store 0.0 for the cell data
                    cell_iter->GetCellData()->SetItem("protrusion notch", 0.0);
                    cell_iter->GetCellData()->SetItem("protrusion delta", 0.0);
                }
            }
            else
            {
                // If this cell has no neighbours, such as an isolated cell in a CaBasedCellPopulation, store 0.0 for the cell data
                cell_iter->GetCellData()->SetItem("protrusion notch", 0.0);
                cell_iter->GetCellData()->SetItem("protrusion delta", 0.0);
            }
        }
    }
}

template<unsigned DIM>
void DeltaNotchReporterProtrusionDifferentiationTrackingModifier<DIM>::SetProtrusionNeighbours(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // this will associate for each cell a set of the indices of all the cells that this cell is in protrusion-mediated
    // contact with, which depends on the length of the protrusion, its angle, its angle of opening,
    // and the activation threshold which we choose

    // get reference to mesh
    AbstractMesh<DIM, DIM>& r_mesh = rCellPopulation.rGetMesh();
    c_vector<double, DIM> first_cell_location = rCellPopulation.GetLocationOfCellCentre(rCellPopulation.GetCellUsingLocationIndex(0));
    c_vector<double, DIM> last_cell_location = rCellPopulation.GetLocationOfCellCentre(rCellPopulation.GetCellUsingLocationIndex(rCellPopulation.GetNumAllCells()-1));
    c_vector<double, DIM> middle_location = 0.5*(last_cell_location - first_cell_location);
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter) 
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
        // get location index of current cell
        unsigned index_A = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        // get location vector of current cell
        c_vector<double, DIM> location_vector_A = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        location_vector_A = location_vector_A - middle_location;

        for (typename AbstractCellPopulation<DIM>::Iterator possible_contact_cell_iter = rCellPopulation.Begin();
         possible_contact_cell_iter != rCellPopulation.End();
         ++possible_contact_cell_iter) // iterate over all cells to find all possible contacts
        {
            // get location index of possible contact
            unsigned index_B = rCellPopulation.GetLocationIndexUsingCell(*possible_contact_cell_iter);

            // possible implementation, although the confusion of elements/nodes in the vertex mesh 
            // means this is confusing
            // c_vector<double, DIM> vector_between_cells = r_mesh.GetVectorFromAtoB(r_mesh.GetNode(index_A)->rGetLocation(), r_mesh.GetNode(index_B)->rGetLocation());
            // c_vector<double, DIM> unit_vector_between_cells = vector_between_cells/norm_2(vector_between_cells);
            
            c_vector<double, DIM> location_vector_B = rCellPopulation.GetLocationOfCellCentre(*possible_contact_cell_iter);
            location_vector_B = location_vector_B - middle_location;
            c_vector<double, DIM> vector_between_cells = r_mesh.GetVectorFromAtoB(location_vector_A, location_vector_B);
            c_vector<double, DIM> unit_vector_between_cells = vector_between_cells/norm_2(vector_between_cells);

            // protrusion contact is dependent on the distance between cells
            double distance_between_cells = norm_2(vector_between_cells);
            // // Baum model; does not incorporate variation in protrusion lengths
            // if ( (1.5 < distance_between_cells) && (distance_between_cells < (2*(this_protrusion_length+this_protrusion_tip_length))) )
            // // Hypothetical signalling/receiving neighborhood model; does not incorporate variation in protrusion lengths
            // if ( ((this_protrusion_length-this_protrusion_tip_length) < distance_between_cells) && (distance_between_cells < (2*(this_protrusion_length+this_protrusion_tip_length))) )
            // Bajpai model
            if ( (2*(this_protrusion_length-this_protrusion_tip_length) < distance_between_cells) && (distance_between_cells < (2*(this_protrusion_length+this_protrusion_tip_length))) )
            {
                if (sin(this_protrusion_angular_opening)*sin(this_protrusion_angular_opening) > angular_activation_threshold)
                {   
                    protrusion_contact_indices.insert(index_B);
                }
                else
                {
                    // calculate the dot product and divide by 2
                    
                    // // polarized projections; projections point away from the centre
                    // c_vector<double, DIM> unit_location_vector_A = location_vector_A/norm_2(location_vector_A);
                    // c_vector<double, DIM> unit_location_vector_B = location_vector_B/norm_2(location_vector_B);
                    // double dot_product_1;
                    // dot_product_1 = unit_vector_between_cells[0]*unit_location_vector_A[0] + unit_vector_between_cells[1]*unit_location_vector_A[1];
                    // double dot_product_2;
                    // dot_product_2 = unit_vector_between_cells[0]*unit_location_vector_B[0] + unit_vector_between_cells[1]*unit_location_vector_B[1];

                    // polarized projections; projections point in a circle
                    c_vector<double, DIM> unit_location_vector_A = location_vector_A/norm_2(location_vector_A);
                    c_vector<double, DIM> unit_location_vector_B = location_vector_B/norm_2(location_vector_B);
                    double dot_product_1;
                    dot_product_1 = unit_vector_between_cells[0]*unit_location_vector_A[1] - unit_vector_between_cells[1]*unit_location_vector_A[0];
                    double dot_product_2;
                    dot_product_2 = unit_vector_between_cells[0]*unit_location_vector_B[1] - unit_vector_between_cells[1]*unit_location_vector_B[0];

                    // // if all cells have same angle of protrusion
                    // double dot_product_1;
                    // dot_product_1 = unit_vector_between_cells[0]*cos(this_protrusion_angle) + unit_vector_between_cells[1]*sin(this_protrusion_angle);
                    // double dot_product_2;
                    // dot_product_2 = unit_vector_between_cells[0]*cos(this_protrusion_angle) + unit_vector_between_cells[1]*sin(this_protrusion_angle);
                    if ( (0.5*(dot_product_1*dot_product_1 + dot_product_2*dot_product_2)) > angular_activation_threshold)
                    {
                        protrusion_contact_indices.insert(index_B);
                    }
                }
            }
        }
        GetCellProtrusionContactsReference(*cell_iter)->SetCellProtrusionContacts(protrusion_contact_indices);
    }
}

template<unsigned DIM>
boost::shared_ptr<CellProtrusionContacts> DeltaNotchReporterProtrusionDifferentiationTrackingModifier<DIM>::GetCellProtrusionContactsReference(CellPtr p_cell)
{
    CellPropertyCollection cell_protrusion_contacts_collection = p_cell->rGetCellPropertyCollection().GetPropertiesType<CellProtrusionContacts>();
    return boost::static_pointer_cast<CellProtrusionContacts>(cell_protrusion_contacts_collection.GetProperty());
}

template<unsigned DIM>
void DeltaNotchReporterProtrusionDifferentiationTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class DeltaNotchReporterProtrusionDifferentiationTrackingModifier<1>;
template class DeltaNotchReporterProtrusionDifferentiationTrackingModifier<2>;
template class DeltaNotchReporterProtrusionDifferentiationTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DeltaNotchReporterProtrusionDifferentiationTrackingModifier)
