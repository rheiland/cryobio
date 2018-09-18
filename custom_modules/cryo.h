/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version 1.3.0) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version 1.3.0) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#    llelized diffusive transport solver for 3-D biological simulations,      #
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730   #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM;
using namespace PhysiCell;

#include <vector>

//extern Cell_Definition oocyte_cell;
//extern Cell_Definition gran_cell;

void create_oocyte_cell_type( void );  //done
void create_gran_cell_type( void );  //done
void create_cell_types( void ); // done
void setup_tissue();
void add_connection(Cell* pCell_1, Cell* pCell_2);
void add_TZP_connection(Cell* pCell_1, Cell* pCell_2);
void add_basement(Cell* pCell_1);
//rwh void find_connections(Cell* pCell);
void find_connections(int idx, Cell* pCell);
void find_basement(Cell* pCell);
void find_TZP(Cell* pCell);

void custom_update_velocity_addhesion(Cell* pCell, Phenotype& phenotype, double dt);
void custom_update_velocity_from_TZP_addhesion(Cell* pCell, Phenotype& phenotype, double dt);
void custom_update_velocity_repulsion(Cell* pCell, Phenotype& phenotype, double dt);
void custom_update_velocity_from_TZP_repulsion(Cell* pCell, Phenotype& phenotype, double dt);
void custom_update_velocity_basement(Cell* pCell, Phenotype& phenotype, double dt);

// set up the microenvironment to include the immunostimulatory factor
void setup_microenvironment( void );   // done

std::vector<std::string> oocyte_coloring_function( Cell* ); // done


void oocyte_cell_phenotype_rule( Cell* pCell , Phenotype& phenotype , double dt ); // done
void oocyte_cell_rule( Cell* pCell, Phenotype& phenotype, double dt ); // done

void gran_cell_phenotype_rule( Cell* pCell , Phenotype& phenotype , double dt ); // done
void gran_cell_rule( Cell* pCell, Phenotype& phenotype, double dt ); // done
// please just work
void attach_cells( Cell* pCell_1, Cell* pCell_2 );
void dettach_cells( Cell* pCell_1 , Cell* pCell_2 );
void add_elastic_velocity( Cell* pActingOn, Cell* pAttachedTo , double elastic_constant );
void extra_elastic_attachment_mechanics( Cell* pCell, Phenotype& phenotype, double dt );
double sample_cell_boundary(Cell* pCell, Phenotype& phenotype, int solute);

void oocyte_volume_function(Cell* pCell, Phenotype& phenotype, double dt);
void gran_volume_function(Cell* pCell, Phenotype& phenotype, double dt);
