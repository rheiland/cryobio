#include "./cryo.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>

#include <cassert>  //rwh

//Establish constants------------------------
static double time_called=0;
#define PI 3.14159265
#define R 8314//granulosa (10^-3 J/mole*k)
constexpr auto R_atm = 0.0821;//oocyte;

//-------------------------- CPAs and Solutes ----------------------------------------------
//---------CPA's and water @ 25 -----------------------
static double V_bar_Water = .018;
static double V_bar_Glycerol = .071;
static double V_bar_EG = .054;//L/mol 
static double V_bar_DMSO = .069;
static double V_bar_xPBS = .070;
//---------Salt----------------------------------
static double V_bar_n = .01;//temporary value->check lit.

std::unordered_map< int, std::vector <Cell* > > cell_id_nbrs; 

Cell_Definition oocyte_cell;
Cell_Definition gran_cell;

// rwh - indices to custom vars
static int oocyte_initial_vol_index=0;
static int gran_initial_vol_index=0;


static int transport_idx;
static int oocyte_connections_idx;
static int oocyte_TZP_status_idx;
static int oocyte_Ps_idx;
static int oocyte_CPA_fluxes_idx;
static int oocyte_Inter_CPA_moles_idx;
static int oocyte_Ext_CPA_concentration_idx;
static int oocyte_Prev_Inter_CPA_moles_idx;

static int gran_connections_idx;
static int gran_TZP_status_idx;
static int gran_Ps_idx;
static int gran_CPA_fluxes_idx;
static int gran_Inter_CPA_moles_idx;
static int gran_Ext_CPA_concentration_idx;
static int gran_Prev_Inter_CPA_moles_idx;



//initial conditions--------------------------
double times=0;


//-------------------------- Construct Cells and Tissue ----------------------------------------------
void create_oocyte_cell_type(void)
{
	oocyte_cell = cell_defaults;
	oocyte_cell.name = "oocyte";
	oocyte_cell.type = 1;

	// turn off proliferation;

	int cycle_start_index = live.find_phase_index(PhysiCell_constants::live);
	int cycle_end_index = live.find_phase_index(PhysiCell_constants::live);
	oocyte_cell.phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = 0.0;

	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);

	// solute
	// rwh - change index to 0
//	oocyte_cell.phenotype.secretion.uptake_rates[1] = 0;
//	oocyte_cell.phenotype.secretion.secretion_rates[1] = 0;
	oocyte_cell.phenotype.secretion.uptake_rates[0] = 0;
	oocyte_cell.phenotype.secretion.secretion_rates[0] = 0;

	//rwh-comment out
	// oocyte_cell.phenotype.secretion.uptake_rates[3] = 0;
	// oocyte_cell.phenotype.secretion.secretion_rates[3] = 0;

	//CPA  (rwh: change index 2 to 1)
	// oocyte_cell.phenotype.secretion.uptake_rates[2] = 0; //
	// oocyte_cell.phenotype.secretion.secretion_rates[2] = 0; //
	// oocyte_cell.phenotype.secretion.saturation_densities[2] = 1;	
	oocyte_cell.phenotype.secretion.uptake_rates[1] = 0; //
	oocyte_cell.phenotype.secretion.secretion_rates[1] = 0; //
	oocyte_cell.phenotype.secretion.saturation_densities[1] = 1;


	// set apoptosis rate to 0
	oocyte_cell.phenotype.death.rates[apoptosis_index] = 0;

	// turn on motility;
	oocyte_cell.phenotype.motility.is_motile = true;
	oocyte_cell.phenotype.motility.restrict_to_2D = true;
	oocyte_cell.phenotype.mechanics.cell_cell_adhesion_strength *= 0;
	oocyte_cell.phenotype.mechanics.cell_cell_repulsion_strength *= 0;

	// set functions - rwh comment out?
	// NOTE: calling user-defined cell functions need to be thread-safe since 
	//       they're invoked inside the OpenMP loop (rf. PhysiCell_cell.cpp)
	// oocyte_cell.functions.update_phenotype = oocyte_cell_phenotype_rule;
	// oocyte_cell.functions.custom_cell_rule = oocyte_cell_rule;

	oocyte_cell.phenotype.motility.migration_bias = 1;


	// set custom data values
	//rwh
//	oocyte_cell.custom_data["initial volume"] = 180000.0;
//rwh	oocyte_cell.custom_data["initial volume"] = 524.0;
	oocyte_cell.custom_data.add_variable( "initial volume" , "dimensionless", 524.0 ); 
//	initial_vol_index = pCell->custom_data.find_variable_index( "activation" );


	Custom_Cell_Data ccd;
	std::vector<int> connected_cells(1000, 0);
//rwh	oocyte_connections_idx = oocyte_cell.custom_data.add_vector_variable("connections", connected_cells);
//	std::cout << "oocyte_connections_idx = " << oocyte_connections_idx << "\n"; //rwh:  "= 1"
	oocyte_cell.custom_data.add_int_vector("connections", connected_cells);
//rwh	oocyte_cell.custom_data.add_variable( "connections_idx" , oocyte_connections_idx );  // crap, only double!

	std::vector<int> connected_TZP(1000, 0);
//	oocyte_TZP_status_idx = oocyte_cell.custom_data.add_vector_variable("TZP_status", connected_TZP);
//	std::cout << "oocyte_TZP_status_idx = " << oocyte_TZP_status_idx<< "\n";
	oocyte_cell.custom_data.add_int_vector("TZP_status", connected_TZP);


	oocyte_cell.custom_data.add_variable("basement", "unitless", -1.0);
	//Parameters

	oocyte_cell.custom_data.add_variable("Ps", "unitless", .18);  //rwh: beware - "Ps" both a scalar and vector??
	oocyte_cell.custom_data.add_variable("Lp", "M/s", 1.08 / 60);
	oocyte_cell.custom_data.add_variable("R", "unitless", R_atm);
	oocyte_cell.custom_data.add_variable("Area", "unitless", 408);
	oocyte_cell.custom_data.add_variable("Temperature", "unitless", 294);
	oocyte_cell.custom_data.add_variable("TZP_length", "unitless", 6);
	
	oocyte_cell.custom_data.add_variable("Inter_salt_moles", "unitless", .2);
	oocyte_cell.custom_data.add_variable("Inter_sucrose_moles", "unitless", 0.0);
	
	oocyte_cell.custom_data.add_variable("Ext_salt_concentration", "unitless", .2);
	oocyte_cell.custom_data.add_variable("Ext_sucrose_concentration", "unitless", 0.0);
	
	oocyte_cell.custom_data.add_variable("Water_volume", "unitless", 1.0);
	oocyte_cell.custom_data.add_variable("Water_flux", "unitless", 0.0);
	oocyte_cell.custom_data.add_variable("Prev_Water_volume", "unitless", 0.0);

	std::vector<double> Ps_Vector(4, 0.0);
//	oocyte_Ps_idx = oocyte_cell.custom_data.add_vector_variable("Ps", Ps_Vector);
	oocyte_cell.custom_data.add_vector("Ps_vec", Ps_Vector);  //rwh: change name, append "_vec"

	std::vector<double> CPA_fluxes(4, 0.0);
	oocyte_cell.custom_data.add_vector("CPA_fluxes", CPA_fluxes);

	std::vector<double> Inter_CPA_moles(4, 0.0);
	oocyte_cell.custom_data.add_vector("Inter_CPA_moles", Inter_CPA_moles);

	// std::cout << "---Inter_CPA_moles[0]=" << oocyte_cell.custom_data["Inter_CPA_moles"[0]] << "\n";
	// std::cout << "Inter_CPA_moles[1]=" << oocyte_cell.custom_data["Inter_CPA_moles"[1]] << "\n";
	// std::cout << "Inter_CPA_moles[2]=" << oocyte_cell.custom_data["Inter_CPA_moles"[2]] << "\n";
	// std::cout << "Inter_CPA_moles[3]=" << oocyte_cell.custom_data["Inter_CPA_moles"[3]] << "\n";

	std::vector<double> Ext_CPA_concentration(4, 0.0);
	oocyte_cell.custom_data.add_vector("Ext_CPA_concentration", Ext_CPA_concentration);

	std::vector<double> Prev_Inter_CPA_moles(4, 0.0);
	oocyte_cell.custom_data.add_vector("Prev_Inter_CPA_moles", Prev_Inter_CPA_moles);

	return;
}
void create_gran_cell_type(void)
{
	//Define cell type
	gran_cell = cell_defaults;
	gran_cell.name = "granulosa cell";
	gran_cell.type = 2;

	// turn on life cycle, make the cells live;
	int cycle_start_index = live.find_phase_index(PhysiCell_constants::live);
	int cycle_end_index = live.find_phase_index(PhysiCell_constants::live);

	// Set life cycle constant and turn off apopotosis
	int apoptosis_index = gran_cell.phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
	gran_cell.phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = 0.0;
	gran_cell.phenotype.death.rates[apoptosis_index] = 0.0;

	//diffusion rates
	//rwh - change index - just 0 and 1 now
	gran_cell.phenotype.secretion.uptake_rates[1] = 0; //
	gran_cell.phenotype.secretion.secretion_rates[1] = 0; //
	gran_cell.phenotype.secretion.uptake_rates[0] = 0; //
	gran_cell.phenotype.secretion.secretion_rates[0] = 0; //
//	gran_cell.phenotype.secretion.uptake_rates[3] = 0; //
//	gran_cell.phenotype.secretion.secretion_rates[3] = 0; //


	//CPA
	//rwh - change index - just 0 and 1 now
//	gran_cell.phenotype.secretion.uptake_rates[2] = 0; //
//	gran_cell.phenotype.secretion.secretion_rates[2] = 0; //
//	gran_cell.phenotype.secretion.saturation_densities[2] = 1;
	gran_cell.phenotype.secretion.uptake_rates[1] = 0; //
	gran_cell.phenotype.secretion.secretion_rates[1] = 0; //
	gran_cell.phenotype.secretion.saturation_densities[1] = 1;


	// turn on motility
	gran_cell.phenotype.motility.is_motile = true;
	gran_cell.phenotype.motility.restrict_to_2D = true;

	// set Mechanics
	gran_cell.phenotype.mechanics.cell_cell_adhesion_strength *= 0;
	gran_cell.phenotype.mechanics.cell_cell_repulsion_strength *= 0;
	gran_cell.phenotype.motility.migration_bias = 1;

	//
	//gran_cell.cell_state.simple_pressure;


	// set functions
	// rwh - comment out functions for now?
	// NOTE: calling user-defined cell functions need to be thread-safe since 
	//       they're invoked inside the OpenMP loop (rf. PhysiCell_cell.cpp)
	gran_cell.functions.update_phenotype = gran_cell_phenotype_rule;
	gran_cell.functions.custom_cell_rule = gran_cell_rule;

	// set custom data values
//	std::vector<double> connected_cells(1000, 0.0);
	std::vector<int> connected_cells(1000, 0.0);
//	gran_connections_idx = gran_cell.custom_data.add_vector_variable("connections", connected_cells);
	gran_cell.custom_data.add_int_vector("connections", connected_cells);

//	std::vector<double> connected_TZP(1000, 0.0);
	std::vector<int> connected_TZP(1000, 0.0);
//	gran_TZP_status_idx = gran_cell.custom_data.add_vector_variable("TZP_status", connected_TZP);
	gran_cell.custom_data.add_int_vector("TZP_status", connected_TZP);
	// set custom data values

	Custom_Cell_Data ccd;
	//rwh
//	gran_cell.custom_data["initial volume"] = 180000.0;
//	gran_cell.custom_data["initial volume"] = 524.0;
	gran_cell.custom_data.add_variable( "initial volume" , "dimensionless", 524.0 ); 

	gran_cell.custom_data.add_variable("basement", "unitless", -1.0);
	

	//Parameters
	gran_cell.custom_data.add_variable("Ps", "unitless", 21 / 25.241);
	gran_cell.custom_data.add_variable("Lp", "unitless", 5.333 / 25.241);
	gran_cell.custom_data.add_variable("Area", "unitless", 0.0);
	gran_cell.custom_data.add_variable("R", "unitless", 1);
	gran_cell.custom_data.add_variable("Temperature", "unitless", 1);
	gran_cell.custom_data.add_variable("connection_length", "unitless", -1.0);

	gran_cell.custom_data.add_variable("Inter_salt_moles", "unitless", .2);
	gran_cell.custom_data.add_variable("Inter_sucrose_moles", "unitless", 0.0);

	gran_cell.custom_data.add_variable("Ext_salt_concentration", "unitless", .2);
	gran_cell.custom_data.add_variable("Ext_sucrose_concentration", "unitless", 0.0);

	gran_cell.custom_data.add_variable("Water_volume", "unitless", 1.0);
	gran_cell.custom_data.add_variable("Water_flux", "unitless", 0.0);
	gran_cell.custom_data.add_variable("Prev_Water_volume", "unitless", 0.0);

	std::vector<double> Ps_Vector(4, 0.0);
//	gran_Ps_idx = gran_cell.custom_data.add_vector_variable("Ps", Ps_Vector);
	gran_cell.custom_data.add_vector("Ps_vec", Ps_Vector);  //rwh: change name "Ps" --> "Ps_vec"

	std::vector<double> CPA_fluxes(4, 0.0);
//	gran_CPA_fluxes_idx = gran_cell.custom_data.add_vector_variable("CPA_fluxes", CPA_fluxes);
	gran_cell.custom_data.add_vector("CPA_fluxes", CPA_fluxes);

	std::vector<double> Inter_CPA_moles(4, 0.0);
	gran_cell.custom_data.add_vector("Inter_CPA_moles", Inter_CPA_moles);

	std::vector<double> Ext_CPA_concentration(4, 0.0);
	gran_cell.custom_data.add_vector("Ext_CPA_concentration", Ext_CPA_concentration);

	std::vector<double> Prev_Inter_CPA_moles(4, 0.0);
	gran_cell.custom_data.add_vector("Prev_Inter_CPA_moles", Prev_Inter_CPA_moles);
	
	return;
}
void create_cell_types(void)
{
	// use the same random seed so that future experiments have the
	// same initial histogram of oncoprotein, even if threading means
	// that future division and other events are still not identical
	// for all runs
	SeedRandom(0);

	// housekeeping

	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment(&microenvironment);

	// turn the default cycle model to live,
	// so it's easier to turn off proliferation

	cell_defaults.phenotype.cycle.sync_to_cycle_model(live);

	// Make sure we're ready for 2D

	cell_defaults.functions.set_orientation = up_orientation;

	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; // true;

	// set to no motility for cancer cells
	cell_defaults.phenotype.motility.is_motile = true;

	// use default proliferation and death

	int cycle_start_index = live.find_phase_index(PhysiCell_constants::live);
	int cycle_end_index = live.find_phase_index(PhysiCell_constants::live);

	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);


	// set default uptake and secretion
	// oxygen
	//rwh - but elsewhere:  use_oxygen_as_first_field = false;  //rwh??
	cell_defaults.phenotype.secretion.secretion_rates[0] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[0] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[0] = 0;

	cell_defaults.name = "default cell";
	cell_defaults.type = 0;

	// add custom data

	// for oocyte
	cell_defaults.custom_data.add_variable("elastic coefficient", "1/sec", 0.05);
	cell_defaults.custom_data.add_variable("receptor", "dimensionless", 0.0);
	cell_defaults.custom_data.add_variable("oocyte release oxygen threshold", "mmHg", 15.0);

	cell_defaults.custom_data.add_variable("damage rate", "1/sec", 1.0 / 30.0);
	cell_defaults.custom_data.add_variable("repair rate", "1/sec", 1.0 / 240.0);
	cell_defaults.custom_data.add_variable("drug death rate", "1/sec", 1.0 / 240.0);
	cell_defaults.custom_data.add_variable("damage", "dimensionless", 0.0);

	std::vector<double> wsvalues = { 0,0,0,0,0,0,0,0,0 };
//	transport_idx = cell_defaults.custom_data.add_vector_variable("transport", wsvalues);
	cell_defaults.custom_data.add_vector("transport", wsvalues);

	// create the oocyte types
	create_oocyte_cell_type();

	// create the gran types
	create_gran_cell_type();

	return;
}
void setup_tissue(void)
{
	static double foll_radius = 60.;
	//place oocyte
	std::cout << "\tPlacing oocyte cells ... " << std::endl;

	// place the oocyte cell at origin
	Cell* pC;
	pC = create_cell(oocyte_cell);

	//rwh
	// pC->custom_data["foobar"] =42.0;
	// std::cout << "-------- foobar=" <<pC->custom_data["foobar"] << "\n";
	// std::cout << "-------- Foo bar=" <<pC->custom_data["Foo bar"] << "\n";
	// std::cout << "-------- wth=" <<pC->custom_data["wth"] << "\n";

	// std::cout << "setup oocyte: Inter_CPA_moles[0]=" << pC->custom_data["Inter_CPA_moles"[0]] << "\n";
	// std::cout << "setup: Inter_CPA_moles[1]=" << pC->custom_data["Inter_CPA_moles"[1]] << "\n";

	pC->assign_position(0, 0, 0); //oocyte position

	static int oocyte_initial_vol_index = pC->custom_data.find_variable_index( "initial volume" );

//rwh	pC->set_total_volume(pC->custom_data["initial volume"]); //sets initial total volume
	pC->set_total_volume(pC->custom_data[oocyte_initial_vol_index]); //sets initial total volume

	//rwh 
//	static double gran_radius = pow((3 / (4 * PI)*(pC->custom_data["initial volume"])), (1 / 3));
//	double gran_radius = pow((3 / (4 * PI)*(pC->custom_data["initial volume"])), (1 / 3));
//	std::cout << "orig gran_radius ="<<gran_radius<<std::endl;
//	gran_radius /= 10.;
//	gran_radius = pow((3. / (4. * PI)*(pC->custom_data["initial volume"])), (1. / 3.));
	double gran_radius = pow((3. / (4. * PI)*(pC->custom_data[oocyte_initial_vol_index])), (1. / 3.));
//	gran_radius /= 10.;
//	gran_radius = 5.;
	std::cout << "new gran_radius ="<<gran_radius<<std::endl;
	
	pC->is_movable = true; //set oocyte to moveable

	// place a cluster of granulosa around the oocyte
	std::cout << "\tPlacing granulosa cells ... " << std::endl;

	double cell_spacing = 0.95 * 2.0 * gran_radius;
	std::cout << "cell_spacing = " << cell_spacing << std::endl;
	Cell* pCell = NULL;
	double x = 0.0;
	double x_outer = foll_radius;
	//rwh 
//	static double oocyte_radius = pow((3 / (4 * PI)*(pCell->custom_data["initial volume"])), (1 / 3));
//	std::cout << PI <<","<<1/3<<std::endl;  // = 0 due to integer division
//	std::cout << pCell->custom_data["initial volume"] <<std::endl;
//	double init_vol = pC->custom_data["initial volume"] ;
	double init_vol = pC->custom_data[oocyte_initial_vol_index] ;
	std::cout << "pC init_vol = "<<init_vol <<std::endl;
//	double oocyte_radius = pow((3 / (4 * PI)*(pCell->custom_data["initial volume"])), (1 / 3));
//	double oocyte_radius = pow((3.0 / (4.0 * PI)*(pC->custom_data["initial volume"])), (1.0 / 3));
	double oocyte_radius = pow((3.0 / (4.0 * PI)*(pC->custom_data[oocyte_initial_vol_index])), (1.0 / 3));
	std::cout << "orig oocyte_radius ="<<oocyte_radius<<std::endl;
//	oocyte_radius /= 10.;
//	oocyte_radius = 5.;
	std::cout << "new oocyte_radius ="<<oocyte_radius<<std::endl;

	double x_inner = oocyte_radius;
	double y_inner = oocyte_radius;
	double y = 0.0;
	int n = 0;

	while (y < foll_radius)
	{
		x = 0.0;
		x_inner = sqrt(oocyte_radius*oocyte_radius - y * y);
		y_inner = sqrt(oocyte_radius*oocyte_radius - x * x);
		x_outer = sqrt(foll_radius*foll_radius - y * y);

		if (n % 2 == 1)
		{
			x = 0.5*cell_spacing;
		}

		while (x < x_outer)
		{
			// if not inside oocyte place granulosa
			if (fabs(x) > (x_inner) || fabs(y) > (y_inner))
			{
				pCell = create_cell(gran_cell); // granulosas in quadrant 1
				pCell->assign_position(x, y, 0.0);
				//rwh
//	std::cout << "setup gran: Inter_CPA_moles[0]=" << pC->custom_data["Inter_CPA_moles"[0]] << "\n";
//	std::cout << "setup gran: Inter_CPA_moles[1]=" << pC->custom_data["Inter_CPA_moles"[1]] << "\n";
//				pCell->set_total_volume(pCell->custom_data["initial volume"]);
				pCell->set_total_volume(pCell->custom_data[oocyte_initial_vol_index]);  //rwh
				
				pCell->is_movable = true;
				if (fabs(y) > (0.01))
				{

					pCell = create_cell(gran_cell); // granulosa in quadrant 4
					pCell->assign_position(x, -y, 0.0);
					pCell->set_total_volume(pCell->custom_data[oocyte_initial_vol_index]);  //rwh
					
					pCell->is_movable = true;
				}

				if (fabs(x) > (0.01))
				{
					pCell = create_cell(gran_cell); // granulosa cell in quadrant 2
					pCell->assign_position(-(x), y, 0.0);
					//pCell->set_total_volume(pCell->custom_data["initial volume"]);  //rwh
					pCell->set_total_volume(pCell->custom_data[oocyte_initial_vol_index]);  //rwh
					
					pCell->is_movable = true;

					if (fabs(y) > (0.01))
					{
						pCell = create_cell(gran_cell); // granulosa cell in quadrant 3
						pCell->assign_position(-(x), -(y), 0.0);
//rwh						pCell->set_total_volume(pCell->custom_data["initial volume"]);
						pCell->set_total_volume(pCell->custom_data[oocyte_initial_vol_index]);
						
						pCell->is_movable = true;
					}
				}
			}
			x += cell_spacing;
		}

		y += cell_spacing * sqrt(3.0) / 2.0;
		n++;
	}
	
	std::cout << "Finding connections/TZP/basement for cells ... " << std::endl;

	//create initial neighbors and connections:


	// Print out absolute positions of nodes (cells) for .dot (graphviz) format
	// double sfact = 4.0;
	// for (int i = 0; i < (*all_cells).size(); i++)
	// {
	// 	// b [pos="100,100"]
	// 	std::cout << i << " [pos=\"" <<  sfact*(*all_cells)[i]->position[0] <<"," << sfact*(*all_cells)[i]->position[1] << "\"]" << std::endl;
	// }
	for (int i = 0; i < (*all_cells).size(); i++)
	{
//		find_connections((*all_cells)[i]);  //rwh: calc all-to-all
		find_connections(i, (*all_cells)[i]);  //rwh: don't calc all-to-all

		find_TZP((*all_cells)[i]);
		find_basement((*all_cells)[i]);
	}


	// int cellid = 0;
 //    std::cout << "cellid= " << cellid << std::endl;
	// for (int idx=0; idx < cell_id_nbrs[cellid].size(); idx++)
 //        std::cout << cell_id_nbrs[cellid][idx] << std::endl;
	// cellid = 84;
 //    std::cout << "cellid= " << cellid << std::endl;
	// for (int idx=0; idx < cell_id_nbrs[cellid].size(); idx++)
 //        std::cout << cell_id_nbrs[cellid][idx] << std::endl;

    //rwh - iterate thru cell_id_nbrs
    int total_count = 0;
    for ( auto it = cell_id_nbrs.begin(); it != cell_id_nbrs.end(); ++it )
    {
//    	std::cout << "cell ID= " << it->first << std::endl;
    	total_count++;
    	//rwh - print out vector of neighbors (as cell ptrs)
//		for (int idx=0; idx < it->second.size(); idx++)
//        	std::cout << "   " << it->second[idx] << std::endl;

//    	std::cout << " " << it->first << ":" << it->second;
    }
    std::cout << "  total (parent) cell (ID) count in cell_id_nbrs = " << total_count << std::endl;


	std::cout << "done!" << std::endl;
	// make a plot 

	PhysiCell_SVG_options.length_bar = 200;
	SVG_plot("initial.svg", microenvironment, 0.0, 0.0, oocyte_coloring_function);

	/*
	//initialize data files
	std::ofstream ofs;
	ofs.open("oocyte_data.txt", std::ofstream::out | std::ofstream::trunc);

	ofs << "time, x, y, z, volume" << std::endl;

	ofs.close();

	ofs.open("granulosa_data.txt", std::ofstream::out | std::ofstream::trunc);

	ofs << "time, x, y, z, volume" << std::endl;;

	ofs.close();
	*/
	return;
}

//-------------------------- Makeing connections ----------------------------------------------

//rwh ------- NOTE:  THIS FUNCTION NOT USED NOW ---------
void add_connection( Cell* pCell_1, Cell* pCell_2 )
{
    if (pCell_1->ID != pCell_2->ID)  //rwh: no longer needed if I don't do N^2 (all-to-all)
    {    
//rwh        pCell_1->custom_data["connections"[pCell_1->ID]]=pCell_2->ID;
//rwh        pCell_2->custom_data["connections"[pCell_2->ID]]=pCell_1->ID;  
        pCell_1->custom_data.user_int_vector["connections"][pCell_1->ID] = pCell_2->ID;
        pCell_2->custom_data.user_int_vector["connections"][pCell_2->ID] = pCell_1->ID;  
    }
    return;
}

void add_TZP_connection( Cell* pCell_1, Cell* pCell_2 )
{
    if (pCell_1->ID != pCell_2->ID)
    {
    	// pCell_1->custom_data["TZP_status"[pCell_2->ID]] = pCell_2->ID;
	    // pCell_2->custom_data["TZP_status"[pCell_1->ID]] = pCell_1->ID;
    	pCell_1->custom_data.user_int_vector["TZP_status"][pCell_2->ID] = pCell_2->ID;
	    pCell_2->custom_data.user_int_vector["TZP_status"][pCell_1->ID] = pCell_1->ID;
    }
    return;
}

void add_basement( Cell* pCell_1 )
{
    pCell_1->custom_data["basement"] = 1;
    return;
}

//-------------------------- Finding connections and neighbors ----------------------------------------------
void find_connections(int idx, Cell* pCell)  //slow way to do this and doesn't allow reconnects
{
//rwh    static double distance = 0;
    double distance = 0.0;
    double diff = 0.0;
    static double turgidity = 1.0;
    static double connection_length = 1.0;
//    static std::vector<double> disp(3,0.0);

//    std::cout << "--------- find_connections:\n";
    //rwh: want to avoid checking myself/pCell
//    for (int i=0; i< (*all_cells).size(); i++)   // loop over all cells
    for (int i=idx+1; i< (*all_cells).size(); i++)   //rwh: "+1" - only check the higher IDs
    {
        distance=0.0;
        for( int j = 0 ; j < 3 ; j++ )
        {
//rwh            disp[j] = (*all_cells)[i]->position[j]-pCell->position[j];
//            distance += disp[j] * disp[j];
            diff = (*all_cells)[i]->position[j] - pCell->position[j];
            distance += diff * diff;
        }
        distance = sqrt(distance);
    
    	//rwh: check if pCell is connected to 
        bool condition1 = fabs(distance) < (pCell->phenotype.geometry.radius +
        		(*all_cells)[i]->phenotype.geometry.radius) *
        		 turgidity + connection_length;

        if (condition1)
        {
//            std::cout << "connect: " << pCell->ID << " -- " << (*all_cells)[i]->ID << std::endl;
            //rwh: for dot/graphviz
            //std::cout << pCell->ID << " -- " << (*all_cells)[i]->ID << std::endl;
            cell_id_nbrs[pCell->ID].push_back((*all_cells)[i]);  //rwh - fill new unordered_map (but never use, yet)
        }

        //original
        // if( fabs(distance) < (pCell->phenotype.geometry.radius +
        // 		(*all_cells)[i]->phenotype.geometry.radius) *
        // 		 turgidity + connection_length 
        if( condition1  &&  (pCell->custom_data.user_int_vector["connections"][(*all_cells)[i]->ID] != -3) )
        {
            add_connection(pCell,(*all_cells)[i]);

//            cell_id_nbrs[pCell->ID].push_back((*all_cells)[i]);  //rwh - using new unordered_map
        }
        else
        {
            pCell->custom_data.user_int_vector["connections"][(*all_cells)[i]->ID] = -3;
            (*all_cells)[i]->custom_data.user_int_vector["connections"][pCell->ID] = -3;
        }
    }
    return;
}

void find_basement(Cell* pCell)//slow way to do this and doesn't allow reconnects
{
	//rwh
	static int initial_vol_index = pCell->custom_data.find_variable_index( "initial volume" );

//	static double gran_radius = pow((3 / (4 * PI)*(pCell->custom_data["initial volume"])), (1 / 3));
//	static double gran_radius = pow((3. / (4. * PI)*(pCell->custom_data["initial volume"])), (1. / 3.));
	static double gran_radius = pow((3. / (4. * PI)*(pCell->custom_data[initial_vol_index])), (1. / 3.));

	static double foll_radius = 60.0;
    static std::vector<double> disp (3, 0.0);
    static double distance = 0.0;

    for (int i=0; i<3; i++)
	    {disp[i] = 0.0 - pCell->position[i];}

    distance = sqrt(disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2]);
    if(distance > (foll_radius-gran_radius))
    {
        add_basement(pCell);
    }
    return;
    
}
void find_TZP(Cell* pCell)//slow way to do this and doesn't allow reconnects
{
    static double distance = 0;
    static double turgidity=1;
    static double TZP_length=10;
    static double min_TZP=.5;
    static std::vector<double> disp(3,0.0);

    if(pCell->type == 1)  // oocyte
    {
        for (int i=0; i< (*all_cells).size(); i++)
        {
            for( int j = 0 ; j < 3 ; j++ )
            {
                disp[j] = (*all_cells)[i]->position[j] - pCell->position[j];
                distance += disp[j] * disp[j];
            }
            distance = sqrt(distance);
            
            if( fabs(distance) < (pCell->phenotype.geometry.radius +
            		(*all_cells)[i]->phenotype.geometry.radius) * turgidity + TZP_length)
            {
                add_TZP_connection(pCell,(*all_cells)[i]);
            }
            else// "disconnect: permanently breaks connection"
            {
                
            }
            distance=0;
        }
    }
    else
    {
    
    }
    return;
}
//-------------------------- Updating Velocity ----------------------------------------------
// take in cell pointer go through other cells if connected calculate dV_x, dV_y, dV_z, do second order adams bashforth to get V_x, V_y, V_z from each connection sum them and hope that is sufficent

void custom_update_velocity_adhesion( Cell* pCell, Phenotype& phenotype, double dt)
{
    // static double connection_length=0;
    // static double prev_velocity[1000][3];
    // static double net_net_force[1000][3];
    // static double net_repulsion_force[1000][3];
    
    // static double force_magnitude=0;
    // static double theta=0;
    // static double phi=0;
    // //std::vector< std::vector<double> > all_the_forces(3, std::vector<double>(100,0));
    
    // static std::vector<double> forcetome (3, 0.0);
    // static std::vector<double> normforce (3, 0.0);
    // static std::vector<double> net_dv (3, 0.0);
    // static std::vector<double> force (3, 0.0);
    // static std::vector<double> dv (3, 0.0);
    // static std::vector<double> dv_net (3, 0.0);
    // static std::vector<double> disp (3, 0.0);
    // //std::vector<double> distance (3, 0.0);
    // static double distance=0;
  
    // static double spring_k=200;
    // static double rel_mass=pCell->phenotype.volume.total;


    double connection_length=0;
    double prev_velocity[1000][3];
    double net_net_force[1000][3];
    double net_repulsion_force[1000][3];
    
    double force_magnitude=0;
    double theta=0;
    double phi=0;
    //std::vector< std::vector<double> > all_the_forces(3, std::vector<double>(100,0));
    
    std::vector<double> forcetome (3, 0.0);
    std::vector<double> normforce (3, 0.0);
    std::vector<double> net_dv (3, 0.0);
    std::vector<double> force (3, 0.0);
    std::vector<double> dv (3, 0.0);
    std::vector<double> dv_net (3, 0.0);
    std::vector<double> disp (3, 0.0);
    //std::vector<double> distance (3, 0.0);
    double distance=0;
  
    double spring_k=200;
    double rel_mass=pCell->phenotype.volume.total;
    

    //rwh: Is the goal to find (only) all cells connected to pCell and recompute pCell's velocity accordingly?

    for(int k; k<(*all_cells).size(); k++)
    {
//        if(pCell->custom_data["connections"[(*all_cells)[k]->ID]] > -1 )  //rwh - wth?
//        if(pCell->custom_data.user_int_vector["connections"][(*all_cells)[k]->ID] >= 0 )  //rwh - wth?
        if ( (pCell->custom_data.user_int_vector["connections"][(*all_cells)[k]->ID] >= 0 ) && 
        	 (pCell != (*all_cells)[k]) )
        {
//rwh            static double r = pCell->phenotype.geometry.radius*2; 
            static double diam = pCell->phenotype.geometry.radius*2;  //rwh: really diameter? ("*2")
            distance = dist( (*all_cells)[k]->position, pCell->position);
            
//rwh            if(fabs(distance)>r)
            if(fabs(distance) > diam)
            {
                for (int i=0;i<3;i++)
                	{disp[i] = (*all_cells)[k]->position[i] - pCell->position[i];}
                
                force_magnitude = distance/3419 * spring_k/rel_mass;
                normforce = normalize(disp);
                forcetome = operator*(force_magnitude,normforce);  //rwh: can just do f*n 
                net_dv += forcetome;
            }
        }
    }
    
    prev_velocity[pCell->ID][0]=pCell->velocity[0];
    prev_velocity[pCell->ID][1]=pCell->velocity[1];
    prev_velocity[pCell->ID][2]=pCell->velocity[2];
    
    
    pCell->velocity[0]=prev_velocity[pCell->ID][0]+(net_dv[0]*dt);
    pCell->velocity[1]=prev_velocity[pCell->ID][1]+(net_dv[1]*dt);
    pCell->velocity[2]=prev_velocity[pCell->ID][2]+(net_dv[2]*dt);
	/*
    std::ofstream ofs;
    ofs.open ("spring_k_data.txt", std::ofstream::out | std::ofstream::app);
    
    ofs << sqrt(net_dv[0]*net_dv[0]+net_dv[1]*net_dv[1]+net_dv[2]*net_dv[2]) <<std::endl;
    
    ofs.close();
	*/
    for (int i=0;i<3;i++)
	    {net_dv[i] = 0;}
    return;
}
void custom_update_velocity_from_TZP_adhesion( Cell* pCell, Phenotype& phenotype, double dt)
{
    // static double antrum_length=2.5;
    // static double prev_velocity[1000][3];
    // static double net_net_force[1000][3];
    // static double net_repulsion_force[1000][3];
    
    // static double force_magnitude=0;
    // //static double theta=0;
    // //static double phi=0;
    // //std::vector< std::vector<double> > all_the_forces(3, std::vector<double>(100,0));
    
    // static std::vector<double> forcetome (3, 0.0);
    // static std::vector<double> normforce (3, 0.0);
    // static std::vector<double> net_oocyte_dv (3, 0.0);
    // static std::vector<double> force (3, 0.0);
    // static std::vector<double> dv (3, 0.0);
    // static std::vector<double> dv_net (3, 0.0);
    // static std::vector<double> disp (3, 0.0);
    // //std::vector<double> distance (3, 0.0);
    // static double distance=0;
   
    // static double spring_k=200;
    // static double rel_mass=pCell->phenotype.volume.total;

    //rwh - remove all 'static's
    double antrum_length=2.5;
    double prev_velocity[1000][3];
    double net_net_force[1000][3];
    double net_repulsion_force[1000][3];
    
    double force_magnitude=0;
    //static double theta=0;
    //static double phi=0;
    //std::vector< std::vector<double> > all_the_forces(3, std::vector<double>(100,0));
    
    std::vector<double> forcetome (3, 0.0);
    std::vector<double> normforce (3, 0.0);
    std::vector<double> net_oocyte_dv (3, 0.0);
    std::vector<double> force (3, 0.0);
    std::vector<double> dv (3, 0.0);
    std::vector<double> dv_net (3, 0.0);
    std::vector<double> disp (3, 0.0);
    //std::vector<double> distance (3, 0.0);
    double distance=0;
   
    double spring_k=200;
    double rel_mass=pCell->phenotype.volume.total;
    
    for(int k; k<(*all_cells).size(); k++)
    {
//        if(pCell->custom_data["TZP_status"[(*all_cells)[k]->ID]]>-1 )
        if(pCell->custom_data.user_int_vector["TZP_status"][(*all_cells)[k]->ID] > -1 )
        {
            static double r= (*all_cells)[k]->phenotype.geometry.radius+antrum_length;
            distance=dist((*all_cells)[k]->position,pCell->position);
            if(fabs(distance)>r)
            {
                for (int i=0;i<3;i++)
                {disp[i] = (*all_cells)[k]->position[i]-pCell->position[i];
                 
                }
                
                force_magnitude = distance/3419*spring_k/rel_mass;
                normforce = normalize(disp);
                forcetome = operator*(force_magnitude,normforce);
                net_oocyte_dv += forcetome;
				/*
                std::ofstream ofs;
                ofs.open ("tzp_cell.txt", std::ofstream::out | std::ofstream::app);
                
                ofs << time_called<<", "<< pCell->position[0]<<", "<<pCell->position[1]<< ", "<< pCell->position[2]<<", "<<pCell->phenotype.geometry.radius<<std::endl;
                
                ofs.close();
                */
            }
        }
    }
    
    prev_velocity[pCell->ID][0] = pCell->velocity[0];
    prev_velocity[pCell->ID][1] = pCell->velocity[1];
    prev_velocity[pCell->ID][2] = pCell->velocity[2];
    
    pCell->velocity[0] = prev_velocity[pCell->ID][0]+(net_oocyte_dv[0]*dt);
    pCell->velocity[1] = prev_velocity[pCell->ID][1]+(net_oocyte_dv[1]*dt);
    pCell->velocity[2] = prev_velocity[pCell->ID][2]+(net_oocyte_dv[2]*dt);
    /*
    std::ofstream ofs;
    ofs.open ("TZP_data.txt", std::ofstream::out | std::ofstream::app);
    
    ofs << sqrt(net_oocyte_dv[0]*net_oocyte_dv[0]+net_oocyte_dv[1]*net_oocyte_dv[1]+net_oocyte_dv[2]*net_oocyte_dv[2]) <<std::endl;
    
    ofs.close();
	*/
    for (int i=0;i<3;i++)
	    {net_oocyte_dv[i] = 0;}
    
    return;
}
void custom_update_velocity_repulsion( Cell* pCell, Phenotype& phenotype, double dt)
{
    
    // static double connection_length=0;
    // static double prev_velocity[1000][3];
    // static double net_net_force[1000][3];
    // static double net_repulsion_force[1000][3];
    
    // static double force_magnitude=0;
    // static double theta=0;
    // static double phi=0;
    // //std::vector< std::vector<double> > all_the_forces(3, std::vector<double>(100,0));
    
    // static std::vector<double> forcetome (3, 0.0);
    // static std::vector<double> normforce (3, 0.0);
    // static std::vector<double> net_dv (3, 0.0);
    // static std::vector<double> force (3, 0.0);
    // static std::vector<double> dv (3, 0.0);
    // static std::vector<double> dv_net (3, 0.0);
    // static std::vector<double> disp (3, 0.0);
    // //std::vector<double> distance (3, 0.0);
    // static double distance=0;
    
    // static double spring_k=500;
    // static double rel_mass=pCell->phenotype.volume.total;


	// rwh NOTE: remove all static declarations!
    double connection_length=0;
    double prev_velocity[1000][3];
    double net_net_force[1000][3];
    double net_repulsion_force[1000][3];
    
    double force_magnitude=0;
    double theta=0;
    double phi=0;
    //std::vector< std::vector<double> > all_the_forces(3, std::vector<double>(100,0));
    
    std::vector<double> forcetome (3, 0.0);
    std::vector<double> normforce (3, 0.0);
    std::vector<double> net_dv (3, 0.0);
    std::vector<double> force (3, 0.0);
    std::vector<double> dv (3, 0.0);
    std::vector<double> dv_net (3, 0.0);
    std::vector<double> disp (3, 0.0);
    //std::vector<double> distance (3, 0.0);
    double distance=0;
    
    double spring_k=500;
    double rel_mass=pCell->phenotype.volume.total;
    
    for(int k; k<(*all_cells).size(); k++)
    {
//rwh        if(pCell->custom_data.user_int_vector["connections"][(*all_cells)[k]->ID] > -1 )
        if ( (pCell->custom_data.user_int_vector["connections"][(*all_cells)[k]->ID] > -1 ) &&
        	 (pCell != (*all_cells)[k]) )
        {
            static double r = pCell->phenotype.geometry.radius+(*all_cells)[k]->phenotype.geometry.radius;
            distance = dist((*all_cells)[k]->position,pCell->position);
            
            if(fabs(distance)<r)
            {
                if ( (*all_cells)[k] == pCell)
                {
                	std::cout << __FUNCTION__ << ": ACK! computing distance from me to me" << "\n";
                	assert(false);
                }
                for (int i=0;i<3;i++)
                {disp[i] = (*all_cells)[k]->position[i] - pCell->position[i];}
                
                force_magnitude = distance/3419 * spring_k/rel_mass;
                normforce = normalize(disp);
                forcetome = operator*(force_magnitude,normforce);
                net_dv -= forcetome;//note the freaking minus sign
            }
        }
    }
    
    prev_velocity[pCell->ID][0]=pCell->velocity[0];
    prev_velocity[pCell->ID][1]=pCell->velocity[1];
    prev_velocity[pCell->ID][2]=pCell->velocity[2];
    
    
    pCell->velocity[0]=prev_velocity[pCell->ID][0]+(net_dv[0]*dt);
    pCell->velocity[1]=prev_velocity[pCell->ID][1]+(net_dv[1]*dt);
    pCell->velocity[2]=prev_velocity[pCell->ID][2]+(net_dv[2]*dt);
    
    for (int i=0;i<3;i++)
    {net_dv[i] = 0;}
    return;
}
void custom_update_velocity_from_TZP_repulsion( Cell* pCell, Phenotype& phenotype, double dt)
{
	//rwh - remove static decls
    // static double connection_length=0;
    // static double prev_velocity[1000][3];
    // static double net_net_force[1000][3];
    // static double net_repulsion_force[1000][3];
    
    // static double force_magnitude=0;
    // //static double theta=0;
    // //static double phi=0;
    // //std::vector< std::vector<double> > all_the_forces(3, std::vector<double>(100,0));
    
    // static std::vector<double> forcetome (3, 0.0);
    // static std::vector<double> normforce (3, 0.0);
    // static std::vector<double> net_oocyte_dv (3, 0.0);
    // static std::vector<double> force (3, 0.0);
    // static std::vector<double> dv (3, 0.0);
    // static std::vector<double> dv_net (3, 0.0);
    // static std::vector<double> disp (3, 0.0);
    // //std::vector<double> distance (3, 0.0);
    // static double distance=0;
    
    // static double spring_k=200;
    // static double rel_mass=pCell->phenotype.volume.total;


    double connection_length=0;
    double prev_velocity[1000][3];
    double net_net_force[1000][3];
    double net_repulsion_force[1000][3];
    
    double force_magnitude=0;
    //static double theta=0;
    //static double phi=0;
    //std::vector< std::vector<double> > all_the_forces(3, std::vector<double>(100,0));
    
    std::vector<double> forcetome (3, 0.0);
    std::vector<double> normforce (3, 0.0);
    std::vector<double> net_oocyte_dv (3, 0.0);
    std::vector<double> force (3, 0.0);
    std::vector<double> dv (3, 0.0);
    std::vector<double> dv_net (3, 0.0);
    std::vector<double> disp (3, 0.0);
    //std::vector<double> distance (3, 0.0);
    double distance=0;
    
    double spring_k=200;
    double rel_mass=pCell->phenotype.volume.total;
    
    for(int k; k<(*all_cells).size(); k++)
    {
//        if(pCell->custom_data["TZP_status"[(*all_cells)[k]->ID]] > -1.0 )
//        if(pCell->custom_data.user_int_vector["TZP_status"][(*all_cells)[k]->ID] > -1.0 )
        if ( (pCell->custom_data.user_int_vector["TZP_status"][(*all_cells)[k]->ID] > -1 ) &&
        	 (pCell != (*all_cells)[k]) )
        {
            static double r= (*all_cells)[k]->phenotype.geometry.radius;
            distance=dist((*all_cells)[k]->position,pCell->position);
            if(fabs(distance)<r)
            {
                for (int i=0;i<3;i++)
                	{disp[i] = (*all_cells)[k]->position[i]-pCell->position[i];
                }
                
                force_magnitude=distance/3419*spring_k/rel_mass;
                normforce=normalize(disp);
                forcetome=operator*(force_magnitude,normforce);
                net_oocyte_dv-=forcetome;
            }
        }
    }
    
    prev_velocity[pCell->ID][0]=pCell->velocity[0];
    prev_velocity[pCell->ID][1]=pCell->velocity[1];
    prev_velocity[pCell->ID][2]=pCell->velocity[2];
    
    pCell->velocity[0]=prev_velocity[pCell->ID][0]+(net_oocyte_dv[0]*dt);
    pCell->velocity[1]=prev_velocity[pCell->ID][1]+(net_oocyte_dv[1]*dt);
    pCell->velocity[2]=prev_velocity[pCell->ID][2]+(net_oocyte_dv[2]*dt);
    
    for (int i=0;i<3;i++)
    	{net_oocyte_dv[i] = 0;}
    
    return;
}
void custom_update_velocity_basement( Cell* pCell, Phenotype& phenotype, double dt)
{
	// static double foll_radius = 60;
 //    static double connection_length=0;
 //    static double prev_velocity[1000][3];
 //    static double net_net_force[1000][3];
 //    static double net_repulsion_force[1000][3];
    
 //    static double force_magnitude=0;
   
 //    //std::vector< std::vector<double> > all_the_forces(3, std::vector<double>(100,0));
    
 //    static std::vector<double> forcetome (3, 0.0);
 //    static std::vector<double> normforce (3, 0.0);
 //    static std::vector<double> net_dv (3, 0.0);
 //    static std::vector<double> force (3, 0.0);
 //    static std::vector<double> dv (3, 0.0);
 //    static std::vector<double> dv_net (3, 0.0);
 //    static std::vector<double> disp (3, 0.0);
 //    //std::vector<double> distance (3, 0.0);
 //    static double distance=0;
    
 //    static double basement_k=200;
 //    static double rel_mass=pCell->phenotype.volume.total;

    
	double foll_radius = 60;
    double connection_length=0;
    double prev_velocity[1000][3];
    double net_net_force[1000][3];
    double net_repulsion_force[1000][3];
    
    double force_magnitude=0;
   
    //std::vector< std::vector<double> > all_the_forces(3, std::vector<double>(100,0));
    
    std::vector<double> forcetome (3, 0.0);
    std::vector<double> normforce (3, 0.0);
    std::vector<double> net_dv (3, 0.0);
    std::vector<double> force (3, 0.0);
    std::vector<double> dv (3, 0.0);
    std::vector<double> dv_net (3, 0.0);
    std::vector<double> disp (3, 0.0);
    //std::vector<double> distance (3, 0.0);
    double distance=0;
    
    double basement_k=200;
    double rel_mass=pCell->phenotype.volume.total;
    
        if(pCell->custom_data["basement"] > -1 )
        {
                for (int i=0;i<3;i++)
                	{disp[i] = 0-pCell->position[i];}   //rwh - check that pCell not at origin!

                distance=sqrt(disp[0]*disp[0]+disp[1]*disp[1]+disp[2]*disp[2]);
                if (distance == 0.0)
                	{return;}
                force_magnitude=distance*basement_k/rel_mass;
                normforce=normalize(disp);
                forcetome=operator*(force_magnitude,normforce);
                net_dv += forcetome;
               
         if (distance > foll_radius)
         {
             net_dv[0]= -1*net_dv[0];
             net_dv[1]= -1*net_dv[1];
             net_dv[2]= -1*net_dv[2];
         }
        }
        
    prev_velocity[pCell->ID][0]=pCell->velocity[0];
    prev_velocity[pCell->ID][1]=pCell->velocity[1];
    prev_velocity[pCell->ID][2]=pCell->velocity[2];
    
    pCell->velocity[0]=prev_velocity[pCell->ID][0]-(net_dv[0]*dt);
    pCell->velocity[1]=prev_velocity[pCell->ID][1]-(net_dv[1]*dt);
    pCell->velocity[2]=prev_velocity[pCell->ID][2]-(net_dv[2]*dt);
    
    for (int i=0;i<3;i++)
    	{net_dv[i] = 0;}
    return;
}

//-------------------------- microenvironment----------------------------------------------
void setup_microenvironment( void )
{
    // set domain parameters
    default_microenvironment_options.use_oxygen_as_first_field = false;  //rwh??
//    default_microenvironment_options.use_oxygen_as_first_field = true;  //rwh??

    //rwh - read from config file now... well, we don't really do that yet. Ugh.
    default_microenvironment_options.X_range = {-150, 150};
    default_microenvironment_options.Y_range = {-150, 150};
    //rwh
//    default_microenvironment_options.X_range = {-250, 250};
//    default_microenvironment_options.Y_range = {-250, 250};

    //default_microenvironment_options.Z_range = {-150, 150};

    default_microenvironment_options.dx = 1;
    default_microenvironment_options.dy = 1;
    default_microenvironment_options.dz = 1;
    default_microenvironment_options.simulate_2D = true;

    // turn on gradients
    default_microenvironment_options.calculate_gradients = true;

    // add the solute

    //rwh          microenvironment.set_density(0, "glucose", "dimensionless" );

//    microenvironment.add_density( "empty", "mol/L" );
//    microenvironment.diffusion_coefficients[1] = 0;
//    microenvironment.decay_rates[1] = 0;
    microenvironment.set_density( 0, "empty", "mol/L" );  //rwh
    microenvironment.diffusion_coefficients[0] = 0;
    microenvironment.decay_rates[0] = 0;
    

    microenvironment.add_density( "CPA", "dimensionless" );
    //rwh - change CPA density index
//    microenvironment.diffusion_coefficients[2] = 1100;
//    microenvironment.decay_rates[2] = 0;
    microenvironment.diffusion_coefficients[1] = 1100;
    microenvironment.decay_rates[1] = 0;
//    microenvironment.diffusion_coefficients[2] = 100;  //rwh
    /*
     microenvironment.add_density( "salt", "dimensionless" );
     microenvironment.diffusion_coefficients[2] = 1e2;
     microenvironment.decay_rates[2] = 0;
     
     microenvironment.add_density( "solute", "mol/L" );
     microenvironment.diffusion_coefficients[1] = 1100;
     microenvironment.decay_rates[1] = 0;
     */
    default_microenvironment_options.outer_Dirichlet_conditions = true;
    
    default_microenvironment_options.time_units = "sec";
    microenvironment.spatial_units = "microns";
    microenvironment.time_units = "seconds";
    microenvironment.mesh.units = "microns";
    
//    default_microenvironment_options.Dirichlet_condition_vector[2] = 0;
    default_microenvironment_options.Dirichlet_condition_vector[1] = 0;   //rwh
  
    initialize_microenvironment();
    
    return;
}

//-------------------------- SVG Coloring function --------------------------------------------
std::vector<std::string> oocyte_coloring_function( Cell* pCell )
{
	std::string color = "black"; 
	std::vector< std::string > output( 4 , color ); 
	
	// black cells if necrotic 
	if( pCell->phenotype.death.dead == true )
		{ return output; }

	output[3] = "none"; // no nuclear outline color 

	if( pCell->type == 1 )
		{ color = "blue"; }

	if( pCell->type == 2 )
		{ color = "limegreen"; }

	output[0] = color; 
	output[2] = color; 
	
	return output; 
}

//---------------------------------------OOCYTE CELL RULES----------------------------------------------------

void oocyte_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
//	std::cout << __FUNCTION__ << ":"<<__LINE__ << "---------- dt= " << dt << "\n";  // = 0.1
    oocyte_volume_function( pCell, phenotype, dt);
//rwh    std::cout<<"oocyte volume: "<< pCell->phenotype.volume.total<<std::endl;

#define dbg_oo_cell_rule
#ifdef dbg_oo_cell_rule
    std::cout<<__FUNCTION__<<":"<<__LINE__ <<": volume: "<< pCell->phenotype.volume.total << std::endl;
#endif
    //custom_update_velocity_from_TZP(pCell, phenotype, dt);
	/*
    std::ofstream ofs;
    ofs.open ("oocyte_data.txt", std::ofstream::out | std::ofstream::app);
    
    ofs << time_called<<", "<< pCell->position[0]<<", "<<pCell->position[1]<< ", "<< pCell->position[2]<<", "<<pCell->phenotype.volume.total<<std::endl;
    
    ofs.close();
    */
    time_called += 0.1;
    return;
}

void oocyte_cell_phenotype_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
    return;
}
/*  GRANULOSA CELL RULES */
// Needs to be thread-safe!
void gran_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
   // std::cout<<"position: "<< pCell->position[0]<<","<< pCell->position[1]<<","<< pCell->position[2]<<std::endl;
    //std::cout<<" monkey cell id"<< pCell->ID<<std::endl;
    oocyte_volume_function(pCell,phenotype,dt);
//#define dbg_gran_cell_rule
#undef dbg_gran_cell_rule
#ifdef dbg_gran_cell_rule
    std::cout<<__FUNCTION__<<":"<<__LINE__ << ": volume: "<< pCell->phenotype.volume.total << std::endl;
#endif
    //std::cout<<"gran volume: "<< pCell->phenotype.volume.total<<std::endl;
    
    custom_update_velocity_adhesion( pCell, phenotype, dt);
    custom_update_velocity_from_TZP_adhesion( pCell,phenotype,dt);
    custom_update_velocity_repulsion( pCell, phenotype, dt);
    custom_update_velocity_from_TZP_repulsion( pCell,phenotype,dt);
    custom_update_velocity_basement( pCell, phenotype, dt);
    //std::cout<<"volume: "<< pCell->phenotype.volume.total<<std::endl;
   // std::cout<<"position: "<< pCell->position[0]<<","<< pCell->position[1]<<","<< pCell->position[2]<<std::endl;
   /*
	std::ofstream ofs;
    ofs.open ("granulosa_data.txt", std::ofstream::out | std::ofstream::app);
    
    ofs << time_called<<", "<< pCell->position[0]<<", "<<pCell->position[1]<< ", "<< pCell->position[2]<<", "<<pCell->phenotype.volume.total<<std::endl;
    
    ofs.close();
	*/
    return;
}

void gran_cell_phenotype_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
    
    
    return;
}

//-------------------------- volume functions ----------------------------------------------
double sample_cell_boundary(Cell* pCell, Phenotype& phenotype, int solute)
{
	std::vector<double> boundary_point(3, 0);

	static int solute_density_index = 0;  //rwh

//rwh	static double solute_average_at_boundary;
	double solute_average_at_boundary;
	// rwh
	solute_average_at_boundary = 0.0;

	//rwh
	static double delta = 2.0 * PI / 8.0;
	for (int i = 0; i < 8; i++)
	{
		//rwh
//		boundary_point[0] = (pCell->position[0]) + pCell->phenotype.geometry.radius*cos(i * 45 * PI / 180.0);
//		boundary_point[1] = (pCell->position[1]) + pCell->phenotype.geometry.radius*sin(i * 45 * PI / 180.0);
		boundary_point[0] = (pCell->position[0]) + pCell->phenotype.geometry.radius * cos(i * delta);
		boundary_point[1] = (pCell->position[1]) + pCell->phenotype.geometry.radius * sin(i * delta);

//rwh		solute_average_at_boundary += microenvironment.nearest_density_vector(boundary_point)[solute];
		solute_average_at_boundary += microenvironment.nearest_density_vector(boundary_point)[solute_density_index];
	}

//rwh - comment out
	// for (int i = 0; i < 8; i++)
	// {
	// 	boundary_point[0] = (pCell->position[0]) + pCell->phenotype.geometry.radius*cos(i * 45 * PI / 180.0);
	// 	boundary_point[2] = (pCell->position[2]) + pCell->phenotype.geometry.radius*sin(i * 45 * PI / 180.0);
	// 	solute_average_at_boundary += microenvironment.nearest_density_vector(boundary_point)[solute];
	// }
	// solute_average_at_boundary = solute_average_at_boundary / 16;//1.5;
	solute_average_at_boundary = solute_average_at_boundary / 8.0;
	return solute_average_at_boundary;
}

void oocyte_volume_function( Cell* pCell, Phenotype& phenotype, double dt)
{
	std::vector<double> center = pCell->position;
	//rwh - drop 'static' on the following
	/*
	static double k1_water=0;
	static double k1_EG=0;
	static double Lp = pCell->custom_data["Lp"];
	static double T = pCell->custom_data["Temperature"];

	static double R1 = pCell->custom_data["R"];
//	static int k= pCell->ID;  //rwh never used?
    */
	double k1_water=0;
	double k1_EG=0;
	double Lp = pCell->custom_data["Lp"];
	double T = pCell->custom_data["Temperature"];
	double R1 = pCell->custom_data["R"];

	// rwh: following line in segfault path (next to last) -- try changing case to match the defined one (next 2 lines)!!
//	pCell->custom_data["ext_salt_concentration"] =.2;//ext_salt_tot/sample_cell_boundary(pCell,phenotype, salt number);
//rwh	pCell->custom_data["Ext_salt_concentration"]=.2;//ext_salt_tot/sample_cell_boundary(pCell,phenotype, salt number);
//	std::cout << "ext_salt_conc=" << pCell->custom_data["Ext_salt_concentration"] << "\n";
//    double Total_external_concentration= pCell->custom_data["ext_salt_concentration"] +sample_cell_boundary(pCell,phenotype, 2);//1.7;
    double Total_external_concentration = pCell->custom_data["Ext_salt_concentration"] +
    										sample_cell_boundary(pCell,phenotype, 2);//1.7;

//rwh - zeros
/*
	std::cout << "Inter_CPA_moles[0]=" << pCell->custom_data["Inter_CPA_moles"[0]] << "\n";
	std::cout << "Inter_CPA_moles[1]=" << pCell->custom_data["Inter_CPA_moles"[1]] << "\n";
	std::cout << "Inter_CPA_moles[2]=" << pCell->custom_data["Inter_CPA_moles"[2]] << "\n";
	std::cout << "Inter_CPA_moles[3]=" << pCell->custom_data["Inter_CPA_moles"[3]] << "\n";
	*/
	

	// double Total_internal_concentration = (pCell->custom_data["Inter_CPA_moles"[0]] + 
	// 	pCell->custom_data["Inter_CPA_moles"[1]] + 
	// 	pCell->custom_data["Inter_CPA_moles"[2]] + 
	// 	pCell->custom_data["Inter_sucrose_moles"] + 
	// 	pCell->custom_data["Inter_salt_moles"]) / pCell->custom_data["Water_volume"];

	double Total_internal_concentration = (pCell->custom_data.user_vector["Inter_CPA_moles"][0] + 
		pCell->custom_data.user_vector["Inter_CPA_moles"][1] + 
		pCell->custom_data.user_vector["Inter_CPA_moles"][2] + 
		pCell->custom_data["Inter_sucrose_moles"] + 
		pCell->custom_data["Inter_salt_moles"]) / pCell->custom_data["Water_volume"];
    
    //rwh - check for negative
//	static int initial_vol_index = pCell->custom_data.find_variable_index( "initial volume" );
    double init_vol = pCell->custom_data["initial volume"];
//    double init_vol = pCell->custom_data[initial_vol_index];
    double water_vol = pCell->custom_data["Water_volume"];
    double vBarEG = V_bar_EG;
//    double interCPAMoles0 = pCell->custom_data["Inter_CPA_moles"[0]];
    double interCPAMoles0 = pCell->custom_data.user_vector["Inter_CPA_moles"][0];

//#define dbg_vol_fn    
#undef dbg_vol_fn    
#ifdef dbg_vol_fn    
    std::cout<<" init_vol, etc= " <<init_vol<<","<<water_vol<<","<<vBarEG<<","<<interCPAMoles0 <<"\n";
#endif

//rwh    pCell->phenotype.volume.total= pCell->custom_data["initial volume"]*
//    	((pCell->custom_data["Water_volume"] + V_bar_EG * pCell->custom_data["Inter_CPA_moles"[0]]));
    pCell->phenotype.volume.total = init_vol * (water_vol + vBarEG * interCPAMoles0);
//    pCell->phenotype.volume.total = init_vol;

//rwh    pCell->phenotype.volume.total += 0.02 * init_vol;
//    std::cout<<" => "<<pCell->phenotype.volume.total<<"\n";

    //rwh: reminder of user/custom vector var names
	// gran_cell.custom_data.add_vector("Ps_vec", Ps_Vector);  //rwh: change name "Ps" --> "Ps_vec"
	// gran_cell.custom_data.add_vector("CPA_fluxes", CPA_fluxes);
	// gran_cell.custom_data.add_vector("Inter_CPA_moles", Inter_CPA_moles);
	// gran_cell.custom_data.add_vector("Ext_CPA_concentration", Ext_CPA_concentration);
	// gran_cell.custom_data.add_vector("Prev_Inter_CPA_moles", Prev_Inter_CPA_moles);

    //toxcity calculation-----------------------------
//	double cost_function_concentration = pCell->custom_data["Inter_CPA_moles"[0]] / pCell->custom_data["Water_volume"];
	double cost_function_concentration = pCell->custom_data.user_vector["Inter_CPA_moles"][0] / pCell->custom_data["Water_volume"];

#ifdef dbg_vol_fn    
    std::cout<<__FUNCTION__<<":"<<__LINE__<<" => "<<pCell-> custom_data.user_vector["Inter_CPA_moles"][0] <<"\n";
#endif   

	pCell->custom_data["Toxicity"] = pCell->custom_data["Toxicity"] + dt*pow(cost_function_concentration,1.6);
	//------------------------------------------------------
        //fluxes
		pCell->custom_data["Water_flux"] = -Lp*R1*T*(Total_external_concentration-Total_internal_concentration);

		for (int idx = 0; idx<4; idx++)
		{
			// pCell->custom_data["CPA_fluxes"[i]] = pCell->custom_data["Ps"[i]] * ((sample_cell_boundary(pCell, phenotype, 2) - 
			// 	(pCell->custom_data["Inter_CPA_moles"[0]] / pCell->custom_data["Water_volume"]) ));

			//rwh: is this what was intended?
			pCell->custom_data.user_vector["CPA_fluxes"][idx] = pCell->custom_data.user_vector["Ps_vec"][idx] * 
				((sample_cell_boundary(pCell, phenotype, 2) - (pCell->custom_data.user_vector["Inter_CPA_moles"][0] / 
					pCell->custom_data["Water_volume"]) ));
				
		}
		
        k1_water = dt * pCell->custom_data["Water_flux"];// k1w= dt* water_flux
        
		pCell->custom_data["Water_volume"] = pCell->custom_data["Prev_Water_volume"] + k1_water;//change in volume= volume change+w_flux*dt//Vn+1=vn+dv*dt
    
//        k1_EG = dt * pCell->custom_data["CPA_fluxes"[0]];
		// pCell->custom_data["Inter_CPA_moles"[0]] = pCell->custom_data["Prev_Inter_CPA_moles"[0]] +k1_EG;
		// pCell->custom_data["Prev_Water_volume"] = pCell->custom_data["Water_volume"];
		// pCell->custom_data["Prev_Inter_CPA_moles"[0]] = pCell->custom_data["Inter_CPA_moles"[0]];

		pCell->custom_data["Prev_Water_volume"] = pCell->custom_data["Water_volume"];

		//rwh: confused - why only use/update 0th element?
        k1_EG = dt * pCell->custom_data.user_vector["CPA_fluxes"][0];
		pCell->custom_data.user_vector["Inter_CPA_moles"][0] = pCell->custom_data.user_vector["Prev_Inter_CPA_moles"][0] + k1_EG;
#ifdef dbg_vol_fn    
		std::cout<<__FUNCTION__<<":"<<__LINE__<<" => "<<pCell-> custom_data.user_vector["Inter_CPA_moles"][0] <<
				", k1_EG=" << k1_EG << "\n";
#endif

		pCell->custom_data.user_vector["Prev_Inter_CPA_moles"][0] = pCell->custom_data.user_vector["Inter_CPA_moles"][0];
      
//   pCell->phenotype.secretion.uptake_rates[2] = pCell->custom_data["CPA_fluxes"[0]];//+1-2*microenvironment.nearest_density_vector(pCell->position)[2]); //
//   pCell->phenotype.secretion.uptake_rates[2] = pCell->custom_data.user_vector["CPA_fluxes"][0];
   pCell->phenotype.secretion.uptake_rates[1] = pCell->custom_data.user_vector["CPA_fluxes"][0];

    //pCell->phenotype.secretion.secretion_rates[2] = (S_flux+1-2*microenvironment.nearest_density_vector(pCell->position)[2]); ////
}

void gran_volume_function( Cell* pCell, Phenotype& phenotype, double dt)
{/*
   // for (int i=1; i< (*all_cells).size();i++)
    //{
        
        std::cout<< "granulosa cell ID: "<< pCell->ID<< std::endl;
        static int k= pCell->ID;
        double N_out=.2;
        double Total_external_concentration=N_out+sample_cell_boundary(pCell,phenotype, 2);//1.7;
        std::cout<< "external C "<< Total_external_concentration << std::endl;
        std::vector<double> center=pCell->position;
        pCell->phenotype.volume.total=pCell->custom_data["initial volume"]*(a[9][k]+V_bar_s*a[0][k]);
        std::cout<<"volume: "<<pCell->phenotype.volume.total<<std::endl;
        //pCell->set_total_volume(pCell->custom_data["initial volume"]*(a[9][pCell->ID]+V_bar_s*a[0][pCell->ID]));
        double cc=a[0][k]/a[9][k];
        a[10][k]=a[10][k]+dt*pow(cc,1.6);
        //fluxes
        W_flux=-(5.333/25.241)*(Total_external_concentration-(a[0][k]+N_mol)/(a[9][k]));
        S_flux=(21/25.241)*(sample_cell_boundary(pCell,phenotype, 2)-(a[0][k]/(a[9][k])));
        N_flux=0;//order of magnitude slower than solute and 2 orders slower than water so approximate as 0
        // k1w= dt* water_flux
        k1w=dt*a[1][k];
        //change in volume= volume change+w_flux*dt//Vn+1=vn+dv*dt
        v_change=a[9][k]+k1w;
        
        k1s=dt*a[5][k];
        s_change=a[0][k]+k1s;
        a[9][k]=v_change;
        std::cout<<"v change: "<<a[9][k]<< std::endl;
        a[0][k]=s_change;
        
        W_fluxiii=W_fluxii;
        W_fluxii=W_fluxi;
        W_fluxi=W_flux;
        
        S_fluxiii=S_fluxii;
        S_fluxii=S_fluxi;
        S_fluxi=S_flux;
        
        
        a[1][k]=W_flux;
        a[2][k]=W_fluxi;
        a[3][k]=W_fluxii;
        a[4][k]=W_fluxiii;
        a[5][k]=S_flux;
        a[6][k]=S_fluxi;
        a[7][k]=S_fluxii;
        a[8][k]=S_fluxiii;
        
      
    pCell->phenotype.secretion.uptake_rates[2] = (S_flux);//+1-2*microenvironment.nearest_density_vector(pCell->position)[2]);///(*all_cells)[i]->phenotype.volume.total; ////
     
   // }
    
    
  */  
	return;
}
//---------------------------------------second order Adams-Bashforth----------------------------------------------------
//---------------------------------------Temperature change----------------------------------------------------

//---------------------------------------Osmoregulation----------------------------------------------------
