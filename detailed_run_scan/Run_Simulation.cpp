#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <chrono>
#include <random>
#include <boost/numeric/odeint.hpp>
#include <stdexcept>
#include <mpi.h>
#include "../code/linterp.h"
#include "../code/Constants.hpp"
#include "../code/Model.hpp"
#include "../code/Simulation_Setups.hpp"

void Read_Setup(std::string filename, std::vector<double> &gene_lengths, std::vector<double> &TSSes, std::vector<int> &Gene_Directions, std::vector<double> &promoters)
{
	std::ifstream f;
	f.open(filename);
	double gene_length, TSS, promoter;
	int direction;
	std::string gene_name;

	while(f >> gene_name >> TSS >> gene_length >> direction >> promoter)
	{
		gene_lengths.push_back(gene_length);
		TSSes.push_back(TSS);
		Gene_Directions.push_back(direction);
		promoters.push_back(promoter);
	}

	f.close();

	return;
}

void Run_Simulation(std::string filename, int plasmid_flag, int promoter_flag, int GQ_on_id, double GQ_on_val, int eR_on_id, double eR_on_val, int topo_id, double topo_val, int eR_off_id, double eR_off_val, int PQS_flag, int world_rank, int start_rank, int restartflag)
{
	double force = 1.0;
	int torque_flag = 0;
	int brute_force_flag = 1;

	sigma_basal = -0.06;

	std::vector<double> gene_lengths, TSSes, promoters;
	std::vector<int> Gene_Directions;
	Read_Setup(filename, gene_lengths, TSSes, Gene_Directions, promoters);

	int clamp0_flag = 1, clamp1_flag = 1;
	double T = 3600.0*2.0;
	int finish_count_limit = -1;
	double topoisomerase = topo_val;
	double rna_degrad = 0.0;
	double barrier_on = 0.0, barrier_off = 0.0; // No barriers needed for this run
	double barrier = 1350.0;
	double nucl_par = 1.0;
	int fileflag = 1;

	int single_pol_II_flag = 0;

	int nucl_file_flag = 1;

//	for(int i = 0; i < promoters.size(); i++)
//	{
//		promoters[i] = k_on_val;
//	}

	if(plasmid_flag == 0)
	{
		clamp0_flag = 1;
		clamp1_flag = 1;

		clamp0_sigma = 0.0;
		clamp1_sigma = clamp0_sigma - (sigma_basal*w0*(3700.0*0.34));

		circular_plasmid_flag = 0;
	}
	else if(plasmid_flag == 1)
	{
		circular_plasmid_flag = 1;

		topoisomerase = 0.0;
	}
	else
	{
		std::cout << "Invalid plasmid parameter" << "\n";
		throw std::invalid_argument("Invalid plasmid parameter");
	}

	if(promoter_flag == 0)
	{
		myong_k_0 = 0.01;
		myong_s_0 = 0.065;
	}
	else if(promoter_flag == 1)
	{
		myong_k_0 = 0.04;
		myong_s_0 = 0.065;
	}
	else if(promoter_flag == 2)
	{
		myong_k_0 = 0.01;
		myong_s_0 = 0.05;
	}

	myong_k_on_G = GQ_on_val / 60.0;
	myong_k_on_eR = eR_on_val / 60.0;
	myong_k_off_eR = eR_off_val / 60.0;

	std::string outputfolder = "outputfiles/RUN_" + std::to_string(plasmid_flag) + "_" + std::to_string(promoter_flag) + "_" + std::to_string(GQ_on_id) + "_" + std::to_string(eR_on_id) + "_" + std::to_string(topo_id) + "_" + std::to_string(eR_off_id) + "_" + std::to_string(PQS_flag);
	std::string inputfolder = "inputfiles/RUN_" + std::to_string(plasmid_flag) + "_" + std::to_string(promoter_flag) + "_" + std::to_string(GQ_on_id) + "_" + std::to_string(eR_on_id) + "_" + std::to_string(topo_id) + "_" + std::to_string(eR_off_id) + "_" + std::to_string(PQS_flag);

	Gillespie_Simulation_Multiple_Genes(force, restartflag, brute_force_flag, torque_flag, gene_lengths, TSSes, Gene_Directions, clamp0_flag, clamp1_flag, T, finish_count_limit, promoters, topoisomerase, rna_degrad, barrier_on, barrier_off, barrier, nucl_par, PQS_flag, fileflag, outputfolder, inputfolder, world_rank + start_rank, single_pol_II_flag, nucl_file_flag);

	return;
}

int main(int argc, char *argv[])
{
	MPI_Init(NULL, NULL);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	std::string filename = argv[1];
	int restartflag = std::stoi(argv[2]);
	int plasmid_flag = std::stoi(argv[3]);
	int promoter_flag = std::stoi(argv[4]);
	int GQ_on_id = std::stoi(argv[5]);
	double GQ_on_val = std::stod(argv[6]);
	int eR_on_id = std::stoi(argv[7]);
	double eR_on_val = std::stod(argv[8]);
	int topo_id = std::stoi(argv[9]);
	double topo_val = std::stod(argv[10]);
	int eR_off_id = std::stoi(argv[11]);
	double eR_off_val = std::stod(argv[12]);
	int PQS_flag = std::stoi(argv[13]);
	int start_rank = std::stoi(argv[14]);

	generator = std::mt19937(std::time(NULL) + plasmid_flag + promoter_flag + GQ_on_id + eR_on_id + topo_id + eR_off_id + PQS_flag + world_rank + start_rank);

	try
	{
		Run_Simulation(filename, plasmid_flag, promoter_flag, GQ_on_id, GQ_on_val, eR_on_id, eR_on_val, topo_id, topo_val, eR_off_id, eR_off_val, PQS_flag, world_rank, start_rank, restartflag);
	}
	catch(const std::invalid_argument &e)
	{
		std::cout << "Exceptional circumstances" << "\n";
		return(0);
	}

	MPI_Finalize();

	return(0);
}
