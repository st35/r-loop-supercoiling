std::mt19937 generator;

void splitstring(std::string line, std::string delim, std::vector<std::string> *l, int popflag)
{
	if(line.length() == 0)
	{
		return;
	}

	if(popflag == 1)
	{
        	line.pop_back();
	}

	if(line.length() == 0)
	{
		return;
	}

        int toklen, pos1 = 0, pos2;
        while(1)
        {
                pos2 = line.find(delim, pos1);
                if(pos2 == std::string::npos)
                {
                        break;
                }
                toklen = pos2 - pos1;
                (*l).push_back(line.substr(pos1, toklen));
                pos1 = pos2 + 1;
        }
        (*l).push_back(line.substr(pos1, std::string::npos));

	return;
}

template<typename T> void Write_Vector_State_To_File(std::ofstream *f, const double t, const std::vector<T> &vec, int num_entries, int shift)
{
	(*f) << t << "\t";
	for(int i = 0; i < num_entries; i++)
	{
		(*f) << vec[i + shift] << "\t";
	}
	(*f) << "\n";
}

void Gillespie_Simulation_Multiple_Genes(double force, int restartflag, int brute_force_flag, int torque_flag, std::vector<double> gene_lengths, std::vector<double> TSSes, const std::vector<int> &Gene_Directions, int clamp0_flag, int clamp1_flag, double T, int finish_count_limit, const std::vector<double> &promoters, double topoisomerase, double rna_degrad, double barrier_on, double barrier_off, double barrier, double nucl_par, int PQS_flag, int fileflag, std::string outputfolder, std::string inputfolder, int world_rank, int single_pol_II_flag, int nucl_file_flag, int probe_gene_bodies_flag = 0, int galactose_switch = -1, double galactose_time = -1.0, double GapR_off_par = 1.0, int brute_gapr_flag = 0, int henikoff_switch_flag = 0, double henikoff_time = -1.0)
{
	InterpMultilinear<2, double> interp_ML = Setup_Interp();
	InterpMultilinear<1, double> interp_ML_Cutoff_0 = Setup_Interp_Cutoffs("../torque_interp/sigma_s.log");
	InterpMultilinear<1, double> interp_ML_Cutoff_1 = Setup_Interp_Cutoffs("../torque_interp/sigma_p.log");

	std::uniform_real_distribution<double> distribution(0.0, 1.0);

	int numgenes = gene_lengths.size();
	if(TSSes.size() != numgenes || Gene_Directions.size() != numgenes || promoters.size() != numgenes)
	{
		std::cout << "Mismatch between gene count and associated vectors" << "\n";
		throw std::invalid_argument("Mismatch between gene count and associated vectors");
	}

	for(int i = 0; i < ((int) TSSes.size()) - 1; i++)
	{
		if(TSSes[i] < TSSes[i + 1])
		{
			std::cout << "TSSes are in the wrong order." << "\n";
			throw std::invalid_argument("TSSes are in the wrong order.");
		}
	}

	std::vector<double> gene_start, gene_end;
	for(int i = 0; i < numgenes; i++)
	{
		TSSes[i] = TSSes[i]*0.34;
		gene_lengths[i] = gene_lengths[i]*0.34;
		if(Gene_Directions[i] == 1)
		{
			gene_start.push_back(TSSes[i]);
			gene_end.push_back(TSSes[i] + gene_lengths[i]);
		}
		else
		{
			gene_start.push_back(TSSes[i] - gene_lengths[i]);
			if(TSSes[i] - gene_lengths[i] < 0.0)
			{
				std::cout << "Possible error in the config file; TSS location not compatible with gene length." << "\n";
				throw std::invalid_argument("Possible error in the config file; TSS location not compatible with gene length.");
			}
			gene_end.push_back(TSSes[i]);
		}
	}
	int overlapflag = 0;
	double overlap = 0.0;
	for(int i = 0; i < gene_start.size(); i++)
	{
		for(int j = i + 1; j < gene_start.size(); j++)
		{
			if(gene_start[i] <= gene_end[j] && gene_start[j] <= gene_end[i] && promoters[i] > 0.0 && promoters[j] > 0.0)
			{
				overlapflag = 1;
				overlap = std::min(gene_end[i], gene_end[j]) - std::max(gene_start[i], gene_start[j]);
			}
		}
	}
	if(overlapflag == 1)
	{
		std::cout << "Gene bodies overlap; code is not configured for scenarios" << "\n";
		throw std::invalid_argument("Gene bodies overlap; code is not configured for scenarios");
	}

	std::vector<double> k_on;
	for(int i = 0; i < numgenes; i++)
	{
//		k_on.push_back((0.5 / 60.0)*promoters[i]);
		k_on.push_back(1.0);
	}
	double k_off = 0.0; // RNAPs don't fall off before crossing the gene end site.
	double k_topo = (1.0 / 60.0)*topoisomerase;
	double k_nucl_on = 12.0*0.1;
	double k_nucl_off = 4.0*0.1*nucl_par;
	double k_rna_degrad = 0.01*rna_degrad;
	double k_prot_barrier_on = (0.5 / 60.0)*barrier_on;
	double k_prot_barrier_off = (0.5 / 60.0)*barrier_off;
	double alpha_marenduzzo = 0.0;
	double GapR_On = (0.5 / 60.0)*0.0;
	double GapR_Off = (0.5 / 60.0)*GapR_off_par*0.0;

	if(torque_flag == 0) // In the case of prokaryotes: no nucleosome binding / unbinding
	{
		k_nucl_on = 0.0;
		k_nucl_off = 0.0;
	}

	int s_R = 0, s_G = 0, s_eR = 0;

	int s_PQS = PQS_flag;

	double k_on_G = myong_k_on_G;
	double k_on_eR = myong_k_on_eR;

	double k_off_R = 1.0;
	double k_off_G = 0.0;
	double k_off_eR = 0.0;

	barrier = 0.34*barrier;
	clamp0 = 0.0;
	if(Gene_Directions[0] == 1)
	{
		clamp1 = TSSes[0] + gene_lengths[0] + barrier;
	}
	else
	{
		clamp1 = TSSes[0] + barrier;
	}

	std::vector<double> GapR_Points;
	if(probe_gene_bodies_flag == 0)
	{
		GapR_Points = Get_GapR_Points(clamp0, clamp1, 100.0*0.34); // 100 bp separation between points where we probe the GapR binding
	}
	else
	{
		GapR_Points = Get_GapR_Points_Gene_Bodies(gene_start, gene_end, gene_lengths);
	}
	std::vector<double> RNAP_Points = Get_GapR_Points(clamp0, clamp1, 100.0*0.34); // 100 bp separation between points for probing RNAP occupancy
	std::vector<double> Nucl = Get_Nucleosome_Array(clamp0, clamp1);
	std::vector<double> Barr = Get_Barrier_Array_Monica(clamp0, clamp1, 1000.0*0.34);
	int orderflag = 0;
	for(int i = 0; i < ((int) Barr.size()) - 1; i++)
	{
		if(Barr[i] < Barr[i + 1])
		{
			orderflag = 1;
		}
	}

	if(orderflag == 1)
	{
		std::cout << "Barrier protein-binding sites are in the wrong order." << "\n";
		throw std::invalid_argument("Barrier protein-binding sites are in the wrong order.");
	}

	std::vector<double> Barr_Phi;
	std::vector<int> Nucl_Status, Alive, Barr_Status, Is_Nucl_Blocked, Is_RNAP_Blocked, Is_TSS_Blocked, Segment_Nucl_Count, GapR_State, Is_GapR_Blocked, Nucl_Type;
	double GQ_spotsigma = 0.0;

	int numevents = 0;
	for(int i = 0; i < numgenes; i++)
	{
		numevents += 1; // Recruitment of RNAPs
	}
	numevents += 1; // Drop-off of RNAPs before end of transcription
	numevents += 1; // Relaxation of DNA supercoiling
	for(int i = 0; i < Nucl.size(); i++)
	{
		numevents += 2; // Binding and unbinding of nucleosomes
		if(torque_flag == 1 && distribution(generator) < (k_nucl_on / (k_nucl_on + k_nucl_off)))
		{
			Nucl_Status.push_back(1);
			Nucl_Type.push_back(1);
		}
		else
		{
			Nucl_Status.push_back(0);
			Nucl_Type.push_back(0);
		}
	}
	for(int i = 0; i < numgenes; i++)
	{
		numevents += 1; // Degradation of RNAs
	}
	for(int i = 0; i < Barr.size(); i++)
	{
		numevents += 2; // Binding and unbinding of proteins at barrier sites
		Barr_Status.push_back(0);
		Barr_Phi.push_back(NAN);
	}
	for(int i = 0; i < GapR_Points.size(); i++)
	{
		numevents += 2; // Binding and unbinding of GapR
		GapR_State.push_back(0);
	}
	numevents += 2; // R-loop formation and dissolution
	numevents += 2; // GQ formation and dissolution
	numevents += 2; // Extended R-loop formation and dissolution

	std::ifstream inputstream;
	int inputindex, inputval0, inputval2;
	double inputval1;
	std::string nanstr;
	if(restartflag == 1)
	{
		inputstream.open(inputfolder + "/nucl_status_" + std::to_string(world_rank) + ".log");
		inputindex = 0;
		while(inputstream >> inputval0)
		{
			Nucl_Status[inputindex] = inputval0;
			inputindex += 1;
		}
		if(inputindex != Nucl_Status.size())
		{
			std::cout << "Error in the nucl_status input file" << "\n";
			throw std::invalid_argument("Error in the nucl_status input file");
		}
		inputstream.close();

		inputstream.open(inputfolder + "/barr_file_" + std::to_string(world_rank) + ".log");
		inputindex = 0;
		while(inputstream >> inputval0 >> nanstr)
		{
			Barr_Status[inputindex] = inputval0;
			if(inputval0 == 0)
			{
				Barr_Phi[inputindex] = NAN;
			}
			else
			{
				Barr_Phi[inputindex] = std::stod(nanstr);
			}
			inputindex += 1;
		}
		if(inputindex != Barr.size())
		{
			std::cout << "Error in the barr_status input file" << "\n";
			throw std::invalid_argument("Error in the barr_status input file");
		}
		inputstream.close();
	}

	std::vector<double> R;
	for(int i = 0; i < numevents; i++)
	{
		R.push_back(0.0);
	}

	std::vector<double> phi_x, dphi_xdt, allx, Segments, Sigma, Torques, TSS_spotsigma, Velocities, Left_Barrier, Right_Barrier, GapR_Spot_Sigma;
	std::vector<std::vector<double>> phi, x, intimes, outtimes, outpos;
	std::vector<int> finish_count, current_rna_count, Gene_ID, Direction, nucl_seg_counts, Segment_RNAP_Count;
	std::vector<std::vector<int>> map, map_gapr;
	for(int i = 0; i < numgenes; i++)
	{
		phi.push_back(std::vector<double>());
		x.push_back(std::vector<double>());
		intimes.push_back(std::vector<double>());
		outtimes.push_back(std::vector<double>());
		outpos.push_back(std::vector<double>());
		finish_count.push_back(0);
		current_rna_count.push_back(0);
		TSS_spotsigma.push_back(0.0);
	}
	for(int i = 0; i < GapR_Points.size(); i++)
	{
		GapR_Spot_Sigma.push_back(0.0);
	}

	std::string line;
	std::vector<std::string> l;
	double last_timepoint = 0.0;
	int check_flag = 0;
	if(restartflag == 1)
	{
		inputstream.open(inputfolder + "/phi_file_" + std::to_string(world_rank) + ".log");
		inputindex = 0;
		while(getline(inputstream, line))
		{
			l.clear();
			splitstring(line, " ", &l, 1);
			for(int i = 0; i < l.size(); i++)
			{
				phi[inputindex].push_back(std::stod(l[i]));
			}
			inputindex += 1;
		}
		inputstream.close();
		if(inputindex != numgenes)
		{
			std::cout << "Error in the phi input file" << "\n";
			throw std::invalid_argument("Error in the phi input file");
		}

		inputstream.open(inputfolder + "/x_file_" + std::to_string(world_rank) + ".log");
		inputindex = 0;
		while(getline(inputstream, line))
		{
			l.clear();
			splitstring(line, " ", &l, 1);
			for(int i = 0; i < l.size(); i++)
			{
				x[inputindex].push_back(std::stod(l[i]));
			}
			inputindex += 1;
		}
		inputstream.close();
		if(inputindex != numgenes)
		{
			std::cout << "Error in the x input file" << "\n";
			throw std::invalid_argument("Error in the x input file");
		}

		inputstream.open(inputfolder + "/intimes_file_" + std::to_string(world_rank) + ".log");
		inputindex = 0;
		while(getline(inputstream, line))
		{
			l.clear();
			splitstring(line, " ", &l, 1);
			for(int i = 0; i < l.size(); i++)
			{
				intimes[inputindex].push_back(std::stod(l[i]));
			}
			inputindex += 1;
		}
		inputstream.close();
		if(inputindex != numgenes)
		{
			std::cout << "Error in the intimes input file" << "\n";
			throw std::invalid_argument("Error in the intimes input file");
		}

		inputstream.open(inputfolder + "/outtimes_file_" + std::to_string(world_rank) + ".log");
		inputindex = 0;
		while(getline(inputstream, line))
		{
			l.clear();
			splitstring(line, " ", &l, 1);
			for(int i = 0; i < l.size(); i++)
			{
				outtimes[inputindex].push_back(std::stod(l[i]));
			}
			inputindex += 1;
		}
		inputstream.close();
		if(inputindex != numgenes)
		{
			std::cout << "Error in the outtimes input file" << "\n";
			throw std::invalid_argument("Error in the outtimes input file");
		}

		inputstream.open(inputfolder + "/outpos_file_" + std::to_string(world_rank) + ".log");
		inputindex = 0;
		while(getline(inputstream, line))
		{
			l.clear();
			splitstring(line, " ", &l, 1);
			for(int i = 0; i < l.size(); i++)
			{
				outpos[inputindex].push_back(std::stod(l[i]));
			}
			finish_count[inputindex] += 1;
			inputindex += 1;
		}
		inputstream.close();
		if(inputindex != numgenes)
		{
			std::cout << "Error in the outpos input file" << "\n";
			throw std::invalid_argument("Error in the outpos input file");
		}

		inputstream.open(inputfolder + "/curr_count_file_" + std::to_string(world_rank) + ".log");
		inputindex = 0;
		while(inputstream >> inputval0)
		{
			current_rna_count[inputindex] = inputval0;
			inputindex += 1;
		}
		inputstream.close();
		if(inputindex != numgenes)
		{
			std::cout << "Error in the current rna input file" << "\n";
			throw std::invalid_argument("Error in the current rna input file");
		}

		inputstream.open(inputfolder + "/time_" + std::to_string(world_rank) + ".log");
		inputstream >> inputval1 >> inputval0 >> inputval2;
		last_timepoint = inputval1;
		if(inputval0 == 0)
		{
			std::cout << "That's all folks!" << "\n";
			return;
		}
		if(inputval2 > 0)
		{
			check_flag = inputval2;
		}
		inputstream.close();

		inputstream.open(inputfolder + "/generator_state_" + std::to_string(world_rank) + ".log");
		inputstream >> generator;
		inputstream.close();
	}

	int PolII_Count = 0, iterid = 0;

	for(int i = 0; i < numgenes; i++)
	{
		PolII_Count += x[i].size();
	}

	odesystem System{clamp0_flag, clamp1_flag, force, TSSes, gene_lengths, Gene_Directions, std::vector<double>(), std::vector<double>(), Nucl_Status, k_on, k_off, k_topo, k_nucl_on, k_nucl_off, Nucl, k_rna_degrad, k_prot_barrier_on, k_prot_barrier_off, Barr, Barr_Status, Barr_Phi, iterid, alpha_marenduzzo, std::vector<int>(), brute_force_flag, torque_flag, current_rna_count, GapR_On, GapR_Off, GapR_Points, GapR_State, galactose_switch, brute_gapr_flag, k_on_G, k_on_eR, k_off_R, k_off_G, k_off_eR, s_R, s_G, s_eR, s_PQS};
	wrappedsystem WS(&System);
	observer Observer{std::vector<double>(), std::vector<double>()};
	wrappedobserver OB(&Observer);

	double t = 0.0, dt = 0.0, a0 = 0.0, p0 = 0.0, p1 = 0.0, sum_new = 0.0, sum_prev = 0.0, epsilon = 1.0e-12, t0 = 0.0, dt0 = 0.0, dt0ini = 1.0, dist = 0.0, newphi = 0.0, hours_lapsed = 0.0, last_write_time = 0.0, write_interval = 1.0, x1 = 0.0, phi1 = 0.0, sigma_dash = 0.0;
	double x2 = 0.0, phi2 = 0.0, LK1, LK2, LK;
	int blocking_barr, henikoff_switch = 0;
	if(restartflag == 1)
	{
		t = last_timepoint;
	}
	last_write_time = t;
	std::vector<int> left_block, right_block;
	for(int i = 0; i < TSSes.size(); i++)
	{
		left_block.push_back(-1);
		right_block.push_back(-1);
	}

	int x_dict_count = 0;

	if(restartflag == 1)
	{
		Gene_ID.clear();
		for(int i = 0; i < numgenes; i++)
		{
			if(Gene_Directions[i] == 1)
			{
				for(int j = 0; j < phi[i].size(); j++)
				{
					phi_x.push_back(phi[i][j]);
					Gene_ID.push_back(i);
				}
			}
			else
			{
				for(int j = phi[i].size() - 1; j >= 0; j--)
				{
					phi_x.push_back(phi[i][j]);
					Gene_ID.push_back(i);
				}
			}
		}
		for(int i = 0; i < numgenes; i++)
		{
			if(Gene_Directions[i] == 1)
			{
				for(int j = 0; j < x[i].size(); j++)
				{
					phi_x.push_back(x[i][j]);
					x_dict_count += 1;
				}
			}
			else
			{
				for(int j = x[i].size() - 1; j >= 0; j--)
				{
					phi_x.push_back(x[i][j]);
					x_dict_count += 1;
				}
			}
		}
		phi_x.push_back(0.0);
		Alive = Is_Alive(phi_x, TSSes, gene_lengths, Gene_ID);
	}
	else
	{
		phi_x.push_back(0.0);
		Alive = std::vector<int>();
	}

	if(phi_x.size() != 2*x_dict_count + 1)
	{
		std::cout << "Mismatch between PolII count and phi_x vector size" << "\n";
		throw std::invalid_argument("Mismatch between PolII count and phi_x vector size");
	}

	Segments = Get_Segments_With_Barriers(phi_x, Alive, Barr, Barr_Status);
        Sigma = Get_Supercoiling_Densities_With_Barriers(clamp0_flag, clamp1_flag, phi_x, Segments, Alive, Barr, Barr_Status, Barr_Phi);

	int countsum = 0, rootindex, event = -1, orig_count, zero_flag, barrindex, barronflag, GapR_Index;

	std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
	std::chrono::steady_clock::time_point end_time;

	int exitflag = 0;

	std::ofstream alive_file, x_file, phi_file, velocity_file, angular_v_file, segments_file, sigma_file, torque_file, nucl_status_file, barr_status_file, barr_phi_file, nucl_pos_file, barr_pos_file, event_file, gapr_pos_file, gapr_file, nucl_type_file, rnap_density_file, rnap_segments_pos_file, spot_sigma_file, myong_file, sua_file;

	event_file.open(outputfolder + "/event_" + std::to_string(world_rank) + ".log");

	if(fileflag == 1)
	{
		alive_file.open(outputfolder + "/alive_" + std::to_string(world_rank) + ".log");
		x_file.open(outputfolder + "/x_" + std::to_string(world_rank) + ".log");
		phi_file.open(outputfolder + "/phi_" + std::to_string(world_rank) + ".log");
		velocity_file.open(outputfolder + "/velocity_" + std::to_string(world_rank) + ".log");
		angular_v_file.open(outputfolder + "/angular_v_" + std::to_string(world_rank) + ".log");
		segments_file.open(outputfolder + "/segments_" + std::to_string(world_rank) + ".log");
		sigma_file.open(outputfolder + "/sigma_" + std::to_string(world_rank) + ".log");
		torque_file.open(outputfolder + "/torque_" + std::to_string(world_rank) + ".log");
		nucl_status_file.open(outputfolder + "/nucl_status_" + std::to_string(world_rank) + ".log");
		barr_status_file.open(outputfolder + "/barr_status_" + std::to_string(world_rank) + ".log");
		barr_phi_file.open(outputfolder + "/barr_phi_" + std::to_string(world_rank) + ".log");
		myong_file.open(outputfolder + "/myong_" + std::to_string(world_rank) + ".log");
	}
	sua_file.open(outputfolder + "/sua_" + std::to_string(world_rank) + ".log");
	if(restartflag == 0 && world_rank == 0)
	{
		nucl_pos_file.open(outputfolder + "/nucl_pos_" + std::to_string(world_rank) + ".log");
		barr_pos_file.open(outputfolder + "/barr_pos_" + std::to_string(world_rank) + ".log");
		for(int i = 0; i < Nucl.size(); i++)
		{
			nucl_pos_file << Nucl[i] << "\n";
		}
		for(int i = 0; i < Barr.size(); i++)
		{
			barr_pos_file << Barr[i] << "\n";
		}
		nucl_pos_file.close();
		barr_pos_file.close();

		gapr_pos_file.open(outputfolder + "/gapr_pos_" + std::to_string(world_rank) + ".log");
		for(int i = 0; i < GapR_Points.size(); i++)
		{
			gapr_pos_file << GapR_Points[i] << "\n";
		}
		gapr_pos_file.close();
	}

	if(galactose_switch > -1)
	{
		gapr_file.open(outputfolder + "/gapr_" + std::to_string(world_rank) + ".log");
	}

	if(nucl_file_flag == 1)
	{
		nucl_status_file.open(outputfolder + "/nucl_status_" + std::to_string(world_rank) + ".log");
		nucl_type_file.open(outputfolder + "/nucl_type_" + std::to_string(world_rank) + ".log");
		rnap_density_file.open(outputfolder + "/rnap_density_" + std::to_string(world_rank) + ".log");
		spot_sigma_file.open(outputfolder + "/spot_sigma_" + std::to_string(world_rank) + ".log");

		if(restartflag == 0 && world_rank == 0)
		{
			rnap_segments_pos_file.open(outputfolder + "/rnap_segments_pos_" + std::to_string(world_rank) + ".log");
			rnap_segments_pos_file << clamp0 << "\n";
			for(int i = 0; i < RNAP_Points.size(); i++)
			{
				rnap_segments_pos_file << RNAP_Points[i] << "\n";
			}
			rnap_segments_pos_file << clamp1 << "\n";
			rnap_segments_pos_file.close();
		}
	}

	while((T > 0.0 && t < T) || (finish_count_limit > -1 && outtimes[0].size() < finish_count_limit))
	{
		orderflag = 0;
		for(int i = 0; i < numgenes; i++)
		{
			for(int j = 0; j < ((int) x[i].size()) - 1; j++)
			{
				if(Gene_Directions[i] == 1 && x[i][j + 1] > x[i][j])
				{
					orderflag = 1;
				}
				if(Gene_Directions[i] == -1 && x[i][j + 1] < x[i][j])
				{
					orderflag = 1;
				}
			}
		}
		if(orderflag == 1)
		{
			std::cout << "Genes are incorrectly ordered in the x vector." << "\n";
			throw std::invalid_argument("Genes are incorrectly ordered in the x vector.");
		}

		countsum = 0;
		for(int i = 0; i < numgenes; i++)
		{
			if(x[i].size() != phi[i].size())
			{
				std::cout << "Mismatch between x and phi vectors" << "\n";
				throw std::invalid_argument("Mismatch between x and phi vectors");
			}
			countsum += x[i].size();
		}
		if(countsum != PolII_Count)
		{
			std::cout << "Mismatch between PolII count and count from vectors" << "\n";
			throw std::invalid_argument("Mismatch between PolII count and count from vectors");
		}

		phi_x.clear();
		Gene_ID.clear();
		Direction.clear();

		System.Direction.clear();
		System.Exit_Times.clear();
		System.Exit_Points.clear();
		System.Gene_ID.clear();

		Observer.T.clear();
		Observer.P.clear();

		System.iterid = iterid;
		System.localiterid = 0;
		System.inipolIIcount = PolII_Count;

		for(int i = 0; i < numgenes; i++)
		{
			if(Gene_Directions[i] == 1)
			{
				for(int j = 0; j < phi[i].size(); j++)
				{
					phi_x.push_back(phi[i][j]);
					Gene_ID.push_back(i);
					Direction.push_back(Gene_Directions[i]);

					System.Direction.push_back(Gene_Directions[i]);
					System.Exit_Times.push_back(-1.0);
					System.Exit_Points.push_back(-1.0);
				}
			}
			else
			{
				for(int j = phi[i].size() - 1; j >= 0; j--)
				{
					phi_x.push_back(phi[i][j]);
					Gene_ID.push_back(i);
					Direction.push_back(Gene_Directions[i]);

					System.Direction.push_back(Gene_Directions[i]);
					System.Exit_Times.push_back(-1.0);
					System.Exit_Points.push_back(-1.0);
				}
			}
		}
		for(int i = 0; i < numgenes; i++)
		{
			if(Gene_Directions[i] == 1)
			{
				for(int j = 0; j < x[i].size(); j++)
				{
					phi_x.push_back(x[i][j]);
				}
			}
			else
			{
				for(int j = x[i].size() - 1; j >= 0; j--)
				{
					phi_x.push_back(x[i][j]);
				}
			}
		}
		phi_x.push_back(0.0);

		for(int i = 0; i < Nucl.size(); i++)
		{
			System.Nucl_Status[i] = Nucl_Status[i];
		}
		for(int i = 0; i < Barr.size(); i++)
		{
			System.Barr_Status[i] = Barr_Status[i];
		}
		for(int i = 0; i < Barr.size(); i++)
		{
			System.Barr_Phi[i] = Barr_Phi[i];
		}
		for(int i = 0; i < Gene_ID.size(); i++)
		{
			System.Gene_ID.push_back(Gene_ID[i]);
		}
		for(int i = 0; i < Nucl.size(); i++)
		{
			System.Nucl_Status[i] = Nucl_Status[i];
		}
		for(int i = 0; i < numgenes; i++)
		{
			System.current_rna_count[i] = current_rna_count[i];
		}
		for(int i = 0; i < GapR_Points.size(); i++)
		{
			System.GapR_State[i] = GapR_State[i];
		}
		System.galactose_switch = galactose_switch;
		System.s_R = s_R;
		System.s_G = s_G;
		System.s_eR = s_eR;

		Alive = Is_Alive(phi_x, TSSes, gene_lengths, Gene_ID);
		Segments = Get_Segments_With_Barriers(phi_x, Alive, Barr, Barr_Status);
		Sigma = Get_Supercoiling_Densities_With_Barriers(clamp0_flag, clamp1_flag, phi_x, Segments, Alive, Barr, Barr_Status, Barr_Phi);
		allx = Get_All_Of_Two_Vectors_Ordered(phi_x, Barr, (phi_x.size() - 1) / 2, 0, phi_x, Alive, Barr, Barr_Status, 7);
		Get_Nucl_RNAP_Map_Pairwise(Nucl, allx, map, delta, brute_force_flag);
		Is_Nucl_Blocked = Find_Blocked_Nucl(map);
		Get_GapR_RNAP_Map_Pairwise(GapR_Points, allx, map_gapr, delta);
		Is_GapR_Blocked = Find_Blocked_GapR(map_gapr, brute_gapr_flag);

		for(int i = 0; i < numgenes; i++)
		{
			TSS_spotsigma[i] = Get_Spot_Sigma(Segments, Sigma, TSSes[i]);
		}
		GQ_spotsigma = Get_Spot_Sigma(Segments, Sigma, TSSes[0] + GQ_spacer);

		for(int i = 0; i < numevents; i++)
		{
			R[i] = 0.0;
		}

		for(int i = 0; i < numgenes; i++)
		{
//			R[i] = k_on[i]*std::max(1.0 - alpha_marenduzzo*TSS_spotsigma[i], 0.0); // RNAP recruitment
			R[i] = k_on[i]*Marenduzzo_Func(TSS_spotsigma[i]); // RNAP recruitment
		}
		R[numgenes + 0] = k_off*PolII_Count; // RNAP drop-off before reaching the end of the gene
		R[numgenes + 1] = k_topo; // Supercoiling relaxation

		for(int i = 0; i < Nucl.size(); i++) // Nucleosome binding / unbinding
		{
			if(Nucl_Status[i] == 0)
			{
				if(Is_Nucl_Blocked[i] == 0)
				{
					R[i + numgenes + 2] = k_nucl_on;
				}
				else
				{
					R[i + numgenes + 2] = 0.0;
				}
				R[i + numgenes + 2 + Nucl.size()] = 0.0;
			}
			else
			{
				R[i + numgenes + 2] = 0.0;
				R[i + numgenes + 2 + Nucl.size()] = k_nucl_off;
			}
		}
		for(int i = 0; i < numgenes; i++) // RNA degradation
		{
			R[i + numgenes + 2 + 2*Nucl.size()] = current_rna_count[i]*k_rna_degrad;
		}
		for(int i = 0; i < Barr.size(); i++) // Barrier protein binding / unbinding
		{
			if(Barr_Status[i] == 0 && Is_Barr_Site_Blocked(Barr[i], phi_x, Alive, delta) == 1)
			{
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes] = Get_Barrier_On_Rate(k_prot_barrier_on, Get_Spot_Sigma(Segments, Sigma, Barr[i]));
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes + Barr.size()] = 0.0;
			}
			else if(Barr_Status[i] == 1)
			{
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes] = 0.0;
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes + Barr.size()] = Get_Barrier_Off_Rate(k_prot_barrier_off, Get_Spot_Sigma(Segments, Sigma, Barr[i]));

			}
		}
		for(int i = 0; i < GapR_Points.size(); i++) // Binding / unbinding of GapR
		{
			if(GapR_State[i] == 0 && Is_GapR_Blocked[i] == 0)
			{
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size()] = Get_GapR_On_Rate(GapR_On, Get_Spot_Sigma(Segments, Sigma, GapR_Points[i]), galactose_switch);
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + GapR_Points.size()] = 0.0;
			}
			else if(GapR_State[i] == 0 && Is_GapR_Blocked[i] == 1)
			{
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size()] = 0.0;
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + GapR_Points.size()] = 0.0;
			}
			else if(GapR_State[i] == 1)
			{
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size()] = 0.0;
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + GapR_Points.size()] = GapR_Off;
			}
		}
		if(s_R == 0)
		{
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size()] = k_on[0]*Marenduzzo_Func(TSS_spotsigma[0])*(1.0 - s_R);
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 1] = 0.0;
		}
		else if(s_R == 1)
		{
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size()] = 0.0;
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 1] = (1.0 - s_eR)*s_R*(k_off_R*(1.0 - s_G) + s_G*k_off_eR);
		}
		else
		{
			std::cout << "Invalid R-loop state (Sim file)" << "\n";
			throw std::invalid_argument("Invalid R-loop state (Sim file)");
		}
		if(s_G == 0)
		{
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 2] = s_PQS*s_R*(1.0 - s_G)*k_on_G;
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 3] = 0.0;
		}
		else if(s_G == 1)
		{
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 2] = 0.0;
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 3] = k_off_G;
		}
		else
		{
			std::cout << "Invalid GQ state (Sim file)" << "\n";
			throw std::invalid_argument("Invalid GQ state (Sim file)");
		}
		if(s_eR == 0)
		{
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 4] = s_R*(1.0 - s_eR)*k_on_eR;
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 5] = 0.0;
		}
		else if(s_eR == 1)
		{
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 4] = 0.0;
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 5] = k_off_eR;
		}
		else
		{
			std::cout << "Invalid eR state (Sim file)" << "\n";
			throw std::invalid_argument("Invalid eR state (Sim file)");
		}

		a0 = 0.0;
		for(int i = 0; i < numevents; i++)
		{
			if(R[i] < 0.0)
			{
				std::cout << "At least one of the rates is invalid." << "\n";
				throw std::invalid_argument("At least one of the rates is invalid.");
			}
			a0 += R[i];
		}

		p0 = distribution(generator);
		boost::numeric::odeint::integrate_const(boost::numeric::odeint::euler<std::vector<double>>(), WS, phi_x, 0.0, (1.0 / a0)*std::log(2.0 / p0), 1.0e-3, OB);
		for(int i = 0; i < Observer.P.size() - 1; i++)
		{
			if(Observer.P[i + 1] < Observer.P[i])
			{
				std::cout << "a0 vector is non-monotonic; this is a mistake." << "\n";
				throw std::invalid_argument("a0 vector is non-monotonic; this is a mistake.");
			}
		}

		rootindex = -1;
		dt = 0.0;
		for(int i = 0; i < Observer.P.size(); i++)
		{
			if(Observer.P[i] > std::log(1.0 / p0))

			{
				break;
			}
			dt = Observer.T[i];
			rootindex = i;
		}
		if(std::log(1.0 / p0) - Observer.P[rootindex] > 1.0)
		{
			approx_count += 1;
			if(std::log(1.0 / p0) - Observer.P[rootindex] > max_approx_diff)
			{
				max_approx_diff = std::log(1.0 / p0) - Observer.P[rootindex];
			}
		}

		t += dt;

		if(PolII_Count > 0) // Integrate only if there are any RNAPs currently transcribing
		{
			phi_x.clear();
			Gene_ID.clear();
			Direction.clear();

			System.Direction.clear();
			System.Exit_Times.clear();
			System.Exit_Points.clear();
			System.Gene_ID.clear();

			Observer.T.clear();
			Observer.P.clear();

			System.iterid = iterid;
			System.localiterid = 0;
			System.inipolIIcount = PolII_Count;

			for(int i = 0; i < numgenes; i++)
			{
				if(Gene_Directions[i] == 1)
				{
					for(int j = 0; j < phi[i].size(); j++)
					{
						phi_x.push_back(phi[i][j]);
						Gene_ID.push_back(i);
						Direction.push_back(Gene_Directions[i]);

						System.Direction.push_back(Gene_Directions[i]);
						System.Exit_Times.push_back(-1.0);
						System.Exit_Points.push_back(-1.0);
					}
				}
				else
				{
					for(int j = phi[i].size() - 1; j >= 0; j--)
					{
						phi_x.push_back(phi[i][j]);
						Gene_ID.push_back(i);
						Direction.push_back(Gene_Directions[i]);

						System.Direction.push_back(Gene_Directions[i]);
						System.Exit_Times.push_back(-1.0);
						System.Exit_Points.push_back(-1.0);
					}
				}
			}
			for(int i = 0; i < numgenes; i++)
			{
				if(Gene_Directions[i] == 1)
				{
					for(int j = 0; j < x[i].size(); j++)
					{
						phi_x.push_back(x[i][j]);
					}
				}
				else
				{
					for(int j = x[i].size() - 1; j >= 0; j--)
					{
						phi_x.push_back(x[i][j]);
					}
				}
			}
			phi_x.push_back(0.0);

			if((phi_x.size() - 1) / 2 != PolII_Count)
			{
				std::cout << "Mismatch between PolII count and phi_x vector size" << "\n";
				throw std::invalid_argument("Mismatch between PolII count and phi_x vector size");
			}

			for(int i = 0; i < Nucl.size(); i++)
			{
				System.Nucl_Status[i] = Nucl_Status[i];
			}
			for(int i = 0; i < Barr.size(); i++)
			{
				System.Barr_Status[i] = Barr_Status[i];
			}
			for(int i = 0; i < Barr.size(); i++)
			{
				System.Barr_Phi[i] = Barr_Phi[i];
			}
			for(int i = 0; i < Gene_ID.size(); i++)
			{
				System.Gene_ID.push_back(Gene_ID[i]);
			}
			for(int i = 0; i < numgenes; i++)
			{
				System.current_rna_count[i] = current_rna_count[i];
			}
			for(int i = 0; i < GapR_Points.size(); i++)
			{
				System.GapR_State[i] = GapR_State[i];
			}
			System.galactose_switch = galactose_switch;
			System.s_R = s_R;
			System.s_G = s_G;
			System.s_eR = s_eR;

			for(int i = 0; i < PolII_Count - 1; i++)
			{
				if(phi_x[i + PolII_Count] < phi_x[i + PolII_Count + 1])
				{
					std::cout << "phi_x vector is incorrectly ordered; must be monotonic" << "\n";
					throw std::invalid_argument("phi_x vector is incorrectly ordered; must be monotonic");
				}
			}

			if(check_flag > 0) // Checking that Get_Equalizing_Phi worked
			{
				Alive = Is_Alive(phi_x, TSSes, gene_lengths, Gene_ID);
				Segments = Get_Segments_With_Barriers(phi_x, Alive, Barr, Barr_Status);
        			Sigma = Get_Supercoiling_Densities_With_Barriers(clamp0_flag, clamp1_flag, phi_x, Segments, Alive, Barr, Barr_Status, Barr_Phi);
				zero_flag = 0;
				for(int i = 0; i < Sigma.size() - 1; i++)
				{
					if(std::abs(Sigma[i] - Sigma[i + 1]) < epsilon)
					{
						zero_flag = 1;
					}
				}
				if(zero_flag == 0)
				{
					std::cout << "At least one segment with 0 sigma was expected; none found." << "\n";
					throw std::invalid_argument("At least one segment with 0 sigma is expected; none found");
				}
				check_flag = 0;
			}

			if(fileflag == 0)
			{
				boost::numeric::odeint::integrate_const(boost::numeric::odeint::euler<std::vector<double>>(), WS, phi_x, t - dt, t, 1.0e-3, OB);
			}
			else // If the trajectory has to be written out, do stepwise integration instead of all at once
			{
				t0 = t - dt;
				dt0 = dt0ini;
				while(t0 < t)
				{
					if(t0 + dt0 > t)
					{
						dt0 = t - t0;
					}

					dphi_xdt.clear();
					for(int i = 0; i < phi_x.size(); i++)
					{
						dphi_xdt.push_back(0.0);
					}

					System(phi_x, dphi_xdt, t0); // Get the derivatives at time t0

					allx = Get_All_Of_Two_Vectors_Ordered(phi_x, Barr, (phi_x.size() - 1) / 2, 0, phi_x, Alive, Barr, Barr_Status, 8);
					Alive = Is_Alive(phi_x, TSSes, gene_lengths, Gene_ID);
					Segments = Get_Segments_With_Barriers(phi_x, Alive, Barr, Barr_Status);
					Sigma = Get_Supercoiling_Densities_With_Barriers(clamp0_flag, clamp1_flag, phi_x, Segments, Alive, Barr, Barr_Status, Barr_Phi);
					nucl_seg_counts = Get_Segment_Nucleosome_Count_Pairwise(allx, Nucl, Nucl_Status);
					if(torque_flag == 0)
					{
						Torques = Get_Prokaryotic_Torques(Sigma, Segments, force);
					}
					else if(torque_flag == 1)
					{
						Torques = Get_Interpolated_Eukaryotic_Torques(interp_ML, interp_ML_Cutoff_0, interp_ML_Cutoff_1, Sigma, Segments, nucl_seg_counts, 1.0);
					}
					Is_RNAP_Blocked = Find_RNAP_Blocked_By_Nucl(Nucl, Nucl_Status, phi_x, Alive, Direction, delta, brute_force_flag);
					Velocities = Get_Velocities_With_Barriers(phi_x, Segments, Torques, System.Direction, Alive, Barr, Barr_Status, Is_RNAP_Blocked, delta);

					// Writing begins

					Write_Vector_State_To_File<int>(&alive_file, t0, Alive, Alive.size(), 0);
					Write_Vector_State_To_File<double>(&x_file, t0, phi_x, PolII_Count, PolII_Count);
					Write_Vector_State_To_File<double>(&phi_file, t0, phi_x, PolII_Count, 0);
					Write_Vector_State_To_File<double>(&velocity_file, t0, Velocities, Velocities.size(), 0);
					Write_Vector_State_To_File<double>(&angular_v_file, t0, dphi_xdt, PolII_Count, 0);
					Write_Vector_State_To_File<double>(&segments_file, t0, Segments, Segments.size(), 0);
					Write_Vector_State_To_File<double>(&sigma_file, t0, Sigma, Sigma.size(), 0);
					Write_Vector_State_To_File<double>(&torque_file, t0, Torques, Torques.size(), 0);
					Write_Vector_State_To_File<int>(&nucl_status_file, t0, Nucl_Status, Nucl.size(), 0);
					Write_Vector_State_To_File<int>(&barr_status_file, t0, Barr_Status, Barr.size(), 0);
					Write_Vector_State_To_File<double>(&barr_phi_file, t0, Barr_Phi, Barr.size(), 0);
					myong_file << t0 << "\t" << Get_Spot_Sigma(Segments, Sigma, TSSes[0] + GQ_spacer) << "\t" << System.s_R << "\t" << System.s_G << "\t" << System.s_eR << "\t" << k_on[0]*Marenduzzo_Func(Get_Spot_Sigma(Segments, Sigma, TSSes[0])) << "\n";

					// Writing ends

					boost::numeric::odeint::integrate_const(boost::numeric::odeint::euler<std::vector<double>>(), WS, phi_x, t0, t0 + dt0, 1.0e-3, OB);
					t0 += dt0;
				}
			} // Stepwise integration ends

			Alive = Is_Alive(phi_x, TSSes, gene_lengths, Gene_ID);
			for(int i = 0; i < Alive.size(); i++)
			{
				if((Alive[i] == 1 && System.Exit_Times[i] > 0.0)) // Alive RNAPs should not have a valid exit time
				{
					std::cout << "An alive RNAP somehow has a valid exit time." << "\n";
					throw std::invalid_argument("An alive RNAP somehow has a valid exit time.");
				}
				if(System.Exit_Times[i] < 0.0 && Alive[i] == 0) // Set exit times for newly dead RNAPs (RNAPs that died after leaving the integrator)
				{
					System.Exit_Times[i] = t;
					System.Exit_Points[i] = phi_x[i + PolII_Count];
					System.current_rna_count[Gene_ID[i]] += 1;
				}
				if(Alive[i] == 0 && System.Exit_Times[i] < 0.0) // Dead RNAPs should no longer have an invalid exit time
				{
					std::cout << "One of the dead RNAPs has an invalid exit time." << "\n";
					throw std::invalid_argument("One of the dead RNAPs has an invalid exit time.");
				}
			}

			for(int i = 0; i < Alive.size(); i++) // Updating the local variables for RNAPs that have died
			{
				if(Alive[i] == 0)
				{
					outtimes[Gene_ID[i]].push_back(System.Exit_Times[i]);
					outpos[Gene_ID[i]].push_back(System.Exit_Points[i]);
					current_rna_count[Gene_ID[i]] += 1;
					finish_count[Gene_ID[i]] += 1;
				}
			}

			// Updating the local state of RNAPs

			for(int i = 0; i < numgenes; i++)
			{
				phi[i].clear();
				x[i].clear();
			}
			PolII_Count = 0;
			for(int i = 0; i < (phi_x.size() - 1) / 2; i++)
			{
				if(Alive[i] == 0)
				{
					continue;
				}
				phi[Gene_ID[i]].push_back(phi_x[i]);
				x[Gene_ID[i]].push_back(phi_x[i + ((phi_x.size() - 1) / 2)]);
				PolII_Count += 1;
			}
			for(int i = 0; i < numgenes; i++) // Reversing for genes transcribing in the reverse direction
			{
				if(Gene_Directions[i] != 1)
				{
					std::reverse(phi[i].begin(), phi[i].end());
					std::reverse(x[i].begin(), x[i].end());
				}
			}
		} // Ending the integration step

		Alive = Is_Alive(phi_x, TSSes, gene_lengths, Gene_ID);
		Segments = Get_Segments_With_Barriers(phi_x, Alive, Barr, Barr_Status);
		Sigma = Get_Supercoiling_Densities_With_Barriers(clamp0_flag, clamp1_flag, phi_x, Segments, Alive, Barr, Barr_Status, Barr_Phi);
		allx = Get_All_Of_Two_Vectors_Ordered(phi_x, Barr, (phi_x.size() - 1) / 2, 0, phi_x, Alive, Barr, Barr_Status, 9);
		Get_Nucl_RNAP_Map_Pairwise(Nucl, allx, map, delta, brute_force_flag);
		Is_Nucl_Blocked = Find_Blocked_Nucl(map);
		Get_GapR_RNAP_Map_Pairwise(GapR_Points, allx, map_gapr, delta);
		Is_GapR_Blocked = Find_Blocked_GapR(map_gapr, brute_gapr_flag);

		Segment_RNAP_Count = Get_Segment_RNAP_Count(RNAP_Points, phi_x, Alive);

		for(int i = 0; i < numgenes; i++)
		{
			TSS_spotsigma[i] = Get_Spot_Sigma(Segments, Sigma, TSSes[i]);
		}
		GQ_spotsigma = Get_Spot_Sigma(Segments, Sigma, TSSes[0] + GQ_spacer);

		for(int i = 0; i < GapR_Points.size(); i++)
		{
			GapR_Spot_Sigma[i] = Get_Spot_Sigma(Segments, Sigma, GapR_Points[i]);
		}

		for(int i = 0; i < numevents; i++)
		{
			R[i] = 0.0;
		}

		for(int i = 0; i < numgenes; i++)
		{
//			R[i] = k_on[i]*std::max(1.0 - alpha_marenduzzo*TSS_spotsigma[i], 0.0); // RNAP recruitment
			R[i] = k_on[i]*Marenduzzo_Func(TSS_spotsigma[i]); // RNAP recruitment
		}
		R[numgenes + 0] = k_off*PolII_Count; // RNAP drop-off before reaching the end of the gene
		R[numgenes + 1] = k_topo; // Supercoiling relaxation

		for(int i = 0; i < Nucl.size(); i++) // Nucleosome binding / unbinding
		{
			if(Nucl_Status[i] == 0)
			{
				if(Is_Nucl_Blocked[i] == 0)
				{
					R[i + numgenes + 2] = k_nucl_on;
				}
				else
				{
					R[i + numgenes + 2] = 0.0;
				}
				R[i + numgenes + 2 + Nucl.size()] = 0.0;
			}
			else
			{
				R[i + numgenes + 2] = 0.0;
				R[i + numgenes + 2 + Nucl.size()] = k_nucl_off;
			}
		}
		for(int i = 0; i < numgenes; i++) // RNA degradation
		{
			R[i + numgenes + 2 + 2*Nucl.size()] = current_rna_count[i]*k_rna_degrad;
		}
		for(int i = 0; i < Barr.size(); i++) // Barrier protein binding / unbinding
		{
			if(Barr_Status[i] == 0 && Is_Barr_Site_Blocked(Barr[i], phi_x, Alive, delta) == 1)
			{
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes] = Get_Barrier_On_Rate(k_prot_barrier_on, Get_Spot_Sigma(Segments, Sigma, Barr[i]));
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes + Barr.size()] = 0.0;
			}
			else if(Barr_Status[i] == 1)
			{
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes] = 0.0;
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes + Barr.size()] = Get_Barrier_Off_Rate(k_prot_barrier_off, Get_Spot_Sigma(Segments, Sigma, Barr[i]));

			}
		}
		for(int i = 0; i < GapR_Points.size(); i++) // Binding / unbinding of GapR
		{
			if(GapR_State[i] == 0 && Is_GapR_Blocked[i] == 0)
			{
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size()] = Get_GapR_On_Rate(GapR_On, Get_Spot_Sigma(Segments, Sigma, GapR_Points[i]), galactose_switch);
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + GapR_Points.size()] = 0.0;
			}
			else if(GapR_State[i] == 0 && Is_GapR_Blocked[i] == 1)
			{
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size()] = 0.0;
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + GapR_Points.size()] = 0.0;
			}
			else if(GapR_State[i] == 1)
			{
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size()] = 0.0;
				R[i + numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + GapR_Points.size()] = GapR_Off;
			}
		}
		if(s_R == 0)
		{
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size()] = k_on[0]*Marenduzzo_Func(TSS_spotsigma[0])*(1.0 - s_R);
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 1] = 0.0;
		}
		else if(s_R == 1)
		{
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size()] = 0.0;
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 1] = (1.0 - s_eR)*s_R*(k_off_R*(1.0 - s_G) + s_G*k_off_eR);
		}
		else
		{
			std::cout << "Invalid R-loop state (Sim file 2)" << "\n";
			throw std::invalid_argument("Invalid R-loop state (Sim file 2)");
		}
		if(s_G == 0)
		{
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 2] = s_PQS*s_R*(1.0 - s_G)*k_on_G;
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 3] = 0.0;
		}
		else if(s_G == 1)
		{
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 2] = 0.0;
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 3] = k_off_G;
		}
		else
		{
			std::cout << "Invalid GQ state (Sim file 2)" << "\n";
			throw std::invalid_argument("Invalid GQ state (Sim file 2)");
		}
		if(s_eR == 0)
		{
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 4] = s_R*(1.0 - s_eR)*k_on_eR;
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 5] = 0.0;
		}
		else if(s_eR == 1)
		{
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 4] = 0.0;
			R[numgenes + 2 + 2*Nucl.size() + numgenes + 2*Barr.size() + 2*GapR_Points.size() + 5] = k_off_eR;
		}
		else
		{
			std::cout << "Invalid eR state (Sim file 2)" << "\n";
			throw std::invalid_argument("Invalid eR state (Sim file 2)");
		}

		a0 = 0.0;
		for(int i = 0; i < numevents; i++)
		{
			if(R[i] < 0.0)
			{
				std::cout << "One of the rates is invalid: once again" << "\n";
				throw std::invalid_argument("One of the rates is invalid: once again");
			}
			a0 += R[i];
		}

		p1 = distribution(generator);
		sum_prev = 0.0;
		sum_new = 0.0;
		event = -1;
		for(int i = 0; i < numevents; i++)
		{
			sum_new = sum_prev + (R[i] / a0);
			if(p1 >= sum_prev && p1 < sum_new)
			{
				event = i;
				break;
			}
			sum_prev = sum_new;
		}

		if(std::abs(R[event]) < epsilon)
		{
			std::cout << "The chosen event has a zero rate." << "\n";
			throw std::invalid_argument("The chosen event has a zero rate.");
		}

		// Event execution begins

		orig_count = PolII_Count;

		if(event < numgenes) // RNAP recruitment
		{
			Is_TSS_Blocked = Find_TSSes_Blocked_By_Nucl(Nucl, Nucl_Status, TSSes, left_block, right_block, brute_force_flag, delta);
			if(left_block[event] > -1 && Nucl[left_block[event]] > TSSes[event])
			{
				std::cout << "Incorrect left blocking nucleosome reported" << "\n";
				throw std::invalid_argument("Incorrect left blocking nucleosome reported");
			}
			if(right_block[event] > -1 && TSSes[event] > Nucl[right_block[event]])
			{
				std::cout << "Incorrect right blocking nucleosome reported" << "\n";
				throw std::invalid_argument("Incorrect right blocking nucleosome reported");
			}
			blocking_barr = -1;
			if(single_pol_II_flag == 0 || (single_pol_II_flag == 1 && intimes[0].size() == 0))
			{
				newphi = Get_Equalizing_Phi(Segments, Barr, TSSes[event], phi_x, Alive, Barr_Status, Barr_Phi, clamp0_flag, clamp1_flag);
				if(x[event].size() == 0)
				{
					phi[event].push_back(newphi);
					x[event].push_back(TSSes[event]);
					intimes[event].push_back(t);
					PolII_Count += 1;
				}
				else
				{
					if(std::abs(TSSes[event] - x[event][x[event].size() - 1]) > delta)
					{
						phi[event].push_back(newphi);
						x[event].push_back(TSSes[event]);
						intimes[event].push_back(t);
						PolII_Count += 1;
					}
				}
			}
			if(PolII_Count > orig_count) // If an RNAP has been recruited, remove the nucleosome blocking the TSS or the barrier protein blocking the TSS (if any)
			{
				if(Is_TSS_Blocked[event] == 1)
				{
					if(left_block[event] > -1)
					{
						if(Nucl_Status[left_block[event]] != 1)
						{
							std::cout << "Left blocking nucleosome does not really exist." << "\n";
							throw std::invalid_argument("Left blocking nucleosome does not really exist.");
						}
						Nucl_Status[left_block[event]] = 0;
					}
					if(right_block[event] > -1)
					{
						if(Nucl_Status[right_block[event]] != 1)
						{
							std::cout << "Right blocking nucleosome does not really exist." << "\n";
							throw std::invalid_argument("Right blocking nucleosome does not really exist.");
						}
						Nucl_Status[right_block[event]] = 0;
					}
				}
				blocking_barr = -1;
				if(Is_TSS_Blocked_By_Barrier(TSSes, Barr, Barr_Status, event, &blocking_barr, delta) == 1)
				{
					if(blocking_barr == -1 || Barr_Status[blocking_barr] == 0)
					{
						std::cout << "TSS blocking barrier does not exist." << "\n";
						throw std::invalid_argument("TSS blocking barrier does not exist.");
					}
					Barr_Status[blocking_barr] = 0;
				}
				check_flag = 1; // Check if the Get_Equalizing_Phi function really works
			}
		}
		else if(event == numgenes) // RNAP drop-off (currently deactivated)
		{
			std::cout << "RNAP dropping-off mid transcription even though the event has been deactivated." << "\n";
			throw std::invalid_argument("RNAP dropping-off mid transcription even though the event has been deactivated.");
		}
		else if(event == numgenes + 1) // Supercoiling relaxation
		{
			if(circular_plasmid_flag == 0)
			{
				if(x[0].size() > 0)
				{
					x1 = x[0][0];
					phi1 = phi[0][0];

					x2 = x[0][x[0].size() - 1];
					phi2 = phi[0][phi[0].size() - 1];

					LK1 = Sigma[0]*(Segments[0] / 3.4);
					LK2 = Sigma[Sigma.size() - 1]*(Segments[Segments.size() - 1] / 3.4);

					LK = LK1 + LK2;
					sigma_dash = (LK) / ((Segments[0] + Segments[Segments.size() - 1]) / 3.4);

					clamp0_sigma = phi2 + (w0*sigma_dash*Segments[Segments.size() - 1]);
					clamp1_sigma = phi1 - (w0*sigma_dash*Segments[0]);
				}
				else
				{
					;
				}
			}

//			for(int i = 0; i < numgenes; i++) // Relaxation at RNAP sites
//			{
//				for(int j = 0; j < phi[i].size(); j++)
//				{
//					phi[i][j] = 0.0;
//				}
//			}
//			for(int i = 0; i < Barr_Status.size(); i++) // Relaxation at barrier protein sites
//			{
//				if(Barr_Status[i] == 1)
//				{
//					Barr_Phi[i] = 0.0;
//				}
//			}
		}
		else if(event > numgenes + 1 && event <= numgenes + 1 + Nucl.size()) // Nucleosome binding
		{
			if(Nucl_Status[event - numgenes - 2] != 0)
			{
				std::cout << "Nucleosome being recruited to already occupied site." << "\n";
				throw std::invalid_argument("Nucleosome being recruited to already occupied site.");
			}
			if(Is_Nucl_Blocked[event - numgenes - 2] != 0)
			{
				std::cout << "Nucleosome being recruited to a blocked site." << "\n";
				throw std::invalid_argument("Nucleosome being recruited to a blocked site.");
			}
			Nucl_Status[event - numgenes - 2] = 1;
			if(henikoff_switch == 0)
			{
				Nucl_Type[event - numgenes - 2] = 1;
			}
			else
			{
				Nucl_Type[event - numgenes - 2] = 2;
			}
		}
		else if(event > numgenes + 1 + Nucl.size() && event <= numgenes + 1 + 2*Nucl.size()) // Nucleosome unbinding
		{
			if(Nucl_Status[event - numgenes - Nucl.size() - 2] != 1)
			{
				std::cout << "Unbinding of a nucleosome that does not exist" << "\n";
				throw std::invalid_argument("Unbinding of a nucleosome that does not exist");
			}
			Nucl_Status[event - numgenes - Nucl.size() - 2] = 0;
			Nucl_Type[event - numgenes - Nucl.size() - 2] = 0;
		}
		else if(event > numgenes + 1 + 2*Nucl.size() && event <= numgenes + 1 + 2*Nucl.size() + numgenes) // RNA degradation
		{
			current_rna_count[event - (numgenes + 1 + 2*Nucl.size() + 1)] -= 1;
		}
		else if(event > 2*numgenes + 1 + 2*Nucl.size() && event <= 2*numgenes + 1 + 2*Nucl.size() + Barr.size()) // Barrier protein binding
		{
			barrindex = event - (2*numgenes + 1 + 2*Nucl.size() + 1);
			barronflag = 1;
			for(int i = 0; i < (phi_x.size() - 1) / 2; i++)
			{
				if(std::abs(phi_x[i + (phi_x.size() - 1) / 2] - Barr[barrindex]) < delta && Alive[i] == 1)
				{
					barronflag = 0;
					break;
				}
			}
			if(barronflag == 1)
			{
				newphi = Get_Equalizing_Phi(Segments, Barr, Barr[event - (2*numgenes + 1 + 2*Nucl.size() + 1)], phi_x, Alive, Barr_Status, Barr_Phi, clamp0_flag, clamp1_flag);
				Barr_Phi[barrindex] = newphi;
				Barr_Status[barrindex] = 1;
				check_flag = 2;
			}
			else
			{
				std::cout << "Barrier protein binding to a site blocked by RNAP" << "\n";
				throw std::invalid_argument("Barrier protein binding to a site blocked by RNAP");
			}
		}
		else if(event > 2*numgenes + 1 + 2*Nucl.size() + Barr.size() && event <= 2*numgenes + 1 + 2*Nucl.size() + 2*Barr.size()) // Barrier protein unbinding
		{
			barrindex = event - (2*numgenes + 1 + 2*Nucl.size() + Barr.size() + 1);
			if(Barr_Status[barrindex] == 0)
			{
				std::cout << "Unbinding of non-existent barrier protein" << "\n";
				throw std::invalid_argument("Unbinding of non-existent barrier protein");
			}
			Barr_Status[barrindex] = 0;
			Barr_Phi[barrindex] = NAN;
		}
		else if(event > 2*numgenes + 1 + 2*Nucl.size() + 2*Barr.size() && event <= 2*numgenes + 1 + 2*Nucl.size() + 2*Barr.size() + GapR_Points.size())
		{
			GapR_Index = event - (2*numgenes + 1 + 2*Nucl.size() + 2*Barr.size() + 1);
			if(GapR_State[GapR_Index] != 0)
			{
				std::cout << "GapR protein binding to an already occupied site" << "\n";
				throw std::invalid_argument("GapR protein binding to an already occupied site");
			}
			if(Is_GapR_Blocked[GapR_Index] == 1)
			{
				std::cout << "GapR protein binding to a blocked site" << "\n";
				throw std::invalid_argument("GapR protein binding to a blocked site");
			}
			GapR_State[GapR_Index] = 1;
		}
		else if(event > 2*numgenes + 1 + 2*Nucl.size() + 2*Barr.size() + GapR_Points.size() && event <= 2*numgenes + 1 + 2*Nucl.size() + 2*Barr.size() + 2*GapR_Points.size())
		{
			GapR_Index = event - (2*numgenes + 1 + 2*Nucl.size() + 2*Barr.size() + GapR_Points.size() + 1);
			if(GapR_State[GapR_Index] != 1)
			{
				std::cout << "GapR protein being unbound does not exist." << "\n";
				throw std::invalid_argument("GapR protein being unbound does not exist.");
			}
			GapR_State[GapR_Index] = 0;
		}
		else if(event == 2*numgenes + 1 + 2*Nucl.size() + 2*Barr.size() + 2*GapR_Points.size() + 1)
		{
			s_R = 1;
		}
		else if(event == 2*numgenes + 1 + 2*Nucl.size() + 2*Barr.size() + 2*GapR_Points.size() + 2)
		{
			s_R = 0;
		}
		else if(event == 2*numgenes + 1 + 2*Nucl.size() + 2*Barr.size() + 2*GapR_Points.size() + 3)
		{
			s_G = 1;
		}
		else if(event == 2*numgenes + 1 + 2*Nucl.size() + 2*Barr.size() + 2*GapR_Points.size() + 4)
		{
			s_G = 0;
		}
		else if(event == 2*numgenes + 1 + 2*Nucl.size() + 2*Barr.size() + 2*GapR_Points.size() + 5)
		{
			s_eR = 1;
			if(circular_plasmid_flag == 0)
			{
				if(x[0].size() > 0)
				{
					sigma_dash = -(10.0) / (Segments[Segments.size() - 1] / 3.4);
					sigma_dash = Sigma[Sigma.size() - 1] - sigma_dash;
					clamp0_sigma = phi[0][phi[0].size() - 1] + (sigma_dash*w0*Segments[Segments.size() - 1]);
				}
				else
				{
					sigma_dash = -(10.0) / (Segments[Segments.size() - 1] / 3.4);
					sigma_dash = Sigma[Sigma.size() - 1] - sigma_dash;
					clamp0_sigma = clamp1_sigma + (sigma_dash*w0*Segments[Segments.size() - 1]);
				}
			}
			else
			{
				sigma_dash = -(10.0) / ((3700.0*0.34) / 3.4);
				sigma_basal = sigma_basal - sigma_dash;
			}
		}
		else if(event == 2*numgenes + 1 + 2*Nucl.size() + 2*Barr.size() + 2*GapR_Points.size() + 6)
		{
			s_eR = 0;
		}
		else
		{
			std::cout << "The chosen event type is invalid." << "\n";
			throw std::invalid_argument("The chosen event type is invalid.");
		}

		if(PolII_Count - orig_count > 1)
		{
			std::cout << "More than one PolII recruited during the simulation: impossible event" << "\n";
			throw std::invalid_argument("More than one PolII recruited during the simulation: impossible event");
		}
		iterid += 1;

//		event_file << t << " "  << event << " " << Nucl.size() << " " << Barr.size() << "\n";
		event_file << t << " "  << event << "\n";

		if(t - last_write_time > write_interval)
		{
			sua_file << t << "\t" << s_R << "\t" << s_G << "\t" << s_eR << "\n";
			if(galactose_switch == 1)
			{
				gapr_file << t << "\t";
				for(int i = 0; i < GapR_Points.size(); i++)
				{
					gapr_file << GapR_State[i] << "\t";
				}
				gapr_file << "\n";
			}
			if(nucl_file_flag == 1)
			{
				spot_sigma_file << t << "\t";
				for(int i = 0; i < GapR_Points.size(); i++)
				{
					spot_sigma_file << GapR_Spot_Sigma[i] << "\t";
				}
				spot_sigma_file << "\n";

				nucl_status_file << t << "\t";
				nucl_type_file << t << "\t";
				for(int i = 0; i < Nucl.size(); i++)
				{
					nucl_status_file << Nucl_Status[i] << "\t";
					nucl_type_file << Nucl_Type[i] << "\t";
				}
				nucl_status_file << "\n";
				nucl_type_file << "\n";
	
				rnap_density_file << t - dt << "\t";
				for(int i = 0; i < Segment_RNAP_Count.size(); i++)
				{
					rnap_density_file << Segment_RNAP_Count[i] << "\t";
				}
				rnap_density_file << "\n";
			}

			last_write_time = t;
		}

		if(t > galactose_time && galactose_switch == 0)
		{
			galactose_switch = 1;
		}

		if(henikoff_switch_flag == 1 && t > henikoff_time && henikoff_switch == 0)
		{
			henikoff_switch = 1;
		}

		end_time = std::chrono::steady_clock::now();
		hours_lapsed = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count() / (3600.0);

		if(hours_lapsed > 18.0)
		{
			exitflag = 1;
			break;
		}
	}

	std::ofstream rates(outputfolder + "/rates_" + std::to_string(world_rank) + ".log"); // Avg. velocity in base pairs / s
	for(int i = 0; i < numgenes; i++)
	{
		if(Gene_Directions[i] == 1)
		{
			for(int j = 0; j < outtimes[i].size(); j++)
			{
				rates << ((outpos[i][j] - TSSes[i]) / 0.34) / (outtimes[i][j] - intimes[i][j]) << "\t";
			}
		}
		else if(Gene_Directions[i] == -1)
		{
			for(int j = 0; j < outtimes[i].size(); j++)
			{
				rates << ((TSSes[i] - outpos[i][j]) / 0.34) / (outtimes[i][j] - intimes[i][j]) << "\t";
			}
		}
		else
		{
			std::cout << "Invalid gene direction while calculating overall RNAP velocity" << "\n";
			throw std::invalid_argument("Invalid gene direction while calculating overall RNAP velocity");
		}
		rates << "\n";
	}
	rates.close();

	// Saving system state for restart

	std::ofstream outputstream;
	outputstream.open(inputfolder + "/nucl_status_" + std::to_string(world_rank) + ".log");
	outputstream << std::setprecision(17);
	for(int i = 0; i < Nucl_Status.size(); i++)
	{
		outputstream << Nucl_Status[i] << "\n";
	}
	outputstream.close();

	outputstream.open(inputfolder + "/barr_file_" + std::to_string(world_rank) + ".log");
	outputstream << std::setprecision(17);
	for(int i = 0; i < Barr.size(); i++)
	{
		outputstream << Barr_Status[i] << " " << Barr_Phi[i] << "\n";
	}
	outputstream.close();

	outputstream.open(inputfolder + "/phi_file_" + std::to_string(world_rank) + ".log");
	outputstream << std::setprecision(17);
	for(int i = 0; i < numgenes; i++)
	{
		for(int j = 0; j < phi[i].size(); j++)
		{
			outputstream << phi[i][j] << " ";
		}
		outputstream << "\n";
	}
	outputstream.close();

	outputstream.open(inputfolder + "/x_file_" + std::to_string(world_rank) + ".log");
	outputstream << std::setprecision(17);
	for(int i = 0; i < numgenes; i++)
	{
		for(int j = 0; j < x[i].size(); j++)
		{
			outputstream << x[i][j] << " ";
		}
		outputstream << "\n";
	}
	outputstream.close();

	outputstream.open(inputfolder + "/intimes_file_" + std::to_string(world_rank) + ".log");
	outputstream << std::setprecision(17);
	for(int i = 0; i < numgenes; i++)
	{
		for(int j = 0; j < intimes[i].size(); j++)
		{
			outputstream << intimes[i][j] << " ";
		}
		outputstream << "\n";
	}
	outputstream.close();

	outputstream.open(inputfolder + "/outtimes_file_" + std::to_string(world_rank) + ".log");
	outputstream << std::setprecision(17);
	for(int i = 0; i < numgenes; i++)
	{
		for(int j = 0; j < outtimes[i].size(); j++)
		{
			outputstream << outtimes[i][j] << " ";
		}
		outputstream << "\n";
	}
	outputstream.close();

	outputstream.open(inputfolder + "/outpos_file_" + std::to_string(world_rank) + ".log");
	outputstream << std::setprecision(17);
	for(int i = 0; i < numgenes; i++)
	{
		for(int j = 0; j < outpos[i].size(); j++)
		{
			outputstream << outpos[i][j] << " ";
		}
		outputstream << "\n";
	}
	outputstream.close();

	outputstream.open(inputfolder + "/curr_count_file_" + std::to_string(world_rank) + ".log");
	outputstream << std::setprecision(17);
	for(int i = 0; i < numgenes; i++)
	{
		outputstream << current_rna_count[i] << "\n";
	}
	outputstream.close();

	outputstream.open(inputfolder + "/time_" + std::to_string(world_rank) + ".log");
	outputstream << std::setprecision(17);
	outputstream << t << " " << exitflag << " " << check_flag << "\n";
	outputstream.close();

	outputstream.open(inputfolder + "/generator_state_" + std::to_string(world_rank) + ".log");
	outputstream << generator;
	outputstream.close();

	event_file.close();

	std::ofstream approx_diff_file(outputfolder + "/approx_diff_" + std::to_string(world_rank) + ".log");
	approx_diff_file << approx_count << "\t" << max_approx_diff << "\n";
	approx_diff_file.close();

	if(fileflag == 1)
	{
		alive_file.close();
		x_file.close();
		phi_file.close();
		velocity_file.close();
		angular_v_file.close();
		segments_file.close();
		sigma_file.close();
		torque_file.close();
		nucl_status_file.close();
		barr_status_file.close();
		barr_phi_file.close();
		myong_file.close();
	}
	if(galactose_switch > -1)
	{
		gapr_file.close();
		spot_sigma_file.close();
	}
	if(nucl_file_flag == 1)
	{
		nucl_status_file.close();
		nucl_type_file.close();
		rnap_density_file.close();
	}
	sua_file.close();

	return;
}
