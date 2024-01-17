std::vector<double> Read_Torque_Calc_Vectors(std::string filename) // Read files containing data for torque interpolation
{
	std::ifstream f;
	f.open(filename);
	if(!f.good())
	{
		std::cout << "One or more files needed for torque calculation is missing." << "\n";
		throw std::invalid_argument("One or more files needed for torque calculation is missing.");
	}
	std::vector<double> vec;
	double val = 0.0;
	while(f >> val)
	{
		vec.push_back(val);
	}
	f.close();

	return(vec);
}

InterpMultilinear<2, double> Setup_Interp() // Return object for 2D (psi, sigma) torque interpolation
{
	std::vector<double> grid1 = Read_Torque_Calc_Vectors("../torque_interp/psi.log");
	std::vector<double> grid2 = Read_Torque_Calc_Vectors("../torque_interp/sigma.log");
	std::vector<double> f_values = Read_Torque_Calc_Vectors("../torque_interp/torque.log");

	std::vector< std::vector<double>::iterator > grid_iter_list;
	grid_iter_list.push_back(grid1.begin());
	grid_iter_list.push_back(grid2.begin());

	array<int,2> grid_sizes;
	grid_sizes[0] = grid1.size();
	grid_sizes[1] = grid2.size();

	int num_elements = grid_sizes[0]*grid_sizes[1];

	InterpMultilinear<2, double> interp_ML(grid_iter_list.begin(), grid_sizes.begin(), f_values.data(), f_values.data() + num_elements);

	return(interp_ML);
}

InterpMultilinear<1, double> Setup_Interp_Cutoffs(std::string filename) // Return object for 1D (psi) sigma_s / sigma_p cutoff interpolation
{
	std::vector<double> grid1 = Read_Torque_Calc_Vectors("../torque_interp/psi.log");
	std::vector<double> f_values = Read_Torque_Calc_Vectors(filename);

	std::vector<std::vector<double>::iterator> grid_iter_list;
	grid_iter_list.push_back(grid1.begin());

	array<int, 1> grid_sizes;
	grid_sizes[0] = grid1.size();

	int num_elements = grid_sizes[0];

	InterpMultilinear<1, double> interp_ML(grid_iter_list.begin(), grid_sizes.begin(), f_values.data(), f_values.data() + num_elements);

	return(interp_ML);
}

double Get_Barrier_On_Rate(double base_rate, double sigma) // Possible dependence of barrier on rate on local sigma
{
	return(base_rate);
}

double Get_Barrier_Off_Rate(double base_rate, double sigma) // Possible dependence of barrier off rate on local sigma
{
	return(base_rate);
}

double Get_GapR_On_Rate(double base_rate, double sigma, int galactose_switch) // GapR binding rate dependent on local sigma
{
	if(galactose_switch == 0)
	{
		return(0.0);
	}
	return(((1.0 + std::tanh((sigma - 0.5) / 0.25)) / 2.0)*base_rate);
}

int Is_Barr_Site_Blocked(double barr_loc, const std::vector<double> &phi_x, const std::vector<int> &Alive, double gap_delta) // Find if a barrier protein binding site is blocked by the presence of an RNAP
{
	int status = 1;
	int PolII_Count = (phi_x.size() - 1) / 2;
	for(int i = 0; i < PolII_Count; i++)
	{
		if(Alive[i] == 0)
		{
			continue;
		}
		if(std::abs(phi_x[i + PolII_Count] - barr_loc) < gap_delta)
		{
			status = 0;
			break;
		}
	}

	return(status);
}

double Get_Prokaryotic_Torque(double sigma, double f, double segment) // Torque(sigma) calculation for the prokaryotic case
{
	double cutoff = 0.34*1000.0;

	double A = 50.0, C = 95.0, P = 24.0;
	double kBT = 4.1; // T = 300K.
	double A_m = 4.0, C_m = 1.75, e_m = 6.0*kBT, sigma_0 = -1.0;

	double c = kBT*C*w0*w0;
	double p = kBT*P*w0*w0;

	double cs = c*(1.0 - ((C / (4.0*A))*std::sqrt(kBT / (A*f))));
	double g = f - std::sqrt((kBT*f / A));

	double factor = std::sqrt((2.0*p*g) / (1.0 - (p / cs)));

	double sigma_s = factor / cs;
	if(finite_size_flag == 1)
	{
		sigma_s = sigma_s*(1.0 + std::pow((cutoff / segment), 2.0));
	}
	double sigma_p = factor / p;

	double g_m = 1.2*(f - std::sqrt(kBT*f / A_m));
	double c_m = kBT*C_m*w0*w0;

	double sigma_sm = (c_m / (cs - c_m))*(-sigma_0 - std::sqrt(sigma_0*sigma_0 + (2.0*(cs - c_m) / (cs*c_m))*(g + e_m - g_m)));
	double sigma_m = sigma_0 + (cs / (cs - c_m))*(-sigma_0 - std::sqrt(sigma_0*sigma_0 + (2.0*(cs - c_m) / (cs*c_m))*(g + e_m - g_m)));

	double torque = 0.0;

	if(sigma <= sigma_m)
	{
		torque = (c_m / w0)*(sigma - sigma_0);
	}
	else if(sigma > sigma_m && sigma <= sigma_sm)
	{
		torque = (c_m / w0)*(sigma_m - sigma_0);
	}
	else if(sigma > sigma_sm && sigma <= sigma_s)
	{
		torque = (cs / w0)*sigma;
	}
	else if(sigma > sigma_s && sigma <= sigma_p)
	{
		torque = factor / w0;
	}
	else if(sigma > sigma_p)
	{
		torque = (p / w0)*sigma;
		if(torque > 40.0)
		{
			torque = 40.0;
		}
	}

	return(torque);
}

double Get_Eukaryotic_Torque(InterpMultilinear<2, double> &interp_ML, InterpMultilinear<1, double> &interp_ML_Cutoff_0, InterpMultilinear<1, double> &interp_ML_Cutoff_1, double nucl_density, double sigma, double force, double segment) // Torque(psi, sigma) calculation for the eukaryotic case
{
	double cutoff = 0.34*1000.0, slope = 747.65; // Slope is for the linear torque in the "DNA being twisted" regime

	array<double, 2> args = {nucl_density, sigma};
	array<double, 2> args_base;
	array<double, 1> args_2 = {nucl_density};

	double sigma_s = 0.0, sigma_p = 0.0, sigma_p_star = 0.0;

	if(finite_size_flag == 1)
	{
		sigma_s = interp_ML_Cutoff_0.interp(args_2.begin());
		args_base = {nucl_density, sigma_s};
		sigma_p = interp_ML_Cutoff_1.interp(args_2.begin());
		sigma_p_star = sigma_p*(1.0 + std::pow((cutoff / segment), 2.0));

		if(sigma > sigma_s && sigma < sigma_p_star)
		{
			return((slope)*(sigma - sigma_s) + interp_ML.interp(args_base.begin()));
		}
	}

	if(sigma < -0.35)
	{
		return(-10.0026);
	}
	if(sigma > 0.35)
	{
		return(40.0);
	}

	return(interp_ML.interp(args.begin()));
}

std::vector<double> Get_GapR_Points(double clamp0, double clamp1, double gap) // Get a list of points at which the supercoiling density will be recorded
{
	std::vector<double> GapR_Points;
	double pointer = clamp0;
	while(pointer < clamp1)
	{
		pointer += gap;
		if(pointer < clamp1)
		{
			GapR_Points.push_back(pointer);
		}
	}

	return(GapR_Points);
}

std::vector<double> Get_GapR_Points_Gene_Bodies(const std::vector<double> &gene_start, const std::vector<double> &gene_end, const std::vector<double> &gene_lengths) // Get a list of points, only within gene bodies, where the supercoiling density will be probed
{
	std::vector<double> GapR_Points;
	double pointer = 0.0, gap = 0.0;
	for(int i = gene_lengths.size() - 1; i >= 0; i--)
	{
		gap = gene_lengths[i] / 50.0;
		if(gene_start[i] < gene_end[i])
		{
			pointer = gene_start[i];
			while(pointer < gene_end[i] + gap)
			{
				GapR_Points.push_back(pointer);
				pointer += gap;
			}
		}
		else if(gene_start[i] > gene_end[i])
		{
			pointer = gene_end[i];
			while(pointer < gene_start[i] + gap)
			{
				GapR_Points.push_back(pointer);
				pointer += gap;
			}
		}
		else
		{
			std::cout << "Spotted a contradiction in the gene start and end sites." << "\n";
			throw std::invalid_argument("Spotted a contradiction in the gene start and end sites.");
		}
	}

	for(int i = 0; i < GapR_Points.size() - 1; i++)
	{
		if(GapR_Points[i + 1] < GapR_Points[i])
		{
			std::cout << "Points chosen in the gene bodies are incorrectly ordered." << "\n";
			throw std::invalid_argument("Points chosen in the gene bodies are incorrectly ordered.");
		}
	}

	return(GapR_Points);
}

int Is_TSS_Blocked_By_Barrier(const std::vector<double> &TSSes, const std::vector<double> &Barr, const std::vector<int> &Barr_Status, int TSS_index, int *blocking_barr, double gap_delta) // Find if the TSS (given by TSS_index) is blocked by the binding of a barrier protein
{
	double dist = 0.0;
	for(int i = 0; i < Barr.size(); i++)
	{
		if(Barr_Status[i] == 0)
		{
			continue;
		}
		dist = std::abs(Barr[i] - TSSes[TSS_index]);
		if(dist < gap_delta)
		{
			(*blocking_barr) = i;
			return(1);
		}
	}

	return(0);
}

void Get_Nucl_RNAP_Map_Pairwise(const std::vector<double> &Nucl, const std::vector<double> &allx, std::vector<std::vector<int>> &map, double gap_delta, int brute_force_flag) // Find, for each nucleosome, the list of RNAPs or barrier proteins within a distance of gap_delta
{
	for(int i = 0; i < map.size(); i++)
	{
		map.clear();
	}
	for(int i = 0; i < Nucl.size(); i++)
	{
		map.push_back(std::vector<int>());
	}
	int flag = 0, mapped_counts = 0;
	for(int i = 0; i < Nucl.size(); i++)
	{
		flag = 0;
		for(int j = 0; j < allx.size(); j++)
		{
			if(Nucl[i] > allx[j] && Nucl[i] - allx[j] < gap_delta)
			{
				if(Nucl[i] - allx[j] < 0.0)
				{
					std::cout << "First check condition failure in Get_Nucl_RNAP_Map_Pairwise" << "\n";
					throw std::invalid_argument("First check condition failure in Get_Nucl_RNAP_Map_Pairwise");
				}
				flag = 1;
				map[i].push_back(j);
				mapped_counts += 1;
			}
			else if(Nucl[i] < allx[j] && allx[j] - Nucl[i] < nucleosome + gap_delta)
			{
				if(allx[j] - Nucl[i] < 0.0)
				{
					std::cout << "Second check condition failure in Get_Nucl_RNAP_Map_Pairwise" << "\n";
					throw std::invalid_argument("Second check condition failure in Get_Nucl_RNAP_Map_Pairwise");
				}
				flag = 1;
				map[i].push_back(j);
				mapped_counts += 1;
			}
		}
	}
}

void Get_GapR_RNAP_Map_Pairwise(const std::vector<double> &GapR_Points, const std::vector<double> &allx, std::vector<std::vector<int>> &map, double gap_delta) // Find, for each GapR site, the list of RNAPs or barrier proteins within a distance of gap_delta
{
	for(int i = 0; i < map.size(); i++)
	{
		map.clear();
	}
	for(int i = 0; i < GapR_Points.size(); i++)
	{
		map.push_back(std::vector<int>());
	}
	int flag = 0, mapped_counts = 0;
	for(int i = 0; i < GapR_Points.size(); i++)
	{
		flag = 0;
		for(int j = 0; j < allx.size(); j++)
		{
			if(GapR_Points[i] > allx[j] && GapR_Points[i] - allx[j] < gap_delta)
			{
				if(GapR_Points[i] - allx[j] < 0.0)
				{
					std::cout << "First check condition failure in Get_GapR_RNAP_Map_Pairwise" << "\n";
					throw std::invalid_argument("First check condition failure in Get_GapR_RNAP_Map_Pairwise");
				}
				flag = 1;
				map[i].push_back(j);
				mapped_counts += 1;
			}
			else if(GapR_Points[i] < allx[j] && allx[j] - GapR_Points[i] < gap_delta)
			{
				if(allx[j] - GapR_Points[i] < 0.0)
				{
					std::cout << "Second check condition failure in Get_GapR_RNAP_Map_Pairwise" << "\n";
					throw std::invalid_argument("Second check condition failure in Get_GapR_RNAP_Map_Pairwise");
				}
				flag = 1;
				map[i].push_back(j);
				mapped_counts += 1;
			}
		}
	}
}

std::vector<int> Find_Blocked_Nucl(const std::vector<std::vector<int>> &map) // Find if a nucleosome binding site is blocked by an RNAP / barrier protein
{
	std::vector<int> Is_Nucl_Blocked;
	for(int i = 0; i < map.size(); i++)
	{
		if(map[i].size() > 0)
		{
			Is_Nucl_Blocked.push_back(1);
		}
		else
		{
			Is_Nucl_Blocked.push_back(0);
		}
	}

	return(Is_Nucl_Blocked);
}

std::vector<int> Find_Blocked_GapR(const std::vector<std::vector<int>> &map, int brute_gapr_flag) // Find if a GapR site is blocked by an RNAP / barrier protein
{
	std::vector<int> Is_GapR_Blocked;
	for(int i = 0; i < map.size(); i++)
	{
		if(map[i].size() > 0 && brute_gapr_flag == 0)
		{
			Is_GapR_Blocked.push_back(1);
		}
		else
		{
			Is_GapR_Blocked.push_back(0);
		}
	}

	return(Is_GapR_Blocked);
}

void Binary_Search_Nucl_Array(const std::vector<double> &Nucl_Bound, double pos, int *left_pos, int *right_pos) // Find the nucleosomes located to the left and right of pos
{
	(*left_pos) = 0;
	(*right_pos) = Nucl_Bound.size() - 1;
	int mid;
	if(Nucl_Bound.size() == 0)
	{
		(*left_pos) = -1;
		(*right_pos) = -1;

		return;
	}
	if(pos < Nucl_Bound[(*left_pos)])
	{
		(*left_pos) = -1;
		(*right_pos) = 0;

		return;
	}
	if(pos > Nucl_Bound[(*right_pos)])
	{
		(*left_pos) = Nucl_Bound.size() - 1;
		(*right_pos) = -1;

		return;
	}
	while((*left_pos) < (*right_pos) && (*left_pos) >= 0 && (*right_pos) < Nucl_Bound.size())
	{
		mid = ((*left_pos) + (*right_pos)) / 2;
		if(mid == (*left_pos) || mid == (*right_pos))
		{
			return;
		}
		if(pos > Nucl_Bound[(*left_pos)] && pos < Nucl_Bound[mid])
		{
			(*right_pos) = mid;
		}
		else if(pos > Nucl_Bound[mid] && pos < Nucl_Bound[(*right_pos)])
		{
			(*left_pos) = mid;
		}
	}

	return;
}

std::vector<int> Find_RNAP_Blocked_By_Nucl(const std::vector<double> &Nucl, const std::vector<int> &Nucl_Status, const std::vector<double> &phi_x, const std::vector<int> &Alive, const std::vector<int> &Direction, double gap_delta, int brute_force_flag) // Find which RNAPs are blocked by a nucleosome
{
	std::vector<int> Is_RNAP_Blocked;
	int PolII_Count = (phi_x.size() - 1) / 2;
	for(int i = 0; i < PolII_Count; i++)
	{
		Is_RNAP_Blocked.push_back(0);
	}

	if(brute_force_flag == 1) // if the nucleosomes do not offer a steric hindrance to RNAP movement
	{
		return(Is_RNAP_Blocked);
	}

	std::vector<double> Nucl_Bound;
	for(int i = 0; i < Nucl.size(); i++)
	{
		if(Nucl_Status[i] == 1)
		{
			Nucl_Bound.push_back(Nucl[i]);
		}
	}
	int left_pos, right_pos, flag;
	for(int i = 0; i < PolII_Count; i++)
	{
		if(Alive[i] == 0)
		{
			Is_RNAP_Blocked[i] = -1;
		}
		Binary_Search_Nucl_Array(Nucl_Bound, phi_x[i + PolII_Count], &left_pos, &right_pos);

		if(left_pos >= right_pos && left_pos != -1 && right_pos != -1)
		{
			std::cout << "First failure of the binary search process" << "\n";
			throw std::invalid_argument("First failure of the binary search process");
		}
		if(left_pos > -1 && right_pos > -1 && !(phi_x[i + PolII_Count] > Nucl_Bound[left_pos] && phi_x[i + PolII_Count] < Nucl_Bound[right_pos]))
		{
			std::cout << "Second failure of the binary search process" << "\n";
			throw std::invalid_argument("Second failure of the binary search process");
		}

		flag = 0;
		if(left_pos == -1 && right_pos == -1)
		{
			flag = 0;
		}
		else if(left_pos == -1)
		{
			if(Nucl_Bound[right_pos] - phi_x[i + PolII_Count] < gap_delta && Direction[i] == 1)
			{
				if(Nucl_Bound[right_pos] - phi_x[i + PolII_Count] < 0.0)
				{
					std::cout << "Third failure of the binary search process" << "\n";
					throw std::invalid_argument("Third failure of the binary search process");
				}
				flag = 1;
			}
		}
		else if(right_pos == -1)
		{
			if(phi_x[i + PolII_Count] - Nucl_Bound[left_pos] < nucleosome + gap_delta && Direction[i] == -1)
			{
				if(phi_x[i + PolII_Count] - Nucl_Bound[left_pos] < 0.0)
				{
					std::cout << "Fourth failure of the binary search process" << "\n";
					throw std::invalid_argument("Fourth failure of the binary search process");
				}
				flag = 1;
			}
		}
		else
		{
			if(!(Nucl_Bound[left_pos] < phi_x[i + PolII_Count] && Nucl_Bound[right_pos] > phi_x[i + PolII_Count]))
			{
				std::cout << "Fifth failure of the binary search process" << "\n";
				throw std::invalid_argument("Fifth failure of the binary search process");
			}
			if(Nucl_Bound[right_pos] - phi_x[i + PolII_Count] < gap_delta && Direction[i] == 1)
			{
				flag = 1;
			}
			if(phi_x[i + PolII_Count] - Nucl_Bound[left_pos] < nucleosome + gap_delta && Direction[i] == -1)
			{
				flag = 1;
			}
		}
		if(flag == 1)
		{
			Is_RNAP_Blocked[i] = 1;
		}
	}

	return(Is_RNAP_Blocked);
}

std::vector<int> Find_TSSes_Blocked_By_Nucl(const std::vector<double> &Nucl, const std::vector<int> &Nucl_Status, const std::vector<double> TSSes, std::vector<int> &left_block, std::vector<int> &right_block, double gap_delta, int brute_force_flag) // Find if any of the TSSes are blocked by a nucleosome
{
	for(int i = 0; i < TSSes.size(); i++)
	{
		left_block[i] = -1;
		right_block[i] = -1;
	}

	std::vector<int> Is_TSS_Blocked;
	for(int i = 0; i < TSSes.size(); i++)
	{
		Is_TSS_Blocked.push_back(0);
	}

	if(brute_force_flag == 1)
	{
		return(Is_TSS_Blocked);
	}

	std::vector<double> Nucl_Bound;
	std::vector<int> Bound_Index;
	for(int i = 0; i < Nucl.size(); i++)
	{
		if(Nucl_Status[i] == 1)
		{
			Nucl_Bound.push_back(Nucl[i]);
			Bound_Index.push_back(i);
		}
	}
	int left_pos = -1, right_pos = -1, flag;
	for(int i = 0; i < TSSes.size(); i++)
	{
		Binary_Search_Nucl_Array(Nucl_Bound, TSSes[i], &left_pos, &right_pos);

		if(left_pos >= right_pos && left_pos != -1 && right_pos != -1)
		{
			std::cout << "First failure of the binary search process in Find_TSSes_Blocked_By_Nucl" << "\n";
			throw std::invalid_argument("First failure of the binary search process in Find_TSSes_Blocked_By_Nucl");
		}
		if(left_pos > -1 && right_pos > -1 && !(TSSes[i] > Nucl_Bound[left_pos] && TSSes[i] < Nucl_Bound[right_pos]))
		{
			std::cout << "Second failure of the binary search process in Find_TSSes_Blocked_By_Nucl" << "\n";
			throw std::invalid_argument("Second failure of the binary search process in Find_TSSes_Blocked_By_Nucl");
		}

		flag = 0;
		if(left_pos == -1 && right_pos == -1)
		{
			flag = 0;
		}
		else if(left_pos == -1)
		{
			if(Nucl_Bound[right_pos] - TSSes[i] < 0.0)
			{
				std::cout << "Third failure of the binary search process in Find_TSSes_Blocked_By_Nucl" << "\n";
				throw std::invalid_argument("Third failure of the binary search process in Find_TSSes_Blocked_By_Nucl");
			}
			if(Nucl_Bound[right_pos] - TSSes[i] < gap_delta)
			{
				flag = 2;
			}
		}
		else if(right_pos == -1)
		{
			if(TSSes[i] - Nucl_Bound[left_pos] < 0.0)
			{
				std::cout << "Fourth failure of the binary search process in Find_TSSes_Blocked_By_Nucl" << "\n";
				throw std::invalid_argument("Fourth failure of the binary search process in Find_TSSes_Blocked_By_Nucl");
			}
			if(TSSes[i] - Nucl_Bound[left_pos] < nucleosome + gap_delta)
			{
				flag = 1;
			}
		}
		else
		{
			if(Nucl_Bound[right_pos] - TSSes[i] < gap_delta && TSSes[i] - Nucl_Bound[left_pos] > nucleosome + gap_delta)
			{
				flag = 2;
			}
			else if(Nucl_Bound[right_pos] - TSSes[i] < gap_delta && TSSes[i] - Nucl_Bound[left_pos] < nucleosome + gap_delta)
			{
				flag = 3;
			}
			else if(Nucl_Bound[right_pos] - TSSes[i] > gap_delta && TSSes[i] - Nucl_Bound[left_pos] > nucleosome + gap_delta)
			{
				flag = 0;
			}
			else if(Nucl_Bound[right_pos] - TSSes[i] > gap_delta && TSSes[i] - Nucl_Bound[left_pos] < nucleosome + gap_delta)
			{
				flag = 1;
			}
		}
		if(flag > 0)
		{
			Is_TSS_Blocked[i] = 1;
			if(flag == 2)
			{
				right_block[i] = Bound_Index[right_pos];
			}
			else if(flag == 1)
			{
				left_block[i] = Bound_Index[left_pos];
			}
			else if(flag == 3)
			{
				right_block[i] = Bound_Index[right_pos];
				left_block[i] = Bound_Index[left_pos];
			}
		}
	}

	return(Is_TSS_Blocked);
}

std::vector<int> Get_Segment_Nucleosome_Count_Pairwise(const std::vector<double> &allx, const std::vector<double> &Nucl, const std::vector<int> &Nucl_Status) // Get the number of nucleosomes in each segment; uses an inefficient pairwise approach
{
	std::vector<int> count, Nucl_ID;
	for(int i = 0; i < allx.size() + 1; i++)
	{
		count.push_back(0);
	}
	for(int i = 0; i < Nucl.size(); i++)
	{
		Nucl_ID.push_back(-1);
	}
	for(int i = 0; i < allx.size(); i++)
	{
		for(int j = 0; j < Nucl.size(); j++)
		{
			if(Nucl_Status[j] == 0)
			{
				continue;
			}
			if(i == 0)
			{
				if(Nucl[j] >= allx[0])
				{
					count[0] += 1;
					Nucl_ID[j] = 0;
				}
			}
			if(i == allx.size() - 1)
			{
				if(Nucl[j] < allx[i])
				{
					count[allx.size()] += 1;
					Nucl_ID[j] = allx.size();
				}
			}
			if(i > 0)
			{
				if(Nucl[j] >= allx[i] && Nucl[j] < allx[i - 1])
				{
					count[i] += 1;
					Nucl_ID[j] = i;
				}
			}
		}
	}
	if(allx.size() == 0)
	{
		for(int i = 0; i < Nucl.size(); i++)
		{
			if(Nucl_Status[i] == 0)
			{
				continue;
			}
			count[0] += 1;
			Nucl_ID[i] = 0;
		}
	}
	int count0 = 0, count1 = 0;
	for(int i = 0; i < Nucl.size(); i++)
	{
		if(Nucl_Status[i] == 1)
		{
			count0 += 1;
		}
	}
	for(int i = 0; i < count.size(); i++)
	{
		count1 += count[i];
	}
	if(count0 != count1)
	{
		std::cout << "Mismatch between the total nucleosome count and the numbers in different segments" << "\n";
//		std::cout << count0 << "\t" << count1 << "\n";
//		for(int i = 0; i < count.size(); i++)
//		{
//			std::cout << count[i] << "\t";
//		}
//		std::cout << "\n";
//		for(int i = 0; i < allx.size(); i++)
//		{
//			std::cout << allx[i] << "\t";
//		}
//		std::cout << "\n";
//		for(int i = 0; i < Nucl.size(); i++)
//		{
//			if(Nucl_Status[i] == 0)
//			{
//				continue;
//			}
//			if(Nucl_ID[i] == -1)
//			{
//				std::cout << "--->\t" << Nucl[i] << "\n";
//			}
//			std::cout << Nucl[i] << "\t";
//		}
//		std::cout << "\n";
//		for(int i = 0; i < Nucl_ID.size(); i++)
//		{
//			if(Nucl_ID[i] == -1)
//			{
//				continue;
//			}
//			std::cout << Nucl_ID[i] << "\t";
//		}
//		std::cout << "\n";
		throw std::invalid_argument("Mismatch between the total nucleosome count and the numbers in different segments");
	}

	return(count);
}

std::vector<int> Get_Segment_RNAP_Count(const std::vector<double> &pos, const std::vector<double> &phi_x, const std::vector<int> &Alive) // Find the number of RNAPs in each segment (segments defined by the vector pos); for the calculation of RNAP density
{
	int PolII_Count = (phi_x.size() - 1) / 2;

	std::vector<int> count;
	for(int i = 0; i < pos.size() + 1; i++)
	{
		count.push_back(0);
	}
	for(int i = 0; i < pos.size() + 1; i++)
	{
		for(int j = 0; j < PolII_Count; j++)
		{
			if(Alive[j] == 0)
			{
				continue;
			}
			if(i == 0)
			{
				if(phi_x[PolII_Count + j] <= pos[0])
				{
					count[0] += 1;
				}
			}
			else if(i == pos.size())
			{
				if(phi_x[PolII_Count + j] > pos[pos.size() - 1])
				{
					count[pos.size()] += 1;
				}
			}
			else
			{
				if(phi_x[PolII_Count + j] > pos[i - 1] && phi_x[PolII_Count + j] <= pos[i])
				{
					count[i] += 1;
				}
			}
		}
	}

	int count0 = 0, count1 = 0;
	for(int i = 0; i < PolII_Count; i++)
	{
		if(Alive[i] == 1)
		{
			count0 += 1;
		}
	}

	for(int i = 0; i < count.size(); i++)
	{
		count1 += count[i];
	}

	if(count0 != count1)
	{
		std::cout << "Mismatch betweeen the total RNAP count and the numbers in different segments" << "\n";
		throw std::invalid_argument("Mismatch betweeen the total RNAP count and the numbers in different segments");
	}

	return(count);
}

std::vector<int> Is_Alive(const std::vector<double> &phi_x, const std::vector<double> &TSSes, const std::vector<double> &gene_lengths, const std::vector<int> &Gene_ID) // Find if an RNAP has finished transcribing
{
	int PolII_Count = (phi_x.size() - 1) / 2;
	std::vector<int> Alive;
	for(int i = 0; i < PolII_Count; i++)
	{
		if(std::abs(phi_x[i + PolII_Count] - TSSes[Gene_ID[i]]) - gene_lengths[Gene_ID[i]] > 0.0)
		{
			Alive.push_back(0); // RNAP has finished transcribing
		}
		else
		{
			Alive.push_back(1); // RNAP is still within the gene body
		}
	}

	return(Alive);
}

template<typename T> std::vector<T> Get_All_Of_Two_Vectors_Ordered(const std::vector<T> vec1, const std::vector<T> vec2, int shift1, int shift2, const std::vector<double> &phi_x, const std::vector<int> &Alive, const std::vector<double> &Barr, const std::vector<int> &Barr_Status, int caller_id) //Put entries from vectors vec1 and vec2 in the order of positions of RNAPs and barrier proteins
{
	int PolII_Count = (phi_x.size() - 1) / 2, Alive_PolII_Count = 0, Alive_Barr_Count = 0;
	for(int i = 0; i < PolII_Count; i++)
	{
		if(Alive[i] == 1)
		{
			Alive_PolII_Count += 1;
		}
	}
	for(int i = 0; i < Barr.size(); i++)
	{
		if(Barr_Status[i] == 1)
		{
			Alive_Barr_Count += 1;
		}
	}
	std::vector<T> allvec;
	int index0 = 0, index1 = 0;
	for(int i = 0; i < Alive_PolII_Count + Alive_Barr_Count; i++)
	{
		while(index0 < Alive.size() && Alive[index0] == 0)
		{
			index0 += 1;
		}
		while(index1 < Barr.size() && Barr_Status[index1] == 0)
		{
			index1 += 1;
		}
		if(index1 >= Barr.size() || (PolII_Count > 0 && index1 < Barr.size() && phi_x[index0 + PolII_Count] > Barr[index1]))
		{
			allvec.push_back(vec1[index0 + shift1]);
			index0 += 1;
		}
		else if(index0 >= Alive.size() || (PolII_Count > 0 && index1 < Barr.size() && Barr[index1] > phi_x[index0 + PolII_Count]))
		{
			allvec.push_back(vec2[index1 + shift2]);
			index1 += 1;
		}
	}
	int errorflag = 0;
	if(allvec.size() != Alive_PolII_Count + Alive_Barr_Count)
	{
		errorflag = 1;
	}
	if(errorflag == 1)
	{
		std::cout << "Mismatch between the sum of vector lengths and the length of the final vector; calling_fn id:" << "\t" << caller_id << "\n";
		throw std::invalid_argument("Mismatch between the sum of vector lengths and the length of the final vector; calling_fn id:" + std::to_string(caller_id));
	}

	return(allvec);
}

std::vector<double> Get_Segments_With_Barriers(const std::vector<double> &phi_x, const std::vector<int> &Alive, const std::vector<double> &Barr, const std::vector<int> &Barr_Status) // Get the lengths of the segments the genomic region is divided into by RNAPs and barrier proteins
{
	int PolII_Count = (phi_x.size() - 1) / 2, Alive_PolII_Count = 0, Alive_Barr_Count = 0;
	if(Alive.size() != PolII_Count || Barr.size() != Barr_Status.size())
	{
		std::cout << "Status vectors are of incorrect size" << "\n";
		throw std::invalid_argument("Status vectors are of incorrect size");
	}
	for(int i = 0; i < PolII_Count; i++)
	{
		if(Alive[i] == 1)
		{
			Alive_PolII_Count += 1;
		}
	}
	for(int i = 0; i < Barr.size(); i++)
	{
		if(Barr_Status[i] == 1)
		{
			Alive_Barr_Count += 1;
		}
	}
	std::vector<double> allx = Get_All_Of_Two_Vectors_Ordered<double>(phi_x, Barr, PolII_Count, 0, phi_x, Alive, Barr, Barr_Status, 0);
	for(int i = 0; i < ((int)(allx.size()) - 1); i++)
	{
		if(allx[i] < allx[i + 1])
		{
			std::cout << "Ordering function failed" << "\n";
			throw std::invalid_argument("Ordering function failed");
		}
	}

	std::vector<double> Segments;
	if(Alive_PolII_Count == 0 && Alive_Barr_Count == 0)
	{
		Segments.push_back(clamp1 - clamp0);
		return(Segments);
	}
	double right_segment, left_segment;
	for(int i = 0; i < allx.size(); i++)
	{
		if(i == 0)
		{
			right_segment = clamp1 - allx[i];
		}
		else
		{
			right_segment = allx[i - 1] - allx[i];
		}
		Segments.push_back(right_segment);
	}
	left_segment = allx[allx.size() - 1] - clamp0;
	Segments.push_back(left_segment);

	if(Segments.size() != allx.size() + 1 && Segments.size() > 0)
	{
		std::cout << "Mismatch between allx vector size and number of segments found" << "\n";
		throw std::invalid_argument("Mismatch between allx vector size and number of segments found");
	}
	double sum = 0.0;
	for(int i = 0; i < Segments.size(); i++)
	{
		sum += Segments[i];
	}
	if(sum - (clamp1 - clamp0) > 1e-3)
	{
		std::cout << "Total segment lengths do not match the overall length of the genomic segment" << "\n";
		throw std::invalid_argument("Total segment lengths do not match the overall length of the genomic segment");
	}
	for(int i = 0; i < Segments.size(); i++)
	{
		if(Segments[i] < 0.0)
		{
			std::cout << "Negative segment length encountered" << "\n";
			throw std::invalid_argument("Negative segment length encountered");
		}
	}

	return(Segments);
}

std::vector<double> Get_Supercoiling_Densities_With_Barriers(int clamp0_flag, int clamp1_flag, const std::vector<double> &phi_x, const std::vector<double> &Segments, const std::vector<int> &Alive, const std::vector<double> &Barr, const std::vector<int> &Barr_Status, const std::vector<double> &Barr_Phi) // Get the supercoiling densities in the segments defined by RNAPs and protein barriers
{
	int PolII_Count = (phi_x.size() - 1) / 2, segmentindex = 0;
	if(Alive.size() != PolII_Count || Barr.size() != Barr_Status.size())
	{
		std::cout << "Status vectors do not have the correct size; in supercoiling density calc" << "\n";
		throw std::invalid_argument("Status vectors do not have the correct size; in supercoiling density calc");
	}

	std::vector<double> allphi = Get_All_Of_Two_Vectors_Ordered<double>(phi_x, Barr_Phi, 0, 0, phi_x, Alive, Barr, Barr_Status, 1);

	if(allphi.size() + 1 != Segments.size())
	{
		std::cout << "Mismatch between allphi vector size and number of segments" << "\n";
		throw std::invalid_argument("Mismatch between allphi vector size and number of segments");
	}

	std::vector<double> Sigma;
	if(allphi.size() == 0)
	{
		if(circular_plasmid_flag == 0)
		{
			Sigma.push_back((clamp0_sigma - clamp1_sigma) / (w0*Segments[0]));
		}
		else if(circular_plasmid_flag == 1)
		{
			Sigma.push_back(0.0 + sigma_basal);
		}
		return(Sigma);
	}
	double right_sigma, left_sigma;
	for(int i = 0; i < allphi.size(); i++)
	{
		if(i == 0)
		{
			if(circular_plasmid_flag == 0)
			{
				if(clamp1_flag == 0)
				{
					right_sigma = 0.0;
				}
				else if(clamp1_flag == 1)
				{
					right_sigma = (allphi[i] - clamp1_sigma) / (w0*Segments[segmentindex]);
				}
			}
			else if(circular_plasmid_flag == 1)
			{
				right_sigma = (allphi[i] - allphi[allphi.size() - 1]) / (w0*(Segments[0] + Segments[Segments.size() - 1]));
			}
			Sigma.push_back(right_sigma);
		}
		else
		{
			right_sigma = (allphi[i] - allphi[i - 1]) / (w0*Segments[segmentindex]);
			Sigma.push_back(right_sigma);
		}
		segmentindex += 1;
	}
	if(allphi.size() == 0)
	{
		;
	}
	else
	{
		if(circular_plasmid_flag == 0)
		{
			if(clamp0_flag == 0)
			{
				left_sigma = 0.0;
			}
			else
			{
				left_sigma = (clamp0_sigma - allphi[allphi.size() - 1]) / (w0*Segments[Segments.size() - 1]);
			}
		}
		else if(circular_plasmid_flag == 1)
		{
			left_sigma = (allphi[0] - allphi[allphi.size() - 1]) / (w0*(Segments[0] + Segments[Segments.size() - 1]));
		}
		Sigma.push_back(left_sigma);
	}

	if(circular_plasmid_flag == 1)
	{
		for(int i = 0; i < Sigma.size(); i++)
		{
			Sigma[i] += sigma_basal;
		}
	}

	if(circular_plasmid_flag == 1 && Sigma.size() > 1)
	{
		if(std::abs(Sigma[0] - Sigma[Sigma.size() - 1]) > 1e-12)
		{
			std::cout << "Moorhead, Minnesota" << "\n";
			throw std::invalid_argument("Moorhead, Minnesota");
		}
	}

	if(Sigma.size() != Segments.size())
	{
		std::cout << "Sigma vector does not match number of segments" << "\n";
		throw std::invalid_argument("Sigma vector does not match number of segments");
	}

	return(Sigma);
}

std::vector<double> Get_Prokaryotic_Torques(const std::vector<double> &Sigma, const std::vector<double> &Segments, double force) // Get torque in each segment for the case of prokaryotes
{
	std::vector<double> Torques;
	for(int i = 0; i < Sigma.size(); i++)
	{
		Torques.push_back(Get_Prokaryotic_Torque(Sigma[i], force, Segments[i]));
	}

	return(Torques);
}

std::vector<double> Get_Interpolated_Eukaryotic_Torques(InterpMultilinear<2, double> &interp_ML, InterpMultilinear<1, double> &interp_ML_Cutoff_0, InterpMultilinear<1, double> &interp_ML_Cutoff_1, const std::vector<double> &Sigma, const std::vector<double> &Segments, const std::vector<int> &nucl_seg_counts, double force) // Get torque in each segment for the case of eukaryotes
{
	std::vector<double> Torques;
	double nucl_density = 0.0;
	for(int i = 0; i < Sigma.size(); i++)
	{
		nucl_density = ((double) nucl_seg_counts[i])*((nucleosome + linker) / Segments[i]);
		Torques.push_back(Get_Eukaryotic_Torque(interp_ML, interp_ML_Cutoff_0, interp_ML_Cutoff_1, nucl_density, Sigma[i], force, Segments[i]));
	}

	return(Torques);
}

std::vector<double> Get_Velocities_With_Barriers(const std::vector<double> &phi_x, const std::vector<double> &Segments, const std::vector<double> &Torques, const std::vector<int> &Direction, const std::vector<int> &Alive, const std::vector<double> &Barr, const std::vector<int> &Barr_Status, const std::vector<int> &Is_RNAP_Blocked, double gap_delta) // Get RNAP velocities
{
	int PolII_Count = (phi_x.size() - 1) / 2, segmentindex = 0;
	if(Alive.size() != PolII_Count || Barr.size() != Barr_Status.size())
	{
		std::cout << "Wrong size of status vectors during velocity calculation" << "\n";
		throw std::invalid_argument("Wrong size of status vectors during velocity calculation");
	}
	std::vector<int> PolII_ID, Barr_ID;
	for(int i = 0; i < PolII_Count; i++)
	{
		PolII_ID.push_back(1);
	}
	for(int i = 0; i < Barr.size(); i++)
	{
		Barr_ID.push_back(2);
	}

	std::vector<int> allid = Get_All_Of_Two_Vectors_Ordered<int>(PolII_ID, Barr_ID, 0, 0, phi_x, Alive, Barr, Barr_Status, 2);

	std::vector<double> Velocities;
	double right_torque, left_torque, v;
	while(allid.size() > 0 && allid[segmentindex] == 2)
	{
		segmentindex += 1;
	}
	for(int i = 0; i < PolII_Count; i++)
	{
		if(Alive[i] == 0)
		{
			Velocities.push_back(0.0);
			continue;
		}
		if(allid[segmentindex] != 1)
		{
			std::cout << "First failure with allid during velocity calc" << "\n";
			throw std::invalid_argument("First failure with allid during velocity calc");
		}
		if(segmentindex + 1 > Torques.size() - 1)
		{
			std::cout << "Segment index out of bounds during velocity calculation" << "\n";
			throw std::invalid_argument("Segment index out of bounds during velocity calculation");
		}
		right_torque = Torques[segmentindex];
		left_torque = Torques[segmentindex + 1];
		if(Direction[i] == 1)
		{
			if(Segments[segmentindex] > gap_delta && Is_RNAP_Blocked[i] == 0)
			{
				v = (v0 / 2.0)*(1.0 - std::tanh((right_torque - left_torque) / tau_c));
			}
			else
			{
				v = 0.0;
			}
		}
		else if(Direction[i] == -1)
		{
			if(Segments[segmentindex + 1] > gap_delta && Is_RNAP_Blocked[i] == 0)
			{
				v = (-v0 / 2.0)*(1.0 - std::tanh(-(right_torque - left_torque) / tau_c));
			}
			else
			{
				v = 0.0;
			}
		}
		else
		{
			std::cout << "Invalid gene direction" << "\n";
			throw std::invalid_argument("Invalid gene direction");
		}
		Velocities.push_back(v);
		segmentindex += 1;
		while(allid[segmentindex] == 2)
		{
			segmentindex += 1;
		}
	}
	if(Velocities.size() != PolII_Count)
	{
		std::cout << "Velocity vector size does not match RNAP count" << "\n";
		throw std::invalid_argument("Velocity vector size does not match RNAP count");
	}

	return(Velocities);
}

std::vector<double> Get_Nucleosome_Array(double clamp0, double clamp1) // Get a list of nucleosome binding sites
{
	std::vector<double> Nucl;
	double startpos = clamp0 + (linker / 2.0);
	while(startpos < clamp1)
	{
		if(startpos + nucleosome > clamp1)
		{
			break;
		}
		Nucl.push_back(startpos);
		startpos += nucleosome + linker;
	}

	return(Nucl);
}

double Get_Equalizing_Phi(const std::vector<double> &Segments, const std::vector<double> &Barr, double eqpos, const std::vector<double> &phi_x, const std::vector<int> &Alive, const std::vector<int> &Barr_Status, const std::vector<double> &Barr_Phi, int clamp0_flag, int clamp1_flag) // Get the value of phi in the middle of a segment such that if an RNAP / barrier protein were placed at that point, the two segments on either side will have the same supercoiling density
{
	int segmentindex = Segments.size() - 1, targetseg = 0, flag = 0;
	double loc = clamp0 + Segments[segmentindex], oldloc = clamp0;
	while(segmentindex >= 0)
	{
		if(eqpos > oldloc && eqpos < loc)
		{
			targetseg = segmentindex;
			flag = 1;
			break;
		}
		oldloc = loc;
		segmentindex -= 1;
		loc += Segments[segmentindex];
	}
	if(flag == 0)
	{
		targetseg = 0;
	}

	int PolII_Count = (phi_x.size() - 1) / 2;
	segmentindex = 0;
	if(Alive.size() != PolII_Count || Barr.size() != Barr_Status.size())
	{
		std::cout << "Wrong status vector size in Get_Equalizing_Phi" << "\n";
		throw std::invalid_argument("Wrong status vector size in Get_Equalizing_Phi");
	}
	std::vector<double> allx = Get_All_Of_Two_Vectors_Ordered<double>(phi_x, Barr, PolII_Count, 0, phi_x, Alive, Barr, Barr_Status, 3);

	for(int i = 0; i < ((int)(allx.size()) - 1); i++)
	{
		if(allx[i] < allx[i + 1])
		{
			std::cout << "Failure in allx vector assembly in Get_Equalizing_Phi" << "\n";
			throw std::invalid_argument("Failure in allx vector assembly in Get_Equalizing_Phi");
		}
	}

	std::vector<double> allphi = Get_All_Of_Two_Vectors_Ordered<double>(phi_x, Barr_Phi, 0, 0, phi_x, Alive, Barr, Barr_Status, 4);

	if(allphi.size() + 1 != Segments.size() || allx.size() + 1 != Segments.size())
	{
		std::cout << "Mismatch between allx / allphi and segment vector size" << "\n";
		throw std::invalid_argument("Mismatch between allx / allphi and segment vector size");
	}

	if(circular_plasmid_flag == 1 && allphi.size() == 1)
	{
		return(allphi[0]);
	}

	double rightphi, leftphi, rightx, leftx;
	if(allphi.size() == 0)
	{
		rightphi = clamp1_sigma;
		leftphi = clamp0_sigma;
		rightx = clamp1;
		leftx = clamp0;
	}
	else
	{
		if(targetseg == 0)
		{
			rightphi = clamp1_sigma;
			leftphi = allphi[0];
			if(circular_plasmid_flag == 0 && clamp1_flag == 0)
			{
				rightphi = leftphi;
			}
			rightx = clamp1;
			leftx = allx[0];
			if(circular_plasmid_flag == 1)
			{
				rightphi = allphi[allphi.size() - 1];
				leftphi = allphi[0];
				rightx = allx[allx.size() - 1];
				leftx = allx[0];

				return(rightphi + ((leftphi - rightphi)*(clamp1 - eqpos + rightx) / (clamp1 - leftx + rightx)));
			}
		}
		else if(targetseg == Segments.size() - 1)
		{
			rightphi = allphi[allphi.size() - 1];
			leftphi = clamp0_sigma;
			if(circular_plasmid_flag == 0 && clamp0_flag == 0)
			{
				leftphi = rightphi;
			}
			rightx = allx[allx.size() - 1];
			leftx = clamp0;
			if(circular_plasmid_flag == 1)
			{
				rightphi = allphi[0];
				leftphi = allphi[allphi.size() - 1];
				rightx = allx[0];
				leftx = allx[allx.size() - 1];

				return(rightphi - ((rightphi - leftphi)*(clamp1 - rightx + eqpos) / (clamp1 - rightx + leftx)));
			}
		}
		else
		{
			rightphi = allphi[targetseg - 1];
			leftphi = allphi[targetseg];
			rightx = allx[targetseg - 1];
			leftx = allx[targetseg];
		}
	}

	double newphi = ((eqpos - rightx) / (leftx - rightx))*(leftphi - rightphi) + rightphi;

	return(newphi);
}

double Get_Spot_Sigma(const std::vector<double> &Segments, const std::vector<double> &Sigma, double spot) // Get sigma at a given spot; for plotting sigma as a function of the genomic location
{
	int segmentindex = Segments.size() - 1, targetseg = 0, flag = 0;
	double loc = clamp0 + Segments[segmentindex], oldloc = clamp0;
	while(segmentindex >= 0)
	{
		if(spot > oldloc && spot < loc)
		{
			targetseg = segmentindex;
			flag = 1;
			break;
		}
		oldloc = loc;
		segmentindex -= 1;
		loc += Segments[segmentindex];
	}
	if(flag == 0)
	{
		targetseg = 0;
	}

	return(Sigma[targetseg]);
}

std::vector<double> Get_Barrier_Array(double clamp0, double clamp1, int barr_count) // Get a list of barrier protein binding sites: fixed number of binding sites
{
	std::vector<double> Barr;
	for(int i = barr_count - 1; i >= 0; i--)
	{
		Barr.push_back(clamp0 + ((clamp1 - clamp0) / (double (barr_count + 1)))*(double (i + 1)));
	}

	return(Barr);
}

std::vector<double> Get_Barrier_Array_Monica(double clamp0, double clamp1, double guo_gap) // Get a list of barrier protein binding sites: fixed gap between sites
{
	std::vector<double> Barr;
	double start = clamp1 - guo_gap;
	while(start >= clamp0 + guo_gap)
	{
		Barr.push_back(start);
		start -= guo_gap;
	}

	return(Barr);
}

double Marenduzzo_Func(double sigma)
{
	return(myong_k_0 + myong_s_0*(1.0 - tanh(60.0*(sigma + 0.04))));
//	return(0.01 + 0.065*(1.0 - tanh(60.0*(sigma + 0.04))));
//	return((1.0 - std::tanh(60.0*(sigma + 0.04))));
}

double GQ_Dependence(double sigma)
{
//	return(0.3*(1.0 - std::tanh(50.0*(sigma + 0.04))) + 0.4);
	return(1.0);
}

struct observer // To store the time and the sum of rates of events (a0) at each point
{
	std::vector<double> T, P;
	observer(const std::vector<double> param0, const std::vector<double> param1)
	{
		for(int i = 0; i < param0.size(); i++)
		{
			T.push_back(param0[i]);
		}
		for(int i = 0; i < param1.size(); i++)
		{
			P.push_back(param1[i]);
		}
	}
	void operator()(const std::vector<double> &phi_x, double t)
	{
		T.push_back(t);
		P.push_back(phi_x[phi_x.size() - 1]);
	}
};

struct odesystem
{
	int clamp0_flag, clamp1_flag;
	double force;
	std::vector<double> TSSes, gene_lengths;
	std::vector<int> Direction;
	std::vector<double> Exit_Times, Exit_Points;
	std::vector<int> Nucl_Status;
	std::vector<int> current_rna_count;
	int iterid, localiterid, inipolIIcount;
	std::vector<double> Barr, Barr_Phi;
	std::vector<int> Barr_Status;
	std::vector<double> k_on;
	double k_off, k_topo, k_nucl_on, k_nucl_off, k_rna_degrad, k_prot_barrier_on, k_prot_barrier_off, alpha_marenduzzo;
	std::vector<double> Nucl;
	std::vector<int> Gene_ID;
	int brute_force_flag = 0;
	int torque_flag = 0;
	double GapR_On, GapR_Off;
	std::vector<double> GapR_Points;
	std::vector<int> GapR_State;
	int galactose_switch, brute_gapr_flag;
	double k_on_G, k_on_eR, k_off_R, k_off_G, k_off_eR;
	int s_R, s_G, s_eR, s_PQS;

	InterpMultilinear<2, double> interp_ML = Setup_Interp();
	InterpMultilinear<1, double> interp_ML_Cutoff_0 = Setup_Interp_Cutoffs("../torque_interp/sigma_s.log");
	InterpMultilinear<1, double> interp_ML_Cutoff_1 = Setup_Interp_Cutoffs("../torque_interp/sigma_p.log");

	odesystem(int par_clamp0_flag, int par_clamp1_flag, double par_force, const std::vector<double> &par_TSSes, const std::vector<double> &par_gene_lengths, const std::vector<int> &par_Direction, const std::vector<double> &par_Exit_Times, const std::vector<double> &par_Exit_Points, const std::vector<int> &par_Nucl_Status, std::vector<double> &par_k_on, double par_k_off, double par_k_topo, double par_k_nucl_on, double par_k_nucl_off, std::vector<double> &par_Nucl, double par_k_rna_degrad, double par_k_prot_barrier_on, double par_k_prot_barrier_off, const std::vector<double> &par_Barr, const std::vector<int> &par_Barr_Status, const std::vector<double> &par_Barr_Phi, int par_iterid, double par_alpha_marenduzzo, const std::vector<int> &par_Gene_ID, int par_brute_force_flag, int par_torque_flag, const std::vector<int> &par_current_rna_count, double par_GapR_On, double par_GapR_Off, const std::vector<double> &par_GapR_Points, const std::vector<int> &par_GapR_State, int par_galactose_switch, int par_brute_gapr_flag, double par_k_on_G, double par_k_on_eR, double par_k_off_R, double par_k_off_G, double par_k_off_eR, int par_s_R, int par_s_G, int par_s_eR, int par_s_PQS)
	{
		clamp0_flag = par_clamp0_flag; // Fixed
		clamp1_flag = par_clamp1_flag; // Fixed
		force = par_force; // Fixed
		for(int i = 0; i < par_TSSes.size(); i++) // Fixed
		{
			TSSes.push_back(par_TSSes[i]);
		}
		for(int i = 0; i < par_gene_lengths.size(); i++) // Fixed
		{
			gene_lengths.push_back(par_gene_lengths[i]);
		}
		for(int i = 0; i < par_Direction.size(); i++) // Updated every time
		{
			Direction.push_back(par_Direction[i]);
		}
		for(int i = 0; i < par_Exit_Times.size(); i++) // Updated every time
		{
			Exit_Times.push_back(par_Exit_Times[i]);
		}
		for(int i = 0; i < par_Exit_Points.size(); i++) // Updated every time
		{
			Exit_Points.push_back(par_Exit_Points[i]);
		}
		for(int i = 0; i < par_Nucl_Status.size(); i++) // Size fixed, value updated every time
		{
			Nucl_Status.push_back(par_Nucl_Status[i]);
		}
		for(int i = 0; i < par_k_on.size(); i++) // Fixed
		{
			k_on.push_back(par_k_on[i]);
		}
		k_off = par_k_off; // Fixed
		k_topo = par_k_topo; // Fixed
		k_nucl_on = par_k_nucl_on; // Fixed
		k_nucl_off = par_k_nucl_off; // Fixed
		for(int i = 0; i < par_Nucl.size(); i++) // Fixed
		{
			Nucl.push_back(par_Nucl[i]);
		}
		k_rna_degrad = par_k_rna_degrad; // Fixed
		k_prot_barrier_on = par_k_prot_barrier_on; // Fixed
		k_prot_barrier_off = par_k_prot_barrier_off; // Fixed
		for(int i = 0; i < par_Barr.size(); i++) // Fixed
		{
			Barr.push_back(par_Barr[i]);
		}
		for(int i = 0; i < par_Barr_Status.size(); i++) // Size fixed, value updated every time
		{
			Barr_Status.push_back(par_Barr_Status[i]);
		}
		for(int i = 0; i < par_Barr_Phi.size(); i++) // Size fixed, value updated every time
		{
			Barr_Phi.push_back(par_Barr_Phi[i]);
		}
		iterid = par_iterid;
		alpha_marenduzzo = par_alpha_marenduzzo; // Fixed
		for(int i = 0; i < par_Gene_ID.size(); i++) // Updated every time
		{
			Gene_ID.push_back(par_Gene_ID[i]);
		}
		brute_force_flag = par_brute_force_flag; // Fixed
		torque_flag = par_torque_flag; // Fixed
		for(int i = 0; i < par_current_rna_count.size(); i++) // Size fixed, value updated every time
		{
			current_rna_count.push_back(par_current_rna_count[i]);
		}
		GapR_On = par_GapR_On; // Fixed
		GapR_Off = par_GapR_Off; // Fixed
		for(int i = 0; i < par_GapR_Points.size(); i++) // Fixed
		{
			GapR_Points.push_back(par_GapR_Points[i]);
		}
		for(int i = 0; i < par_GapR_State.size(); i++) // Size fixed, value updated every time
		{
			GapR_State.push_back(par_GapR_State[i]);
		}
		galactose_switch = par_galactose_switch; // Updated every time
		brute_gapr_flag = par_brute_gapr_flag;
		k_on_G = par_k_on_G; // Fixed
		k_on_eR = par_k_on_eR; // Fixed
		k_off_R = par_k_off_R; // Fixed
		k_off_G = par_k_off_G; // Fixed
		k_off_eR = par_k_off_eR; // Fixed
		s_R = par_s_R; // Updated every time
		s_G = par_s_G; // Updated every time
		s_eR = par_s_eR; // Updated every time
		s_PQS = par_s_PQS; // Fixed
	}

	void operator()(const std::vector<double> &phi_x, std::vector<double> &dphi_xdt, double t)
	{
		localiterid += 1;
		int PolII_Count = (phi_x.size() - 1) / 2, segmentindex = 0, alivecount = 0;
		if(PolII_Count != inipolIIcount)
		{
			std::cout << "Mismatch between PolII count and ini_PolII count" << "\n";
			throw std::invalid_argument("Mismatch between PolII count and ini_PolII count");
		}
		std::vector<int> Alive = Is_Alive(phi_x, TSSes, gene_lengths, Gene_ID);
		if(Alive.size() != Exit_Times.size())
		{
			std::cout << "Current alive Pol II count does not match the size of the exit times array" << "\n";
			throw std::invalid_argument("Current alive Pol II count does not match the size of the exit times array");
		}
		for(int i = 0; i < Alive.size(); i++)
		{
			if(Alive[i] == 1)
			{
				alivecount += 1;
			}
			if((Alive[i] == 1 && Exit_Times[i] > 0.0))
			{
				std::cout << "Alive Pol II has a positive exit time; this is not possible." << "\n";
				throw std::invalid_argument("Alive Pol II has a positive exit time; this is not possible.");
			}
			if(Exit_Times[i] < 0.0 && Alive[i] == 0)
			{
				Exit_Times[i] = t;
				Exit_Points[i] = phi_x[i + PolII_Count];
				current_rna_count[Gene_ID[i]] += 1;
			}
			if(Alive[i] == 0 && Exit_Times[i] < 0.0)
			{
				std::cout << "Dead Pol II has an invalid exit time." << "\n";
				throw std::invalid_argument("Dead Pol II has an invalid exit time.");
			}
		}

		std::vector<int> PolII_ID, Barr_ID;
		for(int i = 0; i < PolII_Count; i++)
		{
			PolII_ID.push_back(1);
		}
		for(int i = 0; i < Barr.size(); i++)
		{
			Barr_ID.push_back(2);
		}
		std::vector<int> allid = Get_All_Of_Two_Vectors_Ordered<int>(PolII_ID, Barr_ID, 0, 0, phi_x, Alive, Barr, Barr_Status, 5);

		std::vector<double> Segments = Get_Segments_With_Barriers(phi_x, Alive, Barr, Barr_Status);
		std::vector<double> Sigma = Get_Supercoiling_Densities_With_Barriers(clamp0_flag, clamp1_flag, phi_x, Segments, Alive, Barr, Barr_Status, Barr_Phi);
		std::vector<double> allx = Get_All_Of_Two_Vectors_Ordered<double>(phi_x, Barr, PolII_Count, 0, phi_x, Alive, Barr, Barr_Status, 6);
		std::vector<int> nucl_seg_counts = Get_Segment_Nucleosome_Count_Pairwise(allx, Nucl, Nucl_Status);
		std::vector<double> Torques;
		if(torque_flag == 0) // Prokaryotic torque
		{
			Torques = Get_Prokaryotic_Torques(Sigma, Segments, force);
		}
		else if(torque_flag == 1) // Eukaryotic torque
		{
			Torques = Get_Interpolated_Eukaryotic_Torques(interp_ML, interp_ML_Cutoff_0, interp_ML_Cutoff_1, Sigma, Segments, nucl_seg_counts, force);
		}
		else
		{
			std::cout << "Invalid torque flag" << "\n";
			throw std::invalid_argument("Invalid torque flag");
		}
		std::vector<int> Is_RNAP_Blocked = Find_RNAP_Blocked_By_Nucl(Nucl, Nucl_Status, phi_x, Alive, Direction, delta, brute_force_flag);
		std::vector<double> Velocities = Get_Velocities_With_Barriers(phi_x, Segments, Torques, Direction, Alive, Barr, Barr_Status, Is_RNAP_Blocked, delta);
		double denom;
		while(allid.size() > 0 && allid[segmentindex] == 2)
		{
			segmentindex += 1;
		}
		for(int i = 0; i < PolII_Count; i++)
		{
			if(allid.size() > 0 && (Alive[i] == 1 && allid[segmentindex] != 1))
			{
				std::cout << "Plausible error in the allid array assembly" << "\n";
				throw std::invalid_argument("Plausible error in the allid array assembly");
			}
			denom = chi + eta*std::pow(std::abs(phi_x[i + PolII_Count] - TSSes[Gene_ID[i]]), alpha);
			dphi_xdt[i] = 0.0;
			dphi_xdt[i + PolII_Count] = 0.0;
			if(std::abs(Velocities[i]) > 0.0)
			{
				dphi_xdt[i] = ((w0*Velocities[i])*(eta*std::pow(std::abs(phi_x[i + PolII_Count] - TSSes[Gene_ID[i]]), alpha)) / denom) - ((Torques[segmentindex] - Torques[segmentindex + 1]) / denom);
			}
			dphi_xdt[i + PolII_Count] = Velocities[i];
			if(Direction[i] == 1 && Segments[segmentindex] < delta && Velocities[i] > 0.0 && Alive[i] == 1)
			{
				std::cout << "Pol II is moving in the wrong direction or despite steric hindrance." << "\n";
				throw std::invalid_argument("Pol II is moving in the wrong direction or despite steric hindrance.");
			}
			if(Direction[i] == -1 && Segments[segmentindex + 1] < delta && std::abs(Velocities[i]) > 0.0 && Alive[i] == 1)
			{
				std::cout << "Pol II is again moving in rhw wrong direction or despite steric hindrance." << "\n";
				throw std::invalid_argument("Pol II is again moving in rhw wrong direction or despite steric hindrance.");
			}
			if(Alive[i] == 1)
			{
				segmentindex += 1;
			}
			while(allid.size() > 0 && allid[segmentindex] == 2)
			{
				segmentindex += 1;
			}
		}
		double insprop = 0.0;
		std::vector<double> TSS_spotsigma;
		std::vector<std::vector<int>> map, map_gapr;
		std::vector<int> Is_Nucl_Blocked, Is_GapR_Blocked;
		for(int i = 0; i < TSSes.size(); i++)
		{
			TSS_spotsigma.push_back(Get_Spot_Sigma(Segments, Sigma, TSSes[i]));
		}
		double GQ_spotsigma = Get_Spot_Sigma(Segments, Sigma, TSSes[0] + GQ_spacer);
		for(int i = 0; i < TSSes.size(); i++)
		{
//			insprop += k_on[i]*std::max(1.0 - alpha_marenduzzo*TSS_spotsigma[i], 0.0); // RNAP recruitment to TSSes
			insprop += k_on[i]*Marenduzzo_Func(TSS_spotsigma[i]); // RNAP recruitment to TSSes
		}
		insprop += k_off*alivecount; // RNAP drop-off before reaching the end of the gene body
		insprop += k_topo; // Supercoiling relaxation
		Get_Nucl_RNAP_Map_Pairwise(Nucl, allx, map, delta, brute_force_flag);
		Is_Nucl_Blocked = Find_Blocked_Nucl(map);
		Get_GapR_RNAP_Map_Pairwise(GapR_Points, allx, map_gapr, delta);
		Is_GapR_Blocked = Find_Blocked_GapR(map_gapr, brute_gapr_flag);
		for(int i = 0; i < Nucl_Status.size(); i++) // Binding / unbinding of nucleosomes
		{
			if(Nucl_Status[i] == 0 && Is_Nucl_Blocked[i] == 0)
			{
				insprop += k_nucl_on;
			}
			else if(Nucl_Status[i] == 1)
			{
				insprop += k_nucl_off;
			}
		}
		for(int i = 0; i < TSSes.size(); i++)
		{
			insprop += current_rna_count[i]*k_rna_degrad; // RNA degradation
		}
		for(int i = 0; i < Barr_Status.size(); i++) // Barrier protein binding / unbinding
		{
			if(Barr_Status[i] == 0 && Is_Barr_Site_Blocked(Barr[i], phi_x, Alive, delta) == 1)
			{
				insprop += Get_Barrier_On_Rate(k_prot_barrier_on, Get_Spot_Sigma(Segments, Sigma, Barr[i]));
			}
			else if(Barr_Status[i] == 1)
			{
				insprop += Get_Barrier_Off_Rate(k_prot_barrier_off, Get_Spot_Sigma(Segments, Sigma, Barr[i]));
			}
		}
		for(int i = 0; i < GapR_State.size(); i++)
		{
			if(GapR_State[i] == 0 && Is_GapR_Blocked[i] == 0)
			{
				insprop += Get_GapR_On_Rate(GapR_On, Get_Spot_Sigma(Segments, Sigma, GapR_Points[i]), galactose_switch);
			}
			else if(GapR_State[i] == 0 && Is_GapR_Blocked[i] == 1)
			{
				;
			}
			else if(GapR_State[i] == 1)
			{
				insprop += GapR_Off;
			}
			else
			{
				std::cout << "Invalid GapR state / condition" << "\n";
				throw std::invalid_argument("Invalid GapR state / condition");
			}
		}
		if(s_R == 0)
		{
			insprop += k_on[0]*Marenduzzo_Func(TSS_spotsigma[0])*(1.0 - s_R); // R-loop formation
		}
		else if(s_R == 1)
		{
			insprop += (1.0 - s_eR)*s_R*(k_off_R*(1.0 - s_G) + s_G*k_off_eR); // R-loop dissolution
		}
		else
		{
			std::cout << "Invalid R-loop state" << "\n";
			throw std::invalid_argument("Invalid R-loop state");
		}
		if(s_G == 0)
		{
			insprop += s_PQS*s_R*(1.0 - s_G)*k_on_G; // GQ formation
		}
		else if(s_G == 1)
		{
			insprop += k_off_G; // GQ dissolution
		}
		else
		{
			std::cout << "Invalid GQ state" << "\n";
			throw std::invalid_argument("Invalid GQ state");
		}
		if(s_eR == 0)
		{
			insprop += s_R*(1.0 - s_eR)*k_on_eR; // Extended R-loop formation
		}
		else if(s_eR == 1)
		{
			insprop += k_off_eR; // Extended R-loop dissolution
		}
		else
		{
			std::cout << "Invalid eR state" << "\n";
			throw std::invalid_argument("Invalid eR state");
		}

		dphi_xdt[2*PolII_Count] = insprop;
	}
};

struct wrappedobserver
{
	observer *Obs;
	wrappedobserver(observer *param0)
	{
		Obs = param0;
	}
	void operator()(const std::vector<double> &phi_x, double t)
	{
		(*Obs)(phi_x, t);
	}
};

struct wrappedsystem
{
	odesystem *System;
	wrappedsystem(odesystem *param0)
	{
		System = param0;
	}
	void operator()(const std::vector<double> &phi_x, std::vector<double> &dphi_xdt, double t)
	{
		(*System)(phi_x, dphi_xdt, t);
	}
};
