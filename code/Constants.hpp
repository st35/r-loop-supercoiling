// Fixed model constants

const double w0 = 1.85; // 1.85 nm-1
const double v0 = 20.0; // 20 nm/s
const double tau_c = 12.0; // 12 pN.nm
const double eta = 0.0005; // 5e-4 pN.s
const double alpha = 1.5;
const double chi = 0.05; // 5e-2 pN.nm.s
const double delta = 15.0; // 15 nm
const double nucleosome = 0.34*147.0; // 147 base pairs
const double linker = 0.34*30.0; // 30 base pairs
const double nalpha = 10000.0; // Marenduzzo parameter; not currently used

// Global variables that may be changed

int finite_size_flag = 1;
int circular_plasmid_flag = 1;
double GQ_spacer = 30.0*0.34; // 30 bp from TSS
double sigma_basal = 0.0;
double myong_k_0 = 0.0, myong_s_0 = 0.0; // Inputs for the k_on function
double myong_k_on_G = 1.0, myong_k_on_eR = 1.0, myong_k_off_eR = 0.0;

double clamp0 = -1.0, clamp1 = -1.0;
double clamp0_sigma = 0.0, clamp1_sigma = 0.0;
int approx_count = 0;
double max_approx_diff = -1.0;
