/*yet to implement:
extension of primers, extension by dimers and monomers and how it affects the rate
how error rate is affected by different types of growth

*/


float rr;
float R; //this will be the total rate
float P_on, P_off, P_ext;
int rand_polymer;
int rand_nt;
float dgpair;

float Av = 6.02e23;
float monomer_conc = 1.0e-3;
float volume = 1.0e-15;
float lam = 0.2;
float num_monomers = monomer_conc*Av*volume;
int max_oligo_len = 10;

int initlen = 30;
float growtime = 0.2;
float t = 0.0;

//rate parameters
float RT = 0.616;
float misstack1 = 0.0;
float misstack2 = 1.0;
float dginit = 4.09;
float dgpartial = 1.5;
float kpri = 1.0e6;
float w = 1.0e5;

float single_mis_eqn = w*monomer_conc*exp(-(dgpartial + misstack1)/RT);
float double_mis_eqn = w*monomer_conc*exp(-(dgpartial + misstack2)/RT);

//for storing monomer and oligomer concentration values
std::vector <float> totconc, meanconc, totnum, meannum;


int primer_len = 4;


std::map<std::string, float> stacking = {
	{"AAUU", -0.93},
	{"UUAA", -0.93},
	{"AUUA", -1.10},
	{"UAAU", -1.33},
	{"CUGA", -2.08},
	{"AGUC", -2.08},
	{"CAGU", -2.11},
	{"UGAC", -2.11},
	{"GUCA", -2.24},
	{"ACUG", -2.24},
	{"GACU", -2.35},
	{"UCAG", -2.35},
	{"CGGC", -2.36},
	{"GGCC", -3.26},
	{"CCGG", -3.26},
	{"GCCG", -3.42},
	{"GUUG", 0.72},
	{"GGUU", -0.25},
	{"UUGG", -0.25},
	{"AGUU", -0.35},
	{"UUGA", -0.35},
	{"UGAU", -0.39},
	{"UAGU", -0.39},
	{"UUAG", -0.51},
	{"GAUU", -0.51},
	{"UGGU", -0.57},
	{"AUUG", -0.90},
	{"GUUA", -0.90},
	{"CGGU", -1.25},
	{"UGGC", -1.25},
	{"CUGG", -1.77},
	{"GGUC", -1.77},
	{"GGCU", -1.80},
	{"UCGG", -1.80},
	{"GUCG", -2.15},
	{"GCUG", -2.15}

};
