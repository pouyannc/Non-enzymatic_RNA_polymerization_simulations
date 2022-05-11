#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <random>
#include <map>
#include "nonE_RNA_poly_header.h" 
#include <ctime>



int totnum_polymers = 1000;
//int max_len = 100;

//random number between 0-1 using time-based seed
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);
std::uniform_real_distribution<double> distribution (0.0,1.0);


char randomMonomer() {
	switch (rand() % 4 +1) {
		case 1: 
			return 'A';
		case 2: 
			return 'G';
		case 3: 
			return 'C';
		case 4: 
			return 'U';
		default:
			break;

		}
}

std::vector <char> randomSeq(int seq_len) {
	std::vector <char> seq;
	for (int i = 0; i<seq_len; ++i){
		seq.push_back(randomMonomer());
	}
	return seq;
}

class Polymer { //creating a user defined type to take incoming polymers
public:
	std::vector <char> tseq, bseq;
	bool full;
	
	Polymer(std::vector <char> strand)
		:tseq(strand), bseq(), full(false) {}
};
std::vector <Polymer> all_polymers;

bool checkMatch (char btop, char bbot) { //given two RNA bases, check if they are a match
	if (btop == 'A' && bbot == 'U') return true;
	else if (btop == 'U' && bbot == 'A') return true;
	else if (btop == 'G' && bbot == 'C') return true;
	else if (btop == 'C' && bbot == 'G') return true;
	else if (btop == 'G' && bbot == 'U') return true;
	else if (btop == 'U' && bbot == 'G') return true;
	else return false;
}

bool checkTemplate (Polymer p) { //check if the polymer is an empty template or not. TRUE means the polymer is a single strand!
	if (p.bseq.size() == 0) return true;
	else return false;
}

char matchingBase(char nt) { //return a matching base
	if (nt == 'A') return 'U';
	else if (nt == 'U') return 'A';
	else if (nt == 'G') return 'C';
	else if (nt == 'C') return 'G';
}

int duplexStartIndex (Polymer d) { //find the start index where the double strand starts
	for (int i = 0; i < d.tseq.size(); ++i) {
		if (d.bseq[i] != '-') return i;
	}
}

int init_seq_length = 30;
void initPolymers () { //initialize set number of 50-mers 
	for (int i = 0; i < totnum_polymers; ++i) {
		all_polymers.push_back(Polymer(randomSeq(initlen)));
	}
}

void initOligomerConc() { //initialize the concentrations for monomers and small oligomers
	totconc.push_back(0);
	meanconc.push_back(0);
	totnum.push_back(0);
	meannum.push_back(0);
	for (int i = 1; i <= max_oligo_len; ++i) {
		totconc.push_back(4*monomer_conc*pow(lam,i-1));
		meanconc.push_back(totconc[i]/pow(4,i)); 
		totnum.push_back(totconc[i]*Av*volume);
		meannum.push_back(meanconc[i]*Av*volume);
	}
}

float getDGStack(Polymer d, int index, char br_base = '-') { 
	if (br_base == '-') {
		br_base = d.bseq[index+1];
	}
	std::string stack;
	stack.push_back(d.tseq[index]);
	stack.push_back(d.tseq[index+1]);
	stack.push_back(d.bseq[index]);
	stack.push_back(br_base);
	return stacking[stack];
	//return stacking[stack.append(d.tseq[index])];
	//,(d.tseq[index+1]),(d.bseq[index]),(d.bseq[index+1])
}

void finalBPMatch (Polymer b, int &BP_match, int &i_BP) { //iterate and check when the index gets to the final base PAIR, 1 == true, 0 ==  false, 2== null
	for (int i = 1; i<b.bseq.size(); ++i) {
		if (b.bseq[i-1] != '-' && b.bseq[i] == '-') {
			i_BP = (i-1);
			if (checkMatch (b.tseq[i-1], b.bseq[i-1])){ // if final bp is match, return 1
				BP_match = 1;
				return;
			}
			else { // if final bp isnt a match
				BP_match = 0;
				return;
			}
		}
	}
	BP_match = 2;
	i_BP = (-1);
	return;
}

float getOnRate() {
	float tot_on_rate = 0;
	int c = 0;
	for (int i = 0; i < all_polymers.size(); ++i){ //count number of templates without primer bound
		if (checkTemplate(all_polymers[i]) == true) ++c;
	}
	//number of available primers x primer on rate x concentration of specific primer * number of possible binding positions
	tot_on_rate = c*kpri*meanconc[primer_len]*(initlen-primer_len);
	return tot_on_rate;
}

float getOffRate() {
	float tot_off_rate = 0;
	int c = 0;
	dgpair = 0;
	for (int i = 0; i < all_polymers.size(); ++i){ //count number of templates with primer bound
		if (checkTemplate(all_polymers[i]) == false) {
			int duplex_i = duplexStartIndex(all_polymers[i]);
			while (all_polymers[i].bseq[duplex_i+1] != '-') {
				dgpair += getDGStack(all_polymers[i], duplex_i);
				++duplex_i;
			}
			tot_off_rate += kpri*exp((dginit + dgpair)/RT);
		}
	}
	return tot_off_rate;
}

float getExtensionRate(Polymer p, char base_request = '-') { //returns total extension rate for all base possibilities, OR returns the extension rate for base_request base ONLY (if specified)
	int final_BP_bool, i_finalBP;
	finalBPMatch(p, final_BP_bool, i_finalBP);
	char open_base = p.tseq[i_finalBP +1];

	if (checkTemplate(p) == false) { //check if there is a primer bound

		if (final_BP_bool == 0) { //if final bp is a mismatch
			if (open_base == 'A') {
				if (base_request == '-') return (3*(double_mis_eqn) + single_mis_eqn);
				else if (base_request == 'A' || base_request == 'G' || base_request == 'C') return double_mis_eqn;
				else if (base_request == 'U') return single_mis_eqn;
			}
			else if (open_base == 'C') {
				if (base_request == '-') return (3*(double_mis_eqn) + single_mis_eqn);
				else if (base_request == 'A' || base_request == 'U' || base_request == 'C') return double_mis_eqn;
				else if (base_request == 'G') return single_mis_eqn;
			}
			else if (open_base == 'U') {
				if (base_request == '-') return (2*(double_mis_eqn) + 2*(single_mis_eqn));
				else if (base_request == 'C' || base_request == 'U') return double_mis_eqn;
				else if (base_request == 'G' || base_request == 'A') return single_mis_eqn;
			}
			else if (open_base == 'G') {
				if (base_request == '-') return (2*(double_mis_eqn) + 2*(single_mis_eqn));
				else if (base_request == 'A' || base_request == 'G') return double_mis_eqn;
				else if (base_request == 'C' || base_request == 'U') return single_mis_eqn;
			}
		}

		else if (final_BP_bool == 1) { //if final bp is a match
			if (open_base == 'A') {
				if (base_request == '-') return (3*(single_mis_eqn) + (w*monomer_conc*exp(-(dgpartial + getDGStack(p, i_finalBP, 'U'))/RT)));
				else if (base_request == 'A' || base_request == 'G' || base_request == 'C') return single_mis_eqn;
				else if (base_request == 'U') return (w*monomer_conc*exp(-(dgpartial + getDGStack(p, i_finalBP, 'U'))/RT));
			}
			else if (open_base == 'C') {
				if (base_request == '-') return (3*(single_mis_eqn) + (w*monomer_conc*exp(-(dgpartial + getDGStack(p, i_finalBP, 'G'))/RT)));
				else if (base_request == 'A' || base_request == 'U' || base_request == 'C') return single_mis_eqn;
				else if (base_request == 'G') return (w*monomer_conc*exp(-(dgpartial + getDGStack(p, i_finalBP, 'G'))/RT));
			}
			else if (open_base == 'U') {
				if (base_request == '-') return (2*(single_mis_eqn) + (w*monomer_conc*exp(-(dgpartial + getDGStack(p, i_finalBP, 'A'))/RT)) + (w*monomer_conc*exp(-(dgpartial + getDGStack(p, i_finalBP, 'G'))/RT)));
				else if (base_request == 'C' || base_request == 'U') return single_mis_eqn;
				else if (base_request == 'G') return (w*monomer_conc*exp(-(dgpartial + getDGStack(p, i_finalBP, 'G'))/RT));
				else if (base_request == 'A') return (w*monomer_conc*exp(-(dgpartial + getDGStack(p, i_finalBP, 'A'))/RT));
			}
			else if (open_base == 'G') {
				if (base_request == '-') return (2*(single_mis_eqn) + (w*monomer_conc*exp(-(dgpartial + getDGStack(p, i_finalBP, 'C'))/RT)) + (w*monomer_conc*exp(-(dgpartial + getDGStack(p, i_finalBP, 'U'))/RT)));
				else if (base_request == 'G' || base_request == 'A') return single_mis_eqn;
				else if (base_request == 'U') return (w*monomer_conc*exp(-(dgpartial + getDGStack(p, i_finalBP, 'U'))/RT));
				else if (base_request == 'C') return (w*monomer_conc*exp(-(dgpartial + getDGStack(p, i_finalBP, 'C'))/RT));
			}
		}

		else {
			p.full = true;
			return 0; //else, the strand cannot extend further.
		}
	}

	else return 0; //else, the polymer is single stranded.
	
}

float getTotExtensionRate() {
	float tot_ext_rate = 0;

	for (int i  = 0; i < all_polymers.size(); ++i) {
		tot_ext_rate += getExtensionRate(all_polymers[i]);

	}
	return tot_ext_rate;
}

char getExtensionBase (Polymer p) { //return a random base to add.
	float random_base = distribution(generator)*getExtensionRate(p);
	float c_rate, g_rate, a_rate, u_rate;
	c_rate = getExtensionRate (p, 'C');
	g_rate = getExtensionRate (p, 'G');
	a_rate = getExtensionRate (p, 'A');
	u_rate = getExtensionRate (p, 'U');
	
	if (random_base < c_rate) return 'C';
	else if (random_base < c_rate + g_rate) return 'C';
	else if (random_base < c_rate + g_rate + a_rate) return 'A';
	else if (random_base < c_rate + g_rate + a_rate + u_rate) return 'U';

}

void printPolymers() { //output the list of Polymers
	int num_mismatches = 0;
	int num_bound_template = 0;

	for (int i = 0; i < all_polymers.size(); ++i) {
		for (int j = 0; j < all_polymers[i].tseq.size(); ++j){
			std::cout << all_polymers[i].tseq[j];
		}
		std::cout << "\n";
		for (int j = 0; j < all_polymers[i].bseq.size(); ++j){
			std::cout << all_polymers[i].bseq[j];
			if (all_polymers[i].bseq[j] != '-') {
				if (checkMatch(all_polymers[i].tseq[j], all_polymers[i].bseq[j]) == false) ++num_mismatches;
			}
		}
		std::cout << "\n \n";

	}
	std::cout << all_polymers.size() << "\n";
	std::cout << "final time: " << t << "\n";


	for (int i = 0; i<all_polymers.size(); ++i) {
		if (checkTemplate(all_polymers[i]) == false) ++num_bound_template;
	}
	std::cout << "number of bound polymers: " << num_bound_template<< "\n";
	std::cout << "number of mismatches: " << num_mismatches<< "\n";
	
}

int main() {

	
	initPolymers();
	initOligomerConc();

	std::cout << "t" <<t<<"\n" << "growtime" << growtime << "\n";

	clock_t t_real;
	t_real = clock();

	while ((clock()-t_real) / CLOCKS_PER_SEC < 30) {

		P_on = getOnRate();
		P_off = getOffRate();
		P_ext = getTotExtensionRate();
		

		//std::cout << "tot ON rate: " << P_on << "\ntot OFF rate: " << P_off << "\ntot EXT rate: " << P_ext << "\n";

		R = P_on + P_off + P_ext;
		t += -log(distribution(generator)) / R;	
		rr = distribution(generator)*R;
		std::cout<<"Current time: " << t << "\n";

		rand_polymer = rand() % totnum_polymers;

		//if (t>growtime) t = growtime;

		if (rr< P_on){ //add primer
			while (checkTemplate(all_polymers[rand_polymer]) == false){ //get a polymer that is just template
				rand_polymer = rand() % totnum_polymers;
			}
			rand_nt = rand() % (initlen-primer_len);
			//set all of bottom sequence to 0 for the chosen polymer
			all_polymers[rand_polymer].bseq.resize(initlen,'-');
			for (int i = rand_nt; i < rand_nt+primer_len; ++i) { //bind a matching primer
				all_polymers[rand_polymer].bseq[i] = matchingBase(all_polymers[rand_polymer].tseq[i]);
			}
		}

		else if (rr < (P_off + P_on) ) { //remove primer
			while (checkTemplate(all_polymers[rand_polymer]) == true){ //get a polymer that is duplex
				rand_polymer = rand() % totnum_polymers;
			}
			all_polymers[rand_polymer].bseq.clear();
		}

		else if (rr < (P_off + P_on + P_ext)) { //extend a primer by a monomer
			int final_BP_bool, i_finalBP;
			i_finalBP = (-1);
			while (checkTemplate(all_polymers[rand_polymer]) == true || i_finalBP == (-1)) { //get a polymer that is duplex and can extend
				rand_polymer = rand() % totnum_polymers;
				finalBPMatch(all_polymers[rand_polymer], final_BP_bool, i_finalBP);
			}
			char ext_base = getExtensionBase(all_polymers[rand_polymer]);
			all_polymers[rand_polymer].bseq[i_finalBP +1] = ext_base;
		}
	}

	printPolymers();
}
	

