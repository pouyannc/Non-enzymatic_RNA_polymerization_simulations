#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <string>

void printVec(const std::vector<std::string>& vec)
{
	//print out the strands and their size distribution
	std::vector <int> size_dist(300, 0);

	for (int i=0; i < vec.size(); ++i)
	{
		std::cout << vec[i] << ", ";
		size_dist[vec[i].size()] +=1;
	}
	std::cout << '\n';
	for (int i=0; i < size_dist.size(); ++i)
	{
		std::cout << i << '\t' << size_dist[i] << '\n';
	}
	std::cout << std::endl;

}


int numBonds(const std::vector<std::string>& strands, const std::vector<std::string>& templates)
{
	//calculate total number of bonds in the set of strands
	int num_bonds = 0;
	for (int i = 0; i < strands.size(); ++i)
	{
		const std::string& v = strands[i];
		num_bonds += (v.size() -1);
	}
	for (int i = 0; i < templates.size(); ++i)
	{
		const std::string& v2 = templates[i];
		num_bonds += (v2.size() -1);
	}
	return num_bonds;
}

float getLigationRate(const std::vector<std::string>& strands, const std::vector<std::string>& templates; float s_rate, float volume)
{
	int total_strands = strands.size() + templates.size();
	return s_rate*(total_strands * (total_strands-1)/volume);
}

float getHydrolysisRate(const std::vector<std::string>& strands, float h_rate)
{
	return (h_rate * numBonds(strands));
}

float getAnnealRate(const std::vector<std::string>& strands, const std::vector<std::string>& templates, float a_rate, float volume)
{
	return (a_rate * (strands.size()*templates.size()/volume));
}

int numDuplexes(const std::vector<std::string>& primers)
{
	int n_duplexes = 0;
	for (int i = 0; i < primers.size(); ++i)
	{
		if (primers[i] != "")
		{
			++n_duplexes;
		}
	}
	return n_duplexes;
}

float getMeltingRate(const std::vector<std::string>& primers, float m_rate, float volume)
{
	int num_duplexes = numDuplexes(primers);
	return (m_rate* (num_duplexes/volume));
}



int main()
{
	//initial conditions
	srand(time(NULL));
	std::vector <std::string> strands;
	std::vector <std::string> templates;
	std::vector <std::string> primers;
	float s_rate = 1;
	float h_rate = 1;
	float a_rate = 1;
	float m_rate = 0.01;
	float volume = 10000;
	float num_nt = 10000;
	//populate strands vector with equal parts A, C, G, U
	for (int i=0; i < num_nt; ++i)
	{
		if (i < 0.25*num_nt)
		{
			strands.push_back("A");
		}
		else if (i < 0.5*num_nt)
		{
			strands.push_back("U");
		}
		else if (i < 0.75*num_nt)
		{
			strands.push_back("C");
		}
		else
		{
			strands.push_back("G");
		}
		
	}
	std::cout << strands.size() << std::endl;

	//loop variables
	float t = 0;
	float k_lig;
	float k_hyd;
	float k_ann;
	float k_melt;
	float R;
	float P_lig;
	float P_ann;
	float P_melt;

	while(t<=15)
	{
		k_lig = getLigationRate(strands, s_rate, volume);
		k_hyd = getHydrolysisRate(strands, h_rate);
		k_ann = getAnnealRate(strands, templates, a_rate, volume);
		k_melt = getMeltingRate(primers, m_rate, volume)
		R = k_lig + k_hyd + k_ann + k_melt;
		t += -(log(float(rand())/RAND_MAX)) / R;
		// std::cout << "Current time: " << t << std::endl;
		std::cout << "Current time: " << t << "                \r";

		P_lig = k_lig / R;
		P_ann = k_ann / R;
		P_melt = k_melt / R;

		//ligation condition
		if (float(rand())/RAND_MAX < P_lig)
		{
			//ligation reaction
			//
			//getting random strand from all strands+templates
			int ri_remove = rand() % (strands.size() + templates.size());
			int ri_extend = rand() % (strands.size() + templates.size());

			//removing from small strands:
			if (ri_remove < strands.size())
			{
				std::string extension = strands[ri_remove];
				strands.erase(strands.begin() + (ri_remove));
				//extending to small strands:
				if (ri_extend < strands.size())
				{
					strands[ri_extend] += extension;
					//moving to templates if strand size is 4 or more:
					if strands[ri_extend].size() >= 4
					{
						templates.push_back(strands[ri_extend]);
						primers.push_back("");
						strands.erase(strands.begin() + (ri_extend))
					}
				}
				//extending to template strands:
				else
				{
					templates[ri_extend - strands.size()] += extension;
				}

			}

			//removing from template strands:
			else
			{
				std::string extension = templates[ri_remove - strands.size()];
				if (ri_extend < strands.size())
					{
						templates[ri_remove - strands.size()] += strands[ri_extend];
						strands.erase(strands.begin() + (ri_extend));
					}
					else
					{
						
					}
			}
			
			strands.erase(strands.begin() + (ri_remove));
			
			strands[ri_extend] += extension;
			// std::cout << "Done Ligating" << std::endl;
			//printVec(strands);
		}

		else if (float(rand())/RAND_MAX < (P_lig + P_ann))
		{

		}

		else if (float(rand())/RAND_MAX < (P_lig + P_ann + P_melt))
		{

		}

		//hydrolysis condition
		else
		{
			//hydrolysis reaction
			// std::cout << "Hydrolysing.. " << std::endl;

			int num_bonds = numBonds(strands);
			int bond = rand() % num_bonds + 1;
			int c = 0;
			for (int i = 0; i < strands.size(); ++i)
			{
				if (strands[i].size() >1)
				{
					for (int j = 0; j < strands[i].size() -1; ++j)
					{
						++c;
						if (c == bond)
						{
							strands.push_back(strands[i].substr(j+1));
							strands[i] = strands[i].substr(0,j+1);
						}
					}
				}
			}
			// std::cout << "Done Hydrolysing" <<  std::endl;
			//printVec(strands);


		}


	}


	printVec(strands);
	std::cout << "done" << '\n';
	//int ri = rand();
	//float rd = float(rand())/RAND_MAX;
	return 0;
}