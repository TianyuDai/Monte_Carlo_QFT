#include <fstream> 
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <random>
std::ofstream fout("block_4d_g0_4-16"); 
std::ofstream fout2("Corr_Func_4d_g0_4-16"); 

const size_t d = 4, LT = 16, LX = 4; 
int dims[d] = {LX, LX, LX, LT}; 
const int nblock = 200, blocks = 10000, nSign = 40; 
double g = 0, kappa = 0.124843945, sigma = 1; 
template<size_t dim, typename T>
struct multidim_vector
{
	typedef std::vector<typename multidim_vector<dim-1, T>::type > type; 
}; 
template<typename T>
struct multidim_vector<0,T>
{
	typedef T type; 
};

multidim_vector<d, double>::type chi; 
multidim_vector<d, bool>::type flag; 
std::default_random_engine generator; 
multidim_vector<2, int>::type visited, nb; 

int Sign_Update()
{
	// Reset the flag grid
	for (size_t i = 0; i < LX; i++)
		for (size_t j = 0; j < LX; j++)
			for (size_t k = 0; k < LX; k++)
				flag[i][j][k] = std::vector<bool> (LT, false); 
	visited.clear(); 
	//visited.reserve(pow(LX, d-1) * LT); 

	// Pick the origin as the first site in cluster
	flag[0][0][0][0] = true; 
	chi[0][0][0][0] *= -1.; 
	visited[0] = std::vector<int> (d, 0); 
	int count = 1; 

	// Going over all visited sites and all neighboring directions
	for (int i = 0; i < count; i++)
	{
		// All neighboring directions
		for (size_t j = 0; j < d; j++)
		{
			int before = (dims[j] + visited[i][j] - 1) % dims[j]; 
			int after = (visited[i][j] + 1) % dims[j]; 
			for (size_t k = 0; k < d; k++)
			{
				nb[2*j][k] = visited[i][k]; 
				nb[2*j+1][k] = visited[i][k]; 
			}
			nb[2*j][j] = before; 
			nb[2*j+1][j] = after; 
		}
		for (size_t j = 0; j < 2 * d; j++)
		{
			double beta = kappa * chi[visited[i][0]][visited[i][1]][visited[i][2]][visited[i][3]] * chi[nb[j][0]][nb[j][1]][nb[j][2]][nb[j][3]]; 
			if (beta < 0 && !flag[nb[j][0]][nb[j][1]][nb[j][2]][nb[j][3]])
			{
				double p = 1. - exp(2. * beta); 
				if (((double) rand() / RAND_MAX) < p)
				{
					visited[count].reserve(d); 
					visited[count][0] = nb[j][0]; 
					visited[count][1] = nb[j][1]; 
					visited[count][2] = nb[j][2]; 
					visited[count][3] = nb[j][3]; 
					flag[nb[j][0]][nb[j][1]][nb[j][2]][nb[j][3]] = true; 
					chi[nb[j][0]][nb[j][1]][nb[j][2]][nb[j][3]] *= -1.; 
					count++; 
				}
			}
		}

	}
	return count; 
}

void Regular_Update()
{
	//Regular update
	//double alpha = 0.; 
	std::vector < std::vector <int> > nb_ij; 
	nb_ij.reserve(2*d); 
	for (size_t i = 0; i < LX; i++)
		for (size_t j = 0; j < LX; j++)
			for (size_t k = 0; k < LX; k++)
				for (size_t l = 0; l < LT; l++)
				{
					// All neighboring directions
					std::vector <int> coor{(int)i, (int)j, (int)k, (int)l}; 
					double alpha = 0.; 
					for (size_t m = 0; m < d; m++)
					{
						nb_ij[2*m].reserve(d); 
						nb_ij[2*m+1].reserve(d); 
						int before = (dims[m] + coor[m] - 1) % dims[m]; 
						int after = (coor[m] + 1) % dims[m]; 
						std::vector<int> nb_ij(d); 
						for (size_t n = 0; n < d; n++)
							nb_ij.at(n) = coor.at(n); 
						nb_ij.at(m) = before; 
						alpha += chi[nb_ij[0]][nb_ij[1]][nb_ij[2]][nb_ij[3]]; 
						nb_ij.at(m) = after; 
						alpha += chi[nb_ij[0]][nb_ij[1]][nb_ij[2]][nb_ij[3]]; 
					}
					alpha *= kappa; 
					double chi_new, chi_old; 
					std::normal_distribution<double> dist_alpha(alpha, sigma); 
					chi_new = dist_alpha(generator); 
					chi_old = chi[i][j][k][l]; 
					double p = exp(-1. * g * pow(chi_new, 4)) / exp(-1. * g * pow(chi_old, 4)); 
					if (((double) rand() / RAND_MAX) < p) chi[i][j][k][l] = chi_new; 
				}
}

int main()
{
	fout2 << "dim = " << d << ", LX = " << LX << ", LT = " << LT <<std::endl; 
	fout2 << "g = " << g << ", kappa = " << kappa << std::endl; 
	fout2 << "number of blocks: " << nblock << " number of steps in each block: " << blocks << std::endl; 

	// Initialize the grid
		
	chi.reserve(LX); 
	flag.reserve(LX); 
	nb.reserve(2*d); 
	visited.reserve(pow(LX, d-1) * LT); 
	for (int i = 0; i < 2 * d; i++)
		nb[i].reserve(d); 
	for (size_t i = 0; i < LX; i++)
	{
		chi[i].reserve(LX); 
		flag[i].reserve(LX); 
		for (size_t j = 0; j < LX; j++)
		{
			chi[i][j].reserve(LX); 
			flag[i][j].reserve(LX); 
			for (size_t k = 0; k < LX; k++)
			{
				chi[i][j][k] = std::vector <double> (LT, 0.); 
				flag[i][j][k].reserve(LT); 
			}
		}
	}

	std::vector <double> corr_block(nblock, 0.); 
	// Calculation
	for (int i = 0; i < nblock; i++)
	{
		for (int j = 0; j < blocks; j++)
		{
			Regular_Update(); 
			for (int k = 0; k < nSign; k++)
			{
				int c = Sign_Update(); 
				for (int l = 0; l < c; l++)
					if (visited[l][d-1] == LT/2) corr_block[i] += chi[0][0][0][0] * chi[visited[l][0]][visited[l][1]][visited[l][2]][visited[l][3]]; 
			}
		}
		corr_block[i] /= blocks * nSign; 
		fout << corr_block[i] << std::endl; 
		std::cout << i << " Block" << std::endl; 
	}

	double corr = 0., corr_sq = 0., corr_err; 
	int thermo = 10, n = nblock - thermo; 
	for (size_t i = 0; i < nblock; i++)
	{
		if (i > thermo - 1)
		{
			corr += corr_block[i]; 
			corr_sq += pow(corr_block[i], 2); 
		}
	}
	corr /= n; 
	corr_err = sqrt((corr_sq / n - pow(corr, 2))/n); 
	fout2 << corr << " " << corr_err << std::endl; 
	return 0; 
}
