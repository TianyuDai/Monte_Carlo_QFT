#include <fstream> 
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <random>
std::ofstream fout("block_2d_g1_k300_64-64"); 
std::ofstream fout2("Corr_2d_g1_k300_64-64"); 

const size_t d = 2, LT = 64, LX = 64; 
int dims[d] = {LX, LT}; 
const int nblock = 1000, blocks = 10000, nSign = 1; 
double g = 0.01, kappa = 0.300, sigma = 1; 
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
std::vector < std::vector<int> > visited; 
std::vector<std::vector <int> > nb; 

int Sign_Update()
{
	// Reset the flag grid
	for (size_t i = 0; i < LX; i++)
		flag[i] = std::vector<bool> (LT, false); 
	visited.clear(); 
	visited.reserve(pow(LX, d-1) * LT); 

	// Pick the origin as the first site in cluster
	flag[0][0] = true; 
	chi[0][0] *= -1.; 
	visited[0] = std::vector<int>(d, 0); 
	int count = 1; 
	
	// Going over all visited sites and all neighboring directions
	for (int i=0; i < count; i++)
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
		for (int j = 0; j < 2 * d; j++)
		{
			double beta = kappa * chi[visited[i][0]][visited[i][1]] * chi[nb[j][0]][nb[j][1]]; 
			if (beta < 0 && !flag[nb[j][0]][nb[j][1]])
			{
				double p = 1. - exp(2. * beta); 
				if (((double) rand() / RAND_MAX) < p)
				{
					visited[count].reserve(d); 
					visited[count][0] = nb[j][0]; 
					visited[count++][1] = nb[j][1]; 
					flag[nb[j][0]][nb[j][1]] = true; 
					chi[nb[j][0]][nb[j][1]] *= -1.; 
				}
			}
		}

	}
	return count; 
}

void Regular_Update()
{
	//Regular update
	for (size_t i=0; i < LX; i++)
	{
		for (size_t j = 0; j < LT; j++)
		{
			// All neighboring directions
			int left = (LT + j - 1) % LT; 
			int right = (j + 1) % LT; 
			int up = (LX + i - 1) % LX; 
			int down = (i + 1) % LX; 
			double alpha = kappa * (chi[i][left] + chi[i][right] + chi[up][j] + chi[down][j]); 
			double chi_new, chi_old; 
			std::normal_distribution<double> dist_alpha(alpha, sigma); 
			chi_new = dist_alpha(generator); 
			chi_old = chi[i][j]; 
			double p = exp(-1. * g * pow(chi_new, 4)) / exp(-1. * g * pow(chi_old, 4)); 
			if (((double) rand() / RAND_MAX) < p) chi[i][j] = chi_new; 
		}
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
	for (int i = 0; i < 2 * d; i++)
		nb[i].reserve(d);  
	for (size_t i = 0; i < LX; i++)
	{
		chi[i] = std::vector <double> (LT, 0.5); 
		flag[i].reserve(LT); 
	}

	std::vector <double> corr_block(nblock, 0.); 
	// Calculation
	for (size_t i = 0; i < nblock; i++)
	{
		for (size_t j = 0; j < blocks; j++)
		{
			Regular_Update(); 
			for (int k = 0; k < nSign; k++)
			{
				int c = Sign_Update(); 
				for (int l = 0; l < c; l++)
					if (visited[l][1] == LT/2) corr_block[i] += chi[0][0] * chi[visited[l][0]][visited[l][1]]; 
			}
		}
		corr_block[i] /= blocks * nSign; 
		fout << corr_block[i] <<std::endl; 
		std::cout << i << " Block" << std::endl; 
	}

	for (size_t i = 0; i < LX; i++)
 	{
		for (size_t j = 0; j < LT; j++)
			fout2 << chi[i][j] << " "; 
		fout2 << std::endl; 
	}

	double corr = 0., corr_sq = 0., corr_err; 
	int thermo = 10, n = nblock - thermo; 
	for (size_t i = 0; i < nblock; i++)
	{
		fout << corr_block[i] << std::endl; 
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
