
# include <iostream>
# include <fstream>
# include <sstream>
# include <stdlib.h>
# include <cmath>
# include <math.h>
# include <random>
# include <cstdlib>
# include <dirent.h>
# include <string>
# include <getopt.h>



using namespace std;



/*******************************************************************************/



// Define global variables
const int _maxsize = 1e5;

// Next-nearest-neighbour neighbourhood -> 8 neighbours 
//const int NEIGHBOURHOOD = 8;

int q, seed;

double radius_double, t, r_birth, r_death, r_birth_normalised, r_death_normalised, rand_double, cell_cell_chord_gradient, cell_cell_chord_yIntercept, min_move_distance, selection_coeff;
double optimal_direction_i,  optimal_direction_j, optimal_vector_norm, vector_norm, rescaled_min_length, scalar_prod, dist, nearest_space_distance, x_boundary_intersect, y_boundary_intersect;
double weibull_shape, weibull_scale, birth_rate_sum;
int radius, Ntot, Nwt, iter, x, y, cell_x, cell_y, dir, queue, ind, length, coordX, coordY, previous_link_direction, chain_length;
int arising_time, x_b, y_b, direction, chosen_direction, min_length, num_mins, chosen_min, doubled_ecDNA_copyNumber, daughter_ecDNA_copyNumber1, daughter_ecDNA_copyNumber2;
int empty_cell_x, empty_cell_y, search_radius, num_nearest_empty_cells, rand_int, vertical_direction, horizontal_direction, i_lowerBound, j_lowerBound, i_upperBound, j_upperBound;
int next_move_direction_x, next_move_direction_y, motherCell_copyNumber, clusterSize, initial_copyNumber;



// Define arrays which will contain relative coordinates of empty neighbours for a chosen cell 
int chainX[(int)_maxsize];
int chainY[(int)_maxsize];
int moves_direction_x[(int)_maxsize];
int moves_direction_y[(int)_maxsize];
double moves_distances[(int)_maxsize];
int previous_moves_direction_x[(int)_maxsize];
int previous_moves_direction_y[(int)_maxsize];
double previous_moves_distances[(int)_maxsize];



bool BIRTH = false;
bool DEATH = false;
bool quiet = true;
bool skip = false;
bool clustering = false;
bool weibull_fitness = false;



/*******************************************************************************/



// Define poisson distributions
default_random_engine generator;

// Declare Poisson distribution
poisson_distribution<int> ecDNA_clusterSizeDistribution(2.0);



/*******************************************************************************/



// Define a cell
class Cell
{
	public:
		int ecDNA;

	// Constructor for Cell object
	Cell(){}

	// Set() and get() methods
	
	void set_ecDNA(int n)
	{
		this->ecDNA = n;
	}
};




// Fitness function cells carrying x copies of ecDNA (either flat, or Weibull function with parameters (shape)=1.8; (scale)=40)
double fitness(int x , bool weibull_fitness)
{
	if (weibull_fitness == true)
	{
		weibull_shape = 1.8;
		weibull_scale = 40.0;

		return 50.0 * (weibull_shape/weibull_scale) * pow((x/weibull_scale) , weibull_shape-1) * exp(-pow((x/weibull_scale) , weibull_shape));
	}

	else
	{
		if (x > 0) return 1;
		else return 0;
	}
}




// Straight line pushing 
void straight_line_division(Cell ** tissue , int cell_x , int cell_y , int *Ntot , int *x_b , int *y_b , double *r_birth , int radius)
{

	//cout << "\nCell at (x , y) = (" << cell_x << " , " << cell_y << ") wants to divide." << endl;

	// Find nearest empty lattice point (this becomes the end of the chain of cells)
	empty_cell_x = -1;
	empty_cell_y = -1;
	nearest_space_distance = 1e10;
	num_nearest_empty_cells = 1;
	search_radius = 0;

	do
	{
		search_radius += 1;
		//cout << "Searching cells within radius of " << search_radius << endl;

		if (search_radius > q) return;

		for (int i = -search_radius; i <= search_radius; ++i)
		{
			for (int j = -search_radius; j <= search_radius; ++j)
			{

				//if ((i == 0) && (j == 0)) continue;

				// Compute distance 
				dist = (i*i) + (j*j);

				// Only check cells within distance search_radius of dividing cell
				if (dist > float(search_radius*search_radius)) continue;

				//cout << "Checking cell status at (x,y) = (" << cell_x + i << " , " << cell_y + j << ") -> " << tissue[cell_x + i][cell_y + j].ecDNA << endl;

				// Check if lattice coordinates are occupied or empty
				if (tissue[cell_x + i][cell_y + j].ecDNA == -1)
				{

					//cout << i << " " << j << endl;

					if (dist < nearest_space_distance)
					{
						nearest_space_distance = dist;
						empty_cell_x = cell_x + i;
						empty_cell_y = cell_y + j;

						num_nearest_empty_cells = 1;
					}

					else if (dist == nearest_space_distance) 
					{
						++num_nearest_empty_cells;
					}
				}
			}
		}
	}
	while((empty_cell_x == -1) && (empty_cell_y == -1));









	// If there is more than 1 empty cell at the nearest distance to the dividing cell, choose one with uniform probability
	if (num_nearest_empty_cells > 1)
	{


		// Create vector of co-ordinates to nearest empty lattice points 
		vector<int> nearest_empty_cells_x;
		vector<int> nearest_empty_cells_y;

		nearest_empty_cells_x.clear();
		nearest_empty_cells_y.clear();

		for (int i = -search_radius; i <= search_radius; ++i)
		{
			for (int j = -search_radius; j <= search_radius; ++j)
			{
				// Compute distance 
				dist = (i*i) + (j*j);

				// Only check cells within 20 cell radius of dividing cell
				if (dist > float(search_radius*search_radius)) continue;

				// Check if lattice coordinates are occupied or empty
				if ((tissue[cell_x + i][cell_y + j].ecDNA == -1) && (dist == nearest_space_distance))
				{
					nearest_empty_cells_x.push_back(cell_x + i);
					nearest_empty_cells_y.push_back(cell_y + j);
					//cout << i << " " << j << endl;
					//cout << "Nearest empty cell at (x , y) = (" << i << " , " << j << ")" << endl;
				}
			}
		}



		// Choose which of the cells will divide
		rand_int = int(round(drand48()*(num_nearest_empty_cells-1)));
		empty_cell_x = nearest_empty_cells_x[rand_int];
		empty_cell_y = nearest_empty_cells_y[rand_int];



	}

	//cout << "Cell at (x , y) = (" << cell_x << " , " << cell_y << ") will push cells towards (" << empty_cell_x << " , " << empty_cell_y << ")" << endl;
	//cout << cell_x - empty_cell_x << " " << cell_y - empty_cell_y << endl;





	// Decide whether horizontal cell movements will be left or right, similarly for vertical movements
	horizontal_direction = 0;
	vertical_direction = 0;

	if (empty_cell_x > cell_x) horizontal_direction = 1;
	else if (cell_x > empty_cell_x) horizontal_direction = -1;

	if (empty_cell_y > cell_y) vertical_direction = 1;
	else if (cell_y > empty_cell_y) vertical_direction = -1;

	//cout << "Horizontal & vertical move directions are: " << horizontal_direction << " & " << vertical_direction << endl;








	// Draw straight line between dividing cell and empty lattice point
	cell_cell_chord_gradient = 0.0;
	cell_cell_chord_yIntercept = 0.0;
	if (empty_cell_x != cell_x)
	{
		cell_cell_chord_gradient = float(empty_cell_y - cell_y)/float(empty_cell_x - cell_x);
		cell_cell_chord_yIntercept = cell_y - (cell_cell_chord_gradient * cell_x);
	}

	//cout << "Cell-cell chord defined as y = (" << cell_cell_chord_gradient << ")x + (" << cell_cell_chord_yIntercept << ")" << endl;
	








	// Construct chain by finding intersection of cell-cell chord with vertical & horizontal lines between (cell_x , cell_y) and (empty_cell_x , empty_cell_y)
	i_lowerBound = cell_x;
	i_upperBound = empty_cell_x;
	j_lowerBound = cell_y;
	j_upperBound = empty_cell_y;

	if (cell_x > empty_cell_x)
	{
		i_lowerBound = empty_cell_x;
		i_upperBound = cell_x;
	}

	if (cell_y > empty_cell_y)
	{
		j_lowerBound = empty_cell_y;
		j_upperBound = cell_y;
	}

	queue = 0;






	for (int i = i_lowerBound+1; i < i_upperBound; ++i)
	{
		//cout << "Row/column interface at x = " << i << endl;

		// Find coordinates of intersection between line x = i and cell-cell chord
		y_boundary_intersect = (cell_cell_chord_gradient * float(i)) + cell_cell_chord_yIntercept;

		// Find (squared) Euclidean distance between (i , y_boundary_intersect) and (cell_x , cell_y)
		dist = pow((cell_x - i) , 2) + pow((cell_y - y_boundary_intersect) , 2);

		// Store move in list
		moves_direction_x[queue] = horizontal_direction;
		moves_direction_y[queue] = 0;

		// Store distance in list
		moves_distances[queue] = dist;

		queue += 1;
	}





	for (int j = j_lowerBound+1; j < j_upperBound; ++j)
	{
		//cout << "Row/column interface at y = " << j << endl;

		// Find coordinates of intersection between line x = i and cell-cell chord
		if (empty_cell_x != cell_x) x_boundary_intersect = (float(j) - cell_cell_chord_yIntercept)/cell_cell_chord_gradient;
		else x_boundary_intersect = cell_x;

		//cout << x_boundary_intersect << endl;

		// Find (squared) Euclidean distance between (i , y_boundary_intersect) and (cell_x , cell_y)
		dist = pow((cell_x - x_boundary_intersect) , 2) + pow((cell_y - j) , 2);

		//cout << dist << endl;

		// Store move in list
		moves_direction_x[queue] = 0;
		moves_direction_y[queue] = vertical_direction;

		// Store distance in list
		moves_distances[queue] = dist;

		queue += 1;
	}




	// cout << "List of moves:" << endl;
	// for (int i = 0; i < queue; ++i)
	// {
	// 	cout << moves_direction_x[i] << " , " <<  moves_direction_y[i] << endl;
	// }








	// Reset lists containing previous move information 
	for (int i = 0; i < queue; ++i)
	{
		previous_moves_distances[i] = 0.0;
		previous_moves_direction_x[i] = -10;
		previous_moves_direction_y[i] = -10;
	}









	// Loop through moves and their distances, and construct chain of cell coordinates.
	chain_length = 0;
	while (chain_length != queue)
	{
		min_move_distance = 1e10;

		for (int i = 0; i < queue; ++i)
		{
			skip = false;

			//cout << "Checking if " << moves_distances[i] << " < " << min_move_distance << endl;
			if (moves_distances[i] <= min_move_distance)
			{
				// Check if distance is in previous_moves list
				for (int j = 0; j < chain_length; ++j)
				{
					if ((previous_moves_distances[j] == moves_distances[i]) && (previous_moves_direction_x[j] == moves_direction_x[i]) && (previous_moves_direction_y[j] == moves_direction_y[i])) skip = true;
				}

				if (skip == true) continue;

				min_move_distance = moves_distances[i];
				next_move_direction_x = moves_direction_x[i];
				next_move_direction_y = moves_direction_y[i];
			}
			//cout << "Adding move (" << next_move_direction_x << " , " << next_move_direction_y << ") to chain." << endl;
		}

		if (chain_length == 0)
		{
			chainX[chain_length] = cell_x + next_move_direction_x;
			chainY[chain_length] = cell_y + next_move_direction_y;

			previous_moves_distances[chain_length] = min_move_distance;
			previous_moves_direction_x[chain_length] = next_move_direction_x;
			previous_moves_direction_y[chain_length] = next_move_direction_y;

			chain_length += 1;

			//cout << "1. Adding cell (" << cell_x + next_move_direction_x << " , " << cell_y + next_move_direction_y << ") to chain." << endl;

		}
		else
		{
			chainX[chain_length] = chainX[chain_length - 1] + next_move_direction_x;
			chainY[chain_length] = chainY[chain_length - 1] + next_move_direction_y;

			previous_moves_distances[chain_length] = min_move_distance;
			previous_moves_direction_x[chain_length] = next_move_direction_x;
			previous_moves_direction_y[chain_length] = next_move_direction_y;

			chain_length += 1;

			//cout << "2. Adding cell (" << chainX[chain_length - 1] + next_move_direction_x << " , " << chainY[chain_length - 1] + next_move_direction_y << ") to chain." << endl;

		}

	}


	// cout << "Cell chain:" << endl;
	// for (int i = 0; i < queue; ++i)
	// {
	// 	cout << chainX[i] << " , " <<  chainY[i] << endl;
	// }









	// Once the chain has been constructed, move all cells along one place
	if (queue > 0)
	{

		// Add two final moves (one vertical and one horizontal) to the end of the chain
		if (drand48() < 0.5)
		{
			chainX[queue] = chainX[queue - 1] + horizontal_direction;
			chainY[queue] = chainY[queue - 1] + 0;

			++queue;

			chainX[queue] = chainX[queue - 1] + 0;
			chainY[queue] = chainY[queue - 1] + vertical_direction;
		}
		else
		{
			chainX[queue] = chainX[queue - 1] + 0;
			chainY[queue] = chainY[queue - 1] + vertical_direction;

			++queue;

			chainX[queue] = chainX[queue - 1] + horizontal_direction;
			chainY[queue] = chainY[queue - 1] + 0;
		}


		for (int i = 0; i < queue; ++i)
		{
			tissue[chainX[queue-i]][chainY[queue-i]].ecDNA = tissue[chainX[queue-i-1]][chainY[queue-i-1]].ecDNA;
			//cout << "(" << chainX[queue-i-1] << " , " << chainY[queue-i-1] << ") -> (" << chainX[queue-i] << " , " << chainY[queue-i] << ")" << endl;


			// Update bounds on tissue size
			if (fabs(chainX[i] + 1 - radius) > *x_b) *x_b = fabs(chainX[i] + 1 - radius);
			if (fabs(chainY[i] + 1 - radius) > *y_b) *y_b = fabs(chainY[i] + 1 - radius);
		}

	}

	else 	// Even if queue=0, check that newly created cell increases any bounds
	{
		// Update bounds on tissue size
		if (fabs(chainX[0] + 1 - radius) > *x_b) *x_b = fabs(chainX[0] + 1 - radius);
		if (fabs(chainY[0] + 1 - radius) > *y_b) *y_b = fabs(chainY[0] + 1 - radius);
	}



	if (queue == 0)
	{
		chainX[0] = empty_cell_x;
		chainY[0] = empty_cell_y;
	}










	// Every ecDNA is copied once
	doubled_ecDNA_copyNumber = tissue[cell_x][cell_y].ecDNA * 2;
	//cout << "Dividing cell ecDNA copy number = " << tissue[cell_x][cell_y].ecDNA << ". Copied to " << doubled_ecDNA_copyNumber << endl;



	if (clustering == false)	// No ecDNA clustering
	{	
		// Distribute ecDNA between two daughter cells according to binomial
		binomial_distribution<int> draw_new_ecDNA_copyNumber(doubled_ecDNA_copyNumber , 0.5);
		daughter_ecDNA_copyNumber1 = draw_new_ecDNA_copyNumber(generator);
		daughter_ecDNA_copyNumber2 = doubled_ecDNA_copyNumber - daughter_ecDNA_copyNumber1;
		//cout << "Daugher cells receive " << daughter_ecDNA_copyNumber1 << " and " << daughter_ecDNA_copyNumber2 << " ecDNA\n" << endl;


		// Create daughter cell
		tissue[cell_x][cell_y].set_ecDNA(daughter_ecDNA_copyNumber1);
		tissue[chainX[0]][chainY[0]].set_ecDNA(daughter_ecDNA_copyNumber2);


		// Update total sum of birth rates in the system after cell division
		*r_birth -= (1.0 + selection_coeff*fitness(int(0.5*doubled_ecDNA_copyNumber) , weibull_fitness));		// Subtract contribution from mother cell
		*r_birth += (1.0 + selection_coeff*fitness(daughter_ecDNA_copyNumber1 , weibull_fitness));		// Add contribution from daughter cell 1
		*r_birth += (1.0 + selection_coeff*fitness(daughter_ecDNA_copyNumber2 , weibull_fitness));		// Add contribution from daughter cell 2

	}


	else 	// ecDNA cluster with Poisson cluster size distribution
	{
  		// Distribute ecDNA in clusters amongst two daughter cells 
  		motherCell_copyNumber = doubled_ecDNA_copyNumber;
  		daughter_ecDNA_copyNumber1 = 0;
  		daughter_ecDNA_copyNumber2 = 0;

  		//cout << "\n\nMother cell (ecDNA = " << motherCell_copyNumber << ") distributes ecDNA hubs across daughter cells" << endl;

  		while(motherCell_copyNumber > 0)
  		{
  			// Draw cluster size from Poisson distribution
  			clusterSize = ecDNA_clusterSizeDistribution(generator);
  			
  			// Skip if cluster size is larger than available remaining ecDNA in mother cell
  			if (clusterSize > motherCell_copyNumber) continue;

  			// Choose which daughter cell to put ecDNA cluster (binomial distribution of ecDNA clusters)
  			rand_double = drand48();

  			if (rand_double < 0.5) daughter_ecDNA_copyNumber1 += clusterSize;

  			// If not chosen to be placed into daughter cell 1, it goes in daughter cell 2
  			else daughter_ecDNA_copyNumber2 += clusterSize;

  			motherCell_copyNumber -= clusterSize;
  		}

  		//cout << "Finished distributing ecDNA hubs. Daughter cell 1 has " << daughter_ecDNA_copyNumber1 << " ecDNA, and daughter cell 2 has " << daughter_ecDNA_copyNumber2 << endl;
	
  		// Create daughter cell
		tissue[cell_x][cell_y].set_ecDNA(daughter_ecDNA_copyNumber1);
		tissue[chainX[0]][chainY[0]].set_ecDNA(daughter_ecDNA_copyNumber2);


		// Update total sum of birth rates in the system after cell division
		*r_birth -= (1.0 + selection_coeff*fitness(int(0.5*doubled_ecDNA_copyNumber) , weibull_fitness));		// Subtract contribution from mother cell
		*r_birth += (1.0 + selection_coeff*fitness(daughter_ecDNA_copyNumber1 , weibull_fitness));		// Add contribution from daughter cell 1
		*r_birth += (1.0 + selection_coeff*fitness(daughter_ecDNA_copyNumber2 , weibull_fitness));		// Add contribution from daughter cell 2

	}



	*Ntot += 1;


	if ((chainX[queue] != empty_cell_x) || (chainY[queue] != empty_cell_y))
	{
		cout << "Error with pushing algorithm!" << endl;
		exit(0);
	}
}






// Parse command line arguments (Flags and numerical arguments)
void parse_command_line_arguments(int argc, char** argv , bool *quiet , int *seed , int *q , bool *clustering , bool *weibull_fitness , int *initial_copyNumber , double *selection_coeff)
{
	int c;

	while ((c = getopt (argc, argv, "vCFx:q:n:s:")) != -1)
	switch (c)
	{
		// Verbose flag
		case 'v':
			*quiet = false;
			break;

		// Random seed
		case 'x':
			*seed = atoi(optarg);		
			break;

		// Set pushing limit for division algorithm
		case 'q':
			*q = atoi(optarg);		
			break;

		// Clustering parameter (C=0 for no ecDNA clustering, C=1 for clustering)
		case 'C':
			*clustering = true;
			break;

		case 'F':
			*weibull_fitness = true;
			break;

		// ecDNA copy number in initial cell
		case 'n':
			*initial_copyNumber = atoi(optarg);
			break;

		// Selection coefficient
		case 's':
			*
			selection_coeff = atof(optarg);
			break;

		case '?':
			if (optopt == 'c')
				fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			else if (isprint (optopt))
				fprintf (stderr, "Unknown option `-%c'.\n", optopt);
			else
				fprintf (stderr,
				"Unknown option character `\\x%x'.\n",
				optopt);
		//return 1;
		default:
		abort ();
	}





	// Checks on input parameter values
	if (*q < 0)
	{
		cout << "Pushing parameter q must be greater than 0. Exiting." << endl;
		exit(0);
	}

	if (*initial_copyNumber < 0)
	{
		cout << "Initial copy number must be 0 or greater. Exiting." << endl;
		exit(0);
	}
}






// Set up tissue (i.e. array of cells)
Cell** initialise_tissue(int _maxsize , int *Ntot , int initial_copyNumber)
{
	// Estimate radius of final system based on _maxsize
	radius_double = pow ( (_maxsize/M_PI) , (1.0/2.0) );


	// Slightly over-estimate this to avoid segmentation errors
	radius = (int)(1.5*radius_double);


	if (!quiet) cout << " " << endl;


	// Set up the system of cells, called tissue
	Cell ** tissue = new Cell*[2*radius];
	for (int i = 0; i < (2*radius); i++)
	{
		tissue[i] = new Cell[2*radius];

		for (int j = 0; j < (2*radius); j++)
		{
			tissue[i][j].set_ecDNA(-1);		// Cells with ecDNA copy number = -1 represent empty lattice points
		}

		if (!quiet) printf(" Initialising tissue... %i%%\r", (int)((i+1)*100.0/(2*radius)));
		if (!quiet) fflush(stdout);
	}
	if (!quiet) printf(" Initialising tissue... Done.\r");
	if (!quiet) cout << " " << endl;
		




	// Seed first tissue cell at (x,y) = (radius , radius)
	tissue[radius][radius].set_ecDNA(initial_copyNumber);

	*Ntot += 1;

	return tissue;
}





// Compute normalised birth and death rates 
void compute_normalised_birth_and_death_rates(int Ntot , double r_birth , double *r_birth_normalised , double *r_death_normalised)
{
	r_death = 0.5*Ntot;

	// Compute normalised reaction rates
	*r_birth_normalised = r_birth/(r_birth + r_death);
	*r_death_normalised = r_death/(r_birth + r_death);
}






// Choose next event in Gillespie algorithm
void choose_next_event(bool *BIRTH , bool *DEATH , double r_birth_normalised , double r_death_normalised)
{
	*BIRTH = false;
	*DEATH = false;

	rand_double = drand48();
	if (rand_double <= r_birth_normalised)
	{
		*BIRTH = true;
		//cout << "BIRTH" << endl;
	}
	else if (rand_double <= r_birth_normalised + r_death_normalised)
	{
		*DEATH = true;
		//cout << "DEATH" << endl;
	}
	else
	{
		cout << "Problem with Gillespie rates..." << endl;
		exit(0);
	}
}






// Select cell in tissue with all cells having equal probability of being chosen
void select_cell_flat_probability(Cell ** tissue , int *cell_x , int *cell_y , int radius , int x_b , int y_b)
{
	// Randomly select one cell to die (all birth rates are equal)
	*cell_x = 0;
	*cell_y = 0;

	do
	{
		*cell_x = (int)((2*(x_b))*drand48()) + radius - x_b;
		*cell_y = (int)((2*(y_b))*drand48()) + radius - y_b;
	}
	while (tissue[*cell_x][*cell_y].ecDNA == -1);
}






// Select cell in tissue when cells have different birth rates
void select_cell_unequal_birth_rates(Cell ** tissue , int *cell_x , int *cell_y , int radius)
{
	// Choose cell to divide based on its division rate
	birth_rate_sum = 0.0;
	*cell_x = -1;
	*cell_y = -1;
	rand_double = drand48();

	// Loop over all cells in system
	for (int i = 0; i < (2*radius); i++)
	{
		for (int j = 0; j < (2*radius); j++)
		{
			if (tissue[i][j].ecDNA == -1) continue; 

			// Add contribution of cell
			birth_rate_sum += (1.0 + selection_coeff*fitness(tissue[i][j].ecDNA , weibull_fitness));
			//cout << "Birth rate sum = " << birth_rate_sum << endl;

			//cout << "Checking if " << rand_double << " < " << birth_rate_sum/r_birth << endl;
			if (rand_double <= birth_rate_sum/r_birth)
			{
				*cell_x = i;
				*cell_y = j;
				break;
			}
		}
		if ((*cell_x != -1) && (*cell_y != -1)) break;
	}
}






// Kill cell and remove from lattice
void kill_cell(Cell ** tissue , double *r_birth , int cell_x , int cell_y , double selection_coeff , int *Ntot)
{
	// Subtract contribution of cell to r_birth
	*r_birth -= (1.0 + selection_coeff*fitness(tissue[cell_x][cell_y].ecDNA , weibull_fitness))
	;

	// Cell dies
	tissue[cell_x][cell_y].set_ecDNA(-1);
	*Ntot -= 1;
}
















/*******************************************************************************/










int main(int argc, char** argv)
{

	// Query number of available cores
	//unsigned concurentThreadsSupported = std::thread::hardware_concurrency();



	// Reset time and tissue size variables
	t = 0.0;
	Ntot = 0;
	selection_coeff = 0.0;
	q = 0;


	//================== Parse command line arguments ====================//
	parse_command_line_arguments(argc , argv , &quiet , &seed , &q , &clustering , &weibull_fitness , &initial_copyNumber , &selection_coeff);




	// Seed random number generator
	srand48(seed);
        default_random_engine generator(seed);





	//================== Initialise tissue ====================//
	Cell ** tissue = initialise_tissue(_maxsize , &Ntot , initial_copyNumber);





	//================== Simulate tissue growth ==================//
	iter = 0;
	x = 0;
	y = 0;
	x_b = 10;
	y_b = 10;

	// Initial total sum of all birth rates in system (initial cell with initial_copyNumber ecDNA)
	r_birth = 1.0 + selection_coeff*fitness(initial_copyNumber , weibull_fitness);


	do
	{
		
		++iter;


		// Re-evaluate birth and death rates
		compute_normalised_birth_and_death_rates(Ntot , r_birth , &r_birth_normalised , &r_death_normalised);



		// After initial expansion, update x- and y-bounds to include entire space
		if (Ntot > 0.1*_maxsize)
		{
			x_b = radius;
			y_b = radius;
		}




		// Choose birth or death based on normalised rates
		choose_next_event(&BIRTH , &DEATH , r_birth_normalised , r_death_normalised);





		// If division:
		if (BIRTH)
		{
			// Randomly select one cell to divide (all birth rates are equal)
			if (selection_coeff == 0) select_cell_flat_probability(tissue , &cell_x , &cell_y , radius , x_b , y_b);

			// Choose cell to divide based on its division rate
			else select_cell_unequal_birth_rates(tissue , &cell_x , &cell_y , radius);

			// Cell divides
			straight_line_division(tissue , cell_x , cell_y , &Ntot , &x_b , &y_b , &r_birth , radius);
		}





		// If death:
		if ((DEATH) && (Ntot > 1))
		{
			// Randomly select one cell to die (all birth rates are equal)
			select_cell_flat_probability(tissue , &cell_x , &cell_y , radius , x_b , y_b);

			// Cell dies
			kill_cell(tissue , &r_birth , cell_x , cell_y , selection_coeff , &Ntot);
		}




		if ((!quiet) && (iter%1000 == 0))
		{
			cout << "Iteration #" << iter << " -- N = " << Ntot << endl;
		}






	} while (Ntot < _maxsize);		// Exit once system has reached total size of _maxsize

	if (!quiet) cout << " " << endl;










	//================== Open data files & write final system data ==================//



	stringstream f;
	f.str("");
	f << "./2D_DATA/Nmax=" << _maxsize << "_initialCopyNumber=" << initial_copyNumber << "_q=" << q << "_s=" << selection_coeff << "_clustering=" << boolalpha << clustering << "_copyNumberDependentSelection=" << boolalpha << weibull_fitness << "/seed=" << seed;
	DIR *dir = opendir(f.str().c_str());
	if(!dir)
	{
		f.str("");
		f << "mkdir -p ./2D_DATA/Nmax=" << _maxsize << "_initialCopyNumber=" << initial_copyNumber << "_q=" << q << "_s=" << selection_coeff << "_clustering=" << boolalpha << clustering << "_copyNumberDependentSelection=" << boolalpha << weibull_fitness << "/seed=" << seed;
		system(f.str().c_str());
	}

	ofstream tissue_file;
	f.str("");
	f << "./2D_DATA/Nmax=" << _maxsize << "_initialCopyNumber=" << initial_copyNumber << "_q=" << q << "_s=" << selection_coeff << "_clustering=" << boolalpha << clustering << "_copyNumberDependentSelection=" << boolalpha << weibull_fitness << "/seed=" << seed << "/tissue.csv";
	tissue_file.open(f.str().c_str());


	if (!quiet) cout << " " << endl;
	if (!quiet) cout << "Created output files..." << endl;





	// Write tissue data to file
	for (int i = 0; i < (2*radius); i++)
	{
		for (int j = 0; j < (2*radius); j++)
		{
			tissue_file << i << "," << j << "," << tissue[i][j].ecDNA << endl;
		}
	}



	tissue_file.close();









	return 0;
}

















