
# include <iostream>
# include <fstream>
# include <sstream>
# include <stdlib.h>
# include <stdio.h>
# include <cmath>
# include <math.h>
# include <random>
# include <cstdlib>
# include <dirent.h>
# include <string>
# include <getopt.h>
# include <algorithm>
# include <initializer_list>



using namespace std;



/*******************************************************************************/



// Declare variables
#define PI 3.14159265


double radius_double, t, r_birth , r_birth_ecDNA_negative , r_birth_ecDNA_positive , r_death, r_birth_normalised , r_birth_ecDNA_negative_normalised , r_birth_ecDNA_positive_normalised , r_death_normalised, rand_double, cell_cell_chord_gradient, cell_cell_chord_yIntercept, min_move_distance, selection_coeff;
double optimal_direction_i,  optimal_direction_j, optimal_vector_norm, vector_norm, rescaled_min_length, scalar_prod, dist, nearest_space_distance, x_boundary_intersect, y_boundary_intersect;
double weibull_shape, weibull_scale, birth_rate_sum;
int Nmax, q, seed, radius, Ntot, N_ecDNA_hot, iter, x, y, cell_x, cell_y, dir, queue_size, ind, length, coordX, coordY, previous_link_direction, chain_length;
int arising_time, x_b, y_b, direction, chosen_direction, min_length, num_mins, chosen_min, doubled_ecDNA_copyNumber, daughter_ecDNA_copyNumber1, daughter_ecDNA_copyNumber2;
int empty_cell_x, empty_cell_y, search_radius, num_nearest_empty_cells, rand_int, vertical_direction, horizontal_direction, i_lowerBound, j_lowerBound, i_upperBound, j_upperBound;
int next_move_direction_x, next_move_direction_y, motherCell_copyNumber, clusterSize, initial_copyNumber, search_radius_incrememnt, minimum_edge_distance;


// Define arrays which will contain relative coordinates of empty neighbours for a chosen cell 
vector<int> chainX(1001);
vector<int> chainY(1001);
vector<int> moves_direction_x(1001);
vector<int> moves_direction_y(1001);
vector<double> moves_distances(1001);
vector<int> previous_moves_direction_x(1001);
vector<int> previous_moves_direction_y(1001);
vector<double> previous_moves_distances(1001);

vector<int> nearest_empty_cells_x;
vector<int> nearest_empty_cells_y;



bool BIRTH_ecDNA_negative = false;
bool BIRTH_ecDNA_positive = false;
bool DEATH = false;
bool skip = false;
bool verbose_flag = false;
bool clustering_flag = false;
bool stop_searching = false;

double search_radius_constant = cos(PI/8)*0.765375;



/*******************************************************************************/



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






//-----------------------






// Find nearest empty lattice points from cell at given coordinates (x,y) = (cell_x,cell_y), and count the number of jointly nearest lattice points if more than one
void find_nearest_empty_lattice_point(Cell** tissue , int cell_x , int cell_y , int* search_radius , int* search_radius_incrememnt , int minimum_edge_distance , double* nearest_space_distance , int* num_nearest_empty_cells , int* empty_cell_x , int* empty_cell_y)
{

	stop_searching = false;
	*empty_cell_x = -1;
	*empty_cell_y = -1;
	*nearest_space_distance = 1e10;
	*num_nearest_empty_cells = 0;

	*search_radius_incrememnt = ceil(q/50.0);
	*search_radius = 1 - *search_radius_incrememnt;

	// Make sure we don't look off the edges of the array (this is necessary for very large increment sizes e.g. 500)
	// Need to find the nearest 'edge' of the tissue array for the given moether cell co-ordinates
	//minimum_edge_distance = std::min({cell_x, (2*radius)-cell_x, cell_y, (2*radius)-cell_y});

	if ((cell_x < (2*radius)-cell_x) && (cell_x < cell_y) && (cell_x < (2*radius)-cell_y)) minimum_edge_distance = cell_x;
	else if (((2*radius)-cell_x < cell_x) && ((2*radius)-cell_x < cell_y) && ((2*radius)-cell_x < (2*radius)-cell_y)) minimum_edge_distance = (2*radius)-cell_x;
	else if ((cell_y < cell_x) && (cell_y < (2*radius)-cell_x) && (cell_y < (2*radius)-cell_y)) minimum_edge_distance = cell_y;
	else minimum_edge_distance = (2*radius)-cell_y;

	while (true)
	{

		// Make sure we don't look further than the max value, q, or off the edges of the array (this is necessary for very large increment sizes e.g. 500)
		if (*search_radius + *search_radius_incrememnt >= minimum_edge_distance)
		{
			*search_radius_incrememnt = minimum_edge_distance - 2;
			stop_searching = true;
		}
		if (*search_radius + *search_radius_incrememnt > q)
		{
			*search_radius_incrememnt = q - *search_radius;
			stop_searching = true;
		}


		*search_radius += *search_radius_incrememnt;


		// Search within a radius of search_radius_incrememnt from dividing cell for nearest empty lattice point 
		// Loop through all cells within box of size 2*search_radius_incrememnt centred around dividing cell
		for (int i = -(*search_radius); i <= *search_radius; ++i)
		{
			for (int j = -(*search_radius); j <= *search_radius; ++j)
			{

				// Compute distance between this cell and dividing cell
				dist = (i*i) + (j*j);


				// Only check cells within distance search_radius of dividing cell
				if (dist > float((*search_radius)*(*search_radius))) continue;


				// Check if lattice coordinates are occupied or empty
				if (tissue[cell_x + i][cell_y + j].ecDNA == -1)
				{
					// If found new nearest empty lattice point, update nearest distance and empty cell coordinates
					if (dist < *nearest_space_distance)
					{
						*nearest_space_distance = dist;
						*empty_cell_x = cell_x + i;
						*empty_cell_y = cell_y + j;

						*num_nearest_empty_cells = 1;
					}

					// If found another empty cell equally close to current nearest empty cell, simply keep track of number of joint nearest cells
					else if (dist == *nearest_space_distance) 
					{
						*num_nearest_empty_cells += 1;
					}
				}
			}
		}

		if ((stop_searching == true) || ((*empty_cell_x != -1) && (*empty_cell_y != -1))) break;

	}
}






//-----------------------






// If >1 joint nearest empty lattice point from cell at (x,y) = (cell_x,cell_y), select one with uniform probability 
void choose_nearest_empty_lattice_point(Cell ** tissue , int cell_x , int cell_y , vector<int>& nearest_empty_cells_x , vector<int>& nearest_empty_cells_y , int search_radius , int nearest_space_distance , int *empty_cell_x , int *empty_cell_y)
{

	// Clear vector of co-ordinates to nearest empty lattice points 
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
			}
		}
	}



	// Choose which of the cells will divide
	rand_int = int(round(drand48()*(num_nearest_empty_cells-1)));
	*empty_cell_x = nearest_empty_cells_x[rand_int];
	*empty_cell_y = nearest_empty_cells_y[rand_int];
}






//-----------------------






void construct_cell_pushing_path(int cell_x , int cell_y , int empty_cell_x , int empty_cell_y , int* queue_size , vector<int>& moves_direction_x , vector<int>& moves_direction_y , vector<double>& moves_distances , vector<int>& previous_moves_direction_x , vector<int>& previous_moves_direction_y , vector<double>& previous_moves_distances , vector<int>& chainX , vector<int>& chainY)
{

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

	*queue_size = 0;



	for (int i = i_lowerBound+1; i < i_upperBound; ++i)
	{

		// Find coordinates of intersection between line x = i and cell-cell chord
		y_boundary_intersect = (cell_cell_chord_gradient * float(i)) + cell_cell_chord_yIntercept;

		// Find (squared) Euclidean distance between (i , y_boundary_intersect) and (cell_x , cell_y)
		dist = pow((cell_x - i) , 2) + pow((cell_y - y_boundary_intersect) , 2);

		// Store move in list
		moves_direction_x[*queue_size] = horizontal_direction;
		moves_direction_y[*queue_size] = 0;

		// Store distance in list
		moves_distances[*queue_size] = dist;

		*queue_size += 1;
	}




	for (int j = j_lowerBound+1; j < j_upperBound; ++j)
	{

		// Find coordinates of intersection between line x = i and cell-cell chord
		if (empty_cell_x != cell_x) x_boundary_intersect = (float(j) - cell_cell_chord_yIntercept)/cell_cell_chord_gradient;
		else x_boundary_intersect = cell_x;

		// Find (squared) Euclidean distance between (i , y_boundary_intersect) and (cell_x , cell_y)
		dist = pow((cell_x - x_boundary_intersect) , 2) + pow((cell_y - j) , 2);

		// Store move in list
		moves_direction_x[*queue_size] = 0;
		moves_direction_y[*queue_size] = vertical_direction;

		// Store distance in list
		moves_distances[*queue_size] = dist;

		*queue_size += 1;
	}





	// Reset lists containing previous move information 
	for (int i = 0; i < *queue_size; ++i)
	{
		previous_moves_distances[i] = 0.0;
		previous_moves_direction_x[i] = -10;
		previous_moves_direction_y[i] = -10;
	}





	// All moves have been determined, but are not yet in order of distance from dividing cell. 
	// So, loop through moves and their distances, and construct chain of cell coordinates.
	chain_length = 0;
	while (chain_length != *queue_size)
	{
		min_move_distance = 1e10;

		for (int i = 0; i < *queue_size; ++i)
		{
			skip = false;

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
		}

		if (chain_length == 0)
		{
			chainX[chain_length] = cell_x + next_move_direction_x;
			chainY[chain_length] = cell_y + next_move_direction_y;

			previous_moves_distances[chain_length] = min_move_distance;
			previous_moves_direction_x[chain_length] = next_move_direction_x;
			previous_moves_direction_y[chain_length] = next_move_direction_y;

			chain_length += 1;

		}
		else
		{
			chainX[chain_length] = chainX[chain_length - 1] + next_move_direction_x;
			chainY[chain_length] = chainY[chain_length - 1] + next_move_direction_y;

			previous_moves_distances[chain_length] = min_move_distance;
			previous_moves_direction_x[chain_length] = next_move_direction_x;
			previous_moves_direction_y[chain_length] = next_move_direction_y;

			chain_length += 1;
		}

	}


}






//-----------------------






void move_cells_along_path(Cell ** tissue , vector<int>& chainX , vector<int>& chainY , int* queue_size , int horizontal_direction , int vertical_direction , int *x_b , int *y_b , int radius , int empty_cell_x , int empty_cell_y)
{
	if (*queue_size > 0)
	{

		//cout << "queue_size = " << queue_size << endl;

		// Add two final moves (one vertical and one horizontal) to the end of the chain
		if (drand48() < 0.5)
		{
			chainX[*queue_size] = chainX[*queue_size - 1] + horizontal_direction;
			chainY[*queue_size] = chainY[*queue_size - 1] + 0;

			*queue_size += 1;

			chainX[*queue_size] = chainX[*queue_size - 1] + 0;
			chainY[*queue_size] = chainY[*queue_size - 1] + vertical_direction;
		}
		else
		{
			chainX[*queue_size] = chainX[*queue_size - 1] + 0;
			chainY[*queue_size] = chainY[*queue_size - 1] + vertical_direction;

			*queue_size += 1;

			chainX[*queue_size] = chainX[*queue_size - 1] + horizontal_direction;
			chainY[*queue_size] = chainY[*queue_size - 1] + 0;
		}


		for (int i = 0; i < *queue_size; ++i)
		{
			tissue[chainX[*queue_size-i]][chainY[*queue_size-i]].ecDNA = tissue[chainX[*queue_size-i-1]][chainY[*queue_size-i-1]].ecDNA;

			// Update bounds on tissue size
			if (fabs(chainX[i] + 1 - radius) > *x_b) *x_b = fabs(chainX[i] + 1 - radius);
			if (fabs(chainY[i] + 1 - radius) > *y_b) *y_b = fabs(chainY[i] + 1 - radius);
		}

	}

	else 	// Even if queue_size=0, check that newly created cell increases any bounds
	{
		// Update bounds on tissue size
		if (fabs(chainX[0] + 1 - radius) > *x_b) *x_b = fabs(chainX[0] + 1 - radius);
		if (fabs(chainY[0] + 1 - radius) > *y_b) *y_b = fabs(chainY[0] + 1 - radius);
	}



	if (*queue_size == 0)
	{
		chainX[0] = empty_cell_x;
		chainY[0] = empty_cell_y;
	}
}






//-----------------------






void choose_horizontal_and_vertical_move_directions(int cell_x , int cell_y , int empty_cell_x , int empty_cell_y , int* horizontal_direction , int* vertical_direction)
{
	*horizontal_direction = 0;
	*vertical_direction = 0;

	if (empty_cell_x > cell_x) *horizontal_direction = 1;
	else if (cell_x > empty_cell_x) *horizontal_direction = -1;

	if (empty_cell_y > cell_y) *vertical_direction = 1;
	else if (cell_y > empty_cell_y) *vertical_direction = -1;
}






//-----------------------






void compute_cell_cell_chord(double* cell_cell_chord_gradient , double* cell_cell_chord_yIntercept , int cell_x , int cell_y , int empty_cell_x , int empty_cell_y)
{
	*cell_cell_chord_gradient = 0.0;
	*cell_cell_chord_yIntercept = 0.0;
	if (empty_cell_x != cell_x)
	{
		*cell_cell_chord_gradient = float(empty_cell_y - cell_y)/float(empty_cell_x - cell_x);
		*cell_cell_chord_yIntercept = cell_y - (*cell_cell_chord_gradient * cell_x);
	}
}






//-----------------------






void distribute_ecDNA(Cell ** tissue , int cell_x , int cell_y , int* daughter_ecDNA_copyNumber1 , int* daughter_ecDNA_copyNumber2 , mt19937_64 *generator)
{
	// Every ecDNA is copied once
	doubled_ecDNA_copyNumber = tissue[cell_x][cell_y].ecDNA * 2;

	// Distribute ecDNA between two daughter cells according to binomial
	binomial_distribution<int> draw_new_ecDNA_copyNumber(doubled_ecDNA_copyNumber , 0.5);
	*daughter_ecDNA_copyNumber1 = draw_new_ecDNA_copyNumber(*generator);
	*daughter_ecDNA_copyNumber2 = doubled_ecDNA_copyNumber - *daughter_ecDNA_copyNumber1;
}






//-----------------------






// Straight line pushing 
void straight_line_division(Cell ** tissue , int cell_x , int cell_y , int *Ntot , int *N_ecDNA_hot , int *x_b , int *y_b , double *r_birth , int radius , mt19937_64 *generator)
{


	// Find nearest (or joint nearest) empty lattice points from dividing cell at (cell_x , cell_y)
	find_nearest_empty_lattice_point(tissue , cell_x , cell_y , &search_radius , &search_radius_incrememnt , minimum_edge_distance , &nearest_space_distance , &num_nearest_empty_cells , &empty_cell_x , &empty_cell_y);



	// If no empty lattice points within max radius of q from dividing cell, cell cannot divide. Exit.
	if ((empty_cell_x == -1) && (empty_cell_y == -1)) return;



	// If there is more than 1 empty cell at the nearest distance to the dividing cell, choose one with uniform probability
	if (num_nearest_empty_cells > 1)
	{
		choose_nearest_empty_lattice_point(tissue , cell_x , cell_y , nearest_empty_cells_x , nearest_empty_cells_y , search_radius , nearest_space_distance , &empty_cell_x , &empty_cell_y);
	}



	// Decide whether horizontal cell movements will be left or right, similarly for vertical movements
	choose_horizontal_and_vertical_move_directions(cell_x , cell_y , empty_cell_x , empty_cell_y , &horizontal_direction , &vertical_direction);



	// Draw straight line between dividing cell and empty lattice point
	compute_cell_cell_chord(&cell_cell_chord_gradient , &cell_cell_chord_yIntercept , cell_x , cell_y , empty_cell_x , empty_cell_y);



	// Construct chain by finding intersection of cell-cell chord with vertical & horizontal lines between (cell_x , cell_y) and (empty_cell_x , empty_cell_y)
	construct_cell_pushing_path(cell_x , cell_y , empty_cell_x , empty_cell_y , &queue_size , moves_direction_x , moves_direction_y , moves_distances , previous_moves_direction_x , previous_moves_direction_y , previous_moves_distances , chainX , chainY);

	

	// Once the chain has been constructed, move all cells along one place
	move_cells_along_path(tissue , chainX , chainY , &queue_size , horizontal_direction , vertical_direction , x_b , y_b , radius , empty_cell_x , empty_cell_y);



	// Nearby cells have been pushed and an empty space created next to dividing cell. Now deal with cell division and ecDNA inheritance
	distribute_ecDNA(tissue , cell_x , cell_y , &daughter_ecDNA_copyNumber1 , &daughter_ecDNA_copyNumber2 , generator);
	


	// Create daughter cell
	tissue[cell_x][cell_y].set_ecDNA(daughter_ecDNA_copyNumber1);
	tissue[chainX[0]][chainY[0]].set_ecDNA(daughter_ecDNA_copyNumber2);


	// Book-keeping
	if ((daughter_ecDNA_copyNumber1 > 0) && (daughter_ecDNA_copyNumber2 > 0)) *N_ecDNA_hot += 1;

	*Ntot += 1;



	if ((chainX[queue_size] != empty_cell_x) || (chainY[queue_size] != empty_cell_y))
	{
		cout << "Error with pushing algorithm!" << endl;
		exit(0);
	}
}






//-----------------------






// Parse command line arguments (Flags and numerical arguments)
void parse_command_line_arguments(int argc, char** argv , bool *verbose_flag , int *seed , int *Nmax , int *q , int *initial_copyNumber , double *selection_coeff)
{
	int c;
	int option_index;
	char* arg_long = nullptr;
	int verbose = 0;

	static struct option long_options[] =
	{
		{"verbose", no_argument, &verbose, 1},
	}; 

	while ((c = getopt_long(argc, argv, "x:q:n:s:N:", long_options, &option_index)) != -1)
	switch (c)
	{
		case 0:
		{
			arg_long = optarg;
			break;
		}

		// Random seed
		case 'x':
			*seed = atoi(optarg);		
			break;

		// Set pushing limit for division algorithm
		case 'q':
			*q = atoi(optarg);		
			break;

		// ecDNA copy number in initial cell
		case 'n':
			*initial_copyNumber = atoi(optarg);
			break;

		// Selection coefficient
		case 's':
			*selection_coeff = atof(optarg);
			break;

		// Selection coefficient
		case 'N':
			*Nmax = atoi(optarg);
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



	// Set boolean values for verbose flag
	if (verbose == 1) *verbose_flag = true;


	// Checks on input parameter values
	if (*q <= 0)
	{
		cout << "Pushing parameter q must be greater than 0. Exiting." << endl;
		exit(0);
	}


	if (*Nmax <= 0)
	{
		cout << "Maximum tumor size Nmax must be greater than 0. Exiting." << endl;
		exit(0);
	}


	if (*initial_copyNumber < 0)
	{
		cout << "Initial copy number must be 0 or greater. Exiting." << endl;
		exit(0);
	}
}






//-----------------------






// Set up tissue (i.e. array of cells)
Cell** initialise_tissue(int Nmax , int *Ntot , int *N_ecDNA_hot , int initial_copyNumber)
{

	
	// Estimate radius of final system based on Nmax
	radius_double = pow ( (Nmax/M_PI) , (1.0/2.0) );


	// Slightly over-estimate this to avoid segmentation errors
	radius = (int)(1.5*radius_double);


	if (verbose_flag) cout << " " << endl;


	// Set up the system of cells, called tissue
	Cell ** tissue = new Cell*[2*radius];
	for (int i = 0; i < (2*radius); i++)
	{
		tissue[i] = new Cell[2*radius];

		for (int j = 0; j < (2*radius); j++)
		{
			tissue[i][j].set_ecDNA(-1);		// Cells with ecDNA copy number = -1 represent empty lattice points
		}

		if (verbose_flag) printf(" Initialising tissue... %i%%\r", (int)((i+1)*100.0/(2*radius)));
		if (verbose_flag) fflush(stdout);
	}
	if (verbose_flag) printf(" Initialising tissue... Done.\r");
	if (verbose_flag) cout << " " << endl;
		


	// Seed first tissue cell at (x,y) = (radius , radius)
	tissue[radius][radius].set_ecDNA(initial_copyNumber);

	// Book-keeping
	*Ntot += 1;
	*N_ecDNA_hot += 1;

	return tissue;
}






//-----------------------






// Compute normalised birth and death rates 
void compute_normalised_birth_and_death_rates(int Ntot , int N_ecDNA_hot , double *r_birth_ecDNA_negative_normalised , double *r_birth_ecDNA_positive_normalised , double *r_death_normalised)
{
	// Compute un-normalised reaction rates
	r_death = 0.5*Ntot;
	r_birth_ecDNA_negative = (double)(Ntot - N_ecDNA_hot) * (1.0);
	r_birth_ecDNA_positive = (double)(N_ecDNA_hot) * (1.0 + selection_coeff);

	// Compute normalised reaction rates
	*r_birth_ecDNA_negative_normalised = r_birth_ecDNA_negative/(r_birth_ecDNA_negative + r_birth_ecDNA_positive + r_death);
	*r_birth_ecDNA_positive_normalised = r_birth_ecDNA_positive/(r_birth_ecDNA_negative + r_birth_ecDNA_positive + r_death);
	*r_death_normalised = r_death/(r_birth_ecDNA_negative + r_birth_ecDNA_positive + r_death);
}






//-----------------------






// Choose next event in Gillespie algorithm
void choose_next_event(bool *BIRTH_ecDNA_negative , bool *BIRTH_ecDNA_positive , bool *DEATH , double r_birth_ecDNA_negative_normalised , double r_birth_ecDNA_positive_normalised , double r_death_normalised)
{
	*BIRTH_ecDNA_negative = false;
	*BIRTH_ecDNA_positive = false;
	*DEATH = false;

	rand_double = drand48();
	if (rand_double <= r_birth_ecDNA_negative_normalised)
	{
		*BIRTH_ecDNA_negative = true;
		//cout << "BIRTH (ecDNA-)" << endl;
	}
	else if (rand_double <= r_birth_ecDNA_negative_normalised + r_birth_ecDNA_positive_normalised)
	{
		*BIRTH_ecDNA_positive = true;
		//cout << "BIRTH (ecDNA+)" << endl;
	}
	else if (rand_double <= r_birth_ecDNA_negative_normalised + r_birth_ecDNA_positive_normalised + r_death_normalised)
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






//-----------------------






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







//-----------------------







// Kill cell and remove from lattice
void kill_cell(Cell ** tissue , int cell_x , int cell_y , double selection_coeff , int *Ntot , int *N_ecDNA_hot)
{

	if (tissue[cell_x][cell_y].ecDNA > 0) *N_ecDNA_hot -= 1;

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
	N_ecDNA_hot = 0;
	selection_coeff = 0.0;
	q = 0;
	Nmax = 0;


	//================== Parse command line arguments ====================//
	parse_command_line_arguments(argc , argv , &verbose_flag , &seed , &Nmax , &q , &initial_copyNumber , &selection_coeff);



	// Seed random number generator
	srand48(seed);
	mt19937_64 generator;
	generator.seed(seed);





	//================== Initialise tissue ====================//
	Cell ** tissue = initialise_tissue(Nmax , &Ntot , &N_ecDNA_hot , initial_copyNumber);





	//================== Simulate tissue growth ==================//
	iter = 0;
	x = 0;
	y = 0;
	x_b = 10;
	y_b = 10;


	do
	{
		
		++iter;


		// Re-evaluate birth and death rates
		compute_normalised_birth_and_death_rates(Ntot , N_ecDNA_hot , &r_birth_ecDNA_negative_normalised , &r_birth_ecDNA_positive_normalised , &r_death_normalised);



		// After initial expansion, update x- and y-bounds to include entire space
		if (Ntot > 0.1*Nmax)
		{
			x_b = radius;
			y_b = radius;
		}




		// Choose birth or death based on normalised rates
		choose_next_event(&BIRTH_ecDNA_negative , &BIRTH_ecDNA_positive , &DEATH , r_birth_ecDNA_negative_normalised , r_birth_ecDNA_positive_normalised , r_death_normalised);





		// If division:
		if (BIRTH_ecDNA_negative)
		{
			// Randomly select one cell to divide (all birth rates are equal)
			while (true)
			{
				select_cell_flat_probability(tissue , &cell_x , &cell_y , radius , x_b , y_b);
				if (tissue[cell_x][cell_y].ecDNA == 0) break;
			}

			// Cell divides
			straight_line_division(tissue , cell_x , cell_y , &Ntot , &N_ecDNA_hot , &x_b , &y_b , &r_birth , radius , &generator);
		}




		if (BIRTH_ecDNA_positive)
		{
			// Randomly select one cell to divide (all birth rates are equal)
			while (true)
			{
				select_cell_flat_probability(tissue , &cell_x , &cell_y , radius , x_b , y_b);
				if (tissue[cell_x][cell_y].ecDNA > 0) break;
			}

			// Cell divides
			straight_line_division(tissue , cell_x , cell_y , &Ntot , &N_ecDNA_hot , &x_b , &y_b , &r_birth , radius , &generator);
		}





		// If death:
		if ((DEATH) && (Ntot > 1))
		{
			// Randomly select one cell to die (all birth rates are equal)
			select_cell_flat_probability(tissue , &cell_x , &cell_y , radius , x_b , y_b);

			// Cell dies
			kill_cell(tissue , cell_x , cell_y , selection_coeff , &Ntot , &N_ecDNA_hot);
		}




		if ((verbose_flag) && (iter%1000 == 0))
		{
			cout << "Iteration #" << iter << " -- N = " << Ntot << " -- N_ecDNA_hot = " << N_ecDNA_hot << endl;
		}






	} while (Ntot < Nmax);		// Exit once system has reached total size of Nmax

	if (verbose_flag) cout << " " << endl;










	//================== Open data files & write final system data ==================//



	stringstream f;
	f.str("");
	f << "./results/Nmax=" << Nmax << "_n=" << initial_copyNumber << "_q=" << q << "_s=" << selection_coeff << "/seed=" << seed;
	DIR *dir = opendir(f.str().c_str());
	if(!dir)
	{
		f.str("");
		f << "mkdir -p ./results/Nmax=" << Nmax << "_n=" << initial_copyNumber << "_q=" << q << "_s=" << selection_coeff << "/seed=" << seed;
		system(f.str().c_str());
	}

	ofstream tissue_file;
	f.str("");
	f << "./results/Nmax=" << Nmax << "_n=" << initial_copyNumber << "_q=" << q << "_s=" << selection_coeff << "/seed=" << seed << "/tissue.csv";
	tissue_file.open(f.str().c_str());


	if (verbose_flag) cout << " " << endl;
	if (verbose_flag) cout << "Created output files..." << endl;





	// Write tissue data to file
	for (int i = 0; i < (2*radius); i++)
	{
		for (int j = 0; j < (2*radius); j++)
		{
			if (tissue[i][j].ecDNA != -1) tissue_file << i << "," << j << "," << tissue[i][j].ecDNA << endl;
		}
	}



	tissue_file.close();









	return 0;
}

















