
#include <iostream>

#define AVAIL_MEM_SIZE 32768
#define AVAIL_MEM (AVAIL_MEM_SIZE / 4)

void inplace_transpose_square(int *x, const int &L, int lda)
{
	for(int l=0; l<L; l++)
		for(int m=0; m<L; m++)
		{
			if(m<l)
				continue;

			int tmp;
			tmp = x[l*lda + m];
			x[l*lda + m] = x[m*lda + l];
			x[m*lda + l] = tmp;
		}
}

void snake(int *x, int *ts, int *td, const int &L, const int &P, int is, int id, int pos)
{
	// pos 0 - start, 1 - inter, 2 - end

	for(int p=0; p<P; p++)
	{
		for(int j=0; j<L; j++)
		{
			if(pos == 0)
			{
				ts[p*L + j] = x[is*P*L + p*L + j];
				td[p*L + j] = x[id*P*L + p*L + j];				
				x[id*P*L + p*L + j] = ts[p*L + j];		
			}
			else if(pos == 1)
			{
				td[p*L + j] = x[id*P*L + p*L + j];				
				x[id*P*L + p*L + j] = ts[p*L + j];		
			}
			else
			{
				x[id*P*L + p*L + j] = ts[p*L + j];
			}
		}
	}
}
/* -> get_cycles funcition gets the swapping logic required for given row x col matrix.
   -> cycle_map[0] holds the total number of cycles required. 
   -> cycles start and end with the same index, hence we can identify individual cycles,
   though we tend to store the cycle index contiguousl.y*/
void get_cycles(int *cycle_map, int num_reduced_row, int num_reduced_col)
{
    int *is_swapped = new int[num_reduced_row * num_reduced_col];
    int i, map_index = 1, num_cycles = 0;
    int swap_id;
    /*initialize swap map*/
    is_swapped[0] = 1;
    is_swapped[num_reduced_row * num_reduced_col - 1] = 1;
    for (i = 1; i < (num_reduced_row * num_reduced_col - 1); i++)
    {
        is_swapped[i] = 0;
    }

    for (i = 1; i < (num_reduced_row * num_reduced_col - 1); i++)
    {
        swap_id = i;
        while (!is_swapped[swap_id])
        {
            is_swapped[swap_id] = 1;
            cycle_map[map_index++] = swap_id;
            swap_id = (num_reduced_row * swap_id) % (num_reduced_row * num_reduced_col - 1);
            if (swap_id == i)
            {
                cycle_map[map_index++] = swap_id;
                num_cycles++;
            }
        }
    }
    cycle_map[0] = num_cycles;
}

/* This function factorizes L and it finds a maximum of
   the 'factors less than max_capacity'*/
int get_num_lines_to_be_loaded(int max_capacity, int L)
{
    if ( L < max_capacity)
    {
        return L;
    }

    int square_root = (int)sqrt(L) + 1;
    int max_factor = 1;
    for (int i = 1; i < square_root; i++)
    {
        if (L % i == 0)
        {
            if (( i > max_factor) && (i <= max_capacity))
            {
                max_factor = i;
            }

            if (((L / i) > max_factor) && ((L / i) <= max_capacity))
            {
                max_factor = L / i;
            }
        }
    }
    return max_factor;
}

void inplace_1_isto_2_transpose_generic(int *x, const int &L, const int &M)
{
    int *te = new int[AVAIL_MEM >> 1];
    int *to = new int[AVAIL_MEM >> 1];
    int max_capacity = (AVAIL_MEM >> 1) / L;
    if (max_capacity <= 0)
    {
        std::cout << "\nIn-place transpose cannot be performed within specified memory constraints.\n";
        exit(1);
    }
    int num_lines_loaded = get_num_lines_to_be_loaded(max_capacity, L);
    int num_reduced_row = std::ceil( (float) M / (float) (2 * num_lines_loaded)); 
    int num_reduced_col = 2;
    /* The reduced row and col comes from the fact that we perform swaps of 'num_lines_loaded' rows first
       and thus reducing the amount of swaps required to be done using the snake function*/
    int i;

    if (num_lines_loaded > 1)
    {
        for (i = 0; i < M; i += 2 * num_lines_loaded)
        {
            // read
            for (int p = 0; p < num_lines_loaded; p++)
            {
                for (int j = 0; j < L; j++)
                {
                    te[p*L + j] = x[i*L + (2 * p + 0)*L + j];
                    to[p*L + j] = x[i*L + (2 * p + 1)*L + j];
                }
            }

            // write
            for (int p = 0; p < num_lines_loaded; p++)
            {
                for (int j = 0; j < L; j++)
                {
                    x[i*L + 0 + p*L + j] = te[p*L + j];
                    x[i*L + num_lines_loaded*L + p*L + j] = to[p*L + j];
                }
            }
        }
    }
    int *cycle_map = new int[num_reduced_row * num_reduced_col * 2];
    /* The memory required by cycle_map canniot exceed 2 times row*col by design*/

    get_cycles(cycle_map, num_reduced_row, num_reduced_col);

    int *ta = te;
    int *tb = to;
    int *ttmp;

    int inx = 0, start_inx;
    for (i = 0; i < cycle_map[0]; i++)
    {
        start_inx = cycle_map[++inx];
        std::cout << "\nCycle:" << (i + 1) <<">\t"<< "("<<start_inx <<","<< cycle_map[inx + 1] <<")";
        snake(x, ta, tb, L, num_lines_loaded, start_inx, cycle_map[inx + 1], 0);

        while (start_inx != cycle_map[++inx])
        {
            ttmp = ta;
            ta = tb;
            tb = ttmp;

            std::cout <<"\t" << "(" <<cycle_map[inx] << "," << cycle_map[inx + 1] << ")";
            int action_var = (cycle_map[inx + 1] == start_inx) ? 2 : 1;
            snake(x, ta, tb, L, num_lines_loaded, cycle_map[inx], cycle_map[inx + 1], action_var);
        }
    }

    delete[] te;
    delete[] to;
}

int main()
{
	const int L = 128;
	const int M = 2*L;

	int *x, *y, *z;

	x = new int[L*M];
	y = new int[L*M];
	z = new int[L*M];

	for(int l=0; l<L; l++)
		for(int m=0; m<M; m++)
		{
			x[l*M + m] = rand()%1000;
		}


	for(int l=0; l<L; l++)
		for(int m=0; m<M; m++)
		{
			y[m*L + l] = x[l*M + m];
		}


	memcpy(z, x, L*M*sizeof(int));

	// for(int i=0; i<L*M; i++)
	// {
	// 	if(!(i%M))
	// 		std::cout << std::endl;
	// 
	// 	std::cout.width(4);
	// 	std::cout << z[i] << " ";
	// }
	// std::cout << std::endl << std::endl;

	inplace_transpose_square(z, L, M);

	// for(int i=0; i<L*M; i++)
	// {
	// 	if(!(i%M))
	// 		std::cout << std::endl;
	// 
	// 	std::cout.width(4);
	// 	std::cout << z[i] << " ";
	// }
	// std::cout << std::endl << std::endl;

	inplace_transpose_square(z+L, L, M);

	// for(int i=0; i<L*M; i++)
	// {
	// 	if(!(i%M))
	// 		std::cout << std::endl;
	// 
	// 	std::cout.width(4);
	// 	std::cout << z[i] << " ";
	// }
	// std::cout << std::endl << std::endl;

    inplace_1_isto_2_transpose_generic(z, L, M);

	for(int l=0; l<L; l++)
		for(int m=0; m<M; m++)
		{
			if(z[m*L + l] != x[l*M + m])
			{
				std::cout << "fail" << std::endl;
				return -1;
			}
		}


// printing
#if 0
	for(int i=0; i<L*M; i++)
	{
		if(!(i%M))
			std::cout << std::endl;

		std::cout.width(4);
		std::cout << x[i] << " ";
	}

	std::cout << std::endl << std::endl;

	for(int i=0; i<L*M; i++)
	{
		if(!(i%L))
			std::cout << std::endl;

		std::cout.width(4);
		std::cout << y[i] << " ";
	}
#endif



	delete[] x;
	delete[] y;
	delete[] z;

	std::cout << std::endl << std::endl;
	return 0;
}
