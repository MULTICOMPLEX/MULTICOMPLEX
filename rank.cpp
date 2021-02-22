#include "stdafx.hpp"


constexpr int64_t R = 6;
constexpr int64_t C = 8;

/* function for exchanging two rows of
	 a matrix */
void swap(int64_t mat[R][C], int64_t row1, int64_t row2,
	int64_t col)
{
	for (int64_t i = 0; i < col; i++)
	{
		int64_t temp = mat[row1][i];
		mat[row1][i] = mat[row2][i];
		mat[row2][i] = temp;
	}
}

// Function to display a matrix 
void display(int64_t mat[R][C], int64_t row, int64_t col);

/* function for finding rank of matrix */
int64_t rankOfMatrix(int64_t mat[R][C])
{
	int64_t rank = C;

	for (int64_t row = 0; row < rank; row++)
	{
		// Before we visit current row 'row', we make 
		// sure that mat[row][0],....mat[row][row-1] 
		// are 0. 

		// Diagonal element is not zero 
		if (mat[row][row])
		{
			for (int col = 0; col < R; col++)
			{
				if (col != row)
				{
					// This makes all entries of current 
					// column as 0 except entry 'mat[row][row]' 
					double mult = (double)mat[col][row] /
						mat[row][row];
					for (int64_t i = 0; i < rank; i++)
						mat[col][i] -= int64_t(mult * mat[row][i]);
				}
			}
		}

		// Diagonal element is already zero. Two cases 
		// arise: 
		// 1) If there is a row below it with non-zero 
		//    entry, then swap this row with that row 
		//    and process that row 
		// 2) If all elements in current column below 
		//    mat[r][row] are 0, then remvoe this column 
		//    by swapping it with last column and 
		//    reducing number of columns by 1. 
		else
		{
			bool reduce = true;

			/* Find the non-zero element in current
					column  */
			for (int64_t i = row + 1; i < R; i++)
			{
				// Swap the row with non-zero element 
				// with this row. 
				if (mat[i][row])
				{
					swap(mat, row, i, rank);
					reduce = false;
					break;
				}
			}

			// If we did not find any row with non-zero 
			// element in current columnm, then all 
			// values in this column are 0. 
			if (reduce)
			{
				// Reduce number of columns 
				rank--;

				// Copy the last column here 
				for (int i = 0; i < R; i++)
					mat[i][row] = mat[i][rank];
			}

			// Process this row again 
			row--;
		}

		// Uncomment these lines to see intermediate results 
		// display(mat, R, C); 
		// printf("\n"); 
	}
	return rank;
}


/* function for displaying the matrix */
void display(int64_t mat[R][C], int64_t row, int64_t col)
{
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
			printf("  %lld", mat[i][j]);
		printf("\n");
	}
}

// Driver program to test above functions 
void rankOfMatrix_driver()
{

	int64_t mat[][8] = {
	{1876629017, 4152359834, 1082407483, 1864343263, 2188214093, 3487775666},
	{2120614597,  865551356, 1003478452, 250111500, 2483185378, 724766492},
	{4184751862, 1804098124, 3193982441, 3950428284, 1237467605, 2094870485},
	{1218239234, 3180983871, 2109061215, 2361771010, 3026984891, 1808271471},
	{3404340681, 3644829873, 3534168138, 1636848830, 2726221497, 1333971010},
	{3871457154, 3728478625, 3718639641, 145579190, 1871660900, 1678297624},
	{2046113202, 4261365484, 1853789709, 2832455776, 742681671, 2374769264},
	{2962937355, 3493917143, 1518025885, 974884375, 1029704752, 3522578350}
	};

	//unsigned int mat[][8] = { 
	/*
	{141315883,434692195,4255249181,4096446407,2361130520,2002268097},
	{1611651047,1408724345,3740769668,2267430008,228843550,1773735796},
	{1139189892,2777050888,659720259,204513614,2169749615,1203581837},
	{2164525295,1939214263,4180449344,3535920116,2094369053,3508832048},
	{4029609629,22138426,144453320,3159060626,379781407,4002658859},
	{311846527,4136615535,3356259691,797419551,494787874,2761167000},
	{2111206118,3742249927,2396576998,428061699,1060996031,2581888523},
	{3561603515,823928132,1881194574,3113384865,3539525447,4091870465}
	};
		*/
		/*
			{3553981139, 3486132163, 562827202, 142729576, 3607851650, 1966831617},
			{334430047, 1982369087, 2784897839, 3161230150, 2248782475, 554594610},
			{1319354616, 4044902126, 4138195424, 2932248383, 2375325699, 199082585},
			{1695128821, 2284370752, 1952634381, 3419985259, 2178796422, 143139904},
			{472097483, 101485014, 4073061488, 2984401253, 3707071390, 1504626319},
			{3381099193, 4261186351, 361177105, 3339717019, 1847623729, 71790825},
			{2019909186, 2328368163, 2716867707, 133090799, 2161043247, 535492356},
			{2705101217, 2828823104, 3791216385, 3677661053, 4271839598, 2318223183}
			};

		*/


		/*
		std::uniform_int_distribution<unsigned int> dist;
		for(int r=0; r < R; r++){
			for(int c=0; c < C; c++){
		mat[r][c] = dist(rng);}}
		*/

	printf("Rank of the matrix is : %lld",
		rankOfMatrix(mat));
}