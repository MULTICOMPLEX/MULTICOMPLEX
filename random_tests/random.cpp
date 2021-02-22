
#include "cn.hpp"

//#include <boost/math/tools/bivariate_statistics.hpp>

//https://blog.demofox.org/2019/06/01/taking-the-max-of-uniform-random-numbers/

void random_dist()
{
  //std::uniform_real_distribution<> dist(0, c_maxValue);
  const size_t c_maxValue = 255;
  const size_t c_numSamples = 1000000;
  
  mxws_64 rng;
 
  std::vector<size_t> counts(c_maxValue + 1);
  std::vector<size_t> countsAverage2(c_maxValue + 1);
  std::vector<size_t> countsAverage3(c_maxValue + 1);
  std::vector<size_t> countsAverage4(c_maxValue + 1);
  std::vector<size_t> countsMax2(c_maxValue + 1);
  std::vector<size_t> countsMax3(c_maxValue + 1);
  std::vector<size_t> countsMax4(c_maxValue + 1);
  std::vector<size_t> countsMin2(c_maxValue + 1);
  std::vector<size_t> countsMin3(c_maxValue + 1);
  std::vector<size_t> countsMin4(c_maxValue + 1);

  std::uniform_int_distribution<size_t> dist(0, c_maxValue);

    // generate data
  for (size_t index = 0; index < c_numSamples; ++index)
  {    
    size_t value1 = dist(rng);
    size_t value2 = dist(rng);
    size_t value3 = dist(rng);
    size_t value4 = dist(rng);

    counts[value1]++;

    countsAverage2[(value1 + value2) / 2]++;
    countsAverage3[(value1 + value2 + value3) / 3]++;
    countsAverage4[(value1 + value2 + value3 + value4) / 4]++;

    countsMax2[std::max(value1, value2)]++;
    countsMax3[std::max(value1, std::max(value2, value3))]++;
    countsMax4[std::max(std::max(value1, value2), std::max(value3, value4))]++;

    countsMin2[std::min(value1, value2)]++;
    countsMin3[std::min(value1, std::min(value2, value3))]++;
    countsMin4[std::min(std::min(value1, value2), std::min(value3, value4))]++;      
   }
    
  // write histograms
    FILE *file = nullptr;
    fopen_s(&file, "histograms.csv", "w+t");

    fprintf(file, "\"Index\",\"Count1\",\"Average2\",\"Average3\",\"Average4\",\"Max2\",\"Max3\",\"Max4\",\"Min2\",\"Min3\",\"Min4\",\"y=x\",\"y=x^2\",\"y=x^3\"\n");

    for (size_t index = 0; index < c_maxValue + 1; ++index)
    {
      fprintf(file, "\"%zu\",", index);
      fprintf(file, "\"%zu\",", counts[index]);
      fprintf(file, "\"%zu\",", countsAverage2[index]);
      fprintf(file, "\"%zu\",", countsAverage3[index]);
      fprintf(file, "\"%zu\",", countsAverage4[index]);
      fprintf(file, "\"%zu\",", countsMax2[index]);
      fprintf(file, "\"%zu\",", countsMax3[index]);
      fprintf(file, "\"%zu\",", countsMax4[index]);
      fprintf(file, "\"%zu\",", countsMin2[index]);
      fprintf(file, "\"%zu\",", countsMin3[index]);
      fprintf(file, "\"%zu\",", countsMin4[index]);

      // make the PDFs be the same scale as the histogram counts.
      float x = float(index) / float(c_maxValue);
      fprintf(file, "\"%zu\",", size_t(float(c_numSamples)*x * 2.0f / float(c_maxValue)));
      fprintf(file, "\"%zu\",", size_t(float(c_numSamples)*x*x * 3.0f / float(c_maxValue)));
      fprintf(file, "\"%zu\"\n", size_t(float(c_numSamples)*x*x*x * 4.0f / float(c_maxValue)));
    }
    
    fclose(file);
}

std::ofstream create_rng_file
(
  std::string s, 
  std::size_t n_samples
)
{
  std::ofstream myfile;
  myfile.open (s);

  std::string ms = 
"#==================================================================\n\
# generator True RNG\n\
#==================================================================\n\
type: d\n\
count: "+std::to_string(n_samples)+"\n\
numbit: 32\n";
  
  myfile << ms;
  return myfile;
}

void random2()
{
  mxws rng;

  std::map<double, double> hist;
  
  //std::uniform_real_distribution<> dist(-10, 11);
  //std::uniform_real_distribution<> dist(0, 2*pi);
  //std::uniform_int_distribution<> dist(0, 255);
  //std::uniform_real_distribution<> dist2(0, pi);
///    
  // values near the mean are the most likely
  // standard deviation affects the dispersion of generated values from the mean
  std::normal_distribution<> dist(0, 5);
///
  //the distribution of a sum of the squares of k independent standard normal random variables
  //std::chi_squared_distribution<> dist(5);
///
  //std::gamma_distribution<> dist(5);
///
  // if an event occurs 4 times a minute on average
  // how often is it that it occurs n times in one minute?
  //std::poisson_distribution<> dist(4);
///  
  // if particles decay once per second on average,
  // how much time, in seconds, until the next one?
  //std::exponential_distribution<> dist(1);
///
  //std::lognormal_distribution<> dist(0, 5);
///  
  double a,b,c;//,d,e;
  
  for (int n = 0; n < 10000000; ++n)
  {
    a = sin(dist(rng));
    b = sin(dist(rng));
    c = cos(dist(rng));
    //d = cos(dist(rng));
    //e = cos(dist(rng));
  
    ++hist[ floor( (a + b + c ) * 5 ) ];
  }

  for (auto& p : hist) 
  {
    // printf("%04d", int((p.first)));
    //std::cout << std::internal << std::setfill('0') << std::setw(4) << int(p.first);
    std::cout << std::internal << std::setw(3) << int(p.first)

    //std::cout << " : ";
      
    //printf("%5.2f", p.second/5000);
    << " : " << std::right << std::setw(6) << std::fixed << std::setprecision(2) << p.second/5000
     
    << " : " << std::string(int(p.second/10000), '*') << std::endl;
  }
  
  std::cout << std::endl << "Prefix sum :" << std::endl;
  
  double i = 0;
  for (auto& p : hist) 
  {
    i+=p.second;
    std::cout << std::internal << std::setw(3) << int(p.first)
    << " : " << std::right << std::setw(6) << std::fixed << std::setprecision(2) << i/200000
    << " : " << std::string(int(i/200000.), '*') << std::endl;
  }

  std::cout << std::endl;
}

void test_rng()
{
  std::map<uint64_t, uint64_t> hist;
 
  uint64_t a,c,x = 1;
  a = uint64_t(pow(2,8))/*256*/;
 
  mxws_64 rng;
  
  for (int n = 0; n < pow(2,23); ++n)
  {
    x = rng();
    c = x % a;
    ++hist[ c ];
  }
  
  for (auto& p : hist) 
  {
    std::cout << std::internal << std::setw(3) << int(p.first)

    << " : " << std::right << std::setw(6) << std::fixed << std::setprecision(2) << p.second
     
    << " : " << std::string((p.second/1000), '*') << std::endl;
  }
  
    std::cout << std::endl;
}

void weyl_sequence()
{
  std::map<uint64_t, uint64_t> hist;
 
  uint64_t a = uint64_t(pow(2, 8))/*256*/;
  
  std::random_device r;

  uint64_t x = 1, w = (uint64_t(r()) << 32) | r();

  for (uint64_t n = 0; n < uint64_t(pow(2,24)); ++n)
  { 
    x *= w;
    x = (x >> 32) | (x << 32);
    w += x;

    ++hist[ x % a ];
  }

  std::cout << "  Weyl sequence :\n\n";
  
  for (auto& p : hist) 
  {
    std::cout << std::internal << std::setw(3) << int(p.first)

    << " : " << std::right << std::setw(6) << std::fixed << std::setprecision(2) << p.second
     
    << " : " << std::string((p.second/1000), '*') << std::endl;
  }
  
    std::cout << std::endl;
}

unsigned int countSetBits(unsigned int n) 
{ 
    unsigned int count = 0; 
    while (n) { 
        count += n & 1; 
        n >>= 1; 
    } 
    return count; 
} 

template <typename T>
double normalize (T& value) {
  return value < 0
    ? -static_cast<double>(value) / std::numeric_limits<T>::min()
    :  static_cast<double>(value) / std::numeric_limits<T>::max()
    ;
}

double stddev
(
  std::vector<double> const & v
)
{
  double sum = std::accumulate(v.begin(), v.end(), 0.0);
  double mean = sum / v.size();

  double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
  double stdev = std::sqrt(sq_sum / v.size() - mean * mean);
  return stdev;
}

void mxms_to_file()
{
  size_t n_samples = 200000000;//51129048;//49999984;
  auto outfile = create_rng_file ("mxms", n_samples);
  
  uint32_t x;
  std::uniform_int_distribution<uint32_t> dist;
  mxws mxws;
  
  for(size_t t=0; t < n_samples;t++)
  {
    x = dist(mxws);
    
    outfile << x << std::endl;
  }
    
  outfile.close();
}

/*  
void ShowConsoleCursor(bool showFlag)
{
    HANDLE out = GetStdHandle(STD_OUTPUT_HANDLE);

    CONSOLE_CURSOR_INFO     cursorInfo;

    GetConsoleCursorInfo(out, &cursorInfo);
    cursorInfo.bVisible = showFlag; // set the cursor visibility
    SetConsoleCursorInfo(out, &cursorInfo);
}

void gotoxy(int x, int y)
{
  COORD coord;
  coord.X = x;
  coord.Y = y;
  SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), coord);
}
*/

int v_bool_to_int
(
  const std::vector<bool>& b
) 
{
  int i;
  i = accumulate(b.rbegin(), b.rend(), 0, [](int x, int y) { return (x << 1) + y; });
  return i;
}

void random_int_to_file(size_t n_samples)
{
  std::ofstream myfile;
  myfile.open ("test_rng");
  
  //rt2 -d 2 -g 64 -f rng003
  
  //example
  //rt2 -o -f example.input -t 10     test output

//int n_samples = 10000000; //10000000 
  std::string ms = 
"#==================================================================\n\
# generator True RNG\n\
#==================================================================\n\
type: d\n\
count: "+std::to_string(n_samples)+"\n\
numbit: 32\n";
  
  mxws rng;
  
  myfile << ms;
    
  for(size_t t=0; t < n_samples; t++)
    myfile << rng() << std::endl;

  myfile.close();
}

// Program to find correlation 
// coefficient 
  
// Utility Function to print 
// a Vector 
template <typename T>
void printVector
(
  const std::vector<T> &X
)  
{ 
    for (auto i: X)  
        std::cout << i << " "; 
      
    std::cout << std::endl; 
} 
  
// Function returns the rank vector 
// of the set of observations 
template <typename T>
std::vector<T> rankify
(
  std::vector<T> & X
)
{   
    int N = X.size(); 
  
    // Rank Vector 
    std::vector<T> Rank_X(N); 
      
    for(int i = 0; i < N; i++)  
    { 
        int r = 1, s = 1; 
          
        // Count no of smaller elements 
        // in 0 to i-1 
        for(int j = 0; j < i; j++) { 
            if (X[j] < X[i] ) r++; 
            if (X[j] == X[i] ) s++; 
        } 
      
        // Count no of smaller elements 
        // in i+1 to N-1 
        for (int j = i+1; j < N; j++) { 
            if (X[j] < X[i] ) r++; 
            if (X[j] == X[i] ) s++; 
        } 
  
        // Use Fractional Rank formula 
        // fractional_rank = r + (n-1)/2 
        Rank_X[i] = r + (s-1) * 0.5;         
    } 
      
    // Return Rank Vector 
    return Rank_X; 
} 
   
// function that returns 
// Pearson correlation coefficient.
template <typename T>
T correlationCoefficient 
(
  std::vector<T> &X, 
  std::vector<T> &Y
) 
{ 
    int n = X.size(); 
    T sum_X = 0, sum_Y = 0,  
                    sum_XY = 0; 
    T squareSum_X = 0,  
        squareSum_Y = 0; 
  
    for (int i = 0; i < n; i++) 
    { 
        // sum of elements of array X. 
        sum_X = sum_X + X[i]; 
  
        // sum of elements of array Y. 
        sum_Y = sum_Y + Y[i]; 
  
        // sum of X[i] * Y[i]. 
        sum_XY = sum_XY + X[i] * Y[i]; 
  
        // sum of square of array elements. 
        squareSum_X = squareSum_X +  
                      X[i] * X[i]; 
        squareSum_Y = squareSum_Y +  
                      Y[i] * Y[i]; 
    } 
  
    // use formula for calculating 
    // correlation coefficient. 
    T corr = (T)(n * sum_XY -  
                  sum_X * sum_Y) /  
                  sqrt((n * squareSum_X - 
                       sum_X * sum_X) *  
                       (n * squareSum_Y - 
                       sum_Y * sum_Y)); 
  
    return corr; 
} 
  
// Driver function 
template <typename T>
double Spearman_Driver
(
  std::vector<T>& x, 
  std::vector<T>& y
) 
{ 
    //X = {15,18,21, 15, 21}; 
    //Y= {25,25,27,27,27}; 
  
    // Get ranks of vector X 
    std::vector<double> X(x.begin(), x.end());
    std::vector<double> rank_x(X.begin(), X.end());
  
    // Get ranks of vector y 
    std::vector<double> Y(y.begin(), y.end());
    std::vector<double> rank_y(Y.begin(), Y.end());
   /*   
    std::cout << "Vector X" << std::endl; 
    printVector(X); 
  
    // Print rank vector of X  
    std::cout << "Rankings of X" << std::endl; 
    printVector(rank_x); 
      
    // Print Vector Y 
    std::cout << "Vector Y" << std::endl; 
    printVector(Y); 
  
    // Print rank vector of Y  
    std::cout << "Rankings of Y" << std::endl; 
    printVector(rank_y); 
  */
    // Print Spearmans coefficient 
    //std::cout << "Spearman's Rank correlation:" << std::endl; 
    
    auto cor = correlationCoefficient(rank_x, rank_y);
    std::cout << std::internal << std::right << std::setw(8) << std::setprecision(2) 
     << cor << std::endl;
     
 return cor;  
} 

/*
// Driver function 
template <typename T>
double Correlation_Driver
(
  std::vector<T>& x, 
  std::vector<T>& y
) 
{ 
    std::vector<double> X(x.begin(), x.end());
    std::vector<double> Y(y.begin(), y.end());
    
    auto cor = boost::math::tools::correlation_coefficient(X,Y);
    std::cout << std::internal << std::right << std::setw(8) << std::setprecision(2) 
     << cor << std::endl;
     
 return cor;  
} 
*/
void vertical_HISTOGRAM()
{
  int c,i,j,arr[10],height=0;
  system("clear");

  for(i=0 ; i<10 ; i++)
    arr[i]=0;

 while( ( c=getchar() ) != EOF)
  {
     if(c >= '0' || c <='9')
     ++arr[c-'0'];
     if( arr[c-'0'] > height )
      {
      height = arr[c-'0'];
      } 
  }
  printf("\n");
  for(j=height ; j>0 ; j--)     // row
  {
    printf("%2d|",j);
    for ( i=0 ; i<=9 ; i++)  // column
    {
    if( j == arr[i] )
     {
      printf(" *|");
      arr[i]--;
     }
     else
        printf("  |");
    }

    printf("\n");
}
  printf("  |");
  
  for(i=0; i<=9; i++) printf(" %d|",i);
    
    printf("\n  ------------DIGITS-------------");
    printf("\n");
} 

unsigned setKthBit(unsigned int n, unsigned int k) 
{ 
    // kth bit of n is being set by this operation 
    return ((1 << k) | n); 
} 
  
// Driver program to test above 
void setKthBit_driver() 
{ 
    unsigned int n = 10, k = 2; 
    std::cout << "Kth bit set number = "
         << setKthBit(n, k); 
} 

void Birthday_Probability()
{
   const int TRIALS = 1500000;
   short int birthdays[365];
   int successfulTrials;
   bool sharedBirthday;

  mxws rng;
  
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    
   for (int people = 2; people < 101; ++people) {

      successfulTrials = 0;
      for (int i = 0; i < TRIALS; ++i) {
         for (int j = 0; j < 365; birthdays[j++] = 0); // set days all to 0
         sharedBirthday = false;
         for (int j = 0; j < people; ++j) {
            // if the given birthday is shared (has more than one person)
            // then we have a shared birthday, stop checking
            if (++birthdays[rng() % 365] > 1){
               sharedBirthday = true;
               break;
            }
         }
         if (sharedBirthday) ++successfulTrials;
      }

      std::cout << "Probability of " << people << " people in a room sharing a birthday is \t"
           << ( double(successfulTrials) / double(TRIALS) ) << std::endl;
  
 }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "\nTime difference = " << 
  std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
  //11815[ms]
}

//38093904702297390785243708291056390518886454060947061/75091883268515350125426207425223147563269805908203125 = 

//23 people = 0.5072972343239854072254172283370325002359718452929878099019740020

inline uint64_t reversebits
(
  uint64_t x
)
{ 
  // swap odd and even bits
  x = ((x >> 1) & 0x5555555555555555) | ((x & 0x5555555555555555) << 1);
  // swap consecutive pairs
  x = ((x >> 2) & 0x3333333333333333) | ((x & 0x3333333333333333) << 2);
  // swap nibbles  
  x = ((x >> 4) & 0x0F0F0F0F0F0F0F0F) | ((x & 0x0F0F0F0F0F0F0F0F) << 4);
  // swap bytes
  x = ((x >> 8) & 0x00FF00FF00FF00FF) | ((x & 0x00FF00FF00FF00FF) << 8);
  // swap 2-byte long pairs
  x = ((x >> 16) & 0x0000FFFF0000FFFF) | ((x & 0x0000FFFF0000FFFF) << 16);
  // swap 4-byte long pairs
  x = ( x >> 32             ) | ( x               << 32);
  return x;
}

void test_Marsaglia_polar()
{
  size_t n_samples = 10000000;
  std::map<double, double> hist{};
  
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  for (size_t n = 0; n < n_samples; ++n)
  {
   // auto v1 = Box_Muller_transform(dist (rng),dist (rng));
    auto v = Box_Muller_transform_uni<float>(0, 1);
    
    ++hist[ floor( ( v[0]  ) * 8 ) ];
    ++hist[ floor( ( v[1]  ) * 8 ) ]; 
  }
  
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "\nTime difference = " << 
  std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
  
  for (auto& p : hist) 
  {
    std::cout << std::internal << std::right  << std::setw(5) << (int)p.first

    << " : " << std::right << std::setw(5) << std::fixed << std::setprecision(0) << p.second
     
    << " : " << std::string(uint64_t(round(p.second)/20000), '*') << std::endl;
  }
  
  hist.clear();
  
  begin = std::chrono::steady_clock::now();
  for (size_t n = 0; n < n_samples; ++n)
  {
    auto v = Marsaglia_polar<float>(0, 1);
    
    ++hist[ floor( ( v[0]  ) * 8 ) ];
    ++hist[ floor( ( v[1]  ) * 8 ) ]; 
  }
  
  end = std::chrono::steady_clock::now();
  std::cout << "\nTime difference = " << 
  std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
  
  for (auto& p : hist) 
  {
    std::cout << std::internal << std::right  << std::setw(5) << (int)p.first

    << " : " << std::right << std::setw(5) << std::fixed << std::setprecision(0) << p.second
     
    << " : " << std::string(uint64_t(round(p.second)/20000), '*') << std::endl;
  }

  std::cout << std::endl;
}

void test_cycle()
{
  mxws mxws;
  std::bitset<64> x(mxws.x);
  std::cout << "x = " << std::hex << x << std::dec << std::endl;
  
  for(int x=0; x < 128; x++)
  {
    std::cout << int(mxws()) << " ";
    if(x%16 == 15)std::cout << std::endl;
  }
}
