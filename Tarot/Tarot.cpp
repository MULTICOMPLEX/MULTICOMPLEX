#include "stdafx.hpp"

std::string read_string_from_file(const std::string& file_path);
template<typename T>
T TRNG();
void create_html_file();
std::vector<std::string> read_lines_from_file(const std::string& FILENAME);
void Write_a_string_to_the_end_of_a_file(const std::string& FILENAME, const std::string& str);
template<typename T>
std::string int_to_hex_string(const T& i);

double ShannonEntropy(const std::string& teststring);
double Entropy(const std::string& teststring);

#define no_console
#ifdef _MSC_VER
#ifdef no_console
#include <Windows.h>
#define App() (APIENTRY WinMain( _In_ HINSTANCE hInstance,  _In_opt_ HINSTANCE hPrevInstance,  _In_ LPSTR  lpCmdLine,  _In_ int nCmdShow)) 
#else 
#define App() (main(int argc, char** argv))
#endif
#else 
#define App() (main(int argc, char** argv))
#endif

int App()
{
  
  Py_Initialize();

  //https://github.com/Yard1/sourcerandom
  PyRun_SimpleStringStd("import wikipedia,webbrowser,os,sourcerandom");

  std::vector<std::string> vs;
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/9/90/RWS_Tarot_00_Fool.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/d/de/RWS_Tarot_01_Magician.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/8/88/RWS_Tarot_02_High_Priestess.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/d/d2/RWS_Tarot_03_Empress.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/c/c3/RWS_Tarot_04_Emperor.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/8/8d/RWS_Tarot_05_Hierophant.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/d/db/RWS_Tarot_06_Lovers.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/9/9b/RWS_Tarot_07_Chariot.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/f/f5/RWS_Tarot_08_Strength.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/4/4d/RWS_Tarot_09_Hermit.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/3/3c/RWS_Tarot_10_Wheel_of_Fortune.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/e/e0/RWS_Tarot_11_Justice.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/2/2b/RWS_Tarot_12_Hanged_Man.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/d/d7/RWS_Tarot_13_Death.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/f/f8/RWS_Tarot_14_Temperance.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/5/55/RWS_Tarot_15_Devil.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/5/53/RWS_Tarot_16_Tower.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/d/db/RWS_Tarot_17_Star.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/7/7f/RWS_Tarot_18_Moon.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/1/17/RWS_Tarot_19_Sun.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/d/dd/RWS_Tarot_20_Judgement.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/f/ff/RWS_Tarot_21_World.jpg");

  vs.push_back("https://upload.wikimedia.org/wikipedia/en/1/11/Wands01.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/0/0f/Wands02.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/f/ff/Wands03.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/a/a4/Wands04.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/9/9d/Wands05.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/3/3b/Wands06.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/e/e4/Wands07.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/6/6b/Wands08.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/e/e7/Wands09.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/0/0b/Wands10.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/6/6a/Wands11.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/1/16/Wands12.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/0/0d/Wands13.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/c/ce/Wands14.jpg");

  vs.push_back("https://upload.wikimedia.org/wikipedia/en/f/fd/Pents01.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/9/9f/Pents02.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/4/42/Pents03.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/3/35/Pents04.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/9/96/Pents05.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/a/a6/Pents06.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/6/6a/Pents07.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/4/49/Pents08.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/f/f0/Pents09.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/4/42/Pents10.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/e/ec/Pents11.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/d/d5/Pents12.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/8/88/Pents13.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/1/1c/Pents14.jpg");

  vs.push_back("https://upload.wikimedia.org/wikipedia/en/3/36/Cups01.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/f/f8/Cups02.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/7/7a/Cups03.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/3/35/Cups04.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/d/d7/Cups05.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/1/17/Cups06.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/a/ae/Cups07.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/6/60/Cups08.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/2/24/Cups09.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/8/84/Cups10.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/a/ad/Cups11.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/f/fa/Cups12.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/6/62/Cups13.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/0/04/Cups14.jpg");

  vs.push_back("https://upload.wikimedia.org/wikipedia/en/1/1a/Swords01.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/9/9e/Swords02.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/0/02/Swords03.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/b/bf/Swords04.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/2/23/Swords05.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/2/29/Swords06.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/3/34/Swords07.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/a/a7/Swords08.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/2/2f/Swords09.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/d/d4/Swords10.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/4/4c/Swords11.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/b/b0/Swords12.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/d/d4/Swords13.jpg");
  vs.push_back("https://upload.wikimedia.org/wikipedia/en/3/33/Swords14.jpg");
  
  //PyRun_SimpleStringStd("print((hex(RAND_GEN.randint(0, x))).upper())");
  
  auto x = Cartomancy<uint256_t>(uint8_t(vs.size()));
  goto_url(vs[uint8_t(x)]);//goto_wiki("Modulo_operation");

  //x = Cartomancy<uint256_t>(uint8_t(vs.size()));
  //goto_url(vs[uint8_t(x)]);

  //geo_wiki();

  //random_wiki();


  Py_Finalize();

  return EXIT_SUCCESS;
}

void PyRun_SimpleStringStd
(
  const std::string& somestring
)
{
  PyRun_SimpleStringFlags(somestring.c_str(), NULL);
}

template<typename T>
T Cartomancy(const T Range)
{
  auto x = TRNG<T>();
  return x % Range;
}

void goto_url
(
  const std::string& url
)
{
  auto str = url;
  str.append("'");
  str.insert(0, "'");
  PyRun_SimpleStringStd("webbrowser.open(" + str + ")");
}

void goto_wiki
(
  const std::string& page
)
{
  auto str = page;
  str.insert(0, "https://en.wikipedia.org/wiki/");
  str.append("'");
  str.insert(0, "'");
  PyRun_SimpleStringStd("webbrowser.open(" + str + ")");
}

void random_wiki()
{
  PyRun_SimpleStringStd("wikipage = wikipedia.random(pages=1)");
  PyRun_SimpleStringStd("wikiload = wikipedia.page(wikipage)");
  PyRun_SimpleStringStd("webbrowser.open(wikiload.url,new=2)");
}

void geo_wiki()
{
  std::ostringstream out1, out2;
  out1.precision(6);
  out2.precision(6);
  auto rlocation = get_random_location(51.924419, 4.477733, 10000);
  out1 << std::fixed << rlocation.first;
  out2 << std::fixed << rlocation.second;
  //std::cout << out1.str() << std::endl;
  //std::cout << out2.str() << std::endl;

  PyRun_SimpleStringStd("wikipage = wikipedia.geosearch(" + out1.str() + ", " + out2.str() + "\
  , title=None, results=10, radius=10000)");
  PyRun_SimpleStringStd("wikiload = wikipedia.page(wikipage)");
  PyRun_SimpleStringStd("webbrowser.open(wikiload.url,new=2)");//open in new tab
}

void open_file_in_browser
(
  const std::string& file
)
{
  auto str = file;
  str.append("'");
  str.insert(0, "'");
  
  // PyRun_SimpleStringStd("os.system(" + str + ")");
  PyRun_SimpleStringStd("webbrowser.open(" + str + ", new=2)");
}

std::pair<double, double> get_random_location
(
  const double& x0,
  const double& y0,
  const int& radius
)
{
  mxws rng;
  std::uniform_real_distribution<double> r(0, 1);

  // Convert Radius from meters to degrees.
  double rd = radius / (double)111300;

  double u = r(rng);
  std::cout << u << std::endl;

  double v = r(rng);

  double w = rd * sqrt(u);
  double t = 2 * pi * v;
  double x = w * cos(t);
  double y = w * sin(t);

  double xp = x / cos(y0);

  // Resulting point.
  return std::make_pair(xp + x0, y + y0);
}

std::string read_string_from_file(const std::string& file_path) 
{
  std::ifstream input_stream(file_path, std::ios_base::binary);

  if (input_stream.fail()) {
    std::cout << "Failed to open file : " << file_path;
  }

  std::stringstream buffer;
  buffer << input_stream.rdbuf();

  input_stream.close();

  return buffer.str();
}

template<typename T>
T TRNG()
{
  PyRun_SimpleStringStd(
    "RAND_GEN1 = sourcerandom.SourceRandom(source = sourcerandom.OnlineRandomnessSource.QRNG_ANU)"
  );

  PyRun_SimpleStringStd(
    "RAND_GEN2 = sourcerandom.SourceRandom(source = sourcerandom.OnlineRandomnessSource.RANDOM_ORG)"
  );

  static int t = 0;
  std::string tel = std::to_string(t++);

  const unsigned n = std::numeric_limits<T>::digits * 2;
  typedef number<cpp_int_backend<n, n, unsigned_magnitude, unchecked, void> > uint;

  PyRun_SimpleStringStd("x = pow(2, "+ std::to_string(n) +")");
  PyRun_SimpleStringStd("f = open('random_number_" + tel + ".txt', 'w')");
  PyRun_SimpleStringStd("rv = RAND_GEN1.randint(0, x)");
  PyRun_SimpleStringStd("print(hex(rv).upper(), file=f)");
  PyRun_SimpleStringStd("rv = RAND_GEN1.randint(0, x)");
  PyRun_SimpleStringStd("print(hex(rv).upper(), file=f)");
  PyRun_SimpleStringStd("f.close()");

  //PyRun_SimpleStringStd("exec(open(\"server.py\").read())");
  //http-server -s
  //http_server localhost 8080 ./

  auto str = read_lines_from_file("./random_number_" + tel + ".txt");

  uint a(str[0]);//std::stoull(str[0], nullptr, 0);
  uint b(str[1]);//std::stoull(str[1], nullptr, 0);

  //std::seed_seq seq = { uint32_t(a >> 32), uint32_t(a) };
  //mxws rng(seq);

  mxws_t<T> rng;

  rng.w = a;
  rng.x = b;
  //std::cout << std::hex << std::uppercase << "0X" << rng.w << std::endl;
  //std::cout << std::hex << std::uppercase << "0X" << rng.x << std::endl;

  T x;

  for (int n = 0; n < 100; n++) //shuffle QRNG_ANU, RANDOM_ORG
    x = rng();

  auto rngs = int_to_hex_string(x);
  
  Write_a_string_to_the_end_of_a_file("./random_number_" + tel + ".txt", rngs);

  //create_html_file();

  //open_file_in_browser("http://localhost:8080/random_number_" + tel + ".html");
  open_file_in_browser("http://localhost:8080/random_number_0.html");
  
  return x;
}

void create_html_file()
{
  std::ofstream myfile;

  static int t = 0;
  std::string tel = std::to_string(t++);
  std::string fs = "random_number_" + tel + ".html";

  myfile.open(fs);

  std::string str_html_1 =
"<!DOCTYPE html>\n\
<html lang = \"en - US\">\n\
<head>\n\
<title>Random.html</title>\n\
<meta charset = \"UTF - 8\">\n\
<meta name = \"description\" content=\"Free Web tutorials\">\n\
<meta name = \"keywords\" content = \"HTML, CSS, JavaScript\">\n\
<meta name = \"author\" content = \"MULTICOMPLEX\">\n\
<style>\n\
div.a {\n\
line-height: normal;\n\
font-family: \"Consolas MS\";\n\
font-size: 22px;\n\
}\n\
h1 {\n\
font-family: \"Consolas MS\";\n\
font-size: 35px;\n\
}\n\
body {background-color: rgb(0, 43, 54);\n\
color: rgb(238, 232, 213);}\n\
</style>\n\
</head>\n\
<body>\n\
\n";
  
  myfile << str_html_1;

  auto str_html_2 = read_lines_from_file("./random_number_" + tel + ".txt");

  auto ds = str_html_2[0];
  auto sh = ShannonEntropy(ds.erase(0, 2));
  myfile << "<h1>Random number from QRNG_ANU  &emsp; &emsp;("+ std::to_string(sh) +") :</h1>" << std::endl;

  myfile << "<div class = \"a\">" << std::endl;
  
  myfile << str_html_2[0] << "<br>" << std::endl;
  
  myfile << "</div>" << std::endl;

  ds = str_html_2[1];
  sh = ShannonEntropy(ds.erase(0, 2));
  myfile << "<h1>Random number from RANDOM_ORG &ensp;("+ std::to_string(sh) +") :</h1>" << std::endl;

  myfile << "<div class = \"a\">" << std::endl;

  myfile << str_html_2[1] << "<br>" << std::endl;

  myfile << "</div>" << std::endl;

  ds = str_html_2[2];
  sh = ShannonEntropy(ds.erase(0, 2));
  myfile << "<h1>Shuffled QRNG_ANU + RANDOM_ORG ("+ std::to_string(sh) +") :</h1>" << std::endl;

  myfile << "<div class = \"a\">" << std::endl;

  myfile << str_html_2[2] << "<br>" << std::endl;

  myfile << "</div>" << std::endl;

  myfile << "<br>" << std::endl;
  

 // myfile << "<object data = \"http://www.shannonentropy.netmark.pl\" width=\"650\" height=\"400\">" << std::endl;
  //myfile << "</object>" << std::endl;
  
  //myfile << "<iframe src = \"https://www.wolframalpha.com/\"  width=\"950\" height=\"400\">" << std::endl;
  //myfile << "</iframe>" << std::endl;
  

  myfile << "<iframe src = \"./webgl/webgl_postprocessing_procedural.html\"  width=\"950\" height=\"510\">" << std::endl;
  myfile << "</iframe>" << std::endl;

  myfile << "<br>" << std::endl;
  myfile << "<br>" << std::endl;

  //myfile << "<iframe src = \"./webgl/clearing-with-colors.html\"  width=\"1000\" height=\"460\">" << std::endl;
  //myfile << "</iframe>" << std::endl;

  myfile << "<iframe src = \"./webgl/webgl_loader_ply.html\"  width=\"660\" height=\"510\">" << std::endl;
  myfile << "</iframe>" << std::endl;

  //myfile << "<iframe src = \"./webgl/webgl_loader_gltf.html\"  width=\"950\" height=\"510\">" << std::endl;
  //myfile << "</iframe>" << std::endl;

  myfile << "<iframe src = \"./webgl/webgl_materials_standard.html\"  width=\"950\" height=\"510\">" << std::endl;
  myfile << "</iframe>" << std::endl;

  myfile << "<br>" << std::endl;
  myfile << "<br>" << std::endl;
  
  myfile << "<iframe src = \"./webgl/webgl_loader_obj_mtl.html\"  width=\"660\" height=\"510\">" << std::endl;
  myfile << "</iframe>" << std::endl;

  //myfile << "<iframe src = \".//webgl/webgl-examples-gh-pages/tutorial/sample8/\"  width=\"660\" height=\"510\">" << std::endl;
  //myfile << "</iframe>" << std::endl;

  
  //myfile << "<iframe src = \"./Sage_Math01.html\"  width=\"660\" height=\"510\">" << std::endl;
  //myfile << "</iframe>" << std::endl;

  std::string str_html_3 =
    "</body>\n\
</html>";

  myfile << str_html_3;

  myfile.close();
}

std::vector<std::string> read_lines_from_file(const std::string& FILENAME)
{
  std::ifstream file(FILENAME);
  std::vector<std::string> linev;
  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      linev.push_back(line);
      // using printf() in all tests for consistency
      //printf("%s", line.c_str());
    }
    file.close();
  }
  else std::cout << "Can't open file :" << FILENAME << std::endl;

  return linev;
}

void Write_a_string_to_the_end_of_a_file
(
  const std::string& FILENAME, 
  const std::string& str
)
{
  std::ofstream out;
  // std::ios::app is the open mode "append" meaning
   // new data will be written to the end of the file.
  out.open(FILENAME, std::ios::app);
  
  if (out.is_open()) {
    out << str;
  }
  else std::cout << "Can't open file :" << FILENAME << std::endl;

  out.close();
}

template<typename T>
std::string int_to_hex_string(const T& i)
{
  std::stringstream sstream;
  sstream << "0X" << std::uppercase << std::hex << i;
  return sstream.str();
}

double Entropy
(
  const std::string& teststring
) 
{
  std::map<char, size_t> frequencies;
  for (auto& c : teststring)
    frequencies[c] ++;
  size_t numlen = teststring.length();
  double infocontent = 0;
  for (auto& p : frequencies) {
    auto freq = (double)(p.second) / numlen;
    infocontent -= freq * log(freq);
  }

  return infocontent;
}

double ShannonEntropy
(
  const std::string& teststring
)
{
  std::map<char, int> frequencies;
  for (auto& c : teststring)
    frequencies[c] ++;
  size_t numlen = teststring.length();
  double infocontent = 0;
  for (auto& p : frequencies) {
    auto freq = (double)(p.second) / numlen;
    infocontent -= freq * log2(freq);
  }

  return infocontent;
}

