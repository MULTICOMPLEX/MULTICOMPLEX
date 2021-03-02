
#include <vector>
#include <iostream>
#include <array>

template <typename elem, int order>
class multicomplex;


template <class T>
class Matrix
{
public:
	std::vector<std::vector<T> > array;
	size_t height;
	size_t width;

	Matrix(size_t height, size_t width);
	Matrix(std::vector<std::vector<T> > const& array);
	Matrix();
	virtual ~Matrix() = default;

	size_t getHeight() const;
	size_t getWidth() const;

	Matrix identity();

	Matrix add(const Matrix& m) const;
	Matrix subtract(const Matrix& m) const;

	Matrix dot(const Matrix<T>& m) const;
	Matrix transpose() const;
	Matrix multiply(const T& value) const;
	Matrix divide(const T& value) const;

	Matrix applyFunction(T(*function)(T)) const;
	Matrix subMat(size_t startH, size_t startW, size_t h, size_t w) const;

	void fill(const T& value);
	void put(size_t h, size_t w, const T& value);
	T get(size_t h, size_t w) const;

	void print(std::ostream&) const;

	bool operator==(const Matrix& m);
	bool operator!=(const Matrix& m);
	Matrix operator+=(const Matrix& m);
	Matrix operator-=(const Matrix& m);
	Matrix operator*=(const Matrix& m);
	Matrix operator*=(const T& m);

	Matrix operator/=(const T& m);
	Matrix operator/(const T& m);

	T& operator()(int y, int x);

	Matrix operator* (const Matrix& b);

	void transpose(Matrix& c, Matrix& d, int n, T determinant);
	void cofactor(Matrix& a, Matrix& d, int n, T determinant);
	void inverse(Matrix& a, Matrix& d, int n, T determinant);

	void minor(Matrix& b, Matrix& a, int i, int n);

	class row
	{
		Matrix& _a;
		int _i;
	public:
		row(Matrix& a, int i) : _a(a), _i(i) {}
		T operator[](int j) { return _a.array[_i][j]; }
	};

	row operator[](int i)
	{
		return row(*this, i);
	}

	T tr();
};


//	calculate minor of matrix OR build new matrix : k-had = minor
template <class T>
void Matrix<T>::minor
(
	Matrix& b,
	Matrix& a,
	int i,
	int n
)
{
	int j, l, h = 0, k = 0;
	for (l = 1; l < n; l++)
		for (j = 0; j < n; j++) {
			if (j == i)
				continue;
			b.array[h][k] = a[l][j];
			k++;
			if (k == (n - 1)) {
				h++;
				k = 0;
			}
		}
}// end function

//---------------------------------------------------

//	calculate determinant of matrix
template <class T>
T det
(
	Matrix<T>& a,
	int n
)
{
	int i;
	auto s = a.height;
	Matrix<T> b(s, s);
	T sum = 0;
	if (n == 1)
		return a[0][0];
	else if (n == 2)
		return (a[0][0] * a[1][1] - a[0][1] * a[1][0]);
	else
		for (i = 0; i < n; i++) {
			a.minor(b, a, i, n);	// read function
			sum = (T)(sum + a[0][i] * std::pow(-1, i) * det(b, (n - 1)));	// read function	// sum = determinte matrix
		}
	return sum;
}// end function

//---------------------------------------------------

//	calculate transpose of matrix
template <class T>
void Matrix<T>::transpose
(
	Matrix& c,
	Matrix& d,
	int n,
	T determinant
)
{
	int i, j;
	auto s = c.height;
	Matrix b(s, s);
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			b.array[i][j] = c[j][i];
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			d.array[i][j] = b[i][j] / determinant;	// array d[][] = inverse matrix
}// end function

//---------------------------------------------------

//	calculate cofactor of matrix
template <class T>
void Matrix<T>::cofactor
(
	Matrix& a,
	Matrix& d,
	int n,
	T determinant
)
{
	auto t = a.height;
	Matrix b(t, t), c(t, t);
	int l, h, m, k, i, j;
	for (h = 0; h < n; h++)
		for (l = 0; l < n; l++) {
			m = 0;
			k = 0;
			for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
					if (i != h && j != l) {
						b.array[m][k] = a[i][j];
						if (k < (n - 2))
							k++;
						else {
							k = 0;
							m++;
						}
					}
			c.array[h][l] = (T)std::pow(-1, (h + l)) * det(b, (n - 1));	// c = cofactor Matrix
		}
	transpose(c, d, n, determinant);	// read function
}// end function

//---------------------------------------------------

//	calculate inverse of matrix
template <class T>
void Matrix<T>::inverse
(
	Matrix& a,
	Matrix& d,
	int n,
	T determinant
)
{
	if (determinant == 0)
	{
		std::cout << "\nInverse of Entered Matrix is not possible\n";
		exit(1);
	}
	else if (n == 1)
		d.array[0][0] = 1;
	else
		cofactor(a, d, n, determinant); // read function
}// end function

//---------------------------------------------------

template <class T>
Matrix<T>::Matrix
(
	size_t height,
	size_t width
)
{
	this->height = height;
	this->width = width;
	this->array = std::vector<std::vector<T> >(height, std::vector<T>(width));
}

//---------------------------------------------------

template <class T>
Matrix<T>::Matrix
(
	std::vector<std::vector<T> > const& array
)
{
	if (array.size() == 0)
		throw std::invalid_argument("Size of array must be greater than 0.");

	this->height = array.size();
	this->width = array[0].size();
	this->array = array;
}

//---------------------------------------------------

template <class T>
Matrix<T>::Matrix
(
)
{
	height = 0;
	width = 0;
}
//---------------------------------------------------

template <class T>
size_t Matrix<T>::getHeight
(
) const
{
	return height;
}

//---------------------------------------------------

template <class T>
size_t Matrix<T>::getWidth
(
) const
{
	return width;
}

//---------------------------------------------------

template <class T>
void Matrix<T>::fill
(
	const T& value
)
{
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			array[i][j] = value;
		}
	}
}

//---------------------------------------------------

template <class T>
T Matrix<T>::tr
(
)
{
	T sum = 0;
	for (size_t i = 0; i < height; i++) {
		for (size_t j = 0; j < width; j++) {
			if (i == j)sum += this->array[i][j];
		}
	}
	return sum;
}

//---------------------------------------------------

template <class T>
void Matrix<T>::put
(
	size_t h,
	size_t w,
	const T& value
)
{
	if (!(h >= 0 && h < height && w >= 0 && w < width))
		throw std::invalid_argument("Index out of bounds.");

	array[h][w] = value;
}

//---------------------------------------------------

template <class T>
T Matrix<T>::get
(
	size_t h,
	size_t w
) const
{
	if (!(h >= 0 && h < height && w >= 0 && w < width))
		throw std::invalid_argument("Index out of bounds.");

	return array[h][w];
}

//---------------------------------------------------

template <class T>
Matrix<T> Matrix<T>::multiply
(
	const T& value
) const
{
	Matrix result(array);
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			result.array[i][j] *= value;
		}
	}

	return result;
}

//---------------------------------------------------

template <class T>
Matrix<T> Matrix<T>::divide
(
	const T& value
) const
{
	Matrix result(array);
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			result.array[i][j] /= value;
		}
	}

	return result;
}

//---------------------------------------------------

template <class T>
Matrix<T> Matrix<T>::add
(
	const Matrix& m
) const
{
	if (!(height == m.height && width == m.width))
		throw std::invalid_argument("Matrix dimension must be the same.");

	Matrix result(height, width);
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			result.array[i][j] = array[i][j] + m.array[i][j];
		}
	}

	return result;
}

//---------------------------------------------------

template <class T>
Matrix<T> Matrix<T>::subtract
(
	const Matrix& m
) const
{
	if (!(height == m.height && width == m.width))
		throw std::invalid_argument("Matrix dimension must be the same.");

	Matrix result(height, width);
	for (size_t i = 0; i < height; i++) {
		for (size_t j = 0; j < width; j++) {
			result.array[i][j] = array[i][j] - m.array[i][j];
		}
	}
	return result;
}

//---------------------------------------------------

template <class T>
Matrix<T> Matrix<T>::dot
(
	const Matrix& m
) const
{
	if (!(width == m.height))
		throw std::invalid_argument("Dot product not compatible.");

	T w = 0;
	int mwidth = m.width;

	Matrix result(height, mwidth);
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < mwidth; j++) {
			for (int h = 0; h < width; h++) {
				w += array[i][h] * m.array[h][j];
			}
			result.array[i][j] = w;
			w = 0;
		}
	}
	return result;
}

//---------------------------------------------------

template <class T>
Matrix<T> Matrix<T>::transpose
(
) const
{
	Matrix result(width, height);

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			result.array[i][j] = array[j][i];
		}
	}
	return result;
}

//---------------------------------------------------

template <class T>
Matrix<T> Matrix<T>::applyFunction
(
	T(*function)(T)
) const
{
	Matrix result(height, width);

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			result.array[i][j] = (*function)(array[i][j]);
		}
	}
	return result;
}

//---------------------------------------------------

template <class T>
Matrix<T> Matrix<T>::identity
(
)
{
	for (size_t i = 0; i < height; i++)
	{
		for (size_t j = 0; j < width; j++)
		{
			if (i == j)array[i][j] = T(1);
		}
	}
	return *this;
}

//---------------------------------------------------

template <class T>
Matrix<T> Matrix<T>::subMat
(
	size_t startH,
	size_t startW,
	size_t h,
	size_t w
) const
{
	if (!(startH >= 0 && startH + h <= height && startW >= 0 && startW + w <= width))
		throw std::invalid_argument("Index out of bounds");

	Matrix result(h, w);
	for (int i = startH; i < startH + h; i++)
	{
		for (int j = startW; j < startW + w; j++)
		{
			result.array[i - startH][j - startW] = array[i][j];
		}
	}
	return result;
}

//---------------------------------------------------

template <class T>
bool Matrix<T>::operator==
(
	const Matrix& m
	)
{
	if (height == m.height && width == m.width)
	{
		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				if (array[i][j] != m.array[i][j])
				{
					return false;
				}
			}
		}
		return true;
	}
	return false;
}

//---------------------------------------------------

template <class T>
bool Matrix<T>::operator!=
(
	const Matrix& m
	)
{
	return !operator==(m);
}

//---------------------------------------------------

template <class T>
Matrix<T> Matrix<T>::operator+=
(
	const Matrix& m
	)
{
	this->array = add(m).array;
	return *this;
}

//---------------------------------------------------

template <class T>
Matrix<T> Matrix<T>::operator-=
(
	const Matrix& m
	)
{
	this->array = subtract(m).array;
	return *this;
}

//---------------------------------------------------

template <class T>
Matrix<T> Matrix<T>::operator*=
(
	const Matrix& m
	)
{
	this->array = multiply(m).array;
	return *this;
}

//---------------------------------------------------

template <class T>
Matrix<T> Matrix<T>::operator*=
(
	const T& m
	)
{
	*this = this->multiply(m);
	return *this;
}

//---------------------------------------------------

template <class T>
Matrix<T> Matrix<T>::operator/=
(
	const T& m
	)
{
	*this = this->divide(m);
	return *this;
}

//---------------------------------------------------

template <class T>
Matrix<T> Matrix<T>::operator/
(
	const T& m
	)
{
	return *this->divide(m);
}

//---------------------------------------------------

template <class T>
T& Matrix<T>::operator()
(
	int y,
	int x
	)
{
	if (!(y >= 0 && y < height && x >= 0 && x < width))
		throw std::invalid_argument("Index out of bounds.");
	return array[y][x];
}

//---------------------------------------------------

template <class T>
Matrix<T> operator+
(
	const Matrix<T>& a,
	const Matrix<T>& b
	)
{
	return a.add(b);
}

//---------------------------------------------------

template <class T>
Matrix<T> operator-
(
	const Matrix<T>& a,
	const Matrix<T>& b
	)
{
	return a.subtract(b);
}

//---------------------------------------------------

template <class elem>
inline const Matrix<elem> operator -
(
	const elem& a,
	const Matrix<elem>& b
	)
{
	const Matrix<elem> n(a);
	return n - b;
}

//---------------------------------------------------

template <class elem>
inline const Matrix<elem> operator -
(
	const Matrix<elem>& a,
	const elem& b
	)
{
	const Matrix<elem> n(b);
	return a - n;
}

//---------------------------------------------------

template <class elem>
Matrix<elem> Matrix<elem>::operator*
(
	const Matrix<elem>& b
	)
{
	size_t h = height;
	size_t w = width;

	Matrix<elem> result(h, w);

	for (size_t i = 0; i < h; i++)
	{
		for (size_t j = 0; j < w; j++)
		{
			result.array[i][j] = 0;
			for (size_t k = 0; k < w; k++)
				result.array[i][j] += this->array[i][k] *
				b.array[k][j];
		}
	}

	return result;
}

//---------------------------------------------------

template <typename T, typename elem>
inline const Matrix<elem> operator*
(
	const T& a,
	const Matrix<elem>& b
	)
{
	auto h = b.height;
	auto w = b.width;

	Matrix<elem> result(b.array);

	for (size_t i = 0; i < h; i++)
	{
		for (size_t j = 0; j < w; j++)
		{
			result.array[i][j] *= elem(a);
		}
	}

	return result;
}

//---------------------------------------------------

template <typename elem, int order>
inline const Matrix<multicomplex<elem, order>> operator*
(
	const multicomplex<elem, order>& a,
	const Matrix<multicomplex<elem, order>>& b
	)
{
	auto h = b.height;
	auto w = b.width;

	Matrix<multicomplex<elem, order>> result(b.array);

	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			result.array[i][j] *= a;
		}
	}

	return result;
}


template <typename elem>
inline const Matrix<multicomplex<elem, 0>> operator*
(
	const multicomplex<elem, 0>& a,
	const Matrix<multicomplex<elem, 0>>& b
	)
{
	auto h = b.height;
	auto w = b.width;

	Matrix<multicomplex<elem, 0>> result(b.array);

	for (size_t i = 0; i < h; i++)
	{
		for (size_t j = 0; j < w; j++)
		{
			result.array[i][j] *= a;
		}
	}

	return result;
}

//---------------------------------------------------

template <typename T, typename elem>
inline const Matrix<elem> operator*
(
	const Matrix<elem>& b,
	const T& a
	)
{
	auto h = b.height;
	auto w = b.width;

	Matrix<elem> result(b.array);

	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			result.array[i][j] *= elem(a);
		}
	}

	return result;
}

//---------------------------------------------------

template <class elem>
inline const std::vector<elem> operator*
(
	const Matrix<elem>& a,
	const std::vector<elem>& b
	)
{
	std::size_t i, j;
	std::size_t m = a.array.size();
	std::size_t n = b.size();

	std::vector<elem> prod(m);

	for (i = 0; i < m; i++) {
		prod[i] = 0.;
		for (j = 0; j < n; j++)
			prod[i] += a.array[i][j] * b[j];
	}
	return prod;
}

//---------------------------------------------------

template <class elem>
inline const std::vector<elem> operator*
(
	const std::vector<elem>& a,
	const Matrix<elem>& b
	)
{
	std::size_t i, j;
	std::size_t m = b.array.size();
	std::size_t n = a.size();

	std::vector<elem> prod(m);

	for (i = 0; i < m; i++) {
		prod[i] = 0.;
		for (j = 0; j < n; j++)
			prod[i] += b.array[i][j] * a[j];
	}
	return prod;
}

//---------------------------------------------------

template <typename T, typename elem>
inline const std::vector<elem> operator*
(
	const std::vector<elem>& a,
	const T& b
	)
{
	std::size_t i;

	std::vector<elem> res = a;

	for (i = 0; i < a.size(); i++)
		res[i] *= elem(b);

	return res;
}

//---------------------------------------------------

template <typename T, typename elem>
inline const std::vector<elem> operator*
(
	const T& a,
	const std::vector<elem>& b
	)
{
	std::size_t i;

	std::vector<elem> res = b;

	for (i = 0; i < b.size(); i++)
		res[i] *= elem(a);

	return res;
}

//---------------------------------------------------

template <class elem>
inline const std::vector<elem> operator/
(
	const std::vector<elem>& a,
	const std::vector<elem>& b
	)
{
	std::size_t j;

	std::vector<elem> prod = a;

	for (j = 0; j < a.size(); j++)
		prod[j] /= elem(b[j]);

	return prod;
}

//---------------------------------------------------

template <class T>
void Matrix<T>::print
(
	std::ostream& flux
) const
{
	for (size_t i = 0; i < height; i++)
	{
		for (size_t j = 0; j < width; j++)
		{
			flux << array[i][j] << " ";
		}
		flux << std::endl;
	}
}

//---------------------------------------------------

template <class elem>
std::ostream& operator <<
(
	std::ostream& flux,
	const std::vector<elem>& a
	)
{
	flux << "\nvector ";

	for (auto& i : a)
		flux << " " << i;

	flux << std::endl;

	return flux;
}

//---------------------------------------------------

template <typename T>
std::ostream& operator <<
(
	std::ostream& flux,
	const Matrix<T>& m
	)
{
	m.print(flux);
	return flux;
}

//---------------------------------------------------

template <class elem>
std::vector<elem> normalize
(
	const std::vector<elem>& arr
)
{
	std::vector<elem> output(arr.size());
	elem mod = 0.0;

	for (std::size_t i = 0; i < arr.size(); ++i) {
		mod += arr[i] * arr[i];
	}

	elem mag = std::sqrt(mod);

	if (mag == 0) {
		throw std::logic_error("The input vector is a zero vector");
	}

	for (std::size_t i = 0; i < arr.size(); ++i) {
		output[i] = arr[i] / mag;
	}

	return output;
}

//---------------------------------------------------

template <class T>
Matrix<T> operator/
(
	const Matrix<T>& a,
	const T& b
	)
{
	return a.divide(b);
}

//---------------------------------------------------