#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <fstream>
#include <vector>
#include "math.h"
#include <cmath>
#define PI 3.14159265359
using namespace std;
string errorHelp = "1) A = [1 2; 3] \nError: Dimensions problem, number of elements in each row and column must be equal. \nCorrect examples: A = [1 2; 3 4] ,  A = [1 2 3; 3 4 3] \n \n2) C = 3 * (2 + 4 +  \nError: Operation must be completed. \nCorrect example: C = 3 * (2 + 4 + 6 ) \n \n3) B = 5 / 0 \nError: Division by zero has no meaning. \nCorrect example: B = 5 / 1 \n \n4) C = [2 3; 4 5; 1 6] * [1 4; 2 4; 5 6] \nError: Multiplication can't be done, number of columns of first matrix must be equal number of rows of second matrix. \nCorrect example: C = [2 3; 4 5; 1 6] * [1 4; 2 4] \n \n5) C = [2 3; 4 5; 1 6] ^ 3    \nError: Matrix power is valid for square matrices (rows = columns) only. \nCorrect example: C = [2 3; 4 5] ^ 3   \n \n6) C = [2 3; 4 5; 1 6] ^ 3.4 \nError: Fractional power is valid for 1*1 matrices only. \nCorrect example: C = 2 ^ 3.4 \n \n7) A = [[1 2;3 4],[5;6;7]] \nError : Horizontal concatenate, dimensions of matrices being concatenated aren't consistent. \nCorrect example: A = [[1 2;3 4],[5;6]] \n \n8) A = [1 2; 3 4; [1 2 3]] \nError : Vertical concatenate, dimensions of matrices being concatenated aren't consistent. \nCorrect example: A = [1 2; 3 4; [1 2]] \n \n9) A = sqrt([1 -2; -3 4]) \nError: sqrt is not valid for negative numbers. \nCorrect example: A = sqrt([1 2; 3 4]) \n \n10) A = log([1 -2; -3 4]) \nError: log is not valid for negative numbers. \nCorrect example: A = log([1 2; 3 4])\n \n11) A = [2 8 9;8 9 7;2 6 6]\nB = [1 2 3;4 5 6]\nC = A / B\nError: to get the Inverse for a matrix it must be a squared Matrix.\nCorrect Example: B = [1 2 3;5 6 8;7 8 9] ";

double StringToDouble(const string &text)
{
	stringstream ss(text);
	double result;
	return ss >> result ? result : 0;
}
int stoi(const string &text)
{
	stringstream ss(text);
	int result;
	return ss >> result ? result : 0;
}

int ValidDimensions(string mString)
{
	int rows = 0, columns = 0, nColumnsOtherRows = 0;
	stringstream ss(mString);
	string token;
	while (getline(ss, token, ';'))
	{
		rows++;
		stringstream sn;
		sn << token;
		double in;
		nColumnsOtherRows = 0;
		while (sn >> in)
		{
			if (rows == 1)
				columns++;
			else if (rows > 1)
				nColumnsOtherRows++;
		}
		if ((nColumnsOtherRows != columns) && (nColumnsOtherRows != 0))
		{
			return 0;
		}
	}
	return 1;
}

//trigonometric functions for doubles NOT matrices
// sin functions --> sind , asind
double sind(double x)
{
	return sin((PI / 180)*x);
}
double asind(double x)
{
	return ((180.0 / PI)*asin(x));
}
//cos functions --> cosd , acosd
double cosd(double x)
{
	return cos((PI / 180)*x);
}
double acosd(double x)
{
	return ((180.0 / PI)*acos(x));
}
//tan functions --> tand , atand ,atan2 , atan2d
double tand(double x)
{
	return tan((PI / 180)*x);
}
double atand(double x)
{
	return ((180.0 / PI)*atan(x));
}
double atan2(double x, double y)
{
	return atan(x / y);
}
double atan2d(double x, double y)
{
	return ((180.0 / PI)*atan(x / y));
}
// csc functions --> csc , cscd ,acsc ,acscd ,csch ,acsch
double csc(double x)
{
	return (1.0 / sin(x));
}
double cscd(double x)
{
	return (1.0 / sind(x));
}
double acsc(double x) // error handling
{
	return asin(1.0 / x);
}
double acscd(double x)
{
	return asind(1.0 / x);
}
double csch(double x)
{
	return (1.0 / sinh(x));
}
double acsch(double x)
{
	return asinh(1.0 / x);
}
//sec functions --> sec ,secd ,asec, asecd ,sech ,asech
double sec(double x)
{
	return (1.0 / cos(x));
}
double secd(double x)
{
	return (1.0 / cosd(x));
}
double asec(double x) // error handling
{
	return acos(1.0 / x);
}
double asecd(double x)
{
	return acosd(1.0 / x);
}
double sech(double x)
{
	return (1.0 / cosh(x));
}
double asech(double x)
{
	return acosh(1.0 / x);
}
// cot functions cot , cotd , acot ,acotd ,coth ,acoth
double cot(double x)
{
	return (1.0 / tan(x));
}
double cotd(double x)
{
	return (1.0 / tand(x));
}
double acot(double x) // error handling
{
	return atan(1.0 / x);
}
double acotd(double x)
{
	return atand(1.0 / x);
}
double coth(double x)
{
	return (1.0 / tanh(x));
}
double acoth(double x)
{
	return atanh(1.0 / x);
}
// Square root of sum of squares (hypotenuse) function
double hypot(double x, double y) // error handling underflow and overflow
{
	return sqrt(abs(x)*abs(x) + abs(y)*abs(y));
}
// angle convertion functions
double deg2rad(double x)
{
	return x * (PI / 180.0);
}
double rad2deg(double x)
{
	return x * (180.0 / PI);
}

struct SElement
{
	double value = 0;
	int isFilled = 0;
};
class matrix
{
public:
	int rows;
	int columns;
	string mName, errorHandler;

	SElement** element;
	int getRows()
	{
		return rows;
	}

	int getColumns()
	{
		return columns;
	}
	void setRows(int r)
	{
		rows = r;
	}
	void setColumns(int c)
	{
		columns = c;
	}

	/*double** getElementArray() //if we make element private
	{
	return element;
	}*/


	matrix()
	{
		rows = 0;
		columns = 0;
		mName = "NULL";
		element = NULL;
	}

	//copy constructor
	matrix(const matrix &mCopied)
	{
		rows = mCopied.rows;
		columns = mCopied.columns;
		mName = mCopied.mName;

		if (rows == 0 && columns == 0)
			element = NULL;
		else
		{
			element = new SElement*[rows];
			for (int i = 0; i<rows; i++)
			{
				element[i] = new SElement[columns];
			}
			for (int i = 0; i<rows; i++)
			{
				for (int j = 0; j<columns; j++)
				{
					element[i][j].value = mCopied.element[i][j].value;
				}
			}
		}

	}

	void setName(string name)
	{
		this->mName = name;
	}

	string getName()
	{
		return mName;
	}



	/* double* getElement(int elementRow, int elementColumn)
	{
	return this->element[elementRow, elementColumn].value;
	}*/
	SElement getElement(int r, int c)
	{
		return this->element[r][c];
	}


	void initialising(string mName, string mString) // give me string wana azzabat isa
	{
		this->mName = mName;
		stringstream ss(mString);
		string token;
		while (getline(ss, token, ';'))
		{
			this->rows++;
			stringstream sn;
			sn << token;
			double in;
			while (sn >> in)
			{
				if (rows == 1)
					this->columns++;
			}
		}
		element = new SElement*[rows];
		for (int i = 0; i < rows; ++i)
			element[i] = new SElement[columns];
		int p = 0, q = 0;
		stringstream ss3(mString);
		string token1;
		while (getline(ss3, token1, ';'))
		{
			stringstream ss1(token1);
			string token2;
			while (ss1 >> token2)
			{
				element[p][q].value = StringToDouble(token2); //p rows --- q columns
				q++;
			}
			q = 0;
			p++;
		}
	}
	void update(string mName, string mString)
	{

		setRows(0);
		setColumns(0);
		for (int i = 0; i<rows; i++)
		{
			delete[] element[i];
		}
		delete[] element;
		initialising(mName, mString);

	}
	void update(string mName, int rows, int columns)
	{

		setRows(0);
		setColumns(0);
		if (element)
		{
			for (int i = 0; i< this->rows; i++)
			{
				delete[] element[i];
			}
			delete[] element;
		}
		
		initialising(mName, rows, columns);

	}
	void initialising(int rows1, int columns1)//give rows and columns wana azzabat isa
	{
		this->columns = columns1;
		this->rows = rows1;
		element = new SElement*[rows1];
		for (int i = 0; i<rows1; i++)
			element[i] = new SElement[columns1];
	}
	void initialising(string mName, int rows1, int columns1)//give rows and columns wana azzabat isa
	{
		this->mName = mName;
		this->columns = columns1;
		this->rows = rows1;
		element = new SElement*[rows1];
		for (int i = 0; i<rows1; i++)
			element[i] = new SElement[columns1];
	}

	void setElement(int elementRow, int elementColumn, double elementValue) // modify the main object's 2d array
	{
		element[elementRow][elementColumn].value = elementValue;
	}

	void getTranspose(matrix &x)
	{

		this->initialising(x.columns, x.rows);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, x.element[n][m].value);
			}
		}

	}

	void subMatrix(matrix& x, int elementRow, int elementColumn)
	{
		this->initialising((x.rows) - 1, (x.columns) - 1);
		for (int i = 0, m = 0; i<x.rows, m<(this->getRows()); i++)
		{
			if (i == elementRow) continue;
			else
			{
				for (int j = 0, n = 0; j<x.columns, n<(this->getColumns()); j++)
				{
					if (j == elementColumn) continue;
					else
					{
						this->setElement(m, n, x.element[i][j].value);
						n++;
					}
				}
				m++;
			}
		}

	}

	//--------------------------------------------------------------------
	double getDeterminant() {
		if (rows == columns) {
			int i, j, k;
			double factor;
			double temp;
			matrix a(*this);
			int counti = 0;
			int m = this->rows;
			for (i = 0; i < m - 1; i++)
			{
				/* Elementary Row Operation I */
				if (a.element[i][i].value == 0)
				{
					for (k = i; k < m; k++)
					{
						if (a.element[k][i].value != 0)
						{
							for (j = 0; j < m; j++)
							{
								temp = a.element[i][j].value;
								a.element[i][j].value = a.element[k][j].value;
								a.element[k][j].value = temp;
							}
							k = m;
						}
					}
					counti++;
				}
				/* Elementary Row Operation III */
				if (a.element[i][i].value != 0)
				{
					for (k = i + 1; k < m; k++)
					{
						factor = -1.0 * a.element[k][i].value / a.element[i][i].value;
						for (j = i; j < m; j++)
						{
							a.element[k][j].value = a.element[k][j].value + (factor * a.element[i][j].value);
						}
					}
				}
			}

			/* Display upper triangular matrix */

			temp = 1.0;


			for (i = 0; i < m; i++)
			{
				temp *= a.element[i][i].value;
			}


			if (counti % 2 == 0)
			{
				return temp;;
			}
			else
			{
				return -1 * temp;
			}
		}
		else return 0;

	}
	//--------------------------------------------------------------------------
	void getInverse(matrix &x)
	{
		matrix z;
		double detObj = x.getDeterminant();
		this->initialising(x.rows, x.columns);
		z.initialising(x.rows, x.columns);
		matrix sub;
		int nozero = 0;
		for (int i = 0; i<x.rows; i++)
		{
			for (int j = 0; j<x.columns; j++)
			{
				sub.subMatrix(x, i, j);
				double minor = sub.getDeterminant();
				if (minor == 0) {
					nozero = nozero / detObj;
					z.setElement(i, j, nozero);
				}
				else {
					if ((i + j) % 2 != 0) minor *= -1;




					minor = minor / detObj;
					z.setElement(i, j, minor);
				}
			}

		}
		this->getTranspose(z);
	}
	void inversePerElement(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				if (x.element[m][n].value != 0)
					this->setElement(m, n, 1 / (x.element[m][n].value));
				else
					this->errorHandler = "Error There's a zero element in the matrix"; //this is to handle 1/0 error me7taga tet3addel tab3an
			}
		}
	}

	void add(matrix& x, matrix& y, int old = 0)
	{

		if (old == 0)
			this->initialising(x.rows, x.columns);
		for (int i = 0; i< (this->rows); i++)
		{
			for (int j = 0; j< (this->columns); j++)
			{
				element[i][j].value = x.element[i][j].value + y.element[i][j].value;
				element[i][j].isFilled = 1;
			}
		}

	}

	void sub(matrix& x, matrix& y, int old = 0)
	{
		if (old == 0)
			this->initialising(x.rows, x.columns);
		for (int i = 0; i< (this->rows); i++)
		{
			for (int j = 0; j< (this->columns); j++)
			{
				element[i][j].value = x.element[i][j].value - y.element[i][j].value;
				element[i][j].isFilled = 1;
			}
		}

	}

	void mult(matrix &x, matrix &y, int asg = 0)//remember to handle errors of dimension
	{
		//declaration the output(returned) matrix
		if (asg == 0)
			this->initialising(x.rows, y.columns);//use constructor instead
		double** temp;
		if (asg == 1)
		{
			temp = new double*[x.rows];
			for (int w = 0; w<x.rows; w++)
				temp[w] = new double[x.columns];
			for (int w = 0; w<x.rows; w++)
				for (int q = 0; q<x.columns; q++)
					temp[w][q] = x.element[w][q].value;
		}
		if (asg == 2)
		{
			temp = new double*[y.rows];
			for (int w = 0; w<y.rows; w++)
				temp[w] = new double[y.columns];
			for (int w = 0; w<y.rows; w++)
				for (int q = 0; q<y.columns; q++)
					temp[w][q] = y.element[w][q].value;
		}

		//filling the matrix with valus
		for (int i = 0; i<rows; i++)
			for (int j = 0; j<columns; j++)
				element[i][j].value = 0;
		for (int i = 0; i < x.rows; i++)
		{
			for (int j = 0; j < y.columns; j++)
			{
				for (int k = 0; k < y.rows; k++)
				{
					if (asg == 0)
						element[i][j].value += x.element[i][k].value * y.element[k][j].value; //use setElement function instead
					else if (asg == 1)
						element[i][j].value += temp[i][k] * y.element[k][j].value;
					else if (asg == 2)
						element[i][j].value += x.element[i][k].value * temp[k][j];
				}
			}
		}

	}

	void div(matrix &x, matrix &y)
	{
		matrix inverseDenom;
		inverseDenom.getInverse(y);
		this->initialising(x.rows, inverseDenom.getColumns());
		this->mult(x, inverseDenom);


	}

	void zeroes(int r, int c)
	{
		rows = r;
		columns = c;
		initialising(r, c);
		for (int rs = 0; rs < r; rs++)
			for (int cs = 0; cs < c; cs++)
				element[rs][cs].value = 0;

	}
	void ones(int r, int c)
	{
		rows = r;
		columns = c;
		initialising(r, c);
		for (int rs = 0; rs < r; rs++)
			for (int cs = 0; cs < c; cs++)
				element[rs][cs].value = 1;

	}
	void eye(int r, int c)
	{
		rows = r;
		columns = c;
		initialising(r, c);
		for (int rs = 0; rs < r; rs++)
			for (int cs = 0; cs < c; cs++)
				if (rs == cs)element[rs][cs].value = 1;

	}
	void randM(int r, int c)
	{
		rows = r;
		columns = c;
		initialising(r, c);
		for (int rs = 0; rs < r; rs++)
			for (int cs = 0; cs < c; cs++)
				element[rs][cs].value = rand() % 10 + 1;

	}
	//trigonometric functions , any trigonometric function for a matrix is started with 'M'
	//sin functions --> sin , sind , asin , asind , sinh , asinh
	void Msin(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, sin(x.element[m][n].value));
			}
		}
	}
	void Msind(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, sind(x.element[m][n].value));
			}
		}
	}
	void Masin(matrix &x)  // error handling
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, asin(x.element[m][n].value));
			}
		}
	}
	void Masind(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (asind(x.element[m][n].value)));
			}
		}
	}
	void Msinh(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, sinh(x.element[m][n].value));
			}
		}
	}
	void Masinh(matrix &x)  // error handling
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, asinh(x.element[m][n].value));
			}
		}
	}
	//cos functions --> cos , cosd, acos , acosd ,cosh, acosh
	void Mcos(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, cos(x.element[m][n].value));
			}
		}
	}
	void Mcosd(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, cosd(x.element[m][n].value));
			}
		}
	}
	void Macos(matrix &x)  // error handling
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, acos(x.element[m][n].value));
			}
		}
	}
	void Macosd(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (acosd(x.element[m][n].value)));
			}
		}
	}
	void Mcosh(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, cosh(x.element[m][n].value));
			}
		}
	}
	void Macosh(matrix &x)  // error handling
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, acosh(x.element[m][n].value));
			}
		}
	}
	// tan functions --> tan,tand,atan,atand,tanh,atanh
	void Mtan(matrix &x) //error handling
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, tan(x.element[m][n].value));
			}
		}
	}
	void Mtand(matrix &x)//error handling
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, tand(x.element[m][n].value));
			}
		}
	}
	void Matan(matrix &x)  // error handling
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, atan(x.element[m][n].value));
			}
		}
	}
	void Matand(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (atand(x.element[m][n].value)));
			}
		}
	}
	void Mtanh(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, tanh(x.element[m][n].value));
			}
		}
	}
	void Matanh(matrix &x)  // error handling
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, atanh(x.element[m][n].value));
			}
		}
	}
	// csc functions --> csc , cscd ,acsc ,acscd ,csch ,acsch
	void Mcsc(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, csc(x.element[m][n].value));
			}
		}
	}
	void Mcscd(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, cscd(x.element[m][n].value));
			}
		}
	}
	void Macsc(matrix &x)  // error handling
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, acsc(x.element[m][n].value));
			}
		}
	}
	void Macscd(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (acscd(x.element[m][n].value)));
			}
		}
	}
	void Mcsch(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, csch(x.element[m][n].value));
			}
		}
	}
	void Macsch(matrix &x)  // error handling
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, acsch(x.element[m][n].value));
			}
		}
	}
	//sec functions --> sec ,secd ,asec, asecd ,sech ,asech
	void Msec(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, sec(x.element[m][n].value));
			}
		}
	}
	void Msecd(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, secd(x.element[m][n].value));
			}
		}
	}
	void Masec(matrix &x)  // error handling
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, asec(x.element[m][n].value));
			}
		}
	}
	void Masecd(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (asecd(x.element[m][n].value)));
			}
		}
	}
	void Msech(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, sech(x.element[m][n].value));
			}
		}
	}
	void Masech(matrix &x)  // error handling
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, asech(x.element[m][n].value));
			}
		}
	}
	// cot functions cot , cotd , acot ,acotd ,coth ,acoth
	void Mcot(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, cot(x.element[m][n].value));
			}
		}
	}
	void Mcotd(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, cotd(x.element[m][n].value));
			}
		}
	}
	void Macot(matrix &x)  // error handling
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, acot(x.element[m][n].value));
			}
		}
	}
	void Macotd(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (acotd(x.element[m][n].value)));
			}
		}
	}
	void Mcoth(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, coth(x.element[m][n].value));
			}
		}
	}
	void Macoth(matrix &x)  // error handling
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, acoth(x.element[m][n].value));
			}
		}
	}

	//end trignometric functions


	//element wise operators
	void addEL(matrix &x, double y)//EL stands for Element Wise  A+2 or A+1
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (x.element[m][n].value) + y);
			}
		}
	}
	void subEL(matrix &x, double y)//EL stands for Element Wise  A-2 or A-1
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (x.element[m][n].value) - y);
			}
		}
	}

	void multEL(matrix &x, double y)//EL stands for Element Wise  A.*2 or A.*1  same as A*2 or A*1
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (x.element[m][n].value)*y);
			}
		}
	}
	void divEL(matrix &x, double y)//EL stands for Element Wise A./2 or A./1  same as A/2 or A/1
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (x.element[m][n].value) / y);
			}
		}
	}
	void multEL(matrix &x, matrix& y)//EL stands for Element Wise A.*B  or  A.* c etc..
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (x.element[m][n].value)*(y.element[m][n].value));
			}
		}
	}
	void divEL(matrix &x, matrix& y)//EL stands for Element Wise A./B  or  A./c etc..
	{
		this->initialising(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (x.element[m][n].value) / (y.element[m][n].value));
			}
		}
	}
	//end element wise operators
	///// farouk
	void multforpower(matrix &x, matrix &y)           //multiplication for square matrices and assigning the result in the first matrix
	{

		this->initialising(x.rows, x.columns);

		for (int i = 0; i<rows; i++)
			for (int j = 0; j<columns; j++)
				element[i][j].value = 0;
		for (int i = 0; i < x.rows; i++)
		{
			for (int j = 0; j < x.columns; j++)
			{
				for (int k = 0; k < x.rows; k++)
				{

					element[i][j].value += x.element[i][k].value * y.element[k][j].value;
				}
			}
		}
		for (int i = 0; i < x.rows; i++)
			for (int j = 0; j < y.columns; j++)
				x.element[i][j].value = element[i][j].value;

	}
	void identityMatrix()                                      //identity matrix I
	{

		for (int i = 0; i<rows; i++)
			for (int j = 0; j<columns; j++)
				element[i][j].value = (i == j);
	}

	void power(matrix &x, double power)                  //matrix power
	{
		if (x.rows == x.columns)                          //must be square matrix
		{
			int mod = power * 10;                          //To check if power is fractional
			if (power == 0)                               //produce identity matrix if power = 0
			{
				this->initialising(x.rows, x.columns);
				this->identityMatrix();
			}
			else if (mod % 10 != 0)                       //To check if power is fractional
			{
				if (x.rows == 1)                          //Support fractional power for 1*1 matrix only
				{
					double y = pow(x.element[0][0].value, power);
					this->initialising(1, 1);
					element[0][0].value = y;
				}
				else
					errorHandler = "Error: Fraction power is supported in 1*1 matrix only.";
			}
			else
			{
				int intPower;
				if (power < 0)
					intPower = (int)abs(power);
				else
					intPower = (int)power;
				matrix y;                                         //Two matrices to be used in calculation
				y.initialising(x.rows, x.columns);
				matrix temp;
				temp.initialising(x.rows, x.columns);
				for (int i = 0; i<temp.rows; i++)
					for (int j = 0; j<temp.columns; j++)
						temp.element[i][j].value = x.element[i][j].value;
				y.identityMatrix();
				while (intPower > 0)                             //matrix power by exponentiation by squaring algorithm
				{
					if (intPower % 2 == 1)
					{
						this->multforpower(y, temp);
					}

					this->multforpower(temp, temp);
					intPower /= 2;


				}
				this->initialising(y.rows, y.columns);              //Get inverse if  power is negative
				if (power < 0)
				{
					this->getInverse(y);

				}
				else
				{
					for (int i = 0; i<rows; i++)
						for (int j = 0; j<columns; j++)
							element[i][j].value = y.element[i][j].value;
				}

			}
		}

		else
			errorHandler = "Error: for A^x, A must be a square matrix.";

	}
	void elementWisePower(matrix &x, double power)                  //matrix power
	{
		int mod = power * 10;                          //To check if power is fractional
		if (power == 0)                               //produce identity matrix if power = 0
		{
			this->initialising(x.rows, x.columns);
			for (int i = 0; i<rows; i++)
				for (int j = 0; j<columns; j++)
					element[i][j].value = 1;
		}
		else if (mod % 10 != 0)                       //To check if power is fractional
		{
			this->initialising(x.rows, x.columns);
			for (int i = 0; i<rows; i++)
				for (int j = 0; j<columns; j++)
					element[i][j].value = pow(x.element[i][j].value, power);
		}
		else
		{
			int intPower;
			if (power < 0)
				intPower = (int)abs(power);
			else
				intPower = (int)power;
			matrix y;                                         //Two matrices to be used in calculation
			y.initialising(x.rows, x.columns);
			for (int i = 0; i<y.rows; i++)
				for (int j = 0; j<y.columns; j++)
					y.element[i][j].value = 1;
			matrix temp;
			temp.initialising(x.rows, x.columns);
			for (int i = 0; i<temp.rows; i++)
				for (int j = 0; j<temp.columns; j++)
					temp.element[i][j].value = x.element[i][j].value;

			while (intPower > 0)                             //matrix power by exponentiation by squaring algorithm
			{
				if (intPower % 2 == 1)
				{
					for (int i = 0; i<temp.rows; i++)
						for (int j = 0; j<temp.columns; j++)
							y.element[i][j].value = y.element[i][j].value * temp.element[i][j].value;
				}

				for (int i = 0; i<temp.rows; i++)
					for (int j = 0; j<temp.columns; j++)
						temp.element[i][j].value = temp.element[i][j].value * temp.element[i][j].value;
				intPower /= 2;


			}
			this->initialising(y.rows, y.columns);              //Get inverse if  power is negative
			if (power < 0)
			{
				this->inversePerElement(y);

			}
			else
			{
				for (int i = 0; i<rows; i++)
					for (int j = 0; j<columns; j++)
						element[i][j].value = y.element[i][j].value;
			}

		}

	}
	void logMatrix(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int i = 0; i<rows; i++)
			for (int j = 0; j<columns; j++)
				element[i][j].value = log(x.element[i][j].value);
	}
	void log10Matrix(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int i = 0; i<rows; i++)
			for (int j = 0; j<columns; j++)
				element[i][j].value = log10(x.element[i][j].value);
	}
	void sqrtMatrix(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int i = 0; i<rows; i++)
			for (int j = 0; j<columns; j++)
				element[i][j].value = sqrt(x.element[i][j].value);
	}
	void expMatrix(matrix &x)
	{
		this->initialising(x.rows, x.columns);
		for (int i = 0; i<rows; i++)
			for (int j = 0; j<columns; j++)
				element[i][j].value = exp(x.element[i][j].value);
	}
	//////

	void print()
	{
		if (mName[0] != '@' && mName[0] != '&' &&
			mName[0] != '#' && mName[0] != '_'  && mName[0] != '!' && mName[0] != '%' && mName[0] != '&')
		{
			if (errorHandler == "Error There's a zero element in the matrix" || errorHandler == "Error The determinant of this matrix is eual to zero")
				cout << errorHandler;
			else
			{
				cout << endl;
				cout << mName << " = " << endl;
				for (int i = 0; i<rows; i++)
				{
					for (int j = 0; j<columns; j++)
					{
						cout << "\t" << element[i][j].value;
					}
					cout << endl;
				}
			}
			cout << endl;

		}
	}

	~matrix()
	{
		if (element)
		{
			for (int i = 0; i<rows; i++)
			{
				delete[] element[i];
			}
			delete[] element;
			element = NULL;
		}
	}

};
class mDynArr
{
	int arraySize;
	int usedSlots;
public:
	matrix *p;

	mDynArr()
	{
		p = new matrix[1];
		usedSlots = 0;
	}
	mDynArr(int s)
	{
		arraySize = s;
		if (arraySize <= 0)
			arraySize = 1;
		p = new matrix[arraySize];
		usedSlots = 0;
	}
	int getArraysize()
	{
		return arraySize;
	}
	int getUsedSlots()
	{
		return usedSlots;
	}
	matrix operator [] (int i)
	{
		return p[i];
	}
	void create(string name)
	{
		p[usedSlots].setName(name);
		refresh();
	}
	void create(string name, string mString)
	{
		p[usedSlots].initialising(name, mString);
		refresh();
	}
	void create(string name, int rows, int columns)
	{
		p[usedSlots].initialising(name, rows, columns);
		refresh();
	}
	void update(int num, string name, string mString)
	{
		p[num].update(name, mString);
	}
	void update(int num, string name, int rows, int columns)
	{
		p[num].update(name, rows, columns);
	}
	void refresh()
	{
		usedSlots++;
		if (usedSlots == arraySize)
		{

			arraySize += 10;
			matrix *t = new matrix[arraySize];
			for (int i = 0; i<usedSlots; i++)
			{
				t[i].setRows(p[i].getRows());
				t[i].setColumns(p[i].getColumns());
				t[i].setName(p[i].getName());

				t[i].element = new SElement*[t[i].getRows()];
				for (int o = 0; o<p[i].getRows(); o++)
				{
					t[i].element[o] = new SElement[t[i].getColumns()];

				}
				for (int k = 0; k<p[i].getRows(); k++)
				{
					for (int j = 0; j<p[i].getColumns(); j++)
					{

						t[i].element[k][j] = p[i].element[k][j];
					}
				}
				for (int g = 0; g<p[i].getRows(); g++)
				{
					delete[] p[i].element[g];
				}
				delete[] p[i].element;



			}

			p = t;
		}
	}

};


//matrix* memory =new matrix[20];
mDynArr memory(10);
int memoryPointer = 0;
int exit1 = 0; //exit is ambiguos to VS so i added the 1 (it's meaningless)

void removeSpaces(string &str)
{
	for (int i = 0; i<str.length(); i++)
	{
		if (str[i] == ' ') { str.erase(i, 1); i--; }//replaced erase fn with =''

	}
}


int memoryCheck(string mName)
{
	for (int i = 0; i<memoryPointer; i++)
	{
		if (memory[i].getName() == mName)
			return i;
	}
	return -1;

}
void memorizeMatrix(int mIndex, string mString, string mName)//0 for all phase1, 1 for sAssignment
{
	if (mIndex == -1)
	{
		// cout<<"mwmory"<<memoryPointer;
		memory.create(mName, mString);
		memory.p[memoryPointer].print();
		memoryPointer++;
	}

	else
	{
		memory.update(mIndex, mName, mString);
		memory.p[mIndex].print();
	}
}
void memorizeMatrix(int mIndex, int rows, int columns, string mName)//0 for all phase1, 1 for sAssignment
{
	if (mIndex == -1)
	{
		// cout<<"mwmory"<<memoryPointer;
		memory.create(mName, rows, columns);
		memoryPointer++;
	}

	else
	{
		memory.update(mIndex, mName, rows, columns);
	}
}
void cut(string &variable1, string &variable2, int &index1, int &index2, char op, string operation)
{
	variable1 = operation.substr(0, operation.find(op));
	variable2 = operation.substr(operation.find(op) + 1, (operation.length() - operation.find(op)) - 1);

	index1 = memoryCheck(variable1);
	index2 = memoryCheck(variable2);
}
void cut(string &variable1, int &index1, char op, string operation)
{
	variable1 = operation.substr(operation.find(op) + 1, (operation.length() - operation.find(op)) - 1);
	index1 = memoryCheck(variable1);
}

void removeSpaces2(string &str)
{
	while (str[str.length() - 1] == ' ')
		str.erase(str.length() - 1, 1);
}
struct sizeValue
{
	int rows;
	int columns;
};

vector <string> separatedString;

void separate(string inputString)
{
	//vector<string> separatedString;
	int flag = 0;
	int beginPostion = 0;
	for (int i = 0; i<inputString.length(); i++)
	{
		if (inputString[i] == ';'&&flag == 0)
		{
			//cout<<inputString.substr(beginPostion,i-beginPostion)<<endl;
			separatedString.push_back(inputString.substr(beginPostion, i - beginPostion));
			beginPostion = i + 1;
		}
		else if (i == inputString.length() - 1)
		{
			//cout<<inputString.substr(beginPostion,inputString.length()-beginPostion);
			separatedString.push_back(inputString.substr(beginPostion, inputString.length() - beginPostion));
		}
		else if (inputString[i] == '[')
			flag++;
		else if (inputString[i] == ']')
			flag--;
	}
	//return separatedString;
}

sizeValue compare(sizeValue m1, sizeValue m2)     //get the total size of 2 concatenated matrices
{
	sizeValue m;
	if (m1.rows == m2.rows)   //if the rows of matrix1 = rows of matrix2 whatever the columns are equal or not
	{
		m.rows = m1.rows;
		m.columns = m1.columns + m2.columns;
		return m;
	}
	else if (m1.columns == m2.columns&&m1.rows != m2.rows)  //if thr columns of matrix1= columns of matrix2 and the rows aren't equal
	{
		m.columns = m1.columns;
		m.rows = m1.rows + m2.rows;
		return m;
	}
}
//============================================================================
sizeValue get_sizeValue(string s)  //s is a string without spaces or ;   ex:   123  or  10*sin(A)
{
	sizeValue n;
	int flag = 0;
	int index;
	for (int i = 0; i<s.length(); i++)
	{
		if ((s[i]>64 && s[i]<91) || (s[i]>96 && s[i]<123))
		{
			flag = 1;
			if (s[i] == 's'&&s[i + 1] == 'i'&&s[i + 2] == 'n')    //if sin
			{
				index = memoryCheck(s.substr(i + 4, s.find(')') - i - 4));
				if (index == -1)
				{
					n.rows = 1;  n.columns = 1;
					return n;
				}
				n.columns = memory.p[index].getColumns();
				n.rows = memory.p[index].getRows();
				return n;
			}
			else if (s[i] == 'c'&&s[i + 1] == 'o'&&s[i + 2] == 's')  // if cos
			{
				index = memoryCheck(s.substr(i + 4, s.find(')') - i - 4));
				if (index == -1)
				{
					n.rows = 1;  n.columns = 1;
					return n;
				}
				n.columns = memory.p[index].getColumns();
				n.rows = memory.p[index].getRows();
				return n;
			}
			else if (s[i] == 't'&&s[i + 1] == 'a'&&s[i + 2] == 'n')   //  if tan
			{
				index = memoryCheck(s.substr(i + 4, s.find(')') - i - 4));
				if (index == -1)
				{
					n.rows = 1;  n.columns = 1;
					return n;
				}
				n.columns = memory.p[index].getColumns();
				n.rows = memory.p[index].getRows();
				return n;
			}
			else if (s[i] == 'r'&&s[i + 1] == 'a'&&s[i + 2] == 'n'&&s[i + 3] == 'd')   // if rand(55,100)
			{
				n.columns = stoi(s.substr(s.find(',') + 1, s.find(')') - s.find(',') - 1));
				n.rows = stoi(s.substr(s.find('(') + 1, s.find(',') - s.find(',') - 1));
				return n;  //return numbers of columns 100 & rows 55
			}
			else if (s[i] == 'e'&&s[i + 1] == 'y'&&s[i + 2] == 'e')     //if eye(55,100)
			{
				n.columns = stoi(s.substr(s.find(',') + 1, s.find(')') - s.find(',') - 1));
				n.rows = stoi(s.substr(s.find('(') + 1, s.find(',') - s.find(',') - 1));
				return n;  //return numbers of columns 100 & rows 55
			}
			else if (s[i] == 'z'&&s[i + 1] == 'e'&&s[i + 2] == 'r'&&s[i + 3] == 'o'&&s[i + 4] == 's')   // if zeros(55,100)
			{
				n.columns = stoi(s.substr(s.find(',') + 1, s.find(')') - s.find(',') - 1));
				n.rows = stoi(s.substr(s.find('(') + 1, s.find(',') - s.find(',') - 1));
				return n;  //return numbers of columns 100 & rows 55
			}
			else if (s[i] == 'o'&&s[i + 1] == 'n'&&s[i + 2] == 'e'&&s[i + 3] == 's')  // if ones(55,100)
			{
				n.columns = stoi(s.substr(s.find(',') + 1, s.find(')') - s.find(',') - 1));
				n.rows = stoi(s.substr(s.find('(') + 1, s.find(',') - s.find(',') - 1));
				return n;  //return numbers of columns 100 & rows 55
			}
			else
			{
				index = memoryCheck(s.substr(i, 1));
				n.columns = memory.p[index].getColumns();
				n.rows = memory.p[index].getRows();
				return n;
			}
		}
	}
	if (flag == 0)
	{
		n.rows = 1;  n.columns = 1;
		return n;
	}
}
//=============================================================================
sizeValue sizing(string matrix)
{
	sizeValue n; n.rows = 0;  n.columns = 0;
	if (matrix[0] == '[')
		matrix = matrix.substr(1, matrix.length() - 2);
	stringstream sMatrix(matrix);
	string token;
	getline(sMatrix, token, ';');
	stringstream sn(token);
	string element;
	while (sn >> element)
	{
		n.columns += get_sizeValue(element).columns;
	}
	//-----------------------------------
	stringstream ssMatrix(matrix);
	while (getline(ssMatrix, token, ';'))
	{
		stringstream sn(token);
		string element;
		sn >> element;
		n.rows += get_sizeValue(element).rows;
	}
	return n;
}
//================================================================================

sizeValue conc(string s)
{
	int i = 0, j = 0, k = 0;
	string miniMatrix;
	sizeValue Vstack[2]; //sizeValue Vstack[2]; supposed
	Vstack[1].rows = 0;
	Vstack[1].columns = 0;
	while (1)
	{
		i = s.rfind('['); //if not found returns -1
		if (i == -1)
		{
			int flag = 0;
			for (int y = 0; y<s.length(); y++)
			{
				if (s[y] != ',' && s[y] != '[' && s[y] != ']')
				{
					flag = 1;
				}
			}
			if (flag == 1)
			{
				Vstack[k] = sizing(s);
				Vstack[0] = compare(Vstack[0], Vstack[1]);
			}
			break;
		}
		for (int o = i; o<s.length(); o++)
		{
			if (s[o] == ']')
			{
				j = o;
				break;
			}
		}
		if (1)//i != -1 )//&& j != -1)
		{
			miniMatrix = s.substr(i, j - i + 1);
			Vstack[k] = sizing(miniMatrix);
			if (Vstack[1].rows != 0)
			{
				Vstack[0] = compare(Vstack[0], Vstack[1]);
				Vstack[1].rows = 0; //it can be overwritten
				Vstack[1].columns = 0;
			}
			s.erase(i, j - i + 1);
			if (s.rfind('[') == -1)
			{
				int flag = 0;
				for (int y = 0; y<s.length(); y++)
				{
					if (s[y] != ',' && s[y] != '[' && s[y] != ']')
					{
						flag = 1;
					}
				}
				if (flag == 1)
				{
					Vstack[1] = sizing(s);
					Vstack[0] = compare(Vstack[0], Vstack[1]);
				}
				break;
			}
			removeSpaces2(s);
			for (int o = i; o <= (s.length()); o++)
			{
				if ((s[o] == ']' && s[o - 1] == '[')
					|| (s[o] == ']' && s[o - 1] == ',' && s[o - 2] == '[')
					|| (s[o] == ']' && s[o - 1] == ',')
					|| (s[o] == '[' && s[o - 1] == ',')
					|| (s[o] == ',' && s[o - 1] == ']')
					|| (s[o] == ',' && s[o - 1] == '[')
					)
				{
					if ((s[o] == ']' && s[o - 1] == '['))
						s.erase(o - 1, 2);
					else if ((s[o] == ']' && s[o - 1] == ',' && s[o - 2] == '['))
						s.erase(o - 2, 3);
					else if ((s[o] == ']' && s[o - 1] == ','))
						s.erase(o - 1, 2);
					else if ((s[o] == '[' && s[o - 1] == ','))
						s.erase(o - 1, 2);
					else if ((s[o] == ',' && s[o - 1] == ']'))
						s.erase(o - 1, 2);
					else if ((s[o] == ',' && s[o - 1] == '['))
						s.erase(o - 1, 2);
				}
			}
			k = 1;
		}
	}
	return Vstack[0];
}

sizeValue calcSize(vector<string>& separatedString)
{
	vector <sizeValue> finStack;
	int k = 0;
	for (int i = 0; i<separatedString.size(); i++)
	{
		if (separatedString[i].find('[') == -1)
			finStack.push_back(sizing(separatedString[i]));
		else
			finStack.push_back(conc(separatedString[i]));
		k++;
	}
	if (finStack.size() == 1)
	{
		separatedString.clear();
		return finStack[0];
	}
	sizeValue sum = compare(finStack[0], finStack[1]);
	for (int i = 2; i<finStack.size(); i++)
	{
		sum = compare(sum, finStack[i]);
	}
	separatedString.clear();
	return sum;
}

string mul_ope_solver(string &ope);
void sFill(matrix &mSoph, string mString, string mName, matrix &trMat)
{
	int lastPos = 0;
	string temp;
	matrix tempM;
	int flag = 0, OBCounter = 0, CBFlag = 0, continued = 0;
	int OBColumnCounter[100] = { 0 }, OBRowCounter[100] = { 0 }, OBColumnReset = 0, OBRowReset = 0;
	for (int r = 0; r < mSoph.rows; r++)
	{
		for (int c = 0; c < mSoph.columns; c++)
		{
			if (mSoph.element[r + OBRowCounter[OBCounter]][c + OBColumnCounter[OBCounter]].isFilled == 1)
			{
				continued = 1;
				continue;
			}
			temp = mString.substr(lastPos, 1);
			//If character is found
			if ((int(temp[0]) >= 65 && int(temp[0]) <= 90) || (int(temp[0]) >= 97 && int(temp[0]) <= 122))
			{
				//If it's a stored matrix
				int index = memoryCheck(temp);
				int spacePos = mString.find(' ', lastPos)
					, semicolumnPos = mString.find(';', lastPos)
					, CBPos = mString.find(']', lastPos);
				if (spacePos == -1) spacePos = 999999;
				if (semicolumnPos == -1) semicolumnPos = 999999;
				if (CBPos == -1) CBPos = 999999;
				if (index != -1)
				{
					//B
					//memory.p[index]
					if (temp == mName)
					{
						for (int rs = 0; rs < trMat.rows; rs++)
						{
							for (int cs = 0; cs < trMat.columns; cs++)
							{
								mSoph.element[r + rs + OBRowCounter[OBCounter]][c + cs + OBColumnCounter[OBCounter]].value = trMat.element[rs][cs].value;
								mSoph.element[r + rs + OBRowCounter[OBCounter]][c + cs + OBColumnCounter[OBCounter]].isFilled = 1;
							}
						}
					}
					else
					{
						for (int rs = 0; rs < memory.p[index].rows; rs++)
						{
							for (int cs = 0; cs < memory.p[index].columns; cs++)
							{
								mSoph.element[r + rs + OBRowCounter[OBCounter]][c + cs + OBColumnCounter[OBCounter]].value = memory.p[index].element[rs][cs].value;
								mSoph.element[r + rs + OBRowCounter[OBCounter]][c + cs + OBColumnCounter[OBCounter]].isFilled = 1;
							}
						}
					}

					lastPos += 1;
					flag = 1;
				}

				else
				{
					temp = mString.substr(lastPos, 3);
					//if trigonometric function
					if (temp == "sin" || temp == "cos" || temp == "tan" || temp == "sec" || temp == "cot"
						|| mString.substr(lastPos, 4) == "cose")
					{
						temp = mString.substr(lastPos, mString.find(')', lastPos) - lastPos + 1);
						//didn't handle if there are brackets inside the asdsin(())
						//tempM = Trig(temp);
						string tempOPE = temp;
						tempOPE = mul_ope_solver(tempOPE);
						if (tempOPE == "Statement or Expression is incorrect" || tempOPE == "Matrix dimensions must agree")
						{
							cout << tempOPE << endl;
							delete &mSoph;
							return;
						}

						for (int rs = 0; rs < memory.p[memoryCheck(tempOPE)].rows; rs++)
						{
							for (int cs = 0; cs < memory.p[memoryCheck(tempOPE)].columns; cs++)
							{
								mSoph.element[r + rs][c + cs].value = memory.p[memoryCheck(tempOPE)].element[rs][cs].value;
								mSoph.element[r + rs][c + cs].isFilled = 1;
							}
						}
						lastPos += temp.length();
						flag = 1;
					}
					//if rand
					else if (temp == "ran")
					{
						temp = mString.substr(lastPos + 4, mString.find(')', lastPos) - lastPos + 1);
						//tempM = rand(temp)
						matrix rMatrix;
						removeSpaces(temp);
						string s1(1, temp[1]);
						string s3(1, temp[3]);
						rMatrix.randM(StringToDouble(s1), StringToDouble(s1));
						for (int rs = 0; rs < rMatrix.rows; rs++)
						{
							for (int cs = 0; cs < rMatrix.columns; cs++)
							{
								mSoph.element[r + rs][c + cs].value = rMatrix.element[rs][cs].value;
								mSoph.element[r + rs][c + cs].isFilled = 1;
							}
						}
						lastPos += temp.length();
						flag = 1;
					}
					//if eye
					else if (temp == "eye")
					{
						temp = mString.substr(lastPos + 3, mString.find(')', lastPos) - lastPos + 1);
						//tempM = eye(temp)
						matrix rMatrix;
						removeSpaces(temp);
						string s1(1, temp[1]);
						string s3(1, temp[3]);
						rMatrix.eye(StringToDouble(s1), StringToDouble(s1));
						for (int rs = 0; rs < rMatrix.rows; rs++)
						{
							for (int cs = 0; cs < rMatrix.columns; cs++)
							{
								mSoph.element[r + rs][c + cs].value = rMatrix.element[rs][cs].value;
								mSoph.element[r + rs][c + cs].isFilled = 1;
							}
						}
						lastPos += temp.length();
						flag = 1;
					}
					//if zeroes
					else if (temp == "zer")
					{
						temp = mString.substr(lastPos + 6, mString.find(')', lastPos) - lastPos + 1);
						//tempM = rand(temp)
						matrix rMatrix;
						removeSpaces(temp);
						string s1(1, temp[1]);
						string s3(1, temp[3]);
						rMatrix.zeroes(StringToDouble(s1), StringToDouble(s1));
						for (int rs = 0; rs < rMatrix.rows; rs++)
						{
							for (int cs = 0; cs < rMatrix.columns; cs++)
							{
								mSoph.element[r + rs][c + cs].value = rMatrix.element[rs][cs].value;
								mSoph.element[r + rs][c + cs].isFilled = 1;
							}
						}
						lastPos += temp.length();
						flag = 1;
					}
					//if ones
					else if (temp == "one")
					{
						temp = mString.substr(lastPos + 4, mString.find(')', lastPos) - lastPos + 1);
						//tempM = rand(temp)
						matrix rMatrix;
						removeSpaces(temp);
						string s1(1, temp[1]);
						string s3(1, temp[3]);
						rMatrix.ones(StringToDouble(s1), StringToDouble(s1));
						for (int rs = 0; rs < rMatrix.rows; rs++)
						{
							for (int cs = 0; cs < rMatrix.columns; cs++)
							{
								mSoph.element[r + rs][c + cs].value = rMatrix.element[rs][cs].value;
								mSoph.element[r + rs][c + cs].isFilled = 1;
							}
						}
						lastPos += temp.length();
						flag = 1;
					}
				}
				if ((spacePos < semicolumnPos) && (spacePos < CBPos))
				{
					if (OBCounter != 0)
					{
						OBColumnCounter[OBCounter]++;//walking through the inner matrix
						c--;//freezing c as it will increase by one the next loop
					}
				}
				else if ((semicolumnPos < spacePos) && (semicolumnPos < CBPos))
				{
					if (OBCounter != 0)
					{
						OBRowCounter[OBCounter]++;//walking through the inner matrix
						OBColumnCounter[OBCounter] = 0;//reseting column to start position
						c--;
					}
				}
				else
				{
					if (OBCounter != 0)
					{
						OBColumnCounter[OBCounter] = 0;//reseting both dimensions to the start position
						OBRowCounter[OBCounter] = 0;
					}
				}
			}
			//if [ 1 2 3]
			else if (temp[0] == '[')
			{
				OBCounter++;
			}
			else if (temp[0] == ']')
			{
				OBColumnCounter[OBCounter] = 0;
				OBRowCounter[OBCounter] = 0;
				OBCounter--;
				if (OBCounter == 0)
				{
					OBColumnCounter[OBCounter] = 0;
					OBRowCounter[OBCounter] = 0;
				}
			}
			//if 2.3
			else if (((int(temp[0]) <= 57) && (int(temp[0]) >= 48)) || temp[0] == '-')
			{

				int spacePos = mString.find(' ', lastPos)
					, semicolumnPos = mString.find(';', lastPos)
					, CBPos = mString.find(']', lastPos);
				int operatorPos = 999999
					, plusPos = mString.find('+', lastPos + 1)
					, minusPos = mString.find('-', lastPos + 1)
					, divPos = mString.find('/', lastPos + 1)
					, multPos = mString.find('*', lastPos + 1)
					, powPos = mString.find('^', lastPos + 1);
				if (plusPos == -1) plusPos = 999999;
				if (minusPos == -1) minusPos = 999999;
				if (divPos == -1) divPos = 999999;
				if (multPos == -1) multPos = 999999;
				if (powPos == -1) powPos = 999999;

				if ((plusPos < minusPos) && (plusPos < divPos) && (plusPos < multPos) && (plusPos < powPos))
					operatorPos = plusPos;
				else if ((minusPos < plusPos) && (minusPos< divPos) && (minusPos< multPos) && (minusPos< powPos))
					operatorPos = minusPos;
				else if ((divPos < plusPos) && (divPos<minusPos) && (divPos< multPos) && (divPos< powPos))
					operatorPos = divPos;
				else if ((multPos < plusPos) && (multPos<minusPos) && (multPos < divPos) && (multPos< powPos))
					operatorPos = multPos;
				else if ((powPos < plusPos) && (powPos<minusPos) && (powPos < divPos) && (powPos< multPos))
					operatorPos = powPos;

				if (spacePos == -1) spacePos = 999999;
				if (semicolumnPos == -1) semicolumnPos = 999999;
				if (CBPos == -1) CBPos = 999999;
				//Because the separators inside [] are numerous
				if ((spacePos < semicolumnPos) && (spacePos < CBPos))
				{
					temp = mString.substr(lastPos, spacePos - lastPos);
					if (operatorPos < spacePos)//if it's not just a number but an operation
					{
						string tempOPE = temp;
						tempOPE = mul_ope_solver(tempOPE);
						if (tempOPE == "Statement or Expression is incorrect" || tempOPE == "Matrix dimensions must agree")
						{
							cout << tempOPE << endl;
							delete &mSoph;
							return;
						}
						for (int rs = 0; rs < memory.p[memoryCheck(tempOPE)].rows; rs++)
						{
							for (int cs = 0; cs < memory.p[memoryCheck(tempOPE)].columns; cs++)
							{
								mSoph.element[r + rs + OBRowCounter[OBCounter]][c + cs + OBColumnCounter[OBCounter]].value = memory.p[memoryCheck(tempOPE)].element[rs][cs].value;
								mSoph.element[r + rs + OBRowCounter[OBCounter]][c + cs + OBColumnCounter[OBCounter]].isFilled = 1;
							}
						}

					}
					else//if it's a number
					{
						stringstream ss;
						ss << temp;
						double value;
						ss >> value;
						mSoph.element[r + OBRowCounter[OBCounter]][c + OBColumnCounter[OBCounter]].value = value;
						mSoph.element[r + OBRowCounter[OBCounter]][c + OBColumnCounter[OBCounter]].isFilled = 1;
					}
					//adjusting the position in the string and in the matrix
					if (OBCounter != 0)
					{
						OBColumnCounter[OBCounter]++;//walking through the inner matrix
						c--;//freezing c as it will increase by one the next loop
					}
					lastPos = spacePos + 1;
				}
				else if ((semicolumnPos < spacePos) && (semicolumnPos < CBPos))
				{
					temp = mString.substr(lastPos, semicolumnPos - lastPos);
					if (operatorPos < semicolumnPos) // if it's not just a number but an operation
					{
						string tempOPE = temp;
						tempOPE = mul_ope_solver(tempOPE);
						if (tempOPE == "Statement or Expression is incorrect" || tempOPE == "Matrix dimensions must agree")
						{
							cout << tempOPE << endl;
							delete &mSoph;
							return;
						}
						for (int rs = 0; rs < memory.p[memoryCheck(tempOPE)].rows; rs++)
						{
							for (int cs = 0; cs < memory.p[memoryCheck(tempOPE)].columns; cs++)
							{
								mSoph.element[r + rs + OBRowCounter[OBCounter]][c + cs + OBColumnCounter[OBCounter]].value = memory.p[memoryCheck(tempOPE)].element[rs][cs].value;
								mSoph.element[r + rs + OBRowCounter[OBCounter]][c + cs + OBColumnCounter[OBCounter]].isFilled = 1;
							}
						}
					}
					else//if it's a number
					{
						stringstream ss;
						ss << temp;
						double value;
						ss >> value;
						mSoph.element[r + OBRowCounter[OBCounter]][c + OBColumnCounter[OBCounter]].value = value;
						mSoph.element[r + OBRowCounter[OBCounter]][c + OBColumnCounter[OBCounter]].isFilled = 1;
					}
					//adjusting the position in the string and in the matrix
					if (OBCounter != 0)
					{
						OBRowCounter[OBCounter]++;//walking through the inner matrix
						OBColumnCounter[OBCounter] = 0;//reseting column to start position
						c--;
					}

					lastPos = semicolumnPos + 1;
				}
				else
				{
					temp = mString.substr(lastPos, CBPos - lastPos);
					if (operatorPos < CBPos)// if it's not just a number but an operation
					{
						string tempOPE = temp;
						tempOPE = mul_ope_solver(tempOPE);
						if (tempOPE == "Statement or Expression is incorrect" || tempOPE == "Matrix dimensions must agree")
						{
							cout << tempOPE << endl;
							delete &mSoph;
							return;
						}
						for (int rs = 0; rs < memory.p[memoryCheck(tempOPE)].rows; rs++)
						{
							for (int cs = 0; cs < memory.p[memoryCheck(tempOPE)].columns; cs++)
							{
								mSoph.element[r + rs + OBRowCounter[OBCounter]][c + cs + OBColumnCounter[OBCounter]].value = memory.p[memoryCheck(tempOPE)].element[rs][cs].value;
								mSoph.element[r + rs + OBRowCounter[OBCounter]][c + cs + OBColumnCounter[OBCounter]].isFilled = 1;
							}
						}
					}
					else//if it's a number
					{
						stringstream ss;
						ss << temp;
						double value;
						ss >> value;
						mSoph.element[r + OBRowCounter[OBCounter]][c + OBColumnCounter[OBCounter]].value = value;
						mSoph.element[r + OBRowCounter[OBCounter]][c + OBColumnCounter[OBCounter]].isFilled = 1;
					}
					//adjusting the position in the string and in the matrix
					if (OBCounter != 0)
					{
						OBColumnCounter[OBCounter] = 0;//reseting both dimensions to the start position
						OBRowCounter[OBCounter] = 0;
					}

					lastPos = CBPos;//I want to start from the ]
				}

				flag = 1;
			}
			if (flag == 0)//if the character im on now is rubish ( space , etc...)
			{
				lastPos += 1;
				c--;
			}
			flag = 0;
			continued = 0;
		}
	}
	mSoph.print();
}
void input_checker(string input) // assignment or operation
{

	short asg = 0;
	//FOB=First Open Bracket, FCB=First Close Bracker, EQ=EQual
	//stop calling find numerous times and use the below for loop
	short FOBPos = input.find('['), SOBPos = input.find('[', FOBPos + 1), FCBPos = input.find(']'), LCBPos = input.rfind(']'), EQPos = input.find('='), plusPos = input.find('+')
		, minusPos = input.find('-'), elemWiseInvPos = input.find("./"), divPos = input.find('/'), multPos = input.find("*")
		, transPos = input.find("'");
	string mName = input.substr(0, EQPos);
	removeSpaces(mName);
	bool assignmentOP = (FOBPos != -1) && (FCBPos != -1);
	bool sAssignmentOP = (SOBPos != -1) && (FOBPos != -1); //sophisticated assignment
	bool mathOP = (plusPos != -1) || (minusPos != -1) || (divPos != -1) || (multPos != -1) || (transPos != -1);

	/*	for (int i = 0; i < input.length(); i++)
	{

	if ((int(input[i])>=65&&int(input[i])<=90) || (int(input[i])>=97&&int(input[i])<=122))//if input contains letters
	{
	sAssignmentOP = true;
	}
	}*/

	int index = memoryCheck(mName);
	if (sAssignmentOP)
	{
		string mString = input.substr(FOBPos + 1, (LCBPos - 2 - FOBPos + 1));// Entering without []
																			 //call concatenation function to calculate the rows and the columns of the matrix and to create it
		separate(mString);
		sizeValue mSize = calcSize(separatedString);
		matrix trickyMat;
		if (index != -1)
		{
			matrix tempM;
			tempM.zeroes(memory.p[index].rows, memory.p[index].columns);
			trickyMat.add(memory.p[index], tempM);
		}
		memorizeMatrix(index, mSize.rows, mSize.columns, mName);
		index = memoryCheck(mName);
		sFill(memory.p[index], mString, mName, trickyMat);
		// memory.p[index]; this will be the matrix im working on
		//memory.p[index].sFill(mString); sophisticated filling
	}
	else if (assignmentOP)
	{
		string mString = input.substr(FOBPos + 1, (LCBPos - 2 - FOBPos + 1));//FOB+1 & LCB-2 to remove braces
		if (!ValidDimensions(mString))
		{
			cout << "Dimensions of matrices being concatenated are not consistent." << endl;
			return;
		}
		else memorizeMatrix(index, mString, mName);
	}

	else if (mathOP)
	{
		string mString = input.substr(EQPos + 1, (input.length() - EQPos + 1));
		string tempOPE = mString;
		matrix zMatrix;
		tempOPE = mul_ope_solver(tempOPE);
		if (tempOPE == "Statement or Expression is incorrect" || tempOPE == "Matrix dimensions must agree")
		{
			cout << tempOPE << endl;
			return;
		}
		memory.create(mName);
		memoryPointer++;
		zMatrix.zeroes(memory.p[memoryCheck(tempOPE)].rows, memory.p[memoryCheck(tempOPE)].columns);
		memory.p[memoryCheck(mName)].add(memory.p[memoryCheck(tempOPE)], zMatrix);
		memory.p[memoryCheck(mName)].print();
		//3shan te3mel copy e3mel add(m,0) , we 3shan tet2akked en law 3amalt call 2 times mayeb2ash 3andak 2 funcs L
		//me7tageen nerga3 nezbot 7ewar el exit we enena lamma nekteb esm matrix ye3mellaha print
	}
	else
	{
		string mString = input.substr(EQPos + 1, (input.length() - EQPos + 1));
		removeSpaces(input);
		removeSpaces(mString);
		string tempTrig = mString.substr(0, 3);
		if ((memoryCheck(input) != -1) &&
			memoryCheck(input)<memoryPointer)
		{
			memory.p[memoryCheck(input)].print();
		}
		else if (tempTrig == "sin" || tempTrig == "cos" || tempTrig == "tan" || tempTrig == "sec" || tempTrig == "cot")
		{
			//basy le yasser
		}
		else if (tempTrig == "ran")
		{
			memory.create(mName);
			memoryPointer++;
			string s1(1, mString[5]);
			string s3(1, mString[7]);
			memory.p[memoryCheck(mName)].randM(StringToDouble(s1), StringToDouble(s1));
			memory.p[memoryCheck(mName)].print();
		}
		else if (tempTrig == "eye")
		{
			memory.create(mName);
			memoryPointer++;
			string s1(1, mString[4]);
			string s3(1, mString[6]);
			memory.p[memoryCheck(mName)].eye(StringToDouble(s1), StringToDouble(s3));
			memory.p[memoryCheck(mName)].print();
		}
		else if (tempTrig == "zer")
		{
			memory.create(mName);
			memoryPointer++;
			string s1(1, mString[6]);
			string s3(1, mString[8]);
			memory.p[memoryCheck(mName)].zeroes(StringToDouble(s1), StringToDouble(s3));
			memory.p[memoryCheck(mName)].print();
		}
		else if (tempTrig == "one")
		{
			memory.create(mName);
			memoryPointer++;
			string s1(1, mString[5]);
			string s3(1, mString[7]);
			memory.p[memoryCheck(mName)].ones(StringToDouble(s1), StringToDouble(s3));
			memory.p[memoryCheck(mName)].print();
		}
		else if (input == "exit")
		{
			exit1 = 1;
			return;
		}
		else
		{
			bool notNumber = false;
			if (mString[mString.length() - 1] == ';') mString = mString.substr(0, mString.length() - 1);
			for (int i = 0; i < mString.length(); i++)
				if (!(mString[i] >= '0' && mString[i] <= '9') && (mString[i] != '.'))
				{
					notNumber = true;
					break;
				}
			if (notNumber)
			{
				cout << "invalid input" << endl;
			}
			else
			{
				memorizeMatrix(index, mString, mName);
			}
		}

	}
}
bool Is_operation(char character) {

	if (character == '(' || character == ')' || character == '^' || character == '*' ||

		character == '/' || character == '+' || character == '-' || character == '~')

		return 1;

	else return 0;

}
string alphanum =
"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
string symb = "@#$_!&%";
int stringLength = alphanum.length();
int Mctr = 0;
string genRandom()  // Random string generator function.
{
	if (Mctr >= stringLength)
	{
		symb.erase(0, 1);
		alphanum =
			"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
		Mctr = 0;
	}
	string ranS = "_a";
	ranS[1] = alphanum[0];
	ranS[0] = symb[0];
	alphanum.erase(0, 1);
	Mctr++;
	return ranS;
}
bool calcAndRep(int i, int j, string  &fullOp, char ch_op) {
	//cout<<i<<" "<<j << " "<<fullOp<<" "<<ch_op<<endl ;

	int opOnNum = 1;
	int pos1 = i, pos2 = j;
	double op1 = 0, op2 = 0, res = 0;

	string s_op1, s_op2;

	stringstream strm, strm2;
	stringstream strm3;
	string result;
	string op = fullOp.substr(pos1, pos2 - pos1 + 1);


	s_op1 = op.substr(0, op.find(ch_op, 1));
	s_op2 = op.substr(op.find(ch_op, 1) + 1, pos2 - op.find(ch_op, 1) + 1);
	if (ch_op == '~')
		s_op2.erase(0, 2);

	//cout<<op<<" opers "<<s_op1<<" "<<s_op2<<endl;
	if (op.find("sqrt") == 0)
	{

		s_op1 = op.substr(4, op.length() - 4);
		if (memoryCheck(s_op1) != -1)
		{
			for (int i = 0; i<memory.p[memoryCheck(s_op1)].getRows(); i++)
			{
				for (int j = 0; j<memory.p[memoryCheck(s_op1)].getColumns(); j++)
				{
					if (memory.p[memoryCheck(s_op1)].element[i][j].value <= 0)
						return 0;
				}
			}
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].sqrtMatrix(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else {
			strm << s_op1;
			strm >> op1;
			if (op1 <= 0)
				return 0;
			res = sqrt(op1);
		}

		// cout<<s_op1<<endl;
	}
	else if (op.find("~P~")<op.length())
	{
      if (memoryCheck(s_op2) != -1 && memoryCheck(s_op1) != -1)
        {
            string nMat = genRandom();
			memory.create(nMat);

			memory.p[memoryPointer].add(memory.p[memoryCheck(s_op2)], memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
        }
		else if (memoryCheck(s_op2) != -1)
		{

			strm << s_op1;
			strm >> op1;
			string nMat = genRandom();
			memory.create(nMat);

			memory.p[memoryPointer].addEL(memory.p[memoryCheck(s_op2)], op1);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else if (memoryCheck(s_op1) != -1)
		{
			strm << s_op2;
			strm >> op2;
			string nMat = genRandom();
			memory.create(nMat);

			memory.p[memoryPointer].addEL(memory.p[memoryCheck(s_op1)], op2);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
	}
	else if (op.find("~S~")<op.length())
	{
		if (memoryCheck(s_op2) != -1 && memoryCheck(s_op1) != -1)
        {
            string nMat = genRandom();
			memory.create(nMat);

			memory.p[memoryPointer].sub(memory.p[memoryCheck(s_op1)], memory.p[memoryCheck(s_op2)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
        }
		else if (memoryCheck(s_op2) != -1)
		{

			strm << s_op1;
			strm >> op1;
			string nMat = genRandom();
			memory.create(nMat);

			memory.p[memoryPointer].subEL(memory.p[memoryCheck(s_op2)], op1);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else if (memoryCheck(s_op1) != -1)
		{
			strm << s_op2;
			strm >> op2;
			string nMat = genRandom();
			memory.create(nMat);

			memory.p[memoryPointer].subEL(memory.p[memoryCheck(s_op1)], op2);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
	}
	else if (op.find("~D~")<op.length())
	{
		if ((memoryCheck(s_op1) != -1) && memoryCheck(s_op2) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);

			memory.p[memoryPointer].divEL(memory.p[memoryCheck(s_op1)],
				memory.p[memoryCheck(s_op2)]);

			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else if (memoryCheck(s_op2) != -1)
		{
			strm << s_op1;
			strm >> op1;
			string nMat = genRandom();
			memory.create(nMat);

			memory.p[memoryPointer].divEL(memory.p[memoryCheck(s_op2)], op1);
			memory.p[memoryPointer].elementWisePower(memory.p[memoryPointer], -1);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else if (memoryCheck(s_op1) != -1)
		{
			strm << s_op2;
			strm >> op2;
			string nMat = genRandom();
			memory.create(nMat);

			memory.p[memoryPointer].divEL(memory.p[memoryCheck(s_op1)], op2);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
	}
	else if (op.find("~M~")<op.length())
	{
		if ((memoryCheck(s_op1) != -1) && memoryCheck(s_op2) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);

			memory.p[memoryPointer].multEL(memory.p[memoryCheck(s_op1)],
				memory.p[memoryCheck(s_op2)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else if (memoryCheck(s_op2) != -1)
		{

			strm << s_op1;
			strm >> op1;
			string nMat = genRandom();
			memory.create(nMat);

			memory.p[memoryPointer].multEL(memory.p[memoryCheck(s_op2)], op1);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else if (memoryCheck(s_op1) != -1)
		{
			strm << s_op2;
			strm >> op2;
			string nMat = genRandom();
			memory.create(nMat);

			memory.p[memoryPointer].multEL(memory.p[memoryCheck(s_op1)], op2);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
	}
	else if (op.find("~W~")<op.length())
	{

		if (memoryCheck(s_op1) != -1)
		{
			strm << s_op2;
			strm >> op2;
			string nMat = genRandom();
			memory.create(nMat);

			memory.p[memoryPointer].elementWisePower(memory.p[memoryCheck(s_op1)], op2);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
	}
	else if (op.find("exp") == 0)
	{
		s_op1 = op.substr(3, op.length() - 3);
		if (memoryCheck(s_op1) != -1)
		{

			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].expMatrix(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		{
			strm << s_op1;
			strm >> op1;
			res = exp(op1);
		}
	}
	else if (op.find("log") == 0)
	{
		s_op1 = op.substr(3, op.length() - 3);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].logMatrix(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		{
			strm << s_op1;
			strm >> op1;
			res = log(op1);
		}
	}
	else if (op.find("tog10") == 0)
	{
		s_op1 = op.substr(5, op.length() - 5);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].log10Matrix(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		{
			strm << s_op1;
			strm >> op1;
			res = log10(op1);
		}
	}
	else if (op.find("sinh") == 0)
	{
		s_op1 = op.substr(4, op.length() - 4);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Msinh(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		strm << s_op1;
		strm >> op1;
		res = sinh(op1);
	}
	else if (op.find("sind") == 0)
	{
		s_op1 = op.substr(4, op.length() - 4);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Msind(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		strm << s_op1;
		strm >> op1;
		res = sind(op1);
	}
	else if (op.find("sin") == 0)
	{
		s_op1 = op.substr(3, op.length() - 3);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Msin(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		{
			strm << s_op1;
			strm >> op1;
			res = sin(op1);
		}
	}

	else if (op.find("asind") == 0)
	{
		s_op1 = op.substr(5, op.length() - 5);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Masind(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = sin(op1);
		}
	}
	else if (op.find("asinh") == 0)
	{
		s_op1 = op.substr(5, op.length() - 5);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Masinh(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = sin(op1);
		}
	}
	else if (op.find("asin") == 0)
	{

		s_op1 = op.substr(4, op.length() - 4);
		if (memoryCheck(s_op1) != -1)
		{

			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Masin(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = sin(op1);
		}
	}

	else if (op.find("cosd") == 0)
	{
		s_op1 = op.substr(4, op.length() - 4);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Mcosd(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = sin(op1);
		}
	}
	else if (op.find("cosh") == 0)
	{
		s_op1 = op.substr(4, op.length() - 4);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Mcosh(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = sin(op1);
		}
	}
	else if (op.find("cos") == 0)
	{
		s_op1 = op.substr(3, op.length() - 3);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Mcos(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = sin(op1);
		}
	}

	else if (op.find("acosd") == 0)
	{
		s_op1 = op.substr(5, op.length() - 5);
		strm << s_op1;
		strm >> op1;
		res = acosd(op1);
		// cout<<s_op1<<endl;
	}
	else if (op.find("acosh") == 0)
	{
		s_op1 = op.substr(5, op.length() - 5);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Macosh(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = sin(op1);
		}
	}
	else if (op.find("acos") == 0)
	{
		s_op1 = op.substr(4, op.length() - 4);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Macos(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = sin(op1);
		}
	}

	else if (op.find("tand") == 0)
	{
		s_op1 = op.substr(4, op.length() - 4);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Mtand(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = sin(op1);
		}
	}
	else if (op.find("tanh") == 0)
	{
		s_op1 = op.substr(4, op.length() - 4);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Mtanh(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = sin(op1);
		}
	}
	else if (op.find("tan") == 0)
	{
		s_op1 = op.substr(3, op.length() - 3);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Mtan(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = sin(op1);
		}
	}

	else if (op.find("atand") == 0)
	{
		s_op1 = op.substr(5, op.length() - 5);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Matand(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = sin(op1);
		}
	}
	else if (op.find("atanh") == 0)
	{
		s_op1 = op.substr(5, op.length() - 5);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Matanh(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = sin(op1);
		}
	}

	else if (op.find("atan2d") == 0) //need some work
	{
		s_op1 = op.substr(6, op.length() - 6);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Matand(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = atanh(op1);
		}
	}
	else if (op.find("atan2") == 0) //need some work
	{
		s_op1 = op.substr(5, op.length() - 5);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Matanh(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = atanh(op1);
		}
	}
	else if (op.find("atan") == 0)
	{
		s_op1 = op.substr(4, op.length() - 4);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Matan(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = atan(op1);
		}
	}

	else if (op.find("cscd") == 0)
	{
		s_op1 = op.substr(4, op.length() - 4);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Mcscd(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = cscd(op1);
		}
	}
	else if (op.find("csch") == 0)
	{
		s_op1 = op.substr(4, op.length() - 4);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Mcsch(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = csch(op1);
		}
	}
	else if (op.find("csc") == 0)
	{
		s_op1 = op.substr(3, op.length() - 3);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Mcsc(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = csc(op1);
		}
	}

	else if (op.find("acscd") == 0)
	{
		s_op1 = op.substr(5, op.length() - 5);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Macscd(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = acscd(op1);
		}
	}
	else if (op.find("acsch") == 0)
	{
		s_op1 = op.substr(5, op.length() - 5);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Macsch(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = acsch(op1);
		}
	}
	else if (op.find("acsc") == 0)
	{
		s_op1 = op.substr(4, op.length() - 4);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Macsc(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = acsc(op1);
		}
	}

	else if (op.find("secd") == 0)
	{
		s_op1 = op.substr(4, op.length() - 4);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Msecd(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = secd(op1);
		}
	}
	else if (op.find("sech") == 0)
	{
		s_op1 = op.substr(4, op.length() - 4);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Msech(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = sech(op1);
		}
	}
	else if (op.find("sec") == 0)
	{
		s_op1 = op.substr(3, op.length() - 3);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Msec(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = sec(op1);
		}
	}

	else if (op.find("asecd") == 0)
	{
		s_op1 = op.substr(5, op.length() - 5);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Masecd(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = asecd(op1);
		}
	}
	else if (op.find("asech") == 0)
	{
		s_op1 = op.substr(5, op.length() - 5);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Masech(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = asech(op1);
		}
	}
	else if (op.find("asec") == 0)
	{
		s_op1 = op.substr(4, op.length() - 4);
		if (memoryCheck(s_op1) != -1)
		{
			string nMat = genRandom();
			memory.create(nMat);
			memory.p[memoryPointer].Masec(memory.p[memoryCheck(s_op1)]);
			opOnNum = 0;
			memoryPointer++;
			result = nMat;
		}
		else
		{
			strm << s_op1;
			strm >> op1;
			res = asec(op1);
			// cout<<"GLLLLLLL"<<" "<<s_op1<< " "<< res<<endl;
		}
	}
	else if (op.find("hypot") == 0) //need some work
	{
		s_op1 = op.substr(5, op.length() - 5);
		strm << s_op1;
		strm >> op1;
		res = atanh(op1);
		//cout<<s_op1<<endl;
	}
	else if (op.find("deg2rad") == 0) //need some work
	{
		s_op1 = op.substr(5, op.length() - 5);
		strm << s_op1;
		strm >> op1;
		res = atanh(op1);
		// cout<<s_op1<<endl;
	}
	else if (op.find("rad2deg") == 0) //need some work
	{
		s_op1 = op.substr(5, op.length() - 5);
		strm << s_op1;
		strm >> op1;
		res = atanh(op1);
		// cout<<s_op1<<endl;
	}
	else
	{
		strm << op;
		strm >> op1;
		op = op.substr(op.find(ch_op, 1), op.length() - op.find(ch_op, 1));
		ch_op = op[0];
		op.erase(0, 1);
		strm2 << op;
		strm2 >> op2;
		//cout << op1 <<"/-/ "<<op2 <<endl<<op<<endl;

		switch (ch_op)
		{
		case '^':


			if (memoryCheck(s_op1) != -1 && memoryCheck(s_op2) != -1)
			{
				if (memory.p[memoryCheck(s_op2)].getColumns() != 1
					|| memory.p[memoryCheck(s_op2)].getRows() != 1)
					return 0;
				if (memory.p[memoryCheck(s_op1)].getColumns() != memory.p[memoryCheck(s_op1)].getRows())
					return 0;

				if (memory.p[memoryCheck(s_op2)].getColumns() == 1 &&
					memory.p[memoryCheck(s_op2)].getRows() == 1)
				{
					string nMat = genRandom();
					memory.create(nMat);
					memory.p[memoryPointer].power(memory.p[memoryCheck(s_op1)],
						memory.p[memoryCheck(s_op2)].getElement(0, 0).value);
					opOnNum = 0;
					memoryPointer++;
					result = nMat;

				}
				else break;
			}
			else if (memoryCheck(s_op1) != -1)
			{
				if (memory.p[memoryCheck(s_op1)].getColumns() != memory.p[memoryCheck(s_op1)].getRows())
					return 0;
				if (op2<1 || op2 - (int)op2 != 0)
					return 0;
				string nMat = genRandom();
				memory.create(nMat);
				memory.p[memoryPointer].power(memory.p[memoryCheck(s_op1)], op2);
				opOnNum = 0;
				memoryPointer++;
				result = nMat;
			}
			res = pow(op1, op2);
			break;
		case '+':
			if (memoryCheck(s_op1) != -1 && memoryCheck(s_op2) != -1)
			{
				if (memory.p[memoryCheck(s_op1)].getColumns() != memory.p[memoryCheck(s_op2)].getColumns()
					|| memory.p[memoryCheck(s_op1)].getRows() != memory.p[memoryCheck(s_op2)].getRows())
					return 0;

				string nMat = genRandom();
				memory.create(nMat);
				memory.p[memoryPointer].add(memory.p[memoryCheck(s_op1)],
					memory.p[memoryCheck(s_op2)]);
				opOnNum = 0;
				memoryPointer++;
				result = nMat;
			}
			else if (memoryCheck(s_op1) != -1 )
            {
                string nMat = genRandom();
				memory.create(nMat);
				memory.p[memoryPointer].addEL(memory.p[memoryCheck(s_op1)],
					op2);
				opOnNum = 0;
				memoryPointer++;
				result = nMat;
            }
            else if (memoryCheck(s_op2) != -1 )
            {

                string nMat = genRandom();
				memory.create(nMat);
				memory.p[memoryPointer].addEL(memory.p[memoryCheck(s_op2)],
					op1);
				opOnNum = 0;
				memoryPointer++;
				result = nMat;
            }
			else
				res = op1 + op2;
			break;
		case '-':
			if (memoryCheck(s_op1) != -1 && memoryCheck(s_op2) != -1)
			{
				if (memory.p[memoryCheck(s_op1)].getColumns() != memory.p[memoryCheck(s_op2)].getColumns()
					|| memory.p[memoryCheck(s_op1)].getRows() != memory.p[memoryCheck(s_op2)].getRows())
					return 0;

				string nMat = genRandom();
				memory.create(nMat);
				memory.p[memoryPointer].sub(memory.p[memoryCheck(s_op1)],
					memory.p[memoryCheck(s_op2)]);
				opOnNum = 0;
				memoryPointer++;
				result = nMat;
			}
			else if (memoryCheck(s_op1) != -1 )
            {
                string nMat = genRandom();
				memory.create(nMat);
				memory.p[memoryPointer].subEL(memory.p[memoryCheck(s_op1)],
					op2);
				opOnNum = 0;
				memoryPointer++;
				result = nMat;
            }
            else if (memoryCheck(s_op2) != -1 )
            {

                string nMat = genRandom();
				memory.create(nMat);
				memory.p[memoryPointer].subEL(memory.p[memoryCheck(s_op2)],
					op1);
				opOnNum = 0;
				memoryPointer++;
				result = nMat;
            }
			else
				res = op1 - op2;
			break;
		case '*':
			if (memoryCheck(s_op1) != -1 && memoryCheck(s_op2) != -1)
			{
				if (memory.p[memoryCheck(s_op1)].getColumns() != memory.p[memoryCheck(s_op2)].getRows()
					|| memory.p[memoryCheck(s_op2)].getColumns() != memory.p[memoryCheck(s_op1)].getRows())
					return 0;
				string nMat = genRandom();
				memory.create(nMat);
				memory.p[memoryPointer].mult(memory.p[memoryCheck(s_op1)],
					memory.p[memoryCheck(s_op2)]);
				opOnNum = 0;
				memoryPointer++;
				result = nMat;
			}
			else if (memoryCheck(s_op1) != -1 )
            {
                string nMat = genRandom();
				memory.create(nMat);
				memory.p[memoryPointer].multEL(memory.p[memoryCheck(s_op1)],
					op2);
				opOnNum = 0;
				memoryPointer++;
				result = nMat;
            }
            else if (memoryCheck(s_op2) != -1 )
            {

                string nMat = genRandom();
				memory.create(nMat);
				memory.p[memoryPointer].multEL(memory.p[memoryCheck(s_op2)],
					op1);
				opOnNum = 0;
				memoryPointer++;
				result = nMat;
            }
			else
				res = op1 * op2;
			break;
		case '/':
			if (memoryCheck(s_op1) != -1 && memoryCheck(s_op2) != -1)
			{
				if (memory.p[memoryCheck(s_op1)].getColumns() != memory.p[memoryCheck(s_op2)].getRows()
					|| memory.p[memoryCheck(s_op2)].getColumns() != memory.p[memoryCheck(s_op1)].getRows())
					return 0;
				string nMat = genRandom();
				memory.create(nMat);
				memory.p[memoryPointer].div(memory.p[memoryCheck(s_op1)],
					memory.p[memoryCheck(s_op2)]);
				opOnNum = 0;
				memoryPointer++;
				result = nMat;
			}
			else if (memoryCheck(s_op1) != -1 )
            {
                string nMat = genRandom();
				memory.create(nMat);
				memory.p[memoryPointer].divEL(memory.p[memoryCheck(s_op1)],op2);

				opOnNum = 0;
				memoryPointer++;
				result = nMat;
            }
            else if (memoryCheck(s_op2) != -1 )
            {
                if(memory.p[memoryCheck(s_op2)].getColumns()!=1 &&
                   memory.p[memoryCheck(s_op2)].getRows()!=1)
                    return 0;
                string nMat = genRandom();
				memory.create(nMat);
				memory.p[memoryPointer].divEL(memory.p[memoryCheck(s_op2)],op1);
				memory.p[memoryPointer].elementWisePower(memory.p[memoryPointer], -1);
				opOnNum = 0;
				memoryPointer++;
				result = nMat;
            }
			else
			{
				if (op2 == 0)
					return 0;
				res = op1 / op2;
			}

			break;
		}
	}
	//if((abs(res-round(res))<=0.00001 && abs(res-round(res))>=-0.00001)||
	//     (abs(round(res)-res )<=0.00001 && abs(round(res) -res)>=-0.00001) )
	//    res=round(res);
	if (opOnNum)
	{
		strm3 << res;
		strm3 >> result;
	}

	fullOp.replace(pos1, pos2 - pos1 + 1, result);
	//cout<<fullOp<<endl;
	return 1;
}

void Operation_solver(string &operation) {
	if (operation == "Statement or Expression is incorrect" || operation == "Matrix dimensions must agree")
		return;
	removeSpaces(operation);
	bool error = 1;
	if (Is_operation(operation[operation.length() - 1]) && operation[operation.length() - 1] != ')') {
		operation = "Statement or Expression is incorrect";
		return;
	}
	if (operation[operation.length() - 1] == '*' && operation[operation.length() - 1] == '.') {
		operation = "Statement or Expression is incorrect";
		return;
	}
	if (operation.find("log10") != -1) {
		int rep = operation.find("log10");
		operation[rep] = 't';
	}
	if (operation.find(".*") != -1) {
		int rep = operation.find(".*");
		operation.replace(rep, 2, "~M~");

	}
	if (operation.find("./") != -1) {
		int rep = operation.find("./");
		operation.replace(rep, 2, "~D~");
	}
	if (operation.find(".+") != -1) {
		int rep = operation.find(".+");
		operation.replace(rep, 2, "~P~");
	}
	if (operation.find(".-") != -1) {
		int rep = operation.find(".-");
		operation.replace(rep, 2, "~S~");
	}
	if (operation.find(".^") != -1) {
		int rep = operation.find(".^");
		operation.replace(rep, 2, "~W~");
	}

	for (int i = operation.length() - 1; i >= 0; i--) {
		if (operation[i] == '(') {
			int j = 0;
			int found_bracket = 0;

			for (j = i; j<operation.length(); j++)
				if (operation[j] == ')') {
					found_bracket = 1;
					break;
				}
			if (found_bracket == 0) {
				operation = "Statement or Expression is incorrect";
				return;

			}



			string temp = operation.substr(i + 1, j - i - 1);
			//cout<<"Temp "<<temp<<endl;
			Operation_solver(temp);
			operation.replace(i, j - i + 1, temp);
			Operation_solver(operation);

			break;
		}
	}

	for (int i = operation.length() - 1; i >= 0; i--) { //sqrt
		if (operation[i] == 's'&&operation[i + 1] == 'q'&&operation[i + 2] == 'r'
			&&operation[i + 3] == 't') {
			int end = 0;
			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			error = calcAndRep(i, end, operation, operation[i]);
			if (error == 0) {
				operation = "Statement or expression is incorrect";
				return;
			}

			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}

	for (int i = operation.length() - 1; i >= 0; i--) { //asind
		if (operation[i] == 'a'&&operation[i + 1] == 's'&&operation[i + 2] == 'i'
			&&operation[i + 3] == 'n'&&operation[i + 4] == 'd') {
			int end = 0;
			for (end = i + 5; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			error = calcAndRep(i, end, operation, operation[i]);
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //asinh
		if (operation[i] == 'a'&&operation[i + 1] == 's'&&operation[i + 2] == 'i'
			&&operation[i + 3] == 'n'&&operation[i + 4] == 'h') {
			int end = 0;
			for (end = i + 5; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			error = calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //asin
		if (operation[i] == 'a'&&operation[i + 1] == 's'&&operation[i + 2] == 'i'
			&&operation[i + 3] == 'n') {
			int end = 0;
			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;

			error = calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //sinh
		if (operation[i] == 's'&&operation[i + 1] == 'i'&&operation[i + 2] == 'n'
			&&operation[i + 3] == 'h') {
			int end = 0;
			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			error = calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //sind
		if (operation[i] == 's'&&operation[i + 1] == 'i'&&operation[i + 2] == 'n'
			&&operation[i + 3] == 'd') {
			int end = 0;
			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			error = calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //sin
		if (operation[i] == 's'&&operation[i + 1] == 'i'&&operation[i + 2] == 'n') {
			int end = 0;
			for (end = i + 3; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			error = calcAndRep(i, end, operation, operation[i]);
			//  cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}

	for (int i = operation.length() - 1; i >= 0; i--) { //acosd
		if (operation[i] == 'a'&&operation[i + 1] == 'c'&&operation[i + 2] == 'o'
			&&operation[i + 3] == 's'&&operation[i + 4] == 'd') {
			int end = 0;
			for (end = i + 5; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			error = calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //acosh
		if (operation[i] == 'a'&&operation[i + 1] == 'c'&&operation[i + 2] == 'o'
			&&operation[i + 3] == 's'&&operation[i + 4] == 'h') {
			int end = 0;
			for (end = i + 5; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) {//acos
		if (operation[i] == 'a'&&operation[i + 1] == 'c'&&operation[i + 2] == 'o'
			&&operation[i + 3] == 's') {
			int end = 0;
			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //cosh
		if (operation[i] == 'c'&&operation[i + 1] == 'o'&&operation[i + 2] == 's'
			&&operation[i + 3] == 'h') {
			int end = 0;
			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //cosd
		if (operation[i] == 'c'&&operation[i + 1] == 'o'&&operation[i + 2] == 's'
			&&operation[i + 3] == 'd') {
			int end = 0;
			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //cos
		if (operation[i] == 'c'&&operation[i + 1] == 'o'&&operation[i + 2] == 's') {
			int end = 0;
			for (end = i + 3; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}



	for (int i = operation.length() - 1; i >= 0; i--) { //atand
		if (operation[i] == 'a'&&operation[i + 1] == 't'&&operation[i + 2] == 'a'
			&&operation[i + 3] == 'n'&&operation[i + 4] == 'd') {
			int end = 0;
			for (end = i + 5; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //atanh
		if (operation[i] == 'a'&&operation[i + 1] == 't'&&operation[i + 2] == 'a'
			&&operation[i + 3] == 'n'&&operation[i + 4] == 'h') {
			int end = 0;
			for (end = i + 5; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //atan2
		if (operation[i] == 'a'&&operation[i + 1] == 't'&&operation[i + 2] == 'a'
			&&operation[i + 3] == 'n'&&operation[i + 4] == '2') {
			int end = 0;
			for (end = i + 5; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //atand
		if (operation[i] == 'a'&&operation[i + 1] == 't'&&operation[i + 2] == 'a'
			&&operation[i + 3] == 'n'&&operation[i + 4] == '2'&&operation[i + 4] == 'd') {
			int end = 0;
			for (end = i + 6; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //atan
		if (operation[i] == 'a'&&operation[i + 1] == 't'&&operation[i + 2] == 'a'
			&&operation[i + 3] == 'n') {
			int end = 0;
			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}

	for (int i = operation.length() - 1; i >= 0; i--) { //tanh
		if (operation[i] == 't'&&operation[i + 1] == 'a'&&operation[i + 2] == 'n'
			&&operation[i + 3] == 'h') {
			int end = 0;
			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //tand
		if (operation[i] == 't'&&operation[i + 1] == 'a'&&operation[i + 2] == 'n'
			&&operation[i + 3] == 'd') {
			int end = 0;
			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //tan
		if (operation[i] == 't'&&operation[i + 1] == 'a'&&operation[i + 2] == 'n') {
			int end = 0;
			for (end = i + 3; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}


	for (int i = operation.length() - 1; i >= 0; i--) { //acscd
		if (operation[i] == 'a'&&operation[i + 1] == 'c'&&operation[i + 2] == 's'
			&&operation[i + 3] == 'c'&&operation[i + 4] == 'd') {
			int end = 0;
			for (end = i + 5; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //acsch
		if (operation[i] == 'a'&&operation[i + 1] == 'c'&&operation[i + 2] == 's'
			&&operation[i + 3] == 'c'&&operation[i + 4] == 'h') {
			int end = 0;
			for (end = i + 5; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //acsc
		if (operation[i] == 'a'&&operation[i + 1] == 'c'&&operation[i + 2] == 's'
			&&operation[i + 3] == 'c') {
			int end = 0;
			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}

	for (int i = operation.length() - 1; i >= 0; i--) { //csch
		if (operation[i] == 'c'&&operation[i + 1] == 's'&&operation[i + 2] == 'c'
			&&operation[i + 3] == 'h') {
			int end = 0;
			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //cscd
		if (operation[i] == 'c'&&operation[i + 1] == 's'&&operation[i + 2] == 'c'
			&&operation[i + 3] == 'd') {
			int end = 0;
			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //csc
		if (operation[i] == 'c'&&operation[i + 1] == 's'&&operation[i + 2] == 'c') {
			int end = 0;
			for (end = i + 3; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}



	for (int i = operation.length() - 1; i >= 0; i--) {  //asecd
		if (operation[i] == 'a'&&operation[i + 1] == 's'&&operation[i + 2] == 'e'
			&&operation[i + 3] == 'c'&&operation[i + 4] == 'd') {
			int end = 0;
			for (end = i + 5; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) {  //asech
		if (operation[i] == 'a'&&operation[i + 1] == 's'&&operation[i + 2] == 'e'
			&&operation[i + 3] == 'c'&&operation[i + 4] == 'h') {
			int end = 0;
			for (end = i + 5; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}

	for (int i = operation.length() - 1; i >= 0; i--) {  //asec
		if (operation[i] == 'a'&&operation[i + 1] == 's'&&operation[i + 2] == 'e'
			&&operation[i + 3] == 'c') {
			int end = 0;
			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //sech
		if (operation[i] == 's'&&operation[i + 1] == 'e'&&operation[i + 2] == 'c'
			&&operation[i + 3] == 'h') {
			int end = 0;
			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //secd
		if (operation[i] == 's'&&operation[i + 1] == 'e'&&operation[i + 2] == 'c'
			&&operation[i + 3] == 'd') {
			int end = 0;
			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //sec
		if (operation[i] == 's'&&operation[i + 1] == 'e'&&operation[i + 2] == 'c') {
			int end = 0;
			for (end = i + 3; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}

	for (int i = operation.length() - 1; i >= 0; i--) {  //acotd
		if (operation[i] == 'a'&&operation[i + 1] == 'c'&&operation[i + 2] == 'o'
			&&operation[i + 3] == 't'&&operation[i + 4] == 'd') {
			int end = 0;
			for (end = i + 5; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) {// acoth
		if (operation[i] == 'a'&&operation[i + 1] == 'c'&&operation[i + 2] == 'o'
			&&operation[i + 3] == 't'&&operation[i + 4] == 'h') {
			int end = 0;
			for (end = i + 5; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //acot
		if (operation[i] == 'a'&&operation[i + 1] == 'c'&&operation[i + 2] == 'o'
			&&operation[i + 3] == 't') {
			int end = 0;
			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //cot
		if (operation[i] == 'c'&&operation[i + 1] == 'o'&&operation[i + 2] == 't') {
			int end = 0;
			for (end = i + 3; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //hypot
		if (operation[i] == 'h'&&operation[i + 1] == 'y'&&operation[i + 2] == 'p'
			&&operation[i + 3] == 'o'&&operation[i + 4] == 't') {
			int end = 0;
			for (end = i + 5; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //deg2rad
		if (operation[i] == 'd'&&operation[i + 1] == 'e'&&operation[i + 2] == 'g'
			&&operation[i + 3] == '2'&&operation[i + 4] == 'r'&&operation[i + 5] == 'a'&&operation[i + 6] == 'd') {
			int end = 0;
			for (end = i + 7; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) { //rad2deg
		if (operation[i] == 'r'&&operation[i + 1] == 'a'&&operation[i + 2] == 'd'
			&&operation[i + 3] == '2'&&operation[i + 4] == 'd'&&operation[i + 5] == 'e'&&operation[i + 6] == 'g') {
			int end = 0;
			for (end = i + 7; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) {
		if (operation[i] == 'r'&&operation[i + 1] == 'a'&&operation[i + 2] == 'n'&&operation[i + 3] == 'd') {
			int start = 0;

			for (start = i - 1; start>0; start--)
				if (Is_operation(operation[start - 1]) && start != 1)
					break;
			int end = 0;

			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);

			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) {
		if (operation[i] == 'e'&&operation[i + 1] == 'y'&&operation[i + 2] == 'e') {
			int start = 0;

			for (start = i - 1; start>0; start--)
				if (Is_operation(operation[start - 1]) && start != 1)
					break;
			int end = 0;

			for (end = i + 3; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) {
		if (operation[i] == 'z'&&operation[i + 1] == 'e'&&operation[i + 2] == 'r'&&operation[i + 3] == 'o'&&operation[i + 4] == 's') {
			int start = 0;

			for (start = i - 1; start>0; start--)
				if (Is_operation(operation[start - 1]) && start != 1)
					break;
			int end = 0;

			for (end = i + 5; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) {
		if (operation[i] == 'o'&&operation[i + 1] == 'n'&&operation[i + 2] == 'e'&&operation[i + 3] == 's') {
			int start = 0;

			for (start = i - 1; start>0; start--)
				if (Is_operation(operation[start - 1]) && start != 1)
					break;
			int end = 0;

			for (end = i + 4; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);

			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}

	for (int i = operation.length() - 1; i >= 0; i--) {
		if (operation[i] == 't'&&operation[i + 1] == 'o'&&operation[i + 2] == 'g'&&operation[i + 3] == '1'&&operation[i + 4] == '0') {
			int start = 0;

			for (start = i - 1; start>0; start--)
				if (Is_operation(operation[start - 1]) && start != 1)
					break;
			int end = 0;

			for (end = i + 5; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) {
		if (operation[i] == 'l'&&operation[i + 1] == 'o'&&operation[i + 2] == 'g') {
			int start = 0;

			for (start = i - 1; start>0; start--)
				if (Is_operation(operation[start - 1]) && start != 1)
					break;
			int end = 0;

			for (end = i + 3; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) {
		if (operation[i] == 'e'&&operation[i + 1] == 'x'&&operation[i + 2] == 'p') {
			int start = 0;

			for (start = i - 1; start>0; start--)
				if (Is_operation(operation[start - 1]) && start != 1)
					break;
			int end = 0;

			for (end = i + 3; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			calcAndRep(i, end, operation, operation[i]);
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = 0; i<operation.length() - 1; i++) {
		if (operation[i] == '~'&&operation[i + 1] == 'W'&&operation[i + 2] == '~') {
			int start = 0;

			for (start = i - 1; start>0; start--)
				if (Is_operation(operation[start - 1]) && start != 1)
					break;
			int end = 0;

			for (end = i + 3; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			error = calcAndRep(start, end, operation, operation[i]);
			if (error == 0) {
				operation = "Matrix dimensions must agree";
				return;
			}
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = 0; i<operation.length() - 1; i++) {
		if (operation[i] == '~'&&operation[i + 1] == 'M'&&operation[i + 2] == '~') {
			int start = 0;

			for (start = i - 1; start>0; start--)
				if (Is_operation(operation[start - 1]) && start != 1)
					break;
			int end = 0;

			for (end = i + 3; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			error = calcAndRep(start, end, operation, operation[i]);
			if (error == 0) {
				operation = "Matrix dimensions must agree";
				return;
			}

			//  cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = 0; i<operation.length() - 1; i++) {
		if (operation[i] == '~'&&operation[i + 1] == 'D'&&operation[i + 2] == '~') {
			int start = 0;

			for (start = i - 1; start>0; start--)
				if (Is_operation(operation[start - 1]) && start != 1)
					break;
			int end = 0;

			for (end = i + 3; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			error = calcAndRep(start, end, operation, operation[i]);
			if (error == 0) {
				operation = "Matrix dimensions must agree";
				return;
			}

			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}


	for (int i = operation.length() - 1; i >= 0; i--) {
		if (operation[i] == '^') {
			int start = 0;

			for (start = i - 1; start>0; start--)
				if (Is_operation(operation[start - 1]) && start != 1)
					break;
			int end = 0;

			for (end = i + 1; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			error = calcAndRep(start, end, operation, operation[i]);
			if (error == 0) {
				operation = "Matrix dimensions must agree";
				return;
			}

			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}

	for (int i = 0; i<operation.length() - 1; i++) {
		if (operation[i] == '*') {
			int start = 0;

			for (start = i - 1; start>0; start--)
				if (Is_operation(operation[start - 1]) && start != 1)
					break;
			int end = 0;

			for (end = i + 1; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			error = calcAndRep(start, end, operation, operation[i]);
			if (error == 0) {
				operation = "Matrix dimensions must agree";
				return;
			}
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}


	for (int i = 0; i<operation.length() - 1; i++) {
		if (operation[i] == '/') {
			int start = 0;

			for (start = i - 1; start>0; start--)
				if (Is_operation(operation[start - 1]) && start != 1)
					break;
			int end = 0;

			for (end = i + 1; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			error = calcAndRep(start, end, operation, operation[i]);
			if (error == 0) {
				operation = "Matrix dimensions must agree";
				return;
			}
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) {
		if (operation[i] == '~'&&operation[i + 1] == 'P'&&operation[i + 2] == '~') {
			int start = 0;

			for (start = i - 1; start>0; start--)
				if (Is_operation(operation[start - 1]) && start != 1)
					break;
			int end = 0;

			for (end = i + 3; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			error = calcAndRep(start, end, operation, operation[i]);
			if (error == 0) {
				operation = "Matrix dimensions must agree";
				return;
			}

			// cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}
	for (int i = operation.length() - 1; i >= 0; i--) {
		if (operation[i] == '~'&&operation[i + 1] == 'S'&&operation[i + 2] == '~') {
			int start = 0;

			for (start = i - 1; start>0; start--)
				if (Is_operation(operation[start - 1]) && start != 1)
					break;
			int end = 0;

			for (end = i + 3; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			error = calcAndRep(start, end, operation, operation[i]);
			if (error == 0) {
				operation = "Matrix dimensions must agree";
				return;
			}


			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}


	for (int i = 0; i<operation.length(); i++) {
		if (operation[i] == '-' && i != 0) {
			int end = 0;

			for (end = i + 1; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			if (Is_operation(operation[i - 1]))
				i--;
			int start = 0;

			for (start = i - 1; start>0; start--)
				if (Is_operation(operation[start - 1]) && start != 1)
					break;

			error = calcAndRep(start, end, operation, operation[i]);
			if (error == 0) {
				operation = "Matrix dimensions must agree";
				return;
			}
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}

	for (int i = operation.length() - 1; i >= 0; i--) {
		if (operation[i] == '+') {
			int start = 0;

			for (start = i - 1; start>0; start--)
				if (Is_operation(operation[start - 1]) && start != 1)
					break;

			int end = 0;

			for (end = i + 1; end<operation.length() - 1; end++)
				if (Is_operation(operation[end + 1]))
					break;
			error = calcAndRep(start, end, operation, operation[i]);
			if (error == 0) {
				operation = "Matrix dimensions must agree";
				return;
			}
			//cout<<operation<<endl;
			Operation_solver(operation);
			break;
		}
	}




}


string mul_ope_solver(string &ope)
{
	Operation_solver(ope);

	if (memoryCheck(ope) == -1 &&
		ope != "Statement or Expression is incorrect"&& ope != "Matrix dimensions must agree")
	{
		string temp2, tName;
		tName = genRandom();
		temp2 = tName + "=[" + ope + "]";
		input_checker(temp2);
		//cout<<"operatin is"<<ope<<endl;
		return tName;
	}
	return ope;

}

/*End Here*/

int main(int argv, char* argc[])

{

	ios_base::sync_with_stdio(false);

	cin.tie(0);



	if (argv>1)

	{

		ifstream infile(argc[1]);

		string sFile="", temp;
		int Ncomplete=0;
		

		while (getline(infile, temp))

		{

			if (temp.find("\r") != -1)
			temp.replace(temp.find("\r"), 2, "");


			
			sFile+=temp;
			if(sFile.find('[')!=-1||sFile.find(']')!=-1){
				for(int i=0; i<sFile.length() ;i++)
				{
					if(sFile[i]=='[')
						Ncomplete++;
					else if(sFile[i]==']')
						Ncomplete--;

				}
			}
			
				if(Ncomplete==0){

			input_checker(sFile);

				sFile = "";
				Ncomplete=0;
				}



			



			/*if (sFile.find("]") != -1 || sFile.find("];") != -1 || sFile.find("+") != -1 || sFile.find("*") != -1 || sFile.find("/") != -1 || sFile.find("'") != -1 || sFile.find("./") != -1 || (sFile.find("-") != -1 && sFile.length() <= 10))

			{*/

				

			//}

		}

		infile.close();

	}

	else

		while (1) {

			exit1 = 0;

			cout << ">> ";



			string ins;

			getline(cin, ins);

			if (ins.find("\r") != -1)

				ins.replace(ins.find("\r"), 2, "");

			input_checker(ins);

			if (exit1 == 1)

				break;





			/*

			cout<<endl<<"memory----------------"<<endl;

			for(int z=0;z<memoryPointer;z++)

			{

			cout<<memory[z].getName()<<" "<<memory[z].getRows()<<" "<<memory[z].getColumns()<<endl;

			memory[z].print();

			cout<<"#############"<<endl;

			}

			cout<<endl<<"--------------------------"<<endl;

			*/

		}



	return 0;

}

