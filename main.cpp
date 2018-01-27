#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <math.h>
#define PI 3.14159265359
using namespace std;
void removeSpaces2(string &str)
{
	while(str[str.length()-1] ==' ')
		str.erase(str.length()-1,1);
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
//=============================================================================
sizeValue sizing(string matrix)  //take any string but without brackets and retur a structure which has the number of rows and columns
{                               // for example    22 45 sin(90)+1;1 cos(0) 5    size will be 2x3
	sizeValue n; n.rows = 0;  n.columns = 0;
	stringstream ss(matrix);
	string token;
	while (getline(ss, token, ';'))
		n.rows++;
	getline(ss, token, ';');
	stringstream sn(token);
	string in;
	while (sn >> in)
		n.columns++;
	return n;
}

sizeValue conc(string s)
{
    int i=0, j=0, k=0;
    string miniMatrix;
    sizeValue Vstack[2]; //sizeValue Vstack[2]; supposed
    Vstack[1].rows = 0;
    Vstack[1].columns = 0;
    while(1)
    {
    i = s.rfind('['); //if not found returns -1
	if(i == -1)
	{   
		int flag=0;
		for(int y=0 ; y<s.length() ; y++)
		{
			if(s[y] != ',' && s[y] != '[' && s[y] != ']')
			{
				flag = 1;
			}
		}
		if(flag == 1)
		{
			Vstack[k] = sizing(s);
			Vstack[0] = compare(Vstack[0], Vstack[1]);
		}
		break;
	}
    for(int o=i ; o<s.length() ; o++)
    {
        if(s[o] == ']')
        {
            j = o;
            break;
        }
    }
    if(1)//i != -1 )//&& j != -1)
    {
        miniMatrix = s.substr(i, j-i+1);
        Vstack[k] = sizing(miniMatrix);
		if(Vstack[1].rows != 0)
        {
            Vstack[0] = compare(Vstack[0], Vstack[1]);
            Vstack[1].rows = 0; //it can be overwritten
            Vstack[1].columns = 0;
        }
        s.erase(i, j-i+1);
	    if(s.rfind('[') == -1)
		{
			int flag=0;
		for(int y=0 ; y<s.length() ; y++)
		{
			if(s[y] != ',' && s[y] != '[' && s[y] != ']')
			{
				flag = 1;
			}
		}
		if(flag == 1)
		{
			Vstack[1] = sizing(s);
			Vstack[0] = compare(Vstack[0], Vstack[1]);
		}
		break;
		}
		removeSpaces2(s);
		for(int o=i ; o<=(s.length()) ; o++)
		{
			if((s[o] == ']' && s[o-1] == '[')
			|| (s[o] == ']' && s[o-1] == ',' && s[o-2] == '[')
			|| (s[o] == ']' && s[o-1] == ',')
			|| (s[o] == '[' && s[o-1] == ',')
			|| (s[o] == ',' && s[o-1] == ']')
			|| (s[o] == ',' && s[o-1] == '[')
			)
			{
				if((s[o] == ']' && s[o-1] == '['))
					s.erase(o-1, 2);
				else if((s[o] == ']' && s[o-1] == ',' && s[o-2] == '['))
					s.erase(o-2, 3);
				else if((s[o] == ']' && s[o-1] == ','))
					s.erase(o-1, 2);
				else if((s[o] == '[' && s[o-1] == ','))
					s.erase(o-1, 2);
				else if((s[o] == ',' && s[o-1] == ']'))
					s.erase(o-1, 2);
				else if((s[o] == ',' && s[o-1] == '['))
					s.erase(o-1, 2);
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
    int k=0;
    for(int i=0 ; i<separatedString.size() ; i++)
    {
        if(separatedString[i].find('[') == -1)
            finStack.push_back(sizing(separatedString[i]));
        else
            finStack.push_back(conc(separatedString[i]));
        k++;
    }
    if(finStack.size() == 1) 
    {
		separatedString.clear(); 
	    return finStack[0];
	}
    sizeValue sum = compare(finStack[0], finStack[1]);
    for(int i=2 ; i<finStack.size() ; i++)
    {
        sum = compare(sum, finStack[i]);
    }
	separatedString.clear();
    return sum;
}

double StringToDouble(const string &text)
{
	stringstream ss(text);
	double result;
	return ss >> result ? result : 0;
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
		//----------------------------------------------------
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
		for (int i = 0; i<rows; i++)
		{
			delete[] element[i];
		}
		delete[] element;
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
void sFill(matrix &mSoph, string mString)
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
				if (index != -1)
				{
					//B
					//memory.p[index]
					for (int rs = 0; rs < memory.p[index].rows; rs++)
					{
						for (int cs = 0; cs < memory.p[index].columns; cs++)
						{
							mSoph.element[r + rs][c + cs].value = memory.p[index].element[rs][cs].value;
							mSoph.element[r + rs][c + cs].isFilled = 1;
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
						for (int rs = 0; rs < tempM.rows; rs++)
						{
							for (int cs = 0; cs < tempM.columns; cs++)
							{
								mSoph.element[r + rs][c + cs].value = tempM.element[rs][cs].value;
								mSoph.element[r + rs][c + cs].isFilled == 1;
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
						for (int rs = 0; rs < tempM.rows; rs++)
						{
							for (int cs = 0; cs < tempM.columns; cs++)
							{
								mSoph.element[r + rs][c + cs].value = tempM.element[rs][cs].value;
								mSoph.element[r + rs][c + cs].isFilled == 1;
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
						for (int rs = 0; rs < tempM.rows; rs++)
						{
							for (int cs = 0; cs < tempM.columns; cs++)
							{
								mSoph.element[r + rs][c + cs].value = tempM.element[rs][cs].value;
								mSoph.element[r + rs][c + cs].isFilled == 1;
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
						for (int rs = 0; rs < tempM.rows; rs++)
						{
							for (int cs = 0; cs < tempM.columns; cs++)
							{
								mSoph.element[r + rs][c + cs].value = tempM.element[rs][cs].value;
								mSoph.element[r + rs][c + cs].isFilled == 1;
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
						for (int rs = 0; rs < tempM.rows; rs++)
						{
							for (int cs = 0; cs < tempM.columns; cs++)
							{
								mSoph.element[r + rs][c + cs].value = tempM.element[rs][cs].value;
								mSoph.element[r + rs][c + cs].isFilled == 1;
							}
						}
						lastPos += temp.length();
						flag = 1;
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
				OBCounter--;
				if (OBCounter == 0)
				{
					OBColumnCounter[OBCounter] = 0;
					OBRowCounter[OBCounter] = 0;
				}
			}
			//if 2.3
			else if ((int(temp[0]) <= 57) && (int(temp[0]) >= 48))
			{

				int spacePos = mString.find(' ', lastPos)
					, semicolumnPos = mString.find(';', lastPos)
					, CBPos = mString.find(']', lastPos);
				if (spacePos == -1) spacePos = 999999;
				if (semicolumnPos == -1) semicolumnPos = 999999;
				//Because the separators inside [] are numerous
				if ((spacePos < semicolumnPos) && (spacePos < CBPos))
				{
					temp = mString.substr(lastPos, spacePos - lastPos);
				}
				else if ((semicolumnPos < spacePos) && (semicolumnPos < CBPos))
				{
					temp = mString.substr(lastPos, semicolumnPos - lastPos);
				}
				else
				{
					temp = mString.substr(lastPos, CBPos - lastPos);
				}

				stringstream ss;
				ss << temp;
				double value;
				ss >> value;
				mSoph.element[r + OBRowCounter[OBCounter]][c + OBColumnCounter[OBCounter]].value = value;
				mSoph.element[r + OBRowCounter[OBCounter]][c + OBColumnCounter[OBCounter]].isFilled = 1;

				if ((spacePos < semicolumnPos) && (spacePos < CBPos))
				{
					if (OBCounter != 0)
					{
						OBColumnCounter[OBCounter]++;//walking through the inner matrix
						c--;//freezing c as it will increase by one the next loop
					}
					lastPos = spacePos + 1;
				}
				else if ((semicolumnPos < spacePos) && (semicolumnPos < CBPos))
				{
					if (OBCounter != 0)
					{
						OBRowCounter[OBCounter]++;//walking through the inner matrix
						OBColumnCounter[OBCounter] = 0;//reseting column to start position
						c--;
					}

					lastPos = semicolumnPos + 1;
				}
				else// ]
				{
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
		string mString = input.substr(FOBPos , (LCBPos - FOBPos + 1));// Entering with []
																			 //call concatenation function to calculate the rows and the columns of the matrix and to create it
		separate(mString);
		sizeValue mSize = calcSize(separatedString);
		memorizeMatrix(index, mSize.rows, mSize.columns, mName);
		index = memoryCheck(mName);
		mString = input.substr(FOBPos + 1, (LCBPos - 2 - FOBPos + 1));
		sFill(memory.p[index], mString);
		// memory.p[index]; this will be the matrix im working on
		//memory.p[index].sFill(mString); sophisticated filling
	}
	else if (assignmentOP)
	{
		string mString = input.substr(FOBPos + 1, (LCBPos - 2 - FOBPos + 1));//FOB+1 & LCB-2 to remove braces
		memorizeMatrix(index, mString, mName); //I think we may handle sin() & 1X1 matrix inside this func
	}

	else if (mathOP)
	{
		string operation = input.substr(EQPos + 1, (input.length() - EQPos + 1));
		removeSpaces(operation);
		string variable1, variable2;
		int index1, index2;


		//memory.push_back(temp);
		if (plusPos != -1)
		{
			cut(variable1, variable2, index1, index2, '+', operation);


			if (index != -1)
			{
				asg = (index1 == index || index2 == index);
				memory.p[index].add(memory.p[index1], memory.p[index2], asg);
				memory.p[index].print();
			}

			else
			{
				memory.create(mName);
				memory.p[memoryPointer].add(memory.p[index1], (memory.p[index2]));
				memory.p[memoryPointer].print();
				memoryPointer++;
			}
		}

		else if (minusPos != -1)
		{
			cut(variable1, variable2, index1, index2, '-', operation);


			if (index != -1)
			{
				asg = (index1 == index || index2 == index);
				memory.p[index].sub(memory.p[index1], (memory.p[index2]), asg);
				memory.p[index].print();
			}

			else
			{
				memory.create(mName);
				memory.p[memoryPointer].sub(memory.p[index1], (memory.p[index2]));
				memory.p[memoryPointer].print();
				memoryPointer++;
			}


		}
		else if (multPos != -1)
		{
			cut(variable1, variable2, index1, index2, '*', operation);

			if (index != -1)
			{
				asg = (index1 == index);
				if (index2 == index) asg = 2;
				memory.p[index].mult(memory.p[index1], (memory.p[index2]), asg);
				memory.p[index].print();
			}

			else
			{
				memory.create(mName);
				memory.p[memoryPointer].mult(memory.p[index1], (memory.p[index2]));
				memory.p[memoryPointer].print();
				memoryPointer++;
			}
		}
		else if (elemWiseInvPos != -1)
		{
			cut(variable1, index1, './', operation);
			if (index != -1)
			{
				memory.p[index].inversePerElement(memory.p[index1]);
				memory.p[index].print();
			}

			else
			{
				memory.create(mName);
				memory.p[memoryPointer].inversePerElement(memory.p[index1]);
				memory.p[memoryPointer].print();
				memoryPointer++;
			}
		}
		else if (divPos != -1)
		{
			cut(variable1, variable2, index1, index2, '/', operation);
			//cout<<memory.p[index2].getDeterminant()<<endl;
			if (memory.p[index2].getDeterminant() == 0)
			{
				cout << "Division cannot done" << endl;
				return;
			}

			else if (index != -1)
			{
				memory.p[index].div(memory.p[index1], (memory.p[index2]));
				memory.p[index].print();
			}

			else
			{
				memory.create(mName);
				memory.p[memoryPointer].div(memory.p[index1], (memory.p[index2]));
				memory.p[memoryPointer].print();
				memoryPointer++;
			}
		}
		else if (transPos != -1)
		{
			string var;
			int ind;
			var = operation.substr(0, operation.find("'"));
			ind = memoryCheck(var);


			if (index != -1)
			{
				memory.p[index].getTranspose(memory.p[ind]);
				memory.p[index].print();
			}

			else
			{
				memory.create(mName);
				memory.p[memoryPointer].getTranspose(memory.p[ind]);
				memory.p[memoryPointer].print();
				memoryPointer++;
			}
		}
	}


	else
	{
		removeSpaces(input);
		if ((memoryCheck(input) != -1) &&
			memoryCheck(input)<memoryPointer)
		{
			memory.p[memoryCheck(input)].print();
		}
		else if (input == "exit")
		{
			exit1 = 1;
			return;
		}
		else
			cout << "invalid input" << endl;
	}


}
#define endl '\n'

int main(int argv, char* argc[])
{
	ios_base::sync_with_stdio(false);
	cin.tie(0);

	if (argv>1)
	{
		ifstream infile(argc[1]);
		string sFile, temp;
		while (getline(infile, temp))
		{
			if (temp.find("\r") != -1)
				temp.replace(temp.find("\r"), 2, "");
			sFile += temp;
			if (sFile.find("]") != -1 || sFile.find("];") != -1 || sFile.find("+") != -1 || sFile.find("*") != -1 || sFile.find("/") != -1 || sFile.find("'") != -1 || sFile.find("./") != -1 || (sFile.find("-") != -1 && sFile.length() <= 10))
			{
				input_checker(sFile);
				sFile = "";
			}
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
