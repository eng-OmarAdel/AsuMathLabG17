#include <iostream>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include<string>
#include<math.h>
#include <fstream>
#include<vector>
#include<math.h>
#define PI 3.14159265359
using namespace std;

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
double atan2(double x,double y)
{
	return atan(x/y);
}
double atan2d(double x,double y)
{
	return ((180.0 / PI)*atan(x/y));
}
// csc functions --> csc , cscd ,acsc ,acscd ,csch ,acsch
double csc(double x)
{
	return ( 1.0 / sin(x) );
}
double cscd(double x)
{
	return ( 1.0 / sind(x) );
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
	return sqrt(abs(x)*abs(x)+abs(y)*abs(y));
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

class matrix
{
    public:
    int rows;
    int columns;
    string mName,errorHandler;

    double** element;
    int getRows()
    {
        return rows;
    }

    int getColumns()
    {
        return columns;
    }
     void setRows (int r)
    {
        rows=r;
    }
    void setColumns (int c)
    {
        columns=c;
    }

    /*double** getElementArray() //if we make element private
    {
        return element;
    }*/


    matrix()
    {
      rows=0;
      columns=0;
	  mName="NULL" ;
    }

    //copy constructor
    matrix(const matrix &mCopied)
    {
        rows=mCopied.rows;
        columns=mCopied.columns;
        mName=mCopied.mName;

        if(rows==0&&columns==0)
            element=NULL;
        else
        {
        element=new double*[rows];
        for(int i=0;i<rows;i++)
        {
            element[i]=new double[columns];
        }
            for(int i=0;i<rows;i++)
            {
                for(int j=0;j<columns;j++)
                {
                    element[i][j]=mCopied.element[i][j];
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
    double* getElement(int elementRow, int elementColumn)
    {
        return this->element[elementRow, elementColumn];
    }

    void initialling(string mName,string mString) // give me string wana azzabat isa
    {
        this->mName=mName;
        stringstream ss(mString);
        string token;
        while(getline(ss, token, ';'))
        {
            this->rows++;
            stringstream sn;
            sn<<token;
            double in;
            while( sn >> in)
            {
                if(rows==1)
                    this->columns++;
            }

        }
    //----------------------------------------------------
        element = new double*[rows];
        for(int i = 0; i < rows; ++i)
        element[i] = new double[columns];
        int p=0, q=0;
        stringstream ss3(mString);
        string token1;
        while(getline(ss3, token1, ';'))
        {
            stringstream ss1(token1);
            string token2;
            while(ss1>>token2)
            {
                element[p][q] = StringToDouble(token2); //p rows --- q columns
                q++;
            }
            q=0;
            p++;
        }
    }
    void update(string mName,string mString)
    {

        setRows(0);
        setColumns(0);
        for(int i=0;i<rows;i++)
		{
			delete[] element[i];
		}
		delete[] element;
		initialling(mName,mString);

    }

    void initialling(int rows1, int columns1)//give rows and columns wana azzabat isa
    {
        this->columns = columns1;
        this->rows = rows1;
        element=new double* [rows1];
        for(int i=0; i<rows1;i++)
            element[i]=new double[columns1];
    }

    void setElement(int elementRow, int elementColumn, double elementValue) // modify the main object's 2d array
    {
        element[elementRow][elementColumn] = elementValue;
    }

     void getTranspose(matrix &x)
    {

        this->initialling(x.columns, x.rows);
        for(int m=0 ; m<(this->getRows()) ; m++)
        {
            for(int n=0 ; n<(this->getColumns()) ; n++)
            {
                this->setElement(m,n,x.element[n][m]);
            }
        }

    }

    void subMatrix(matrix& x,int elementRow, int elementColumn)
    {
        this->initialling((x.rows)-1, (x.columns)-1);
        for(int i=0,m=0 ; i<x.rows,m<(this->getRows()) ; i++)
        {
            if(i == elementRow) continue;
            else
            {
                for(int j=0,n=0 ; j<x.columns,n<(this->getColumns()) ; j++)
                {
                    if(j == elementColumn) continue;
                    else
                    {
                        this->setElement(m,n,x.element[i][j]);
                        n++;
                    }
                }
                m++;
            }
        }

    }

    //--------------------------------------------------------------------
	double getDeterminant() {
		if(rows==columns){
		int i, j, k;
		double factor;
		double temp;
		matrix a(*this) ;
		int counti=0;
			int m=this->rows ;
	for(i = 0; i < m - 1; i++)
	{
		/* Elementary Row Operation I */
		if(a.element[i][i] == 0)
		{
			for(k = i; k < m; k++)
			{
				if(a.element[k][i] != 0)
				{
					for(j = 0; j < m; j++)
					{
						temp = a.element[i][j];
						a.element[i][j] = a.element[k][j];
						a.element[k][j] = temp;
					}
				k = m;
				}
			}
			counti++;
		}
		/* Elementary Row Operation III */
		if(a.element[i][i] != 0)
		{
			for(k = i + 1; k < m; k++)
			{
				factor = -1.0 * a.element[k][i] /  a.element[i][i];
				for(j = i; j < m; j++)
				{
					a.element[k][j] = a.element[k][j] + (factor * a.element[i][j]);
				}
			}
		}
	}

	/* Display upper triangular matrix */

	temp = 1.0;


	for(i = 0; i < m; i++)
	{
		temp *= a.element[i][i];
	}


	if(counti % 2 == 0)
	{
		return temp;;
	}
	else
	{
		return -1*temp ;
	}
		}
		else return 0;

}
		//--------------------------------------------------------------------------
    void getInverse(matrix &x)
    {
	matrix z;
        double detObj = x.getDeterminant();
        this->initialling(x.rows, x.columns);
	z.initialling(x.rows, x.columns);
	matrix sub;
	int nozero=0;
        for(int i=0 ; i<x.rows ; i++)
        {
            for(int j=0 ; j<x.columns ; j++)
            {
		sub.subMatrix(x,i,j);
                double minor = sub.getDeterminant();
				if(minor==0){
					nozero=nozero/ detObj;
					z.setElement(i,j,nozero);
				}
				else{
			if((i+j)%2!=0) minor *= -1;




                minor = minor / detObj;
                z.setElement(i,j,minor);
				}
            }

        }
        this->getTranspose(z);
    }
	void inversePerElement(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				if(x.element[m][n]!=0)
					this->setElement(m, n, 1/(x.element[m][n]));
				else
					this->errorHandler="Error There's a zero element in the matrix"; //this is to handle 1/0 error me7taga tet3addel tab3an
			}
		}
	}

    void add(matrix& x, matrix& y,int old=0)
    {

         if(old==0)
        this->initialling(x.rows, x.columns);
        for(int i=0 ; i< (this->rows) ; i++)
        {
            for(int j=0 ; j< (this->columns) ; j++)
            {
                element[i][j] = x.element[i][j] + y.element[i][j];
            }
        }

    }

     void sub(matrix& x, matrix& y,int old=0)
    {
        if(old==0)
        this->initialling(x.rows, x.columns);
        for(int i=0 ; i< (this->rows) ; i++)
        {
            for(int j=0 ; j< (this->columns) ; j++)
            {
                element[i][j] = x.element[i][j] - y.element[i][j];
            }
        }

    }

    void mult(matrix &x,matrix &y,int asg=0)//remember to handle errors of dimension
	{
		//declaration the output(returned) matrix
		if(asg==0)
        this->initialling(x.rows, y.columns);//use constructor instead
        double** temp;
        if(asg==1)
        {
            temp=new double*[x.rows];
            for(int w=0;w<x.rows;w++)
                temp[w]=new double [x.columns];
            for(int w=0;w<x.rows;w++)
                for(int q=0;q<x.columns;q++)
                temp[w][q]=x.element[w][q];
        }
        if(asg==2)
        {
            temp=new double*[y.rows];
            for(int w=0;w<y.rows;w++)
                temp[w]=new double [y.columns];
            for(int w=0;w<y.rows;w++)
                for(int q=0;q<y.columns;q++)
                temp[w][q]=y.element[w][q];
        }

        //filling the matrix with valus
		for(int i=0;i<rows;i++)
			for(int j=0; j<columns ; j++)
				element[i][j]=0;
		for (int i = 0; i < x.rows; i++)
		{
			for (int j = 0; j < y.columns; j++)
			{
				for (int k = 0; k < y.rows; k++)
				{
				    if(asg==0)
					    element[i][j] += x.element[i][k] * y.element[k][j]; //use setElement function instead
					else if(asg==1)
                        element[i][j] += temp[i][k] * y.element[k][j];
                    else if(asg==2)
                        element[i][j] += x.element[i][k] * temp[k][j];
				}
			}
		}

	}

	void div(matrix &x, matrix &y)
	{
	    matrix inverseDenom;
        inverseDenom.getInverse(y);
        this->initialling(x.rows, inverseDenom.getColumns());
        this->mult(x,inverseDenom);


	}
	//trigonometric functions , any trigonometric function for a matrix is started with 'M'
	//sin functions --> sin , sind , asin , asind , sinh , asinh
	void Msin(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, sin(x.element[m][n]));
			}
		}
	}
	void Msind(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, sind(x.element[m][n]));
			}
		}
	}
	void Masin(matrix &x)  // error handling
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, asin(x.element[m][n]));
			}
		}
	}
	void Masind(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (asind(x.element[m][n])));
			}
		}
	}
	void Msinh(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, sinh(x.element[m][n]));
			}
		}
	}
	void Masinh(matrix &x)  // error handling
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, asinh(x.element[m][n]));
			}
		}
	}
	//cos functions --> cos , cosd, acos , acosd ,cosh, acosh
	void Mcos(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, cos(x.element[m][n]));
			}
		}
	}
	void Mcosd(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, cosd(x.element[m][n]));
			}
		}
	}
	void Macos(matrix &x)  // error handling
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, acos(x.element[m][n]));
			}
		}
	}
	void Macosd(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (acosd(x.element[m][n])));
			}
		}
	}
	void Mcosh(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, cosh(x.element[m][n]));
			}
		}
	}
	void Macosh(matrix &x)  // error handling
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, acosh(x.element[m][n]));
			}
		}
	}
	// tan functions --> tan,tand,atan,atand,tanh,atanh
	void Mtan(matrix &x) //error handling
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, tan(x.element[m][n]));
			}
		}
	}
	void Mtand(matrix &x)//error handling
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, tand(x.element[m][n]));
			}
		}
	}
	void Matan(matrix &x)  // error handling
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, atan(x.element[m][n]));
			}
		}
	}
	void Matand(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (atand(x.element[m][n])));
			}
		}
	}
	void Mtanh(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, tanh(x.element[m][n]));
			}
		}
	}
	void Matanh(matrix &x)  // error handling
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, atanh(x.element[m][n]));
			}
		}
	}
	// csc functions --> csc , cscd ,acsc ,acscd ,csch ,acsch
	void Mcsc(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, csc(x.element[m][n]));
			}
		}
	}
	void Mcscd(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, cscd(x.element[m][n]));
			}
		}
	}
	void Macsc(matrix &x)  // error handling
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, acsc(x.element[m][n]));
			}
		}
	}
	void Macscd(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (acscd(x.element[m][n])));
			}
		}
	}
	void Mcsch(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, csch(x.element[m][n]));
			}
		}
	}
	void Macsch(matrix &x)  // error handling
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, acsch(x.element[m][n]));
			}
		}
	}
	//sec functions --> sec ,secd ,asec, asecd ,sech ,asech
	void Msec(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, sec(x.element[m][n]));
			}
		}
	}
	void Msecd(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, secd(x.element[m][n]));
			}
		}
	}
	void Masec(matrix &x)  // error handling
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, asec(x.element[m][n]));
			}
		}
	}
	void Masecd(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (asecd(x.element[m][n])));
			}
		}
	}
	void Msech(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, sech(x.element[m][n]));
			}
		}
	}
	void Masech(matrix &x)  // error handling
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, asech(x.element[m][n]));
			}
		}
	}
	// cot functions cot , cotd , acot ,acotd ,coth ,acoth
	void Mcot(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, cot(x.element[m][n]));
			}
		}
	}
	void Mcotd(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, cotd(x.element[m][n]));
			}
		}
	}
	void Macot(matrix &x)  // error handling
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, acot(x.element[m][n]));
			}
		}
	}
	void Macotd(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (acotd(x.element[m][n])));
			}
		}
	}
	void Mcoth(matrix &x)
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, coth(x.element[m][n]));
			}
		}
	}
	void Macoth(matrix &x)  // error handling
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, acoth(x.element[m][n]));
			}
		}
	}

	//end trignometric functions


	//element wise operators
	void addEL(matrix &x, double y)//EL stands for Element Wise  A+2 or A+1
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (x.element[m][n])+y);
			}
		}
	}
	void subEL(matrix &x, double y)//EL stands for Element Wise  A-2 or A-1
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (x.element[m][n])-y);
			}
		}
	}

	void multEL(matrix &x, double y)//EL stands for Element Wise  A.*2 or A.*1  same as A*2 or A*1
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n,(x.element[m][n])*y);
			}
		}
	}
	void divEL(matrix &x, double y)//EL stands for Element Wise A./2 or A./1  same as A/2 or A/1
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (x.element[m][n])/y);
			}
		}
	}
	void multEL(matrix &x, matrix& y)//EL stands for Element Wise A.*B  or  A.* c etc..
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (x.element[m][n])*(y.element[m][n]));
			}
		}
	}
	void divEL(matrix &x, matrix& y)//EL stands for Element Wise A./B  or  A./c etc..
	{
		this->initialling(x.rows, x.columns);
		for (int m = 0; m<(this->getRows()); m++)
		{
			for (int n = 0; n<(this->getColumns()); n++)
			{
				this->setElement(m, n, (x.element[m][n])/(y.element[m][n]));
			}
		}
	}
	//end element wise operators
	///// farouk
	    void multforpower(matrix &x,matrix &y)           //multiplication for square matrices and assigning the result in the first matrix
	{

        this->initialling(x.rows, x.columns);

		for(int i=0;i<rows;i++)
			for(int j=0; j<columns ; j++)
				element[i][j]=0;
		for (int i = 0; i < x.rows; i++)
		{
			for (int j = 0; j < x.columns; j++)
			{
				for (int k = 0; k < x.rows; k++)
				{

                    element[i][j] += x.element[i][k] * y.element[k][j];
				}
			}
		}
		for (int i = 0; i < x.rows; i++)
            for (int j = 0; j < y.columns; j++)
                x.element[i][j] = element[i][j];

	}
	void identityMatrix()                                      //identity matrix I
	{

		for(int i=0;i<rows;i++)
			for(int j=0; j<columns ; j++)
				element[i][j]= (i == j);
	}

    void power(matrix &x, double power)                  //matrix power
	{
	    if(x.rows == x.columns)                          //must be square matrix
        {
            int mod = power*10;                          //To check if power is fractional
            if(power == 0)                               //produce identity matrix if power = 0
            {
                this->initialling(x.rows,x.columns);
                this->identityMatrix();
            }
            else if(mod % 10 != 0)                       //To check if power is fractional
            {
                if(x.rows == 1)                          //Support fractional power for 1*1 matrix only
                {
                    double y = pow(x.element[0][0], power);
                    this->initialling(1,1);
                    element[0][0] = y;
                }
                else
                    errorHandler = "Error: Fraction power is supported in 1*1 matrix only.";
            }
            else
            {
                int intPower;
                if(power < 0)
                  intPower = (int)abs(power);
                else
                    intPower = (int)power;
                matrix y;                                         //Two matrices to be used in calculation
                y.initialling(x.rows, x.columns);
                matrix temp;
                temp.initialling(x.rows,x.columns);
                for(int i=0;i<temp.rows;i++)
                        for(int j=0; j<temp.columns ; j++)
                            temp.element[i][j] = x.element[i][j];
                y.identityMatrix();
                while (intPower > 0)                             //matrix power by exponentiation by squaring algorithm
                {
                    if (intPower % 2 == 1)
                    {
                    this->multforpower(y,temp);
                    }

                    this->multforpower(temp, temp);
                    intPower /= 2;


                }
                this->initialling(y.rows,y.columns);              //Get inverse if  power is negative
                if(power < 0)
                {
                    this->getInverse(y);

                }
                else
                {
                    for(int i=0;i<rows;i++)
                        for(int j=0; j<columns ; j++)
                            element[i][j] = y.element[i][j];
                }

            }
        }

            else
                errorHandler = "Error: for A^x, A must be a square matrix.";

	}
	void elementWisePower(matrix &x, double power)                  //matrix power
	{
            int mod = power*10;                          //To check if power is fractional
            if(power == 0)                               //produce identity matrix if power = 0
            {
                this->initialling(x.rows,x.columns);
                 for(int i=0;i<rows;i++)
                        for(int j=0; j<columns ; j++)
                            element[i][j] = 1;
            }
            else if(mod % 10 != 0)                       //To check if power is fractional
            {
                this->initialling(x.rows,x.columns);
                 for(int i=0;i<rows;i++)
                        for(int j=0; j<columns ; j++)
                            element[i][j] = pow(x.element[i][j],power);
            }
            else
            {
                int intPower;
                if(power < 0)
                  intPower = (int)abs(power);
                else
                  intPower = (int)power;
                matrix y;                                         //Two matrices to be used in calculation
                y.initialling(x.rows, x.columns);
                for(int i=0;i<y.rows;i++)
                        for(int j=0; j<y.columns ; j++)
                            y.element[i][j] = 1;
                matrix temp;
                temp.initialling(x.rows,x.columns);
                for(int i=0;i<temp.rows;i++)
                        for(int j=0; j<temp.columns ; j++)
                            temp.element[i][j] = x.element[i][j];

                while (intPower > 0)                             //matrix power by exponentiation by squaring algorithm
                {
                    if (intPower % 2 == 1)
                    {
                    for(int i=0;i<temp.rows;i++)
                        for(int j=0; j<temp.columns ; j++)
                            y.element[i][j] = y.element[i][j]*temp.element[i][j];
                    }

                    for(int i=0;i<temp.rows;i++)
                        for(int j=0; j<temp.columns ; j++)
                            temp.element[i][j] = temp.element[i][j]*temp.element[i][j];
                    intPower /= 2;


                }
                this->initialling(y.rows,y.columns);              //Get inverse if  power is negative
                if(power < 0)
                {
                    this->inversePerElement(y);

                }
                else
                {
                    for(int i=0;i<rows;i++)
                        for(int j=0; j<columns ; j++)
                            element[i][j] = y.element[i][j];
                }

            }

	}
	void logMatrix(matrix &x)
	{
	    this->initialling(x.rows,x.columns);
	    for(int i=0;i<rows;i++)
                        for(int j=0; j<columns ; j++)
                            element[i][j] = log(x.element[i][j]);
	}
	void log10Matrix(matrix &x)
	{
	    this->initialling(x.rows,x.columns);
	    for(int i=0;i<rows;i++)
                        for(int j=0; j<columns ; j++)
                            element[i][j] = log10(x.element[i][j]);
	}
	void sqrtMatrix(matrix &x)
	{
	    this->initialling(x.rows,x.columns);
	    for(int i=0;i<rows;i++)
                        for(int j=0; j<columns ; j++)
                            element[i][j] = sqrt(x.element[i][j]);
	}
	void expMatrix(matrix &x)
	{
	    this->initialling(x.rows,x.columns);
	    for(int i=0;i<rows;i++)
                        for(int j=0; j<columns ; j++)
                            element[i][j] = exp(x.element[i][j]);
	}
	//////
	void print()
	{
		if(errorHandler=="Error There's a zero element in the matrix" ||
     errorHandler=="Error The determinant of this matrix is eual to zero")
			cout<<errorHandler;
		else
		{
		   // if(mName[0]!='@' && mName[0]!='&' &&
    //mName[0]!='#' && mName[0]!='_' )
          {
             cout << endl;
			cout << mName << " = " << endl;
			for (int i = 0; i<rows; i++)
			{
				for (int j = 0; j<columns; j++)
				{
					cout << "\t" << element[i][j];
				}
				cout << endl;
			}
          }

		}
		cout<<endl;
	}

	~matrix()
	{
		for(int i=0;i<rows;i++)
		{
			delete[] element[i];
		}
		delete[] element;
		element = NULL;

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
        p=new matrix[1];
        usedSlots=0;
    }
    mDynArr(int s)
    {
        arraySize=s;
        if(arraySize<=0)
            arraySize=1;
        p=new matrix[arraySize];
        usedSlots=0;
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
    void create(string name,string mString)
    {
        p[usedSlots].initialling(name,mString);
        refresh();
    }
    void update(int num,string name,string mString)
    {
        p[num].update(name,mString);
    }
    void refresh()
    {
        usedSlots++;
        if(usedSlots==arraySize)
        {

            arraySize+=10;
            matrix *t=new matrix [arraySize];
            for(int i=0;i<usedSlots;i++)
            {
                t[i].setRows(p[i].getRows());
                t[i].setColumns(p[i].getColumns());
                t[i].setName(p[i].getName());

                t[i].element=new double*[t[i].getRows()];
                for(int o=0;o<p[i].getRows();o++)
                {
                   t[i].element[o]=new double[t[i].getColumns()];

                }
                for(int k=0;k<p[i].getRows();k++)
                {
                    for(int j=0;j<p[i].getColumns();j++)
                    {

                        t[i].element[k][j]=p[i].element[k][j];
                    }
                }
                 for(int g=0;g<p[i].getRows();g++)
                 {
                     delete[] p[i].element[g];
                 }
		          delete[] p[i].element;



            }

		    p=t;
        }
    }

};


//matrix* memory =new matrix[20];
int tcntr=0;
mDynArr memory (10);
int memoryPointer=0;
int exit1=0; //exit is ambiguos to VS so i added the 1 (it's meaningless)

void removeSpaces(string &str)
{
    for(int i=0; i<str.length(); i++)
     {
         if(str[i] == ' ') {str.erase(i,1);i--;}//replaced erase fn with =''

     }
}


int memoryCheck(string mName)
{
    for(int i=0;i<memoryPointer;i++)
    {
        if(memory[i].getName()==mName)
            return i;
    }
    return -1;

}
void memorizeMatrix(int mIndex, string mString,string mName)
{
    if(mIndex==-1)
    {
        //cout<<"mwmory"<<memoryPointer;
        memory.create(mName,mString);
        memory.p[memoryPointer].print();
        memoryPointer++;
    }

    else
    {
        memory.update(mIndex,mName,mString);
        memory.p[mIndex].print();
    }
}
void cut(string &variable1,string &variable2,int &index1,int &index2,char op,string operation)
{
    variable1= operation.substr(0,operation.find(op)) ;
    variable2= operation.substr(operation.find(op)+1,(operation.length()-operation.find(op))-1);

    index1=memoryCheck(variable1);
    index2=memoryCheck(variable2);
}
void cut(string &variable1, int &index1, char op, string operation)
{
	variable1 = operation.substr(operation.find(op) + 1, (operation.length() - operation.find(op)) - 1);
	index1 = memoryCheck(variable1);
}
void input_checker(string input) // assignment or operation
{

	short asg=0;
	//FOB=First Open Bracket, FCB=First Close Bracker, EQ=EQual
	short FOBPos = input.find('['), FCBPos = input.find(']'), LCBPos = input.rfind(']'), EQPos = input.find('='), plusPos = input.find('+')
		,minusPos= input.find('-'),elemWiseInvPos= input.find("./"),divPos= input.find('/'),multPos= input.find("*")
		,transPos= input.find("'");
    string mName=input.substr(0, EQPos);
	removeSpaces(mName);

    int index = memoryCheck(mName);
	bool assignmentOP = (FOBPos != -1) && (FCBPos != -1);
	bool mathOP = (plusPos != -1) || (minusPos != -1) || (divPos != -1) || (multPos != -1) || (transPos != -1);
    if(assignmentOP)
    {
		string mString = input.substr(FOBPos+1, (LCBPos - 2 - FOBPos + 1));//FOB+1 & LCB-2 to remove braces
		if (!mathOP && mString.find("[")==-1 )//not complete yet. We need to check for letters
		{
			memorizeMatrix(index, mString, mName); //I think we may handle sin() & 1X1 matrix inside this func
		}
		else
		{

		}
    }

    else if(mathOP)
    {
		    string operation = input.substr(EQPos +1,(input.length()- EQPos +1));
            removeSpaces(operation);
            string variable1,variable2;
            int index1,index2;


				//memory.push_back(temp);
			if(plusPos !=-1)
			{
			    cut(variable1,variable2,index1,index2,'+',operation);


                if(index!=-1)
                {
                asg=(index1==index || index2==index);
                memory.p[index].add(memory.p[index1],memory.p[index2],asg) ;
                memory.p[index].print();
                }

                else
                {
                memory.create(mName);
                memory.p[memoryPointer].add(memory.p[index1],(memory.p[index2])) ;
                memory.p[memoryPointer].print();
                memoryPointer++;
                }
            }

			else if(minusPos !=-1)
			{
			    cut(variable1,variable2,index1,index2,'-',operation);

                if(index!=-1)
                {
                asg=(index1==index || index2==index);
                memory.p[index].sub(memory.p[index1],(memory.p[index2]),asg);
                memory.p[index].print();
                }

                else
                {
                memory.create(mName);
                memory.p[memoryPointer].sub(memory.p[index1],(memory.p[index2])) ;
                memory.p[memoryPointer].print();
                memoryPointer++;
                }
			}
			else if(multPos !=-1)
			{
			    cut(variable1,variable2,index1,index2,'*',operation);

                if(index!=-1)
                {
                    asg=(index1==index);
                    if( index2==index ) asg=2;
                    memory.p[index].mult(memory.p[index1],(memory.p[index2]),asg) ;
                    memory.p[index].print();
                }

                else
                {
                memory.create(mName);
                memory.p[memoryPointer].mult(memory.p[index1],(memory.p[index2])) ;
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
			else if(divPos !=-1)
			{
			    cut(variable1,variable2,index1,index2,'/',operation);
//cout<<memory.p[index2].getDeterminant()<<endl;
				if (memory.p[index2].getDeterminant() == 0)
				{
					cout << "Division cannot done" << endl;
					return;
				}

                else if(index!=-1)
                {
                memory.p[index].div(memory.p[index1],(memory.p[index2])) ;
                memory.p[index].print();
                }

                else
                {
                memory.create(mName);
                memory.p[memoryPointer].div(memory.p[index1],(memory.p[index2])) ;
                memory.p[memoryPointer].print();
                memoryPointer++;
                }
            }
			else if(transPos !=-1)
			{
			    string var;
			    int ind;
			    var= operation.substr(0, operation.find("'")) ;
                ind=memoryCheck(var);


                if(index!=-1)
                {
                memory.p[index].getTranspose(memory.p[ind]) ;
                memory.p[index].print();
                }

                else
                {
                memory.create(mName);
                memory.p[memoryPointer].getTranspose(memory.p[ind]) ;
                memory.p[memoryPointer].print();
                memoryPointer++;
                }
			}
}

		else
        {
            removeSpaces(input);
             if((memoryCheck(input)!=-1) &&
                 memoryCheck(input)<memoryPointer)
             {
                 memory.p[memoryCheck(input)].print();
             }
             else if(input=="exit")
             {
                 exit1=1;
                 return;
             }
             else
            cout<<"invalid input"<<endl ;
        }


}
//==========================================================================
struct size { int rows,columns; };   //structure has the number of rows&columns for any string matrix
//===========================================================================
vector<string> separate(string inputString)
{
    vector<string> separatedString;
    int flag=0;
    int beginPostion=0;
    for(int i=0;i<inputString.length();i++)
    {
        if(inputString[i]==';'&&flag==0)
        {
            cout<<inputString.substr(beginPostion,i-beginPostion)<<endl;
            separatedString.push_back(inputString.substr(beginPostion,i-beginPostion));
            beginPostion=i+1;
        }
        else if(i==inputString.length()-1)
        {
            cout<<inputString.substr(beginPostion,inputString.length()-beginPostion);
            separatedString.push_back(inputString.substr(beginPostion,inputString.length()-beginPostion));
        }
        else if(inputString[i]=='[')
            flag++;
        else if(inputString[i]==']')
            flag--;
    }
    return separatedString;
}
//===========================================================================
size compare(size m1, size m2)     //get the total size of 2 concatenated matrices
{
    size m;
    if(m1.rows==m2.rows)   //if the rows of matrix1 = rows of matrix2 whatever the columns are equal or not
    {
        m.rows=m1.rows;
        m.columns=m1.columns+m2.columns;
        return m;
    }
    else if(m1.columns==m2.columns&&m1.rows!=m2.rows)  //if thr columns of matrix1= columns of matrix2 and the rows aren't equal
    {
        m.columns=m1.columns;
        m.rows=m1.rows+m2.rows;
        return m;
    }
}
//=============================================================================
size sizing(string matrix)  //take any string but without brackets and retur a structure which has the number of rows and columns
{                               // for example    22 45 sin(90)+1;1 cos(0) 5    size will be 2x3
	size n; n.rows=0;  n.columns=0;
	 stringstream ss(matrix);
     string token;
     while(getline(ss, token, ';'))
		 n.rows++;
	 getline(ss, token, ';');
	 stringstream sn(token);
	 string in;
     while( sn >> in)
		 n.columns++;
	 return n;
}
//=================================================================================
#define endl '\n'
/* Omar Yasser */


bool Is_operation(char character) {
	if(character=='('||character==')'||character=='^'||character=='*'||
    character=='/'||character=='+'||character=='-' ||character=='~')
		return 1;
	else return 0;
}

string alphanum =
"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
string symb="@#$_";
int stringLength = alphanum.length();
int Mctr=0;
string genRandom()  // Random string generator function.
{
    if(Mctr>=stringLength)
    {
        symb.erase(0,1);
        alphanum =
         "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
         Mctr=0;
    }
    string ranS="_a";
    ranS[1]=alphanum[0];
    ranS[0]=symb[0];
    alphanum.erase(0,1);
    Mctr++;
    return ranS;
}
void calcAndRep(int i, int j,string  &fullOp,char ch_op){
	//cout<<i<<" "<<j << " "<<fullOp<<" "<<ch_op<<endl ;

    int opOnNum=1;
	 int pos1=i,pos2=j;
    double op1=0,op2=0 ,res=0;

    string s_op1,s_op2;

     stringstream strm,strm2;
     stringstream strm3;
     string result;
    string op=fullOp.substr(pos1,pos2-pos1+1);


    s_op1=op.substr(0,op.find(ch_op,1));
    s_op2=op.substr(op.find(ch_op,1)+1,pos2-op.find(ch_op,1)+1);
    if(ch_op=='~')
        s_op2.erase(0,2);

    //cout<<op<<" opers "<<s_op1<<" "<<s_op2<<endl;
    if(op.find("sqrt")==0)
    {
       s_op1=op.substr(4,op.length()-4);
        if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].sqrtMatrix(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       strm<<s_op1;
       strm>> op1;
       res=sqrt(op1);
      // cout<<s_op1<<endl;
    }
     else if(op.find("~P~")<op.length())
    {
        if(memoryCheck(s_op2)!=-1)
        {

            strm<<s_op1;
            strm>> op1;
            string nMat=genRandom();
           memory.create(nMat);

           memory.p[memoryPointer].addEL(memory.p[memoryCheck(s_op2)],op1);
           opOnNum=0;
           memoryPointer++;
          result=nMat;
        }
        if(memoryCheck(s_op1)!=-1)
        {
            strm<<s_op2;
            strm>> op2;
            string nMat=genRandom();
           memory.create(nMat);

           memory.p[memoryPointer].addEL(memory.p[memoryCheck(s_op1)],op2);
           opOnNum=0;
           memoryPointer++;
          result=nMat;
        }
    }
     else if(op.find("~S~")<op.length())
    {
        if(memoryCheck(s_op2)!=-1)
        {

            strm<<s_op1;
            strm>> op1;
            string nMat=genRandom();
           memory.create(nMat);

           memory.p[memoryPointer].subEL(memory.p[memoryCheck(s_op2)],op1);
           opOnNum=0;
           memoryPointer++;
          result=nMat;
        }
        if(memoryCheck(s_op1)!=-1)
        {
            strm<<s_op2;
            strm>> op2;
            string nMat=genRandom();
           memory.create(nMat);

           memory.p[memoryPointer].subEL(memory.p[memoryCheck(s_op1)],op2);
           opOnNum=0;
           memoryPointer++;
          result=nMat;
        }
    }
         else if(op.find("~D~")<op.length())
    {
        if((memoryCheck(s_op1)!=-1) && memoryCheck(s_op2)!=-1)
        {
            string nMat=genRandom();
           memory.create(nMat);

           memory.p[memoryPointer].divEL(memory.p[memoryCheck(s_op1)],
                                         memory.p[memoryCheck(s_op2)]);

           opOnNum=0;
           memoryPointer++;
          result=nMat;
        }
        else if(memoryCheck(s_op2)!=-1)
        {
            strm<<s_op1;
            strm>> op1;
            string nMat=genRandom();
           memory.create(nMat);

           memory.p[memoryPointer].divEL(memory.p[memoryCheck(s_op2)],op1);

          memory.p[memoryPointer].elementWisePower(memory.p[memoryPointer],-1);
           opOnNum=0;
           memoryPointer++;
          result=nMat;
        }
        else if(memoryCheck(s_op1)!=-1)
        {
            strm<<s_op2;
            strm>> op2;
            string nMat=genRandom();
           memory.create(nMat);

           memory.p[memoryPointer].divEL(memory.p[memoryCheck(s_op1)],op2);
           opOnNum=0;
           memoryPointer++;
          result=nMat;
        }
    }
    else if(op.find("~M~")<op.length())
    {
         if((memoryCheck(s_op1)!=-1) && memoryCheck(s_op2)!=-1)
        {
            string nMat=genRandom();
           memory.create(nMat);

           memory.p[memoryPointer].multEL(memory.p[memoryCheck(s_op1)],
                                         memory.p[memoryCheck(s_op2)]);
           opOnNum=0;
           memoryPointer++;
          result=nMat;
        }
       else if(memoryCheck(s_op2)!=-1)
        {

            strm<<s_op1;
            strm>> op1;
            string nMat=genRandom();
           memory.create(nMat);

           memory.p[memoryPointer].multEL(memory.p[memoryCheck(s_op2)],op1);
           opOnNum=0;
           memoryPointer++;
          result=nMat;
        }
        else if(memoryCheck(s_op1)!=-1)
        {
            strm<<s_op2;
            strm>> op2;
            string nMat=genRandom();
           memory.create(nMat);

           memory.p[memoryPointer].multEL(memory.p[memoryCheck(s_op1)],op2);
           opOnNum=0;
           memoryPointer++;
          result=nMat;
        }
    }
        else if(op.find("~W~")<op.length())
    {

        if(memoryCheck(s_op1)!=-1)
        {
            strm<<s_op2;
            strm>> op2;
            string nMat=genRandom();
           memory.create(nMat);

           memory.p[memoryPointer].elementWisePower(memory.p[memoryCheck(s_op1)],op2);
           opOnNum=0;
           memoryPointer++;
          result=nMat;
        }
    }
     else if(op.find("exp")==0)
    {
       s_op1=op.substr(3,op.length()-3);
       if(memoryCheck(s_op1)!=-1)
       {

           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].expMatrix(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       {
         strm<<s_op1;
         strm>> op1;
           res=exp(op1);
       }
    }
    else if(op.find("log")==0)
    {
       s_op1=op.substr(3,op.length()-3);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].logMatrix(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       {
         strm<<s_op1;
         strm>> op1;
           res=log(op1);
       }
    }
    else if(op.find("tog10")==0)
    {
       s_op1=op.substr(5,op.length()-5);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].log10Matrix(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       {
         strm<<s_op1;
         strm>> op1;
           res=log10(op1);
       }
    }
    else if(op.find("sinh")==0)
    {
       s_op1=op.substr(4,op.length()-4);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Msinh(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       strm<<s_op1;
       strm>> op1;
       res=sinh(op1);
    }
    else if(op.find("sind")==0)
    {
       s_op1=op.substr(4,op.length()-4);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Msind(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       strm<<s_op1;
       strm>> op1;
       res=sind(op1);
    }
    else if(op.find("sin")==0)
    {
       s_op1=op.substr(3,op.length()-3);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Msin(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       {
         strm<<s_op1;
         strm>> op1;
           res=sin(op1);
       }
    }

	else if(op.find("asind")==0)
    {
       s_op1=op.substr(5,op.length()-5);
      if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Masind(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
         strm<<s_op1;
         strm>> op1;
           res=sin(op1);
       }
    }
	else if(op.find("asinh")==0)
    {
       s_op1=op.substr(5,op.length()-5);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Masinh(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
         strm<<s_op1;
         strm>> op1;
           res=sin(op1);
       }
    }
    else if(op.find("asin")==0)
    {

       s_op1=op.substr(4,op.length()-4);
       if(memoryCheck(s_op1)!=-1)
       {

           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Masin(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
         strm<<s_op1;
         strm>> op1;
           res=sin(op1);
       }
    }

	else if(op.find("cosd")==0)
    {
       s_op1=op.substr(4,op.length()-4);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Mcosd(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
         strm<<s_op1;
         strm>> op1;
           res=sin(op1);
       }
    }
	else if(op.find("cosh")==0)
    {
       s_op1=op.substr(4,op.length()-4);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Mcosh(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
         strm<<s_op1;
         strm>> op1;
           res=sin(op1);
       }
    }
     else if(op.find("cos")==0)
    {
       s_op1=op.substr(3,op.length()-3);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Mcos(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
         strm<<s_op1;
         strm>> op1;
           res=sin(op1);
       }
    }

	else if(op.find("acosd")==0)
    {
       s_op1=op.substr(5,op.length()-5);
       strm<<s_op1;
       strm>> op1;
       res=acosd(op1);
      // cout<<s_op1<<endl;
    }
	else if(op.find("acosh")==0)
    {
       s_op1=op.substr(5,op.length()-5);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Macosh(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
         strm<<s_op1;
         strm>> op1;
           res=sin(op1);
       }
    }
    else if(op.find("acos")==0)
    {
       s_op1=op.substr(4,op.length()-4);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Macos(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
         strm<<s_op1;
         strm>> op1;
           res=sin(op1);
       }
    }

   else if(op.find("tand")==0)
    {
       s_op1=op.substr(4,op.length()-4);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Mtand(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
         strm<<s_op1;
         strm>> op1;
           res=sin(op1);
       }
    }
	else if(op.find("tanh")==0)
    {
       s_op1=op.substr(4,op.length()-4);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Mtanh(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
         strm<<s_op1;
         strm>> op1;
           res=sin(op1);
       }
    }
    else if(op.find("tan")==0)
    {
       s_op1=op.substr(3,op.length()-3);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Mtan(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
         strm<<s_op1;
         strm>> op1;
           res=sin(op1);
       }
    }

	else if(op.find("atand")==0)
    {
       s_op1=op.substr(5,op.length()-5);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Matand(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
         strm<<s_op1;
         strm>> op1;
           res=sin(op1);
       }
    }
	else if(op.find("atanh")==0)
    {
       s_op1=op.substr(5,op.length()-5);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Matanh(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
         strm<<s_op1;
         strm>> op1;
           res=sin(op1);
       }
    }

	else if(op.find("atan2d")==0) //need some work
    {
       s_op1=op.substr(6,op.length()-6);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Matand(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
       strm<<s_op1;
       strm>> op1;
       res=atanh(op1);
       }
    }
    else if(op.find("atan2")==0) //need some work
    {
       s_op1=op.substr(5,op.length()-5);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Matanh(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
       strm<<s_op1;
       strm>> op1;
       res=atanh(op1);
       }
    }
    else if(op.find("atan")==0)
    {
       s_op1=op.substr(4,op.length()-4);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Matan(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
       strm<<s_op1;
       strm>> op1;
       res=atan(op1);
       }
    }

	else if(op.find("cscd")==0)
    {
       s_op1=op.substr(4,op.length()-4);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Mcscd(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
       strm<<s_op1;
       strm>> op1;
       res=cscd(op1);
       }
    }
	else if(op.find("csch")==0)
    {
       s_op1=op.substr(4,op.length()-4);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Mcsch(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
       strm<<s_op1;
       strm>> op1;
       res=csch(op1);
       }
    }
    else if(op.find("csc")==0)
    {
       s_op1=op.substr(3,op.length()-3);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Mcsc(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
       strm<<s_op1;
       strm>> op1;
       res=csc(op1);
       }
    }

	else if(op.find("acscd")==0)
    {
       s_op1=op.substr(5,op.length()-5);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Macscd(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
       strm<<s_op1;
       strm>> op1;
       res=acscd(op1);
       }
    }
	else if(op.find("acsch")==0)
    {
       s_op1=op.substr(5,op.length()-5);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Macsch(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
       strm<<s_op1;
       strm>> op1;
       res=acsch(op1);
       }
    }
    else if(op.find("acsc")==0)
    {
       s_op1=op.substr(4,op.length()-4);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Macsc(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
       strm<<s_op1;
       strm>> op1;
       res=acsc(op1);
       }
    }

	else if(op.find("secd")==0)
    {
       s_op1=op.substr(4,op.length()-4);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Msecd(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
       strm<<s_op1;
       strm>> op1;
       res=secd(op1);
       }
    }
	else if(op.find("sech")==0)
    {
       s_op1=op.substr(4,op.length()-4);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Msech(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
       strm<<s_op1;
       strm>> op1;
       res=sech(op1);
       }
    }
    else if(op.find("sec")==0)
    {
       s_op1=op.substr(3,op.length()-3);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Msec(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
       strm<<s_op1;
       strm>> op1;
       res=sec(op1);
       }
    }

	else if(op.find("asecd")==0)
    {
       s_op1=op.substr(5,op.length()-5);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Masecd(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
       strm<<s_op1;
       strm>> op1;
       res=asecd(op1);
       }
    }
	else if(op.find("asech")==0)
    {
       s_op1=op.substr(5,op.length()-5);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Masech(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
       strm<<s_op1;
       strm>> op1;
       res=asech(op1);
       }
    }
    else if(op.find("asec")==0)
    {
       s_op1=op.substr(4,op.length()-4);
       if(memoryCheck(s_op1)!=-1)
       {
           string nMat=genRandom();
           memory.create(nMat);
           memory.p[memoryPointer].Masec(memory.p[memoryCheck(s_op1)]);
           opOnNum=0;
           memoryPointer++;
           result=nMat;
       }
       else
       {
       strm<<s_op1;
       strm>> op1;
       res=asec(op1);
      // cout<<"GLLLLLLL"<<" "<<s_op1<< " "<< res<<endl;
       }
    }
	else if(op.find("hypot")==0) //need some work
    {
       s_op1=op.substr(5,op.length()-5);
       strm<<s_op1;
       strm>> op1;
       res=atanh(op1);
       //cout<<s_op1<<endl;
    }
	else if(op.find("deg2rad")==0) //need some work
    {
       s_op1=op.substr(5,op.length()-5);
       strm<<s_op1;
       strm>> op1;
       res=atanh(op1);
      // cout<<s_op1<<endl;
    }
	else if(op.find("rad2deg")==0) //need some work
    {
       s_op1=op.substr(5,op.length()-5);
       strm<<s_op1;
       strm>> op1;
       res=atanh(op1);
      // cout<<s_op1<<endl;
    }
    else
    {
    strm << op;
    strm >> op1 ;
    op=op.substr(op.find(ch_op,1),op.length()-op.find(ch_op,1));
    ch_op=op[0];
    op.erase(0,1);
    strm2 << op;
    strm2 >> op2;
    //cout << op1 <<"/-/ "<<op2 <<endl<<op<<endl;

    switch (ch_op)
       {
       case '^':

     if(memoryCheck(s_op1)!=-1 && memoryCheck(s_op2)!=-1)
    {

         if(memory.p[memoryCheck(s_op2)].getColumns()==1&&
            memory.p[memoryCheck(s_op2)].getRows()==1)
         {
        string nMat=genRandom();
         memory.create(nMat);
           memory.p[memoryPointer].power(memory.p[memoryCheck(s_op1)],
                                memory.p[memoryCheck(s_op2)].element[0][0]);
           opOnNum=0;
           memoryPointer++;
        result=nMat;

         }
         else break;
    }
    else if(memoryCheck(s_op1)!=-1)
    {

        string nMat=genRandom();
         memory.create(nMat);
           memory.p[memoryPointer].power(memory.p[memoryCheck(s_op1)],op2);
           opOnNum=0;
           memoryPointer++;
        result=nMat;
    }
    res = pow(op1,op2);
      break;
  case '+':
    if(memoryCheck(s_op1)!=-1 && memoryCheck(s_op2)!=-1)
    {

        string nMat=genRandom();
         memory.create(nMat);
           memory.p[memoryPointer].add(memory.p[memoryCheck(s_op1)],
                                memory.p[memoryCheck(s_op2)]);
           opOnNum=0;
           memoryPointer++;
        result=nMat;
    }
    else
  res = op1 + op2;
  break;
case '-':
    if(memoryCheck(s_op1)!=-1 && memoryCheck(s_op2)!=-1)
    {
        string nMat=genRandom();
         memory.create(nMat);
           memory.p[memoryPointer].sub(memory.p[memoryCheck(s_op1)],
                                memory.p[memoryCheck(s_op2)]);
           opOnNum=0;
           memoryPointer++;
        result=nMat;
    }
    else
  res = op1 - op2;
  break;
case '*':
        if(memoryCheck(s_op1)!=-1 && memoryCheck(s_op2)!=-1)
    {
        string nMat=genRandom();
         memory.create(nMat);
           memory.p[memoryPointer].mult(memory.p[memoryCheck(s_op1)],
                                memory.p[memoryCheck(s_op2)]);
           opOnNum=0;
           memoryPointer++;
        result=nMat;
    }
    else
  res = op1 * op2;
  break;
  case '/':
   if(memoryCheck(s_op1)!=-1 && memoryCheck(s_op2)!=-1)
     {
        string nMat=genRandom();
         memory.create(nMat);
           memory.p[memoryPointer].div(memory.p[memoryCheck(s_op1)],
                                memory.p[memoryCheck(s_op2)]);
           opOnNum=0;
           memoryPointer++;
        result=nMat;
     }
    else
  res = op1 / op2;
  break;
}
    }
//if((abs(res-round(res))<=0.00001 && abs(res-round(res))>=-0.00001)||
  //     (abs(round(res)-res )<=0.00001 && abs(round(res) -res)>=-0.00001) )
    //    res=round(res);
if(opOnNum)
{
    strm3 << res;
    strm3 >> result;
}

fullOp.replace(pos1,pos2-pos1+1,result);
//cout<<fullOp<<endl;
}

void Operation_solver(string &operation){
	removeSpaces(operation) ;
	if(operation.find("log10")!=-1){
        int rep= operation.find("log10");
        operation[rep]='t' ;
	}
	if(operation.find(".*")!=-1){
        int rep= operation.find(".*");
        operation.replace(rep,2,"~M~");

	}
	if(operation.find("./")!=-1){
        int rep= operation.find("./");
        operation.replace(rep,2,"~D~");
	}
	if(operation.find(".+")!=-1){
        int rep= operation.find(".+");
        operation.replace(rep,2,"~P~");
	}
	if(operation.find(".-")!=-1){
        int rep= operation.find(".-");
        operation.replace(rep,2,"~S~");
	}
	if(operation.find(".^")!=-1){
        int rep= operation.find(".^");
        operation.replace(rep,2,"~W~");
	}

	for(int i=operation.length()-1 ; i>=0 ; i--){
		if(operation[i]=='('){
			int j=0;

			for(j=i ; j<operation.length(); j++)
				if(operation[j]==')')
			break;

			string temp=operation.substr(i+1,j-i-1);
			//cout<<"Temp "<<temp<<endl;
                Operation_solver(temp);
		    	operation.replace(i,j-i+1,temp);
                Operation_solver(operation) ;

			break;
		 }
		}

		for(int i=operation.length()-1 ; i>=0 ; i--){ //sqrt
			if(operation[i]=='s'&&operation[i+1]=='q'&&operation[i+2]=='r'
                &&operation[i+3]=='t'){
				int end=0 ;
				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}

		for(int i=operation.length()-1 ; i>=0 ; i--){ //asind
			if(operation[i]=='a'&&operation[i+1]=='s'&&operation[i+2]=='i'
                &&operation[i+3]=='n'&&operation[i+4]=='d'){
				int end=0 ;
				for( end=i+5 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){ //asinh
			if(operation[i]=='a'&&operation[i+1]=='s'&&operation[i+2]=='i'
                &&operation[i+3]=='n'&&operation[i+4]=='h'){
				int end=0 ;
				for( end=i+5 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){ //asin
			if(operation[i]=='a'&&operation[i+1]=='s'&&operation[i+2]=='i'
                &&operation[i+3]=='n'){
				int end=0 ;
				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;

				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){ //sinh
			if(operation[i]=='s'&&operation[i+1]=='i'&&operation[i+2]=='n'
                &&operation[i+3]=='h'){
				int end=0 ;
				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){ //sind
			if(operation[i]=='s'&&operation[i+1]=='i'&&operation[i+2]=='n'
                &&operation[i+3]=='d'){
				int end=0 ;
				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){ //sin
			if(operation[i]=='s'&&operation[i+1]=='i'&&operation[i+2]=='n'){
				int end=0 ;
				for( end=i+3 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
              //  cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}

		for(int i=operation.length()-1 ; i>=0 ; i--){ //acosd
			if(operation[i]=='a'&&operation[i+1]=='c'&&operation[i+2]=='o'
                &&operation[i+3]=='s'&&operation[i+4]=='d'){
				int end=0 ;
				for( end=i+5 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){ //acosh
			if(operation[i]=='a'&&operation[i+1]=='c'&&operation[i+2]=='o'
                &&operation[i+3]=='s'&&operation[i+4]=='h'){
				int end=0 ;
				for( end=i+5 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){//acos
			if(operation[i]=='a'&&operation[i+1]=='c'&&operation[i+2]=='o'
                &&operation[i+3]=='s'){
				int end=0 ;
				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		  }
		  for(int i=operation.length()-1 ; i>=0 ; i--){ //cosh
			if(operation[i]=='c'&&operation[i+1]=='o'&&operation[i+2]=='s'
                &&operation[i+3]=='h'){
				int end=0 ;
				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){ //cosd
			if(operation[i]=='c'&&operation[i+1]=='o'&&operation[i+2]=='s'
                &&operation[i+3]=='d'){
				int end=0 ;
				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		   for(int i=operation.length()-1 ; i>=0 ; i--){ //cos
			if(operation[i]=='c'&&operation[i+1]=='o'&&operation[i+2]=='s'){
				int end=0 ;
				for( end=i+3 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}



		for(int i=operation.length()-1 ; i>=0 ; i--){ //atand
			if(operation[i]=='a'&&operation[i+1]=='t'&&operation[i+2]=='a'
                &&operation[i+3]=='n'&&operation[i+4]=='d'){
				int end=0 ;
				for( end=i+5 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){ //atanh
			if(operation[i]=='a'&&operation[i+1]=='t'&&operation[i+2]=='a'
                &&operation[i+3]=='n'&&operation[i+4]=='h'){
				int end=0 ;
				for( end=i+5 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){ //atan2
			if(operation[i]=='a'&&operation[i+1]=='t'&&operation[i+2]=='a'
                &&operation[i+3]=='n'&&operation[i+4]=='2'){
				int end=0 ;
				for( end=i+5 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){ //atand
			if(operation[i]=='a'&&operation[i+1]=='t'&&operation[i+2]=='a'
                &&operation[i+3]=='n'&&operation[i+4]=='2'&&operation[i+4]=='d'){
				int end=0 ;
				for( end=i+6 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){ //atan
			if(operation[i]=='a'&&operation[i+1]=='t'&&operation[i+2]=='a'
                &&operation[i+3]=='n'){
				int end=0 ;
				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}

			for(int i=operation.length()-1 ; i>=0 ; i--){ //tanh
			if(operation[i]=='t'&&operation[i+1]=='a'&&operation[i+2]=='n'
                &&operation[i+3]=='h'){
				int end=0 ;
				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
       for(int i=operation.length()-1 ; i>=0 ; i--){ //tand
			if(operation[i]=='t'&&operation[i+1]=='a'&&operation[i+2]=='n'
                &&operation[i+3]=='d'){
				int end=0 ;
				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){ //tan
			if(operation[i]=='t'&&operation[i+1]=='a'&&operation[i+2]=='n'){
				int end=0 ;
				for( end=i+3 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}


for(int i=operation.length()-1 ; i>=0 ; i--){ //acscd
			if(operation[i]=='a'&&operation[i+1]=='c'&&operation[i+2]=='s'
                &&operation[i+3]=='c'&&operation[i+4]=='d'){
				int end=0 ;
				for( end=i+5 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
for(int i=operation.length()-1 ; i>=0 ; i--){ //acsch
			if(operation[i]=='a'&&operation[i+1]=='c'&&operation[i+2]=='s'
                &&operation[i+3]=='c'&&operation[i+4]=='h'){
				int end=0 ;
				for( end=i+5 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
			for(int i=operation.length()-1 ; i>=0 ; i--){ //acsc
			if(operation[i]=='a'&&operation[i+1]=='c'&&operation[i+2]=='s'
                &&operation[i+3]=='c'){
				int end=0 ;
				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}

		for(int i=operation.length()-1 ; i>=0 ; i--){ //csch
			if(operation[i]=='c'&&operation[i+1]=='s'&&operation[i+2]=='c'
                &&operation[i+3]=='h'){
				int end=0 ;
				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
	for(int i=operation.length()-1 ; i>=0 ; i--){ //cscd
			if(operation[i]=='c'&&operation[i+1]=='s'&&operation[i+2]=='c'
                &&operation[i+3]=='d'){
				int end=0 ;
				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){ //csc
			if(operation[i]=='c'&&operation[i+1]=='s'&&operation[i+2]=='c'){
				int end=0 ;
				for( end=i+3 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}



for(int i=operation.length()-1 ; i>=0 ; i--){  //asecd
			if(operation[i]=='a'&&operation[i+1]=='s'&&operation[i+2]=='e'
                &&operation[i+3]=='c'&&operation[i+4]=='d'){
				int end=0 ;
				for( end=i+5 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
for(int i=operation.length()-1 ; i>=0 ; i--){  //asech
			if(operation[i]=='a'&&operation[i+1]=='s'&&operation[i+2]=='e'
                &&operation[i+3]=='c'&&operation[i+4]=='h'){
				int end=0 ;
				for( end=i+5 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}

		for(int i=operation.length()-1 ; i>=0 ; i--){  //asec
			if(operation[i]=='a'&&operation[i+1]=='s'&&operation[i+2]=='e'
                &&operation[i+3]=='c'){
				int end=0 ;
				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){ //sech
			if(operation[i]=='s'&&operation[i+1]=='e'&&operation[i+2]=='c'
                &&operation[i+3]=='h'){
				int end=0 ;
				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
	for(int i=operation.length()-1 ; i>=0 ; i--){ //secd
			if(operation[i]=='s'&&operation[i+1]=='e'&&operation[i+2]=='c'
                &&operation[i+3]=='d'){
				int end=0 ;
				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){ //sec
			if(operation[i]=='s'&&operation[i+1]=='e'&&operation[i+2]=='c'){
				int end=0 ;
				for( end=i+3 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}

for(int i=operation.length()-1 ; i>=0 ; i--){  //acotd
			if(operation[i]=='a'&&operation[i+1]=='c'&&operation[i+2]=='o'
                &&operation[i+3]=='t'&&operation[i+4]=='d'){
				int end=0 ;
				for( end=i+5 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
for(int i=operation.length()-1 ; i>=0 ; i--){// acoth
			if(operation[i]=='a'&&operation[i+1]=='c'&&operation[i+2]=='o'
                &&operation[i+3]=='t'&&operation[i+4]=='h'){
				int end=0 ;
				for( end=i+5 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){ //acot
			if(operation[i]=='a'&&operation[i+1]=='c'&&operation[i+2]=='o'
                &&operation[i+3]=='t'){
				int end=0 ;
				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){ //cot
			if(operation[i]=='c'&&operation[i+1]=='o'&&operation[i+2]=='t'){
				int end=0 ;
				for( end=i+3 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
for(int i=operation.length()-1 ; i>=0 ; i--){ //hypot
			if(operation[i]=='h'&&operation[i+1]=='y'&&operation[i+2]=='p'
                &&operation[i+3]=='o'&&operation[i+4]=='t'){
				int end=0 ;
				for( end=i+5 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
for(int i=operation.length()-1 ; i>=0 ; i--){ //deg2rad
			if(operation[i]=='d'&&operation[i+1]=='e'&&operation[i+2]=='g'
                &&operation[i+3]=='2'&&operation[i+4]=='r'&&operation[i+5]=='a'&&operation[i+6]=='d'){
				int end=0 ;
				for( end=i+7 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
for(int i=operation.length()-1 ; i>=0 ; i--){ //rad2deg
			if(operation[i]=='r'&&operation[i+1]=='a'&&operation[i+2]=='d'
                &&operation[i+3]=='2'&&operation[i+4]=='d'&&operation[i+5]=='e'&&operation[i+6]=='g'){
				int end=0 ;
				for( end=i+7 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){
			if(operation[i]=='r'&&operation[i+1]=='a'&&operation[i+2]=='n'&&operation[i+3]=='d'){
				int start=0;

				for(start=i-1 ; start>0 ; start--)
					if(Is_operation(operation[start-1])&&start!=1)
						break;
				int end=0;

				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);

               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){
			if(operation[i]=='e'&&operation[i+1]=='y'&&operation[i+2]=='e'){
				int start=0;

				for(start=i-1 ; start>0 ; start--)
					if(Is_operation(operation[start-1])&&start!=1)
						break;
				int end=0;

				for( end=i+3 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){
			if(operation[i]=='z'&&operation[i+1]=='e'&&operation[i+2]=='r'&&operation[i+3]=='o'&&operation[i+4]=='s'){
				int start=0;

				for(start=i-1 ; start>0 ; start--)
					if(Is_operation(operation[start-1])&&start!=1)
						break;
				int end=0;

				for( end=i+5 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
for(int i=operation.length()-1 ; i>=0 ; i--){
			if(operation[i]=='o'&&operation[i+1]=='n'&&operation[i+2]=='e'&&operation[i+3]=='s'){
				int start=0;

				for(start=i-1 ; start>0 ; start--)
					if(Is_operation(operation[start-1])&&start!=1)
						break;
				int end=0;

				for( end=i+4 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);

                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}

		for(int i=operation.length()-1 ; i>=0 ; i--){
			if(operation[i]=='t'&&operation[i+1]=='o'&&operation[i+2]=='g'&&operation[i+3]=='1'&&operation[i+4]=='0'){
				int start=0;

				for(start=i-1 ; start>0 ; start--)
					if(Is_operation(operation[start-1])&&start!=1)
						break;
				int end=0;

				for( end=i+5 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){
			if(operation[i]=='l'&&operation[i+1]=='o'&&operation[i+2]=='g'){
				int start=0;

				for(start=i-1 ; start>0 ; start--)
					if(Is_operation(operation[start-1])&&start!=1)
						break;
				int end=0;

				for( end=i+3 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){
			if(operation[i]=='e'&&operation[i+1]=='x'&&operation[i+2]=='p'){
				int start=0;

				for(start=i-1 ; start>0 ; start--)
					if(Is_operation(operation[start-1])&&start!=1)
						break;
				int end=0;

				for( end=i+3 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(i,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}


		for(int i=operation.length()-1 ; i>=0 ; i--){
			if(operation[i]=='^'){
				int start=0;

				for(start=i-1 ; start>0 ; start--)
					if(Is_operation(operation[start-1])&&start!=1)
						break;
				int end=0;

				for( end=i+1 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(start,end,operation,operation[i]);

               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}

for(int i=0;  i<operation.length()-1 ; i++){
			if(operation[i]=='*'){
				int start=0;

				for(start=i-1 ; start>0 ; start--)
					if(Is_operation(operation[start-1])&&start!=1)
						break;
				int end=0;

				for( end=i+1 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(start,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}


		for(int i=0;  i<operation.length()-1 ; i++){
			if(operation[i]=='/'){
				int start=0;

				for(start=i-1 ; start>0 ; start--)
					if(Is_operation(operation[start-1])&&start!=1)
						break;
				int end=0;

				for( end=i+1 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(start,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=0;  i<operation.length()-1 ; i++){
			if(operation[i]=='~'&&operation[i+1]=='W'&&operation[i+2]=='~'){
				int start=0;

				for(start=i-1 ; start>0 ; start--)
					if(Is_operation(operation[start-1])&&start!=1)
						break;
				int end=0;

				for( end=i+3 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(start,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=0;  i<operation.length()-1 ; i++){
			if(operation[i]=='~'&&operation[i+1]=='M'&&operation[i+2]=='~'){
				int start=0;

				for(start=i-1 ; start>0 ; start--)
					if(Is_operation(operation[start-1])&&start!=1)
						break;
				int end=0;

				for( end=i+3 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(start,end,operation,operation[i]);

              //  cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=0;  i<operation.length()-1 ; i++){
			if(operation[i]=='~'&&operation[i+1]=='D'&&operation[i+2]=='~'){
				int start=0;

				for(start=i-1 ; start>0 ; start--)
					if(Is_operation(operation[start-1])&&start!=1)
						break;
				int end=0;

				for( end=i+3 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(start,end,operation,operation[i]);

                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
			for(int i=0 ; i<operation.length() ; i++){
			if(operation[i]=='-' && i!=0){
            int end=0;

				for( end=i+1 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
						if(Is_operation(operation[i-1]))
                            i--;
				int start=0;

				for(start=i-1 ; start>0 ; start--)
					if(Is_operation(operation[start-1])&&start!=1)
						break;

				calcAndRep(start,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}

			for(int i=operation.length()-1 ; i>=0 ; i--){
			if(operation[i]=='+'){
				int start=0;

				for(start=i-1 ; start>0 ; start--)
					if(Is_operation(operation[start-1])&&start!=1)
						break;

				int end=0;

				for( end=i+1 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(start,end,operation,operation[i]);
                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){
			if(operation[i]=='~'&&operation[i+1]=='P'&&operation[i+2]=='~'){
				int start=0;

				for(start=i-1 ; start>0 ; start--)
					if(Is_operation(operation[start-1])&&start!=1)
						break;
				int end=0;

				for( end=i+3 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(start,end,operation,operation[i]);

               // cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}
		for(int i=operation.length()-1 ; i>=0 ; i--){
			if(operation[i]=='~'&&operation[i+1]=='S'&&operation[i+2]=='~'){
				int start=0;

				for(start=i-1 ; start>0 ; start--)
					if(Is_operation(operation[start-1])&&start!=1)
						break;
				int end=0;

				for( end=i+3 ; end<operation.length()-1; end++)
					if(Is_operation(operation[end+1])  )
						break;
				calcAndRep(start,end,operation,operation[i]);

                //cout<<operation<<endl;
				Operation_solver(operation) ;
				break;
			}
		}



	}

string mul_ope_solver(string &ope)
	{
	    Operation_solver(ope);

	    if(memoryCheck(ope)==-1 )
         {
              string temp2,tName;
             tName=genRandom();
              temp2=tName+"=["+ope+"]";
              input_checker(temp2);
              //cout<<"operatin is"<<ope<<endl;
              return tName;
           }
        return ope;

	}

/*End Here*/
int main(int argv, char* argc[])
{
   /* matrix a; a.initialling("a","1 2 3");
    matrix b; b.initialling("b","1 2 3");
    b.Msind(a);
    b.print();*/

    string ts;
    string inp="A=[2 1 1]";
    input_checker(inp);
    inp="D=[20 25;30 35]";
    input_checker(inp);
    inp="C=[1 2;3 4]";
    input_checker(inp);
  ts="((C*D .+ 4)./2.1 + sqrt(D))./C.^2";

//ts="(4+1)*2^3+1";
//ts="((1+1)-2^2^2)+1";

   cout<<endl<<"----******--"<<endl;
   string mainRes=mul_ope_solver(ts);
    cout<<mainRes<<endl;
    cout<<"--******--"<<endl;

   /* ios_base::sync_with_stdio(false);
    cin.tie(0);
	if (argv>1)
	{
		ifstream infile(argc[1]);
		string sFile, temp;
		while (getline(infile, temp))
		{
			if(temp.find("\r")!=-1)
				temp.replace(temp.find("\r"),2,"");
			sFile += temp;
			if (sFile.find("]") != -1 || sFile.find("];") != -1 || sFile.find("+") != -1 || sFile.find("*") != -1 || sFile.find("/") != -1 || sFile.find("'") != -1 || sFile.find("./") != -1 || (sFile.find("-") != -1 && sFile.length() <= 10))
			{
				input_checker(sFile);
				sFile = "";
			}
		}
		infile.close();
	}else
    while(1){
            exit1=0;
            cout<<">> ";
     string ins;
     getline (cin,ins);
     if(ins.find("\r")!=-1)
				ins.replace(ins.find("\r"),2,"");
     input_checker(ins);
     if(exit1==1)
        break;
*/

cout<<endl<<"memory----------------"<<endl;
for(int z=0;z<memoryPointer;z++)
   {
    cout<<memory[z].getName()<<" "<<memory[z].getRows()<<" "<<memory[z].getColumns()<<endl;
    memory[z].print();
    cout<<"#############"<<endl;
   }
cout<<endl<<"--------------------------"<<endl;

//}

    return 0;
}
