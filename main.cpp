#include <iostream>
#include <string>
#include<vector>
#include <sstream>
#include <fstream>

using namespace std;

double StringToDouble(const string &text)
{
	stringstream ss(text);
	double result;
	return ss >> result ? result : 0;
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

    void initialling(string matrix_name,string matrix) // give me string wana azzabat isa
    {
        mName=matrix_name;
        string matrix2 = matrix.substr(1,matrix.length()-2); // to remove curly brackets
        stringstream ss(matrix2);
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
        stringstream ss3(matrix2);
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
    void update(string matrix_name,string matrix)
    {

        setRows(0);
        setColumns(0);
        for(int i=0;i<rows;i++)
		{
			delete[] element[i];
		}
		delete[] element;
		initialling(matrix_name,matrix);

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

    void swapRows(int ft,int sc)
        {
            double *temp=this->element[ft];
            this->element[ft]=this->element[sc];
            this->element[sc]=temp;

        }
     double getDeterminant()
    {

        if(this->columns==this->rows)
        {
            int n= this->rows;

            for(int col = 0; col < n; ++col) {


      bool found = false;
      for(int row = col; row < n; ++row) {
         if(this->element[row][col]) {
                    if ( row != col )
            {
               this->swapRows( row, col);

            }
            else
            {

            }
            found = true;
            break;
         }
      }
      if(!found) {

         return 0;
      }

      for(int row = col + 1; row < n; ++row) {
         while(true) {
            int del = this->element[row][col] / this->element[col][col];

            for (int j = col; j < n; ++j) {
               this->element[row][j] -= del * this->element[col][j];
            }
            if (this->element[row][col] == 0)
            {
               break;
            }
            else
            {
               this->swapRows(col,row);

            }
         }
      }

        }
        double res = 1;

   for(int i = 0; i < n; ++i) {
      res *= this->element[i][i];
   }
   return res;

    }
    else
            return 0;
    }
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
	void print()
	{
		if(errorHandler=="Error There's a zero element in the matrix" || errorHandler=="Error The determinant of this matrix is eual to zero")
			cout<<errorHandler;
		else
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
		cout<<endl;
	}

	~matrix()
	{
		for(int i=0;i<rows;i++)
		{
			delete[] element[i];
		}
		delete[] element;

	}

};
class ram
{
    int mSize;
    int used;
public:
    matrix *p;

    ram()
    {
        p=new matrix[1];
        used=0;
    }
    ram(int s)
    {
        mSize=s;
        if(mSize<=0)
            mSize=1;
        p=new matrix[mSize];
        used=0;
    }
    int size()
    {
        return mSize;
    }
    int getUsed()
    {
        return used;
    }

    matrix get(int i)
    {
        return p[i];
    }
    matrix operator [] (int i)
    {
        return get(i);
    }
    void creat(string name)
    {
        p[used].setName(name);
        refresh();
    }
    void creat(string name,string matrix)
    {
        p[used].initialling(name,matrix);
        refresh();
    }
    void update(int num,string name,string matrix)
    {
        p[num].update(name,matrix);
    }
    void refresh()
    {
        used++;
        if(used==mSize)
        {

            mSize+=10;
            matrix *t=new matrix [mSize];
            for(int i=0;i<used;i++)
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
ram memory (10);
int memory_arrow=0;
int exit1=0;

void removeSpaces(string &str)
{
    for(int i=0; i<str.length(); i++)
     {
         if(str[i] == ' ') {str.erase(i,1);i--;}

     }
}


int memory_check(string matrix_name)
{
    for(int i=0;i<memory_arrow;i++)
    {
        if(memory[i].getName()==matrix_name)
            return i;
    }
    return -1;

}
void memorize_matrix(int index, string matrix,string matrix_name)
{
    if(index==-1)
    {
        // cout<<"mwmory"<<memory_arrow;
        memory.creat(matrix_name,matrix);
        memory.p[memory_arrow].print();
        memory_arrow++;
    }

    else
    {
        memory.update(index,matrix_name,matrix);
        memory.p[index].print();
    }
}
void cut(string &variable1,string &variable2,int &index1,int &index2,char op,string operation)
{
    variable1= operation.substr(0,operation.find(op)) ;
    variable2= operation.substr(operation.find(op)+1,(operation.length()-operation.find(op))-1);

    index1=memory_check(variable1);
    index2=memory_check(variable2);
}
void cut(string &variable1, int &index1, char op, string operation)
{
	variable1 = operation.substr(operation.find(op) + 1, (operation.length() - operation.find(op)) - 1);
	index1 = memory_check(variable1);
}
void input_checker(string input) // assignment or operation
{

	int flag=0;
	int asg=0;

    string matrix_name;
	matrix_name=input.substr(0,input.find('='));
	removeSpaces(matrix_name);


    int index = memory_check(matrix_name);

    if((input.find('[')!=-1)&&(input.find(']')!=-1))
    {

        flag=0;
        string matrix;
	    matrix=input.substr(input.find('['),(input.find(']')-input.find('[')+1));
        memorize_matrix(index,matrix,matrix_name);


    }

    else if((input.find("+")!=-1)||
    (input.find("-")!=-1)||(input.find("/")!=-1)||(input.find("*")!=-1)||(input.find("'")!=-1))
    {
		    string operation;
		    operation=input.substr(input.find('=')+1,(input.length()-input.find('=')+1));
            removeSpaces(operation);
            string variable1,variable2;
            int index1,index2;

			flag=1;

				//memory.push_back(temp);
			if(operation.find("+")!=-1)
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
                memory.creat(matrix_name);
                memory.p[memory_arrow].add(memory.p[index1],(memory.p[index2])) ;
                memory.p[memory_arrow].print();
                memory_arrow++;
                }
            }

			else if(input.find("-")!=-1)
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
                memory.creat(matrix_name);
                memory.p[memory_arrow].sub(memory.p[index1],(memory.p[index2])) ;
                memory.p[memory_arrow].print();
                memory_arrow++;
                }


			}
			else if(input.find("*")!=-1)
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
                memory.creat(matrix_name);
                memory.p[memory_arrow].mult(memory.p[index1],(memory.p[index2])) ;
                memory.p[memory_arrow].print();
                memory_arrow++;
                }
			}
			else if (input.find("./") != -1)
			{
				cut(variable1, index1, './', operation);
				if (index != -1)
				{
					memory.p[index].inversePerElement(memory.p[index1]);
					memory.p[index].print();
				}

				else
				{
					memory.creat(matrix_name);
					memory.p[memory_arrow].inversePerElement(memory.p[index1]);
					memory.p[memory_arrow].print();
					memory_arrow++;
				}
			}
			else if(input.find("/")!=-1)
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
                memory.creat(matrix_name);
                memory.p[memory_arrow].div(memory.p[index1],(memory.p[index2])) ;
                memory.p[memory_arrow].print();
                memory_arrow++;
                }
            }
			else if(operation.find("'")!=-1)
			{
			    string var;
			    int ind;
			    var= operation.substr(0,operation.find("'")) ;
                ind=memory_check(var);


                if(index!=-1)
                {
                memory.p[index].getTranspose(memory.p[ind]) ;
                memory.p[index].print();
                }

                else
                {
                memory.creat(matrix_name);
                memory.p[memory_arrow].getTranspose(memory.p[ind]) ;
                memory.p[memory_arrow].print();
                memory_arrow++;
                }
			}


}

		else
        {
            removeSpaces(input);
             if((memory_check(input)!=-1) &&
                 memory_check(input)<memory_arrow)
             {
                 memory.p[memory_check(input)].print();
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

int main(int argv,char* argc[])
{
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

/*
cout<<endl<<"memory----------------"<<endl;
for(int z=0;z<memory_arrow;z++)
   {
    cout<<memory[z].getName()<<" "<<memory[z].getRows()<<" "<<memory[z].getColumns()<<endl;
    memory[z].print();
    cout<<"#############"<<endl;

   }
cout<<endl<<"--------------------------"<<endl;

   */ }
    return 0;
}
