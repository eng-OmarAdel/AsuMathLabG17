#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

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

    void initialising(string mName,string mString) // give me string wana azzabat isa
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
		initialising(mName,mString);

    }

    void initialising(int rows1, int columns1)//give rows and columns wana azzabat isa
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

        this->initialising(x.columns, x.rows);
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
        this->initialising((x.rows)-1, (x.columns)-1);
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
        this->initialising(x.rows, x.columns);
	z.initialising(x.rows, x.columns);
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
		this->initialising(x.rows, x.columns);
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
        this->initialising(x.rows, x.columns);
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
        this->initialising(x.rows, x.columns);
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
        this->initialising(x.rows, y.columns);//use constructor instead
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
        this->initialising(x.rows, inverseDenom.getColumns());
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
        p[usedSlots].initialising(name,mString);
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
        // cout<<"mwmory"<<memoryPointer;
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
struct sizeValue { int rows,columns; };   //structure has the number of rows&columns for any string matrix
//===========================================================================

vector <string> separatedString;

void separate(string inputString)
{
    //vector<string> separatedString;
    int flag=0;
    int beginPostion=0;
    for(int i=0;i<inputString.length();i++)
    {
        if(inputString[i]==';'&&flag==0)
        {
            //cout<<inputString.substr(beginPostion,i-beginPostion)<<endl;
            separatedString.push_back(inputString.substr(beginPostion,i-beginPostion));
            beginPostion=i+1;
        }
        else if(i==inputString.length()-1)
        {
            //cout<<inputString.substr(beginPostion,inputString.length()-beginPostion);
            separatedString.push_back(inputString.substr(beginPostion,inputString.length()-beginPostion));
        }
        else if(inputString[i]=='[')
            flag++;
        else if(inputString[i]==']')
            flag--;
    }
}

//==================================================================
sizeValue get_sizeValue (string s)  //s is a string without spaces or ;   ex:   123  or  10*sin(A)
{
	sizeValue n;
    int flag=0;
    int index;
    for(int i=0;i<s.length();i++)
    {
        if((s[i]>64&&s[i]<91)||(s[i]>96&&s[i]<123))
        {
			flag=1;
            if(s[i]=='s'&&s[i+1]=='i'&&s[i+2]=='n')    //if sin
			{
                index=memoryCheck(s.substr(i+4,s.find(')')-i-4));
				if(index==-1)
				{
					n.rows=1;  n.columns=1;
					return n;
				}
				n.columns=memory.p[index].getColumns();
				n.rows=memory.p[index].getRows();
				return n;
			}
            else if(s[i]=='c'&&s[i+1]=='o'&&s[i+2]=='s')  // if cos
            {   
		        index=memoryCheck(s.substr(i+4,s.find(')')-i-4));
				if(index==-1)
				{
					n.rows=1;  n.columns=1;
					return n;
				}
				n.columns=memory.p[index].getColumns();
				n.rows=memory.p[index].getRows();
				return n;
			}
            else if(s[i]=='t'&&s[i+1]=='a'&&s[i+2]=='n')   //  if tan
			{
                index=memoryCheck(s.substr(i+4,s.find(')')-i-4));
				if(index==-1)
                {
					n.rows=1;  n.columns=1;
					return n;
				}
				n.columns=memory.p[index].getColumns();
				n.rows=memory.p[index].getRows();
				return n;
			}
            else if(s[i]=='r'&&s[i+1]=='a'&&s[i+2]=='n'&&s[i+3]=='d')   // if rand(55,100)
			{
				n.columns=stoi(s.substr(s.find(',')+1,s.find(')')-s.find(',')-1));
				n.rows=stoi(s.substr(s.find('(')+1,s.find(',')-s.find(',')-1));
                return n;  //return numbers of columns 100 & rows 55
			}
            else if(s[i]=='e'&&s[i+1]=='y'&&s[i+2]=='e')     //if eye(55,100)
            {
				n.columns=stoi(s.substr(s.find(',')+1,s.find(')')-s.find(',')-1));
				n.rows=stoi(s.substr(s.find('(')+1,s.find(',')-s.find(',')-1));
                return n;  //return numbers of columns 100 & rows 55
			}
            else if(s[i]=='z'&&s[i+1]=='e'&&s[i+2]=='r'&&s[i+3]=='o'&&s[i+4]=='s')   // if zeros(55,100)
            {
				n.columns=stoi(s.substr(s.find(',')+1,s.find(')')-s.find(',')-1));
				n.rows=stoi(s.substr(s.find('(')+1,s.find(',')-s.find(',')-1));
                return n;  //return numbers of columns 100 & rows 55
			}
            else if(s[i]=='o'&&s[i+1]=='n'&&s[i+2]=='e'&&s[i+3]=='s')  // if ones(55,100)
            {
				n.columns=stoi(s.substr(s.find(',')+1,s.find(')')-s.find(',')-1));
				n.rows=stoi(s.substr(s.find('(')+1,s.find(',')-s.find(',')-1));
                return n;  //return numbers of columns 100 & rows 55
			}
            else
			{
				index=memoryCheck(s.substr(i,1));
				n.columns=memory.p[index].getColumns();
				n.rows=memory.p[index].getRows();
				return n;
			}
        }
	}
	if (flag==0)
	{
		n.rows=1;  n.columns=1;
		return n;
	}
}
//======================================================================================================

sizeValue sizing (string matrix)
{
    sizeValue n; n.rows=0;  n.columns=0;
	if(matrix[0]=='[')
		matrix=matrix.substr(1,matrix.length()-2);
    stringstream sMatrix(matrix);
    string token;
    getline(sMatrix,token,';');
    stringstream sn(token);
    string element;
    while(sn>>element)
    {
		n.columns+=get_sizeValue(element).columns;
    }
	//-----------------------------------
	stringstream ssMatrix(matrix);
	while(getline(ssMatrix,token,';'))
	{
		stringstream sn(token);
        string element;
		sn>>element;
		n.rows+=get_sizeValue(element).rows;
	}
	return n;
}
//=============================Concatenation==================================================
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
		for(int o=i ; o<(s.length()) ; o++)
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
	int temp=Vstack[0].rows;               //these 3 lines are because the function reverse rows with columns
	Vstack[0].rows=Vstack[0].columns;      // so i exchang them
	Vstack[0].columns=temp;
    return Vstack[0];
}
//===============================================================================================
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
    if(finStack.size() == 1) separatedString.clear(); return finStack[0];
    sizeValue sum = compare(finStack[0], finStack[1]);
    for(int i=2 ; i<finStack.size() ; i++)
    {
        sum = compare(sum, finStack[i]);
    }
	separatedString.clear();
    return sum;
}
//=================================================================================
#define endl '\n'

int main(int argv,char* argc[])
{
    ios_base::sync_with_stdio(false);
    cin.tie(0);

	if (1)//(argv>1)
	{
		ifstream infile("D://example.txt");
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
for(int z=0;z<memoryPointer;z++)
   {
    cout<<memory[z].getName()<<" "<<memory[z].getRows()<<" "<<memory[z].getColumns()<<endl;
    memory[z].print();
    cout<<"#############"<<endl;
   }
cout<<endl<<"--------------------------"<<endl;
   */ }

    return 0;
}
