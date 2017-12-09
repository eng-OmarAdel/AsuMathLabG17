#include <iostream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

struct sizeValue
{
    int rows;
    int columns;
};

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
    //return separatedString;
}

sizeValue compare(sizeValue m1, sizeValue m2)     //get the total size of 2 concatenated matrices
{
    sizeValue m;
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
sizeValue sizing(string matrix)  //take any string but without brackets and retur a structure which has the number of rows and columns
{                               // for example    22 45 sin(90)+1;1 cos(0) 5    size will be 2x3
	sizeValue n; n.rows=0;  n.columns=0;
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
    if(finStack.size() == 1) return finStack[0];
    sizeValue sum = compare(finStack[0], finStack[1]);
    for(int i=2 ; i<finStack.size() ; i++)
    {
        sum = compare(sum, finStack[i]);
    }
    return sum;
}

int main()
{
    string s = "[1 2 3 4 5],[[1 2;1 2],[[1 2;1 2],[1;2]]]";//"[1 2 3;1 2 3],[1;2]";
    separate(s);
    sizeValue x = calcSize(separatedString);
    cout << x.rows << " * " <<x.columns;
    return 0;
}
