#include <cstdlib>
#include <string>
#include <iostream>


class LinearAlgebraFile
{
protected:
  static void fromString(const string& s, double& out)
  {
    out=std::atof(s.c_str());
  }
  static void fromString(const string& s, long double& out)
  {
    out=std::atof(s.c_str());
  }
  static void fromString(const string& s, float& out)
  {
    out=std::atof(s.c_str());
  }
  static void fromString(const string& s, int& out)
  {
    out=std::atoi(s.c_str());
  }
  static void  fromString(const string& s, long int& out)
  {
    out=std::atol(s.c_str());
  }

public:
  template<class ComponentSType>  static Array < Array < ComponentSType > > readMatrixFromFile(const char* fileName,const Count skipN=0)
  {
    Array < Array < ComponentSType > > matrix;
    Count rowSN=0,colSN=0;
    ifstream read(fileName);
    string lineString;

    for(Count skippedN=0;skippedN<skipN;++skippedN)
      getline(read,lineString);

    while(!read.eof()&&!read.fail()) //get matrix size
    {
      getline(read,lineString);
      if (lineString.c_str()!=NULL&&lineString!=""&&lineString.c_str()[0]!='#')
      {
        if(rowSN==0)
        {
          while (lineString!="")
          {
            int position = lineString.find_first_of(" ");
            if(position==0)
              lineString.erase(0,1); // remove spaces at the beginning
            else
            {
              if (position==(int)string::npos) // check if a space was found at all
                position = lineString.length();
              lineString.erase(0,position); // erase the processed part from the string
              ++colSN; // update the column counter to the current position
            }
          }
        }
        ++rowSN;
      }
    }
    read.close();
    if(rowSN==0||colSN==0)
      return matrix;
    matrix=LinearAlgebra::createMatrix< ComponentSType >(rowSN,colSN);

    Count i=0; //row
    read.open(fileName);
    for(Count skippedN=0;skippedN<skipN;++skippedN)
      getline(read,lineString);
    while(!read.eof()) //get matrix size
    {
      getline(read,lineString);
      if (lineString.c_str()!=NULL&&lineString!=""&&lineString.c_str()[0]!='#')
      {
        Count j=0;//column

        while (lineString!=""&&j<colSN)
        {
          int position = lineString.find_first_of(" ");
          if(position==0)
            lineString.erase(0,1); // remove spaces at the beginning
          else
          {
            if (position==(int)string::npos) // check if a space was found at all
              position = lineString.length();
            fromString(lineString.substr(0,position),matrix[i][j]);
            lineString.erase(0,position); // erase the processed part from the string
            ++j; // update the column counter to the current position
          }
        }
        for(;j<colSN;j++) //initialise missing data
          matrix[i][j]=0;

        ++i;
      }
    }
    read.close();

    return matrix;
  }


  template<class ComponentSType>static void writeMatrix(const char* fileName, const Matrix(ComponentSType) matrix, bool append=false)
  {
    if(matrix.pData==NULL)
      return;
    std::ofstream write(fileName, ( append ? std::ios::out|std::ios::app : std::ios::out ));

    for (Count i=0;i<matrix.pData->componentSCount;++i)
    {
      for (Count j=0;j<matrix[i].pData->componentSCount;++j)
        write<<matrix[i][j]<<" ";
      write<<"\n";
    }

    write.close();
  }


  template<class ComponentSType>static void writeVector(const char* fileName, const Vector<ComponentSType> vector, bool writeIndex=false, bool append=false)
  {
    if(vector.pData==NULL)
      return;
    std::ofstream write(fileName, ( append ? std::ios::out|std::ios::app : std::ios::out ));

    for (Count i=0;i<vector.pData->componentSCount;++i)
    {
      if(writeIndex)
        write<<i<<" ";
      write<<vector[i]<<"\n";
    }

    write.close();
  }


  template<class ScalarUser=void,class K,class Scalar=first_nonVoid_t<ScalarUser,fundamental_t<K>>> static void writeTensor(const char* fileName, const K tensor, bool append=false) // writes tensor componentwise w leading indices.
  //actually this function is able to replace all write[LinearAlgebraElemente] functions: Mat(Mat(K)) m; writeTensor_wIndices<Mat(K)>("...",m); writes the first two dimensions w leading indices and the remaining Mat(K) as record
  {
    std::ofstream write(fileName, ( append ? std::ios::out|std::ios::app : std::ios::out ));

    if constexpr(std::is_same_v<Scalar,K>)
      write<<"0 "<<tensor<<"\n";
    else
      ArrayTools::map(tensor,[&](const Scalar &s,const Array<Count> &index) {  write<<index<<"  "<<s<<"\n";  });

    write.close();
  }


  template<class ComponentSType>static void write(const char* fileName, const Vector<ComponentSType> vector, bool writeIndex=false, bool append=false)
  {  writeVector( fileName, vector, writeIndex );  }
  template<class ComponentSType>static void write(const char* fileName, const Matrix(ComponentSType) matrix, bool writeIndex=false, bool append=false)
  {  writeMatrix( fileName, matrix, writeIndex );  }
};
