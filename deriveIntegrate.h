class Derivative
{
public:
  virtual Field f(Field argument)=0;

  Field twoPoint(Field argument,Field distance)  //O(distance)
  {  return (f(argument+distance)-f(argument))/distance;  }

  Field threePoint(Field argument,Field distance) //O(distance^2)
  {  return (f(argument+distance)-f(argument-distance))/(2*distance);  }
};


class Integrate
{
public:
  virtual Field f(Field argument)=0;

  Field trapezium(Field lowerBound,Field upperBound,Count n)
  {
    Field distance=(upperBound-lowerBound)/Field(n),
          result=(f(upperBound)+f(lowerBound))*distance/2.;

    --n;for(int i=1;i<n;++i)
      result+=f(Field(i)*distance+lowerBound)*distance;
   
    return result;
  }

  Field simpson(Field lowerBound,Field upperBound,Count n)
  {
    n-=n%2;
    Field distance=(upperBound-lowerBound)/Field(n),
          norm    =2.*distance/3.,
          result  =(f(upperBound)+f(lowerBound))*norm/2.;

    --n;for(int i=1;i<n;++i)
      result+=(1<<(i%2))*f(Field(i)*distance+lowerBound)*norm;
   
    return result;
  }

};
