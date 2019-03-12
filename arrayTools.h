#pragma once
#include<tuple>

template<class Tt,class Ts> struct copy_referenceType //similar to forward, can also handle non-reference types
{
  template<class=typename std::enable_if< std::is_reference<Ts>::value>::type > static typename remove_reference<Tt>::type &&  helper(typename remove_reference<Ts>::type &&);
  template<class=typename std::enable_if< std::is_reference<Ts>::value>::type > static typename remove_reference<Tt>::type  &  helper(typename remove_reference<Ts>::type & );
  template<class=typename std::enable_if<!std::is_reference<Ts>::value>::type > static typename remove_reference<Tt>::type     helper(                          Ts          );
  //typename remove_reference<Tt>::type     helper(typename remove_reference<Ts>::type   );

  typedef decltype(helper(std::declval<Ts>())) type;
};


class ArrayTools
{
protected:

  template<class ComponentSType> struct SubArraySArrayLoop
  {
    Array < ComponentSType > result;
    const Count& dimensionIdRemoveSIndex;

    template<class ScalarType>void process(const Array< ScalarType > & array __attribute__ ((unused)),double result __attribute__ ((unused)),const Count dimensionIdRemaining __attribute__ ((unused)))
    {  }//template recursion anchor

    template<class SubComponentSType>void process(const Array< Array< SubComponentSType > >& array,Array< SubComponentSType >& result,const Count dimensionIdRemaining)
    {
      result.resize(array.pData->componentSCount);
      if(dimensionIdRemaining==1)
      {
        for(Count i=0;i<array.pData->componentSCount;++i)
          result[i]=array[i][dimensionIdRemoveSIndex];
      }
      else
        for(Count i=0;i<array.pData->componentSCount;++i)
          process(array[i],result[i],dimensionIdRemaining-1);
    }


    SubArraySArrayLoop(const Array< Array<ComponentSType> >& array,const int& dimensionIdRemove,const Count& i)
      :dimensionIdRemoveSIndex(i)
    {
      process (array,result,dimensionIdRemove);
    }
  };

public:
  template<class SubSettable,class Scalar=fundamental_t<SubSettable>> inline static constexpr Count rank=Array<bool>::rank<SubSettable,Scalar>;

  template<class Invocable> struct get_args
  {
    template<class R,class... Args>               static std::tuple<Args...>                    helper( R(*)(Args...) );                   // ptr to func
    template<class R,class Closure,class... Args> static std::tuple<Args...>                    helper( R(Closure::*)(Args...) );          // ptr to member func
    template<class R,class Closure,class... Args> static std::tuple<Args...>                    helper( R(Closure::*)(Args...) const);     // ptr to member func const
    template<class Closure>                       static decltype(helper(&Closure::operator())) helper( Closure );                         // lambda
    //template<class Closure,class=typename enable_if<is_class<Closure>::value>::type> static auto                helperA( Closure ) -> decltype(helperA(&Closure::operator()));  // lambda

    typedef decltype(helper(std::declval<Invocable>())) type;
  };
  template<class Invocable> using get_args_t = typename get_args<Invocable>::type;


  template< class Array, class Invocable, class... Parameter > struct Map
  {
  public:
    typedef remove_cr_t < typename std::tuple_element<0,typename get_args<Invocable>::type>::type> Scalar;
    typedef remove_cr_t < typename std::tuple_element<1,typename get_args<Invocable>::type>::type> Index;

  private:
    Invocable                    &&invocable;
    std::tuple<Parameter &&...>    parameter;
    static const bool   subsettableIndex = is_subsettable<Index>::value;

    //template<std::enable_if_t<!subsettableIndex>* =nullptr,class ArrayL=Array> struct IndexContainer             {  };
    template<class ArrayL=Array, bool=subsettableIndex> struct IndexContainer {};
    template<class ArrayL>                              struct IndexContainer<ArrayL,true>
    {
      Index v;
      IndexContainer()
      {
        if constexpr ( inheritsArray<Index>::value )
          v.n(DimensionCount<ArrayL,Scalar>::value);
      }
    };

    IndexContainer<> indexContainer;

    template< int d=0, std::size_t... PI, class ArrayL > // PI:<=>ParameterIndices: needed to transform tupel parameter into a varic list via "std::get<PI>(param)...
       inline void map (ArrayL &&array, std::index_sequence< PI... >* = NULL )
    {
      auto indexLambda = [this]()-> decltype(auto) // is either Index{0} or (Index[]&){indexContainer.v[d]}  => we loop over a plain scalar or use the previously
                         {  if constexpr ( subsettableIndex ) return indexContainer.v[d]=0; else return std::decay_t<Index>{0};  };

      for( decltype(auto) i=indexLambda() ;i<array.size(); ++i ) //declt(auto) : replace auto w initialization list & uses rules of decltype => reference detection.
        if constexpr( ! std::is_same
                      <
                        std::remove_cv_t< std::remove_reference_t<decltype(array[0])> >,
                        Scalar
                      >::value )
          map< d+1, PI... >(array[i]);
        else
          if constexpr( subsettableIndex )
            invocable( array[i],indexContainer.v, get<PI> (parameter)... );
          else
            invocable( array[i],i               , get<PI> (parameter)... );
      //return std::forward<Array>(array);
    }


  public:

    Map(Array &&array, Invocable &&invocable, Parameter &&... parameter)
      : invocable(std::forward<Invocable>(invocable)), parameter(std::forward<Parameter>(parameter)...)
    {
      map(std::forward<Array>(array), static_cast< std::make_index_sequence<std::tuple_size_v<std::tuple<Parameter &&...>>>* >(NULL) );
    }

  };


  template<class Array, class Invocable, class... Parameter>
     inline static decltype(auto)
       map (Array &&array,Invocable &&invocable,Parameter &&... parameter)
  // The invocable should have the following signature void (ComponentSType&(&), IntegralType OR Subsettable Integraltype, Parameter...)
  //   whereas Array is resolved into Elements recursively until its ComponentSType equals ComponentSType from Invocable up to cv qualification
  //   Simple one-parameter-Invocables can be mapped directly using the suceeding wrapper.
  // map could be realised with a different-efficient signature: Array map(Array, Invocable&&, Parameters...): for_array(array,i) array[i]=invocable(array[i],i,parameter...);
  //   return value optimisation & parameter replacement theoretically enable copy elusion, st. array=map(array, . , . ); does not result in an array copy.
  //   this ansatz allows might in some cases allow for more copy elusion: pointer-alias-problems.
  {
    if constexpr(!inheritsArray<Array>::value)
    {
       typename Map< Array, Invocable, Parameter... >::Index i {0};
       invocable(array,i,parameter...);
    }
    else
      Map< Array, Invocable, Parameter... > map(std::forward<Array>(array),std::forward<Invocable>(invocable),std::forward<Parameter>(parameter)...);

    return std::forward<Array>(array);
  }

  template<class Array, class Invocable, std::enable_if_t< std::tuple_size<get_args_t<Invocable>>::value==1 >* =nullptr >
    inline static decltype(auto)
      map (Array &&array,Invocable &&func)
  //Invocables w only one parameter have a different calling signature than Lambda(F,Integral,Parameters...), that is used in the more general "map" above.a
  //  For these simple Invocables, this wrapper of "map" can be invoked directly:
  //   if func has a return value: x=func(x) ; x is an element of array matching the type of f's argument
  //   if func returns void, only func(x) is invoced.
  //     func should than a have a reference typed argument or side effects in order to have a non-trivial map expression
  {
    typedef typename std::tuple_element<0,get_args_t<Invocable>>::type FuncSArg;
    return map( std::forward<Array>(array), [&func](FuncSArg x,Count) // FuncSArg already may contain && or &
                                            {
                                              if constexpr ( std::is_same<void,decltype(func(std::declval<FuncSArg>()))>::value )
                                                func(std::forward<FuncSArg>(x));
                                              else
                                                x=func(std::forward<FuncSArg>(x));
                                            } );
  }


  template<class ComponentSType> static Array < ComponentSType > subArray(const Array < Array < ComponentSType > > array,const int  dimensionIdRemove,const Count i)
  {
    if(dimensionIdRemove==0)
      return array[i];
    SubArraySArrayLoop<ComponentSType> arrayLoop(array,dimensionIdRemove,i);

    return arrayLoop.result;
  }


  template<class ComponentSType> static Array < ComponentSType > interval(const Array < ComponentSType > array,Count min,Count max)
  {
    if(min<0)  min=0;  if(max>=array.n())  max=array.n()-1;
    Array < ComponentSType > result(max-min+1,NULL);

    for(Count i=0;i<result.n();++i)
      result.atFast(i)=array[i+min];

    return result;
  }


  template<bool maximum,bool returnIndex=false, class ScalarUser=void,class SubSettable,class Scalar=remove_cr_t<first_nonVoid_t<ScalarUser,decltype(std::declval<SubSettable>()[0])>>> static decltype(auto) extremum(SubSettable &&array)
  {
    typedef std::array<Count,rank<SubSettable,Scalar>> I;
    I r{0};
    Scalar       extrV=std::as_const(array).template at<Scalar>(r);
    map(array, [&r,&extrV](const Scalar &s, const I &i) // later, if opr[] passes rvalues (or: map uses a global at(Array,Index) is used that passes rvalues) Scalar &s might need to become an rvalue if array is one.
    {
      if( (maximum&&s>extrV) || (!maximum&&extrV>s) )
        { extrV=s; for(Count d=0; d< rank<SubSettable,Scalar>; ++d) r[d]=i[d]; }
    } );

    if constexpr (returnIndex)
      if constexpr (rank<SubSettable,Scalar> ==1)
        return r[0]; else return Array<Count>(r);
    else
      return forwardMember<SubSettable>(array.template at<Scalar>(r));
  }

  template<class ComponentSType> static Count maximalComponentId(const Array< ComponentSType > array) //simple but fast.
  {
    Count result=0;
    for(Count i=1;i<array.pData->componentSCount;++i)
      if(array[i]>array[result])
    result=i;

    return result;
  }

  template<class ComponentSType> static Count minimalComponentId(const Array< ComponentSType > array)
  {
    Count result=0;
    for(Count i=1;i<array.pData->componentSCount;++i)
      if(array[i]<array[result])
    result=i;

    return result;
  }

  template<class ScalarUser=void,class K> static decltype(auto) max(K &&array)
  {  if constexpr (!is_subsettable_v<K>) return array; else return extremum<true ,false,ScalarUser>(std::forward<K>(array));  }
  template<class ScalarUser=void,class K> static decltype(auto) min(K &&array)
  {  if constexpr (!is_subsettable_v<K>) return array; else return extremum<false,false,ScalarUser>(std::forward<K>(array));  }


  template<class ComponentSType> static Array< Array <Count> > histogramClassesSIndexArray(const Array< ComponentSType >& array, Count classN, ComponentSType minimum, ComponentSType maximum)
  {
    Array< Array <Count> > classesSIndexArray(classN,NULL);
    Array <Count> classesSSize(classN,0,NULL);

    struct
    {
      Count get(Count classN,ComponentSType& minimum, ComponentSType& maximum,ComponentSType value)
      {	return Count(double(value-minimum)/double(maximum-minimum)*(double)(classN)); }
    } classId;

    for(Count i=0;i<array.pData->componentSCount;++i)
    {
      Count c=classId.get(classN,minimum,maximum,array[i]);
      if(c>=0&&c<classN&&array[i]>=minimum&&array[i]<=maximum)
	      classesSSize[c]++;
    }
    for(Count c=0;c<classN;++c)
      classesSIndexArray[c].resize(classesSSize[c]);
    {
      Array <Count>& classesSFillId=classesSSize; //reuse memory
        classesSFillId=0;

      for(Count i=0;i<array.pData->componentSCount;++i)
      {
        Count c=classId.get(classN,minimum,maximum,array[i]);
        if(c>=0&&c<classN&&array[i]>=minimum&&array[i]<=maximum)
        {
      	  classesSIndexArray[c][classesSFillId[c]]=i;
	        ++classesSFillId[c];
        }
      }
    }

    return classesSIndexArray;
  }


  template<class ComponentSType> static void quickSortIndexArray(const Array< ComponentSType >& array,Array< Count >& indexArray,Count componentIdMinimal,Count componentIdMaximal)
  {
    std::cout<<"ERROR: template<class ComponentSType> static void quickSortIndexArray(const Array< ComponentSType >& array,Array< Count >& indexArray,Count componentIdMinimal,Count componentIdMaximal) contains a bug in the recursive calls.  -- limits? See. Cf. wikipedia -> quicksort\n";
    ComponentSType pivot   (array[indexArray[componentIdMaximal]]);
    Count          leftId =componentIdMinimal;
    Count          rightId=componentIdMaximal-1;

    std::cout<<"A_1\n";
    do
    {
      std::cout<<"A_2 "<<leftId<<" "<<rightId<<"\n";
      for( ;array[indexArray[leftId]]<=pivot && leftId<componentIdMaximal; ++leftId);
      for( ;array[indexArray[rightId]]>=pivot && rightId>componentIdMinimal; --rightId);

      if(leftId<rightId) //swap ids
      {
        Count tmpId=indexArray[leftId];
        indexArray[leftId]=indexArray[rightId];
        indexArray[rightId]=tmpId;
      }
    }
    while(leftId<rightId);

    if( array[indexArray[leftId]] > pivot )
    {
      Count tmpId=indexArray[leftId];
      indexArray[leftId]=indexArray[componentIdMaximal];
      indexArray[componentIdMaximal]=tmpId;
    }// => all elements on the left are smaller than pivot...

    if(leftId<componentIdMaximal)
      quickSortIndexArray(array,indexArray,leftId,componentIdMaximal); //Error here?
    if(leftId>0&&--leftId>componentIdMinimal)
      quickSortIndexArray(array,indexArray,componentIdMinimal,leftId-1);
  }
  template<class ComponentSType> static Array<Count> quickSortIndexArray(Array< ComponentSType >& array,Count componentIdMinimal=0,Count componentIdMaximal=0)
  {
    if(componentIdMaximal==0)
      componentIdMaximal=array.pData->componentSCount-1;
    if(componentIdMinimal>=componentIdMaximal)
    {
      Array<Count> indexArrayEmpty;
      return indexArrayEmpty;
    }
    Array<Count> indexArray(componentIdMaximal-componentIdMinimal+1,NULL);
    for(Count i=componentIdMinimal;i<=componentIdMaximal;++i)
      indexArray[i]=i;

    quickSortIndexArray(array,indexArray,componentIdMinimal,componentIdMaximal);

    return indexArray;
  }
};


template<class ComponentSType> class MathematicaArray: public Array<ComponentSType>
{
public:
  using Array<ComponentSType>::Array;
  template<class K>MathematicaArray(Array<K> a)
  {
    this->resize(a.size());
    for(Count i=0;i<a.size();++i)
      this->operator[](i)=a[i];  //this->operator=(a); // results in recursive construction.
  }
};

template<class K> MathematicaArray(Array<K>) -> MathematicaArray< nested_t<DimensionCount<K,fundamental_t<K>>::value,MathematicaArray,fundamental_t<K>> >; // deduction guide


template<class ComponentSType> std::ostream& operator<<(std::ostream& ostream,const MathematicaArray<ComponentSType>& array)
{
  ostream<<"{";
  if(array.pData!=NULL)
  {
    for(Count i=0;i<array.size();++i)
    {
                            ostream << (array[i]);
      if(i+1!=array.size()) ostream << ",";
      if(arraySStreamOprOutSWidth>0 && ostream.tellp()<int(i)*(arraySStreamOprOutSWidth+1))
      {
        Count fillN=Count(int(i)*(arraySStreamOprOutSWidth+1)-ostream.tellp());
        for(Count i=0;i<fillN;++i)
          ostream<<" ";
      }
    }
  }
  ostream<<"}";
  return ostream;
}

/* C++11 implemantation without "if constexpr" but SFINAE
     template<class Array, class Lambda, class Parameter...>
       inline static typename enable_if< inheritsArray<Array>::value && !is_same< Array, std::tuple_elements<0,lambda_args<Lambda>::type>::type> >::value > >::value, decltype(auto)>::type
         map (Array&& array,Lambda&& lambda,Parameter... parameter)
 {
      for(auto i=0;i<array.size();++i)
        map(array[i],i,parameter...);
        lambda(array[i],i,parameter...);
    }
*/

    /*struct ArrayLoop//bzw. mapping-template
    {
      virtual void evaluate(Real& sum,Real& squareSum,Real& count,const int& dimensionId) = 0;
      virtual void evaluate(Real& real,const int& dimensionId) = 0;

      void dimensionSLoop(Real& sum,Real& squareSum,Real& count,const int& dimensionId)
      {  evaluate(sum,squareSum,count,dimensionId);  }
      void dimensionSLoop(Real& real,const int dimensionId)
      {  evaluate(real,dimensionId);  }

      template<class ComponentSType>void dimensionSLoop(Array<ComponentSType>& sumArray,Array<ComponentSType>& squareSumArray,Array<ComponentSType>& countArray,const int dimensionId=0)
      {
        for(int i=0;i<sumArray.pData->componentSCount;++i)
          dimensionSLoop(sumArray[i],squareSumArray[i],countArray[i],dimensionId+1);
      }
      template<class ComponentSType>void dimensionSLoop(Array<ComponentSType>& array,const int dimensionId=0)
      {
        for(int i=0;i<array.pData->componentSCount;++i)
          dimensionSLoop(array[i],dimensionId+1);
      }
    };
    typedef ArrayLoop Map;*/
