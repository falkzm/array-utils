 #include <boost/type_traits/is_complex.hpp>



class DataAnalysis
{
public:

  struct ArrayLoop//bzw. mapping-template
  {
    virtual void evaluate(Real& sum,Real& squareSum,Real& count,const int& dimensionId __attribute__((unused)) )
    {//noetig, da direkte ueberladung von dimensionSLoop zum templateFkt ueber abgeleitete Klasse unansprechbar macht.
      sum=0; squareSum=0; count=0;
      cerr<<"ERROR: dataAnalysis::ValueLoop::dimensionSLoop\n";
    }
    virtual void evaluate(Real& real,const int& dimensionId __attribute__((unused)) )
    {//noetig, da direkte ueberladung von dimensionSLoop zum templateFkt ueber abgeleitete Klasse unansprechbar macht.
      real=0;
      cerr<<"ERROR: dataAnalysis::ValueLoop::dimensionSLoop\n";
    }

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


  template<class Type>static void mapSqrt(Type& type) // deprecated.
  // -only useful to achieve mapping inside a scope, so if recursive template functions with elementwise treatment should not be declared.
  // -in newer C++ versions, which support lambda function, mapping functionality can be achieved more elegant. See ArrayTools.h
  {
    struct :ArrayLoop
    {
      virtual void evaluate(Real& real,const int& dimensionId __attribute__ ((unused)) )
      {  real=std::sqrt(real);  }
    } arrayLoop;

    arrayLoop.dimensionSLoop(type);
  }


 template< class ScalarUser=void, class ComponentSType, keyClasses (Key,Scalar,ScalarUser,ComponentSType) > static Scalar componentsSum(const Array<ComponentSType>& array,const bool considerComponentIdIntervall=false,int componentIdMinimal=0,int componentIdMaximal=0, Key key=0)//stepsObservablesSMean sollte auf Std-Typ aufbauen
    //the const-length key access provides a more general access sheme that might become std for similar functions
    //  provide a key (Index-Array) and desired return type to enable summing over a selected subArray
  {
    if(!considerComponentIdIntervall||componentIdMinimal<0)                              componentIdMinimal=0;
    if(!considerComponentIdIntervall||componentIdMaximal>=array.pData->componentSCount)  componentIdMaximal=array.pData->componentSCount-1;

    Scalar componentsSum = array.template at< const Scalar >( array[componentIdMinimal], key );//this function is static & becomes the identity if Scalar is ComponentSType, furthermore, key is a nameless arg and we will be optimized away.
    for(int componentId=componentIdMinimal+1;componentId<=componentIdMaximal;++componentId)
      componentsSum += array.template at< Scalar >( array[componentId], key );

    return componentsSum;
  }
  /*
  template<class ComponentSType>static ComponentSType componentsSum(const Array<ComponentSType>& array,const bool considerComponentIdIntervall=false,int componentIdMinimal=0,int componentIdMaximal=0)//stepsObservablesSMean sollte auf Std-Typ aufbauen
  {
    if(!considerComponentIdIntervall||componentIdMinimal<0)
      componentIdMinimal=0;
    if(!considerComponentIdIntervall||componentIdMaximal>=array.pData->componentSCount)
      componentIdMaximal=array.pData->componentSCount-1;

    ComponentSType componentsSum=array[componentIdMinimal];
    for(int componentId=componentIdMinimal+1;componentId<=componentIdMaximal;++componentId)
      componentsSum+=array[componentId];

    return componentsSum;
  }*/


  struct ResizeValues
  {//kann leider nicht durch ValueLoop dargestellt werden, da Koerper (s.u.) virtuell+template waere, was sich gegenseitig ausschliesst (v-table-Groesse unbestimmbar)
    Array<Count> dimensionsSSize;

    void dimensionSLoop(Real& sum,Real& squareSum,Real& count,const int dimensionId __attribute__((unused))=0 )
    {
      sum       =0;
      squareSum =0;
      count     =0;
    }//RekursionsAnker

    template<class ComponentSType>void dimensionSLoop(Array<ComponentSType>& sumArray,Array<ComponentSType>& squareSumArray,Array<ComponentSType>& countArray,const int dimensionId=0)
    {
      {//Koerper
        sumArray.resize(dimensionsSSize[dimensionId]); squareSumArray.resize(dimensionsSSize[dimensionId]); countArray.resize(dimensionsSSize[dimensionId]);
      }
      for(int v=0;v<sumArray.pData->componentSCount;++v)
        dimensionSLoop(sumArray[v],squareSumArray[v],countArray[v],dimensionId+1);
    }
  };
  template<class ComponentSType>static void resize(const Count stepN,const Count observableN,Array<Count> valueSdimensionsSSize,Array< Array<ComponentSType> >& stepsObservablesSValueSum,Array< Array<ComponentSType> >& stepsObservablesSValueSquareSum,Array< Array<ComponentSType> >& stepsObservablesSValueCount)
  {
    ResizeValues resizeValues;
      resizeValues.dimensionsSSize=valueSdimensionsSSize;

    stepsObservablesSValueSum.resize(stepN); stepsObservablesSValueSquareSum.resize(stepN); stepsObservablesSValueCount.resize(stepN);
    for(int stepId=0;stepId<stepN;++stepId)
    {
      stepsObservablesSValueSum[stepId].resize(observableN); stepsObservablesSValueSquareSum[stepId].resize(observableN); stepsObservablesSValueCount[stepId].resize(observableN);
      for(int o=0;o<observableN;++o)
        resizeValues.dimensionSLoop(stepsObservablesSValueSum[stepId][o],stepsObservablesSValueSquareSum[stepId][o],stepsObservablesSValueCount[stepId][o]);
    }
  }


  /*template<class ComponentSType>void reset(Array<ComponentSType>& resultArray)
  {
    struct :ArrayLoop
    {
      virtual void dimensionSLoop(Result& result,const int dimensionId)
      {
        result.sum=0;
        result.squareSum=0;
        result.count=0;
      }
    } arrayLoop;
    arrayLoop.dimensionSLoop(resultArray);
  }*/


  //im Gegensatz zu autoCorrelation(rekursiv bis skalarTyp) wir hier jeweils nur das oberste array aufgeloest
  template<class Type>static Type calculateMean(Type valueSum,Type valueCount)//stepsObservablesSMean sollte auf Std-Typ aufbauen
  {
    valueSum/=valueCount;
    return valueSum;
  }
  template<class ComponentSType>static ComponentSType calculateMean(const Array<ComponentSType>& stepsSValueSum,const Array<ComponentSType>& stepsSValueCount,const bool considerStepIntervall=false,const int stepIdMinimal=0,const int stepIdMaximal=0)//stepsObservablesSMean sollte auf Std-Typ aufbauen
  {
    return calculateMean(componentsSum(stepsSValueSum,considerStepIntervall,stepIdMinimal,stepIdMaximal),componentsSum(stepsSValueCount,considerStepIntervall,stepIdMinimal,stepIdMaximal));
  }


  template<class ComponentSType>static void binning(Array<Count> valueSdimensionsSSize,Array< Array<ComponentSType> >& stepsObservablesSValueSum,Array< Array<ComponentSType> >& stepsObservablesSValueSquareSum,Array< Array<ComponentSType> >& stepsObservablesSValueCount,Count binCount=10)//stepsObservablesSMean sollte auf Std-Typ aufbauen
  {
    Count stepN=stepsObservablesSValueSum.pData->componentSCount;
      if(stepN<=binCount)
        return;
    Array< Array<ComponentSType> > stepsObservablesSValueSumBinned,stepsObservablesSValueSquareSumBinned,stepsObservablesSValueCountBinned;
      resize(binCount,stepsObservablesSValueSum[0].pData->componentSCount,valueSdimensionsSSize,stepsObservablesSValueSumBinned,stepsObservablesSValueSquareSumBinned,stepsObservablesSValueCountBinned);

    for(int binId=0;binId<binCount;++binId)
      stepsObservablesSValueSumBinned[binId]=calculateMean(stepsObservablesSValueSum,stepsObservablesSValueCount,true,(binId*int(stepN))/int(binCount),((binId+1)*int(stepN))/int(binCount)-1);
    //stepsObservablesSValueCountBinned=1;//wenn man statt einem mittelwert sum; sum�; count exakt zusammenf�hren wuerde, dann wuerde das ergebnis exakt identsch sein. (mittelwert, stdAbw)
    stepsObservablesSValueSquareSumBinned=stepsObservablesSValueSumBinned;
      stepsObservablesSValueSquareSumBinned*=stepsObservablesSValueSquareSumBinned;

    struct ArrayLoopRemoveNAN:ArrayLoop
    {
      virtual void evaluate(Real& sum,Real& squareSum,Real& count,const int& dimensionId __attribute__ ((unused)) )
      {
        if(count==0)
        {  sum=0; squareSum=0;  }//no data
        else
          count=1;//binned
      }
    } arrayLoop;
    arrayLoop.dimensionSLoop(stepsObservablesSValueSumBinned,stepsObservablesSValueSquareSumBinned,stepsObservablesSValueCount);
    stepsObservablesSValueSum       =stepsObservablesSValueSumBinned;
    stepsObservablesSValueSquareSum =stepsObservablesSValueSquareSumBinned;
  }


  template<class ComponentSType>static Array<ComponentSType> calculateStandardDeviation(const Array< Array<ComponentSType> >& stepsObservablesSValueSum,const Array< Array<ComponentSType> >& stepsObservablesSValueSquareSum,const Array< Array<ComponentSType> >& stepsObservablesSValueCount,Array<ComponentSType> observablesSMean)//stepsObservablesSMean sollte auf Std-Typ aufbauen
  {
    Array<ComponentSType> observablesSStandardDeviation,
                          observablesSValueSum      =componentsSum(stepsObservablesSValueSum),
                          observablesSValueSquareSum=componentsSum(stepsObservablesSValueSquareSum),
                          observablesSValueCount    =componentsSum(stepsObservablesSValueCount);

    observablesSStandardDeviation=observablesSValueSquareSum;//2. binomische Formel (fuer jeden Messwert extra);1. quadratischer Term
    observablesSValueSum*=2;//Mischterm
      observablesSValueSum*=observablesSMean;//* ist als skalarProdukt def.; *= komponentenweise
      observablesSStandardDeviation-=observablesSValueSum;
    observablesSMean*=observablesSMean;//2. quadratischer Term
      observablesSMean*=observablesSValueCount;//geht sooft wie der messwert selbst ein
      observablesSStandardDeviation+=observablesSMean;

    observablesSValueCount-=1;//Normierung der QuadratSumme mit N(N-1)
      observablesSStandardDeviation/=observablesSValueCount;

    //WurzelZiehen
    mapSqrt(observablesSStandardDeviation);

    return observablesSStandardDeviation;
  }


  template<class ComponentSType>static Array< Array<ComponentSType> > calculateStepsObservablesSJackknifeMean(const Array< Array<ComponentSType> >& stepsObservablesSValueSum,const Array< Array<ComponentSType> >& stepsObservablesSValueCount)//stepsObservablesSMean sollte auf Std-Typ aufbauen
  {//der Fehler bei N Messwerten wird abgeschaetzt, indem der Fehler aller mgl N-1 MesswertMengen (Sample) durch Differenzbildung ihres Mittelwertes zum N-Sample-Mittelwerte berechnet wird. Anschliesend quadratische Mittelwertbildung wie von StdAbw gewohnt
    Count stepN=stepsObservablesSValueSum.pData->componentSCount;

    Array< Array<ComponentSType> > stepsObservablesSJackknifeMean(stepN,NULL);
    {
      //solange hier nur rein statistische mittelwerte benutzt werden, bleibt jackknife identisch zum Fehler aus der StdAbweichung.
      //vorschlag: MittelwertBildung (Jackknife, Global auslagern: fkt(messwerte)->mittelWerte abgeleiteter Groesse
      //virtual nutzen; diese Fkt kann dann per system("command.sh") externes Skript aufrufen(, um zb Gnuplot zum fitten zu benutzen).
      //  template virtual ...  nicht mgl.
      //   =>nutze Spaltung der der JackknifeError-Funktion um samples zu manipulieren
      Array<ComponentSType> observablesSValueSum  =componentsSum(stepsObservablesSValueSum),
                            observablesSValueCount=componentsSum(stepsObservablesSValueCount);

      for(int stepId=0;stepId<stepN;++stepId)
      {
        stepsObservablesSJackknifeMean[stepId]=observablesSValueSum;
          stepsObservablesSJackknifeMean[stepId]-=stepsObservablesSValueSum[stepId];

        Array<ComponentSType> observablesSValueNormalisation=observablesSValueCount;
          observablesSValueNormalisation-=stepsObservablesSValueCount[stepId];

        stepsObservablesSJackknifeMean[stepId]/=observablesSValueNormalisation;
      }
    }

    return stepsObservablesSJackknifeMean;
  }

  template<class ComponentSType>static Array<ComponentSType> calculateJackknifeError(Array< Array<ComponentSType> > stepsObservablesSJackknifeMean,const Array<ComponentSType> observablesSMean)//stepsObservablesSMean sollte auf Std-Typ aufbauen
  {
    Count stepN=stepsObservablesSJackknifeMean.pData->componentSCount;

    Array<ComponentSType> observablesSJackknifeError;
    {
      for(int stepId=0;stepId<stepN;++stepId)
        stepsObservablesSJackknifeMean[stepId]-=observablesSMean;
      stepsObservablesSJackknifeMean*=stepsObservablesSJackknifeMean;

      observablesSJackknifeError=componentsSum(stepsObservablesSJackknifeMean);
      observablesSJackknifeError*=Real(stepN-1)/Real(stepN);
      mapSqrt(observablesSJackknifeError);
    }

    return observablesSJackknifeError;
  }

  template<class ComponentSType>static Array<ComponentSType> calculateJackknifeError(const Array< Array<ComponentSType> >& stepsObservablesSValueSum,const Array< Array<ComponentSType> >& stepsObservablesSValueCount,const Array<ComponentSType> observablesSMean)//stepsObservablesSMean sollte auf Std-Typ aufbauen
  {
    return calculateJackknifeError(calculateStepsObservablesSJackknifeMean(stepsObservablesSValueSum,stepsObservablesSValueCount),observablesSMean);//stepsObservablesSMean sollte auf Std-Typ aufbauen
  }


  template<class ComponentSType>static Array< Array< Array<ComponentSType> > > calculateBootstrapSamplesStepsObservablesSValue(const Array< Array<ComponentSType> >& stepsObservablesSValueSum,const Count bootstrapSampleN,const Count bootstrapSamplesSStepN,Array< Array< Array<ComponentSType> > >* bootstrapSamplesStepsObservablesSValueCountOutputP=NULL,const Array< Array<ComponentSType> >* stepsObservablesSValueCountP=NULL)
  {//if bootstrapSamplesStepsObservablesSValueCountOutputP && stepsObservablesSValueCountP !=NULL then the method will generate samples which represents values as sum and Count instead a of a single value
    Array< Array< Array<ComponentSType> > > bootstrapSamplesStepsObservablesSValueSum          (bootstrapSampleN,NULL);
    Count&                                  stepN                                             =stepsObservablesSValueSum.pData->componentSCount;
    bool                                    generateBootstrapSamplesStepsObservablesSValueCount=bootstrapSamplesStepsObservablesSValueCountOutputP!=NULL&&stepsObservablesSValueCountP!=NULL;
      if(generateBootstrapSamplesStepsObservablesSValueCount)
        bootstrapSamplesStepsObservablesSValueCountOutputP->resize(bootstrapSampleN);

    for(int bootstrapSampleId=0;bootstrapSampleId<bootstrapSampleN;++bootstrapSampleId)
    {
      bootstrapSamplesStepsObservablesSValueSum              [bootstrapSampleId].resize(bootstrapSamplesSStepN);
      if(generateBootstrapSamplesStepsObservablesSValueCount)
        (*bootstrapSamplesStepsObservablesSValueCountOutputP)[bootstrapSampleId].resize(bootstrapSamplesSStepN);

      for(int stepId=0;stepId<bootstrapSamplesSStepN;++stepId)
      {
        int sourceStepId=int( bootstrapSampleId==0 ? stepId : CMonteCarlo::random(stepN-1));

        bootstrapSamplesStepsObservablesSValueSum              [bootstrapSampleId][stepId]=stepsObservablesSValueSum      [sourceStepId];
        if(generateBootstrapSamplesStepsObservablesSValueCount)
          (*bootstrapSamplesStepsObservablesSValueCountOutputP)[bootstrapSampleId][stepId]=(*stepsObservablesSValueCountP)[sourceStepId];
      }
    }

    return bootstrapSamplesStepsObservablesSValueSum;
  }

  template<class ComponentSType>static Array< Array<ComponentSType> > calculateSamplesObservablesSValueMean(Array< Array< Array<ComponentSType> > > samplesStepsObservablesSValueSum,Array< Array< Array<ComponentSType> > >* samplesStepsObservablesSValueCountP=NULL)
  {
    Array< Array<ComponentSType> > samplesObservablesSValueMean(samplesStepsObservablesSValueSum.pData->componentSCount);

    for(int sampleId=0;sampleId<samplesStepsObservablesSValueSum.pData->componentSCount;++sampleId)
    {
      if(samplesStepsObservablesSValueCountP!=NULL)
        samplesObservablesSValueMean[sampleId]=calculateMean(samplesStepsObservablesSValueSum[sampleId],(*samplesStepsObservablesSValueCountP)[sampleId]);
      else
      {
        samplesObservablesSValueMean[sampleId]= componentsSum(samplesStepsObservablesSValueSum[sampleId]);
        samplesObservablesSValueMean[sampleId]/=Real(samplesStepsObservablesSValueSum[sampleId].pData->componentSCount);
        /**///say<<samplesObservablesSValueMean[sampleId][0]<<"\n";
      }
    }

    return samplesObservablesSValueMean;
  }

  template<class ComponentSType>static Array<ComponentSType> calculateObservablesSSampleError(Array< Array<ComponentSType> > samplesObservablesSValueMean,const Array<ComponentSType> observablesSValueMean)//samplesObservablesSMean sollte auf Std-Typ aufbauen
  {
    Count sampleN=samplesObservablesSValueMean.pData->componentSCount;

    Array<ComponentSType> observablesSSampleError;
    {
      for(int sampleId=0;sampleId<sampleN;++sampleId)
        samplesObservablesSValueMean[sampleId]-=observablesSValueMean;
      samplesObservablesSValueMean*=samplesObservablesSValueMean;

      observablesSSampleError=componentsSum(samplesObservablesSValueMean);
      observablesSSampleError/=Real(sampleN-1);// only difference between bootstrap and jackknife equivalent! (JackknifeSSampleN==stepN and samples are called steps therefore)
      mapSqrt(observablesSSampleError);
    }

    return observablesSSampleError;
  }

  template<class ComponentSType>static Array<ComponentSType> calculateObservablesSSampleError(Array< Array<ComponentSType> > samplesObservablesSValueMean)//samplesObservablesSMean sollte auf Std-Typ aufbauen
  {
    return calculateObservablesSSampleError(samplesObservablesSValueMean,componentsSum(samplesObservablesSValueMean)/samplesObservablesSValueMean.pData->componentSCount);
  }


  template<class ComponentSType>static void exportToFile(const Count stepN,const Count observableN,Array<Count> valueSdimensionsSSize,Array< Array<ComponentSType> >& stepsObservablesSValueSum,Array< Array<ComponentSType> >& stepsObservablesSValueSquareSum,Array< Array<ComponentSType> >& stepsObservablesSValueCount,const char* fileName="results/genericData",const char* fileNameAppendix=".dat",const char* fileAppendix=NULL)
  {//naturaly need same structure information as resize to export dataStructure beside the data

    std::stringstream fileNameStringstream;
      fileNameStringstream<<fileName<<fileNameAppendix;

    std::ofstream write(fileNameStringstream.str().c_str());
    write<<"#format 0: stepN observableN valueSdimensionN  valueSdimensionsSSize\n";

    write <<stepN<<" "
          <<observableN<<" "
          <<(valueSdimensionsSSize.pData->componentSCount)<<" ";
    for(int componentId=0;componentId<valueSdimensionsSSize.pData->componentSCount;++componentId)
      write<<" "<<(valueSdimensionsSSize[componentId]);
    write<<"\n\n";

    cout<<"dataAnalysis::export("<<fileNameStringstream.str().c_str()<<")...\n";
    struct ArrayLoopWrite:ArrayLoop
    {
      std::ofstream& write;

      virtual void evaluate(Real& sum,Real& squareSum,Real& count,const int& dimensionId)
      {
        write<<sum<<" "<<squareSum<<" "<<count<<"\n";
      }

      ArrayLoopWrite(std::ofstream& writeD)
        :write(writeD)
      {}
    } arrayLoop(write);
    arrayLoop.dimensionSLoop(stepsObservablesSValueSum,stepsObservablesSValueSquareSum,stepsObservablesSValueCount);

    if(fileAppendix!=NULL)
      write<<fileAppendix;
    write.close();
    cout<<"\ndataAnalysis::export: finished\n";
  }

  template<class ComponentSType>static void importFromFile(Count& stepN,bool readStepNAndResize,Count& observableN,Array<Count>& valueSdimensionsSSize,Array< Array<ComponentSType> >& stepsObservablesSValueSum,Array< Array<ComponentSType> >& stepsObservablesSValueSquareSum,Array< Array<ComponentSType> >& stepsObservablesSValueCount,const char* fileName="results/genericData",const char* fileNameAppendix=".dat")
  {//naturaly need same structure information as resize to export dataStructure beside the data
    std::stringstream fileNameStringstream;
      fileNameStringstream<<fileName<<fileNameAppendix;
    cout<<"dataAnalysis::import("<<fileNameStringstream.str().c_str()<<")...\n";

    std::ifstream read(fileNameStringstream.str().c_str());
    {//KommentarZeile wegwerfen
      stringbuf stringbufTrash;
      read.get(stringbufTrash);// den Rest der Zeile weglesen
    }
    {
      Count stepND;Count valueSdimensionsSSizeSComponentSCount;

      read>>stepND
          >>observableN
          >>valueSdimensionsSSizeSComponentSCount;
      valueSdimensionsSSize.resize(valueSdimensionsSSizeSComponentSCount);
      cout<<"stepND="<<stepND<<" observableN="<<observableN<<" valueSdimensionsSSizeSComponentSCount="<<valueSdimensionsSSizeSComponentSCount<<": ";
      for(int componentId=0;componentId<valueSdimensionsSSize.pData->componentSCount;++componentId)
      {
        read>>(valueSdimensionsSSize[componentId]);
        cout<<valueSdimensionsSSize[componentId]<<" ";
      }
      cout<<"\n";

      if(readStepNAndResize)
      {
        cout<<"resize...\n";
        stepN=stepND;
        resize(stepN,observableN,valueSdimensionsSSize,stepsObservablesSValueSum,stepsObservablesSValueSquareSum,stepsObservablesSValueCount);
      }
      else
        if(stepN<stepND)
        {
          cerr<<"ERROR: dataAnalysis::importFromFile: stepN<stepND\n";
          return;
        }
    }

    cout<<"read data...\n";
    struct ArrayLoopRead:ArrayLoop
    {
      std::ifstream& read;

      virtual void evaluate(Real& sum,Real& squareSum,Real& count,const int& dimensionId)
      {
        read>>sum>>squareSum>>count;
        //cout<<sum<<" "<<squareSum<<" "<<count<<"\n";
      }

      ArrayLoopRead(std::ifstream& readD)
        :read(readD)
      {}
    } arrayLoop(read);
    arrayLoop.dimensionSLoop(stepsObservablesSValueSum,stepsObservablesSValueSquareSum,stepsObservablesSValueCount);
    //cin>>trash;

    read.close();
    cout<<"\ndataAnalysis::import: finished\n";
  }


  static void check(Count stepN=1000000,Count observableN=3)
  {
    cout<<"DataAnalysis::check(): random distribution [0,1] ...\n";
    Array<Count> valueSdimensionsSSize=2;
    Array< Array< Array<Real> > > stepsObservablesSValueSum,stepsObservablesSValueSquareSum,stepsObservablesSValueCount;
      resize(stepN,observableN,valueSdimensionsSSize, stepsObservablesSValueSum,stepsObservablesSValueSquareSum,stepsObservablesSValueCount);
      for(int s=0;s<stepN;++s)
        for(int o=0;o<observableN;++o)
          for(int v=0;v<valueSdimensionsSSize[0];++v)
      {
        Real a=CMonteCarlo::random(), b=CMonteCarlo::random();
          //a=( a>0.5 ? 1 : -1);b=( b>0.5 ? 1 : -1);
        stepsObservablesSValueSum       [s][o][v]=a+b;
        stepsObservablesSValueSquareSum [s][o][v]=a*a+b*b;
        stepsObservablesSValueCount     [s][o][v]=2;
      }

    //mean
    cout<<"  mean:\n    ";
    Array< Array<Real> > observablesSMean=calculateMean(stepsObservablesSValueSum,stepsObservablesSValueCount);
    for(int o=0;o<observableN;++o)
      for(int v=0;v<valueSdimensionsSSize[0];++v)
        cout<<observablesSMean[o][v]<<" ";
    cout<<"\n";
    //stdAbw
    cout<<"  stdDeviation:\n    ";
    Array< Array<Real> > observablesSStandardDeviation=calculateStandardDeviation(stepsObservablesSValueSum,stepsObservablesSValueSquareSum,stepsObservablesSValueCount,observablesSMean);
    for(int o=0;o<observableN;++o)
      for(int v=0;v<valueSdimensionsSSize[0];++v)
        cout<<observablesSStandardDeviation[o][v]<<" ";//soll (N->inf): 1/sqrt(12) approx 0,288675
    cout<<"\n";
    cout<<"  stdDeviation:\n    ";
    Array< Array<Real> > observablesSJackknifeError=calculateJackknifeError(stepsObservablesSValueSum,stepsObservablesSValueCount,observablesSMean);
    for(int o=0;o<observableN;++o)
      for(int v=0;v<valueSdimensionsSSize[0];++v)
        cout<<observablesSJackknifeError[o][v]<<" ";//soll: statistischer Fehler
    cout<<"\n";

    //binning
    cout<<"  binning: "<<stepsObservablesSValueSum.pData->componentSCount<<" ~> ...";
    binning(valueSdimensionsSSize,stepsObservablesSValueSum,stepsObservablesSValueSquareSum,stepsObservablesSValueCount);
    cout<<stepsObservablesSValueSum.pData->componentSCount<<"\n";

    //mean
    cout<<"  mean:\n    ";
    observablesSMean=calculateMean(stepsObservablesSValueSum,stepsObservablesSValueCount);
    for(int o=0;o<observableN;++o)
      for(int v=0;v<valueSdimensionsSSize[0];++v)
        cout<<observablesSMean[o][v]<<" ";
    cout<<"\n";
    //stdAbw
    cout<<"  stdDeviation:\n    ";
    observablesSStandardDeviation=calculateStandardDeviation(stepsObservablesSValueSum,stepsObservablesSValueSquareSum,stepsObservablesSValueCount,observablesSMean);
    for(int o=0;o<observableN;++o)
      for(int v=0;v<valueSdimensionsSSize[0];++v)
        cout<<observablesSStandardDeviation[o][v]<<" ";//soll (N->inf): 1/sqrt(12) approx 0,288675
    cout<<"\n";
    cout<<"  stdDeviation:\n    ";
    observablesSJackknifeError=calculateJackknifeError(stepsObservablesSValueSum,stepsObservablesSValueCount,observablesSMean);
    for(int o=0;o<observableN;++o)
      for(int v=0;v<valueSdimensionsSSize[0];++v)
        cout<<observablesSJackknifeError[o][v]<<" ";//soll: statistischer Fehler
    cout<<"\n";
  }



 // standard template algorithm section: these are implementations for data analysis standard task.  They do not use the mutual data representation/ordering via arrays used above.
 //
 // This class might be splitted into a Tool/standard algorithm file and the analysis framework above.
 // sqrt(Array<.>) and conj(Array<.>) might be moved to a more basic math namespace.



  template<class RUser=void, class F, keyClasses(K,R,RUser,F)> static R sum(const Array<F> a,Count componentIdMinimal=-1,Count componentIdMaximal=-1, K   key=0)
  {  return componentsSum<R>(a,componentIdMinimal>=0&&componentIdMaximal>=0,componentIdMinimal,componentIdMaximal,key);  }
  template<class R>                                            static R sum(const R        a,Count componentIdMinimal=-1,Count componentIdMaximal=-1, int key=0)
  {  return a;  }


  template<class RUser=void, class F, keyClasses(K,R,RUser,F)> static R mean(const Array<F> a,Count componentIdMinimal,Count componentIdMaximal, K key=0)
  {  return sum<R>(a,componentIdMinimal,componentIdMaximal,key)/(componentIdMaximal-componentIdMinimal+1);  }

  template<class RUser=void, class F, keyClasses(K,R,RUser,F)> static R mean(const Array<F> a, K   key=0)
  {  return sum<R>(a,-1,-1,key)/double(a.n());  }
  template<class R>                                            static R mean(const R        a, int key=0)
  {  return a;  }

  template<class RUser=void, class F, keyClasses(K,R,RUser,F)> static R corr(const Array<F> a,const Array<F> b, K key=0)
  {
    auto n=min(a.n(),b.n());
    R d=0;
    for(auto i=0;i<n;++i)
      d+=a.template at<R>(a[i],key)*a.template at<R>(b[i],key);

    return d/R(n);
  }


  template<class F> static Array<F> autoCorr(const Array<F> a,Array<F> *pSD=NULL,bool periodic=false,bool normaliseViaSd=true)
  {
    const Count n(a.n());
    Array<F> ac(n,0,NULL);

    for(auto t=0;t<n;++t)
      for(auto i=0;i<(periodic ? n : n-t );++i)
        ac[t]+=a[i]*a[(i+t)%n]/F( periodic ? n : n-t );
    if(pSD!=NULL)
    {
      Array<F> &sd=*pSD;
        sd=Array<F>(n,0,NULL);

      for(auto t=0;t<n;++t)
        for(auto i=0;i<(periodic ? n : n-t );++i)
          sd[t]+=pow((a[i]*a[(i+t)%n]-ac[t])/F( periodic ? n : n-t ),2);

      sd=sqrt(std::move(sd));
      if(normaliseViaSd)
        sd/=ac[0];
    }

    if(normaliseViaSd)
      return ac/=ac[0];
    return ac;
  }


  template<class F> static Array<F> conj(Array<F> a)
  {
    if constexpr(!is_complex_v<std::decay_t<fundamental_t<F>>>)
      for(Count i=0;i<a.size();++i)
        a[i]=conj(a[i]);
    return a;
  }
  template<class F, class=enable_if_notInheritsArray_t<F>> static complex<F> conj(complex<F> f)
  {  return std::conj(f);  }
  template<class F, class=enable_if_notInheritsArray_t<F>> static F conj(F f)
  {  return f; }

  template<class F> static Array<F> imag_complex(Array<F> a)
  {
    if constexpr(!is_complex_v<std::decay_t<fundamental_t<F>>>)
      for(Count i=0;i<a.size();++i)
        a[i]=imag_complex(a[i]);
    return a;
  }
  template<class F, class=enable_if_notInheritsArray_t<F>> static F imag_complex(F f)
  {  return std::imag(f);  }

  template<class F> static auto imag(Array<F> a)
  {
    if constexpr(!is_complex_v<std::decay_t<fundamental_t<F>>>)
      return a;
    else
    {
      Array<decltype(imag(std::declval<F>()))> r(a.n(),NULL);
      for(Count i=0;i<a.size();++i)
        r[i]=imag(a[i]);
      return r;
    }
  }
  template<class F, class=enable_if_notInheritsArray_t<F>> static auto imag(F f)
  {  return std::imag(f);  }

  template<class F> static auto real(Array<F> a)
  {
    if constexpr(!is_complex_v<std::decay_t<fundamental_t<F>>>)
      return a;
    else
    {
      Array<decltype(real(std::declval<F>()))> r(a.n(),NULL);
      for(Count i=0;i<a.size();++i)
        r[i]=real(a[i]);
      return r;
    }
  }
  template<class F, class=enable_if_notInheritsArray_t<F>> static auto real(F f)
  {  return std::real(f);  }

  template<class F> static auto real_complex(Array<F> a)
  {
    if constexpr(!is_complex_v<std::decay_t<fundamental_t<F>>>)
      for(Count i=0;i<a.size();++i)
        a[i]=real_complex(a[i]);
    return a;
  }
  template<class F, class=enable_if_notInheritsArray_t<F>> static auto real_complex(F f)
  {  return std::real(f);  }

  template< class RFundamentalUser=void, class F, class RFundamental=nonVoid_t< RFundamentalUser, std::decay_t< decltype(std::abs(std::declval< fundamental_t<F> >())) > > >
       static auto abs(F f)                                                       // this abs template for tensors automatically replaces the underlying type. If the type remains unchanged, no memory reallocation will occur, ie. the array is reused.
  {                                                                               // Should serve as prototype for remaining funcs here.
    if constexpr(!inheritsArray_v<F>)
      return RFundamental(std::abs(f));
    else
    {
      if constexpr( is_same_v< fundamental_t<F>, std::remove_reference_t<RFundamental> > )
      { // same type
        for(Count i=0;i<f.size();++i)
        {
          auto &f_i=f[i]; //saves one non-const-[]-access (which is slower). abs(static_cast<const F&>(f)[i]) is no solution since this function is abs and the recursive call will fail in this line: the then constant subfieldelement is going to be reassigned.
          f_i=abs<RFundamental>(f_i);
        }
        return f;
      }
      else
      { // replace type;
        Array< decltype(abs<RFundamental>(f[0])) > r(f.n(),NULL);
        for(Count i=0;i<f.size();++i)
          r[i]=abs<RFundamental>(static_cast<const F&>(f)[i]);
        return r;
      }
    }
  }

  template<class F> static F abs_complex(F f)
  {
    return abs< fundamental_t<F> >(f);
    /*if constexpr(!inheritsArray<A>::value)
      return std::abs(a);
    for(Count i=0;i<a.size();++i)
      a[i]=abs(a[i]);
    return a;*/
  }

  template<class F> static Array<F> sqrt(Array<F> a)
  {
    for(Count i=0;i<a.size();++i)
      a[i]=sqrt(a[i]);
    return a;
  }
  template<class F, class=enable_if_notInheritsArray_t<F>> static F sqrt(F f)
  {  return std::sqrt(f);  }

  template<class F> static Array<F> pow(Array<F> a,auto p)
  {
    for(Count i=0;i<a.size();++i)
      a[i]=pow(a[i],p);
    return a;
  }
  template<class F, class=enable_if_notInheritsArray_t<F>> static F pow(F f,auto p)
  {  return std::pow(f,p);  }


  template<class F=double> static F sd(const Array<F> a) // standard deviation for algebraic field F
                                                         //   F could be an Array<.>, too.
  {
    F e=componentsSum(a);
      e*=conj(e)/-a.n();
      for(auto i=0;i<a.n();++i)
      {
        F a_i_conj=conj(a[i]);
        e+=(a_i_conj*=a[i]); // due to * beeing the scalar product (should be changed), using the non-const return of conj and the elementwise *= instead, enables elementwise support for non-scalar F
      }

    return sqrt(e/a.n());
  }
  template<class F,class T> static F sd(const Array<T> a,Array<Count> indexArray,Count dimId,bool resetIndexArrayAtDim=true) // In general, subset arrays might be cheaper and there is no need of a 2nd implementation of an algorithm.
                                                                                                                              // If components are (large?) subarrays by themself: even wo. reference typed subarrays, constness prevents copying.
  {
    //auto lambda = [&dimId,&indexArray,&lambda](auto array,Count i=0) {  return ( i<dimId-1 ? 1 /*lambda(array[indexArray[i]],i-1)*/ : array.n() );  } -> Count;
    Count &i = indexArray[dimId],
          i0 = ( resetIndexArrayAtDim ? 0 : i ),
           n = a.getSize(indexArray)[dimId];

    F e=0;
      for( i=i0;i<n;++i )
        e += a[indexArray];
      e *= -conj(e)/(n-i0);
      for( i=i0;i<n;++i )
        e += a[indexArray]*conj(a[indexArray]);

    return sqrt(e/(n-i0));
  }


  template<class T> static Count median(const Array<T> a)
  {
    T sumAll=DataAnalysis::sum(a),
      sum_i =0;
    Count i=0;
    for(;i<=a.n()&&(sum_i+=a[i])<=sumAll/2; ++i);
    return std::max(0,i-1);
  }

};
