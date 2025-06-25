
//static double RateMovement = 1;
//static double RateMovement_toN = 1;
static double RateInteraction = 1;
static double RateDissociation = 1;
//static double RateInhibition = 1;
//static double RateInput = 1;
static double CostTrascriptMAPK14 = 1;
static double RateTrascriptMAPK = 1;
static double CostTrascriptMAPK13 = 1;
static double max_Atac = 1;
static double n = 1;
static double constant = -1;

void read_constant(string fname, double& Ka)
{
  ifstream f (fname);
  string line;
  if(f.is_open())
  {
    int i = 1;
    while (getline(f,line))
    {
      switch(i)
      {
      case 1:
        Ka = stod(line);
        //cout << "p" << i << ": " << line << "\t" << p1 << endl;
        break;
      }
      ++i;
    }
    f.close();
  }
  else
  {
    std::cerr<<"\nUnable to open " << fname << ": file do not exists\": file do not exists\n";
    exit(EXIT_FAILURE);
  }
}

void init_data_structures()
{
  read_constant("./DissGeneration", RateDissociation);
  read_constant("./InterGeneration", RateInteraction );
  read_constant("./RateMAPK", RateTrascriptMAPK);
  read_constant("./CostMAPK14", CostTrascriptMAPK14);
  read_constant("./CostMAPK13", CostTrascriptMAPK13);
  read_constant("./maxAtac", max_Atac);
  //read_constant("./MovementParam", RateMovement);
  //read_constant("./MovementParam_toN", RateMovement_toN);
  //read_constant("./RateInhibition", RateInhibition );
  //read_constant("./InputATAC", RateInput );
  
  constant = 1; 
}

double DissRate(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{
  if( constant == -1)   init_data_structures();
  
  double intensity = 1.0;
  
  for (unsigned int k=0; k<Trans[T].InPlaces.size(); k++)
  {
    intensity *= pow(Value[Trans[T].InPlaces[k].Id],Trans[T].InPlaces[k].Card);
  }
  
  double rate = RateDissociation * intensity;
  
  return(rate);
}

double IntRate(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{
  if( constant == -1)   init_data_structures();
  
  double intensity = 1.0;
  
  for (unsigned int k=0; k<Trans[T].InPlaces.size(); k++)
  {
    intensity *= pow(Value[Trans[T].InPlaces[k].Id],Trans[T].InPlaces[k].Card);
  }
  
  double rate = RateInteraction * intensity;
  
  return(rate);
}

double TrascriptMAPK14(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{
  if( constant == -1)   init_data_structures();
  
  double intensity = 1.0;
  
  int idxATAC = NumPlaces.find("Atac") -> second ;
  
  double ATAC = Value[idxATAC];
  //double rate = RateTrascriptMAPK - (RateTrascriptMAPK* 0.9 * pow(ATAC,2)/(pow(max_Atac/2,2) + pow(ATAC,2) ) );
  double rate = 0.0022 - CostTrascriptMAPK14 * ATAC;
  if( rate < 0.0 ) rate = 0.0;
  
  return(rate);
}

double TrascriptMAPK13(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{
  if( constant == -1)   init_data_structures();
  
  double intensity = 1.0;
  
  int idxATAC = NumPlaces.find("Atac") -> second ;
  
  double ATAC = Value[idxATAC];
  // double rate = RateTrascriptMAPK + (RateTrascriptMAPK *10* pow(ATAC,2)/(pow(max_Atac/2,2) + pow(ATAC,2) ) );
  double rate = RateTrascriptMAPK + CostTrascriptMAPK13 * ATAC;
  
  return(rate);
}

/*double Movement(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
 {
 if( constant == -1)   init_data_structures();
 
 double intensity = 1.0;
 
 for (unsigned int k=0; k<Trans[T].InPlaces.size(); k++)
 {
 intensity *= pow(Value[Trans[T].InPlaces[k].Id],Trans[T].InPlaces[k].Card);
 }
 
 double rate = RateMovement * intensity;
 
 
 return(rate);
 }
 
 double Movement_toN(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
 {
 if( constant == -1)   init_data_structures();
 
 double intensity = 1.0;
 
 for (unsigned int k=0; k<Trans[T].InPlaces.size(); k++)
 {
 intensity *= pow(Value[Trans[T].InPlaces[k].Id],Trans[T].InPlaces[k].Card);
 }
 
 double rate = RateMovement_toN * intensity;
 
 
 return(rate);
 }
 
 double Activation(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
 {
 if( constant == -1)   init_data_structures();
 
 
 
 map<string, int>::iterator itr;
 for (itr = NumPlaces.begin(); itr != NumPlaces.end(); ++itr) {
 if( Value[itr->second] < 0 ) {
 cout << "Value[" << itr->first  << " ]: " << Value[itr->second] << endl;
 Value[itr->second] = 0;
 } 
 }
 
 
 int idxMapk = NumPlaces.find("Mapk13acH3K14") -> second ;
 int idxPrkd1X = NumPlaces.find("Prkd1x") -> second ;
 
 double Mapk = Value[idxMapk];
 double Prkd1X = Value[idxPrkd1X];
 
 double rate = RateInhibition * Prkd1X / ( 1 + Mapk );
 
 return(rate);
 }
 
 
 double Input_ATAC(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
 {
 if( constant == -1)   init_data_structures();
 
 double rate = RateInput;
 
 return(rate);
 }
 
 */