gsl_rng *random_setup(void)
{
const gsl_rng_type *TT;
long seedy;
srand((unsigned) time(NULL));
//seedy = -rand();
seedy = time(NULL)*(int)getpid();  //JL: Including the process id to avoid the situation where the jobs submitted at the same time point have identical initial random numbers drawn.
gsl_rng_env_setup ();
gsl_rng_default_seed = seedy;
TT = gsl_rng_default;

return gsl_rng_alloc(TT);
}
//Joe's version
//gsl_rng *rand1;  //This has to be global, to ensure that the generator doesn't start over again.

//int main(){
//  long seed;
//  srand((unsigned) time(NULL));
//  seed = time(NULL)*(int)getpid();

//  gsl_rng_env_setup ();
//  gsl_rng_default_seed = seed;
//  const gsl_rng_type *T1;
//  T1 = gsl_rng_default;
//  rand1 = gsl_rng_alloc(T1);

//}

// -------------------------------------------------------------------- //

int ipow(int x, int p)
{
  if (p == 0) return 1;
  if (p == 1) return x;

  int tmp = ipow(x, p/2);
  if (p%2 == 0)             return tmp * tmp;
  else                      return x * tmp * tmp;
}
