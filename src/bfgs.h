/*
  Copyright (C) 2008 Thorfinn Sand Korneliussen  thorfinn@binf.ku.dk
  
  This file is part of HMMld 0.801
  
  HMMld is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  Foobar is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with HMMld.  If not, see <http://www.gnu.org/licenses/>.
*/  
/*
  Copyright (C) 2008 Thorfinn Sand Korneliussen  thorfinn@binf.ku.dk
  
  This file is part of HMMld 0.987
  
  HMMld is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  Foobar is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/  
/*
  these three parameters affect precision of optimization.
  m is not recommended to be higher than 20
  decreasing factr, pgtol will increase precision
  factr is the multiple of machine precision that result will be
  pgtol is size of gradient on exit
*/
#define MVAL 30
#define FACTR 1.0e6
#define PGTOL 1.0e-4

/*
  nbd is a vector of integers of dimension numpars.
  nbd[i]=0 if there are no bounds for parameter i, 
  =1 if there are only lower bounds
  =2 if there are both lower/upper bounds
  =3 if there is only upper bound                
  or send nbd=NULL is equivalent to nbd=2 for all parameters
  
  noisy=0 => no output, noisy=1 => one line of output, noisy < 99 some output,
  noisy>=100 probably too much
  
  dfun is derivative function or send NULL to use numerical derivative
  (getgradient function)
  
*/
double findmax_bfgs(int numpars, double *invec, double (*fun)( double x[]),
		    void (*dfun)( double x[], double y[]),
		    double *lowbound, double *upbound,
		    int *nbd, int noisy);

void getgradient(int npar,  double invec[],  int need_gradient[], 
		 double outvec[], double(*func)( double []), 
		  double* lowbound,  double* upbound);
