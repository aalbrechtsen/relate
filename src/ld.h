/*
  Copyright (C) 2008 Thorfinn Sand Korneliussen  thorfinn@binf.ku.dk
  
  This file is part of HMMld 0.83
  
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
  
typedef struct {
  int count[9];
  double *expt;
  int total;
  int a;
  int b;
  int c;
  int d;
  int e;
  int f;
  int g;
  int h;
  int i;
  int j;
  int k;
  int l;
  int m;
  int a0;
  int a1;
  int a2;
  int a3;
  double p;
  double u;
  double v;
  double w;
  double x;
  double bigD;
  double rsq2;
  double dprime;
  double lod;
  /* */
  double post_p;
  signed char sign_of_r; /* sign_of_r is really just sign of bigD, but eventually
			 we might want to separate/hide all the internal values
		      */
} geno_count;

typedef geno_count *geno_cptr;  

/*
  given two snp.matrix position and length, read along the vector
  and populates the initial inputs in a new geno_count object
  and returns it.

  caller is responsible for free()
 */
geno_cptr get_geno_count(int *pos1,int *pos2, int length);
