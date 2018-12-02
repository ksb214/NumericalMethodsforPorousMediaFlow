#include<iostream>
#include<math.h>
using namespace std;

double zf( double, double );

const double Pc = 672.;
const double Tc = 383.;


int main(){
  double p, t, z; 
  
  p = 4000; 
  t = 679.67; 
 
  z = zf(p,t);  
  cout<< " z is "<< z ;
}


double zf(double p1, double t1){

//  p, psia - pressure
//  t, Rankine - temperature
//  zf - gas compressibility factor

  double *a;
  double dr, pr, tr, z, c1, c2, c3, c4, crit;
  double fun, dfun, del, dc4dr, drn;  //dr is rhor 
  a = new double [11];
  double z1;
  int i;
  if(p1 < 14.0) z1 = 1.0; //At low pressure z=1;
      
  else{

    tr = t1/Tc;
    pr = p1/Pc;
    cout<<" pr is "<<pr << '\n';
    

    // Dranchuk & Abou-Kassem Eq.
    a[1] = 0.3265 ;
    a[2] =-1.07;
    a[3] =-0.5339;
    a[4] = 0.01569;
    a[5] =-0.05165;
    a[6] = 0.5475;
    a[7] =-0.7361;
    a[8] = 0.1844;
    a[9] = 0.1056;
    a[10]= 0.6134;
    a[11]= 0.7210;

    c1 = a[1] + a[2]/tr + a[3]/pow(tr,3) + a[4]/pow(tr,4) + a[5]/pow(tr,5);
    c2 = a[6] + a[7]/tr + a[8]/pow(tr,2);
    c3 = a[9]*( a[7]/tr + a[8]/pow(tr,2) );

    crit = 1.0e-6;
    z1 = 1.0; // Initialising z1 to 1
    
  for(i = 1; i=10; i++){
  
    //  dr = 0.27*pr/(z1*tr);
    // c4 = a[10]*( 1.0 + a[11]*pow(dr,2))*(pow(dr,2)/pow(tr,3))*exp(-a[11]*pow(dr,2));
  
    // dc4dr = ( 2.*a[10]*dr/pow(tr,3) )*( 1. + a[11]*pow(dr,2) - pow( a[11]*pow(dr,2),2 )* exp( -a[11]*pow(dr,2) )); //d(c4)/dr 
    
    // fun = -z + (1. + c1*dr + c2*( pow(dr,2)) - c3*pow(dr,5) + c4);      // Function 
    // dfun = 0.27*pr/(tr*pow(dr,2)) + ( c1 + 2.*c2*dr - 5.*c3*a[9]*( pow(dr,4) ) + dc4dr) ;  // Differentiated function
    // del = -(fun/dfun);
    // drn = dr + del;
   
    cout<<i<<'\n';

    // if(fabs(del)< crit) break;
    
    // z1 = 0.27*pr/(drn*tr);
  }

  z1 = 0.27*pr/(dr*tr);
  
}

  return z1;
}
