// Code written by Kamalesh Bhambare
#include<iostream>
#include<fstream>
#include<math.h>
using namespace std;

// Declaring constant data
const int NX=5;
const int NY=5;
const float PI=3.14159265;
const int MAX_IT=100000;
const float TOL=1e-8;

// Defining the structure for storing variables at each grid point
struct Element
{
	double x, y, dx, dy; //x, y corresponds to CV x and y (All are in Ft)
	double struc, phi, depth, kx, ky; //kx and ky are in md
	int domain, well;
	double Ax, Ay;
	double G; // Net depth
	double P; //Pressure
	double bc; //To specify flow/pressure BC at well Flow in STB/D, Psf in psi
	double rw; //radius of well, 0 for no well
} ;

// Well = 0 for no well, -1 for producer, 1 for injector
// domain = 1 inside domain, 0 for outside domain

// Number of nodes
int nx, ny, n;

//Fluid Properties
double myu, B, rho, g, gc;

void props(void);
void display_props(void);
void init (Element *);
void print_node (Element *);
void read_input(Element *);
int check_conv(Element *, double *);
void thomas (double *, double * , double *, double *, double *, int );

int main(){
	
	int i,j,k;	
	int left, right, top, bot;
	int no_node, it, flag;

	double G_term, sigma, k_bar, re;
	double s; //slip factor
	double residual, max_R, maxP;

	Element *node;
	double *C, *W, *N, *E, *S, *Q;
	double *as, *bs, *cs, *ds, *x;

	double *newP;

	//Setting the fluid properties
	props();
	display_props();
	
	// Defining array of structure and pointing starting address of array to the pointer node
	node= new Element[n];
	init(node);
	
	C = new double[n+1];
	W = new double[n+1];
	N = new double[n+1];
	E = new double[n+1];
	S = new double[n+1];
	Q = new double[n+1];
	
	// Initialising arrays to zero
	for( i=0; i<=nx; i++)
		C[i] = W[i] = N[i] = E[i] = S[i] =Q[i]= 0.0;
			
	//Reading input from the file (input.dat)
	read_input(node);

	//Print node to file (for testing & validating);
	print_node(node);
	
	ofstream outfile1;
  	outfile1.open("Coefficients.dat",ios::out);
	
	outfile1<<"Cell     E        W        N      S      C      Q\n";
//	outfile1.precision(3);
	//Calculating all terms in the equation for each cell
	for (j=1; j<=ny; j++)
	{
		for (i=1; i<=nx; i++)
		{
			k = i + (j-1)*nx;
			left = (k-1);
			right = (k+1);
			top = i + j*nx;
			bot = i + (j-2)*nx;

			//Grid Boundary
			if(i==1) left = 0;
			if(i==nx) right = 0;
			if(j==1) bot = 0;
			if(j==ny) top = 0;

			//Calculating transmissivity term in each direction

			// For non existant node
			if(node[k].domain == 0 ) {
				C[k] = S[k] = W[k] = E[k] = N[k] = Q[k] = 0;
				outfile1<<k<<'\t'<<E[k]<<'\t'<<W[k]<<'\t'<< N[k]<<'\t'<< S[k]<<'\t'<< C[k]<<'\t'<<Q[k]<<'\n';
				continue;
				}

			//For existant node
			//printf("%d node exist\n",k);

			//Checking of any nbhd node doesn't exist
			if(node[left].domain==0) W[k] = 0;
				else W[k] = (2*node[k].Ax*node[left].Ax*node[k].kx*node[left].kx )/( (myu*B)*((node[k].Ax*node[k].kx*node[left].dx) + (node[left].Ax*node[left].kx*node[k].dx)) );

			if(node[top].domain==0) N[k] = 0;
				else N[k] = (2*node[k].Ay*node[top].Ay*node[k].ky*node[top].ky )/( (myu*B)*((node[k].Ay*node[k].ky*node[top].dy) + (node[top].Ay*node[top].ky*node[k].dy)) );

			if(node[right].domain==0) E[k] = 0;
				else E[k] = (2*node[k].Ax*node[right].Ax*node[k].kx*node[right].kx )/( (myu*B)*((node[k].Ax*node[k].kx*node[right].dx) + (node[right].Ax*node[right].kx*node[k].dx)) );

			if(node[bot].domain==0) S[k]=0;
				else S[k] = (2*node[k].Ay*node[bot].Ay*node[k].ky*node[bot].ky )/( (myu*B)*((node[k].Ay*node[k].ky*node[bot].dy) + (node[bot].Ay*node[bot].ky*node[k].dy)) );


			//Calculation of 'C' and 'Q' term

			//curly bracket term in the notes
			G_term = ( (g*rho)/(144*gc) )*( E[k]*(node[right].G - node[k].G) - W[k]*(node[k].G - node[left].G) + N[k]*(node[top].G - node[k].G) - S[k]*(node[k].G - node[bot].G) );
			if(node[k].domain == 0) { C[k] = 0.0; Q[k] = 0.0; }
				else {
						if(node[k].well == 0) { // No well
						C[k] = -(W[k] + N[k] + E[k] + S[k]);
						Q[k] = G_term;
						}
						else {
							if( (node[k].well == 11) || (node[k].well == -11) ) { //checking for flow rate BC
								C[k] = -(W[k] + N[k] + E[k] + S[k]);
								if (node[k].well == 11 ) Q[k] = - node[k].bc + G_term; //injector
									else Q[k] = node[k].bc + G_term; //producer
							}//End of if block (flow rate BC)
							if( (node[k].well == 12) || (node[k].well == -12) ) { //checking for sand face pressure
								k_bar = sqrt( node[k].kx*node[k].ky );
								re = 0.28*sqrt( sqrt(node[k].kx/node[k].ky)*(node[k].dx*node[k].dx) + sqrt(node[k].ky/node[k].kx)*(node[k].dy*node[k].dy) );
								re = re/( pow( node[k].ky/node[k].kx, 0.25) + pow( node[k].kx/node[k].ky, 0.25 ) );
								s = -0.5;
								sigma = 2*PI*k_bar*node[k].depth/( myu*B*( log(re/node[k].rw) + s ) ) ;

								C[k] = -(W[k] + N[k] + E[k] + S[k] + sigma);
								Q[k] = -sigma*node[k].bc + G_term;

							}

						}//End of else of Well BC
					}//End of else of domain check if loop



			outfile1<<k<<'\t'<<E[k]<<'\t'<<W[k]<<'\t'<<N[k]<<'\t'<<S[k]<<'\t'<<C[k]<<'\t'<<Q[k]<<'\n';

		}// End of i loop
	}// End of j loop

	outfile1.close();

//Solving equations using LSOR Method (Thomas algorithm will be called at each j)

	//swapping in y-direction (solving for constant 'j'
	// Creating one d arrays to store digonals of TDMA
	as = new double[nx+1];
	bs = new double[nx+1];
	cs = new double[nx+1];
	ds = new double[nx+1];
	x = new double[nx+1];
	newP = new double[n+1];	

// Initialising arrays to zero
for( i=0; i<=nx; i++)
		as[i] = bs[i] = cs[i] = ds[i] = x[i] = 0.0;	
	
	it = 0;
	while(1)
	{
		it++;
		//cout<<"Iteration number is: "<< it<<'\n';
		for( i=0; i<=n; i++)
			newP[i] = 0.0;

		for (j=1; j<=ny; j++)
		{
			no_node = 0;
			for (i=1; i<=nx; i++)
			{
				k = i + (j-1)*nx;

				if (node[k].domain == 0) continue;

				left = (k-1);
				right = (k+1);
				top = i + j*nx;
				bot = i + (j-2)*nx;

				//Grid Boundary
				if(i==1) left = 0;
				if(i==nx) right = 0;
				if(j==1) bot = 0;
				if(j==ny) top = 0;

				as[no_node+1] = W[k];
				bs[no_node+1] = C[k];
				cs[no_node+1] = E[k];
				ds[no_node+1] = Q[k] - S[k]*newP[bot] - N[k]*node[top].P ;

				no_node++;

			}// end of i loop

		thomas(as, bs, cs, ds, x, no_node);

			
			for (i=nx; i>=1; i--)
			{
				k = i + (j-1)*nx;
				if (node[k].domain == 0) { newP[k] = 0; }
					else {
						newP[k] = x[no_node];
						no_node--;
					}// end of else
			}

		
			if(no_node != 0) cout<<"Some issues in solver, data may be incorrectly mapped\n";

		}//end of j loop
	
	flag = check_conv( node, newP); 

		if( flag==1 || it == MAX_IT) break; 

		for( i=0; i<=n; i++)
			node[i].P = newP[i];

	}//End of while loop (solver loop)

	for( i=0; i<=n; i++)
			node[i].P = newP[i];
	
	
	//  Finding maximum pressure in the domain
	ofstream outfile2;
  	outfile2.open("MaxP.dat",ios::out);
	for( i=1; i<=n; i++)
	{
		//outfile2 << node[i].x <<'\t'<<node[i].y<<'\t' <<node[i].P<<'\t'<< i<<'\n';
		if( (maxP - node[i].P) < 0 ) maxP = node[i].P;
	}
	outfile2<<"Maximum Pressure in the domain is: "<<maxP;
	outfile2.close();
	
	
	//Writting data for Gnuplot 
	ofstream outfile3;
	outfile3.open("gnu_result.dat",ios::out);
	for (j=1; j<=ny; j++)
	{
		for (i=1; i<=nx; i++)
		{
			k = i + (j-1)*nx;
			outfile3<<node[k].x<<'\t'<<node[k].y<<'\t'<<node[k].P<<'\n';
		}
		outfile3<<'\n';
	}
	outfile3.close();


//Material Balance check (Residuals)
	ofstream outfile4;
	outfile4.open("Residual.dat",ios::out);
	max_R = 0;
	for (j=1; j<=ny; j++)
	{
		for (i=1; i<=nx; i++)
		{
			k = i + (j-1)*nx;
			left = (k-1);
			right = (k+1);
			top = i + j*nx;
			bot = i + (j-2)*nx;

			//Grid Boundary
			if(i==1) left = 0;
			if(i==nx) right = 0;
			if(j==1) bot = 0;
			if(j==ny) top = 0;

			if(node[k].domain !=0) {
				residual = fabs(Q[k] - S[k]*node[bot].P - W[k]*node[left].P - C[k]*node[k].P - E[k]*node[right].P - N[k]*node[top].P);
				outfile4<<k<<'\t'<<residual<<'\n';
				if(residual > max_R ) max_R = residual;
				}
		}//End of i loop
	}//End of j loop
	outfile4<<"Max Residual is;"<<max_R;
	outfile4.close();
	
	
	//Writing pressure at each cell 
	
	ofstream outfile6, outfile7,outfile8, outfile5; 
	outfile5.open("Pdata.dat",ios::out);
	outfile5.precision(4);
	//outfile5.width(5);
	outfile5<<"Cell"<<"Pressure"<<'\n';
	for( i=1; i<=n; i++)
		if( node[i].domain == 0 ) outfile5<<i<<'\t'<<0<<'\n';
			else outfile5<<i<<'\t'<<node[i].P<<'\n';
	outfile5.close();

	//Writing Production and Injection Summary
	outfile6.open("WellSummary.dat",ios::out);
		outfile6.precision(4);
		outfile6.width(30);
	outfile7.open("Injection.dat",ios::out);
		outfile7.precision(4);
		outfile7.width(30);
	outfile8.open("Production.dat",ios::out);
		outfile8.precision(4);
		outfile8.width(30);
	
	outfile6<<"Type        rw(ft)    k(md)    re(ft)    h(feet)       Sigma(bbl/day-psi)\n";
	outfile7<< "Injection(STB/day)    P(psia)     Psf(psia)\n";
	outfile8<<"Production(STB/day)   P(psia)     Psf(psia)\n";

	for( k=1; k<=n; k++)
	{
		if( node[k].well != 0) {
			k_bar = sqrt( node[k].kx*node[k].ky );
			re = 0.28*sqrt( sqrt(node[k].kx/node[k].ky)*(node[k].dx*node[k].dx) + sqrt(node[k].ky/node[k].kx)*(node[k].dy*node[k].dy) );
			re = re/( pow( node[k].ky/node[k].kx, 0.25) + pow( node[k].kx/node[k].ky, 0.25 ) );
			s = -0.5;
			sigma = 2*PI*k_bar*node[k].depth/( myu*B*( log(re/node[k].rw) + s ) ) ;

			//Injector with flow rate specified
			if (node[k].well == 11) {
				outfile6<<"Injector"<<'\t'<<node[k].rw<<'\t'<< k_bar/.001127<<'\t'<<re<<'\t'<< node[k].depth<<'\t'<<sigma<<'\n';
				outfile7<<node[k].bc<<'\t'<<node[k].P<<'\t'<<(node[k].P + node[k].bc/sigma)<<'\n';
			}
			//Injector with Psf specified
			if (node[k].well == 12) {
				outfile6<<"Injector"<<'\t'<<node[k].rw<<'\t'<< k_bar/.001127<<'\t'<< re<<'\t'<<node[k].depth<<'\t'<< sigma<<'\n';
				outfile7<<-sigma*(node[k].P - node[k].bc)<<'\t'<< node[k].P<<'\t'<<node[k].bc<<'\n';
			}

			//Producer with flow rate specified
			if (node[k].well == -11) {
				outfile6<<"Producer"<<'\t'<<node[k].rw<<k_bar/.001127<<'\t'<<re<<'\t'<< node[k].depth<<'\t'<<sigma<<'\n';
				outfile8<<node[k].bc<<'\t'<<node[k].P<<'\t'<<(node[k].P - node[k].bc/sigma)<<'\n';
				}
			
			//producer with Psf specified
			if (node[k].well == -12) {
			outfile6<<"Producer"<<'\t'<< node[k].rw<<'\t'<< k_bar/.001127<<'\t'<< re<<'\t'<< node[k].depth<<'\t'<< sigma<<'\n';
			outfile8<<sigma*(node[k].P - node[k].bc)<<'\t'<< node[k].P<<'\t'<<node[k].bc<<'\n';
				}

		}//End of if for node[k].well
	}
	outfile6.close();
	outfile7.close();
	outfile8.close();
}


// Function definition starts
void props(void)
{
	//Property Value for HW
	nx = NX;
	ny = NY;
	n = nx*ny;
	myu = 0.93; //cp
	rho = 62; //lb/ft^3
	B = 1.0; //RB/STB
	g = 32.2; //Gravity
	gc = 32.2;
}


void display_props()
{
	cout<<"No. of elements in x and y direction are "<< nx<<" and "<< ny <<'\n';
	cout<<"Total number of elements are "<< n<<'\n';
	cout<<"Fluid Properties are, myu: "<<myu <<" cp"<<", rho: "<<rho<< " lb/ft3"<< ", B: "<< B<<" RB/STB"<<'\n';
}


void init (Element *node)
{
	int i;

	for ( i=0; i<=n; i++ )
	{
		node[i].x = node[i].y = 0;
		node[i].dx = node[i].dy = 0;
		node[i].struc = node[i].phi = node[i].depth = 0;
		node[i].kx = node[i].ky = 0;
		node[i].domain = node[i].well = 0;
		node[i].Ax = node[i].Ay = 0;
		node[i].G = node[i].P = 0;
	}
}


void read_input(Element *node)
{
	int i,j,k;

	ifstream myfile;
  	myfile.open("input_validation.dat",ios::in);

	for (j=1; j<=ny; j++)
	{
		for (i=1; i<=nx; i++)
		{
			//Calculating the node number corresponding to i,j
			k = i + (j-1)*nx;
			//kx and ky are in md
			myfile >> node[k].x >> node[k].y >> node[k].dx >> node[k].dy >> node[k].domain >> node[k].well >> node[k].struc >> node[k].depth >> node[k].phi >> node[k].kx >> node[k].ky >> node[k].bc >> node[k].rw;

			if (node[k].domain == 1) {
				//Converting md to perm
				node[k].kx = node[k].kx*1.127*0.001;
				node[k].ky = node[k].ky*1.127*0.001;

				node[k].phi = node[k].phi*.01;

				//Calculating other parameters related to cell
				node[k].G = node[k].struc + node[k].depth/2.0 ;
				node[k].Ax = node[k].dy * node[k].depth;
				node[k].Ay = node[k].dx * node[k].depth;
			}//End of if loop
		}// End of i loop
	}//End of j loop


	myfile.close();
}

void print_node(Element *node)
{
	int i;
	ofstream outfile;
  	outfile.open("FormatedDataInput.txt",ios::out);
	
	outfile << "X         Y       dx    dy   domain     well   struc  depth    phi      kx       ky      Ax      Ay     G     P    WellBCs     rw"<<'\n';
	for (i=1; i<=n; i++ )
	{
		outfile<<node[i].x<<'\t'<<node[i].y<<'\t'<< node[i].dx<<'\t'<< node[i].dy<<'\t';
		outfile<<node[i].domain<<'\t'<<node[i].well<<'\t'<< node[i].struc<<'\t'<< node[i].depth<<'\t'<< node[i].phi<<'\t';
		outfile<<node[i].kx<<'\t'<< node[i].ky<<'\t';
		outfile<<node[i].Ax<<'\t'<<node[i].Ay<<'\t';
		outfile<<node[i].G<<'\t'<<node[i].P<<'\t'<< node[i].bc<<'\t'<< node[i].rw<<'\n';

		if (node[i].well < 0 ) cout<<"Cell number "<<i<<" is a producer well"<<'\n';
		if (node[i].well > 0 ) cout<<"Cell number "<<i<<" is an injector well"<<'\n';

	}

	outfile.close();
}


int check_conv(Element *node, double *newP)
{
	for(int i=1; i<=n; i++)
	if( fabs(node[i].P - newP[i]) > TOL )return 0;
	return 1;

}

void thomas (double *a, double *b, double *c, double *d, double *x, int size)
{
	int i;

	double *w, *q;

	w = new double [size+1];
	q = new double [size+1];

	//Initializating the result
	for( i=0; i<=size; i++)
	{
		w[i] = q[i] = x[i] = 0;
	
	}

	w[1] = c[1]/b[1];
	q[1] = d[1]/b[1];

	for( i=2; i<=(size-1); i++)
	{
		w[i] = c[i]/(b[i] - a[i]*w[i-1]) ;
		q[i] = (d[i] - a[i]*q[i-1])/(b[i] - a[i]*w[i-1]);
	}

	q[size] = (d[size] - a[size]*q[size-1])/(b[size] - a[size]*w[size-1]) ;

	//Solution
	x[size] = q[size];
	for( i=size-1; i>0; i-- )
	{
		x[i] = q[i] - w[i]*x[i+1] ;
	
	}

	
	delete [] q;
	delete [] w;


}
