// Code written by Kamalesh Bhambare, ksb214@psu.edu

#include<iostream>
#include<fstream>
using namespace std;
int main(){
  	float pft=0, pp, sum=0, counter=0;
  	bool flag=false;
  	double A[30];                  // Remember in cpp first element of array is A[0] and not A[1]. 
    for(int i=0;i<30;i++)
  	A[i]=0;
  // Purchase price of the oil 
  	pp = 75 ;               
  

  //Reading the price data from file data.txt
  	ifstream myfile;
  	myfile.open("data.txt",ios::in);
  	for(int i=0;i<30;i++)
    {
      myfile >> A[i];
    }
  	myfile.close();

  // Writting the program output to a file called Output.txt	

	ofstream outfile;
	outfile.open("output.txt",ios::out);
	outfile<<"Selling start day"<<'\t'<<'\t'<<'\t'<<"Selling end day"<<'\t'<<'\t'<<'\t'<<'\t'<<"Profit Earned (US $)"<<'\n'; 
		
  	for(int i=0;i<30;i++){
    pft=(A[i]-pp)*1000;
    
    
    if (pft<0){
    	if(i==0)
    	continue;
       flag=true;	
    }
   // Writting the profit earned and last day of selling    
    if (flag==true && sum!=0){
    outfile<<counter<<'\t'<<'\t'<<'\t'<<'\t'<<'\t'<<sum<<'\n';
    pft=0;	
    sum=0; 
    flag=false;
    }
   
    // Summing the profit when it is positive
    else if (pft>0) {
 	if(sum==0){
	outfile<<i+1<<'\t'<<'\t'<<'\t'<<'\t'<<'\t'<<'\t';
    	  }
      sum=pft+sum;
      counter=i+1;
      flag=false;
          }
    
  }	
  outfile.close();
}



