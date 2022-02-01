//====================================================================================
//
//   Description: A C++ program to solve the Grad-Shafranov equation to find equilibria models
//                                  
//   
//   Version: 1.0
//   Date: Nov 2020
//  
//
//   Equation: GS A = -r**2*sin(theta)*(electron denisty) + Beta*d(Beta)/dA
//               B = strength * (A - A(R,theta))**(power)
//              
//   The code requires you to put an input strength, power, and threshold error
//          
//   Author: Ankan Sur
//   Affiliation: Nicolaus Copernicus Astronomical Center, Warsaw, Poland
//
//====================================================================================



#include <iostream>
#include <math.h>
#include <typeinfo>
#include <fstream>
#include <omp.h>
#include <chrono>

using namespace std;


int main() {


    int Nr = 128;
    int Nu = 128;
    double u[Nu], r[Nr];
    double rc = 0.1*100.0;
    double rout = 0.001*1e5;
    double Rs = 1e6;
    double R = 1e5;
    double Mdot = 6e-13;
    double Bs = 1e8;
    double rA = 3.53*1000*1e5*pow(Bs/1e12,4.0/7.0)*pow(Mdot/1e-9,-2.0/7.0);
    double u_in = 1.0;
    double u_out = sqrt(1-Rs/rA);
    double du = (u_in-u_out)/(Nu-1);
    double dr = rout/(Nr-1);
    int i,j;
    double A[Nr][Nu], Ac[Nr][Nu], S[Nr][Nu],rho[Nr][Nu], Q[Nr][Nu];
    int error = 1.0;
    double den,Qp,Qc,Source,w,max1,min1,counter,threshold,dAdr,dAdu;
    int rid,uid;
    ofstream outputfile1,outputfile2,outputfile3;
    double Pi = atan(1)*4;
    double g = 1.86e14;
    double gamma=5.0/3.0;
    double Kappa = 5.4e9;
    double f,r0,t1,t2;
    double Rp = pow(Rs/rA,0.5)*Rs;
    
    //required input parameters
    threshold = 1e-4;

    //extra parameters
    counter=0;

    
    //under-relaxation parameter
    w = 0.01;
    
    // initialize grid
    for (j=0; j<Nu; j++){ 
        u[j] = u_out + j*du;
        }
    
    //cout<<acos(u[0])*180/Pi;
    
    for (i=0; i<Nr; i++){ 
        r[i] = (Rs + i*dr)/R;
        }
        
    cout<<"Domain above star = "<<rout/100<<" m\n";
    cout<<"printing colatitude range = "<<acos(u[Nu-1])*180/Pi<<" to "<<acos(u[0])*180/Pi<<" degrees"<<"\n";
    cout<<"printing radial range = "<<r[0]<<" to "<<r[Nr-1]<<" km"<<"\n";
    cout<<"Alfven radius = "<<rA/1e5<<" km \n";
    cout<<"Polar cap radius = "<<Rp/1e5<<" km \n";
    cout<<"Threshold height = "<<54*pow(Bs/1e12,0.41)*pow(Rp/1e5,0.42)<<" m \n";

    // initialize source
    for (i=0; i<Nr; i++){ 
        for(j=0;j<Nu;j++){
                S[i][j] = 0.0;
        }
    }
    
    // initiliaze A
    for (i = 0; i < Nr; i++) {
        for (j=0; j < Nu; j++){
            A[i][j] = (1-pow(u[j],2))/r[i];
       }
   }
   
   // set boundaries
    for (i=0;i<Nr;i++){
         A[i][0] = 0.0;
         A[i][Nu-1] = Rs*(1-pow(u[Nu-1],2))/r[i]/R;
         }
    for (j=0;j<Nu;j++){
        A[0][j] = Rs*(1-pow(u[j],2))/r[0]/R;
        A[Nr-1][j] = Rs*(1-pow(u[j],2))/r[Nr-1]/R;
    }
    
    outputfile3.open("A0.txt");

    for (int count = 0; count < Nr; count ++)
        {
            for (int index= 0; index < Nu; index++)
               
                outputfile3<<A[count][index]<<" ";  
            outputfile3<<endl;                          
        }
    outputfile3.close(); 
 
    auto t_start = std::chrono::high_resolution_clock::now();

    while (counter<10000000.0){
        
    std::copy(&A[0][0], &A[0][0]+Nr*Nu,&Ac[0][0]);


    //#pragma omp parallel for collapse(2)
    for (i=1;i<Nr-1;i++){
        for (j=1;j<Nu-1;j++){
            r0 = Rs + rc*(1-pow(A[i][j]*rA/Rs,2));
            //cout<<r0-r[i]*R<<"\t";
            //rho[i][j] = pow(g*(gamma-1.0)/(gamma*Kappa),1.0/(gamma-1.0));//*pow(r0-r[i]*R,1.0/(gamma-1));
            if (r0-r[i]*R>0){
                //cout<<"entered";
                rho[i][j] = pow(g*2.0/5.0/Kappa,1.5)*pow(r0-r[i]*R,1.5);
                }
            else{ 
                rho[i][j] = 0.0;
                }
            Q[i][j] = -16.0*Pi*pow(R/Rs,4)*(rho[i][j]*g*rc/Bs/Bs)*(pow(r[i]*rA/Rs,2))*Ac[i][j];
            //dAdr = (Ac[i+1][j]-Ac[i-1][j])/2/dr;
            //dAdu = (Ac[i][j+1]-Ac[i][j-1])/2/du;
            //Qp = (Q[i+1][j]-Q[i-1][j])/dr/dAdr/2.0 + (Q[i][j+1]-Q[i][j-1])/du/dAdu;
            Qp = (Q[i+1][j] - Q[i-1][j])/(Ac[i+1][j]-Ac[i-1][j]) - (Q[i][j+1]-Q[i][j-1])/(Ac[i][j+1]-Ac[i][j-1]);
            Qc = Q[i][j] - Qp*Ac[i][j];
            den = (2.0/dr/dr + 2.0*(1-u[j]*u[j])/r[i]/r[i]/du/du) - std::min(0.0,Qp);
            Source = Q[i][j] - Qp*Ac[i][j] + std::max(0.0,Qp)*Ac[i][j];
            A[i][j] = (1-w)*Ac[i][j] + w*((A[i+1][j]+A[i-1][j])/dr/dr + (1-u[j]*u[j])*(A[i][j+1]+A[i][j-1])/r[i]/r[i]/du/du + Source)/den;
       }
    }

    double e = 0;
    for (int i=1; i< Nr-1; i++) {
      for (int j=1; j< Nu-1; j++){
         //cout<<e;
         e += pow((A[i][j]-Ac[i][j]),2)/pow(Ac[i][j],2);
      }
    }

    cout<<"error= "<<sqrt(e)<<" at step= "<<counter<<"\n"; 
    
    if (sqrt(e)<threshold){
        std::cout<<"Convergence reached: exiting computations"<<"\n";
        break;
        }

    if (sqrt(e)>100.0 || isnan(sqrt(e))){
        std::cout<<"Exiting loop, solution diverged"<<"\n";
        break;
        }
          
    
    counter+=1;
    }


    outputfile2.open("A.txt");
    cout<<"Produced output file"<<"\n";

    for (int count = 0; count < Nr; count ++)
        {
            for (int index= 0; index < Nu; index++)
               
                outputfile2<<A[count][index]<<" ";  
            outputfile2<<endl;                          
        }
    outputfile2.close(); 
    
    outputfile1.open("rho.txt");
    
    for (int count = 0; count < Nr; count ++)
        {
            for (int index= 0; index < Nu; index++)
               
                outputfile1<<rho[count][index]<<" ";  
            outputfile1<<endl;                          
        }
    outputfile1.close(); 
   
   auto t_end = std::chrono::high_resolution_clock::now();
   double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
   cout<<"time elapsed = "<<elapsed_time_ms/1000.0;
   
    return 0;
}
