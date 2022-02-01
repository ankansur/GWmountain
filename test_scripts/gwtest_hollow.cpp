//========================================================================================================================================
//
//   Description: A C++ program to solve the Grad-Shafranov equation to find equilibria models
//                                  
//   
//   Version: 3.0 (working version)
//   Date: April 2021
//  
//   Author: Ankan Sur
//   Affiliation: Nicolaus Copernicus Astronomical Center, Warsaw, Poland
//
//========================================================================================================================================


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
    double u[Nu], x[Nr],r[Nr];
    double rc = 0.04*100.0; //(convert m to cm)
    double rout = 100*rc; //convert into cm
    double Rs = 1e6;
    double R = 20*rc;
    double Mdot = 6e-13;
    double Bs = 1e8;
    double rA = 3.53*1000*1e5*pow(Bs/1e12,4.0/7.0)*pow(Mdot/1e-9,-2.0/7.0);
    double M = 1.4*1.989e33;
    double G = 6.67e-8;
    double cs = 1e8;
    double x0 = R;
    double a = Rs/R;
    double u_in = 1.0;
    double u_out = sqrt(1-Rs/rA);
    double du = (u_out-u_in)/(Nu-1);
    double dr1 = rc/(Nr/2.0);
    double dr2 = (rout-rc)/(Nr/2.0 -1);
    double dr = rout/(Nr-1);
    int i,j;
    double A[Nr][Nu],Ac[Nr][Nu],S[Nr][Nu],rho[Nr][Nu],Q[Nr][Nu];
    double e = 100.0;
    double den,Qp,Qc,Source,w,max1,min1,counter,threshold,dAdr,dAdu,ncounter,diffr;
    int rid,uid;
    ofstream outputfile1,outputfile2,outputfile3,outputfile4;
    double Pi = atan(1)*4;
    double g = G*M/Rs/Rs;
    double gamma=5.0/3.0;
    double Kappa = 5.4e9;
    double f,r0,t1,t2,sinthp,PsiA;
    double Rp = pow(Rs/rA,0.5)*Rs;
    double dx;
    
    //required input parameters
    cout<<"enter threshold =";
    cin>>threshold ;

    //extra parameters
    counter=0;

    sinthp = sin(Rp/Rs);
    
    //under-relaxation parameter
    w = 0.1;
    
    // initialize grid
    for (int j=0; j<Nu; j++){ 
        u[j] = u_in + j*du;
        }
    
    //cout<<acos(u[0])*180/Pi;
    
    for (int i=0; i<Nr; i++){ 
        if (i<=Nr/2.0){
            r[i] = (i*dr1+Rs);
            }
        else{
            r[i] = (rc + i*dr2+Rs);
            }
        x[i] = (r[i]-Rs)/R;
        }

    /*for (int i=0; i<Nr; i++){ 
        r[i] = (i*dr+Rs);
        x[i] = (r[i]-Rs)/R;
        }*/

    
    outputfile4.open("Source.txt");
    /*
    outputfile4<<"Domain above star = "<<rout/100.0<<" m\n";
    outputfile4<<"printing colatitude range = "<<acos(u[0])*180/Pi<<" to "<<acos(u[Nu-1])*180/Pi<<" degrees"<<"\n";
    outputfile4<<"printing radial range = "<<r[0]<<" to "<<r[Nr-1]<<" cm"<<"\n";
    outputfile4<<"printing x range = "<<x[0]<<" to "<<x[Nr-1]<<"\n"; 
    outputfile4<<"printing a = "<<a<<"\n";
    outputfile4<<"printing x0 = "<<x0<<" cm\n";
    outputfile4<<"Alfven radius = "<<rA/1e5<<" km \n";
    outputfile4<<"Polar cap radius = "<<Rp/1e5<<" km \n";
    outputfile4<<"Threshold height = "<<54*pow(Bs/1e12,0.41)*pow(Rp/1e5,0.42)<<" m \n";
    */
    outputfile4<<"Rout\tColatitude_min\tColatitude_max\trc\tRad_min\tRad_max\tAlfven_radius\tPolar_cap_rad\tTh_height\n";
    outputfile4<<rout<<"\t"<<acos(u[0])*180/Pi<<"\t"<<acos(u[Nu-1])*180/Pi<<"\t"<<rc<<"\t"<<r[0]<<"\t"<<r[Nr-1]<<"\t"<<rA<<"\t"<<Rp<<"\t"<<54*pow(Bs/1e12,0.41)*pow(Rp/1e5,0.42)<<"\n";
    outputfile4.close();


    // initiliaze A
    for (i = 0; i < Nr; i++) {
        for (j=0; j < Nu; j++){
            A[i][j] = pow(Rs,3)*(1-pow(u[j],2))/Rp/Rp/r[i];
       }
   }
   
   // set boundaries
    for (i=0;i<Nr;i++){
         A[i][0] = 0.0;
         A[i][Nu-1] = Rs*Rs*Rs*(1-pow(u[Nu-1],2))/Rp/Rp/r[i]; 
         }
    for (j=0;j<Nu;j++){
        A[0][j] =  Rs*Rs*(1-pow(u[j],2))/Rp/Rp; 
        A[Nr-1][j] = Rs*Rs*Rs*(1-pow(u[j],2))/Rp/Rp/r[Nr-1]; 
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

    while (e>threshold){


    std::copy(&A[0][0], &A[0][0]+Nr*Nu,&Ac[0][0]);

    for (i=1;i<Nr-1;i++){
        dx = x[i+1]-x[i];
        for (j=1;j<Nu-1;j++){
            diffr = Rs + rc*(1.0-4.0*A[i][j]*A[i][j])-r[i];
            if (diffr>0.0){
                rho[i][j] = pow(g*2.0/5.0/Kappa,1.5)*pow(diffr,1.5);
                Q[i][j] = 64.0*Pi*pow(R,4)*pow(x[i]+a,2)*(1.0-u[j]*u[j])*rho[i][j]*g*rc/pow(Bs*Rp*Rp,2);
                }
            else{ 
                rho[i][j] = 0.0;
                Q[i][j] = 0.0;
                }
            den = (2.0/dx/dx + 2.0*(1-u[j]*u[j])/pow(x[i]+a,2.0)/du/du + Q[i][j]);
            A[i][j] = (1-w)*Ac[i][j] + w*((A[i+1][j]+A[i-1][j])/dx/dx + (1-u[j]*u[j])*(A[i][j+1]+A[i][j-1])/pow(x[i]+a,2)/du/du)/den;
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

    if (sqrt(e)>10000.0 || isnan(sqrt(e))){
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
               
                outputfile2<<(A[count][index]+0.5)<<" ";  
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
   cout<<"time elapsed = "<<elapsed_time_ms/1000.0<<"\n";
   
    return 0;
}
