/*----------------------------------------

Major Assignment CFD 1
Member : Syahrir Ginanjar (23623005) 
		 Muhammad Farras Arira (23623006)
		 Bryan (23623011)	 
Case : Inviscid Burger solver (Steady condition, Explicit Scheme)

----------------------------------------*/
/*----------------------------------------
Notes:

CFL is Courant-Friedrich-Lewy number


----------------------------------------*/



#include<iostream>
#include<cmath>
#include<string>
#include<fstream>
#include "../../../OTODIDAK/SIM-Tool/Core/eigen-3.4.0/Eigen/Dense" //Sesuaikan aj rir utk bagian ni
using namespace std;
using namespace Eigen;

//Solver :
//#define Lax
//#define MacCormack
#define RungeKutta
//#define LW

//Geometric, simulation and flow properties
double L = 1., W = 1.; 
int Nx = 41, Ny = Nx;	
int N = Nx*Ny;
double dx = L/(Nx-1), dy = dx;

#if defined Lax
	double CFL = 0.1;
	double u_max = 1.5;
	double c = 340;
	double dt = sqrt( CFL/( pow(u_max+c,2)/pow(dx,2) + pow(c,2)/pow(dy,2)) );
	double tau_x = dt/dx, tau_y = dt/dy;
#elif defined MacCormack
	double CFL = 1.;
	double u_max = 1.5;
	double c = 340;
	double dt = CFL*dx*dy/(u_max*dy + c*sqrt(pow(dx,2) + pow(dy,2)));
	double tau_x = dt/dx, tau_y = dt/dy;
	
	#define WithAV
	#if defined WithAV
		double omega_x = 2e-3;
		double omega_y = 1e-9;
	#endif
#elif defined RungeKutta
	double CFL = 1;
	double u_max = 1.5;
	double c = 340;
	double dt = CFL*dx*dy/(u_max*dy + c*sqrt(pow(dx,2) + pow(dy,2)));
	double tau_x = dt/dx, tau_y = dt/dy;
	double epsilon_e = 0.01;
	
	int num_stage = 10;
#elif defined LW
	double CFL = 0.5;
	double u_max = 1.5;
	double c = 340;
	double dt = sqrt( CFL/( pow(u_max+c,2)/pow(dx,2) + pow(c,2)/pow(dy,2)) );
	double tau_x = dt/dx, tau_y = dt/dy;
#endif



double T = 1.8;
int Nt = T/dt;

//Print simulation info 
void Print_Info() {
	cout<<"---------------------------------------------------------------------------\n";
	cout<<"Case                                     : 2D Inviscid Burger\n";
	cout<<"Size of domain (L x W)                   : "<<L<<"x"<<W<<endl;
	cout<<"Number of nodes                          : "<<N<<endl;
	cout<<"Number of nodes x                        : "<<Nx<<endl;
	cout<<"Number of nodes y                        : "<<Ny<<endl;
	cout<<"Dirichlet BC                             : \n";
	cout<<"           -left                         : 1.5\n";
	cout<<"           -right                        : -0.5\n";
	cout<<"           -bottom                       : 1.5-2x \n";
	cout<<"CFL                                      : "<<CFL<<"\n";
	cout<<"Time step                                : "<<dt<<"\n";
	cout<<"Total time                               : "<<T<<"\n";
	cout<<"Total time step                          : "<<Nt<<"\n";
	
	#if defined MacCormack
		cout<<"Numerical Method                         : MacCormack\n";
	#elif defined RungeKutta
		cout<<"Numerical Method                         : Runge Kutta\n";
	#elif defined LW
		cout<<"Numerical Method                         : Lax Wendroff\n";
	#elif defined Lax
		cout<<"Numerical Method                         : Lax Method\n";
	#endif
	cout<<"---------------------------------------------------------------------------\n";
}

//Numerical Initialization for BC ()
MatrixXd Initialization() {
//	MatrixXd T_num_p = MatrixXd::Zero(Ny,Nx);
	MatrixXd T_num_p = MatrixXd::Constant(Ny,Nx, 1.473169786);
	
	//At left and right of boundary
	for(int i = 0; i<Ny; i++) {
		T_num_p(i,0) = 1.5;		//left
		T_num_p(i,Nx-1) = -0.5;	//right
	}
	
	
	//At top and bottom of boundary
	for(int j = 0; j<Nx; j++) {
		T_num_p(0,j) = 1.5;		//top
		T_num_p(Ny-1,j) = 1.5-2*j*dx;		//top
	}
	
	return T_num_p;
}

//Analytical solution : 
MatrixXd Analytic_Solution(MatrixXd& x, MatrixXd& y) {
	MatrixXd u_analytic = MatrixXd::Zero(Ny,Nx);
	
	for (int j = 0; j<Ny; j++) {
		for (int i= 0; i<Nx; i++) {
			if(y(j,i) <= 0.5) {
				if(x(j,i) <= 1.5*y(j,i)) {u_analytic(j,i) = 1.5;}
				else if (x(j,i)>=1-0.5*y(j,i)) { u_analytic(j,i) = -0.5;}
				else { u_analytic(j,i) = (1.5 - 2*x(j,i) )/( 1-2*y(j,i) ); }
			}
			else {
				if(x(j,i) <= 0.5+0.5*y(j,i)) {u_analytic(j,i) = 1.5;}
				else { u_analytic(j,i) = -0.5;}
			}
		}
	}
	
	return u_analytic;
}

//Print CSV result :
void OutputCSV_Error(int &nout, double *error, string name) {
        ofstream ofs;
		ofs.open(name);
		ofs << "time,error\n";
		for (int i = 0; i<nout; i++) {
			ofs << (i+1) * dt << "," << *(error+i)<<"\n";
		}
		ofs.close();
}
// For extracting data
void OutputCSV_Result(int &nx, int &ny, MatrixXd& result, string name) {
        ofstream ofs;
		ofs.open(name);
		ofs << "Position_x,Position_y,value\n";
		for (int j = 0; j<ny; j++) {
			for(int i = 0; i<nx; i++) {
				ofs << (i) * dx << "," << (Ny-j-1) * dy<< ","<< result(j,i)<<"\n";	
			}
			
		}
		ofs.close();
}


//Main function :
int main() {
	Print_Info();
	
	//Analytical and numerical initialization
	MatrixXd x = MatrixXd::Zero(Ny,Nx);	//physical format
	MatrixXd y = MatrixXd::Zero(Ny,Nx); //physical format
	for (int i = 0; i<Nx; i++) {
		for (int j = 0; j<Ny; j++) {
			x(j,i) = i*dx;
			y(j,i) = (Ny-j-1)*dy;
		}
	}

	
	//The format that are being used here is the physical format
	MatrixXd u_numeric = Initialization();	//For presentation purpose physical format
	
	//Analytical solution
	MatrixXd u_analytic = Analytic_Solution(x,y);	//For presentation purpose physical format
	
	//Numerical solution
	double error[Nt];
	int stop_time = Nt;
	MatrixXd u_new = u_numeric;
	MatrixXd u_before = u_numeric;
	for (int t = 0; t<Nt; t++) {
		//Updating E
		MatrixXd E = pow(u_numeric.array(),2)/2;
		
		//Employing numerical method
		#if defined Lax
			for(int j = 1; j<Ny-1; j++) {
				for(int i = 1; i<Nx-1; i++) {
					u_new(j,i) = 0.25*(u_numeric(j,i+1) + u_numeric(j,i-1) + u_numeric(j+1,i) + u_numeric(j-1,i) ) - 0.5*tau_x*(E(j,i+1) - E(j,i-1)) - 0.5*tau_y*( u_numeric(j-1,i) - u_numeric(j+1,i) );
				}
			}
			
		#elif defined MacCormack
			MatrixXd u_star = u_numeric;
			MatrixXd E_star = E;
			
			//Predictor
			for (int j = 1; j<Ny-1; j++) {
				for (int i= 1; i<Nx-1; i++) {
					u_star(j,i) = u_numeric(j,i) - tau_x*(E(j,i+1) - E(j,i-1)) - tau_y*(u_numeric(j-1,i) - u_numeric(j+1,i));
					
					#if defined WithAV
						if(i != 1 and i != Nx-2 and j != 1 and j != Ny-2) {
							u_star(j,i) += -omega_x*(u_numeric(j,i+2) - 4*u_numeric(j,i+1) + 6*u_numeric(j,i) - 4*u_numeric(j,i-1) + u_numeric(j,i-2))/8;	
							u_star(j,i) += -omega_y*(u_numeric(j+2,i) - 4*u_numeric(j+1,i) + 6*u_numeric(j,i) - 4*u_numeric(j-1,i) + u_numeric(j-2,i))/8;	
						}
					#endif
					
					E_star(j,i) = pow(u_star(j,i),2)/2;
				}
			}
			
			
			//Corrector and updating value
			for(int j = 1; j<Ny-1; j++) {
				for(int i = 1; i<Nx-1; i++) {
					double u_bb = u_numeric(j,i) - tau_x*(E_star(j,i) - E_star(j,i-1)) - tau_y*(u_star(j,i) - u_star(j+1,i));
					#if defined WithAV
						if (i != 1 and i != Nx-2 and j != 1 and j != Ny-2) {
							u_bb += 0.5*omega_x*(u_star(j,i+2) - 4*u_star(j,i+1) + 6*u_star(j,i) - 4*u_star(j,i-1) + u_star(j,i-2))/8;		
							u_bb += 0.5*omega_y*(u_star(j+2,i) - 4*u_star(j+1,i) + 6*u_star(j,i) - 4*u_star(j-1,i) + u_star(j-2,i))/8;		
						}
					#endif
					
					u_new(j,i) = 0.5*(u_star(j,i) + u_bb);
				}
			}
			
		
		#elif defined RungeKutta
			//Step 1 : 
			MatrixXd u_transition = u_numeric;	//u1
			MatrixXd E_transition = E;
			MatrixXd R_transition = E;
			
			/*
				The inviscid burger equation in this problem is 
				du/dt + u du/dx + du/dy = 0 -> du/dt = - u du/dx - du/dy = R(u)
				
				In other word : 
				R(u) = - (dE/dx + du/dy),  E = u^2/2
			*/
			
			//Calculation process
			for(int stage = num_stage; stage>=1; stage--) {
				for(int j = 1; j<Ny-1; j++) {
					for(int i = 1; i<Nx-1; i++) {
						//Calcuate R
						R_transition(j,i) = - (E_transition(j,i+1) - E_transition(j,i-1)/(2*dx)) - (u_transition(j-1,i) - u_transition(j+1,i))/(2*dy);
						
						//Calculate u
						if(stage != 1) {
							u_transition(j,i) = u_numeric(j,i) + R_transition(j,i)*dt/stage;
							E_transition(j,i) = pow(u_transition(j,i),2)/2;	
						}
						else {
							u_new(j,i) = u_numeric(j,i) + R_transition(j,i)*dt;
						}
					}
				}	
			}
			
			//Augmenting the damping term : 
			for(int j = 1; j<Ny-1; j++) {
				for(int i = 1; i<Nx-1; i++) {
					//Damping x : 
					double Dx = 0;
					if(i != 1 and i != Nx-2) Dx = u_numeric(j,i-2) - 4*u_numeric(j,i-1) + 6*u_numeric(j,i) - 4*u_numeric(j,i+1) + u_numeric(j,i+2);
					 
					
					//Damping y : 
					double Dy = 0;
					if(j != 1 and j !=Ny-2) Dy = u_numeric(j-2,i) - 4*u_numeric(j-1,i) + 6*u_numeric(j,i) - 4*u_numeric(j+1,i) + u_numeric(j+2,i);
					
					//Augmenting to the solution : 
					u_new(j,i) += -epsilon_e*(Dx + Dy);
				}
			}
			
		#elif defined LW
			MatrixXd u_star = u_numeric;
			MatrixXd u_dstar = u_numeric;
			MatrixXd F_star = E;

			for(int j = 1; j<Ny-1; j++) {
				for(int i = 1; i<Nx-1; i++) {
					u_star(j,i) = u_numeric(j,i)-tau_x*(E(j,i+1)-E(j,i))-tau_y*(u_numeric(j+1,i)-u_numeric(j,i));
					F_star(j,i) = pow(u_star(j,i),2)/2;
				}
			}

			for(int j = 1; j<Ny-1; j++) {
				for(int i = 1; i<Nx-1; i++) {
					u_dstar(j,i) = u_numeric(j,i)-tau_x*(F_star(j,i)-F_star(j,i-1))-tau_y*(u_star(j,i)-u_star(j-1,i));
					u_new(j,i) = 0.5*(u_star(j,i)+u_dstar(j,i));
				}
			}

		#endif
		
		u_before = u_numeric;
		u_numeric = u_new;
		
		//Error calculation for all numerical methods
		double sum = 0;
		double sum_analytic = 0;
		for (int i = 0; i<Nx; i++) {
			for(int j = 0; j<Ny; j++) {
				sum += pow((u_numeric(j,i) - u_analytic(j,i)),2) ;
				sum_analytic += pow(u_analytic(j,i),2);	
			}
		}
		error[t] = sqrt(sum/sum_analytic)*100;
		cout<<"t = "<<t<<"\nError = "<<error[t]<<"\n\n";
		
		if(error[t]>error[t-1] and t != 0)  {
			u_numeric = u_before;
			stop_time = t;
			break;
		}
		
	}
	#if defined Lax
		string Method = "Lax";
	#elif defined MacCormack
		string Method = "MacCormack";
	#elif defined RungeKutta
		string Method = "RungeKutta";
	#elif defined LW
		string Method = "Lax_Wendroff";
	#endif
	
	string name_error = "Error for the " + Method + " method.csv";
	OutputCSV_Error(stop_time, error, name_error);
	
	string name_numeric = "Numerical result for " + Method + " method.csv";
	OutputCSV_Result(Nx, Ny, u_numeric, name_numeric);
	
	string name_analytic = "Analytical result for " + Method + " method.csv";
	OutputCSV_Result(Nx, Ny, u_analytic, name_analytic);
	

}