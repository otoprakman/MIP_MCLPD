#include <time.h>
#include <ctime>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include <chrono> 
#include <ilcplex/ilocplex.h>

using namespace std;

//define global constants
#define infinity 10000000000
#define MODLUS 2147483647
#define MULT1      24112
#define MULT2      26143
ILOSTLBEGIN
IloEnv env;

///////////////////////////////////////////////////////
//////*** DEFINE PARAMETERS FOR THE MODEL ***//////////
///////////////////////////////////////////////////////

string DataSetType = "Cluster";  // Type of the DataSet

const int N = 100;				// Number of nodes for recharge stations
const int B = 100;				// Number of nodes for base stations
const int C = 600;				// Number of customers
const int rd = 1;				// Coverage Distance
const int rc = 2;				// Connectivity Distance
const int M = N + B;
const float c_dc = 4.0;
const float c_rs = 1.0;
const float total_cost = 14.0;
int P = 5;						// Number of recharge stations to locate
int F = 1;						// Number of facilities to locate
bool Plot = 0;					// Plot at R environment
double Obj;						// Assign Objective Value to a variable
string RPath;					// Path of R.exe
int Tilim= 10800;				// Time Limit in terms of the seconds

IloRangeArray ConstraintSet(env);
IloExprArray expression(env);
IloExpr expr(env);
IloExpr MaxSatisfiedDemand(env);

/////**** OUTPUT FILES ****/////

ofstream output("C:\\Users\\...\\1-Results.txt");
ofstream outFacilities("\\Users\\...\\Selected_Facilities.txt");
ofstream outCustomers("\\Users\\...\\Selected_Customers.txt");
ofstream outFlows("\\Users\\...\\Flows.txt");
ofstream outFromFlows("\\Users\\...\\FromFlows.txt");
ofstream outInFlows("\\Users\\...\\ToFlows.txt");

/////**** TO CREATE VARIABLE NAMES ****/////

string getvarname_x(int i) //** VARIABLE X: FACILITY LOCATIONS **//
{
	char ch[40];

	sprintf_s(ch, "x%d", i + 1);
	return(ch);
}
string getvarname_y(int w) //** VARIABLE Y: DEMAND LOCATIONS **//
{
	char ch[40];

	sprintf_s(ch, "y%d", w + 1);
	return(ch);
}
string getvarname_z(int w, int i) //** VARIABLE Z: FLOW btw FACILITIES **//
{
	char ch[40];

	sprintf_s(ch, "z%d,%d", w + 1, i + 1);
	return(ch);
}
/////**** FUNCTIONS FOR DEFINING DOUBLE INDICES VARIABLES ****//////

int index(int i, int j)
{
	return(i * M + j);
}


/////**** OBJECTIVE FUNCTION****/////

void ObjectiveFunction(IloModel model, IloIntVarArray y, IloIntArray a, IloIntVarArray x)
{
	for (int k = 0; k < C; k++)
	{
		expr += a[k] * y[k];

	}

	MaxSatisfiedDemand += expr;
	expr.clear();
}


/////**** CONSTRAINT SETS ****/////
/*------------------------------------*/
/////*** LIMIT NUMBER OF RS TO P ****/////

void const_NUM_RS(IloModel model, IloIntVarArray x) 
{
	for (int i = 0; i < N; i++)
	{
		expr += x[i];
	}
	ConstraintSet.add(expr <= P);
	expr.clear();
}

/////*** LIMIT NUMBER OF DC TO F ****/////

void const_NUM_DC(IloModel model, IloIntVarArray x) 
{
	for (int i = N; i < M; i++)
	{
		expr += x[i];
	}
	ConstraintSet.add(expr <= F);
	expr.clear();
}

void const_COST(IloModel model, IloIntVarArray x) {
	for (int i = 0; i < N; i++)
	{
		expr += x[i] * c_rs;
	}

	for (int i = N; i < M; i++)
	{
		expr += x[i] * c_dc;
	}
	ConstraintSet.add(expr == total_cost);
	expr.clear();
}

/////*** DEMANDS CAN BE SATISFIED ONLY IF A FACILITY IS OPENED IN A SPECIFIC DISTANCE ****/////

void const_SAT_DEMAND(IloModel model, IloIntVarArray x, IloIntVarArray y, IloNumArray cd)
{
	for (int k = 0; k < C; k++) 
	{
		expr = expr - y[k];
		for (int j = 0; j < M; j++)
		{

			if (cd[index(k, j)] <= rd)
			{		
				expr = expr + x[j];
			}
		}
		ConstraintSet.add(expr >= 0);
		expr.clear();
	}
}

/////*** THERE SHOULD BE OUTFLOW IF A RS IS LOCATED ****/////

void const_RS_OUTFLOW(IloModel model, IloIntVarArray z, IloIntVarArray x, IloNumArray d)
{
	for (int i = 0; i < N; i++) 
	{
		expr = expr - M * x[i];

		for (int j = 0; j < M; j++)
		{
			if (i != j & d[index(i, j)] <= rc )
			{
				expr = expr + z[index(i, j)];
			}
		}
		ConstraintSet.add(expr <= 0);
		expr.clear();
	}
}

/////*** THERE SHOULD BE INFLOW IF A DC IS LOCATED ****/////

void const_DC_INFLOW(IloModel model, IloIntVarArray z, IloIntVarArray x, IloNumArray d)
{
	for (int k = N; k < M; k++) 
	{
		expr = expr - M * x[k];

		for (int i = 0; i < M; i++)
		{
			if (k != i & d[index(i, k)] <= rc)
			{
				expr = expr + z[index(i, k)];
			}
		}
		ConstraintSet.add(expr <= 0);
		expr.clear();
	}
}

/////*** FLOW BALANCE OF EACH LOCATED RS ****/////

void const_RS_FLOW_BALANCE(IloModel model, IloIntVarArray z, IloIntVarArray x, IloNumArray d)
{
	for (int i = 0; i < N; i++)
	{
		expr = expr - x[i];

		for (int k = 0; k < M; k++)
		{
			if (i != k & d[index(i, k)] <= rc)
			{
				expr = expr + z[index(i, k)] - z[index(k, i)];
			}
		}

		ConstraintSet.add(expr >= 0);
		expr.clear();
	}
}


/////**** GET SOLUTION VALUES ****/////
/*------------------------------------*/

void solution_values(IloCplex cplex, IloIntVarArray x, IloIntVarArray y, IloIntVarArray z, 
	IloNumArray d, IloNum elapsed_time, IloNumArray cX, IloNumArray cY)
{
	int x_var[M],
		y_var[C],
		z_var[M*M];

	env.out() << "Solution status = " << cplex.getStatus() << endl;
	env.out() << "Solution value = " << cplex.getObjValue() << endl;
	env.out() << "Number of Alternative Solutions = " << cplex.getSolnPoolNsolns() << endl;
	Obj = cplex.getObjValue();

	outFacilities << "Type\tCoordX\tCoordY\tCost\tCov\ttime\tRS\tDC" << endl;

	for (int i = 0; i < M; i++)
	{
		x_var[i] = cplex.getValue(x[i]);

		output << "X" << i + 1 << " " << "=" << " " << x_var[i] << endl;

		if (x_var[i] == 1)
		{
			if (i < N);
			{
				outFacilities << "RS\t" << cX[i] << "\t" << cY[i] << "\t" << P * 1 + F * 4 << "\t" << Obj << "\t" << elapsed_time << 
					"\t" << P << "\t" << F << endl;
			}
			if (i>=N)
			{
				outFacilities << "DC\t" << cX[i] << "\t" << cY[i] << "\t" << P * 1 + F * 4 << "\t" << Obj << "\t" << elapsed_time <<
					"\t" << P << "\t" << F << endl;
			}
		}

		for (int j = 0; j < M; j++)
		{
			if (i != j & d[index(i, j)] <= rc)
			{
				z_var[index(i, j)] = cplex.getValue(z[index(i, j)]);
				output << "Z" << i + 1 << "," << j + 1 << " " << "=" << " " << z_var[index(i, j)] << endl;
				if (z_var[index(i, j)] > 0)
				{
					outFromFlows << i + 1 << endl;
					outInFlows << j + 1 << endl;
					outFlows << z_var[index(i, j)] << endl;
				}

			}
		}
	}

	for (int i = 0; i < C; i++)
	{
		y_var[i] = cplex.getValue(y[i]);
		output << "Y" << i + 1 << " " << "=" << " " << y_var[i] << endl;
		outCustomers << y_var[i] << endl;
	}
}

/////**** MODELLING USING CPLEX ****/////
/*------------------------------------*/

void model(IloNumArray d, IloNumArray cd, IloIntArray a, IloNumArray cX, IloNumArray cY)
{

	string var_name;
	auto start = std::chrono::high_resolution_clock::now();

	
	// Initiate Variables //

	IloIntVarArray x(env, M, 0, 1);
	IloIntVarArray y(env, C, 0, 1);
	IloIntVarArray z(env, M * M, 0, 999);

	// Naming Variables //

	for (int w = 0; w < M; w++)
	{
		var_name = getvarname_x(w);
		x[w].setName(var_name.c_str());
		for (int i = 0; i < M; i++)
		{
			var_name = getvarname_z(w, i);
			z[index(w, i)].setName(var_name.c_str());
		}
	}

	for (int r = 0; r < C; r++)
	{
		var_name = getvarname_y(r);
		y[r].setName(var_name.c_str());
	}

	IloModel model(env);   // Initiate Model Environment

	// Initiate Objective Function and Constraint Set //

	ObjectiveFunction(model, y, a, x);
	/*const_NUM_RS(model, x);
	const_NUM_DC(model, x);*/
	const_COST(model, x);
	const_SAT_DEMAND(model, x, y, cd);
	const_RS_OUTFLOW(model, z, x, d);
	const_DC_INFLOW(model, z, x, d);
	const_RS_FLOW_BALANCE(model, z, x, d);

	expr.clear();

	model.add(ConstraintSet);
	IloObjective objective = IloMaximize(env, MaxSatisfiedDemand);
	model.add(objective);
	IloCplex cplex(model);
	//cplex.setOut(env.getNullStream());			// Turn off the cplex screen outputs
	cplex.exportModel("1-My_Model.lp");
	//  cplex.setParam(cplex.EpGap, 0.0008);		// for testing purposes
	cplex.setParam(cplex.TiLim, Tilim);				// Time Limit
	cplex.solve();


	//////TO GET NUMBER OF ALTERNATIVE SOLUTIONS/////
	/*cplex.setParam(IloCplex::Param::MIP::Pool::AbsGap, 0.0);
	cplex.setParam(IloCplex::Param::MIP::Pool::Intensity, 4);
	cplex.setParam(IloCplex::Param::MIP::Limits::Populate, 21000000);
	cplex.solve();
	cplex.populate();*/

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	float elapsed_time = (float)duration.count() / 1000000.0;

	cout << "elapsed time: "<<elapsed_time<< endl;
	solution_values(cplex, x, y, z, d, elapsed_time, cX, cY);

	ConstraintSet.end();
	MaxSatisfiedDemand.end();
	cplex.end();
	objective.end();
	model.end();
}


int main(int, char**)
{
	int STOP;
	string distance,
		cust_dist,
		demand,
		coordX,
		coordY;

	distance = "C:\\Users\\...\\distance.txt";
	cust_dist = "C:\\Users\\...\\cust_dist.txt";
	demand = "C:\\Users\\...\\demand.txt";
	coordX = "C:\\Users\\...\\coordx(Facility).txt";
	coordY = "C:\\Users\\...\\coordy(Facility).txt";

	IloNumArray d(env, 1);       // Distances Between Facilities
	IloNumArray cd(env, 1);		 // Distances btw Customers and Facilities
	IloIntArray a(env, 1);		 // Demand amount
	IloNumArray cX(env, 1);
	IloNumArray cY(env, 1);

	d.setSize(M * M);
	cd.setSize(C * M);
	a.setSize(C);
	cX.setSize(M);
	cY.setSize(M);

	/////**** GET INPUT VALUES FROM FILES ****/////
	/*-------------------------------------------*/

	ifstream in;
	in.open(demand);				//Get Demands
	for (int i = 0; i < C; i++)
	{
		in >> a[i];
	}
	in.close();

	in.open(distance);				//Get Distances

	for (int j = 0; j < M; j++)
	{
		for (int i = 0; i < M; i++)
		{
			in >> d[index(j, i)];
		}
	}
	in.close();

	in.open(cust_dist);				//Get Customer Distances

	for (int j = 0; j < C; j++)
	{
		for (int i = 0; i < M; i++)
		{
			in >> cd[index(j, i)];			
		}
	}
	
	in.close();

	in.open(coordX);				

	for (int i = 0; i < M; i++)
	{
			in >> cX[i];
	}
	in.close();

	in.open(coordY);				

	for (int i = 0; i < M; i++)
	{
		in >> cY[i];
	}
	in.close();

	model(d, cd, a, cX, cY);				// Initiate Model

	int counter;
	int TotalReachableDemand = 0;
	int TotalDemand = 0;
	for (int i = 0; i < C; i++)
	{
		counter = 0;
		for (int j = 0; j < M; j++)
		{
			if (cd[index(i, j)] <= rd)
			{
				counter += 1;
			}
		}
		TotalDemand += a[i];
		if (counter == 0)
		{
			a[i] = 0;
		}
		TotalReachableDemand += a[i];
	}
	float Coverage;
	Coverage = Obj / TotalReachableDemand;
	cout << "Total Coverage: %" << Coverage * 100 << endl <<
		"Total Reachable Demand: " << TotalReachableDemand << endl << "Total Demand: " << TotalDemand << endl;
	env.end();

	if (Plot)
	{
		system(RPath.c_str());
	}

	cin >> STOP;
	return 0;
}