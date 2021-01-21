#include <time.h>
#include <ctime>
#include <sstream>
#include <math.h>
#include <stdio.h>
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

string DataSetType = "Cluster";
string RPath;

const int N = 100;  // Number of nodes for recharge stations
const int B = 100;  // Number of nodes for base stations
const int C = 600;  //Number of customers
const double rd = 1;
const double rc = 2;
int MinCoverage = 401;
const float C_dc = 4;
const float C_rs = 1;
int M = N + B;
bool HinderSameLocation = 0;
bool Plot = 0;
double Obj;
int TotalRS=0;
int TotalDC=0;
int TotalSatDemand;
int Tilim = 10800;				// Time Limit in terms of the seconds


IloRangeArray ConstraintSet(env);
IloExprArray expression(env);
IloExpr expr(env);
IloExpr MinimizeCost(env);
ofstream output("C:\\Users\\...\\1-Results.txt");
ofstream outFacilities("C:\\Users\\...\\Selected_Facilities.txt");
ofstream outCustomers("C:\\Users\\...\\Selected_Customers.txt");
ofstream outFlows("C:\\Users\\...\\Flows.txt");
ofstream outFromFlows("C:\\Users\\...\\FromFlows.txt");
ofstream outInFlows("C:\\Users\\...\\ToFlows.txt");

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

void ObjectiveFunction(IloModel model, IloIntVarArray x)
{
	for (int i = 0; i < N; i++)
	{
		expr += C_rs * x[i];
	}

	for (int i = N; i < M; i++)
	{
		expr += C_dc * x[i];
	}

	MinimizeCost += expr;
	expr.clear();
}

/////**** CONSTRAINT SETS ****/////
/*------------------------------------*/
/////*** DEMANDS CAN BE SATISFIED ONLY IF A FACILITY IS OPENED IN A SPECIFIC DISTANCE ****/////

void const_SAT_DEMAND(IloModel model, IloIntVarArray x, IloIntVarArray y, IloNumArray cd)
{
	//Demand can be satisfied from any reachable existing facility

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
			if (d[index(i, j)] <= rc & i != j)
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
	for (int k = N; k < M; k++) //candidate DC sayýsý kadar
	{
		expr = expr - M * x[k];

		for (int i = 0; i < M; i++)
		{
			if (d[index(i, k)] <= rc & k != i)
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
			if (d[index(i, k)] <= rc & i != k)
			{
				expr = expr + z[index(i, k)] - z[index(k, i)];
			}
		}

		ConstraintSet.add(expr >= 0);
		expr.clear();
	}
}

/////*** MIN COVERAGE CONSTRAINT****/////

void const_MIN_COV(IloModel model, IloIntVarArray y, IloIntArray a)
{
	for (int k = 0; k < C; k++)
	{
		expr += a[k] * y[k];

	}
	ConstraintSet.add(expr >= MinCoverage);
	expr.clear();
}

/////*** HINDER SAME LOCATION ****/////

void const_HIND_SAME_LOC(IloModel model, IloIntVarArray x)
{
	for (int i = 0; i < N; i++)
	{
		expr = x[i] + x[i + B];
		ConstraintSet.add(expr <= 1);
		expr.clear();
	}
}


void const5(IloModel model, IloIntVarArray z, IloIntVarArray x, IloNumArray d)
{  /////NOT USING/////
	for (int i = 0; i < N; i++)
	{
		expr += x[i];
	}
	for (int i = N; i < M; i++) //candidate DC sayýsý kadar
	{
		for (int j = 0; j < M; j++)
		{
			if (d[index(i, j)] <= rc & i != j)
			{
				expr += z[index(i, j)];
				expr -= z[index(j, i)];
			}
		}
	}
	ConstraintSet.add(expr <= 0);
	expr.clear();
}


/////**** GET SOLUTION VALUES ****/////
/*------------------------------------*/

void solution_values(IloCplex cplex, IloIntVarArray x, IloIntArray x_var, IloIntVarArray y, IloIntArray y_var,
	IloIntVarArray z, IloIntArray z_var, IloNumArray d, IloNum elapsed_time, IloNumArray cX, IloNumArray cY)
{
	env.out() << "Solution status = " << cplex.getStatus() << endl;
	env.out() << "Solution value = " << cplex.getObjValue() << endl;
	env.out() << "Number of Alternative Solutions = " << cplex.getSolnPoolNsolns() << endl;

	Obj = cplex.getObjValue();

	outFacilities << "Type\tCoordX\tCoordY\tCost\tCov\ttime\tRS\tDC" << endl;
	
	for (int i = 0; i < C; i++)
	{
		y_var[i] = cplex.getValue(y[i]);
		TotalSatDemand += y_var[i];
		output << "Y" << i + 1 << " " << "=" << " " << y_var[i] << endl;
		outCustomers << y_var[i] << endl;
	}

	for (int i = 0; i < M; i++)
	{
		x_var[i] = cplex.getValue(x[i]);

		if (x_var[i] == 1)
		{
			if (i < N)
			{
				TotalRS += x_var[i];

				outFacilities << "RS\t" << cX[i] << "\t" << cY[i] << "\t" << Obj << 
					"\t" << TotalSatDemand << "\t" << elapsed_time << "\t" << TotalRS << "\t" << TotalDC << endl;
			}
			if (i >= N)
			{
				TotalDC += x_var[i];

				outFacilities << "DC\t" << cX[i] << "\t" << cY[i] << "\t" << Obj << 
					"\t" << TotalSatDemand << "\t" << elapsed_time << "\t" << TotalRS << "\t" << TotalDC << endl;
			}
		}

		output << "X" << i + 1 << " " << "=" << " " << x_var[i] << endl;
		
		for (int j = 0; j < M; j++)
		{
			if (d[index(i, j)] <= rc & i != j)
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


}

/////**** MODELLING USING CPLEX ****/////
/*------------------------------------*/

void model(IloNumArray d, IloNumArray cd, IloIntArray a, IloNumArray cX, IloNumArray cY)
{

	string var_name;
	IloTimer elapsed_time(env);

	IloIntArray x_var(env, M);
	IloIntArray z_var(env, M * M);
	IloIntArray y_var(env, C);

	elapsed_time.start();	// Start Timer 

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

	IloModel model(env);	// Initiate Model Environment

	// Initiate Objective Function and Constraint Set //

	ObjectiveFunction(model, x);
	const_SAT_DEMAND(model, x, y, cd);
	const_RS_OUTFLOW(model, z, x, d);
	const_RS_FLOW_BALANCE(model, z, x, d);
	//const5(model, z, x, d);
	const_DC_INFLOW(model, z, x, d);
	const_MIN_COV(model, y, a);
	if (HinderSameLocation)
	{
		const_HIND_SAME_LOC(model, x);
	}

	expr.clear();

	model.add(ConstraintSet);
	IloObjective objective = IloMinimize(env, MinimizeCost);
	model.add(objective);
	IloCplex cplex(model);
	//cplex.setOut(env.getNullStream()); // turn off the cplex screen outputs
	cplex.exportModel("1-My_Model.lp");
	//  cplex.setParam(cplex.EpGap, 0.0008); // for testing purposes
	cplex.setParam(cplex.TiLim, Tilim);
	cplex.solve();


	//////TO GET NUMBER OF ALTERNATIVE SOLUTIONS/////
	/*cplex.setParam(IloCplex::Param::MIP::Pool::AbsGap, 0.0);
	cplex.setParam(IloCplex::Param::MIP::Pool::Intensity, 4);
	cplex.setParam(IloCplex::Param::MIP::Limits::Populate, 21000000);
	cplex.solve();
	cplex.populate();*/


	cout << "elapsed time: "<< elapsed_time.getTime() << endl;

	solution_values(cplex, x, x_var, y, y_var, z, z_var, d, elapsed_time.getTime(), cX, cY);
	elapsed_time.reset();
	ConstraintSet.end();
	MinimizeCost.end();
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

	//////////////////////////////////////////////////////////////////////////////
	// now define all parameters necessary for the model

	IloNumArray d(env, 1);       //distance btw i and j
	IloNumArray cd(env, 1);    //customer distances
	IloIntArray a(env, 1);  // demand amount
	IloNumArray cX(env, 1);
	IloNumArray cY(env, 1);

	d.setSize(M * M);
	cd.setSize(C * M);
	a.setSize(C);
	cX.setSize(M);
	cY.setSize(M);

	///////////////////////////////////////////////////////
	///ADD INPUT DATA!!!! ///////
	//now input data
	ifstream in;
	in.open(demand); //Get Demands
	for (int i = 0; i < C; i++)
	{

		in >> a[i];
	}
	in.close();

	in.open(distance); //Get Distances

	for (int j = 0; j < M; j++)
	{
		for (int i = 0; i < M; i++)
		{
			in >> d[index(j, i)];
		}
	}
	in.close();

	in.open(cust_dist); //Get Customer Distances

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


	model(d, cd, a, cX, cY);		// Initiate Model


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
	
	cout << "Total RS: " << TotalRS << endl <<
		"Total DC: "<<TotalDC<<endl<<
		"Total Satisfied Demand: "<< TotalSatDemand<<endl<<
		"Total Reachable Demand: " << TotalReachableDemand << endl << "Total Demand: " << TotalDemand << endl;
	env.end();

	if (Plot)
	{
		system(RPath.c_str());
	}

	cin >> STOP;
	return 0;
}