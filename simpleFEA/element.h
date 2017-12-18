#pragma once
#include<string>
#include<vector>
#include <gsl\gsl_vector.h>
#include <gsl\gsl_matrix.h>
#include <map>
#include <vector>

//---------
//Root element class
//---------

class Element {
public:
	Element();
	int numNodes;
	std::string elementTypeName;
	int spaceDimensions;
	int numIntegrationPts;
	int numDOFS;
	int numDOFSperNode;
	int stateVarsPerIntPoint;
	bool get3DIntegrationPoints(const int n_points, const int n_nodes, gsl_matrix * xi, gsl_vector * w);
	bool get3DShapeFunctions(const gsl_matrix * xi, const int n_nodes, gsl_vector * f, gsl_matrix * df);
	bool get2DIntegrationPoints(const int n_points, const int n_nodes, gsl_matrix **xi, gsl_vector **w);
	bool get2DShapeFunctions(const std::vector<double> xi, const int n_nodes, gsl_vector **f, gsl_matrix **df);
	virtual bool getElementStiffnessAndResidualForce(gsl_vector * RHS, gsl_matrix * Kelm, double Time, double Dtime, gsl_matrix * Coords, 
		gsl_vector * U, int NNODE, std::map<std::string, std::vector<double>> &fields,int sVarsStart) = 0;
	~Element();
		
};


//------------------
//Create element library
//-------------------
class Element_CPE8 : public Element
{
public:
	Element_CPE8();	
	bool getElementStiffnessAndResidualForce(gsl_vector * RHS, gsl_matrix * Kelm,double Time,double Dtime,
		gsl_matrix * Coords, gsl_vector * U,int NNODE, std::map<std::string, std::vector<double>> &fields, int sVarsStart);


	~Element_CPE8();

};


//Qudratic Tri
class Element_CPE6 : public Element
{
public:
	Element_CPE6();
	bool getElementStiffnessAndResidualForce(gsl_vector * RHS, gsl_matrix * Kelm, double Time, double Dtime,
		gsl_matrix * Coords, gsl_vector * U, int NNODE, std::map<std::string, std::vector<double>> &fields, int sVarsStart);
	
	~Element_CPE6();
};




//Qudratic Tri
class Element_CPE3 : public Element
{
public:
	Element_CPE3();
	bool getElementStiffnessAndResidualForce(gsl_vector * RHS, gsl_matrix * Kelm, double Time, double Dtime,
		gsl_matrix * Coords, gsl_vector * U, int NNODE, std::map<std::string, std::vector<double>> &fields, int sVarsStart);

	~Element_CPE3();
};



//-----------------
//2nd level element class - ELement 2D,3D 
//Define common functions - like shape funcs etc
//-----------------
/*

class Element2D : public ElementBase {
public:
	bool get2DIntegrationPoints(const int n_points,const int n_nodes, gsl_matrix * xi, gsl_vector * w);
	bool get2DShapeFunctions(const gsl_matrix * xi, const int n_nodes, gsl_vector * f, gsl_matrix * df);
	void setElementTypeName(std::string);	
};

class Element3D : public ElementBase {
public:
	bool get3DIntegrationPoints(const int n_points, const int n_nodes, gsl_matrix * xi, gsl_vector * w);
	bool get3DShapeFunctions(const gsl_matrix * xi, const int n_nodes, gsl_vector * f, gsl_matrix * df);
	void setElementTypeName(std::string);
};

*/

class Material {
public:
	int numParameters;
};
