#include<string>
#include<vector>
#include <gsl\gsl_vector.h>
#include <gsl\gsl_matrix.h>
#include <map>
#include "element.h"
#include "stdafx.h"
#include <iostream>
#include <sstream>
#include <regex>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <gsl\gsl_test.h>
#include <gsl\gsl_ieee_utils.h>


#ifndef __KARTHIK_simModel_CPP__
#define __KARHTIK_simModel_CPP__




//----------
//Model Data
//----------
class ModelData 
{
public:
	ModelData();	
	gsl_matrix * nodeCoordArray;
	int numNodes;
	int numElements;
	std::string elementType;
	gsl_matrix_int * elementConnectivityArray;		
	std::map<std::string, gsl_vector_int *> nodeSetMap;
	std::map<std::string, gsl_vector_int *> elementSetMap;	
    std::map<std::string, Element*> elementLibrary;
    Element *currElement;
	int iDimension;
	Material materialObj; //hard coded for now
	bool addElementSetToMap(std::string);
	bool addNodeSetToMap(std::string);
	bool addAmplitudesToMap(std::string);
	std::map<std::string, gsl_matrix *> amplitudeMap;	
	~ModelData();
};


//----------
//History Data
//----------

class LoadData
{
public:
	bool isNodeSetUsed;
	std::string nodeSet;
	int nodeId;
	int dof;
	double mag = 0.0;
};



class Load
{
public:
	Load();
	std::string loadType;
	bool isAmplitudeON;
	std::string amplitudeName;			
	std::vector<LoadData> loadDataVec;
	bool addLoadDataToSequence(std::string);
	~Load();
};








class BCData
{
public:
	bool isNodeSetUsed;	
	std::string nodeSet;
	int nodeId;
	int startDOF;
	int endDOF ;
	double mag = 0.0;
};


class BC
{
public:
	BC();
	std::string BCType;
	//std::string nodeSetName;
	bool isAmplitudeON;
	std::string amplitudeName;
	std::vector<BCData> bcData;
	bool addtobcData(std::string);

	~BC();
};


class OutputRequest
{
public:
	std::string OutputRequestName;
	std::string nodeSetName;
};

class Step
{
public:
	Step(std::string);
	std::string stepName;
	std::string stepType;
	std::string timeIncScheme;
	gsl_vector * timeParameters;
	double timeIncrement;
	std::vector<Load> loads;
	std::vector<BC> BCs;	
	bool addBCtoSequence(std::string);
	bool addLoadtoSequence(std::string);
	


	~Step();


};


class HistoryData
{

public:
	HistoryData();
	std::vector<Step> stepSequence;	
	//std::vector<OutputRequest> outputRequests;
	bool addSteps(std::string);
	~HistoryData();
};

class ResultData
{
public:
	ResultData();
	double time;	
	std::map<std::string, std::vector<double>> fields;
	bool setSVec(int);
	~ResultData();
};

//------------
//Simulation Model
//-----------
class SimModel
{
public:
	SimModel(std::string);
	std::string simModelName;
	ModelData modelData;
	HistoryData historyData;
	bool createSimModelFromInp(std::string strInpFileIn);
	int getNumNodes(const std::string& nodeBlockString);
	int getNumElements(const std::string& elementBlockString, int& numNodesPerElement);
	bool populateNodeArray(const std::string& nodeBlockString, gsl_matrix * nodeArray, int iDimension);
	bool populateELementArray(const std::string& elementBlockString, gsl_matrix_int * elementArray, const int numElements, const int numNodesPerElement);
	std::vector<ResultData> results;

	~SimModel();
};



#endif