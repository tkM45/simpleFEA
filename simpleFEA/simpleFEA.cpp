
//
// simpleFEA.cpp : Defines the entry point for the console application.
//GSL compiled in VS2015 using - https://www.gnu.org/software/gsl/extras/native_win_builds.html
//

#include "stdafx.h"
#include <iostream>
#include <sstream>
#include <string>
#include <regex>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <gsl\gsl_math.h>
#include <gsl\gsl_vector.h>
#include <gsl\gsl_test.h>
#include <gsl\gsl_ieee_utils.h>
#include <gsl\gsl_sf_bessel.h>
#include <gsl\gsl_matrix.h>
#include <gsl\gsl_splinalg.h>
#include "SimModel.h"
#include "myGSLMod.h"



bool computeProblemSize(const SimModel, gsl_spmatrix * globalStiffnes, gsl_vector * RHS, gsl_vector * U);
bool createGlobalStiffness(SimModel &simModel,gsl_spmatrix * globalStiffnes, gsl_vector * globalRHS,gsl_vector * U,double TIME,double DTIME, bool isCutback);
bool applyBC(SimModel, gsl_spmatrix * globalStiffnes, gsl_vector * globalRHS, gsl_vector * U, double TIME, double DTIME,int stepNumber);
double getScaleFromAmp(gsl_matrix *, double d);
bool writeOutputToTxt(SimModel);

int main()

{
	std::string inpFileName = "PlateTrianglesLinear";
	//std::string inpFileName = "PlateTriangles";
	//std::string inpFileName = "singlePlate_Force";
	//std::string inpFileName = "plateCoarse";

	SimModel mySimModel(inpFileName);

	mySimModel.createSimModelFromInp(inpFileName+".inp");

	
	gsl_spmatrix * globalStiffness;
	gsl_vector * globalRHS;
	gsl_vector * U;
	


	int numDofPerNode = mySimModel.modelData.currElement->numDOFSperNode;
	int numNodesInModel = mySimModel.modelData.numNodes;
	int totalDOFs = numNodesInModel*numDofPerNode;
	int numElements = mySimModel.modelData.numElements;	
	int totalIntPts = mySimModel.modelData.currElement->numIntegrationPts;
	int stateVarsPerIntPoint = mySimModel.modelData.currElement->stateVarsPerIntPoint;

	globalStiffness = gsl_spmatrix_alloc(totalDOFs, totalDOFs);
	globalRHS = gsl_vector_calloc(totalDOFs);
	U = gsl_vector_calloc(totalDOFs);
	std::vector<double> UVector;


	double TIME=0.0, DTIME = 0.0, initTimeInc,totalStepTime,minAllowedTimeInc, maxAllowedTimeInc;
	bool isCutBack = false;
	

	int numSteps = mySimModel.historyData.stepSequence.size();
	for (int i = 0; i < numSteps; i++) 
	{
		//Process steps
		Step currStep = mySimModel.historyData.stepSequence[i];

		if (currStep.stepType == "Static")
		{
			
			initTimeInc = gsl_vector_get(currStep.timeParameters,0);
			totalStepTime = gsl_vector_get(currStep.timeParameters, 1); 
			minAllowedTimeInc = gsl_vector_get(currStep.timeParameters, 2); 
			maxAllowedTimeInc = gsl_vector_get(currStep.timeParameters, 3);

			TIME = 0.0;
			DTIME = initTimeInc;
			
			while (true)
			{

				//Add frame to model 
				UVector.clear();
				if (mySimModel.results.size() == 0)
				{
					ResultData currResult;
					currResult.time = TIME;
					currResult.setSVec(numElements*totalIntPts*stateVarsPerIntPoint);
					mySimModel.results.push_back(currResult);
				}
				
				gsl_spmatrix_set_zero(globalStiffness);
				gsl_vector_set_zero(globalRHS);
				//gsl_vector_set_zero(U);

				createGlobalStiffness(mySimModel, globalStiffness, globalRHS, U, TIME, DTIME, isCutBack);

				if (isCutBack)
				{					
					DTIME = 0.5*DTIME;
					if (DTIME<= minAllowedTimeInc)
					{
						std::cout << "Time increment cutback below minimum allowed value";
						std::cout << "Analysis failed";
						return 0;
					}
				}

				if (TIME > totalStepTime) break;


				
				applyBC(mySimModel, globalStiffness, globalRHS, U, TIME, DTIME, i);



				//solve 
				const double tol = 1.0e-6;  /* solution relative tolerance */
				const size_t max_iter = 1000; /* maximum iterations */
				const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
				gsl_splinalg_itersolve *work =
					gsl_splinalg_itersolve_alloc(T, totalDOFs, 0);
				size_t iter = 0;
				double residual;
				int status;
				do
				{
					status = gsl_splinalg_itersolve_iterate(globalStiffness, globalRHS, tol, U, work);
					/* print out residual norm ||A*u - f|| */
					residual = gsl_splinalg_itersolve_normr(work);
					fprintf(stderr, "iter %zu residual = %.12e\n", iter, residual);
					if (status == GSL_SUCCESS)
						fprintf(stderr, "Converged\n");
				} while (status == GSL_CONTINUE && ++iter < max_iter);

				//writeGSLMatrix(globalStiffness, "globalStiffnessBeforeSolve.csv");
				////writeGSLMatrix(globalRHS, "globalRHSBeforeSolve.csv");

				//Add U to results

				if (mySimModel.results.back().time < TIME)
				{
					ResultData currResult;
					currResult.time = TIME;
					currResult.setSVec(numElements*totalIntPts*stateVarsPerIntPoint);
					mySimModel.results.push_back(currResult);
				}


				for (int i = 0; i < U->size; i++)
				{
					UVector.push_back(gsl_vector_get(U, i));
				}
				mySimModel.results.back().fields["U"] = UVector;


				//Update time steps
				TIME = TIME + DTIME;
				

			}
		}
		
		////writeGSLMatrix(U, "U_results.txt");

	}

	//-----------------------------
	//Write the odb
	//-----------------------------	
	


	writeOutputToTxt(mySimModel);
	

	std::cout << "\n\nCreating odb";

	std::string commandStr = "abq2017 python outToODB.py " + inpFileName;
	const char *commandCstr = commandStr.c_str();
	std::system(commandCstr);

	//gsl_vector_free(v);
	std::cout << "\n\n\DONE";
    return 0;
}


bool computeProblemSize(const SimModel simModel, gsl_spmatrix * globalStiffnes, gsl_vector * RHS, gsl_vector * U)
{
	int numNodes = simModel.modelData.numNodes;
	int numDOFsPerNode = simModel.modelData.currElement->numDOFSperNode;
	int totalDOFs = numNodes*numDOFsPerNode;
	int totalIntPts = simModel.modelData.currElement->numIntegrationPts;
	int numElements = simModel.modelData.numElements;

	
	globalStiffnes = gsl_spmatrix_alloc(totalDOFs, totalDOFs);
	RHS = gsl_vector_calloc(totalDOFs);
	U = gsl_vector_calloc(totalDOFs);	
	
	return true;
}



bool createGlobalStiffness(SimModel &simModel, gsl_spmatrix * globalStiffness, gsl_vector * globalRHS, gsl_vector * U, double TIME,double DTIME,bool isCutback)
{
	int numElements = simModel.modelData.numElements;
	int numDOFs = simModel.modelData.currElement->numDOFS;
	int numDOFperNode = simModel.modelData.currElement->numDOFSperNode;	
	int numDimension = simModel.modelData.currElement->spaceDimensions;
	int numNodesPerElement = simModel.modelData.currElement->numNodes;
	int numNodesInModel = simModel.modelData.numNodes;
	int totalIntPts = simModel.modelData.currElement->numIntegrationPts;
	int stateVarsPerIntPoint = simModel.modelData.currElement->stateVarsPerIntPoint;
	




	gsl_matrix * Kelm = gsl_matrix_calloc(numDOFs, numDOFs);
	gsl_vector * RHS_elem = gsl_vector_calloc(numDOFs);
	gsl_vector * U_elem = gsl_vector_calloc(numDOFs);
	gsl_matrix * Coords = gsl_matrix_calloc(numNodesPerElement, numDimension);

	std::string elementType=simModel.modelData.elementType ;


	const gsl_matrix_int * elemArray = simModel.modelData.elementConnectivityArray;
	const gsl_matrix * nodeArray = simModel.modelData.nodeCoordArray;


	////writeGSLMatrix(nodeArray, "NodeArray.txt");
	//writeGSLMatrix(elemArray, "ElemArray.txt");

	double currUx, currUy, nodeCoordX, nodeCoordY,tempData,Ux,Uy;
	int nodeId, dofNumber,globI,globJ,sVarsStart;
	

	for (int i = 0; i < numElements; i++)
	{	
		sVarsStart = i*totalIntPts*stateVarsPerIntPoint;


		//Create the Uelem and RHS_elem vectors
		for (int nodeIndex = 0; nodeIndex < numNodesPerElement; nodeIndex++)
		{
			nodeId = gsl_matrix_int_get(elemArray, i, nodeIndex);
			gsl_vector_set(U_elem, nodeIndex*2, gsl_vector_get(U,(nodeId-1)*2));			
			gsl_vector_set(U_elem, nodeIndex*2+1, gsl_vector_get(U,((nodeId-1)*2)+1));
		}
		
		//writeGSLMatrix(U_elem, "UElemAfterCreation.txt");
		//writeGSLMatrix(U, "UGlobalatUelem.txt");

		//Create the Coord Matrix for the element
		
		for (int j = 0; j < numNodesPerElement; j++)
		{			
			//currUx = gsl_vector_get(U, numDOFs*i + 2 * j );
			//currUy = gsl_vector_get(U, numDOFs*i + 2 * j +1);
			nodeCoordX = gsl_matrix_get(nodeArray, gsl_matrix_int_get(elemArray, i, j)-1, 0);
			nodeCoordY = gsl_matrix_get(nodeArray, gsl_matrix_int_get(elemArray, i, j)-1, 1);			
			gsl_matrix_set(Coords, j, 0, nodeCoordX);
			gsl_matrix_set(Coords, j, 1, nodeCoordY);
		}


		//std::cout << i<<std::endl;
		simModel.modelData.currElement->getElementStiffnessAndResidualForce(RHS_elem, Kelm, 0.0, .1, Coords, U_elem, numNodesPerElement, simModel.results.back().fields,sVarsStart);
		
		


		for (int i1 = 0; i1 < numNodesPerElement * numDOFperNode; i1++)
		{
			for (int j1 = 0; j1 < numNodesPerElement * numDOFperNode; j1++)
			{					

				
				globI = (gsl_matrix_int_get(elemArray, i,((i1 - i1%numDOFperNode) / numDOFperNode))-1)*numDOFperNode + i1%numDOFperNode;				
				globJ = (gsl_matrix_int_get(elemArray, i, ((j1 - j1%numDOFperNode) / numDOFperNode))-1)*numDOFperNode + j1%numDOFperNode;
				//std::cout << globI << " , " << globJ << " , "<<i1<<","<<j1<<std::endl;



				tempData = gsl_spmatrix_get(globalStiffness, globI, globJ);
				
				gsl_spmatrix_set(globalStiffness, globI, globJ, tempData+gsl_matrix_get(Kelm, i1, j1));
				
				
			}
		}
		

		//writeGSLMatrix(globalStiffness, "globalMatrixInitial.csv");
	}
	
	//Check for NAN is rhs and force a cutback
	for (int i1 = 0; i1 < numNodesPerElement * numDOFperNode; i1++)
	{
		tempData = gsl_vector_get(globalRHS, i1);
		if (isnan(tempData))
		{
			isCutback = true;
			break;
		}
	}

	

	gsl_matrix_free(Kelm);
	gsl_vector_free(RHS_elem); 
	gsl_vector_free(U_elem  ); 
	gsl_matrix_free(Coords  );


	return true;
}

bool applyBC(SimModel mySimModel, gsl_spmatrix * globalStiffness, gsl_vector * globalRHS, gsl_vector * U, double TIME, double DTIME,int stepNumber)
{
	int nodeId;
	int globI, globJ;
	int numDofPerNode = mySimModel.modelData.currElement->numDOFSperNode;
	int numNodesInModel = mySimModel.modelData.numNodes;
	int totalDOFinModel = numNodesInModel*numDofPerNode;

	//Check for loads
	int numLoads = mySimModel.historyData.stepSequence[stepNumber].loads.size();

	if (numLoads > 0)
	{
		Load currLoad;
		LoadData currLoadData;
		std::vector <LoadData> currLoadDataSeq;
		gsl_vector_int *nodesInLoadSet;
		double loadMag,loadAmpScale;

		for (int i = 0; i < numLoads; i++)
		{
			currLoad = mySimModel.historyData.stepSequence[stepNumber].loads[i];
			if (currLoad.loadType == "cload")
			{
				currLoadDataSeq = currLoad.loadDataVec;
				for (std::vector<LoadData>::iterator loadIter = currLoadDataSeq.begin(); loadIter != currLoadDataSeq.end(); ++loadIter)
				{
					currLoadData = *loadIter;
					loadMag = currLoadData.mag;
					if (currLoad.isAmplitudeON)
					{
						loadAmpScale = getScaleFromAmp(mySimModel.modelData.amplitudeMap[currLoad.amplitudeName], TIME);
						loadMag = loadMag*loadAmpScale;
					}

					if (currLoadData.isNodeSetUsed)
					{
						nodesInLoadSet = mySimModel.modelData.nodeSetMap[currLoadData.nodeSet];
						for (int i1 = 0; i1 < nodesInLoadSet->size; i1++)
						{
							nodeId = gsl_vector_int_get(nodesInLoadSet, i1);
							globI = (nodeId - 1)*numDofPerNode + currLoadData.dof - 1;
							gsl_vector_set(globalRHS, globI, loadMag);
						}
					}
					else
					{
						nodeId = currLoadData.nodeId;
						globI = (nodeId - 1)*numDofPerNode + currLoadData.dof - 1;
						gsl_vector_set(globalRHS, globI, loadMag);

					}

				}
			}
		}
	}

	//writeGSLMatrix(globalRHS, "RHSAfterLoad.txt");
	//Check for BCs
	int numBCs = mySimModel.historyData.stepSequence[stepNumber].BCs.size();


	if (numBCs > 0)
	{

		std::vector<BCData> currBCDataSeq;
		BC currBC;
		BCData currBCData;
		double bcMag, ampScale, tempVal;
		gsl_vector_int *nodesInSet;
		std::vector<int> constrainedglobEqIds;

		for (int i = 0; i < numBCs; i++)
		{
			currBC = mySimModel.historyData.stepSequence[stepNumber].BCs[i];
			if (currBC.BCType == "boundary") // only supported BC
			{
				currBCDataSeq = currBC.bcData;
				for (std::vector<BCData>::iterator bcIter = currBCDataSeq.begin(); bcIter != currBCDataSeq.end(); ++bcIter)
				{
					currBCData = *bcIter;
					bcMag = currBCData.mag;
					if (currBC.isAmplitudeON)
					{
						ampScale = getScaleFromAmp(mySimModel.modelData.amplitudeMap[currBC.amplitudeName], TIME);
						bcMag = bcMag*ampScale;
					}

					//if node set is used
					if (currBCData.isNodeSetUsed)
					{
						nodesInSet = mySimModel.modelData.nodeSetMap[currBCData.nodeSet];
						for (int i1 = 0; i1 < nodesInSet->size; i1++)
						{
							nodeId = gsl_vector_int_get(nodesInSet, i1);
							for (int dofId = currBCData.startDOF; dofId <= currBCData.endDOF; dofId++)
							{
								globI = (nodeId - 1)*numDofPerNode + dofId - 1;

								//Modify RHS and GlobalK
								gsl_vector_set(globalRHS, globI, bcMag);
								for (int i2 = 0; i2 < totalDOFinModel; i2++)
								{
									gsl_spmatrix_set(globalStiffness, globI, i2, 0.0);
								}
								gsl_spmatrix_set(globalStiffness, globI, globI, 1.0);
								constrainedglobEqIds.push_back(globI);

							}
						}
					}

					//if displacements are specified directly on nodes
					else
					{
						nodeId = currBCData.nodeId;
						for (int dofId = currBCData.startDOF; dofId <= currBCData.endDOF; dofId++)
						{
							globI = (nodeId - 1)*numDofPerNode + dofId - 1;

							//Modify RHS and GlobalK
							gsl_vector_set(globalRHS, globI, bcMag);
							for (int i2 = 0; i2 < totalDOFinModel; i2++)
							{
								gsl_spmatrix_set(globalStiffness, globI, i2, 0.0);
							}
							gsl_spmatrix_set(globalStiffness, globI, globI, 1.0);
							constrainedglobEqIds.push_back(globI);
						}
					}
				}
			}
		}
		//writeGSLMatrix(globalStiffness, "stiffnessAfterBC.csv");
		//writeGSLMatrix(globalRHS, "RHSAfterBC.txt");


		//symmetrize the global stiffness	
		double tempK;
		for (int i = 0; i < constrainedglobEqIds.size(); i++)
		{
			globI = constrainedglobEqIds[i];
			tempVal = gsl_vector_get(globalRHS, globI);

			for (int j = 0; j < totalDOFinModel; j++)
			{
				if (j == globI) continue;
				tempK = gsl_spmatrix_get(globalStiffness, j, globI);
				gsl_vector_set(globalRHS, j, gsl_vector_get(globalRHS, j) - tempVal*tempK);
				gsl_spmatrix_set(globalStiffness, j, globI, 0.0);
			}
		}

		//writeGSLMatrix(globalStiffness, "stiffnessAfterSymmetrize.csv");
		//writeGSLMatrix(globalRHS, "RHSAfterSymmetry.txt");
	}

	//gsl_vector_int_free(nodesInSet);
	return true;
}

double getScaleFromAmp(gsl_matrix * ampArray, double lhs)
{
	//Linear interpolation for amplitudes
	int numRows = ampArray->size1;
	double x1,x2, y1,y2,yReq;
	double tol = 1e-10;
	
	if (lhs < gsl_matrix_get(ampArray, 0, 0)) return gsl_matrix_get(ampArray, 0, 1);
	else if (lhs > gsl_matrix_get(ampArray, numRows - 1, 0)) return gsl_matrix_get(ampArray, numRows - 1, 1);
	else
	{
		for (int i = 0; i < numRows-1; i++)
		{
			x1 = gsl_matrix_get(ampArray, i, 0);
			y1 = gsl_matrix_get(ampArray, i, 1);
			x2 = gsl_matrix_get(ampArray, i+1, 0);
			y2 = gsl_matrix_get(ampArray, i+1, 1);
			
			if ((x1-lhs)>tol) continue;
			else if ((x1 - lhs) < tol && (x2 - lhs) < tol) continue;
			else if (abs(lhs - x1) <= tol) return y1;
			else if (abs(lhs - x2) <= tol) return y2;
			else if ((lhs - x1)>tol && (x2-lhs)>tol)
			{
				yReq = y1 + (y2 - y1)*((lhs - x1) / (x2 - x1));
				return yReq;
			}
		}
	}

	


	
}


bool writeOutputToTxt(SimModel simModel)
{
	std::ofstream outFile;

	std::string modelName = simModel.simModelName;

	outFile.open(modelName + "_Out.txt", std::ios::out | std::ios::trunc);
	//myfile << "Num Rows = " << std::to_string(A_numRows) << std::endl;
	
	//Write Node Coords
	
	outFile << "*NODE START\n";

	int numRows = simModel.modelData.nodeCoordArray->size1;
	int numCols = simModel.modelData.nodeCoordArray->size2;

	for (int i = 0; i < numRows; i++)
	{
		for (int j = 0; j < numCols; j++)
		{
			outFile << gsl_matrix_get(simModel.modelData.nodeCoordArray, i, j)<<",";
		}
		outFile << "\n";
	}
	outFile << "*NODE END\n";

	outFile << "*INTEGRATION POINTS = " << simModel.modelData.currElement->numIntegrationPts<<"\n";
	outFile << "*ELEMENT START\n";

	numRows = simModel.modelData.elementConnectivityArray->size1;
	numCols = simModel.modelData.elementConnectivityArray->size2;

	for (int i = 0; i < numRows; i++)
	{
		outFile << i + 1 << ",";
		for (int j = 0; j < numCols; j++)
		{
			outFile << gsl_matrix_int_get(simModel.modelData.elementConnectivityArray, i, j) << ",";
		}
		outFile << "\n";
	}
	outFile << "*ELEMENT END\n";


	//Results
	
	std::string res;	

	res = "U";
	outFile << "*RESULT START\n";

	
	int numFrames = simModel.results.size();
	for (int frameId = 0; frameId < numFrames; frameId++)
	{
		outFile << "*FRAME START\n";
		outFile << "*TIME =" << simModel.results[frameId].time;

		//Write U
		outFile << "\n*U START\n";		
		int numData = simModel.results[frameId].fields["U"].size();
		for (int i = 0; i < numData; i++)
		{
			outFile << simModel.results[frameId].fields["U"][i]<<",";
		}
		outFile << "\n*U END\n";

		//Write S
		outFile << "*S START\n";
		numData = simModel.results[frameId].fields["S"].size();
		for (int i = 0; i < numData; i++)
		{
			outFile << simModel.results[frameId].fields["S"][i] << ",";
		}
		outFile << "\n*S END\n";
		outFile << "*FRAME END\n";

	}
	
	outFile << "*RESULT END\n";


	outFile.close();

	return true;


}

