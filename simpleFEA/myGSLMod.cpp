#pragma once
#include "myGSLMod.h"
#include <iostream>
#include <fstream>
#include <string>
//Get transpose of a GSL matrix
gsl_matrix * transpose_GSL(const gsl_matrix * A)
{
	int numRows = A->size1;
	int numColumns = A->size2;
	double temp;
	gsl_matrix * transA = gsl_matrix_alloc(numColumns, numRows);

	for (int i = 0; i < numRows; i++)
	{
		for (int j = 0; j < numColumns; j++)
		{
			temp = gsl_matrix_get(A, i, j);
			gsl_matrix_set(transA, j, i,temp);
		}
	}
	return transA;
}

gsl_matrix * multiply_GSL(const gsl_matrix * A, const gsl_matrix * B)
{
	int A_numRows = A->size1;
	int A_numCols = A->size2;

	int B_numRows = B->size1;
	int B_numCols = B->size2;


	if (A_numCols != B_numRows)
	{
		std::cout << "Dimensions for mutliplying matrices are not compatible";
		return nullptr;
	}

	gsl_matrix * C = gsl_matrix_alloc(A_numRows, B_numCols);
	double sum = 0.0;

	for (int i = 0; i < A_numRows; i++)
	{
		for (int j = 0; j <B_numCols; j++)
		{
			sum = 0;
			for (int k = 0; k < A_numCols; k++)
			{
				sum += gsl_matrix_get(A, i, k)*gsl_matrix_get(B, k, j);								
			}
			gsl_matrix_set(C,i, j,sum);
		}
	}

	return C;
}

gsl_vector * multiply_GSL(const gsl_matrix * A, const gsl_vector * B)
{

	int A_numRows = A->size1;
	int A_numCols = A->size2;

	int B_numRows = B->size;
	


	if (A_numCols != B_numRows)
	{
		std::cout << "Dimensions for mutliplying matrices are not compatible";
		return nullptr;
	}

	gsl_vector * C = gsl_vector_alloc(A_numRows);
	double sum = 0.0;

	for (int i = 0; i < A_numRows; i++)
	{
		sum = 0;
		for (int k = 0; k < A_numCols; k++)
		{
			sum += gsl_matrix_get(A, i, k)*gsl_vector_get(B, k);
		}
		gsl_vector_set(C, i ,sum);
		
	}

	return C;
}

bool writeGSLMatrix(const gsl_matrix * A, std::string fileName)
{

	std::ofstream myfile;
	myfile.open(fileName, std::ios::out | std::ios::trunc);
	int A_numRows = A->size1;
	int A_numCols = A->size2;

	//myfile << "Num Rows = "<< std::to_string(A_numRows)<<std::endl;
	//myfile << "Num Columns = " << std::to_string(A_numCols) << std::endl<< std::endl<< std::endl;
	myfile.precision(15);
	myfile << std::showpoint;
	

	for (int i = 0; i < A_numRows; i++)
	{
		for (int j = 0; j < A_numCols; j++)
		{
			myfile << gsl_matrix_get(A,i, j) << ",";
		}
		myfile << std::endl;
	}
	
	myfile << std::endl;

	myfile.close();
	return true;
}

bool writeGSLMatrix(const gsl_matrix_int * A, std::string fileName)
{
	std::ofstream myfile;
	myfile.open(fileName, std::ios::out | std::ios::trunc);
	int A_numRows = A->size1;
	int A_numCols = A->size2;

	//myfile << "Num Rows = " << std::to_string(A_numRows) << std::endl;
	//myfile << "Num Columns = " << std::to_string(A_numCols) << std::endl << std::endl << std::endl;
	
	for (int i = 0; i < A_numRows; i++)
	{
		for (int j = 0; j < A_numCols; j++)
		{
			myfile << gsl_matrix_int_get(A, i, j) << ",";
		}
		myfile << std::endl;
	}

	myfile << std::endl;

	myfile.close();
	return true;
}

bool writeGSLMatrix(const gsl_vector * A, std::string fileName)
{
	std::ofstream myfile;
	myfile.open(fileName, std::ios::out | std::ios::trunc);
	int A_numRows = A->size;
	

	myfile << "Num Rows = " << std::to_string(A_numRows) << std::endl;	

	for (int i = 0; i < A_numRows; i++)
	{
		
		myfile << gsl_vector_get(A, i);		
		myfile << std::endl;
	}

	myfile << std::endl;

	myfile.close();
	return true;
}

bool writeGSLMatrix(const gsl_vector_int * A, std::string fileName)
{
	std::ofstream myfile;
	myfile.open(fileName, std::ios::out | std::ios::trunc);
	int A_numRows = A->size;

	myfile << "Num Rows = " << std::to_string(A_numRows) << std::endl;

	for (int i = 0; i < A_numRows; i++)
	{

		myfile << gsl_vector_int_get(A, i);
		myfile << std::endl;
	}

	myfile << std::endl;

	myfile.close();
	return true;
}

bool writeGSLMatrix(const gsl_spmatrix * A, std::string fileName)
{

	std::ofstream myfile;
	myfile.open(fileName, std::ios::out | std::ios::trunc);
	int A_numRows = A->size1;
	int A_numCols = A->size2;

	//myfile << "Num Rows = " << std::to_string(A_numRows) << std::endl;
	//myfile << "Num Columns = " << std::to_string(A_numCols) << std::endl << std::endl << std::endl;

	for (int i = 0; i < A_numRows; i++)
	{
		for (int j = 0; j < A_numCols; j++)
		{
			myfile << gsl_spmatrix_get(A, i, j) <<",";
		}
		myfile << std::endl;
	}

	myfile << std::endl;

	myfile.close();
	return true;

	
}




