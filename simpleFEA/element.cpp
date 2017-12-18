#include "element.h"
#include <iostream>
#include "myGSLMod.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
Element::Element()
{
	elementTypeName = "";
}

bool Element::get3DIntegrationPoints(const int n_points, const int n_nodes, gsl_matrix * xi, gsl_vector * w)
{
	
	return false;
}

bool Element::get3DShapeFunctions(const gsl_matrix * xi, const int n_nodes, gsl_vector * f, gsl_matrix * df)
{
	return false;
}

bool Element::get2DIntegrationPoints(const int n_points, const int n_nodes, gsl_matrix **xi_p, gsl_vector **w_p)
{
	
	
	int i, j, k, n;
	double cn, w1, w2, w11, w12, w22;
	//n_points -- No. of integration points
	//n_nodes  -- No. of nodes per element
	//xi  -- Returns the isoparametric coordinates of the integration points
	//w -- Returns the weight factor

	//Defines integration points and weights for 2D continuum elements
	//gsl_matrix * xi = *xi_p;
	//gsl_vector * w = *w_p;

	*xi_p = gsl_matrix_calloc(n_points, 2);
	*w_p = gsl_vector_calloc(n_points);



	if (n_points == 1)
	{

		if (n_nodes == 4 || n_nodes == 9)   //4 or 9 noded quad
		{
			gsl_matrix_set(*xi_p, 0, 0, 0.0);
			gsl_matrix_set(*xi_p, 0, 1, 0.0);
			gsl_vector_set(*w_p, 0, 4.0);
		}
		else if (n_nodes == 3 || n_nodes == 6)
		{// 3 or 6 noded triangle
			gsl_matrix_set(*xi_p, 0, 0, 1.0 / 3.0);
			gsl_matrix_set(*xi_p, 0, 1, 1.0 / 3.0);
			gsl_vector_set(*w_p, 0, 0.5);
		}
	}
	else if (n_points == 3)
	{
		//*xi_p
		gsl_matrix_set(*xi_p, 0, 0, 0.50);
		gsl_matrix_set(*xi_p, 0, 1, 0.50);
		gsl_matrix_set(*xi_p, 1, 0, 0.0);
		gsl_matrix_set(*xi_p, 1, 1, 0.50);
		gsl_matrix_set(*xi_p, 2, 0, 0.50);
		gsl_matrix_set(*xi_p, 2, 1, 0.0);

		//*w_p
		gsl_vector_set(*w_p, 0, 1.0 / 6.0);
		gsl_vector_set(*w_p, 1, 1.0 / 6.0);
		gsl_vector_set(*w_p, 2, 1.0 / 6.0);

	}
	else if (n_points == 4)
	{
		if (n_nodes == 4 || n_nodes == 8 || n_nodes == 9)
		{

			cn = 0.57735026918962600;
			gsl_matrix_set(*xi_p, 0, 0, -cn);
			gsl_matrix_set(*xi_p, 1, 0, cn);
			gsl_matrix_set(*xi_p, 2, 0, cn);
			gsl_matrix_set(*xi_p, 3, 0, -cn);
			gsl_matrix_set(*xi_p, 0, 1, -cn);
			gsl_matrix_set(*xi_p, 1, 1, -cn);
			gsl_matrix_set(*xi_p, 2, 1, cn);
			gsl_matrix_set(*xi_p, 3, 1, cn);
			gsl_vector_set(*w_p, 0, 1.0);
			gsl_vector_set(*w_p, 1, 1.0);
			gsl_vector_set(*w_p, 2, 1.0);
			gsl_vector_set(*w_p, 3, 1.0);
		}

		else if (n_nodes == 3 || n_nodes == 6)
		{
			gsl_matrix_set(*xi_p, 0, 0, 1.0 / 3.0);
			gsl_matrix_set(*xi_p, 0, 1, 1.0 / 3.0);
			gsl_matrix_set(*xi_p, 1, 0, 0.60);
			gsl_matrix_set(*xi_p, 1, 1, 0.20);
			gsl_matrix_set(*xi_p, 2, 0, 0.20);
			gsl_matrix_set(*xi_p, 2, 1, 0.60);
			gsl_matrix_set(*xi_p, 3, 0, 0.20);
			gsl_matrix_set(*xi_p, 3, 1, 0.20);

			gsl_vector_set(*w_p, 0, -27.0 / 96.0);
			gsl_vector_set(*w_p, 1, 25.0 / 96.0);
			gsl_vector_set(*w_p, 2, 25.0 / 96.0);
			gsl_vector_set(*w_p, 3, 25.0 / 96.0);
		}
	}


	else if (n_points == 7)
	{
		//Quintic integration for triangle
		gsl_matrix_set(*xi_p, 0, 0, 1.0 / 3.0);
		gsl_matrix_set(*xi_p, 0, 1, 1.0 / 3.0);
		gsl_matrix_set(*xi_p, 1, 0, 0.05971587170);
		gsl_matrix_set(*xi_p, 1, 1, 0.47014206410);
		gsl_matrix_set(*xi_p, 2, 0, 0.47014206410);
		gsl_matrix_set(*xi_p, 2, 1, 0.05971587170);
		gsl_matrix_set(*xi_p, 3, 0, 0.47014206410);
		gsl_matrix_set(*xi_p, 3, 1, 0.47014206410);
		gsl_matrix_set(*xi_p, 4, 0, 0.79742698530);
		gsl_matrix_set(*xi_p, 4, 1, 0.10128650730);
		gsl_matrix_set(*xi_p, 5, 0, 0.10128650730);
		gsl_matrix_set(*xi_p, 5, 1, 0.79742698530);
		gsl_matrix_set(*xi_p, 6, 0, 0.10128650730);
		gsl_matrix_set(*xi_p, 6, 1, 0.10128650730);

		gsl_vector_set(*w_p, 0, 0.11250);
		gsl_vector_set(*w_p, 1, 0.06619707630);
		gsl_vector_set(*w_p, 2, 0.06619707630);
		gsl_vector_set(*w_p, 3, 0.06619707630);
		gsl_vector_set(*w_p, 4, 0.06296959020);
		gsl_vector_set(*w_p, 5, 0.06296959020);
		gsl_vector_set(*w_p, 6, 0.06296959020);
	}


	else if (n_points == 9)
	{
		cn = 0.77459666924148300;
		gsl_matrix_set(*xi_p, 0, 0, -cn);
		gsl_matrix_set(*xi_p, 1, 0, 0.0);
		gsl_matrix_set(*xi_p, 2, 0, cn);
		gsl_matrix_set(*xi_p, 3, 0, -cn);
		gsl_matrix_set(*xi_p, 4, 0, 0.0);
		gsl_matrix_set(*xi_p, 5, 0, cn);
		gsl_matrix_set(*xi_p, 6, 0, -cn);
		gsl_matrix_set(*xi_p, 7, 0, 0.0);
		gsl_matrix_set(*xi_p, 8, 0, cn);
		gsl_matrix_set(*xi_p, 0, 1, -cn);
		gsl_matrix_set(*xi_p, 1, 1, -cn);
		gsl_matrix_set(*xi_p, 2, 1, -cn);
		gsl_matrix_set(*xi_p, 3, 1, 0.0);
		gsl_matrix_set(*xi_p, 4, 1, 0.0);
		gsl_matrix_set(*xi_p, 5, 1, 0.0);
		gsl_matrix_set(*xi_p, 6, 1, cn);
		gsl_matrix_set(*xi_p, 7, 1, cn);
		gsl_matrix_set(*xi_p, 8, 1, cn);


		w1 = 0.55555555555555600;
		w2 = 0.88888888888888900;
		w11 = w1*w1;
		w12 = w1*w2;
		w22 = w2*w2;
		gsl_vector_set(*w_p, 0, w11);
		gsl_vector_set(*w_p, 1, w12);
		gsl_vector_set(*w_p, 2, w11);
		gsl_vector_set(*w_p, 3, w12);
		gsl_vector_set(*w_p, 4, w22);
		gsl_vector_set(*w_p, 5, w12);
		gsl_vector_set(*w_p, 6, w11);
		gsl_vector_set(*w_p, 7, w12);
		gsl_vector_set(*w_p, 8, w11);
	}
	return true;

}

bool Element::get2DShapeFunctions(const std::vector<double> xi, const int n_nodes, gsl_vector **f, gsl_matrix **df)
{

	double g1, g2, g3, dg1, dg2, dg3;
	double h1, h2, h3, dh1, dh2, dh3;
	double z, dzdp, dzdq;
	double xi1, xi2;
	xi1 = xi[0];
	xi2 = xi[1];
	
	*f = gsl_vector_calloc(n_nodes);
	*df = gsl_matrix_calloc(n_nodes, 2);
	if (n_nodes == 3)
	{
		gsl_vector_set(*f, 0, xi1);
		gsl_vector_set(*f, 1, xi2);
		gsl_vector_set(*f, 2, 1.0 - xi1 - xi2);

		gsl_matrix_set(*df , 0, 0 , 1.0 );
		gsl_matrix_set(*df , 0, 1 , 0.0 );
		gsl_matrix_set(*df , 1, 0, 0.0);
		gsl_matrix_set(*df , 1, 1, 1.0 );
		gsl_matrix_set(*df , 2, 0, -1.0);
		gsl_matrix_set(*df , 2, 1, -1.0);
	}
	else if (n_nodes == 4)
	{
		g1 = 0.50*(1.0 - xi1);
		g2 = 0.50*(1.0 + xi1);
		h1 = 0.50*(1.0 - xi2);
		h2 = 0.50*(1.0 + xi2);

		gsl_vector_set(*f, 0, g1*h1);
		gsl_vector_set(*f, 1, g2*h1);
		gsl_vector_set(*f, 2, g2*h2);
		gsl_vector_set(*f, 3, g1*h2);

		dg1 = -0.50;
		dg2 = 0.50;
		dh1 = -0.50;
		dh2 = 0.50;
		gsl_matrix_set(*df, 0, 0, dg1*h1);
		gsl_matrix_set(*df, 1, 0, dg2*h1);
		gsl_matrix_set(*df, 2, 0, dg2*h2);
		gsl_matrix_set(*df, 3, 0, dg1*h2);
		gsl_matrix_set(*df, 0, 1, g1*dh1);
		gsl_matrix_set(*df, 1, 1, g2*dh1);
		gsl_matrix_set(*df, 2, 1, g2*dh2);
		gsl_matrix_set(*df, 3, 1, g1*dh2);
	}

	else if (n_nodes == 6)
	{
		z = 1.0 - xi1 - xi2;
		gsl_vector_set(*f,0,(2.0*xi1 - 1.0)*xi1);
		gsl_vector_set(*f, 1, (2.0*xi2 - 1.0)*xi2);
		gsl_vector_set(*f, 2, (2.0*z - 1.0)*z);
		gsl_vector_set(*f, 3, 4.0*xi1*xi2);
		gsl_vector_set(*f, 4, 4.0*xi2*z);
		gsl_vector_set(*f, 5,4.0*xi1*z		  );

		dzdp = -1.0;
		dzdq = -1.0;
		gsl_matrix_set(*df, 0, 0,4.0*xi1 - 1.0	   );
		gsl_matrix_set(*df, 1, 0, 0.0);
		gsl_matrix_set(*df, 2, 0, 4.0*z*dzdp - dzdp);
		gsl_matrix_set(*df, 3, 0, 4.0*xi2);
		gsl_matrix_set(*df, 4, 0, 4.0*xi2*dzdp);
		gsl_matrix_set(*df, 5, 0,4.0*z + 4.0*xi1*dzdp);
		gsl_matrix_set(*df, 0, 1, 0.0);
		gsl_matrix_set(*df, 1, 1, 4.0*xi2 - 1.0);
		gsl_matrix_set(*df, 2, 1, 4.0*z*dzdq - dzdq);
		gsl_matrix_set(*df, 3, 1, 4.0*xi1);
		gsl_matrix_set(*df, 4, 1, 4.0*z + 4.0*xi2*dzdq);
		gsl_matrix_set(*df, 5, 1,4.0*xi1*dzdq        );
	}

	else if (n_nodes == 8)
	{
		gsl_vector_set(*f, 0,-0.25*(1. - xi1)*(1. - xi2)*(1. + xi1 + xi2));
		gsl_vector_set(*f, 1, 0.25*(1. + xi1)*(1. - xi2)*(xi1 - xi2 - 1.));
		gsl_vector_set(*f, 2, 0.25*(1. + xi1)*(1. + xi2)*(xi1 + xi2 - 1.));
		gsl_vector_set(*f, 3, 0.25*(1. - xi1)*(1. + xi2)*(xi2 - xi1 - 1.));
		gsl_vector_set(*f, 4, 0.5*(1. - xi1*xi1)*(1. - xi2));
		gsl_vector_set(*f, 5, 0.5*(1. + xi1)*(1. - xi2*xi2));
		gsl_vector_set(*f, 6, 0.5*(1. - xi1*xi1)*(1. + xi2));
		gsl_vector_set(*f, 7,0.5*(1. - xi1)*(1. - xi2*xi2));

		gsl_matrix_set(*df, 0, 0, 0.25*(1. - xi2)*(2.*xi1 + xi2));
		gsl_matrix_set(*df, 0, 1, 0.25*(1. - xi1)*(xi1 + 2.*xi2));
		gsl_matrix_set(*df, 1, 0, 0.25*(1. - xi2)*(2.*xi1 - xi2));
		gsl_matrix_set(*df, 1, 1, 0.25*(1. + xi1)*(2.*xi2 - xi1));
		gsl_matrix_set(*df, 2, 0, 0.25*(1. + xi2)*(2.*xi1 + xi2));
		gsl_matrix_set(*df, 2, 1, 0.25*(1. + xi1)*(2.*xi2 + xi1));
		gsl_matrix_set(*df, 3, 0, 0.25*(1. + xi2)*(2.*xi1 - xi2));
		gsl_matrix_set(*df, 3, 1, 0.25*(1. - xi1)*(2.*xi2 - xi1));
		gsl_matrix_set(*df, 4, 0, -xi1*(1. - xi2));
		gsl_matrix_set(*df, 4, 1, -0.5*(1. - xi1*xi1));
		gsl_matrix_set(*df, 5, 0, 0.5*(1. - xi2*xi2));
		gsl_matrix_set(*df, 5, 1, -(1. + xi1)*xi2);
		gsl_matrix_set(*df, 6, 0, -xi1*(1. + xi2));
		gsl_matrix_set(*df, 6, 1, 0.5*(1. - xi1*xi1));
		gsl_matrix_set(*df, 7, 0, -0.5*(1. - xi2*xi2));
		gsl_matrix_set(*df, 7, 1, -(1. - xi1)*xi2);
	}
	else if (n_nodes == 9)
	{
		g1 = -.50*xi1*(1.0 - xi1);
		g2 = (1.0 - xi1)*(1.0 + xi1);
		g3 = .50*xi1*(1.0 + xi1);
		h1 = -.50*xi2*(1.0 - xi2);
		h2 = (1.0 - xi2)*(1.0 + xi2);
		h3 = .50*xi2*(1.0 + xi2);
		dg1 = xi1 - 0.50;
		dg2 = -2.0*xi1;
		dg3 = xi1 + 0.50;
		dh1 = xi2 - 0.50;
		dh2 = -2.0*xi2;
		dh3 = xi2 + 0.50;


		gsl_vector_set(*f,0,g1*h1)  ;
		gsl_vector_set(*f, 1, g2*h1);
		gsl_vector_set(*f, 2, g3*h1);
		gsl_vector_set(*f, 3, g1*h2);
		gsl_vector_set(*f, 4, g2*h2);
		gsl_vector_set(*f, 5, g3*h2);
		gsl_vector_set(*f, 6, g1*h3);
		gsl_vector_set(*f, 7, g2*h3);
		gsl_vector_set(*f, 8,g3*h3) ;

		gsl_matrix_set(*df,0, 0, dg1*h1);
		gsl_matrix_set(*df,0, 1, g1*dh1);
		gsl_matrix_set(*df, 1, 0, dg2*h1);
		gsl_matrix_set(*df, 1, 1, g2*dh1);
		gsl_matrix_set(*df, 2, 0, dg3*h1);
		gsl_matrix_set(*df, 2, 1, g3*dh1);
		gsl_matrix_set(*df, 3, 0, dg1*h2);
		gsl_matrix_set(*df, 3, 1, g1*dh2);
		gsl_matrix_set(*df, 4, 0, dg2*h2);
		gsl_matrix_set(*df, 4, 1, g2*dh2);
		gsl_matrix_set(*df, 5, 0, dg3*h2);
		gsl_matrix_set(*df, 5, 1, g3*dh2);
		gsl_matrix_set(*df, 6, 0, dg1*h3);
		gsl_matrix_set(*df, 6, 0, g1*dh3);
		gsl_matrix_set(*df, 7, 1, dg2*h3);
		gsl_matrix_set(*df, 7, 0, g2*dh3);
		gsl_matrix_set(*df, 8, 1, dg3*h3);
		gsl_matrix_set(*df, 8, 0, g3*dh3);
	}



	return false;
}

bool Element::getElementStiffnessAndResidualForce(gsl_vector * RHS, gsl_matrix * Kelm, double Time, double Dtime, gsl_matrix * Coords, gsl_vector * U, int NNODE, std::map<std::string, std::vector<double>> &fields,int sVarsStart)
{
	
	return false;
}

Element::~Element()
{
}

Element_CPE8::Element_CPE8()
{
	elementTypeName = "CPE8";
	numDOFS = 16;
	numDOFSperNode = 2;
	spaceDimensions = 2;
	numIntegrationPts = 9;
	numNodes = 8;
	stateVarsPerIntPoint = 4;

}

bool Element_CPE8::getElementStiffnessAndResidualForce(gsl_vector * RHS, gsl_matrix * Kelm, double Time, double Dtime, gsl_matrix * Coords, gsl_vector * Uelm, int NNODE, std::map<std::string, std::vector<double>> &fields,int sVarsStart)
{


	int n_points;

	//writeGSLMatrix(Kelm, "KelmBefore.txt");
	
	

	n_points = numIntegrationPts;

	gsl_matrix * D = gsl_matrix_calloc(4, 4);
	gsl_vector * tempVec;
	gsl_vector * strain;
	gsl_vector * stress;
	gsl_matrix * tempMat;
	gsl_matrix * xi;
	gsl_vector * w;
	gsl_vector * N;
	gsl_matrix * dNdxi;
	std::vector<double> stressVector;

	//if (fields.find("S") != fields.end())
	//{
	//	stressVector = fields["S"];
	//}
	//


	//Properties are hard coded for now
	double xnu =.3 ;
	double E = 100.0;
	double d44, d11, d12;
	



	gsl_matrix_set_zero(Kelm);
	gsl_vector_set_zero(RHS);


	d44 = 0.50*E / (1 + xnu);
	d11 = (1.0 - xnu)*E / ((1.0 + xnu)*(1.0 - 2.0*xnu));
	d12 = xnu*E / ((1.0 + xnu)*(1.0 - 2.0*xnu));

	gsl_matrix_set (D,0,1,d12);
	gsl_matrix_set (D,0,2,d12);
	gsl_matrix_set (D,1,0,d12);
	gsl_matrix_set (D,1,2,d12);
	gsl_matrix_set (D,2,0,d12);
	gsl_matrix_set (D,2,1,d12);
	gsl_matrix_set(D, 0, 0, d11);
	gsl_matrix_set(D, 1, 1, d11);
	gsl_matrix_set(D, 2, 2, d11);
	gsl_matrix_set(D, 3, 3, d44);


	//writeGSLMatrix(D, "DMatrix.txt");
	//writeGSLMatrix(Uelm, "U_element.txt");
	//gsl_matrix * xi = nullptr;
	//gsl_vector * w = nullptr;


	bool getIntegrationPts = get2DIntegrationPoints(n_points, NNODE,&xi, &w);

	//std::cout << xi->size1;

	//writeGSLMatrix(xi, "xi.text");
	//writeGSLMatrix(w, "xi.text");

	
	if (!getIntegrationPts) 
	{
		std::cout << "the integration points could not be obtained";
		return false;
	}
	
	std::vector<double> xiTemp;
	xiTemp.push_back(0.0);
	xiTemp.push_back(0.0);

	int sVarIter = 0;
	

	for (int kint=0;kint<n_points;kint++)
	{					
		xiTemp[0] = gsl_matrix_get(xi, kint, 0);
		xiTemp[1] = gsl_matrix_get(xi, kint, 1);

		get2DShapeFunctions(xiTemp, NNODE, &N, &dNdxi);

		gsl_matrix * dxdxi = multiply_GSL(transpose_GSL(Coords), dNdxi);
		//writeGSLMatrix(Coords, "coords.txt");
		//writeGSLMatrix(dNdxi, "dNdxi.txt");	
		//writeGSLMatrix(dxdxi, "dxdxi.txt");
		//writeGSLMatrix(N, "N.txt");



		double determinant = gsl_matrix_get(dxdxi, 0, 0)*gsl_matrix_get(dxdxi, 1, 1) - gsl_matrix_get(dxdxi, 0, 1)*gsl_matrix_get(dxdxi, 1, 0);

		gsl_matrix * dxidx = gsl_matrix_alloc(2, 2);

		gsl_matrix_set(dxidx, 0, 0,gsl_matrix_get(dxdxi, 1, 1) / determinant);
		gsl_matrix_set(dxidx, 0, 1, -1.0*gsl_matrix_get(dxdxi, 0, 1) / determinant);
		gsl_matrix_set(dxidx, 1, 0, -1.0*gsl_matrix_get(dxdxi, 1, 0) / determinant);
		gsl_matrix_set(dxidx, 1, 1, gsl_matrix_get(dxdxi, 0, 0) / determinant);

		//writeGSLMatrix(dxidx, "dxidx.txt");

		gsl_matrix * dNdx = multiply_GSL(dNdxi, dxidx);
		//writeGSLMatrix(dNdx, "dNdx.txt");


		//----------------
		//Define the B Matrix
		//--------------------
		gsl_matrix * B = gsl_matrix_calloc(4,2 * NNODE);		

		for (int i = 0; i < NNODE; i++)
		{
			gsl_matrix_set(B, 0, 2 * i, gsl_matrix_get(dNdx, i, 0));
			gsl_matrix_set(B, 3, 2 * i, gsl_matrix_get(dNdx, i, 1));
		}

		for (int i = 1; i <= NNODE; i++)
		{				
			gsl_matrix_set(B, 1, (i * 2) - 1, gsl_matrix_get(dNdx, i-1, 1));
			gsl_matrix_set(B, 3, (i*2)-1, gsl_matrix_get(dNdx, i-1, 0));						
		}

		//writeGSLMatrix(B, "BMatrix_"+std::to_string(kint)+".csv");


		//Update RHS
		strain = multiply_GSL(B, Uelm);

		//writeGSLMatrix(strain, "strain_" + std::to_string(kint) + ".txt");
		

		stress = multiply_GSL(D, strain);

		//writeGSLMatrix(stress, "stress_" + std::to_string(kint) + ".txt");


		tempVec = multiply_GSL(transpose_GSL(B), stress);
		gsl_vector_scale(tempVec, gsl_vector_get(w, kint)*determinant*-1.0);
		gsl_vector_add(RHS, tempVec);


		//Update element stiffness
		tempMat = multiply_GSL(transpose_GSL(B), multiply_GSL(D, B));
		gsl_matrix_scale(tempMat, gsl_vector_get(w, kint)*determinant);
		gsl_matrix_add(Kelm, tempMat);

		//write stress to element stress vector
		for (int i = 0; i < stress->size; i++)
		{
			//stressVector.push_back(gsl_vector_get(stress, i));
			fields["S"][sVarsStart + sVarIter] = gsl_vector_get(stress, i);
			sVarIter++;
			//std::cout << sVarIter<<" , " << sVarsStart <<std::endl;
		}
			
	}
	
	//writeGSLMatrix(Kelm, "Kelement.txt");
	//writeGSLMatrix(RHS, "RHSelement.txt");

	//Write the results to map.

	
	

	gsl_matrix_free(D);
	gsl_vector_free(tempVec);
	gsl_vector_free(strain);
	gsl_vector_free(stress);
	gsl_matrix_free(tempMat);
	gsl_matrix_free(xi);
	gsl_vector_free(w);
	gsl_vector_free(N);
	gsl_matrix_free(dNdxi);

	return true;
}

Element_CPE8::~Element_CPE8()
{
}


//
//Triangle elements
//

Element_CPE6::Element_CPE6()
{
	elementTypeName = "CPE6";
	numDOFS = 12;
	numDOFSperNode = 2;
	spaceDimensions = 2;
	numIntegrationPts = 3;
	numNodes = 6;
	stateVarsPerIntPoint = 4;

}

bool Element_CPE6::getElementStiffnessAndResidualForce(gsl_vector * RHS, gsl_matrix * Kelm, double Time, double Dtime, gsl_matrix * Coords, gsl_vector * Uelm, int NNODE, std::map<std::string, std::vector<double>> &fields, int sVarsStart)
{


	int n_points;

	//writeGSLMatrix(Kelm, "KelmBefore.txt");

	//if (NNODE == 3) n_points = 1;
	//if (NNODE == 4) n_points = 4;
	//if (NNODE == 6) n_points = 4;
	//if (NNODE == 8) n_points = 9;
	//if (NNODE == 9) n_points = 9;

	n_points = numIntegrationPts;

	gsl_matrix * D = gsl_matrix_calloc(4, 4);
	gsl_vector * tempVec;
	gsl_vector * strain;
	gsl_vector * stress;
	gsl_matrix * tempMat;
	gsl_matrix * xi;
	gsl_vector * w;
	gsl_vector * N;
	gsl_matrix * dNdxi;
	std::vector<double> stressVector;

	//if (fields.find("S") != fields.end())
	//{
	//	stressVector = fields["S"];
	//}
	//


	//Properties are hard coded for now
	double xnu = .3;
	double E = 100.0;
	double d44, d11, d12;


	gsl_matrix_set_zero(Kelm);
	gsl_vector_set_zero(RHS);


	d44 = 0.50*E / (1 + xnu);
	d11 = (1.0 - xnu)*E / ((1.0 + xnu)*(1.0 - 2.0*xnu));
	d12 = xnu*E / ((1.0 + xnu)*(1.0 - 2.0*xnu));

	gsl_matrix_set(D, 0, 1, d12);
	gsl_matrix_set(D, 0, 2, d12);
	gsl_matrix_set(D, 1, 0, d12);
	gsl_matrix_set(D, 1, 2, d12);
	gsl_matrix_set(D, 2, 0, d12);
	gsl_matrix_set(D, 2, 1, d12);
	gsl_matrix_set(D, 0, 0, d11);
	gsl_matrix_set(D, 1, 1, d11);
	gsl_matrix_set(D, 2, 2, d11);
	gsl_matrix_set(D, 3, 3, d44);


	//writeGSLMatrix(D, "DMatrix.txt");
	//writeGSLMatrix(Uelm, "U_element.txt");
	//gsl_matrix * xi = nullptr;
	//gsl_vector * w = nullptr;


	bool getIntegrationPts = get2DIntegrationPoints(n_points, NNODE, &xi, &w);

	//std::cout << xi->size1;

	//writeGSLMatrix(xi, "xi.text");
	//writeGSLMatrix(w, "xi.text");


	if (!getIntegrationPts)
	{
		std::cout << "the integration points could not be obtained";
		return false;
	}

	std::vector<double> xiTemp;
	xiTemp.push_back(0.0);
	xiTemp.push_back(0.0);

	int sVarIter = 0;


	for (int kint = 0; kint<n_points; kint++)
	{
		xiTemp[0] = gsl_matrix_get(xi, kint, 0);
		xiTemp[1] = gsl_matrix_get(xi, kint, 1);

		get2DShapeFunctions(xiTemp, NNODE, &N, &dNdxi);

		gsl_matrix * dxdxi = multiply_GSL(transpose_GSL(Coords), dNdxi);
		//writeGSLMatrix(Coords, "coords.txt");
		//writeGSLMatrix(dNdxi, "dNdxi.txt");	
		//writeGSLMatrix(dxdxi, "dxdxi.txt");
		//writeGSLMatrix(N, "N.txt");



		double determinant = gsl_matrix_get(dxdxi, 0, 0)*gsl_matrix_get(dxdxi, 1, 1) - gsl_matrix_get(dxdxi, 0, 1)*gsl_matrix_get(dxdxi, 1, 0);

		gsl_matrix * dxidx = gsl_matrix_alloc(2, 2);

		gsl_matrix_set(dxidx, 0, 0, gsl_matrix_get(dxdxi, 1, 1) / determinant);
		gsl_matrix_set(dxidx, 0, 1, -1.0*gsl_matrix_get(dxdxi, 0, 1) / determinant);
		gsl_matrix_set(dxidx, 1, 0, -1.0*gsl_matrix_get(dxdxi, 1, 0) / determinant);
		gsl_matrix_set(dxidx, 1, 1, gsl_matrix_get(dxdxi, 0, 0) / determinant);

		//writeGSLMatrix(dxidx, "dxidx.txt");

		gsl_matrix * dNdx = multiply_GSL(dNdxi, dxidx);
		//writeGSLMatrix(dNdx, "dNdx.txt");


		//----------------
		//Define the B Matrix
		//--------------------
		gsl_matrix * B = gsl_matrix_calloc(4, 2 * NNODE);

		for (int i = 0; i < NNODE; i++)
		{
			gsl_matrix_set(B, 0, 2 * i, gsl_matrix_get(dNdx, i, 0));
			gsl_matrix_set(B, 3, 2 * i, gsl_matrix_get(dNdx, i, 1));
		}

		for (int i = 1; i <= NNODE; i++)
		{
			gsl_matrix_set(B, 1, (i * 2) - 1, gsl_matrix_get(dNdx, i - 1, 1));
			gsl_matrix_set(B, 3, (i * 2) - 1, gsl_matrix_get(dNdx, i - 1, 0));
		}

		//writeGSLMatrix(B, "BMatrix_"+std::to_string(kint)+".csv");


		//Update RHS
		strain = multiply_GSL(B, Uelm);

		//writeGSLMatrix(strain, "strain_" + std::to_string(kint) + ".txt");


		stress = multiply_GSL(D, strain);

		//writeGSLMatrix(stress, "stress_" + std::to_string(kint) + ".txt");


		tempVec = multiply_GSL(transpose_GSL(B), stress);
		gsl_vector_scale(tempVec, gsl_vector_get(w, kint)*determinant*-1.0);
		gsl_vector_add(RHS, tempVec);


		//Update element stiffness
		tempMat = multiply_GSL(transpose_GSL(B), multiply_GSL(D, B));
		gsl_matrix_scale(tempMat, gsl_vector_get(w, kint)*determinant);
		gsl_matrix_add(Kelm, tempMat);

		//write stress to element stress vector
		for (int i = 0; i < stress->size; i++)
		{
			//stressVector.push_back(gsl_vector_get(stress, i));
			fields["S"][sVarsStart + sVarIter] = gsl_vector_get(stress, i);
			sVarIter++;
			//std::cout << sVarIter<<" , " << sVarsStart <<std::endl;
		}

	}

	//writeGSLMatrix(Kelm, "Kelement.txt");
	//writeGSLMatrix(RHS, "RHSelement.txt");

	//Write the results to map.

	gsl_matrix_free(D);
	gsl_vector_free(tempVec);
	gsl_vector_free(strain);
	gsl_vector_free(stress);
	gsl_matrix_free(tempMat);
	gsl_matrix_free(xi);
	gsl_vector_free(w);
	gsl_vector_free(N);
	gsl_matrix_free(dNdxi);

	return true;
}

Element_CPE6::~Element_CPE6()
{
}


//3 NODED TRI - Linear
Element_CPE3::Element_CPE3()
{
	elementTypeName = "CPE3";
	numDOFS = 6;
	numDOFSperNode = 2;
	spaceDimensions = 2;
	numIntegrationPts = 1;
	numNodes = 3;
	stateVarsPerIntPoint = 4;

}

bool Element_CPE3::getElementStiffnessAndResidualForce(gsl_vector * RHS, gsl_matrix * Kelm, double Time, double Dtime, gsl_matrix * Coords, gsl_vector * Uelm, int NNODE, std::map<std::string, std::vector<double>> &fields, int sVarsStart)
{


	int n_points = numIntegrationPts;;

	gsl_matrix * D = gsl_matrix_calloc(4, 4);
	gsl_vector * tempVec;
	gsl_vector * strain;
	gsl_vector * stress;
	gsl_matrix * tempMat;
	gsl_matrix * xi;
	gsl_vector * w;
	gsl_vector * N;
	gsl_matrix * dNdxi;
	std::vector<double> stressVector;

	//if (fields.find("S") != fields.end())
	//{
	//	stressVector = fields["S"];
	//}
	//


	//Properties are hard coded for now
	double xnu = .3;
	double E = 100.0;
	double d44, d11, d12;


	gsl_matrix_set_zero(Kelm);
	gsl_vector_set_zero(RHS);


	d44 = 0.50*E / (1 + xnu);
	d11 = (1.0 - xnu)*E / ((1.0 + xnu)*(1.0 - 2.0*xnu));
	d12 = xnu*E / ((1.0 + xnu)*(1.0 - 2.0*xnu));

	gsl_matrix_set(D, 0, 1, d12);
	gsl_matrix_set(D, 0, 2, d12);
	gsl_matrix_set(D, 1, 0, d12);
	gsl_matrix_set(D, 1, 2, d12);
	gsl_matrix_set(D, 2, 0, d12);
	gsl_matrix_set(D, 2, 1, d12);
	gsl_matrix_set(D, 0, 0, d11);
	gsl_matrix_set(D, 1, 1, d11);
	gsl_matrix_set(D, 2, 2, d11);
	gsl_matrix_set(D, 3, 3, d44);


	//writeGSLMatrix(D, "DMatrix.txt");
	//writeGSLMatrix(Uelm, "U_element.txt");
	//gsl_matrix * xi = nullptr;
	//gsl_vector * w = nullptr;


	bool getIntegrationPts = get2DIntegrationPoints(n_points, NNODE, &xi, &w);

	//std::cout << xi->size1;

	//writeGSLMatrix(xi, "xi.text");
	//writeGSLMatrix(w, "xi.text");


	if (!getIntegrationPts)
	{
		std::cout << "the integration points could not be obtained";
		return false;
	}

	std::vector<double> xiTemp;
	xiTemp.push_back(0.0);
	xiTemp.push_back(0.0);

	int sVarIter = 0;


	for (int kint = 0; kint<n_points; kint++)
	{
		xiTemp[0] = gsl_matrix_get(xi, kint, 0);
		xiTemp[1] = gsl_matrix_get(xi, kint, 1);

		get2DShapeFunctions(xiTemp, NNODE, &N, &dNdxi);

		gsl_matrix * dxdxi = multiply_GSL(transpose_GSL(Coords), dNdxi);
		//writeGSLMatrix(Coords, "coords.txt");
		//writeGSLMatrix(dNdxi, "dNdxi.txt");	
		//writeGSLMatrix(dxdxi, "dxdxi.txt");
		//writeGSLMatrix(N, "N.txt");



		double determinant = gsl_matrix_get(dxdxi, 0, 0)*gsl_matrix_get(dxdxi, 1, 1) - gsl_matrix_get(dxdxi, 0, 1)*gsl_matrix_get(dxdxi, 1, 0);

		gsl_matrix * dxidx = gsl_matrix_alloc(2, 2);

		gsl_matrix_set(dxidx, 0, 0, gsl_matrix_get(dxdxi, 1, 1) / determinant);
		gsl_matrix_set(dxidx, 0, 1, -1.0*gsl_matrix_get(dxdxi, 0, 1) / determinant);
		gsl_matrix_set(dxidx, 1, 0, -1.0*gsl_matrix_get(dxdxi, 1, 0) / determinant);
		gsl_matrix_set(dxidx, 1, 1, gsl_matrix_get(dxdxi, 0, 0) / determinant);

		//writeGSLMatrix(dxidx, "dxidx.txt");

		gsl_matrix * dNdx = multiply_GSL(dNdxi, dxidx);
		//writeGSLMatrix(dNdx, "dNdx.txt");


		//----------------
		//Define the B Matrix
		//--------------------
		gsl_matrix * B = gsl_matrix_calloc(4, 2 * NNODE);

		for (int i = 0; i < NNODE; i++)
		{
			gsl_matrix_set(B, 0, 2 * i, gsl_matrix_get(dNdx, i, 0));
			gsl_matrix_set(B, 3, 2 * i, gsl_matrix_get(dNdx, i, 1));
		}

		for (int i = 1; i <= NNODE; i++)
		{
			gsl_matrix_set(B, 1, (i * 2) - 1, gsl_matrix_get(dNdx, i - 1, 1));
			gsl_matrix_set(B, 3, (i * 2) - 1, gsl_matrix_get(dNdx, i - 1, 0));
		}

		//writeGSLMatrix(B, "BMatrix_"+std::to_string(kint)+".csv");


		//Update RHS
		strain = multiply_GSL(B, Uelm);

		//writeGSLMatrix(strain, "strain_" + std::to_string(kint) + ".txt");


		stress = multiply_GSL(D, strain);

		//writeGSLMatrix(stress, "stress_" + std::to_string(kint) + ".txt");


		tempVec = multiply_GSL(transpose_GSL(B), stress);
		gsl_vector_scale(tempVec, gsl_vector_get(w, kint)*determinant*-1.0);
		gsl_vector_add(RHS, tempVec);


		//Update element stiffness
		tempMat = multiply_GSL(transpose_GSL(B), multiply_GSL(D, B));
		gsl_matrix_scale(tempMat, gsl_vector_get(w, kint)*determinant);
		gsl_matrix_add(Kelm, tempMat);

		//write stress to element stress vector
		for (int i = 0; i < stress->size; i++)
		{
			//stressVector.push_back(gsl_vector_get(stress, i));
			fields["S"][sVarsStart + sVarIter] = gsl_vector_get(stress, i);
			sVarIter++;
			//std::cout << sVarIter<<" , " << sVarsStart <<std::endl;
		}

	}

	//writeGSLMatrix(Kelm, "Kelement.txt");
	//writeGSLMatrix(RHS, "RHSelement.txt");

	//Write the results to map.

	gsl_matrix_free(D);
	gsl_vector_free(tempVec);
	gsl_vector_free(strain);
	gsl_vector_free(stress);
	gsl_matrix_free(tempMat);
	gsl_matrix_free(xi);
	gsl_vector_free(w);
	gsl_vector_free(N);
	gsl_matrix_free(dNdxi);

	return true;
}

Element_CPE3::~Element_CPE3()
{
}