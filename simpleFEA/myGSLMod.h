#pragma once
#include <gsl\gsl_math.h>
#include <gsl\gsl_vector.h>
#include <gsl\gsl_test.h>
#include <gsl\gsl_ieee_utils.h>
#include <gsl\gsl_sf_bessel.h>
#include <gsl\gsl_matrix.h>
#include <gsl\gsl_splinalg.h>
#include <iostream>


#ifndef __KARTHIK_GSLMods_H__
#define __KARTHIK_GSLMods_H__


gsl_matrix * transpose_GSL(const gsl_matrix * A);
gsl_matrix * multiply_GSL(const gsl_matrix * A, const gsl_matrix * B);
gsl_vector * multiply_GSL(const gsl_matrix * A, const gsl_vector * B);
bool writeGSLMatrix(const gsl_matrix * A, std::string fileName);
bool writeGSLMatrix(const gsl_matrix_int * A, std::string fileName);
bool writeGSLMatrix(const gsl_vector * A, std::string fileName);
bool writeGSLMatrix(const gsl_vector_int * A, std::string fileName);
bool writeGSLMatrix(const gsl_spmatrix * A, std::string fileName);





#endif