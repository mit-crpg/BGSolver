// mexUtilities.h
/* Header file for a collection of utility functions for MEX file development.
 *
 * Package:    BGSolver v1.03
 * Subpackage: Utilities
 * Date:       November 8, 2012
 * Author:     Eugeny Sosnovsky
 *             esos@mit.edu
 */

#ifndef MEXUTILITIES_H
#define MEXUTILITIES_H

#include "matrix.h"
#include "mex.h"

/*
 * MATLAB Data Structure Handling Functions
 */

// Array size scalar structure field retrieval function
// REQUIREMENTS
// - The data structure MATLAB array must be a 1x1 array.
mwSize mxGetDStructFieldSize(const mxArray* pmDStruct, const char* fldName);

// Array of array sizes structure field retrieval function
// REQUIREMENTS
// - The data structure MATLAB array must be a 1x1 array.
// - The array of array sizes field must be of an unsigned integer MATLAB type.
mwSize* mxGetDStructFieldSizes(const mxArray* pmDStruct, const char* fldName);

// Double scalar structure field retrieval function
// REQUIREMENTS
// - The data structure MATLAB array must be a 1x1 array.
double mxGetDStructFieldDouble(const mxArray* pmDStruct, const char* fldName);

// Integer scalar structure field retrieval function
// REQUIREMENTS
// - The data structure MATLAB array must be a 1x1 array.
int mxGetDStructFieldInt(const mxArray* pmDStruct, const char* fldName);

// Logical scalar structure field retrieval function
// REQUIREMENTS
// - The data structure MATLAB array must be a 1x1 array.
// DETAILS
// - Returns the specified default value if the field does not exist.
mxLogical mxGetDStructFieldLogical(const mxArray* pmDStruct,
      const char* fldName, const mxLogical defValue);

// Logical array pointer structure field retrieval function
// REQUIREMENTS
// - The data structure MATLAB array must be a 1x1 array.
// DETAILS
// - Returns a NULL pointer if the field does not exist.
mxLogical* mxGetDStructFieldLogicalPtr(const mxArray* pmDStruct,
      const char* fldName);

// Array of array indices structure field retrieval function
// REQUIREMENTS
// - The data structure MATLAB array must be a 1x1 array.
// - The indices field must be of an unsigned integer MATLAB type.
// - The indices must be one-based.
// DETAILS
// - If the index array is empty, a null pointer is returned.
// - The indices are returned as zero-based.
mwIndex* mxGetDStructFieldIndices(const mxArray* pmDStruct,
      const char* fldName);

// Array of array indices structure field element retrieval function
// REQUIREMENTS
// - The data structure MATLAB array must be a 1x1 array.
// - The indices element of the cell array field must be of an unsigned integer
//   MATLAB type.
// - The indices must be one-based, the element index must be zero-based.
// DETAILS
// - If the index array is empty, a null pointer is returned.
// - The indices are returned as zero-based.
mwIndex* mxGetDStructFieldElementIndices(const mxArray* pmDStruct,
      const char* fldName, const mwIndex elemID);

// Array size from array of indices structure field retrieval function
// REQUIREMENTS
// - The data structure MATLAB array must be a 1x1 array.
mwSize mxGetDStructFieldSizeFromIndices(const mxArray* pmDStruct,
      const char* fldName);

// Array size from array of indices structure field element retrieval function
// REQUIREMENTS
// - The data structure MATLAB array must be a 1x1 array.
// - The field array containing the indices array must be a cell array.
mwSize mxGetDStructFieldElementSizeFromIndices(const mxArray* pmDStruct,
      const char* fldName, const mwIndex elemID);

// Array of array sizes from array of indices structure field retrieval function
// REQUIREMENTS
// - The data structure MATLAB array must be a 1x1 array.
// - The field array containing the indices must be a cell array.
mwSize* mxGetDStructFieldSizesFromIndicesArray(const mxArray* pmDStruct,
      const char* fldName);

// Cell array element pointer structure field element retrieval function
// REQUIREMENTS
// - The data structure MATLAB array must be a 1x1 array.
// - The field array containing the object must be a cell array.
mxArray* mxGetDStructCellElement(const mxArray* pmDStruct, const char* fldName,
      const mwIndex elemID);

/*
 * MATLAB Array Unitialized Preallocation Functions
 *
 * The idea for this came out of this document:
 * http://www.mathworks.com/matlabcentral/fileexchange/
 *    27151-writing-matlab-cmex-code
 */

// Double MATLAB Vertical Vector Uninitialized Preallocation Function
mxArray* mxPreallocDoubleVector(const mwSize nElems);

// Double MATLAB Full Matrix Unitialized Preallocation Function
mxArray* mxPreallocDoubleMatrix(const mwSize nRows, const mwSize nCols);

/*
 * Size and Index Array Conversion Functions
 */

// Array of array sizes type conversion function
// REQUIREMENTS
// - The array of array sizes must be of an unsigned integer MATLAB type.
int mxConvertSizeArray(const mxArray* inpArray, mwSize nSizes,
      mwSize* outArray);

// Array of array indices type conversion function
// REQUIREMENTS
// - The array of array indices must be of an unsigned integer MATLAB type.
// - The indices must be one-based.
// DETAILS
// - The indices are returned as zero-based.
int mxConvertIndexArray(const mxArray* inpArray, mwSize nIndices,
      mwIndex* outArray);

#endif // MEXUTILITIES_H