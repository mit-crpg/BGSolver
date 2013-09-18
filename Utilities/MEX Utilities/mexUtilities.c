// mexUtilities.c
/* Source file for a collection of utility functions for MEX file development.
 *
 * Package:    BGSolver v1.03
 * Subpackage: Utilities
 * Date:       November 8, 2012
 * Author:     Eugeny Sosnovsky
 *             esos@mit.edu
 */

#include "mexUtilities.h"
#include "matrix.h"
#include "mex.h"

/*
 * MATLAB Data Structure Handling Functions
 */

// Array size scalar structure field retrieval function
mwSize mxGetDStructFieldSize(const mxArray* pmDStruct, const char* fldName)
{
   return (mwSize)mxGetScalar(mxGetField(pmDStruct,0,fldName));
}

// Array of array sizes structure field retrieval function
mwSize* mxGetDStructFieldSizes(const mxArray* pmDStruct, const char* fldName)
{
   // Declarations
   mwSize nSizes;
   mwSize* arrSizes;
   const mxArray* fld;
   
   // Extracting the field
   fld = mxGetField(pmDStruct,0,fldName);
   
   // Sizing the field
   nSizes = (mwSize)mxGetNumberOfElements(fld);
   
   // Preallocating array of array sizes
   arrSizes = (mwSize*)mxMalloc(sizeof(mwSize)*nSizes);
   
   // Converting data type
   mxConvertSizeArray(fld,nSizes,arrSizes);
   
   return arrSizes;
}

// Double scalar structure field retrieval function
double mxGetDStructFieldDouble(const mxArray* pmDStruct, const char* fldName)
{
   return mxGetScalar(mxGetField(pmDStruct,0,fldName));
}

// Integer scalar structure field retrieval function
int mxGetDStructFieldInt(const mxArray* pmDStruct, const char* fldName)
{
   // Declarations
   mxArray* pmFld;
   mxClassID fldClassID;
   int fldVal;
   
   pmFld = mxGetField(pmDStruct,0,fldName);
   fldClassID = mxGetClassID(pmFld);
   
   switch (fldClassID)
   {
      case mxINT8_CLASS:
      {
         fldVal = (int)(*(int8_T*)mxGetData(pmFld));
         break;
      }
      case mxINT16_CLASS:
      {
         fldVal = (int)(*(int16_T*)mxGetData(pmFld));
         break;
      }
      case mxINT32_CLASS:
      {
         fldVal = (int)(*(int32_T*)mxGetData(pmFld));
         break;
      }
      case mxINT64_CLASS:
      {
         fldVal = (int)(*(int64_T*)mxGetData(pmFld));
         break;
      }
      default:
      {
         mexErrMsgTxt("Error! Signed integer field type unknown!!!");
      }
   }
   
   return fldVal;
}

// Logical scalar structure field retrieval function
mxLogical mxGetDStructFieldLogical(const mxArray* pmDStruct,
      const char* fldName, const mxLogical defValue)
{
   // Declarations
   mxLogical* bPtr;
   mxLogical fldValue;
   
   bPtr = mxGetDStructFieldLogicalPtr(pmDStruct,fldName);
   
   if (bPtr != (mxLogical*)0)
   {
      fldValue = *bPtr;
   }
   else
   {
      fldValue = defValue;
   }
   
   return fldValue;
}

// Logical array pointer structure field retrieval function
mxLogical* mxGetDStructFieldLogicalPtr(const mxArray* pmDStruct,
      const char* fldName)
{
   // Declarations
   mxArray* fld;
   mxLogical* bPtr;
   
   fld = mxGetField(pmDStruct,0,fldName);
   if (fld != (mxArray*)0)
   {
      bPtr = mxGetLogicals(fld);
   }
   else
   {
      bPtr = (mxLogical*)0;
   }
   
   return bPtr;
}

// Array of array indices structure field retrieval function
mwIndex* mxGetDStructFieldIndices(const mxArray* pmDStruct,
      const char* fldName)
{
   // Declarations
   mwSize nIndices;
   mwIndex* arrIndices;
   const mxArray* fld;
   
   // Extracting the field
   fld = mxGetField(pmDStruct,0,fldName);
   
   // Sizing the field
   nIndices = (mwSize)mxGetNumberOfElements(fld);
   
   // Checking for empty indices
   if (nIndices == 0)
   {
      arrIndices = (mwIndex*)0;
   }
   else
   {
      // Preallocating array of array indices
      arrIndices = (mwIndex*)mxMalloc(sizeof(mwIndex)*nIndices);
      
      // Converting data type
      mxConvertIndexArray(fld,nIndices,arrIndices);
   }
   
   return arrIndices;
}

// Array of array indices structure field element retrieval function
mwIndex* mxGetDStructFieldElementIndices(const mxArray* pmDStruct,
      const char* fldName, const mwIndex elemID)
{
   // Declarations
   mwSize nIndices;
   mwIndex* arrIndices;
   const mxArray* fldCellArray;
   const mxArray* fldIndices;
   
   // Extracting the indices
   fldCellArray = mxGetField(pmDStruct,0,fldName);
   fldIndices = mxGetCell(fldCellArray,elemID);
   
   // Sizing the field
   nIndices = (mwSize)mxGetNumberOfElements(fldIndices);
   
   // Checking for empty indices
   if (nIndices == 0)
   {
      arrIndices = 0;
   }
   else
   {
      // Preallocating array of array indices
      arrIndices = (mwIndex*)mxMalloc(sizeof(mwIndex)*nIndices);
      
      // Converting data type
      mxConvertIndexArray(fldIndices,nIndices,arrIndices);
   }
   
   return arrIndices;
}

// Array size from array of indices structure field retrieval function
mwSize mxGetDStructFieldSizeFromIndices(const mxArray* pmDStruct,
      const char* fldName)
{
   return (mwSize)mxGetNumberOfElements(mxGetField(pmDStruct,0,fldName));
}

// Array size from array of indices structure field element retrieval function
mwSize mxGetDStructFieldElementSizeFromIndices(const mxArray* pmDStruct,
      const char* fldName, const mwIndex elemID)
{
   // Declarations
   const mxArray* fldCellArray;
   const mxArray* fldIndices;
   mwSize nIndices;
   
   // Extracting the indices
   fldCellArray = mxGetField(pmDStruct,0,fldName);
   fldIndices = mxGetCell(fldCellArray,elemID);
   
   // Extracting the number of indices
   nIndices = (mwSize)mxGetNumberOfElements(fldIndices);
   
   return nIndices;
}

// Array of array sizes from array of indices structure field retrieval function
mwSize* mxGetDStructFieldSizesFromIndicesArray(const mxArray* pmDStruct,
      const char* fldName)
{
   // Declarations
   mwSize nIndexArrays;
   mwSize* nIndicesArr;
   mwIndex i;
   
   // Extracting the sizes
   nIndexArrays = mxGetDStructFieldSizeFromIndices(pmDStruct,fldName);
   nIndicesArr = (mwSize*)mxMalloc(sizeof(mwSize)*nIndexArrays);
   for (i = 0; i < nIndexArrays; i++)
   {
      nIndicesArr[i] = \
      	mxGetDStructFieldElementSizeFromIndices(pmDStruct,fldName,i);
   }
   
   return nIndicesArr;
}

// Cell array element pointer structure field element retrieval function
mxArray* mxGetDStructCellElement(const mxArray* pmDStruct, const char* fldName,
      const mwIndex elemID)
{
   return mxGetCell(mxGetField(pmDStruct,0,fldName),elemID);
}

/*
 * MATLAB Array Unitialized Preallocation Functions
 */

// Double MATLAB Vertical Vector Unitialized Preallocation Function
mxArray* mxPreallocDoubleVector(const mwSize nElems)
{
   // Declaring output
   mxArray* pmV;
   
   // Preallocating vector
   pmV = mxCreateDoubleMatrix(0,0,mxREAL);
   mxSetM(pmV,nElems);
   mxSetN(pmV,1);
   mxSetData(pmV,mxMalloc((mwSize)sizeof(double)*nElems*1));
   
   return pmV;
}

// Double MATLAB Full Matrix Unitialized Preallocation Function
mxArray* mxPreallocDoubleMatrix(const mwSize nRows, const mwSize nCols)
{
   // Declaring output
   mxArray* pmV;
   
   // Preallocating vector
   pmV = mxCreateDoubleMatrix(0,0,mxREAL);
   mxSetM(pmV,nRows);
   mxSetN(pmV,nCols);
   mxSetData(pmV,mxMalloc((mwSize)sizeof(double)*nRows*nCols));
   
   return pmV;
}

/*
 * Size and Index Array Conversion Functions
 */

// Array of array sizes type conversion function
int mxConvertSizeArray(const mxArray* inpArray, mwSize nSizes, mwSize* outArray)
{
   // Declarations
   mxClassID inpClassID;
   mwIndex sid;
   
   // Examining the input
   inpClassID = mxGetClassID(inpArray);
   
   // Extracting correctly converted data pointer
   switch (inpClassID)
   {
      case mxUINT8_CLASS:
      {
         uint8_T* pArrSizes;
         pArrSizes = (uint8_T*)mxGetData(inpArray);
         for (sid = 0; sid < nSizes; sid++)
         {
            outArray[sid] = (mwSize)pArrSizes[sid];
         }
         break;
      }
      case mxUINT16_CLASS:
      {
         uint16_T* pArrSizes;
         pArrSizes = (uint16_T*)mxGetData(inpArray);
         for (sid = 0; sid < nSizes; sid++)
         {
            outArray[sid] = (mwSize)pArrSizes[sid];
         }
         break;
      }
      case mxUINT32_CLASS:
      {
         uint32_T* pArrSizes;
         pArrSizes = (uint32_T*)mxGetData(inpArray);
         for (sid = 0; sid < nSizes; sid++)
         {
            outArray[sid] = (mwSize)pArrSizes[sid];
         }
         break;
      }
      case mxUINT64_CLASS:
      {
         uint64_T* pArrSizes;
         pArrSizes = (uint64_T*)mxGetData(inpArray);
         for (sid = 0; sid < nSizes; sid++)
         {
            outArray[sid] = (mwSize)pArrSizes[sid];
         }
         break;
      }
      default:
      {
         mexErrMsgTxt("Error! Array of array sizes field type unknown!!!");
      }
   }
   
   return 0;
}

// Array of array indices type conversion function
int mxConvertIndexArray(const mxArray* inpArray, mwSize nIndices,
      mwIndex* outArray)
{
   // Declarations
   mxClassID inpClassID;
   mwIndex iid;
   
   // Examining the input
   inpClassID = mxGetClassID(inpArray);
   
   // Extracting correctly converted data pointer
   switch (inpClassID)
   {
      case mxUINT8_CLASS:
      {
         uint8_T* pArrIndices;
         pArrIndices = (uint8_T*)mxGetData(inpArray);
         for (iid = 0; iid < nIndices; iid++)
         {
            outArray[iid] = (mwIndex)pArrIndices[iid]-1;
         }
         break;
      }
      case mxUINT16_CLASS:
      {
         uint16_T* pArrIndices;
         pArrIndices = (uint16_T*)mxGetData(inpArray);
         for (iid = 0; iid < nIndices; iid++)
         {
            outArray[iid] = (mwIndex)pArrIndices[iid]-1;
         }
         break;
      }
      case mxUINT32_CLASS:
      {
         uint32_T* pArrIndices;
         pArrIndices = (uint32_T*)mxGetData(inpArray);
         for (iid = 0; iid < nIndices; iid++)
         {
            outArray[iid] = (mwIndex)pArrIndices[iid]-1;
         }
         break;
      }
      case mxUINT64_CLASS:
      {
         uint64_T* pArrIndices;
         pArrIndices = (uint64_T*)mxGetData(inpArray);
         for (iid = 0; iid < nIndices; iid++)
         {
            outArray[iid] = (mwIndex)pArrIndices[iid]-1;
         }
         break;
      }
      default:
      {
         mexErrMsgTxt("Error! Array of array sizes field type unknown!!!");
      }
   }
   
   return 0;
}