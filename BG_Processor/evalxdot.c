// evalxdot.c
/* State derivative vector evaluation mexFunction.
 *
 * xdot = evalxdot(t,x,xdotsys)
 *
 * This function evaluates the state derivative vector at a given time, with a
 * given (possibly vectorized) x vector, based on a give state derivative vector
 * system data structure.
 *
 * Package:    BGSolver v1.03
 * Subpackage: BG_Processor
 * Date:       November 8, 2012
 * Author:     Eugeny Sosnovsky
 *             esos@mit.edu
 */

#include "mexUtilities.h"
#include "matrix.h"
#include "mex.h"
#include <string.h>
#include <omp.h>

/*
 * Forward declarations. Preallocation, filling and indexing functions.
 */

// Number of input arrays for a given NV type calculation function
// (forward declaration)
mwSize getNVInpNum(const mwIndex nvtype);

// State and BV combination into the C-array function (forward declaration)
double* combineXandBVs(mxArray* pmx, const mwSize nVect);

// Array-filling function with source-based indices (forward declaration)
void fillArrayWithSrcIDs(double* pdest, const double* psrc, const mwSize nElems,
      const mwIndex* srcIDs, const mwSize nStep, const mwSize nVect);

// Array-filling function with destination-based indices (forward declaration)
void fillArrayWithDestIDs(double* pdest, const double* psrc,
      const mwSize nElems, const mwIndex* destIDs, const mwSize nStep,
      const mwSize nVect);

// Array-filling single vector input function with destination-based indices
// (forward declaration)
void fillArrayWithDestIDsScalInp(double* pdest, const double* psrc,
      const mwSize nElems, const mwIndex* destIDs, const mwSize nStep,
      const mwSize nVect);

/*
 * Forward declarations. NV Evaluation functions.
 */

// Numeric function types 1 and 3 evaluation function (forward declaration)
void evalNfh13(double* pNVrfx, const mwSize nElems, const mwIndex* tnidsrfx,
      const double* pC, const mwSize nVect, mxArray*** pmInputArrays,
      const mwSize* inpSizes, const mwIndex** inpCIDs);

// Numeric function type 2 evaluation function (forward declaration)
void evalNfh2(double* pNVrfx, const mwSize nElems, const mwIndex* tnidsrfx,
      const mwSize nVect, mxArray*** pmInputArrays);

// Numeric function type 4 evaluation function (forward declaration)
void evalNfh4(double* pNVrfx, const mwSize nElems, const mwIndex* tnidsrfx,
      const double* pC, const mwSize nVect, mxArray*** pmInputArrays,
      const mwSize* inpSizes, const mwIndex** inpCIDs);

// Numeric function type 5 evaluation function (forward declaration)
void evalNfh5(double* pNVrfx, const mwSize nElems, const mwIndex* tnidsrfx,
      const double* pC, const mwSize nVect, mxArray*** pmInputArrays,
      const mwSize* inpSizes, const mwIndex** inpCIDs, const mwSize* minpSizes,
      const mwIndex** minpCIDs);

// Numeric function type 6 evaluation function (forward declaration)
void evalNfh6(double* pNVrfx, const mwSize nElems, const mwIndex* tnidsrfx,
      const double* pC, const mwSize nVect, mxArray*** pmInputArrays,
      const mwSize* inpSizes, const mwIndex** inpCIDs);

// Numeric function type 7 evaluation function (forward declaration)
void evalNfh7(double* pNVrfx, const mwSize nElems, const mwIndex* tnidsrfx,
      const double* pC, const mwSize nVect, mxArray*** pmInputArrays,
      const mwSize* inpSizes, const mwIndex** inpCIDs, const mwSize* minpSizes,
      const mwIndex** minpCIDs);

/*
 * Forward declarations. Cleanup function.
 */

// Cleanup function (forward declaration)
static void evalxdotCleanup(void);

/* 
 * Static declarations to be used in subsequent evaluations.
 */
// Parallel-to-serial loop size barrier
static const int ser2par = 10;
// 1st Evaluation Status
static mxLogical ifFirstEval = 1;
// Vectorization Requirement Status
static mxLogical ifVectorReq;
// Layer 1 BV vector scalar status
static mxLogical ifBV1Scal;
// Array sizes
static mwSize nX;
static mwSize bvnumrfx, nvnumrfx;
static mwSize lnumrfx;
static mwSize nC;
// Static time array
static mxArray* st_pmt;
static double* st_pt;
// Intermediate BV sizes
static mwSize* nBVs; // nBVs[lid]
// Intermediate BV function handles
static const mxArray** pmBVfhs; // pmBVfhs[lid]
// Intermediate evaluated BV indices
static const mwIndex** pBVlbvids; // pBVlbvids[lid]
// Intermediate BV input statuses
static const mxLogical* pIfBVInpT; // pIfBVInpT[lid]
static const mxLogical* pIfBVInpX; // pIfBVInpX[lid]
static const mxLogical* pIfBVInpN; // pIfBVInpN[lid]
static mxLogical* pIfBVInpIDt; // pIfBVInpIDt[bviid]
// Intermediate BV input sizes
static mwSize nBVInps;
static mwSize* pBVInpSiz; // pBVInpSiz[lid]
static mwSize* pBVInpXSiz; // pBVInpXSiz[lid]
static mwSize* pBVInpNSiz; // pBVInpNSiz[lid]
// Intermediate BV input indices
static mwIndex** pBVInpXIDs; // pBVInpXIDs[lid][bvixid]
static mwIndex** pBVInpNIDs; // pBVInpNIDs[lid][bvinid]
// Intermediate BV inputs
static mxArray** pmBVInps; // pmBVInps[bviid] = mxArray* (input array)
static double** pBVInps; // pBVInps[bviid] = double* (input array)
// Intermediate MATLAB BV I/O arrays
static const mxArray*** pmBVrhs; // pmBVrhs[lid][bvliid] = mxArray*
                                 // (input array)
// Intermediate NV sizes
static mwSize nNVtypes;
static mwSize* nNVs; // nNVs[nvtyp+7*lid]
// Intermediate evaluated NV indices
static const mwIndex** pNVltnids; // pNVltnids[nvtyp+7*lid][ltnid]
// Intermediate NV input sizes
static const mwSize** pNVInpSiz; // pNVInpSiz[nvtyp+7*lid][ltnid]
static const mwSize** pNVMInpSiz; // pNVMInpSiz[nvtyp+7*lid][ltnid]
// Intermediate NV input indices
static const mwIndex*** pNVInpIDs; // pNVInpIDs[nvtyp+7*lid][ltnid][nvicid]
static const mwIndex*** pNVMInpIDs; // pNVMInpIDs[nvtyp+7*lid][ltnid][nvicid]
// Intermediate NV inputs
static mxArray**** pmNVInps; // pmNVInps[nvtyp+7*lid][ltnid][nviid] = mxArray*
                             // (input array)
// Final xdot indices
static mwIndex* pXDotIDs; // pXDotIDs[xid]

// Gateway function
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
   // Evaluation-local declarations
   // Array sizes
   mwSize nVect; // Number of input x-vectors (1 if not vectorized)
   // Indices
   mwIndex lid; // Layer
   mwIndex bviid; // BV input
   mwIndex lbviid; // First layer-specific BV input 
   mwIndex nvtyp; // NV type
   mwIndex ltnid; // Layer- and type-local NV
   mwIndex nidrfx; // NV (RFX)
   mwIndex nviid; // NV Input
   // Input data arrays
   const mxArray* pmt; // Time (mxArray)
   const mxArray* pmx; // State (mxArray)
   const double* px; // State (double)
   const mxArray* pmxdotsys; // xdot vector system data structure (mxArray)
   // Intermediate (RFX) data arrays
   double* pBVrfx; // BVs (RFX) (double) (within C-array, vector 0)
   double* pNVrfx; // NVs (RFX) (double)
   double* pC; // Combined C-array (double)
   // Intermediate MATLAB BV I/O arrays
   mxArray* pmBVlhs[1]; // pmBVlhs[0] = mxArray* (BV output array)
   double* pBVlhs; // BV output array (double)
   // Output data arrays
   mxArray* pmxdot; // xdot (mxArray)
   double* pxdot; // xdot (double)
   
   // Input check (to prevent accidental execution, remove before release)
   if (nrhs != 3)
   {
      mexErrMsgTxt("nargin wrong, evalxdotMEX needs the same inputs as \
         evalxdot!");
   }
   
   // Retrieving inputs
   pmt = prhs[0];
   pmx = prhs[1];
   pmxdotsys = prhs[2];
   
   // Sizing the problem (evaluation-variable)
   nVect = (mwSize)mxGetN(pmx);
   
   // 1st evaluation procedures
   if (ifFirstEval)
   {
      // Declarations (temporary)
      const mxArray** pmNVInpTrackers;
      int bv1vectstat;
      
      // Extracting vectorization requirement
      ifVectorReq = mxGetDStructFieldLogical(pmxdotsys,"ifvectreq",1);
      
      // Extracting BV layer 1 scalar status
      bv1vectstat = mxGetDStructFieldInt(pmxdotsys,"bv1vectstat");
      if (bv1vectstat == 2)
      {
         ifBV1Scal = (mxLogical)1;
      }
      else
      {
         ifBV1Scal = (mxLogical)0;
      }
      
      // Sizing the problem
      nX = mxGetDStructFieldSize(pmxdotsys,"xnum");
      bvnumrfx = mxGetDStructFieldSize(pmxdotsys,"bvnumrfx");
      nvnumrfx = mxGetDStructFieldSize(pmxdotsys,"nvnumrfx");
      lnumrfx = mxGetDStructFieldSize(pmxdotsys,"lnumrfx");
      nC = nX + bvnumrfx;
      nBVs = mxGetDStructFieldSizesFromIndicesArray(pmxdotsys,"lbvidsrfx");
      mexMakeMemoryPersistent(nBVs);
      nNVtypes = 7*(lnumrfx-1);
      
      // Initializing static time array
      st_pmt = mxPreallocDoubleVector(1);
      mexMakeArrayPersistent(st_pmt);
      st_pt = mxGetPr(st_pmt);
      
      // Extracting intermediate BV input statuses
      pIfBVInpT = mxGetDStructFieldLogicalPtr(pmxdotsys,"bvifit");
      pIfBVInpX = mxGetDStructFieldLogicalPtr(pmxdotsys,"bvifix");
      pIfBVInpN = mxGetDStructFieldLogicalPtr(pmxdotsys,"bvifin");
      
      // Extracting intermediate BV input sizes
      pBVInpSiz = mxGetDStructFieldSizes(pmxdotsys,"bvInpSiz");
      mexMakeMemoryPersistent(pBVInpSiz);
      nBVInps = 0;
      for (lid = 0; lid < lnumrfx; lid++)
      {
         nBVInps += pBVInpSiz[lid];
      }
      
      // Preallocating intermediate BV input statuses, BVIDs, BV inputs,
      // input sizes, input IDs, BV function handles and BV I/O arrays
      pIfBVInpIDt = (mxLogical*)mxCalloc(nBVInps,sizeof(mxLogical));
      mexMakeMemoryPersistent(pIfBVInpIDt);
      pBVlbvids = (mwIndex**)mxCalloc(lnumrfx,sizeof(mwIndex*));
      mexMakeMemoryPersistent(pBVlbvids);
      pmBVInps = (mxArray**)mxCalloc(nBVInps,sizeof(mxArray*));
      mexMakeMemoryPersistent(pmBVInps);
      pBVInps = (double**)mxCalloc(nBVInps,sizeof(double*));
      mexMakeMemoryPersistent(pBVInps);
      pBVInpXSiz = (mwSize*)mxCalloc(lnumrfx,sizeof(mwSize));
      mexMakeMemoryPersistent(pBVInpXSiz);
      pBVInpNSiz = (mwSize*)mxCalloc(lnumrfx,sizeof(mwSize));
      mexMakeMemoryPersistent(pBVInpNSiz);
      pBVInpXIDs = (mwIndex**)mxCalloc(lnumrfx,sizeof(mwIndex*));
      mexMakeMemoryPersistent(pBVInpXIDs);
      pBVInpNIDs = (mwIndex**)mxCalloc(lnumrfx,sizeof(mwIndex*));
      mexMakeMemoryPersistent(pBVInpNIDs);
      pmBVfhs = (mxArray**)mxMalloc(lnumrfx*sizeof(mxArray*));
      mexMakeMemoryPersistent(pmBVfhs);
      pmBVrhs = (mxArray***)mxMalloc(lnumrfx*sizeof(mxArray**));
      mexMakeMemoryPersistent(pmBVrhs);
      
      // Preallocating intermediate NV sizes
      nNVs = (mwSize*)mxMalloc(nNVtypes*sizeof(mwSize));
      mexMakeMemoryPersistent(nNVs);
      
      // Preallocating intermediate evaluated NIDs, NV input sizes,
      // NV input indices and NV inputs
      pNVltnids = (mwIndex**)mxCalloc(nNVtypes,sizeof(mwIndex*));
      mexMakeMemoryPersistent(pNVltnids);
      pNVInpSiz = (mwSize**)mxCalloc(nNVtypes,sizeof(mwSize*));
      mexMakeMemoryPersistent(pNVInpSiz);
      pNVMInpSiz = (mwSize**)mxCalloc(nNVtypes,sizeof(mwSize*));
      mexMakeMemoryPersistent(pNVMInpSiz);
      pmNVInpTrackers = (mxArray**)mxCalloc(nNVtypes,sizeof(mxArray*));
      pNVInpIDs = (mwIndex***)mxCalloc(nNVtypes,sizeof(mwIndex**));
      mexMakeMemoryPersistent(pNVInpIDs);
      pNVMInpIDs = (mwIndex***)mxCalloc(nNVtypes,sizeof(mwIndex**));
      mexMakeMemoryPersistent(pNVMInpIDs);
      pmNVInps = (mxArray****)mxCalloc(nNVtypes,sizeof(mxArray***));
      mexMakeMemoryPersistent(pmNVInps);
      
      // Looping through levels
      bviid = 0; // Initializing BV input counter
      for (lid = 0; lid < lnumrfx; lid++)
      {
         // Preparing level-specific BV evaluations
         if (nBVs[lid] > 0)
         {
            lbviid = bviid;
            if (pIfBVInpT[lid])
            {
               pmBVInps[bviid] = st_pmt;
               pIfBVInpIDt[bviid] = (mxLogical)1;
               
               bviid++;
            }
            if (pIfBVInpX[lid])
            {
               // Extracting X input size
               pBVInpXSiz[lid] = mxGetDStructFieldElementSizeFromIndices(
                  pmxdotsys,"lbvinpxids",lid);
               
               // Extracting X input indices
               pBVInpXIDs[lid] = mxGetDStructFieldElementIndices(pmxdotsys,
                  "lbvinpxids",lid);
               mexMakeMemoryPersistent(pBVInpXIDs[lid]);
               
               // Preallocating X input array, if vectorization not required
               if (!ifVectorReq)
               {
                  pmBVInps[bviid] = mxPreallocDoubleVector(pBVInpXSiz[lid]);
                  mexMakeArrayPersistent(pmBVInps[bviid]);
                  pBVInps[bviid] = mxGetPr(pmBVInps[bviid]);
               }
               
               bviid++;
            }
            if (pIfBVInpN[lid])
            {
               // Extracting NV input size
               pBVInpNSiz[lid] = mxGetDStructFieldElementSizeFromIndices(
                  pmxdotsys,"lbvinpnidsrfx",lid);
               
               // Extracting NV input indices
               pBVInpNIDs[lid] = mxGetDStructFieldElementIndices(
                  pmxdotsys,"lbvinpnidsrfx",lid);
               mexMakeMemoryPersistent(pBVInpNIDs[lid]);
               
               // Preallocating NV input array, if vectorization not required
               if (!ifVectorReq)
               {
                  pmBVInps[bviid] = mxPreallocDoubleVector(pBVInpNSiz[lid]);
                  mexMakeArrayPersistent(pmBVInps[bviid]);
                  pBVInps[bviid] = mxGetPr(pmBVInps[bviid]);
               }
               
               bviid++;
            }
            pBVlbvids[lid] = mxGetDStructFieldElementIndices(pmxdotsys,
               "lbvidsrfx",lid);
            mexMakeMemoryPersistent(pBVlbvids[lid]);
            // Retrieving level-specific BV function handle
            pmBVfhs[lid] = mxGetDStructCellElement(pmxdotsys,"b",lid);
            pmBVrhs[lid] = \
            	(mxArray**)mxMalloc((1+pBVInpSiz[lid])*sizeof(mxArray*));
            mexMakeMemoryPersistent(pmBVrhs[lid]);
            pmBVrhs[lid][0] = pmBVfhs[lid];
            if (!ifVectorReq && pBVInpSiz[lid] > 0)
            {
               memcpy(&(pmBVrhs[lid][1]),&(pmBVInps[lbviid]),
               	pBVInpSiz[lid]*sizeof(mxArray*));
            }
         }
      }
      
      // Preparing level-specific NV evaluations
      for (lid = 0; lid < (lnumrfx-1); lid++)
      {
         // Preparing NV type-specific NV evaluations
         for (nvtyp = 0; nvtyp < 7; nvtyp++)
         {
            // Extracting number of NVs at this lid and type
            nNVs[nvtyp+7*lid] = mxGetDStructFieldElementSizeFromIndices(
               pmxdotsys,"lnvidsrfxByType",nvtyp+7*lid);
            
            // Proceding if NVs are present at this level and type
            if (nNVs[nvtyp+7*lid] > 0)
            {
               // Retrieving level- and type-specific NV indices
               pNVltnids[nvtyp+7*lid] = mxGetDStructFieldElementIndices(
               	pmxdotsys,"lnvidsrfxByType",nvtyp+7*lid);
               mexMakeMemoryPersistent(pNVltnids[nvtyp+7*lid]);
               // Retrieving level- and type-specific NV input objects
               if (nvtyp != 1)
               {
                  pmNVInpTrackers[nvtyp+7*lid] = mxGetDStructCellElement(
                     pmxdotsys,"lnvinpsByType",nvtyp+7*lid);
               }
               // Retrieving level= and type-specific NV input indices
               // Regular inputs
               if (nvtyp == 0 || nvtyp == 4 || nvtyp == 5 || nvtyp == 6)
               {
                  pNVInpIDs[nvtyp+7*lid] = \
                  	(mwIndex**)mxCalloc(nNVs[nvtyp+7*lid],sizeof(mwIndex*));
                  mexMakeMemoryPersistent(pNVInpIDs[nvtyp+7*lid]);
                  pNVInpSiz[nvtyp+7*lid] = \
                     mxGetDStructFieldSizesFromIndicesArray(
                     pmNVInpTrackers[nvtyp+7*lid],"nvinpids");
                  mexMakeMemoryPersistent(pNVInpSiz[nvtyp+7*lid]);
                  for (ltnid = 0; ltnid < nNVs[nvtyp+7*lid]; ltnid++)
                  {
                     pNVInpIDs[nvtyp+7*lid][ltnid] = \
                        mxGetDStructFieldElementIndices(
                        pmNVInpTrackers[nvtyp+7*lid],"nvinpids",ltnid);
                     mexMakeMemoryPersistent(pNVInpIDs[nvtyp+7*lid][ltnid]);
                  }
               }
               // Modulating inputs
               if (nvtyp == 2 || nvtyp == 3 || nvtyp == 4 || nvtyp == 6)
               {
                  pNVMInpIDs[nvtyp+7*lid] = \
                  	(mwIndex**)mxCalloc(nNVs[nvtyp+7*lid],sizeof(mwIndex*));
                  mexMakeMemoryPersistent(pNVMInpIDs[nvtyp+7*lid]);
                  pNVMInpSiz[nvtyp+7*lid] = \
                     mxGetDStructFieldSizesFromIndicesArray(
                     pmNVInpTrackers[nvtyp+7*lid],"nvminpids");
                  mexMakeMemoryPersistent(pNVMInpSiz[nvtyp+7*lid]);
                  for (ltnid = 0; ltnid < nNVs[nvtyp+7*lid]; ltnid++)
                  {
                     pNVMInpIDs[nvtyp+7*lid][ltnid] = \
                     	mxGetDStructFieldElementIndices(
                        pmNVInpTrackers[nvtyp+7*lid],"nvminpids",ltnid);
                     mexMakeMemoryPersistent(pNVMInpIDs[nvtyp+7*lid][ltnid]);
                  }
               }
               
               // Preallocating level- and type-specific NV function handles
               // and inputs
               pmNVInps[nvtyp+7*lid] = \
                  (mxArray***)mxCalloc(nNVs[nvtyp+7*lid],sizeof(mxArray**));
               mexMakeMemoryPersistent(pmNVInps[nvtyp+7*lid]);
               
               // Preallocating arrays of NV input arrays and extracting NV
               // numeric function handles
               for (ltnid = 0; ltnid < nNVs[nvtyp+7*lid]; ltnid++)
               {
                  pmNVInps[nvtyp+7*lid][ltnid] = \
                     (mxArray**)mxCalloc(getNVInpNum(nvtyp)+1,sizeof(mxArray*));
                  mexMakeMemoryPersistent(pmNVInps[nvtyp+7*lid][ltnid]);
                  
                  nidrfx = pNVltnids[nvtyp+7*lid][ltnid];
                  pmNVInps[nvtyp+7*lid][ltnid][0] = \
                     mxGetDStructCellElement(pmxdotsys,"nhs",nidrfx);
               }
               
               // Preallocating NV inputs
               // Initializing input counter
               nviid = 1;
               
               // Regular inputs (if vectorization is not required)
               if (!ifVectorReq && nvtyp == 0 || nvtyp == 4 || nvtyp == 5 ||
                   nvtyp == 6)
               {
                  for (ltnid = 0; ltnid < nNVs[nvtyp+7*lid]; ltnid++)
                  {
                     pmNVInps[nvtyp+7*lid][ltnid][nviid] = \
                     	mxPreallocDoubleVector(pNVInpSiz[nvtyp+7*lid][ltnid]);
                     mexMakeArrayPersistent(
                        pmNVInps[nvtyp+7*lid][ltnid][nviid]);
                  }
                  nviid++;
               }
               
               // Time inputs (whether or not vectorization is required)
               if (nvtyp == 1 || nvtyp == 3 || nvtyp == 5 || nvtyp == 6)
               {
                  for (ltnid = 0; ltnid < nNVs[nvtyp+7*lid]; ltnid++)
                  {
                     pmNVInps[nvtyp+7*lid][ltnid][nviid] = st_pmt;
                  }
                  nviid++;
               }
               
               // Modulating inputs (if vectorization is not required)
               if (!ifVectorReq && nvtyp == 2 || nvtyp == 3 || nvtyp == 4 ||
                   nvtyp == 6)
               {
                  for (ltnid = 0; ltnid < nNVs[nvtyp+7*lid]; ltnid++)
                  {
                     pmNVInps[nvtyp+7*lid][ltnid][nviid] = \
                        mxPreallocDoubleVector(pNVMInpSiz[nvtyp+7*lid][ltnid]);
                     mexMakeArrayPersistent(
                        pmNVInps[nvtyp+7*lid][ltnid][nviid]);
                  }
               }
            }
         }
      }
      
      // Extracting the final xdot indices
      pXDotIDs = mxGetDStructFieldIndices(pmxdotsys,"xdotIDsrfx");
      mexMakeMemoryPersistent(pXDotIDs);
      
      // Releasing temporary dynamic memory
      mxFree(pmNVInpTrackers);
      
      // Flip the flag
      ifFirstEval = 0;
   }
   
   // Registering cleanup function
   mexAtExit(&evalxdotCleanup);
   
   // Reallocating X and BVs (RFX) into contiguous memory
   pC = combineXandBVs(pmx,nVect);
   px = mxGetPr(pmx);
   pBVrfx = &(pC[nX]);
   
   // Preallocating intermediate NV (RFX) array
   pNVrfx = (double*)mxMalloc(nvnumrfx*nVect*sizeof(double));
   
   // Updating static time value
   st_pt[0] = mxGetScalar(pmt);
   
   // Initializing counters
   bviid = 0;
   
   // Looping through layers and evaluating bond and numeric variables
   for (lid = 0; lid < lnumrfx; lid++)
   {
      // Evaluating layer-specific BVs
      if (nBVs[lid] > 0)
      {
         lbviid = bviid;
         if (pIfBVInpT[lid])
         {
            bviid++;
         }
         if (pBVInpXSiz[lid] > 0)
         {
            // Preallocating BV input array, if vectorization required
            // (otherwise it's persistently preallocated)
            if (ifVectorReq)
            {
               pmBVInps[bviid] = mxPreallocDoubleMatrix(pBVInpXSiz[lid],nVect);
               pBVInps[bviid] = mxGetPr(pmBVInps[bviid]);
            }
            // Filling BV X-input array
            fillArrayWithSrcIDs(pBVInps[bviid],pC,pBVInpXSiz[lid],
               pBVInpXIDs[lid],nC,nVect);
            // Incrementing counter
            bviid++;
         }
         if (pBVInpNSiz[lid] > 0)
         {
            // Preallocating NV input array, if vectorization required
            // (otherwise it's persistently preallocated)
            if (ifVectorReq)
            {
               pmBVInps[bviid] = mxPreallocDoubleMatrix(pBVInpNSiz[lid],nVect);
               pBVInps[bviid] = mxGetPr(pmBVInps[bviid]);
            }
            // Filling BV NV-input array
            fillArrayWithSrcIDs(pBVInps[bviid],pNVrfx,pBVInpNSiz[lid],
               pBVInpNIDs[lid],nvnumrfx,nVect);
            // Incrementing counter
            bviid++;
         }
         
         // Updating MATLAB BV I/O arrays, if vectorization required
         // (otherwise they are persistently present)
         if (ifVectorReq && pBVInpSiz[lid] > 0)
         {
            memcpy(&(pmBVrhs[lid][1]),&(pmBVInps[lbviid]),
            	pBVInpSiz[lid]*sizeof(mxArray*));
         }
         
         // Evaluating BVs
         mexCallMATLAB(1,pmBVlhs,1+pBVInpSiz[lid],pmBVrhs[lid],"feval");
         pBVlhs = mxGetPr(pmBVlhs[0]);
         if (lid == 0 && ifBV1Scal)
         {
            fillArrayWithDestIDsScalInp(pBVrfx,pBVlhs,nBVs[lid],pBVlbvids[lid],
               nC,nVect);
         }
         else
         {
            fillArrayWithDestIDs(pBVrfx,pBVlhs,nBVs[lid],pBVlbvids[lid],nC,
            	nVect);
         }
         // Clearing memory
         mxDestroyArray(pmBVlhs[0]);
      }
      
      // Evaluating layer-specific NVs
      if (lid < (lnumrfx-1))
      {
         // NV type 1
         nvtyp = 0;
         if (nNVs[nvtyp+7*lid] > 0)
         {
            evalNfh13(pNVrfx,nNVs[nvtyp+7*lid],pNVltnids[nvtyp+7*lid],pC,nVect,
            	pmNVInps[nvtyp+7*lid],pNVInpSiz[nvtyp+7*lid],
               pNVInpIDs[nvtyp+7*lid]);
         }
         
         // NV type 2
         nvtyp = 1;
         if (nNVs[nvtyp+7*lid] > 0)
         {
            evalNfh2(pNVrfx,nNVs[nvtyp+7*lid],pNVltnids[nvtyp+7*lid],nVect,
               pmNVInps[nvtyp+7*lid]);
         }
         
         // NV type 3
         nvtyp = 2;
         if (nNVs[nvtyp+7*lid] > 0)
         {
            evalNfh13(pNVrfx,nNVs[nvtyp+7*lid],pNVltnids[nvtyp+7*lid],pC,nVect,
            	pmNVInps[nvtyp+7*lid],pNVMInpSiz[nvtyp+7*lid],
               pNVMInpIDs[nvtyp+7*lid]);
         }
         
         // NV type 4
         nvtyp = 3;
         if (nNVs[nvtyp+7*lid] > 0)
         {
            evalNfh4(pNVrfx,nNVs[nvtyp+7*lid],pNVltnids[nvtyp+7*lid],pC,nVect,
               pmNVInps[nvtyp+7*lid],pNVMInpSiz[nvtyp+7*lid],
               pNVMInpIDs[nvtyp+7*lid]);
         }
         
         // NV type 5
         nvtyp = 4;
         if (nNVs[nvtyp+7*lid] > 0)
         {
            evalNfh5(pNVrfx,nNVs[nvtyp+7*lid],pNVltnids[nvtyp+7*lid],pC,nVect,
               pmNVInps[nvtyp+7*lid],pNVInpSiz[nvtyp+7*lid],
               pNVInpIDs[nvtyp+7*lid],pNVMInpSiz[nvtyp+7*lid],
               pNVMInpIDs[nvtyp+7*lid]);
         }
         
         // NV type 6
         nvtyp = 5;
         if (nNVs[nvtyp+7*lid] > 0)
         {
            evalNfh6(pNVrfx,nNVs[nvtyp+7*lid],pNVltnids[nvtyp+7*lid],pC,nVect,
               pmNVInps[nvtyp+7*lid],pNVInpSiz[nvtyp+7*lid],
               pNVInpIDs[nvtyp+7*lid]);
         }
         
         // NV type 7
         nvtyp = 6;
         if (nNVs[nvtyp+7*lid] > 0)
         {
            evalNfh7(pNVrfx,nNVs[nvtyp+7*lid],pNVltnids[nvtyp+7*lid],pC,nVect,
               pmNVInps[nvtyp+7*lid],pNVInpSiz[nvtyp+7*lid],
               pNVInpIDs[nvtyp+7*lid],pNVMInpSiz[nvtyp+7*lid],
               pNVMInpIDs[nvtyp+7*lid]);
         }
      }
   }
   
   // Filling xdot array
   pmxdot = mxPreallocDoubleMatrix(nX,nVect);
   pxdot = mxGetPr(pmxdot);
   fillArrayWithSrcIDs(pxdot,pBVrfx,nX,pXDotIDs,nC,nVect);
   plhs[0] = pmxdot;
   
   // Freeing input and intermediate arrays, if vectorization required
   // (otherwise they are persistently preallocated)
   if (ifVectorReq)
   {
      // BV inputs
      bviid = 0;
      for (lid = 0; lid < lnumrfx; lid++)
      {
         if (nBVs[lid] > 0)
         {
            if (pIfBVInpT[lid])
            {
               bviid++;
            }
            if (pBVInpXSiz[lid] > 0)
            {
               mxDestroyArray(pmBVInps[bviid]);
               bviid++;
            }
            if (pBVInpNSiz[lid] > 0)
            {
               mxDestroyArray(pmBVInps[bviid]);
               bviid++;
            }
         }
      }
      
      // NV inputs
      for (lid = 0; lid < (lnumrfx-1); lid++)
      {
         for (nvtyp = 0; nvtyp < 7; nvtyp++)
         {
            if (nNVs[nvtyp+7*lid] > 0)
            {
               // Initializing counter
               nviid = 1;
               
               // Regular inputs
               if (nvtyp == 0 || nvtyp == 4 || nvtyp == 5 || nvtyp == 6)
               {
                  for (ltnid = 0; ltnid < nNVs[nvtyp+7*lid]; ltnid++)
                  {
                     mxDestroyArray(pmNVInps[nvtyp+7*lid][ltnid][nviid]);
                  }
                  
                  nviid++;
               }
               
               // Time inputs
               if (nvtyp == 1 || nvtyp == 3 || nvtyp == 5 || nvtyp == 6)
               {
                  nviid++;
               }
               
               // Modulating inputs
               if (nvtyp == 2 || nvtyp == 3 || nvtyp == 4 || nvtyp == 6)
               {
                  for (ltnid = 0; ltnid < nNVs[nvtyp+7*lid]; ltnid++)
                  {
                     mxDestroyArray(pmNVInps[nvtyp+7*lid][ltnid][nviid]);
                  }
               }
            }
         }
      }
   }
   // Freeing intermediate (RFX) arrays
   mxFree(pNVrfx);
   mxFree(pC);
   
   return;
}

// Number of input arrays for a given NV type calculation function
mwSize getNVInpNum(const mwIndex nvtype)
{
   // Declarations
   mwSize nInputs;
   
   switch (nvtype)
   {
      case 0:
      case 1:
      case 2:
         nInputs = 1;
         break;
      case 3:
      case 4:
      case 5:
         nInputs = 2;
         break;
      case 6:
         nInputs = 3;
         break;
   }
   
   return nInputs;
}

// State and BV combination into the C-array function
double* combineXandBVs(mxArray* pmx, const mwSize nVect)
{
   // Declarations
   double* pC; // Contiguous combined memory data vector
   double* px; // Pointer to x (double)
   mwIndex vid;
   mwIndex cid;
   mwIndex xid;
   
   // Retrieving X data pointer
   px = mxGetPr(pmx);
   
   // Making deep copy of X into the C-array
   pC = (double*)mxMalloc(nVect*nC*sizeof(double));
   for (vid = 0; vid < nVect; vid++)
   {
      cid = nC*vid;
      xid = nX*vid;
      
      memcpy(&(pC[cid]),&(px[xid]),nX*sizeof(double));
   }
   
   return pC;
}

// Array-filling function with source-based indices
void fillArrayWithSrcIDs(double* pdest, const double* psrc, const mwSize nElems,
      const mwIndex* srcIDs, const mwSize nStep, const mwSize nVect)
{
   // Declarations
   mwIndex vid; // Vector index
   mwIndex idInSrcID; // Index in source indices
   mwIndex srcID; // Source index
   mwIndex destID; // Destination index
   
   // Looping through input vectors
   for (vid = 0; vid < nVect; vid++)
   {
      // Looping through input variables in this input vector
      for (idInSrcID = 0; idInSrcID < nElems; idInSrcID++)
      {
         srcID = vid*(mwIndex)nStep + srcIDs[idInSrcID];
         destID = vid*(mwIndex)nElems + idInSrcID;
         pdest[destID] = psrc[srcID];
      }
   }
   
   return;
}

// Array-filling function with destination-based indices
void fillArrayWithDestIDs(double* pdest, const double* psrc,
      const mwSize nElems, const mwIndex* destIDs, const mwSize nStep,
      const mwSize nVect)
{
   // Declarations
   mwIndex vid; // Vector index
   mwIndex idInDestID; // Index in destination indices
   mwIndex srcID; // Source index
   mwIndex destID; // Destination index
   
   // Looping through input vectors
   for (vid = 0; vid < nVect; vid++)
   {
      // Looping through input variables in this input vector
      for (idInDestID = 0; idInDestID < nElems; idInDestID++)
      {
         destID = vid*(mwIndex)nStep + destIDs[idInDestID];
         srcID = vid*(mwIndex)nElems + idInDestID;
         pdest[destID] = psrc[srcID];
      }
   }
   
   return;
}

// Array-filling single vector input function with destination-based indices
void fillArrayWithDestIDsScalInp(double* pdest, const double* psrc,
      const mwSize nElems, const mwIndex* destIDs, const mwSize nStep,
      const mwSize nVect)
{
   // Declarations
   mwIndex vid; // Vector index
   mwIndex srcID; // Source index
   mwIndex destID; // Destination index
   
   // Looping through input vectors
   for (vid = 0; vid < nVect; vid++)
   {
      // Looping through input variables in the input vector
      for (srcID = 0; srcID < nElems; srcID++)
      {
         destID = vid*(mwIndex)nStep + destIDs[srcID];
         pdest[destID] = psrc[srcID];
      }
   }
   
   return;
}

// Numeric function types 1 and 3 evaluation function
void evalNfh13(double* pNVrfx, const mwSize nElems, const mwIndex* tnidsrfx,
      const double* pC, const mwSize nVect, mxArray*** pmInputArrays,
      const mwSize* inpSizes, const mwIndex** inpCIDs)
{
   // Declarations
   mwSignedIndex idInGivenNFs; // Index in given numeric functions
   mwIndex idInSpecInput; // Index in a specific input array
   mwIndex inpCID; // Index in C-array
   mwIndex idInNVrfx; // Index in destination NV values
   const mwIndex* inpCIDsForGivenNF; // Array of inpCIDs for a specific NV
   mwSize inpSizForGivenNF; // Number of inputs for a specific NV
   double** pInputs; // pInputs[idInGivenNFs] - Input array (double*)
   mxArray** prhs; // Function handle, Regular/Modulating input
   mxArray* plhs[1]; // NV evaluation (mxArray)
   
   // Checking if input preallocation is required
   if (ifVectorReq)
   {
      // Declarations
      mwIndex idInSpecInput_0; // idInSpecInput, vector 0
      mwIndex inpCID_0; // inpCID, vector 0
      mwIndex idInNVrfx_0; // idInNVrfx, vector 0
      mwIndex vid; // Vector index
      double* pResult; // NV evaluation (double)
      
      // Preallocating
      pInputs = (double**)mxMalloc(nElems*sizeof(double*));
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         // Preallocating the input array
         pmInputArrays[idInGivenNFs][1] = \
            mxPreallocDoubleMatrix(inpSizes[idInGivenNFs],nVect);
         pInputs[idInGivenNFs] = mxGetPr(pmInputArrays[idInGivenNFs][1]);
      }
      
      // Filling
      #pragma omp parallel for private(idInGivenNFs,inpCIDsForGivenNF, \
              inpSizForGivenNF,idInSpecInput_0,inpCID_0,vid,inpCID, \
              idInSpecInput) if(nElems >= ser2par)
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         inpCIDsForGivenNF = inpCIDs[idInGivenNFs];
         inpSizForGivenNF = inpSizes[idInGivenNFs];
         
         for (idInSpecInput_0 = 0; idInSpecInput_0 < inpSizForGivenNF;
              idInSpecInput_0++)
         {
            inpCID_0 = inpCIDsForGivenNF[idInSpecInput_0];
            
            for (vid = 0; vid < nVect; vid++)
            {
               inpCID = inpCID_0 + vid*nC;
               idInSpecInput = idInSpecInput_0 + vid*inpSizForGivenNF;
               pInputs[idInGivenNFs][idInSpecInput] = pC[inpCID];
            }
         }
      }
      
      // Evaluating and distributing to output
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         prhs = pmInputArrays[idInGivenNFs];
         mexCallMATLAB(1,plhs,2,prhs,"feval");
         pResult = mxGetPr(plhs[0]);
         idInNVrfx_0 = tnidsrfx[idInGivenNFs];
         
         for (vid = 0; vid < nVect; vid++)
         {
            idInNVrfx = idInNVrfx_0 + vid*nvnumrfx;
            pNVrfx[idInNVrfx] = pResult[vid];
         }
         
         mxDestroyArray(plhs[0]);
      }
      
      // Freeing input array pointers
      mxFree(pInputs);
   }
   else
   {
      // Preallocating
      pInputs = (double**)mxMalloc(nElems*sizeof(double*));
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         pInputs[idInGivenNFs] = mxGetPr(pmInputArrays[idInGivenNFs][1]);
      }
      
      // Filling
      #pragma omp parallel for private(idInGivenNFs,inpCIDsForGivenNF, \
              inpSizForGivenNF,idInSpecInput,inpCID) if(nElems >= ser2par)
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         inpCIDsForGivenNF = inpCIDs[idInGivenNFs];
         inpSizForGivenNF = inpSizes[idInGivenNFs];
         
         for (idInSpecInput = 0; idInSpecInput < inpSizForGivenNF;
              idInSpecInput++)
         {
            inpCID = inpCIDsForGivenNF[idInSpecInput];
            pInputs[idInGivenNFs][idInSpecInput] = pC[inpCID];
         }
      }
      
      // Evaluating and distributing to output
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         prhs = pmInputArrays[idInGivenNFs];
         mexCallMATLAB(1,plhs,2,prhs,"feval");
         idInNVrfx = tnidsrfx[idInGivenNFs];
         pNVrfx[idInNVrfx] = mxGetScalar(plhs[0]);
         
         mxDestroyArray(plhs[0]);
      }
      
      // Freeing input array pointers
      mxFree(pInputs);
   }
}

// Numeric function type 2 evaluation function
void evalNfh2(double* pNVrfx, const mwSize nElems, const mwIndex* tnidsrfx,
      const mwSize nVect, mxArray*** pmInputArrays)
{
   // Declarations
   mwIndex idInGivenNFs; // Index in given numeric functions
   mwIndex idInNVrfx; // Index in destination NV values
   mwIndex vid; // Vector index
   mxArray* plhs[1]; // NV evaluation
   double result;
   
   // Looping through numeric functions
   for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
   {
      mexCallMATLAB(1,plhs,2,pmInputArrays[idInGivenNFs],"feval");
      result = mxGetScalar(plhs[0]);
      
      for (vid = 0; vid < nVect; vid++)
      {
         idInNVrfx = vid*(mwIndex)nvnumrfx + tnidsrfx[idInGivenNFs];
         pNVrfx[idInNVrfx] = result;
      }
   }
   
   mxDestroyArray(plhs[0]);
}

// Numeric function type 4 evaluation function
void evalNfh4(double* pNVrfx, const mwSize nElems, const mwIndex* tnidsrfx,
      const double* pC, const mwSize nVect, mxArray*** pmInputArrays,
      const mwSize* inpSizes, const mwIndex** inpCIDs)
{
   // Declarations
   mwSignedIndex idInGivenNFs; // Index in given numeric functions
   mwIndex idInSpecInput; // Index in a specific input array
   mwIndex inpCID; // Index in C-array
   mwIndex idInNVrfx; // Index in destination NV values
   const mwIndex* inpCIDsForGivenNF; // Array of inpCIDs for a specific NV
   mwSize inpSizForGivenNF; // Number of inputs for a specific NV
   double** pInputs; // pInputs[idInGivenNFs] - Input array (double*)
   mxArray** prhs; // Function handle, Time, Modulating input
   mxArray* plhs[1]; // NV evaluation (mxArray)
   
   // Checking if input preallocation is required
   if (ifVectorReq)
   {
      // Declarations
      mwIndex idInSpecInput_0; // idInSpecInput, vector 0
      mwIndex inpCID_0; // inpCID, vector 0
      mwIndex idInNVrfx_0; // idInNVrfx, vector 0
      mwIndex vid; // Vector index
      double* pResult; // NV evaluation (double)
      
      // Preallocating
      pInputs = (double**)mxMalloc(nElems*sizeof(double*));
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         // Preallocating the modulating input array
         pmInputArrays[idInGivenNFs][2] = \
            mxPreallocDoubleMatrix(inpSizes[idInGivenNFs],nVect);
         pInputs[idInGivenNFs] = mxGetPr(pmInputArrays[idInGivenNFs][2]);
      }
      
      // Filling
      #pragma omp parallel for private(idInGivenNFs,inpCIDsForGivenNF, \
              inpSizForGivenNF,idInSpecInput_0,inpCID_0,vid,inpCID, \
              idInSpecInput) if(nElems >= ser2par)
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         inpCIDsForGivenNF = inpCIDs[idInGivenNFs];
         inpSizForGivenNF = inpSizes[idInGivenNFs];
         
         for (idInSpecInput_0 = 0; idInSpecInput_0 < inpSizForGivenNF;
              idInSpecInput_0++)
         {
            inpCID_0 = inpCIDsForGivenNF[idInSpecInput_0];
            
            for (vid = 0; vid < nVect; vid++)
            {
               inpCID = inpCID_0 + vid*nC;
               idInSpecInput = idInSpecInput_0 + vid*inpSizForGivenNF;
               pInputs[idInGivenNFs][idInSpecInput] = pC[inpCID];
            }
         }
      }
      
      // Evaluating and distributing to output
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         prhs = pmInputArrays[idInGivenNFs];
         mexCallMATLAB(1,plhs,3,prhs,"feval");
         pResult = mxGetPr(plhs[0]);
         idInNVrfx_0 = tnidsrfx[idInGivenNFs];
         
         for (vid = 0; vid < nVect; vid++)
         {
            idInNVrfx = idInNVrfx_0 + vid*nvnumrfx;
            pNVrfx[idInNVrfx] = pResult[vid];
         }
         
         mxDestroyArray(plhs[0]);
      }
      
      // Freeing input array pointers
      mxFree(pInputs);
   }
   else
   {
      // Preallocating
      pInputs = (double**)mxMalloc(nElems*sizeof(double*));
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         pInputs[idInGivenNFs] = mxGetPr(pmInputArrays[idInGivenNFs][2]);
      }
      
      // Filling
      #pragma omp parallel for private(idInGivenNFs,inpCIDsForGivenNF, \
              inpSizForGivenNF,idInSpecInput,inpCID) if(nElems >= ser2par)
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         inpCIDsForGivenNF = inpCIDs[idInGivenNFs];
         inpSizForGivenNF = inpSizes[idInGivenNFs];
         
         for (idInSpecInput = 0; idInSpecInput < inpSizForGivenNF;
              idInSpecInput++)
         {
            inpCID = inpCIDsForGivenNF[idInSpecInput];
            pInputs[idInGivenNFs][idInSpecInput] = pC[inpCID];
         }
      }
      
      // Evaluating and distributing to output
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         prhs = pmInputArrays[idInGivenNFs];
         mexCallMATLAB(1,plhs,3,prhs,"feval");
         idInNVrfx = tnidsrfx[idInGivenNFs];
         pNVrfx[idInNVrfx] = mxGetScalar(plhs[0]);
         
         mxDestroyArray(plhs[0]);
      }
      
      // Freeing input array pointers
      mxFree(pInputs);
   }
}

// Numeric function type 5 evaluation function
void evalNfh5(double* pNVrfx, const mwSize nElems, const mwIndex* tnidsrfx,
      const double* pC, const mwSize nVect, mxArray*** pmInputArrays,
      const mwSize* inpSizes, const mwIndex** inpCIDs, const mwSize* minpSizes,
      const mwIndex** minpCIDs)
{
   // Declarations
   mwSignedIndex idInGivenNFs; // Index in given numeric functions
   mwIndex idInSpecInput; // Index in a specific input array
   mwIndex inpCID; // Index in C-array
   mwIndex idInNVrfx; // Index in destination NV values
   const mwIndex* inpCIDsForGivenNF; // Array of inpCIDs for a specific NV
   mwSize inpSizForGivenNF; // Number of inputs for a specific NV
   double** pInputs; // pInputs[idInGivenNFs] - Reg. input array (double*)
   double** pMInputs; // pMInputs[idInGivenNFs] - Mod. input array (double*)
   mxArray** prhs; // Function handle, Regular input, Modulating input
   mxArray* plhs[1]; // NV evaluation (mxArray)
   
   // Checking if input preallocation is required
   if (ifVectorReq)
   {
      // Declarations
      mwIndex idInSpecInput_0; // idInSpecInput, vector 0
      mwIndex inpCID_0; // inpCID, vector 0
      mwIndex idInNVrfx_0; // idInNVrfx, vector 0
      mwIndex vid; // Vector index
      double* pResult; // NV evaluation (double)
      
      // Preallocating
      pInputs = (double**)mxMalloc(nElems*sizeof(double*));
      pMInputs = (double**)mxMalloc(nElems*sizeof(double*));
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         // Preallocating the regular input array
         pmInputArrays[idInGivenNFs][1] = \
            mxPreallocDoubleMatrix(inpSizes[idInGivenNFs],nVect);
         pInputs[idInGivenNFs] = mxGetPr(pmInputArrays[idInGivenNFs][1]);
         
         // Preallocating the modulating input array
         pmInputArrays[idInGivenNFs][2] = \
            mxPreallocDoubleMatrix(minpSizes[idInGivenNFs],nVect);
         pMInputs[idInGivenNFs] = mxGetPr(pmInputArrays[idInGivenNFs][2]);
      }
      
      // Filling
      #pragma omp parallel private(idInGivenNFs,inpCIDsForGivenNF, \
              inpSizForGivenNF,idInSpecInput_0,inpCID_0,vid,inpCID, \
              idInSpecInput) if(nElems >= ser2par)
      {
         #pragma omp for nowait
         for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
         {
            inpCIDsForGivenNF = inpCIDs[idInGivenNFs];
            inpSizForGivenNF = inpSizes[idInGivenNFs];
            
            for (idInSpecInput_0 = 0; idInSpecInput_0 < inpSizForGivenNF;
                 idInSpecInput_0++)
            {
               inpCID_0 = inpCIDsForGivenNF[idInSpecInput_0];
               
               for (vid = 0; vid < nVect; vid++)
               {
                  inpCID = inpCID_0 + vid*nC;
                  idInSpecInput = idInSpecInput_0 + vid*inpSizForGivenNF;
                  pInputs[idInGivenNFs][idInSpecInput] = pC[inpCID];
               }
            }
         }
         
         #pragma omp for
         for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
         {
            inpCIDsForGivenNF = minpCIDs[idInGivenNFs];
            inpSizForGivenNF = minpSizes[idInGivenNFs];
            
            for (idInSpecInput_0 = 0; idInSpecInput_0 < inpSizForGivenNF;
                 idInSpecInput_0++)
            {
               inpCID_0 = inpCIDsForGivenNF[idInSpecInput_0];
               
               for (vid = 0; vid < nVect; vid++)
               {
                  inpCID = inpCID_0 + vid*nC;
                  idInSpecInput = idInSpecInput_0 + vid*inpSizForGivenNF;
                  pMInputs[idInGivenNFs][idInSpecInput] = pC[inpCID];
               }
            }
         }
      }
      
      // Evaluating and distributing to output
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         prhs = pmInputArrays[idInGivenNFs];
         mexCallMATLAB(1,plhs,3,prhs,"feval");
         pResult = mxGetPr(plhs[0]);
         idInNVrfx_0 = tnidsrfx[idInGivenNFs];
         
         for (vid = 0; vid < nVect; vid++)
         {
            idInNVrfx = idInNVrfx_0 + vid*nvnumrfx;
            pNVrfx[idInNVrfx] = pResult[vid];
         }
         
         mxDestroyArray(plhs[0]);
      }
      
      // Freeing input array pointers
      mxFree(pInputs);
      mxFree(pMInputs);
   }
   else
   {
      // Preallocating
      pInputs = (double**)mxMalloc(nElems*sizeof(double*));
      pMInputs = (double**)mxMalloc(nElems*sizeof(double*));
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         pInputs[idInGivenNFs] = mxGetPr(pmInputArrays[idInGivenNFs][1]);
         pMInputs[idInGivenNFs] = mxGetPr(pmInputArrays[idInGivenNFs][2]);
      }
      
      // Filling
      #pragma omp parallel private(idInGivenNFs,inpCIDsForGivenNF, \
              inpSizForGivenNF,inpCID,idInSpecInput) if(nElems >= ser2par)
      {
         #pragma omp for nowait
         for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
         {
            inpCIDsForGivenNF = inpCIDs[idInGivenNFs];
            inpSizForGivenNF = inpSizes[idInGivenNFs];
            
            for (idInSpecInput = 0; idInSpecInput < inpSizForGivenNF;
                 idInSpecInput++)
            {
               inpCID = inpCIDsForGivenNF[idInSpecInput];
               pInputs[idInGivenNFs][idInSpecInput] = pC[inpCID];
            }
         }
         
         #pragma omp for
         for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
         {
            inpCIDsForGivenNF = minpCIDs[idInGivenNFs];
            inpSizForGivenNF = minpSizes[idInGivenNFs];
            
            for (idInSpecInput = 0; idInSpecInput < inpSizForGivenNF;
                 idInSpecInput++)
            {
               inpCID = inpCIDsForGivenNF[idInSpecInput];
               pMInputs[idInGivenNFs][idInSpecInput] = pC[inpCID];
            }
         }
      }
      
      // Evaluating and distributing to output
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         prhs = pmInputArrays[idInGivenNFs];
         mexCallMATLAB(1,plhs,3,prhs,"feval");
         idInNVrfx = tnidsrfx[idInGivenNFs];
         pNVrfx[idInNVrfx] = mxGetScalar(plhs[0]);
         
         mxDestroyArray(plhs[0]);
      }
      
      // Freeing input array pointers
      mxFree(pInputs);
      mxFree(pMInputs);
   }
}

// Numeric function type 6 evaluation function
void evalNfh6(double* pNVrfx, const mwSize nElems, const mwIndex* tnidsrfx,
      const double* pC, const mwSize nVect, mxArray*** pmInputArrays,
      const mwSize* inpSizes, const mwIndex** inpCIDs)
{
   // Declarations
   mwSignedIndex idInGivenNFs; // Index in given numeric functions
   mwIndex idInSpecInput; // Index in a specific input array
   mwIndex inpCID; // Index in C-array
   mwIndex idInNVrfx; // Index in destination NV values
   const mwIndex* inpCIDsForGivenNF; // Array of inpCIDs for a specific NV
   mwSize inpSizForGivenNF; // Number of inputs for a specific NV
   double** pInputs; // pInputs[idInGivenNFs] - Input array (double*)
   mxArray** prhs; // Function handle, Regular input, Time
   mxArray* plhs[1]; // NV evaluation (mxArray)
   
   // Checking if input preallocation is required
   if (ifVectorReq)
   {
      // Declarations
      mwIndex idInSpecInput_0; // idInSpecInput, vector 0
      mwIndex inpCID_0; // inpCID, vector 0
      mwIndex idInNVrfx_0; // idInNVrfx, vector 0
      mwIndex vid; // Vector index
      double* pResult; // NV evaluation (double)
      
      // Preallocating
      pInputs = (double**)mxMalloc(nElems*sizeof(double*));
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         // Preallocating the regular input array
         pmInputArrays[idInGivenNFs][1] = \
            mxPreallocDoubleMatrix(inpSizes[idInGivenNFs],nVect);
         pInputs[idInGivenNFs] = mxGetPr(pmInputArrays[idInGivenNFs][1]);
      }
      
      // Filling
      #pragma omp parallel for private(idInGivenNFs,inpCIDsForGivenNF, \
              inpSizForGivenNF,idInSpecInput_0,inpCID_0,vid,inpCID, \
              idInSpecInput) if(nElems >= ser2par)
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         inpCIDsForGivenNF = inpCIDs[idInGivenNFs];
         inpSizForGivenNF = inpSizes[idInGivenNFs];
         
         for (idInSpecInput_0 = 0; idInSpecInput_0 < inpSizForGivenNF;
              idInSpecInput_0++)
         {
            inpCID_0 = inpCIDsForGivenNF[idInSpecInput_0];
            
            for (vid = 0; vid < nVect; vid++)
            {
               inpCID = inpCID_0 + vid*nC;
               idInSpecInput = idInSpecInput_0 + vid*inpSizForGivenNF;
               pInputs[idInGivenNFs][idInSpecInput] = pC[inpCID];
            }
         }
      }
      
      // Evaluating and distributing to output
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         prhs = pmInputArrays[idInGivenNFs];
         mexCallMATLAB(1,plhs,3,prhs,"feval");
         pResult = mxGetPr(plhs[0]);
         idInNVrfx_0 = tnidsrfx[idInGivenNFs];
         
         for (vid = 0; vid < nVect; vid++)
         {
            idInNVrfx = idInNVrfx_0 + vid*nvnumrfx;
            pNVrfx[idInNVrfx] = pResult[vid];
         }
         
         mxDestroyArray(plhs[0]);
      }
      
      // Freeing input array pointers
      mxFree(pInputs);
   }
   else
   {
      // Preallocating
      pInputs = (double**)mxMalloc(nElems*sizeof(double*));
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         pInputs[idInGivenNFs] = mxGetPr(pmInputArrays[idInGivenNFs][1]);
      }
      
      // Filling
      #pragma omp parallel for private(idInGivenNFs,inpCIDsForGivenNF, \
              inpSizForGivenNF,idInSpecInput,inpCID) if(nElems >= ser2par)
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         inpCIDsForGivenNF = inpCIDs[idInGivenNFs];
         inpSizForGivenNF = inpSizes[idInGivenNFs];
         
         for (idInSpecInput = 0; idInSpecInput < inpSizForGivenNF;
              idInSpecInput++)
         {
            inpCID = inpCIDsForGivenNF[idInSpecInput];
            pInputs[idInGivenNFs][idInSpecInput] = pC[inpCID];
         }
      }
      
      // Evaluating and distributing to output
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         prhs = pmInputArrays[idInGivenNFs];
         mexCallMATLAB(1,plhs,3,prhs,"feval");
         idInNVrfx = tnidsrfx[idInGivenNFs];
         pNVrfx[idInNVrfx] = mxGetScalar(plhs[0]);
         
         mxDestroyArray(plhs[0]);
      }
      
      // Freeing input array pointers
      mxFree(pInputs);
   }
}

// Numeric function type 7 evaluation function
void evalNfh7(double* pNVrfx, const mwSize nElems, const mwIndex* tnidsrfx,
      const double* pC, const mwSize nVect, mxArray*** pmInputArrays,
      const mwSize* inpSizes, const mwIndex** inpCIDs, const mwSize* minpSizes,
      const mwIndex** minpCIDs)
{
   // Declarations
   mwSignedIndex idInGivenNFs; // Index in given numeric functions
   mwIndex idInSpecInput; // Index in a specific input array
   mwIndex inpCID; // Index in C-array
   mwIndex idInNVrfx; // Index in destination NV values
   const mwIndex* inpCIDsForGivenNF; // Array of inpCIDs for a specific NV
   mwSize inpSizForGivenNF; // Number of inputs for a specific NV
   double** pInputs; // pInputs[idInGivenNFs] - Reg. input array (double*)
   double** pMInputs; // pMInputs[idInGivenNFs] - Mod. input array (double*)
   mxArray** prhs; // Function handle, Regular input, Time, Modulating input
   mxArray* plhs[1]; // NV evaluation (mxArray)
   
   // Checking if input preallocation is required
   if (ifVectorReq)
   {
      // Declarations
      mwIndex idInSpecInput_0; // idInSpecInput, vector 0
      mwIndex inpCID_0; // inpCID, vector 0
      mwIndex idInNVrfx_0; // idInNVrfx, vector 0
      mwIndex vid; // Vector index
      double* pResult; // NV evaluation (double)
      
      // Preallocating
      pInputs = (double**)mxMalloc(nElems*sizeof(double*));
      pMInputs = (double**)mxMalloc(nElems*sizeof(double*));
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         // Preallocating the regular input array
         pmInputArrays[idInGivenNFs][1] = \
            mxPreallocDoubleMatrix(inpSizes[idInGivenNFs],nVect);
         pInputs[idInGivenNFs] = mxGetPr(pmInputArrays[idInGivenNFs][1]);
         
         // Preallocating the modulating input array
         pmInputArrays[idInGivenNFs][3] = \
            mxPreallocDoubleMatrix(minpSizes[idInGivenNFs],nVect);
         pMInputs[idInGivenNFs] = mxGetPr(pmInputArrays[idInGivenNFs][3]);
      }
      
      // Filling
      #pragma omp parallel private(idInGivenNFs,inpCIDsForGivenNF, \
              inpSizForGivenNF,idInSpecInput_0,inpCID_0,vid,inpCID, \
              idInSpecInput) if(nElems >= ser2par)
      {
         #pragma omp for nowait
         for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
         {
            inpCIDsForGivenNF = inpCIDs[idInGivenNFs];
            inpSizForGivenNF = inpSizes[idInGivenNFs];
            
            for (idInSpecInput_0 = 0; idInSpecInput_0 < inpSizForGivenNF;
                 idInSpecInput_0++)
            {
               inpCID_0 = inpCIDsForGivenNF[idInSpecInput_0];
               
               for (vid = 0; vid < nVect; vid++)
               {
                  inpCID = inpCID_0 + vid*nC;
                  idInSpecInput = idInSpecInput_0 + vid*inpSizForGivenNF;
                  pInputs[idInGivenNFs][idInSpecInput] = pC[inpCID];
               }
            }
         }
         
         #pragma omp for
         for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
         {
            inpCIDsForGivenNF = minpCIDs[idInGivenNFs];
            inpSizForGivenNF = minpSizes[idInGivenNFs];
            
            for (idInSpecInput_0 = 0; idInSpecInput_0 < inpSizForGivenNF;
                 idInSpecInput_0++)
            {
               inpCID_0 = inpCIDsForGivenNF[idInSpecInput_0];
               
               for (vid = 0; vid < nVect; vid++)
               {
                  inpCID = inpCID_0 + vid*nC;
                  idInSpecInput = idInSpecInput_0 + vid*inpSizForGivenNF;
                  pMInputs[idInGivenNFs][idInSpecInput] = pC[inpCID];
               }
            }
         }
      }
      
      // Evaluating and distributing to output
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         prhs = pmInputArrays[idInGivenNFs];
         mexCallMATLAB(1,plhs,4,prhs,"feval");
         pResult = mxGetPr(plhs[0]);
         idInNVrfx_0 = tnidsrfx[idInGivenNFs];
         
         for (vid = 0; vid < nVect; vid++)
         {
            idInNVrfx = idInNVrfx_0 + vid*nvnumrfx;
            pNVrfx[idInNVrfx] = pResult[vid];
         }
         
         mxDestroyArray(plhs[0]);
      }
      
      // Freeing input array pointers
      mxFree(pInputs);
      mxFree(pMInputs);
   }
   else
   {
      // Preallocating
      pInputs = (double**)mxMalloc(nElems*sizeof(double*));
      pMInputs = (double**)mxMalloc(nElems*sizeof(double*));
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         pInputs[idInGivenNFs] = mxGetPr(pmInputArrays[idInGivenNFs][1]);
         pMInputs[idInGivenNFs] = mxGetPr(pmInputArrays[idInGivenNFs][3]);
      }
      
      // Filling
      #pragma omp parallel private(idInGivenNFs,inpCIDsForGivenNF, \
              inpSizForGivenNF,inpCID,idInSpecInput) if(nElems >= ser2par)
      {
         #pragma omp for nowait
         for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
         {
            inpCIDsForGivenNF = inpCIDs[idInGivenNFs];
            inpSizForGivenNF = inpSizes[idInGivenNFs];
            
            for (idInSpecInput = 0; idInSpecInput < inpSizForGivenNF;
                 idInSpecInput++)
            {
               inpCID = inpCIDsForGivenNF[idInSpecInput];
               pInputs[idInGivenNFs][idInSpecInput] = pC[inpCID];
            }
         }
         
         #pragma omp for
         for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
         {
            inpCIDsForGivenNF = minpCIDs[idInGivenNFs];
            inpSizForGivenNF = minpSizes[idInGivenNFs];
            
            for (idInSpecInput = 0; idInSpecInput < inpSizForGivenNF;
                 idInSpecInput++)
            {
               inpCID = inpCIDsForGivenNF[idInSpecInput];
               pMInputs[idInGivenNFs][idInSpecInput] = pC[inpCID];
            }
         }
      }
      
      // Evaluating and distributing to output
      for (idInGivenNFs = 0; idInGivenNFs < nElems; idInGivenNFs++)
      {
         prhs = pmInputArrays[idInGivenNFs];
         mexCallMATLAB(1,plhs,4,prhs,"feval");
         idInNVrfx = tnidsrfx[idInGivenNFs];
         pNVrfx[idInNVrfx] = mxGetScalar(plhs[0]);
         
         mxDestroyArray(plhs[0]);
      }
      
      // Freeing input array pointers
      mxFree(pInputs);
      mxFree(pMInputs);
   }
}

// Cleanup function
static void evalxdotCleanup(void)
{
   // Declarations
   // Indices
   mwIndex lid; // Layer
   mwIndex bviid; // BV input
   mwIndex nvtyp; // NV type
   mwIndex ltnid; // Layer- and type-local NV
   
   // Freeing static time array
   mxDestroyArray(st_pmt);
   
   // Freeing BV inputs
   if (!ifVectorReq)
   {
      // Freeing all but time BV inputs
      for (bviid = 0; bviid < nBVInps; bviid++)
      {
         if (!pIfBVInpIDt[bviid])
         {
            mxDestroyArray(pmBVInps[bviid]);
         }
      }
   }
   mxFree(pmBVInps);
   mxFree(pBVInps);
   mxFree(pIfBVInpIDt);
   
   // Freeing evaluated BV and BV input indices
   for (lid = 0; lid < lnumrfx; lid++)
   {
      mxFree(pBVlbvids[lid]);
      mxFree(pBVInpXIDs[lid]);
      mxFree(pBVInpNIDs[lid]);
   }
   mxFree(pBVlbvids);
   mxFree(pBVInpXIDs);
   mxFree(pBVInpNIDs);
   
   // Freeing BV input sizes
   mxFree(pBVInpSiz);
   mxFree(pBVInpXSiz);
   mxFree(pBVInpNSiz);
   
   // Freeing BV function handles
   mxFree(pmBVfhs);
   
   // Freeing BV I/O arrays
   for (lid = 0; lid < lnumrfx; lid++)
   {
      if (nBVs[lid] > 0)
      {
         mxFree(pmBVrhs[lid]);
      }
   }
   mxFree(pmBVrhs);
   
   // Freeing BV sizes
   mxFree(nBVs);
   
   // Freeing NV inputs and function handles
   for (lid = 0; lid < (lnumrfx-1); lid++)
   {
      for (nvtyp = 0; nvtyp < 7; nvtyp++)
      {
         for (ltnid = 0; ltnid < nNVs[nvtyp+7*lid]; ltnid++)
         {
            // Freeing NV inputs
            // (if not vectorized, otherwise they are not persistent)
            if (!ifVectorReq)
            {
               switch (nvtyp)
               {
                  case 0:
                  case 2:
                  case 5:
                     mxDestroyArray(pmNVInps[nvtyp+7*lid][ltnid][1]);
                     break;
                  case 3:
                     mxDestroyArray(pmNVInps[nvtyp+7*lid][ltnid][2]);
                     break;
                  case 4:
                     mxDestroyArray(pmNVInps[nvtyp+7*lid][ltnid][1]);
                     mxDestroyArray(pmNVInps[nvtyp+7*lid][ltnid][2]);
                     break;
                  case 6:
                     mxDestroyArray(pmNVInps[nvtyp+7*lid][ltnid][1]);
                     mxDestroyArray(pmNVInps[nvtyp+7*lid][ltnid][3]);
                     break;
               }
            }
            mxFree(pmNVInps[nvtyp+7*lid][ltnid]);
         }
         mxFree(pmNVInps[nvtyp+7*lid]);
      }
   }
   mxFree(pmNVInps);
   
   // Freeing NV input indices
   for (lid = 0; lid < (lnumrfx-1); lid++)
   {
      for (nvtyp = 0; nvtyp < 7; nvtyp++)
      {
         // Freeing regular NV input indices
         if (nvtyp == 0 || nvtyp == 4 || nvtyp == 5 || nvtyp == 6)
         {
            for (ltnid = 0; ltnid < nNVs[nvtyp+7*lid]; ltnid++)
            {
               mxFree(pNVInpIDs[nvtyp+7*lid][ltnid]);
            }
            mxFree(pNVInpIDs[nvtyp+7*lid]);
         }
         // Freeing modulating NV input indices
         if (nvtyp == 2 || nvtyp == 3 || nvtyp == 4 || nvtyp == 6)
         {
            for (ltnid = 0; ltnid < nNVs[nvtyp+7*lid]; ltnid++)
            {
               mxFree(pNVMInpIDs[nvtyp+7*lid][ltnid]);
            }
            mxFree(pNVMInpIDs[nvtyp+7*lid]);
         }
      }
   }
   mxFree(pNVInpIDs);
   mxFree(pNVMInpIDs);
   
   // Freeing evaluated NV indices
   for (lid = 0; lid < (lnumrfx-1); lid++)
   {
      for (nvtyp = 0; nvtyp < 7; nvtyp++)
      {
         mxFree(pNVltnids[nvtyp+7*lid]);
      }
   }
   mxFree(pNVltnids);
   
   // Freeing NV input sizes
   for (lid = 0; lid < (lnumrfx-1); lid++)
   {
      for (nvtyp = 0; nvtyp < 7; nvtyp++)
      {
         // Freeing regular NV input sizes
         if (nvtyp == 0 || nvtyp == 4 || nvtyp == 5 || nvtyp == 6)
         {
            mxFree(pNVInpSiz[nvtyp+7*lid]);
         }
         // Freeing modulating NV input sizes
         if (nvtyp == 2 || nvtyp == 3 || nvtyp == 4 || nvtyp == 6)
         {
            mxFree(pNVMInpSiz[nvtyp+7*lid]);
         }
      }
   }
   mxFree(pNVInpSiz);
   mxFree(pNVMInpSiz);
   
   // Freeing NV sizes
   mxFree(nNVs);
   
   // Freeing the final xdot indices
   mxFree(pXDotIDs);
   
   return;
}