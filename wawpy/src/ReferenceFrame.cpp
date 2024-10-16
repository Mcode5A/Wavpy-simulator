#ifndef     REFERENCEFRAME_CPP_
#define     REFERENCEFRAME_CPP_

/**
    \file   ReferenceFrame.cpp
    \brief  Implementation of the ReferenceFrame class.
*/


#include "ReferenceFrame.hpp"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <string>
#include <iostream>
#include <iomanip>

//--------------------------- CheckPositionDimensions
bool ReferenceFrame::CheckPositionDimensions(gsl_matrix *Pos, unsigned int L)
{
  bool retval=true;
  // Check that it has x,y,z coordinates
  if (REFFRAME_DIMS != Pos->size1)
  {
    retval = false;
  }
  if (L != Pos->size2)
  {
    retval = false;
  }
  return retval;
}

//--------------------------- CheckAttitudeDimensions
bool ReferenceFrame::CheckAttitudeDimensions(gsl_matrix *Att, unsigned int L)
{
  bool retval=true;
  // Check that it has x,y,z coordinates
  if (REFFRAME_DIMS != Att->size1)
  {
    retval = false;
  }
  if (REFFRAME_DIMS*L != Att->size2)
  {
    retval = false;
  }
  return retval;
}

//--------------------------- CheckVelocityDimensions
bool ReferenceFrame::CheckVelocityDimensions(gsl_matrix *Vel, unsigned int L)
{
  bool retval=true;
  // Check that it has x,y,z coordinates
  if (REFFRAME_DIMS != Vel->size1)
  {
    retval = false;
  }
  if (L != Vel->size2)
  {
    retval = false;
  }
  return retval;
}

//--------------------------- CheckVectorLength
bool ReferenceFrame::CheckVectorLength(gsl_vector* v, unsigned int len)
{
  if (NULL==v)
  {
    return false;
  }
  if (v->size != len)
  {
    return false;
  }
  return true;
}

//--------------------------- CheckTimeInstant
bool ReferenceFrame::CheckTimeInstant(const ReferenceFrame *rf, unsigned int L)
{
  if (NULL==rf)
  {
    return false;
  }
  if ( ! (rf->mSolidary && rf->mAllReferredSolidary) ) // there is some non-solidary ReferenceFrame
  {
    if (L>=rf->mL) // L is larger than the number of elements in the ReferenceFrame
    { 
      return false;
    }
  }
  return true;
}

//--------------------------- CheckEndsAtECEF
bool ReferenceFrame::CheckEndsAtECEF(const ReferenceFrame *rf, std::list<const ReferenceFrame*> *mylist)
{
  //std::cout << "Checking " << rf << std::endl;
  // if list does not exists, create list
  if (NULL==mylist) 
  {
    //std::cout << "  Creating list" << std::endl;
    mylist = new std::list<const ReferenceFrame*>();
  }
  if (!mylist->empty()) // there are already some reference frames in the list
  {
    //std::cout << "  Checking if " << rf << " already in list." << std::endl;
    // Check if current rf is already in list
    for (std::list<const ReferenceFrame*>::iterator it=mylist->begin(); it !=mylist->end(); it++)
    {
      //std::cout << "  iterator: " << *it << "  rf pointer: " << rf << std::endl;
      if (*it == rf)
      {
        //std::cout << "  Reference frame: " << rf << " is already in list." << std::endl;
        //std::cout << "  Removing " << rf << " from list." << std::endl;
        mylist->pop_back();
        return false;
      }
    }
  }

  // we reach this point if rf is not in list
  {
    //std::cout << "  " << rf << " not in list." << std::endl;
    // add current rf to list
    mylist->push_back(rf);
    //std::cout << "  Adding ReferenceFrame to list: " << rf << std::endl;
    
    // check if it is ECEF
    if (NULL==rf->mMyRefFrame) // ECEF Reference frame found
    {
      //std::cout << "  ECEF found. " << std::endl;
      return true;
    }
    else // iterate recursively (with mMyRefFrame) 
    {
      const ReferenceFrame *myrf = rf->GetMyReferenceFrame();
      bool b = CheckEndsAtECEF(myrf, mylist);
      //std::cout << "  Removing one frame." << std::endl;
      mylist->pop_back();
      if (mylist->empty())
      {
        //std::cout << "  Deleting list" << std::endl; // potser posar això al final de tot.
        delete(mylist);
      }
      return b;
    }
  }
}


//--------------------------- ReferenceFrame
ReferenceFrame::ReferenceFrame():ReferenceFrame(NULL)
{
  return;
}


//--------------------------- ReferenceFrame
ReferenceFrame::ReferenceFrame(ReferenceFrame *RefFrame, bool Solidary, gsl_matrix *Pos, gsl_matrix *Att, gsl_matrix *Vel, std::string Name)
{
  if (NULL==RefFrame)  // The global reference frame
  {
    mMyRefFrame = NULL;
    mSolidary = true;
    mAllReferredSolidary = true;
    //mAllSolidary=mSolidary;
    mL = 1;
    mName = std::string("ECEF");
    //mLmax = mL;
    
    // gsl library already includes error handling for malloc & calloc
    mPos = gsl_matrix_calloc(REFFRAME_DIMS,mL);
    mAtt = gsl_matrix_alloc(REFFRAME_DIMS,mL*REFFRAME_DIMS);
    gsl_matrix_set_identity(mAtt);
    mVel = gsl_matrix_calloc(REFFRAME_DIMS,mL);
    return;
  }
  else
  {
    mName = Name;
    mMyRefFrame = RefFrame;
  }

  if (true==Solidary)
  {
    mSolidary = true;
    mAllReferredSolidary = mMyRefFrame->mAllReferredSolidary && mMyRefFrame->mSolidary; // inherit if all other are also solidary.
    //mL=1; // solidary frames are constant with time, respect to their reference frame, so mL=1.
    //mLmax = mMyRefFrame->mLmax; // mLmax inherited from MyRefFrame
    mL = mMyRefFrame->mL;    // We take the length of our referred reference 
                             // frame, so we will be able to know what is the
                             // maximum lenght of all the tree.
    //mAllSolidary = mMyRefFrame->mAllSolidary && mMyRefFrame->mSolidary; // inherit if all other are also solidary.
    //std::cout << "solidary: " << mSolidary << " mAllSolidary:" << mAllSolidary << std::endl;
    if (NULL==Pos) // No position defined => we use the origin
    {
      mPos = gsl_matrix_calloc(REFFRAME_DIMS,mL);
    }
    else // Assign position
    {
      if (CheckPositionDimensions(Pos,mL))  // Check that the position is correct
      {
        mPos = gsl_matrix_calloc(REFFRAME_DIMS, mL);
        gsl_matrix_memcpy(mPos,Pos);
      }
      else
      {
        throw("Position matrix has wrong size.");
      }
    }
    if (NULL==Att) // Attitude not defined => frame aligned with its reference frame
    {
      mAtt = gsl_matrix_alloc(REFFRAME_DIMS,mL*REFFRAME_DIMS);
      gsl_matrix_set_identity(mAtt);
    }
    else // Assign attitude
    {
      if (CheckAttitudeDimensions(Att, mL))
      {
        mAtt = gsl_matrix_alloc(REFFRAME_DIMS, mL*REFFRAME_DIMS);
        gsl_matrix_memcpy(mAtt, Att);
      }
      else
      {
        throw("Attitude matrix has wrong size.");
      }
    }
    if (NULL==mVel) // Velocity not defined => we use zero velocity
    {
      mVel = gsl_matrix_calloc(REFFRAME_DIMS, mL);
    }
    else // Assign velocity
    {
      //if (CheckVelocityDimensions(Vel,mL))
      //{
        mVel = gsl_matrix_calloc(REFFRAME_DIMS, mL); // velocity is zero for solidary frames
                                                     // We ignore the provided velocity.
        //gsl_matrix_memcpy(mVel,Vel);
      //}
      //else
      //{
      //  throw("Velocity matrix has wrong size.");
      //}
    }
  }
  else // non solidary frame
  {
    // Check that none of the vectors (position, attitiude, velocity) is null
    if (NULL==Pos || NULL==Att || NULL==Vel)
    {
      throw("Non-solidary ReferenceFrames must explicitly state their position, attitude, and velocity arrays.");
    }
    // Check that the position, attitude and velocity dimensions are consistent
    if ( !(CheckPositionDimensions(Pos,Pos->size2) && 
           CheckAttitudeDimensions(Att,Pos->size2) && 
           CheckVelocityDimensions(Vel,Pos->size2))   )
    {
      throw("Size of Position, Attitude and/or Velocity matrices are wrong.");
    }

    mSolidary = false;
    mAllReferredSolidary = mMyRefFrame->mAllReferredSolidary && mMyRefFrame->mSolidary; // inherit from referred frame
    if (true == mAllReferredSolidary) // all refered solidary
    {
      mL = Pos->size2; // we take the length of the Position vector as the correct one.
    }
    else
    {
      // Check that the provided length is consistent witht the length of the referred frames
      if (Pos->size2 != mMyRefFrame->mL)
      {
        throw("Length of the Reference Frame is different from length of the referred Reference Frame.");
      }
      mL = Pos->size2;
    }
    // Assign position, atttitude and velocity
    mPos = gsl_matrix_alloc(REFFRAME_DIMS,mL);
    mAtt = gsl_matrix_alloc(REFFRAME_DIMS,mL*REFFRAME_DIMS);
    mVel = gsl_matrix_alloc(REFFRAME_DIMS,mL);

    gsl_matrix_memcpy(mPos, Pos);
    gsl_matrix_memcpy(mAtt, Att);
    gsl_matrix_memcpy(mVel, Vel);
  }
  return;
}


//--------------------------- ~ReferenceFrame
ReferenceFrame::~ReferenceFrame()
{
  mMyRefFrame = NULL;
  mL = 0;
  gsl_matrix_free(mPos);
  gsl_matrix_free(mAtt);
  gsl_matrix_free(mVel);
}

//--------------------------- IsSolidary
bool ReferenceFrame::IsSolidary() const
{
  return mSolidary;
}

//--------------------------- IsMainReferenceFrame
bool ReferenceFrame::IsMainReferenceFrame() const
{
  bool retval=true;
  // Check if it refers to nobody
  if (NULL!=mMyRefFrame)
  {
    retval = false;
    return retval;
  }
  // Check if is solidary
  if (!IsSolidary()) // is solidary
  {
    retval=false;
    return retval;
  }
  // Check if is centred at the origin
  if (!gsl_matrix_isnull(mPos))
  {
    retval=false;
    return retval;
  }
  // Check if it is oriented as ECEF
  for (int k=0;k<REFFRAME_DIMS;k++)
  {
    for (int l=0;l<REFFRAME_DIMS;l++)
    {
      double v = gsl_matrix_get(mAtt,k,l);
      if (k==l) // diagonal terms shoud be =1
      {
        if (v!=1) { retval = false; return retval; }
      }
      else // off-diagonal terms should be !=0
      {
        if (v!=0) { retval = false; return retval; }
      }
    }
  }
  // Check the speed
  if (!gsl_matrix_isnull(mVel))
  {
    retval=false;
    return retval;
  }
  return retval;
}

//--------------------------- GetMyReferenceFrame
const ReferenceFrame* ReferenceFrame::GetMyReferenceFrame() const
{
  return mMyRefFrame;
}

//--------------------------- GetLength
unsigned int ReferenceFrame::GetLength() const
{
  return mL;
}

//--------------------------- GetPosition
gsl_vector* ReferenceFrame::GetPosition(unsigned int L) const
{
  unsigned int col;
  if (true==mSolidary) // ignores the instant L, as it is always the same.
  {
    col=0;
  }
  else
  {
    if (L>=mL) //
    {
      throw("Requesting position at non-existing instant. (Instant larger than elements in array).");
    }
    else
    {
      col = L;
    }
  }
  gsl_vector* pos;
  pos = gsl_vector_alloc(REFFRAME_DIMS);
  gsl_matrix_get_col(pos,mPos,col);
  return pos;
}

//--------------------------- GetAttitude
gsl_matrix* ReferenceFrame::GetAttitude(unsigned int L) const
{
  unsigned int col;
  if (true==mSolidary) // ignores instant L, as it is always the same
  {
    col=0;
  }
  else
  {
    if (L>=mL)
    {
      throw("Requesting attitude at non-existing instant. (Instant larger than elements in array).");
    }
    else
    {
      col = L*REFFRAME_DIMS;
    }
  }
  gsl_matrix_view v = gsl_matrix_submatrix(mAtt,0,col,REFFRAME_DIMS,REFFRAME_DIMS);
  gsl_matrix* m;
  m = gsl_matrix_alloc(REFFRAME_DIMS,REFFRAME_DIMS);
  gsl_matrix_memcpy(m,&v.matrix);
  return m;
}

//--------------------------- GetVelocity
gsl_vector* ReferenceFrame::GetVelocity(unsigned int L) const
{
  unsigned int col;
  if (true==mSolidary) // ignores the instant L, as it is always the same.
  {
    col=0;
  }
  else
  {
    if (L>=mL) //
    {
      throw("Requesting velocity at non-existing instant. (Instant larger than elements in array).");
    }
    else
    {
      col = L;
    }
  }
  gsl_vector* vel;
  vel = gsl_vector_alloc(REFFRAME_DIMS);
  gsl_matrix_get_col(vel,mVel,col);
  return vel;
}


//--------------------------- Print
void ReferenceFrame::Print() const
{
  std::cout << "-------- ReferenceFrame" << std::endl;
  std::cout << "  this:        " << this        << std::endl;
  std::cout << "  mName        " << mName       << std::endl;
  std::cout << "  mMyRefFrame: " << mMyRefFrame << std::endl;
  std::cout << "  mSolidary:   " << mSolidary   << std::endl;
  std::cout << "  mL:          " << mL          << std::endl;
  std::cout << "  Position:    "                << std::endl;
  for (unsigned int k=0;k<mL;k++)
  {
    std::cout <<   "  mPos[:," << k << "]=" 
              <<   "  [" << gsl_matrix_get(mPos,0,k) << "; "
              <<            gsl_matrix_get(mPos,1,k) << "; "
              <<            gsl_matrix_get(mPos,2,k) << "]'  " << std::endl;
  }
  std::cout << "  Attitude:    "                << std::endl;
  for (unsigned int k=0;k<mL;k++)
  {
    const int  w = 8;
    const char sep = ' ';
    std::cout <<   "  mAtt[:,:]_" << k << "=" << std::endl;
    for (unsigned int l=0;l<REFFRAME_DIMS;l++)
    {
      std::cout << "    " << std::left << std::setw(w) << std::setfill(sep) << gsl_matrix_get(mAtt,l,k*REFFRAME_DIMS)
                          << std::left << std::setw(w) << std::setfill(sep) << gsl_matrix_get(mAtt,l,k*REFFRAME_DIMS+1)
                          << std::left << std::setw(w) << std::setfill(sep) << gsl_matrix_get(mAtt,l,k*REFFRAME_DIMS+2) << std::endl;
    }
  }
  std::cout << "  Velocity:    " << std::endl;
  for (unsigned int k=0;k<mL;k++)
  {
    std::cout <<   "  mVel[:," << k << "]=" 
              <<   "  [" << gsl_matrix_get(mVel,0,k) << "; "
              <<            gsl_matrix_get(mVel,1,k) << "; "
              <<            gsl_matrix_get(mVel,2,k) << "]'  " << std::endl;
  }
  std::cout << std::endl;
}


//--------------------------- PosInParentFrame2PosInBF
gsl_vector* ReferenceFrame::PosInParentFrame2PosInBF(gsl_vector *p0, unsigned int L) const
{
  // Check the correctness of the input parameters
  // Check dimension of p0
  if (!CheckVectorLength(p0, REFFRAME_DIMS))
  {
    throw("Position vector with wrong size or undefined in PosInParentFrame2PosInBF.");
  }
  // Check the correctness of L
  if (mSolidary)  // Solidary frame
  {
    L=0; // to ignore the instant position
  }
  else
  {
    if (!CheckTimeInstant(this, L))
    {
      throw("Requesting to compute position at non-existent instant of time. L too large."); 
    }
  }
  gsl_vector *p1 = gsl_vector_alloc(REFFRAME_DIMS);

  // Get the attitude for the current time instant
  //gsl_matrix_view Att = gsl_matrix_submatrix(mAtt,0,L*REFFRAME_DIMS, REFFRAME_DIMS, REFFRAME_DIMS);
  gsl_matrix *Att = GetAttitude(L);
  //std::cout << "Showing attitude matrix: " << std::endl;
  //for (unsigned int row=0;row<REFFRAME_DIMS;row++)
  //{
  //  for (unsigned int col=0;col<REFFRAME_DIMS;col++)
  //  {
  //    std::cout << Att->data[row*REFFRAME_DIMS+col] << "  ";
  //  }
  //  std::cout << std::endl;
  //}
  // Prepare elements for matrix computation
  gsl_vector *p0_T;
  p0_T = gsl_vector_alloc(REFFRAME_DIMS);
  gsl_vector_memcpy(p0_T, p0);

  gsl_vector *T;
  T = GetPosition(L);           // Position of current instant
  gsl_vector_sub(p0_T,T);

  //std::cout << "- - - - " << std::endl;
  //Print();
  //std::cout << "p0= (" << p0->data[0] << ", " <<
  //                        p0->data[1] << ", " << 
  //                        p0->data[2] << ")" << std::endl;
  //std::cout << "T= (" << T->data[0] << ", " << 
  //                       T->data[1] << ", " <<
  //                       T->data[2] << ")" << std::endl;
  //std::cout << "p0_T= (" << p0_T->data[0] << ", " << 
  //                          p0_T->data[1] << ", " <<
  //                          p0_T->data[2] << ")" << std::endl;
  double alpha = 1;
  double beta =  0;

  gsl_blas_dgemv(CblasNoTrans, alpha, Att, p0_T, beta, p1);
  //std::cout << "p1= (" << p1->data[0] << ", " << 
  //                        p1->data[1] << ", " <<
  //                        p1->data[2] << ")" << std::endl;
  gsl_vector_free(T);
  gsl_vector_free(p0_T); 
  return p1;
}

//--------------------------- PosInBF2PosInParentFrame
gsl_vector* ReferenceFrame::PosInBF2PosInParentFrame(gsl_vector *p1, unsigned int L) const
{
  // Check the correctness of the input parameters
  // Check dimension of p1
  if (!CheckVectorLength(p1, REFFRAME_DIMS))
  {
    throw("Position vector with wrong size or undefined in PosInBF2PosInParentFrame.");
  }
  // Check the correctness of L
  if (mSolidary)  // Solidary frame
  {
    L=0; // to ignore the instant position
  }
  else
  {
    if (!CheckTimeInstant(this, L))
    {
      throw("Requesting to compute position at non-existent instant of time. L too large."); 
    }
  }
  gsl_vector *p0 = gsl_vector_alloc(REFFRAME_DIMS);
  // copy mPos vector into PosInParent
  gsl_matrix_get_col(p0, mPos, L);
  gsl_matrix_view Att = gsl_matrix_submatrix(mAtt,0,L*REFFRAME_DIMS, REFFRAME_DIMS, REFFRAME_DIMS); 
  double alpha = 1;
  double beta  = 0;
  //std::cout << "L=" << L << std::endl;
  //std::cout << "p1: (" << p1->data[0] << "; " <<
  //                        p1->data[1] << "; " <<
  //                        p1->data[2] << ")" << std::endl;
  //std::cout << "p0: (" << p0->data[0] << "; " <<
  //                        p0->data[1] << "; " <<
  //                        p0->data[2] << ")" << std::endl;
  gsl_vector *x = gsl_vector_alloc(REFFRAME_DIMS);
  gsl_vector_memcpy(x,p1);
  gsl_vector *aux;
  aux = GetPosition(L);
  /*std::cout << "::::: (" << aux->data[0] << ", " 
                         << aux->data[1] << ", "
                         << aux->data[2] << ")" << std::endl;*/
  gsl_vector_add(x,aux);
  gsl_vector_free(aux);
  gsl_blas_dgemv(CblasTrans, alpha, &Att.matrix, x, beta, p0);
  /*std::cout << "p0: (" << p0->data[0] << "; " <<
                          p0->data[1] << "; " <<
                          p0->data[2] << ")" << std::endl;*/
  gsl_vector_free(x);

  return p0;
}

//---------------------------PosInBF2PosInECEF
gsl_vector* ReferenceFrame::PosInBF2PosInECEF(gsl_vector *posBF, unsigned int L) const
{
  bool b=false;
  
  // Check input parameters
  // Check that the position vector has the correct size
  b = (REFFRAME_DIMS == posBF->size);
  if (!b)
  {
    throw("PosInBF2PosInECEF: Cannot compute position. Position vector has wrong size.");
  }

  // Check that the reference frame can be evaluated at L
  b = CheckTimeInstant(this,L);
  if (!b)
  {
    throw("PosInBF2PosInECEF: Cannot compute position. Unexisting time instant.");
  }

  // Check that the reference ends at ECEF
  b = CheckEndsAtECEF(this);
  if (!b)
  {
    throw("PosInBF2PosInECEF: Cannot compute position. ECEF is not the main reference frame.");
  }

  gsl_vector *pos0;  // The position in the parent frame (towards ECEF)
  gsl_vector *pos1;  // The position in the current frame
  pos1 = gsl_vector_alloc(REFFRAME_DIMS);
  gsl_vector_memcpy(pos1,posBF);

  const ReferenceFrame *rf;
  rf = this;
  bool IsECEF=false;
  IsECEF = ( NULL == rf->mMyRefFrame );

  
  while (!IsECEF)
  {
    //rf->Print();
    /*std::cout << std::endl;
    std::cout << "Frame: " << rf << "   " << rf->mName << std::endl;
    std::cout << " referred: " << rf->mMyRefFrame << std::endl;*/
    /*std::cout << "input position= (" << pos1->data[0] << ", " << 
                                        pos1->data[1] << ", " <<
                                        pos1->data[2] << ")" << std::endl;*/
    pos0 = rf->PosInBF2PosInParentFrame(pos1,L);
    /*std::cout << "output position= (" << pos0->data[0] << ", " << 
                                         pos0->data[1] << ", " <<
                                         pos0->data[2] << ")" << std::endl;*/
   
    //std::cout << "copying" << std::endl; 
    gsl_vector_memcpy(pos1,pos0);
    gsl_vector_free(pos0);
    //std::cout << "getting next reference frame" << std::endl;
    rf = rf->mMyRefFrame;
    IsECEF = (NULL == rf->mMyRefFrame);
  }
  return pos1;
}

//---------------------------PosInECEF2PosInBF
gsl_vector* ReferenceFrame::PosInECEF2PosInBF(gsl_vector *posECEF, unsigned int L) const
{
  bool b=false;
  
  // Check input parameters
  // Check that the position vector has the correct size
  b = (REFFRAME_DIMS == posECEF->size);
  if (!b)
  {
    throw("PosInECEF2PosInBF: Cannot compute position. Position vector has wrong size.");
  }

  // Check that the reference frame can be evaluated at L
  b = CheckTimeInstant(this,L);
  if (!b)
  {
    throw("PosInECEF2PosInBF: Cannot compute position. Unexisting time instant.");
  }

  // Check that the reference ends at ECEF
  b = CheckEndsAtECEF(this);
  if (!b)
  {
    throw("PosInECEF2PosInBF: Cannot compute position. ECEF is not the main reference frame.");
  }

  //gsl_vector *pos0;  // The position in the parent frame (towards ECEF)
  //gsl_vector *pos1;  // The position in the current frame
  //pos1 = gsl_vector_alloc(REFFRAME_DIMS);
  //gsl_vector_memcpy(pos1,posBF);

  const ReferenceFrame *rf;
  rf = this;
  bool IsECEF=false;
  IsECEF = ( NULL == rf->mMyRefFrame );

  gsl_vector *pos0; // Position expressed in the parent frame
  gsl_vector *pos1; // Position expressed in the current frame

  if (!IsECEF)
  {
    //std::cout << "It is not ECEF." << std::endl;
    //this->Print();
    /*std::cout << std::endl;
    std::cout << "Frame: " << rf << "   " << rf->mName << std::endl;
    std::cout << " referred: " << rf->mMyRefFrame << std::endl;*/
    /*std::cout << "input position= (" << pos1->data[0] << ", " << 
                                        pos1->data[1] << ", " <<
                                        pos1->data[2] << ")" << std::endl;*/
    // Get the position expressed in the parent frame
    pos0 = (rf->mMyRefFrame)->PosInECEF2PosInBF(posECEF,L);
    //std::cout << "position in parent frame= (" << pos0->data[0] << ", " << 
    //                                              pos0->data[1] << ", " <<
    //                                              pos0->data[2] << ")" << std::endl;
    pos1 = PosInParentFrame2PosInBF(pos0,L);
    //std::cout << "position in my frame= (" << pos1->data[0] << ", " << 
    //                                          pos1->data[1] << ", " <<
    //                                          pos1->data[2] << ")" << std::endl;
    gsl_vector_free(pos0);
    return pos1;
  }
  else
  {
    //std::cout << "I am ECEF !" << std::endl;
    //this->Print();
    pos1 = gsl_vector_alloc(REFFRAME_DIMS);
    gsl_vector_memcpy(pos1,posECEF);
    //std::cout << "position in ECEF= (" << pos1->data[0] << ", " << 
    //                                      pos1->data[1] << ", " <<
    //                                      pos1->data[2] << ")" << std::endl;
    
    return pos1;
  }
}


int ReferenceFrame::testfun(int a)
{
  return a*2;
}

//---------------------------VelInBF2VelInECEF

//---------------------------VelInBF2VelInParentFrame
//---------------------------VelInECEF2VelInBF














#endif  // RERERENCEFRAME_CPP_

