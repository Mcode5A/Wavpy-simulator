#ifndef     REFERENCEFRAME_HPP_
#define     REFERENCEFRAME_HPP_

/**
    \file   ReferenceFrame.hpp
    \brief  Header file for the ReferenceFrame class.
*/

/**
    \brief  The ReferenceFrame class allows to define positions, orientation
            and velocities of reference frames respect to other reference
            frames.
*/

//#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

//#include <iostream>
#include <list>
#include <string>

#define REFFRAME_DIMS 3


class ReferenceFrame
{
#ifdef TEST_REFERENCEFRAME
  public:
#else
  private:
#endif
    std::string     mName;
    ReferenceFrame  *mMyRefFrame=NULL;
    bool            mSolidary;
    bool            mAllReferredSolidary;  // si tots els seus referenciats son solidaris
    unsigned int    mL;                    // nombre de punts
    //unsigned int    mLmax;                 // nombre de punts maxim dels RefFrame referenciats
    gsl_matrix      *mPos=NULL;            // posició respecte el marc original (posicio)
    gsl_matrix      *mAtt=NULL;            // rotació respecte el marc original (attitude)
    gsl_matrix      *mVel=NULL;            // velocitat

    static bool CheckPositionDimensions(gsl_matrix *Pos, unsigned int L);
    static bool CheckAttitudeDimensions(gsl_matrix *Pos, unsigned int L);
    static bool CheckVelocityDimensions(gsl_matrix *Pos, unsigned int L);
    static bool CheckVectorLength(gsl_vector *v, unsigned int len);
    static bool CheckTimeInstant(const ReferenceFrame *rf, unsigned int L);
    static bool CheckEndsAtECEF(const ReferenceFrame *rf, std::list<const ReferenceFrame*> *mylist=NULL);
  public:
    ReferenceFrame();  // The main reference frame
    ReferenceFrame(ReferenceFrame *RefFrame, bool Solidary=true, gsl_matrix *Pos=NULL, gsl_matrix *Att=NULL, gsl_matrix *Vel=NULL, std::string Name="");
    ~ReferenceFrame();
    bool                  IsSolidary() const;
    bool                  IsMainReferenceFrame() const;
    const ReferenceFrame* GetMyReferenceFrame() const;
    unsigned int          GetLength() const;
    gsl_vector*           GetPosition(unsigned int L=0) const;
    gsl_matrix*           GetAttitude(unsigned int L=0) const;
    gsl_vector*           GetVelocity(unsigned int L=0) const;

    void                  Print() const;

    gsl_vector*  PosInParentFrame2PosInBF(gsl_vector *PosParent, unsigned int L=0) const;
    gsl_vector*  PosInBF2PosInParentFrame(gsl_vector *PosBF, unsigned int L=0) const;
    gsl_vector*  PosInBF2PosInECEF(gsl_vector *PosBF, unsigned int L=0) const;
    gsl_vector*  PosInECEF2PosInBF(gsl_vector *PosECEF, unsigned int L=0) const;
    int          testfun(int a);
    //gsl_vector*     VelInBF2VelInECEF(gsl_vector *VelBF, unsigned int L=0) const;
    //gsl_vector*     VelInBF2VelInParentFrame(gsl_vector *VelBF, unsigned int L=0) const;
    //gsl_vector*     VelInECEF2VelInBF(gsl_vector *PosECEF, unsigned int L=0) const

    /*
    */
};


#endif  // RERERENCEFRAME_HPP_

