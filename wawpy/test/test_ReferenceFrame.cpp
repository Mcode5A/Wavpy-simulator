#define BOOST_TEST_MODULE test_ReferenceFrame
#include <boost/test/unit_test.hpp>

#define TEST_REFERENCEFRAME
#include <ReferenceFrame.hpp>
#include <iostream>


BOOST_AUTO_TEST_SUITE (test_ReferenceFrame)


BOOST_AUTO_TEST_CASE (test_ReferenceFrame_constructor)
{

  unsigned int L=2;
  bool val;
  // =================== Check Position Size
  {
    gsl_matrix *Pos = gsl_matrix_alloc(REFFRAME_DIMS, L);
    //--------------------- When everything works fine
    val = ReferenceFrame::CheckPositionDimensions(Pos, L);
    BOOST_CHECK_MESSAGE(val == true, "Result 1.1: " << val);
    //--------------------- When the number of columns is wrong
    val = ReferenceFrame::CheckPositionDimensions(Pos, L+1);
    BOOST_CHECK_MESSAGE( val == false, "Result 1.2: " << val);
    gsl_matrix_free(Pos);
    //--------------------- When the number of rows is wrong
    Pos = gsl_matrix_alloc(REFFRAME_DIMS+1, L);
    val = ReferenceFrame::CheckPositionDimensions(Pos, L);
    BOOST_CHECK_MESSAGE(val == false, "Result 1.3:" << val);
    gsl_matrix_free(Pos);
    Pos = NULL;
  }
  // =================== Check Attitude Size
  {
    gsl_matrix *Att = gsl_matrix_alloc(REFFRAME_DIMS, L*REFFRAME_DIMS);
    //--------------------- When everything works fine
    val = ReferenceFrame::CheckAttitudeDimensions(Att, L);
    BOOST_CHECK_MESSAGE(val == true, "Result 1.4: " << val);
    //--------------------- When the number of columns is wrong
    val = ReferenceFrame::CheckPositionDimensions(Att, L+1);
    BOOST_CHECK_MESSAGE( val == false, "Result 1.5: " << val);
    gsl_matrix_free(Att);
    //--------------------- When the number of rows is wrong
    Att = gsl_matrix_alloc(REFFRAME_DIMS+1, L*REFFRAME_DIMS);
    val = ReferenceFrame::CheckAttitudeDimensions(Att, L);
    BOOST_CHECK_MESSAGE(val == false, "Result 1.6:" << val);
    gsl_matrix_free(Att);
    Att = NULL;
  }
  // =================== Check Velocity Size
  {
    gsl_matrix *Vel = gsl_matrix_alloc(REFFRAME_DIMS, L);
    //--------------------- When everything works fine
    val = ReferenceFrame::CheckVelocityDimensions(Vel, L);
    BOOST_CHECK_MESSAGE(val == true, "Result 1.7: " << val);
    //--------------------- When the number of columns is wrong
    val = ReferenceFrame::CheckVelocityDimensions(Vel, L+1);
    BOOST_CHECK_MESSAGE( val == false, "Result 1.8: " << val);
    gsl_matrix_free(Vel);
    //--------------------- When the number of rows is wrong
    Vel = gsl_matrix_alloc(REFFRAME_DIMS+1, L);
    val = ReferenceFrame::CheckVelocityDimensions(Vel, L);
    BOOST_CHECK_MESSAGE(val == false, "Result 1.9:" << val);
    gsl_matrix_free(Vel);
    Vel = NULL;
  }

  // =================== Check constructor without parameters
  // Should construct the ECEF reference frame
  {
    ReferenceFrame *ECEF = new ReferenceFrame();
    const ReferenceFrame *aux;
    bool b=false;
    // ---------------- GetMyReferenceFrame()
    aux = ECEF->GetMyReferenceFrame();
    BOOST_CHECK_MESSAGE ( NULL == aux, "Result 1.10:" << aux );
    // ---------------- IsSolidary()
    b = ECEF->IsSolidary();
    BOOST_CHECK_MESSAGE ( true == b, "Result 1.11:" << b );
    // ---------------- IsMainReferenceFrame()
    b = ECEF->IsMainReferenceFrame();
    BOOST_CHECK_MESSAGE ( true == b, "Result 1.12:" << b );
    delete ECEF;
  }

  // =================== Check Long constructor 1
  // One input parameter = NULL
  // Should construct the ECEF frame 
  {
    ReferenceFrame *ECEF = new ReferenceFrame(NULL);
    const ReferenceFrame *aux;
    bool b=false;
    // ---------------- GetMyReferenceFrame()
    aux = ECEF->GetMyReferenceFrame();
    BOOST_CHECK_MESSAGE ( NULL == aux, "Result 1.13:" << aux );
    // ---------------- IsSolidary()
    b = ECEF->IsSolidary();
    BOOST_CHECK_MESSAGE ( true == b, "Result 1.14:" << b );
    // ---------------- IsMainReferenceFrame()
    b = ECEF->IsMainReferenceFrame();
    BOOST_CHECK_MESSAGE ( true == b, "Result 1.15:" << b );
    delete ECEF;
  }
  // One input parameter != NULL
  {
    ReferenceFrame *ECEF = new ReferenceFrame();
    ReferenceFrame *bf   = new ReferenceFrame(ECEF);
    const ReferenceFrame *aux;
    bool b;
    // ---------------- GetMyReferenceFrame()
    aux = bf->GetMyReferenceFrame();
    BOOST_CHECK_MESSAGE (aux == ECEF, "Result 1.16: EFEC:" << ECEF << "  aux:" << aux );
    // ---------------- IsSolidary()
    b = bf->IsSolidary();
    BOOST_CHECK_MESSAGE ( true==b, "Result 1.17:" << b );
    // ---------------- IsMainReferenceFrame()
    b = bf->IsMainReferenceFrame();

    delete ECEF;
    delete bf;
  }

  // Solidary frame with specified position and attitude
  {
    ReferenceFrame *ECEF = new ReferenceFrame();
    L = 1;
    gsl_matrix *Pos = gsl_matrix_calloc(REFFRAME_DIMS,L);
    {
      gsl_matrix_set(Pos,0,0,6371e3);
      gsl_matrix_set(Pos,1,0,0);
      gsl_matrix_set(Pos,2,0,0);
    }
    gsl_matrix *Att = gsl_matrix_calloc(REFFRAME_DIMS,L*REFFRAME_DIMS);
    {
      gsl_matrix_set(Att,0,1,1);
      gsl_matrix_set(Att,1,0,-1);
      gsl_matrix_set(Att,2,2,1);
    }
    ReferenceFrame *bf   = new ReferenceFrame(ECEF,true,Pos,Att);
    const ReferenceFrame *aux;
    bool b;
    // ---------------- GetMyReferenceFrame()
    aux = bf->GetMyReferenceFrame();
    BOOST_CHECK_MESSAGE (aux == ECEF, "Result 1.18: EFEC:" << ECEF << "  aux:" << aux );
    // ---------------- IsSolidary()
    b = bf->IsSolidary();
    BOOST_CHECK_MESSAGE ( true==b, "Result 1.19:" << b );
    // ---------------- IsMainReferenceFrame()
    b = bf->IsMainReferenceFrame();
    BOOST_CHECK_MESSAGE ( false == b, "Result 1.20:" << b );
    // ---------------- Check values of the position matrix
    for (int k=0;k<REFFRAME_DIMS;k++)
    {
      b = (gsl_matrix_get(Pos,k,0) == gsl_matrix_get(bf->mPos,k,0));
      BOOST_CHECK_MESSAGE (b, "Result 1.21: k=" << k << ": " << b );
    }
    // ---------------- Check values of the attitude matrix
    for (int k=0;k<REFFRAME_DIMS;k++)
    {
      for (int l=0;l<REFFRAME_DIMS;l++)
      {
        b = (gsl_matrix_get(Att,k,l) == gsl_matrix_get(bf->mAtt,k,l));
        BOOST_CHECK_MESSAGE ( b, "Result 1.22: k=" << k << ", l=" << l << ": " << b );
      }
    }
    gsl_matrix_free(Pos);
    gsl_matrix_free(Att);
    delete ECEF;
    delete bf;  
  }

  // Non-solidary frame with specified position, attitude and velocity
  {
    ReferenceFrame *ECEF = new ReferenceFrame();
    L = 2;
    gsl_matrix *Pos = gsl_matrix_calloc(REFFRAME_DIMS,L);
    {
      for (unsigned int l=0;l<L;l++)
      {
        gsl_matrix_set(Pos,0,l,6371e3+l);
        gsl_matrix_set(Pos,1,l,0+l*2);
        gsl_matrix_set(Pos,2,l,0-l);
      }
    }
    gsl_matrix *Att = gsl_matrix_calloc(REFFRAME_DIMS,L*REFFRAME_DIMS);
    {
      for (unsigned int l=0;l<L;l++)
      { // first col
        gsl_matrix_set(Att,0,l*REFFRAME_DIMS,l);
        gsl_matrix_set(Att,1,l*REFFRAME_DIMS,l*2);
        gsl_matrix_set(Att,2,l*REFFRAME_DIMS,l*3);
        // second col
        gsl_matrix_set(Att,0,1+l*REFFRAME_DIMS,-l);
        gsl_matrix_set(Att,1,1+l*REFFRAME_DIMS,-l*2);
        gsl_matrix_set(Att,2,1+l*REFFRAME_DIMS,-l*3);
        // third col
        gsl_matrix_set(Att,0,2+l*REFFRAME_DIMS,10-l);
        gsl_matrix_set(Att,1,2+l*REFFRAME_DIMS,10-l*2);
        gsl_matrix_set(Att,2,2+l*REFFRAME_DIMS,10-l*3);
      }
    }
    gsl_matrix *Vel = gsl_matrix_calloc(REFFRAME_DIMS,L);
    {
      for (unsigned int l=0;l<L;l++)
      {
        gsl_matrix_set(Vel,0,l,l);
        gsl_matrix_set(Vel,1,l,100+l);
        gsl_matrix_set(Vel,2,l,50-l);
      }
    }
    ReferenceFrame *bf   = new ReferenceFrame(ECEF,false,Pos,Att,Vel);
    bool b;
    // ---------------- Check values of the position matrix
    for (unsigned int l=0;l<L;l++)
    {
      for (int k=0;k<REFFRAME_DIMS;k++)
      {
        b = (gsl_matrix_get(Pos,k,l) == gsl_matrix_get(bf->mPos,k,l));
        BOOST_CHECK_MESSAGE (b, "Result 1.23: k=" << k << ", l=" << l << ": " << b );
      }
    }
    // ---------------- Check values of the attitude matrix
    for (unsigned int l=0;l<L;l++)
    {
      for (int k=0;k<REFFRAME_DIMS;k++)
      {
        for (int m=0;m<REFFRAME_DIMS;m++)
        {
          b = (gsl_matrix_get(Att,k,m+l*REFFRAME_DIMS) == gsl_matrix_get(bf->mAtt,k,m+l*REFFRAME_DIMS));
          BOOST_CHECK_MESSAGE ( b, "Result 1.24: l=" << l << ", k=" << k << ", m=" << m << ": " << b );
        }
      }
    }
    gsl_matrix_free(Pos);
    gsl_matrix_free(Att);
    gsl_matrix_free(Vel);
    delete ECEF;
    delete bf;  
  }

  // Non-solidary frame with  specified position, attitude and velocity (wrong-sized matrices)
  {
    ReferenceFrame *ECEF = new ReferenceFrame();
    L = 2;
    gsl_matrix *Pos   = gsl_matrix_calloc(REFFRAME_DIMS,L);
    gsl_matrix *Att   = gsl_matrix_calloc(REFFRAME_DIMS,L*REFFRAME_DIMS);
    gsl_matrix *Att_w = gsl_matrix_calloc(REFFRAME_DIMS,L*REFFRAME_DIMS+1);
    gsl_matrix *Vel   = gsl_matrix_calloc(REFFRAME_DIMS,L);
    gsl_matrix *Vel_w = gsl_matrix_calloc(REFFRAME_DIMS,L+1);

    // Matrices with correct sizes
    BOOST_CHECK_NO_THROW (ReferenceFrame *bf   = new ReferenceFrame(ECEF,false,Pos,Att,Vel); delete bf;);
    BOOST_CHECK_THROW (ReferenceFrame(ECEF,false,Pos,Att_w,Vel), const char*);
    BOOST_CHECK_THROW (ReferenceFrame(ECEF,false,Pos,Att,Vel_w), const char*);
    BOOST_CHECK_THROW (ReferenceFrame(ECEF,false,Pos,Att_w,Vel_w), const char*);
    gsl_matrix_free(Pos);
    gsl_matrix_free(Att);
    gsl_matrix_free(Att_w);
    gsl_matrix_free(Vel);
    gsl_matrix_free(Vel_w);
    delete ECEF;
  }

  // Check that mAllSolidary is set up correctly
  {
    //std::cout << "ECEF" << std::endl;
    ReferenceFrame *ECEF = new ReferenceFrame();
    L=2;
    gsl_matrix *Pos   = gsl_matrix_calloc(REFFRAME_DIMS,L);
    gsl_matrix *Att   = gsl_matrix_calloc(REFFRAME_DIMS,L*REFFRAME_DIMS);
    gsl_matrix *Vel   = gsl_matrix_calloc(REFFRAME_DIMS,L);
    //std::cout << "RF1" << std::endl;
    ReferenceFrame *RF1 = new ReferenceFrame(ECEF,false,Pos,Att,Vel);
    //std::cout << "RF2" << std::endl;
    ReferenceFrame *RF2 = new ReferenceFrame(RF1,true);
    //std::cout << "RF3" << std::endl;
    ReferenceFrame *RF3 = new ReferenceFrame(RF2,true);

    bool b;
    // ECEF
    b = (ECEF->mSolidary == true);
    BOOST_CHECK_MESSAGE ( b, "Result 1.25: " << b );
    b = (ECEF->mAllReferredSolidary == true);
    BOOST_CHECK_MESSAGE ( b, "Result 1.26: " << b );

    // RF1 -> ECEF
    b = (RF1->mSolidary == false);
    BOOST_CHECK_MESSAGE ( b, "Result 1.27: " << b );
    b = (RF1->mAllReferredSolidary == true);
    BOOST_CHECK_MESSAGE ( b, "Result 1.28: " << b );
    
    // RF2 -> RF1 -> ECEF
    b = (RF2->mSolidary == true);
    BOOST_CHECK_MESSAGE ( b, "Result 1.29: " << b );
    b = (RF2->mAllReferredSolidary == false);
    BOOST_CHECK_MESSAGE ( b, "Result 1.30: " << b );
    
    delete ECEF;
    delete RF1;
    delete RF2;
    delete RF3;
  }

  // Check that mLmax, mSolidary, mAllReferredSolidary are set up correcty
  {
    //std::cout << "ECEF" << std::endl;
    ReferenceFrame *ECEF = new ReferenceFrame();
    L=2;
    gsl_matrix *Pos   = gsl_matrix_calloc(REFFRAME_DIMS,L);
    gsl_matrix *Att   = gsl_matrix_calloc(REFFRAME_DIMS,L*REFFRAME_DIMS);
    gsl_matrix *Vel   = gsl_matrix_calloc(REFFRAME_DIMS,L);
    gsl_matrix *Pos_  = gsl_matrix_calloc(REFFRAME_DIMS,L+1);
    gsl_matrix *Att_  = gsl_matrix_calloc(REFFRAME_DIMS,(L+1)*REFFRAME_DIMS);
    gsl_matrix *Vel_  = gsl_matrix_calloc(REFFRAME_DIMS,L+1);
    bool b;

    // ECEF
    b = (ECEF->mL == 1);
    BOOST_CHECK_MESSAGE ( b, "Result 1.31: " << b );
    b = (true == ECEF->mSolidary);
    BOOST_CHECK_MESSAGE ( b, "Result 1.32: " << b );
    b = (true == ECEF->mAllReferredSolidary);
    BOOST_CHECK_MESSAGE ( b, "Result 1.33: " << b );

    // RF1 -> ECEF
    //std::cout << "RF1" << std::endl;
    ReferenceFrame *RF1 = new ReferenceFrame(ECEF,true);
    b = (RF1->mL == 1);
    BOOST_CHECK_MESSAGE ( b, "Result 1.34: " << b );
    b = (true == RF1->mSolidary);
    BOOST_CHECK_MESSAGE ( b, "Result 1.35: " << b );
    b = (true == RF1->mAllReferredSolidary);
    BOOST_CHECK_MESSAGE ( b, "Result 1.36: " << b );

    // RF2 -> RF1 -> ECEF
    //std::cout << "RF2" << std::endl;
    ReferenceFrame *RF2 = new ReferenceFrame(RF1,false,Pos,Att,Vel);
    b = (RF2->mL == L);
    BOOST_CHECK_MESSAGE ( b, "Result 1.37: " << b );
    b = (false == RF2->mSolidary);
    BOOST_CHECK_MESSAGE ( b, "Result 1.38: " << b );
    b = (true == RF2->mAllReferredSolidary);
    BOOST_CHECK_MESSAGE ( b, "Result 1.39: " << b );


    // RF3 -> RF2 -> RF1 -> ECEF
    //std::cout << "RF3" << std::endl;
    // aqui hauria de petar
    BOOST_CHECK_THROW ( ReferenceFrame *RF3 = new ReferenceFrame(RF2,false,Pos_,Att_,Vel_); delete RF3, const char*);


    // RF4 -> -> RF1 -> ECEF
    //std::cout << "RF4" << std::endl;
    ReferenceFrame *RF4 = new ReferenceFrame(RF1,false,Pos_,Att_,Vel_);
    b = (RF4->mL == (L+1));
    BOOST_CHECK_MESSAGE ( b, "Result 1.40: " << b );
    b = (RF4->mL == (L+1));
    BOOST_CHECK_MESSAGE ( b, "Result 1.41: " << b );


    gsl_matrix_free(Pos);    
    gsl_matrix_free(Att);    
    gsl_matrix_free(Vel);    
    gsl_matrix_free(Pos_);    
    gsl_matrix_free(Att_);    
    gsl_matrix_free(Vel_);    
    delete ECEF;
    delete RF1;
    delete RF2;
    delete RF4;

  }
}




BOOST_AUTO_TEST_CASE (test_ReferenceFrame_get_methods)
{
  {
    ReferenceFrame *ECEF = new ReferenceFrame();

    unsigned int L = 3;
    gsl_matrix *Pos = gsl_matrix_calloc(REFFRAME_DIMS,L);
    gsl_matrix *Att = gsl_matrix_calloc(REFFRAME_DIMS,L*REFFRAME_DIMS);
    gsl_matrix *Vel = gsl_matrix_calloc(REFFRAME_DIMS,L);

    for (unsigned int r=0;r<REFFRAME_DIMS;r++)
    {
      for (unsigned int l=0;l<L;l++)
      {
        gsl_matrix_set(Pos,r,l,r+l);
        gsl_matrix_set(Att,r,l,r+l);
        gsl_matrix_set(Att,r,l+1,r+l+1);
        gsl_matrix_set(Att,r,l+2,r+l+2);
        gsl_matrix_set(Vel,r,l, (r+l)*2);
      }
    }
    ReferenceFrame *Bf = new ReferenceFrame(ECEF, false, Pos, Att, Vel);

    bool b;
    // Extract Position from solidary frame
    gsl_vector *v;
    v = ECEF->GetPosition();
    b = gsl_vector_isnull(v);
    //std::cout << "ECEF[0]:" << gsl_vector_get(v,0) << std::endl <<
    //             "ECEF[1]:" << gsl_vector_get(v,1) << std::endl <<
    //             "ECEF[2]:" << gsl_vector_get(v,2) << std::endl << std::endl;
    BOOST_CHECK_MESSAGE( b, "Result 2.1: " << b);
    // Extract Attitude from solidary frame
    gsl_matrix *m;
    m = ECEF->GetAttitude();
    //std::cout << "m: " << m << std::endl;
    //std::cout << "m: " << m->size1 << "m" << m->size2 << std::endl;

    gsl_matrix *I = gsl_matrix_alloc(REFFRAME_DIMS, REFFRAME_DIMS);
    gsl_matrix_set_identity(I); 
    b = gsl_matrix_equal(I,m);
    BOOST_CHECK_MESSAGE( b, "Result 2.2:" << b);
    // Extract Velocity from solidary frame
    v = ECEF->GetVelocity();
    b = gsl_vector_isnull(v);
    BOOST_CHECK_MESSAGE( b, "Result 2.3: " << b);

    // Extract MyReferenceFrame
    const ReferenceFrame *aux;
    aux = Bf->GetMyReferenceFrame();
    b = (ECEF == aux);
    BOOST_CHECK_MESSAGE( b, "Result 2.4: " << b );
    // Extract Length
    unsigned int l;
    l = Bf->GetLength();
    b = (L == l);
    BOOST_CHECK_MESSAGE( b, "Result 2.5: " << b );    
    // Extract Position from non-solidary frame
    for (unsigned int l=0;l<L; l++)
    {
      gsl_vector *p = Bf->GetPosition(l);
      gsl_matrix *a = Bf->GetAttitude(l);
      gsl_vector *v = Bf->GetVelocity(l);
    for (unsigned int r=0;r<REFFRAME_DIMS;r++)
      { // Position
        b = (gsl_vector_get(p,r) == gsl_matrix_get(Pos,r,l));
        //std::cout << "v(" << r << ")=" << gsl_vector_get(v,r) << std::endl;  
        BOOST_CHECK_MESSAGE (b, "Result 2.5 (Position): " << b << " r=" << r << ",  l=" << l);
        // Attitude
        for (unsigned int j=0;j<REFFRAME_DIMS;j++)
        {
          b = (gsl_matrix_get(a,r,j) == gsl_matrix_get(Att,r,l*REFFRAME_DIMS + j));
          BOOST_CHECK_MESSAGE (b, "Result 2.6 (Attitude): " << b << " r=" << r << ", l=" << l << "  j=" << j );
        }
        // Velocity
        b = (gsl_vector_get(v,r) == gsl_matrix_get(Vel,r,l));
        BOOST_CHECK_MESSAGE (b, "Result 2.7 (Velocity): " << b << " r=" << r << ",  l=" << l);
      }
    gsl_vector_free(p);    
    gsl_matrix_free(a);
    gsl_vector_free(v);    
    }
  }
}


BOOST_AUTO_TEST_CASE (test_ReferenceFrame_CheckVectorLength)
{
  gsl_vector *v=NULL;
  bool b;
  // Check when providing a NULL pointer as a gsl_vector
  for (unsigned int l=0;l<4;l++)
  {
    b = ReferenceFrame::CheckVectorLength(v,l);
    BOOST_CHECK_MESSAGE( !b, "Result 3.1 (CheckVectorLength): " << b );
  }
  v = gsl_vector_alloc(REFFRAME_DIMS);
  for (unsigned int l=0;l<REFFRAME_DIMS+2;l++)
  {
    b = ReferenceFrame::CheckVectorLength(v,l);
    if (REFFRAME_DIMS==l)
    {
      BOOST_CHECK_MESSAGE( b, "Result 3.2 (CheckVectorLength): " << b << " l=" << l);
    }
    else
    {
      BOOST_CHECK_MESSAGE( !b, "Result 3.3 (CheckVectorLength): " << b << " l=" << l );
    }
  }
  gsl_vector_free(v);
}

BOOST_AUTO_TEST_CASE (test_ReferenceFrame_CheckTimeInstant)
{
  ReferenceFrame *rf=NULL;
  bool b;
  // Check when providing a NULL pointer as a ReferenceFrame
  for (unsigned int l=0;l<4;l++)
  {
    b = ReferenceFrame::CheckTimeInstant(rf,l);
    BOOST_CHECK_MESSAGE( !b, "Result 4.1 (CheckTimeInstant): " << b );
  }

  // Checking ECEF
  ReferenceFrame *ECEF;
  ECEF = new ReferenceFrame();
  for (unsigned int l=0;l<3;l++)
  {
    b = ReferenceFrame::CheckTimeInstant(ECEF,l);
    BOOST_CHECK_MESSAGE( b, "Result 4.2 (CheckTimeInstant): " << b << " l=" << l);
  }

  // Checking RF1 (solidary) -> ECEF
  ReferenceFrame *RF1;
  double p1[] = {1, -3, 4};
  double a1[] = { 0, 1, 0,
                -1, 0, 0,
                 0, 0, 1};
  double v1[] = {0, 0, 0};
  gsl_matrix_view pos1 = gsl_matrix_view_array(p1,REFFRAME_DIMS,1);
  gsl_matrix_view att1 = gsl_matrix_view_array(a1,REFFRAME_DIMS,1*3);
  gsl_matrix_view vel1 = gsl_matrix_view_array(v1,REFFRAME_DIMS,1);

  RF1 = new ReferenceFrame(ECEF,true,&pos1.matrix,&att1.matrix,&vel1.matrix);
  for (unsigned int l=0;l<3;l++)
  {
    b = ReferenceFrame::CheckTimeInstant(RF1,l);
    BOOST_CHECK_MESSAGE( b, "Result 4.3 (CheckTimeInstant): " << b << " l=" << l);
  }
 
  // Checking RF2 (non-solidary) -> RF1 (solidary) -> ECEF
  ReferenceFrame *RF2;
  double p2[] = {1, -3, 4, 2, -3, 4, 3, -3, 4};
  double a2[] = { 0, 1, 0,  0, 1, 0,  0, 1, 0,
                 -1, 0, 0, -1, 0, 0, -1, 0, 0,
                  0, 0, 1,  0, 0, 1,  0, 0, 1 };
  double v2[] = {1, 0, 0, 1, 0, 0, 1, 0, 0};
  gsl_matrix_view pos2 = gsl_matrix_view_array(p2,REFFRAME_DIMS,3);
  gsl_matrix_view att2 = gsl_matrix_view_array(a2,REFFRAME_DIMS,3*3);
  gsl_matrix_view vel2 = gsl_matrix_view_array(v2,REFFRAME_DIMS,3);

  RF2 = new ReferenceFrame(RF1,false,&pos2.matrix,&att2.matrix,&vel2.matrix);
  for (unsigned int l=0;l<5;l++)
  {
    b = ReferenceFrame::CheckTimeInstant(RF2,l);
    if (l<3)
    {
      BOOST_CHECK_MESSAGE( b, "Result 4.4 (CheckTimeInstant): " << b << " l=" << l);
    }
    else
    {
      BOOST_CHECK_MESSAGE( !b, "Result 4.5 (CheckTimeInstant): " << b << " l=" << l);
    }
  }
 
  delete(RF2);
  delete(RF1);
  delete(ECEF);
}


BOOST_AUTO_TEST_CASE (test_ReferenceFrame_PosInBF2PosInParentFrame)
{
  unsigned int L = 3;
  ReferenceFrame *ECEF = new ReferenceFrame();
  gsl_matrix *Pos = gsl_matrix_calloc(REFFRAME_DIMS,L);
  for (unsigned int l=0;l<L;l++)
  {
    gsl_matrix_set(Pos,0,l,1*l);
    gsl_matrix_set(Pos,1,l,2*l);
    gsl_matrix_set(Pos,2,l,3*l);
  }
  gsl_matrix *Att = gsl_matrix_calloc(REFFRAME_DIMS,L*REFFRAME_DIMS);
  for (unsigned int l=0;l<L;l++) // Attitude the same as ECEF
  {
    // first row
    gsl_matrix_set(Att,0,l*REFFRAME_DIMS+0, 1);
    gsl_matrix_set(Att,0,l*REFFRAME_DIMS+1, 0);
    gsl_matrix_set(Att,0,l*REFFRAME_DIMS+2, 0);
    // second row
    gsl_matrix_set(Att,1,l*REFFRAME_DIMS+0, 0);
    gsl_matrix_set(Att,1,l*REFFRAME_DIMS+1, 1);
    gsl_matrix_set(Att,1,l*REFFRAME_DIMS+2, 0);
    // second row
    gsl_matrix_set(Att,2,l*REFFRAME_DIMS+0, 0);
    gsl_matrix_set(Att,2,l*REFFRAME_DIMS+1, 0);
    gsl_matrix_set(Att,2,l*REFFRAME_DIMS+2, 1);
  }
  gsl_matrix *Vel = gsl_matrix_calloc(REFFRAME_DIMS,L);
  ReferenceFrame *BF  = new ReferenceFrame(ECEF,false,Pos,Att,Vel);


  // Chosing p1=(0,0,0) at the BF
  {
    gsl_vector *p1 = gsl_vector_calloc(REFFRAME_DIMS);
    double result[] = { 0 ,1, 2,
                        0, 2, 4,
                        0, 3, 6};
    for (unsigned int l=0;l<L;l++)
    {
      //std::cout << std::endl << "l:" << l << std::endl;
      gsl_vector *p0 = BF->PosInBF2PosInParentFrame(p1,l);
      //std::cout << "p0: ";
      for (unsigned int r=0;r<REFFRAME_DIMS;r++)
      {
        bool b;
        b = ( gsl_vector_get(p0,r) == result[r*REFFRAME_DIMS+l] );
        BOOST_CHECK_MESSAGE ( b, "Result 5.1 (PosInBF2PosInParentFrame): " << b << " (r, l)=(" << r << ", " << l << ")" );
        //std::cout  << gsl_vector_get(p0,r) << ", ";      //result = result + gsl_vector_get(PosInParent,r);
      }
      //std::cout << std::endl;
      gsl_vector_free(p0);
    }
    gsl_vector_free(p1);
  }

  // Chosing p1=(1,1,1) at the BF
  {
    gsl_vector *p1 = gsl_vector_calloc(REFFRAME_DIMS);
    gsl_vector_set_all(p1,1);
    double result[] = { 1, 2, 3,
                        1, 3, 5,
                        1, 4, 7};
    for (unsigned int l=0;l<L;l++)
    {
      //std::cout << std::endl << "l:" << l << std::endl;
      gsl_vector *p0 = BF->PosInBF2PosInParentFrame(p1,l);
      //std::cout << "p0: ";
      for (unsigned int r=0;r<REFFRAME_DIMS;r++)
      {
        bool b;
        b = ( gsl_vector_get(p0,r) == result[r*REFFRAME_DIMS+l] );
        BOOST_CHECK_MESSAGE ( b, "Result 5.2 (PosInBF2PosInParentFrame): " << b << " (r, l)=(" << r << ", " << l << ") Result: " << gsl_vector_get(p0,r) << ", but expected: " << result[r*REFFRAME_DIMS+l] );
        //std::cout  << gsl_vector_get(p0,r) << ", ";      //result = result + gsl_vector_get(PosInParent,r);
      }
      //std::cout << std::endl;
      gsl_vector_free(p0);
    }
    gsl_vector_free(p1);
  }

  // Chosing (0,0,0) at the BF
  // And a rotated Attitude matrix 
  //   [0   1  0]
  //   [-1  0  0]
  //   [0   0  1]
  for (unsigned int l=0;l<L;l++)
  {
    // first row
    gsl_matrix_set(Att,0,l*REFFRAME_DIMS+0, 0);
    gsl_matrix_set(Att,0,l*REFFRAME_DIMS+1, 1);
    gsl_matrix_set(Att,0,l*REFFRAME_DIMS+2, 0);
    // second row
    gsl_matrix_set(Att,1,l*REFFRAME_DIMS+0,-1);
    gsl_matrix_set(Att,1,l*REFFRAME_DIMS+1, 0);
    gsl_matrix_set(Att,1,l*REFFRAME_DIMS+2, 0);
    // second row
    gsl_matrix_set(Att,2,l*REFFRAME_DIMS+0, 0);
    gsl_matrix_set(Att,2,l*REFFRAME_DIMS+1, 0);
    gsl_matrix_set(Att,2,l*REFFRAME_DIMS+2, 1);
  }
  delete BF;
  BF = new ReferenceFrame(ECEF,false,Pos,Att,Vel);
  {
    gsl_vector *p1 = gsl_vector_calloc(REFFRAME_DIMS);
    gsl_vector_set_all(p1,1);
    double result[] = {-1,-3,-5,
                        1, 2, 3,
                        1, 4, 7};
    for (unsigned int l=0;l<L;l++)
    {
      //std::cout << std::endl << "l:" << l << std::endl;
      gsl_vector *p0 = BF->PosInBF2PosInParentFrame(p1,l);
      //std::cout << "p0: ";
      for (unsigned int r=0;r<REFFRAME_DIMS;r++)
      {
        bool b;
        b = ( gsl_vector_get(p0,r) == result[r*REFFRAME_DIMS+l] );
        BOOST_CHECK_MESSAGE ( b, "Result 5.3 (PosInBF2PosInParentFrame): " << b << " (r, l)=(" << r << ", " << l << ") Result: " << gsl_vector_get(p0,r) << ", but expected: " << result[r*REFFRAME_DIMS+l] );
        //std::cout  << gsl_vector_get(p0,r) << ", ";      //result = result + gsl_vector_get(PosInParent,r);
      }
      //std::cout << std::endl;
      gsl_vector_free(p0);
    }
    gsl_vector_free(p1);
  }

  gsl_matrix_free(Pos);
  gsl_matrix_free(Att);
  gsl_matrix_free(Vel);
  delete(BF);
  delete(ECEF);
}


BOOST_AUTO_TEST_CASE (test_ReferenceFrame_PosInParentFrame2PosInBF)
{
  unsigned int L = 3;
  ReferenceFrame *ECEF = new ReferenceFrame();
  // RF1
  double pos[] = {0, 1, 2,
                  0, 2, 4,
                  0, 3, 6 };
  double att[] = {1, 0, 0,   0, 1, 0,  1, 0, 0,
                  0, 1, 0,  -1, 0, 0,  0, 0, 1,
                  0, 0, 1,   0, 0, 1,  0,-1, 0 };
  double vel[] = {0, 0, 0,
                  0, 0, 0,
                  0, 0, 0 };
  gsl_matrix_view Pos = gsl_matrix_view_array(pos, REFFRAME_DIMS, L);
  gsl_matrix_view Att = gsl_matrix_view_array(att, REFFRAME_DIMS, L*REFFRAME_DIMS);
  gsl_matrix_view Vel = gsl_matrix_view_array(vel, REFFRAME_DIMS, L);
  ReferenceFrame *BF = new ReferenceFrame(ECEF,false,&Pos.matrix, &Att.matrix, &Vel.matrix);

  // Chosing p0=(0,0,0) at the parent frame
  {
    double p[] = {0,
                  0,
                  0 };
    gsl_vector_view p0 = gsl_vector_view_array(p, REFFRAME_DIMS);
    double result[] = { 0, -2, -2,
                        0,  1, -6,
                        0, -3,  4};

    for (unsigned int l=0;l<L;l++)
    {
      gsl_vector *p1 = BF->PosInParentFrame2PosInBF(&p0.vector,l);
      for (unsigned int r=0;r<REFFRAME_DIMS;r++)
      {
        bool b;
        b = ( gsl_vector_get(p1,r) == result[r*REFFRAME_DIMS+l] );
        BOOST_CHECK_MESSAGE ( b, "Result 6.1 (PosInParentFrame2PosInBF): " << b << " (r, l)=(" << r << ", " << l << ")" );
      }
      gsl_vector_free(p1);
    }
  }

  // Chosing p0=(1,1,1) at the Parent frame
  {
    double p[] = {1,
                  1,
                  1 };
    gsl_vector_view p0 = gsl_vector_view_array(p, REFFRAME_DIMS);  
    double result[] = { 1,-1,-1,
                        1, 0,-5,
                        1,-2, 3};
    for (unsigned int l=0;l<L;l++)
    {
      gsl_vector *p1 = BF->PosInParentFrame2PosInBF(&p0.vector,l);
      for (unsigned int r=0;r<REFFRAME_DIMS;r++)
      {
        bool b;
        b = ( gsl_vector_get(p1,r) == result[r*REFFRAME_DIMS+l] );
        BOOST_CHECK_MESSAGE ( b, "Result 6.2 (PosInParentFrame2PosInBF): " << b << " (r, l)=(" << r << ", " << l << ") Result: " << gsl_vector_get(p1,r) << ", but expected: " << result[r*REFFRAME_DIMS+l] );
      }
      gsl_vector_free(p1);
    }
  }
  delete(BF);
  delete(ECEF);
}


BOOST_AUTO_TEST_CASE (test_ReferenceFrame_CheckEndsAtECEF)
{
  bool b;
  unsigned int L = 3;
  ReferenceFrame *ECEF = new ReferenceFrame();

  gsl_matrix *Pos1 = gsl_matrix_calloc(REFFRAME_DIMS,L);
  gsl_matrix *Att1 = gsl_matrix_calloc(REFFRAME_DIMS,L*REFFRAME_DIMS);
  gsl_matrix *Vel1 = gsl_matrix_calloc(REFFRAME_DIMS,L);
  ReferenceFrame *RF1 = new ReferenceFrame(ECEF,false,Pos1,Att1,Vel1);

  gsl_matrix *Pos2 = gsl_matrix_calloc(REFFRAME_DIMS,L);
  gsl_matrix *Att2 = gsl_matrix_calloc(REFFRAME_DIMS,L*REFFRAME_DIMS);
  gsl_matrix *Vel2 = gsl_matrix_calloc(REFFRAME_DIMS,L);
  ReferenceFrame *RF2 = new ReferenceFrame(RF1,false,Pos2,Att2,Vel2);

  //ReferenceFrame *RF3 = new ReferenceFrame((ReferenceFrame*)0x1234,false,Pos2,Att2,Vel2);
  // the acutal content for position, attitudes and velocities 
  // is non-important for the test case
  // Check ECEF
  //std::cout << std::endl << std::endl;
  //std::cout << "The ReferenceFrame chain is:" << std::endl;
  //std::cout << "ECEF: " << ECEF << std::endl;
  //std::cout << "------" << std::endl;
  b = ReferenceFrame::CheckEndsAtECEF(ECEF, NULL);
  BOOST_CHECK_MESSAGE( b, "Result 7.1:" << b);

  // Check RF1 -> ECEF
  //std::cout << std::endl << std::endl;
  //std::cout << "The ReferenceFrame chain is:" << std::endl;
  //std::cout << "RF1 : " << RF1 << std::endl;
  //std::cout << "ECEF: " << ECEF << std::endl;
  //std::cout << "------" << std::endl;
  b = ReferenceFrame::CheckEndsAtECEF(RF1, NULL);
  BOOST_CHECK_MESSAGE( b, "Result 7.2:" << b);

  // Check RF2-> RF1 -> ECEF
  //std::cout << std::endl << std::endl;
  //std::cout << "The ReferenceFrame chain is:" << std::endl;
  //std::cout << "RF2 : " << RF2 << std::endl;
  //std::cout << "RF1 : " << RF1 << std::endl;
  //std::cout << "ECEF: " << ECEF << std::endl;
  //std::cout << "------" << std::endl;
  b = ReferenceFrame::CheckEndsAtECEF(RF2, NULL);
  BOOST_CHECK_MESSAGE( b, "Result 7.3:" << b);

  // Check RF3 (wrong ReferenceFrame)
  //std::cout << std::endl << std::endl;
  //std::cout << "The ReferenceFrame chain is:" << std::endl;
  //std::cout << "RF3 : " << RF3 << std::endl;
  //std::cout << "------" << std::endl;
  //b = ReferenceFrame::CheckEndsAtECEF(RF3, NULL);
  //BOOST_CHECK_MESSAGE( !b, "Result 51:" << b);
   

  gsl_matrix_free(Pos1);
  gsl_matrix_free(Pos2);
  gsl_matrix_free(Att1);
  gsl_matrix_free(Att2);
  gsl_matrix_free(Vel1);
  gsl_matrix_free(Vel2);
  //delete (RF3);
  delete (RF2);
  delete (RF1);
  delete (ECEF);

}


BOOST_AUTO_TEST_CASE (test_ReferenceFrame_PosInECEF)
{
  bool b;
  ReferenceFrame *ECEF;
  ECEF = new ReferenceFrame();
  { // ECEF -------------------------------------
    { // point at the origin
      double p[]       = {0, 0, 0};
      double presult[] = {0, 0, 0};
      gsl_vector_view pos = gsl_vector_view_array(p, REFFRAME_DIMS);
      
      gsl_vector *pECEF;
      pECEF = ECEF->PosInBF2PosInECEF(&pos.vector,0);
     
      for (unsigned int k=0; k<REFFRAME_DIMS;k++)
      {
        b = (presult[k] == pECEF->data[k]);
        BOOST_CHECK_MESSAGE ( b, "Result 8.1: k=" << k << " b=" << b );
      } 
      gsl_vector_free(pECEF);
    }
    { // point different from the origin
      double p[]       = {100, -2, 33};
      double presult[] = {100, -2, 33};
      gsl_vector_view pos = gsl_vector_view_array(p, REFFRAME_DIMS);
      gsl_vector *pECEF;
      pECEF = ECEF->PosInBF2PosInECEF(&pos.vector,0);
      for (unsigned int k=0; k<REFFRAME_DIMS;k++)
      {
        b = (presult[k] == pECEF->data[k]);
        BOOST_CHECK_MESSAGE ( b, "Result 8.2: k=" << k << " b=" << b );
      } 
      gsl_vector_free(pECEF);
    }
      
    { // RF1 (solidary) -> ECEF -----------------
      double pRF1[] = { 5,
                        3,
                        2};
      double aRF1[] = { 0, 1, 0,
                       -1, 0, 0,
                        0, 0, 1};
      double vRF1[] = { 0,
                        0,
                        0};
      gsl_matrix_view vpRF1 = gsl_matrix_view_array(pRF1, REFFRAME_DIMS, 1);
      gsl_matrix_view vaRF1 = gsl_matrix_view_array(aRF1, REFFRAME_DIMS, 1*REFFRAME_DIMS);
      gsl_matrix_view vvRF1 = gsl_matrix_view_array(vRF1, REFFRAME_DIMS, 1);
      ReferenceFrame *RF1 = new ReferenceFrame(ECEF, true, &vpRF1.matrix, &vaRF1.matrix, &vvRF1.matrix, "RF1");
      { // point at the origin
        double p[]       = {0, 0, 0};
        double presult[] = {-3, 5, 2};
        gsl_vector_view pos = gsl_vector_view_array(p, REFFRAME_DIMS);
    
        gsl_vector *pECEF;
        pECEF = RF1->PosInBF2PosInECEF(&pos.vector,0);
       
        for (unsigned int k=0; k<REFFRAME_DIMS;k++)
        {
          b = (presult[k] == pECEF->data[k]);
          BOOST_CHECK_MESSAGE ( b, "Result 8.3: k=" << k << " b=" << b );
        } 
        gsl_vector_free(pECEF);
      }
      { // point different from  origin
        double p[]       = {6,  -2, -4};
        double presult[] = {-1, 11, -2};
        gsl_vector_view pos = gsl_vector_view_array(p, REFFRAME_DIMS);
    
        gsl_vector *pECEF;
        pECEF = RF1->PosInBF2PosInECEF(&pos.vector,0);
       
        for (unsigned int k=0; k<REFFRAME_DIMS;k++)
        {
          b = (presult[k] == pECEF->data[k]);
          BOOST_CHECK_MESSAGE ( b, "Result 8.4: k=" << k << " b=" << b );
        } 
        gsl_vector_free(pECEF);
      }
      { // RF2 (non solidary) -> RF1 (solidary) -> ECEF
        double pRF2[] = { 2, 1, 0,
                         -1, 0, 1,
                          3, 4, 5};
        double aRF2[] = { 0, 1, 0, 1, 0, 0, 1, 0, 0, 
                         -1, 0, 0, 0, 1, 0, 0, 0,-1,
                          0, 0, 1, 0, 0, 1, 0, 1, 0};
        double vRF2[] = {-1,-1,-1,
                          1, 1, 1,
                          1 ,1 ,1};
        gsl_matrix_view vpRF2 = gsl_matrix_view_array(pRF2, REFFRAME_DIMS, 3);
        gsl_matrix_view vaRF2 = gsl_matrix_view_array(aRF2, REFFRAME_DIMS, 3*REFFRAME_DIMS);
        gsl_matrix_view vvRF2 = gsl_matrix_view_array(vRF2, REFFRAME_DIMS, 3);
        ReferenceFrame *RF2 = new ReferenceFrame(RF1, false, &vpRF2.matrix, &vaRF2.matrix, &vvRF2.matrix,"RF2");
        { // point at the origin
          double p2[]       = { 0,
                                0,
                                0};
          double presult2[] = {-5,-3,-8,
                                6, 6, 5,
                                5, 6, 1};

          
          gsl_vector_view pos2 = gsl_vector_view_array(p2, REFFRAME_DIMS);
     
          for (unsigned int l=0;l<3;l++) 
          {
            gsl_vector *pECEF;
            pECEF = RF2->PosInBF2PosInECEF(&pos2.vector,l);
            for (unsigned int k=0; k<REFFRAME_DIMS;k++)
            {
              b = (presult2[k*REFFRAME_DIMS+l] == pECEF->data[k]);
              BOOST_CHECK_MESSAGE ( b, "Result 8.5: l=" << l << " k=" << k << " b=" << b );
            }
          gsl_vector_free(pECEF);
          }
          { // point different from  origin
            double p3[]      = {6,
                                2,
                                4};
            double presult3[] = {-11, -5,-12,
                                   4, 12, 11,
                                   9, 10, -1};
            gsl_vector_view pos3 = gsl_vector_view_array(p3, REFFRAME_DIMS);
            
            for (unsigned int l=0;l<3;l++)
            {
              //std::cout << "------ L=" << l << std::endl;
              gsl_vector *pECEF;
              //gsl_vector_view aux;
              //aux = gsl_matrix_column(&pos3.matrix,l);
              pECEF = RF2->PosInBF2PosInECEF(&pos3.vector,l);
           
              for (unsigned int k=0; k<REFFRAME_DIMS;k++)
              {
                //std::cout << "p[" << k << "]=" << pECEF->data[k] << std::endl;
                b = (presult3[k*REFFRAME_DIMS+l] == pECEF->data[k]);
                BOOST_CHECK_MESSAGE ( b, "Result 8.6: l=" << l << " k=" << k << " b=" << b );
              } 
              gsl_vector_free(pECEF);
            }
          }
        }
        delete (RF2);
      }
      delete(RF1);
    }
  }
  delete (ECEF);
}


BOOST_AUTO_TEST_CASE (test_ReferenceFrame_PosInECEF2PosInBF)
{
  bool b;
  ReferenceFrame *ECEF;
  ECEF = new ReferenceFrame();
  { // ECEF -------------------------------------
    { // point at the origin
      double p[]       = {0, 0, 0};
      double presult[] = {0, 0, 0};
      gsl_vector_view pos = gsl_vector_view_array(p, REFFRAME_DIMS);
      
      gsl_vector *pBF;
      pBF = ECEF->PosInECEF2PosInBF(&pos.vector,0);
     
      for (unsigned int k=0; k<REFFRAME_DIMS;k++)
      {
        b = (presult[k] == pBF->data[k]);
        BOOST_CHECK_MESSAGE ( b, "Result 9.1: k=" << k << " b=" << b );
      } 
      gsl_vector_free(pBF);
    }
    { // point different from the origin
      double p[]       = {100, -2, 33};
      double presult[] = {100, -2, 33};
      gsl_vector_view pos = gsl_vector_view_array(p, REFFRAME_DIMS);
      gsl_vector *pBF;
      pBF = ECEF->PosInECEF2PosInBF(&pos.vector,0);
      for (unsigned int k=0; k<REFFRAME_DIMS;k++)
      {
        b = (presult[k] == pBF->data[k]);
        BOOST_CHECK_MESSAGE ( b, "Result 9.2: k=" << k << " b=" << b );
      } 
      gsl_vector_free(pBF);
    }
      
    { // RF1 (solidary) -> ECEF -----------------
      double pRF1[] = { 5,
                        3,
                        2};
      double aRF1[] = { 0, 1, 0,
                       -1, 0, 0,
                        0, 0, 1};
      double vRF1[] = { 0,
                        0,
                        0};
      gsl_matrix_view vpRF1 = gsl_matrix_view_array(pRF1, REFFRAME_DIMS, 1);
      gsl_matrix_view vaRF1 = gsl_matrix_view_array(aRF1, REFFRAME_DIMS, 1*REFFRAME_DIMS);
      gsl_matrix_view vvRF1 = gsl_matrix_view_array(vRF1, REFFRAME_DIMS, 1);
      ReferenceFrame *RF1 = new ReferenceFrame(ECEF, true, &vpRF1.matrix, &vaRF1.matrix, &vvRF1.matrix, "RF1");
      { // point at the origin
        double p[]       = { 0,
                             0,
                             0 };
        double presult[] = {-3,
                             5,
                            -2 };
        gsl_vector_view pos = gsl_vector_view_array(p, REFFRAME_DIMS);
        gsl_vector *pBF;
        pBF = RF1->PosInECEF2PosInBF(&pos.vector,0);
       
        for (unsigned int k=0; k<REFFRAME_DIMS;k++)
        {
          b = (presult[k] == pBF->data[k]);
          BOOST_CHECK_MESSAGE ( b, "Result 9.3: k=" << k << " b=" << b );
        } 
        gsl_vector_free(pBF);
      }
      { // point different from  origin
        double p[]       = { 6,
                            -2,
                            -4 };
        double presult[] = { -5,
                             -1,
                             -6 };
        gsl_vector_view pos = gsl_vector_view_array(p, REFFRAME_DIMS);
    
        gsl_vector *pBF;
        pBF = RF1->PosInECEF2PosInBF(&pos.vector,0);
       
        for (unsigned int k=0; k<REFFRAME_DIMS;k++)
        {
          b = (presult[k] == pBF->data[k]);
          BOOST_CHECK_MESSAGE ( b, "Result 9.4: k=" << k << " b=" << b );
        } 
        gsl_vector_free(pBF);
      }
      { // RF2 (non solidary) -> RF1 (solidary) -> ECEF
        double pRF2[] = { 2, 1, 0,
                         -1, 0, 1,
                          3, 4, 5};
        double aRF2[] = { 0, 1, 0, 1, 0, 0, 1, 0, 0, 
                         -1, 0, 0, 0, 1, 0, 0, 0,-1,
                          0, 0, 1, 0, 0, 1, 0, 1, 0};
        double vRF2[] = {-1,-1,-1,
                          1, 1, 1,
                          1 ,1 ,1};
        gsl_matrix_view vpRF2 = gsl_matrix_view_array(pRF2, REFFRAME_DIMS, 3);
        gsl_matrix_view vaRF2 = gsl_matrix_view_array(aRF2, REFFRAME_DIMS, 3*REFFRAME_DIMS);
        gsl_matrix_view vvRF2 = gsl_matrix_view_array(vRF2, REFFRAME_DIMS, 3);
        ReferenceFrame *RF2 = new ReferenceFrame(RF1, false, &vpRF2.matrix, &vaRF2.matrix, &vvRF2.matrix,"RF2");
        { // point at the origin
          double p2[]       = { 0,
                                0,
                                0};
          double presult2[] = { 6,-4,-3,
                                5, 5, 7,
                               -5,-6, 4};

          
          gsl_vector_view pos2 = gsl_vector_view_array(p2, REFFRAME_DIMS);
     
          for (unsigned int l=0;l<3;l++) 
          {
            gsl_vector *pBF2;
            pBF2 = RF2->PosInECEF2PosInBF(&pos2.vector,l);
            for (unsigned int k=0; k<REFFRAME_DIMS;k++)
            {
              b = (presult2[k*REFFRAME_DIMS+l] == pBF2->data[k]);
              BOOST_CHECK_MESSAGE ( b, "Result 9.5: l=" << l << " k=" << k << " b=" << b );
            }
            gsl_vector_free(pBF2);
          }
          { // point different from  origin
            double p3[]      = {6,
                                2,
                                4};
            double presult3[] = { 0, -2, -1,
                                  3, -1,  3,
                                 -1, -2, -2 };
            gsl_vector_view pos3 = gsl_vector_view_array(p3, REFFRAME_DIMS);
            for (unsigned int l=0;l<3;l++)
            {
              gsl_vector *pBF3;
              pBF3 = RF2->PosInECEF2PosInBF(&pos3.vector,l);
           
              for (unsigned int k=0; k<REFFRAME_DIMS;k++)
              {
                //std::cout << "p[" << k << "]=" << pECEF->data[k] << std::endl;
                b = (presult3[k*REFFRAME_DIMS+l] == pBF3->data[k]);
                BOOST_CHECK_MESSAGE ( b, "Result 9.6: l=" << l << " k=" << k << " b=" << b );
              } 
              gsl_vector_free(pBF3);
            }
          }
        } 
        delete (RF2);
      }
      delete(RF1);
    }
  }
  delete (ECEF);
}










BOOST_AUTO_TEST_SUITE_END()






  // six ways to detect and report the same error
  /*
  BOOST_CHECK( add(2,2)==4);            // continues on error

  BOOST_REQUIRE( add(2,2)==4);          // throws on error

  if (add(2,2) !=4)
    BOOST_ERROR("Here is an error..."); // continues on error
  
  if (add(2,2) !=4)
    BOOST_FAIL("Here is an error...");  // throws on error

  BOOST_CHECK_MESSAGE( add(2,2) ==4,    // continues on error  
                       "add(...) result: " << add(2,2)); 

  BOOST_CHECK_EQUAL( add(2,2), 4);      // continues on error
    */


