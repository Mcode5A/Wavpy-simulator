#ifndef     GEOMETRY_HPP_
#define     GEOMETRY_HPP_

/**
    \file   Geometry.hpp 
    \brief  Header file of the Geometry class.
*/

/**
    \brief  The Geometry class implements the geometry computations needed in
            a GNSS-R scenario.

*/
class Geometry
{
    public:
        Geometry();
        ~Geometry();
        void Compute();
        void Run(class ObservableProcessor* Proc);
        void GetTx2RxDelay();
        void GetTx2Scatterer2RxDelay();
        void GetTx2RxDelayRate();      
        void GetTx2Scatterer2RxDelayRate();
        void GetScatterer2TxUnitVector();
        void GetRx2ScattererUnitVector();
        void GetTx2RxUnitVector();
        void GetTxAntennaPosition();
        void GetRxAntennaPosition();
        void GetTxAntennaVelocity();
        void GetRxAntennaVelocity();
        void GetTxAntennaAcceleration();
        void GetRxAntennaAcceleration();
        void GetTxDirectionInRxAntennaFrame();
        void GetScattererDirectionInRxAntennaFrame();
        void GetIsoDelayMapOverSurfaceGrid();
        void GetIsoDopplerMapOverSurfaceGrid();

};

#endif
