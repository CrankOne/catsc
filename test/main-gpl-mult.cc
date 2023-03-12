/* Testing fixture for TrickTrack track finding lib
 *
 * Randomly distributes points on set of layers, fills the TT graph, evolves it
 * and provides track candidates. */

#include <TVector3.h>
#include <TRandom.h>

#include "catsc/cats.hh"

#include <iostream>
#include <fstream>
#include <cstdio>

/// Base class generating random spatial point
class iLayerPointGenerator {
public:
    /// Samples point
    virtual TVector3 get( TRandom & ) const = 0;
    virtual ~iLayerPointGenerator(){}
};

/// Base class for placed/oriented generators
class iSpatialGenerator : public iLayerPointGenerator {
protected:
    TVector3 _placement;
    TVector3 _rotation;
protected:
    iSpatialGenerator( const TVector3 & position
                     , const TVector3 & rotation )
        : _placement(position)
        , _rotation(rotation)
        {}
    /// Returns random point in local space
    virtual TVector3 get_local(TRandom &) const = 0;
public:
    /// Returns random point in global space
    TVector3 get( TRandom & rg ) const override {
        TVector3 r = get_local(rg);
        r.RotateZ(_rotation(2));
        r.RotateY(_rotation(1));
        r.RotateX(_rotation(0));
        r += _placement;
        return r;
    }
};

class GausPlaneGenerator : public iSpatialGenerator {
private:
    const Float_t _p[2][2];  // {x, y}*{mean, sigma}
public:
    GausPlaneGenerator( const TVector3 & p, const TVector3 & r
                      , Float_t xMean, Float_t xSigma
                      , Float_t yMean, Float_t ySigma )
        : iSpatialGenerator(p, r)
        , _p{{xMean, xSigma}, {yMean, ySigma}}
        {}
protected:
    TVector3 get_local(TRandom & rg) const override {
        TVector3 r;
        r.SetX(rg.Gaus(_p[0][0], _p[0][1]));
        r.SetY(rg.Gaus(_p[1][0], _p[1][1]));
        r.SetZ(0);
        return r;
    }
};

typedef catsc::TrackFinder<std::array<float, 3>> CATSTrackFinder;

/** Filter callback type */
struct ByAngleFilter : public CATSTrackFinder::iTripletFilter {
    float threshold;
    ByAngleFilter(float a) : threshold(a) {}
    bool matches( const std::array<float, 3> & a
                , const std::array<float, 3> & b
                , const std::array<float, 3> & c ) const override {
        float p[3][3] = {
            { a[0], a[1], a[2] },
            { b[0], b[1], b[2] },
            { c[0], c[1], c[2] }
        };
        float scProduct = 0, mod1 = 0, mod2 = 0;
        for( int i = 0; i < 3; ++i ) {
            p[0][i] -= p[1][i];
            mod1 += p[0][i]*p[0][i];
            p[2][i] -= p[1][i];
            mod2 += p[2][i]*p[2][i];

            scProduct += p[0][i]*p[2][i];
        }
        if( 0. == mod1 || 0. == mod2 || std::isnan(scProduct)
         || std::isnan(mod1)|| std::isnan(mod2) ) {
            return 0;
        }
        const float angle = acos(scProduct/sqrt(mod1*mod2));
        return angle > threshold ? 1 : 0;
    }
};

struct Collector : public CATSTrackFinder::iTrackCandidateCollector {
    void collect(const cats_HitData_t *, size_t ) override {}
};

int
main(int argc, char * argv[]) {
    // Init generators
    struct {
        iLayerPointGenerator * l;
        size_t nHitsMin, nHitsMax;
    } layers[] = {
        { new GausPlaneGenerator( {0, 0, 0}, {0, 0, 0}
                              ,  0.23, 1.2
                              , -0.12, 1.43 )
        , 15, 50 },
        { new GausPlaneGenerator( {-0.17, 0.23, 5}, {M_PI/6, M_PI/12, 0}
                              , -0.12, 3.2
                              ,  1.2,  0.75 )
        , 15, 100 },
        { new GausPlaneGenerator( {0.7, 0.25, 12}, {M_PI/8, -M_PI/8, M_PI/4}
                              ,  0.07, 1.5
                              ,  0.2,  1.75 )
        , 10, 10 },
        { new GausPlaneGenerator( {-0.02, 0.035, 17}, {-M_PI/9, -M_PI/67, M_PI/3}
                              ,  2.23, 1.5
                              ,  0.2,  0.26 )
        , 15, 50 },
        { nullptr, 0 }
    };
    // Other parameters
    size_t nIts = 10;
    unsigned int nMissingLayers = 1
               , minLength = 4;
    // Create reentrant track finder instance, filter and collector
    ByAngleFilter filter(.97*M_PI);
    Collector collector;
    CATSTrackFinder cats( sizeof(layers)/sizeof(*layers) - 1
                        , 10  // cells
                        , 10  // hits
                        , 10  // refs
                        );
    // Generate points
    TRandom gr;
    for(size_t nIt = 0; nIt < nIts; ++nIt) {
        cats_LayerNo_t nLayer = 0;
        for( auto p = layers; p->l; ++p, ++nLayer ) {
            size_t nSamples = p->nHitsMin + gr.Integer(p->nHitsMax - p->nHitsMin);
            for( size_t nPoint = 0; nPoint < nSamples; ++nPoint ) {
                TVector3 point = p->l->get(gr);
                cats.add(nLayer, std::array<float, 3>{ (float) point(0)
                                                     , (float) point(1)
                                                     , (float) point(2)} );
            }
        }

        for( cats_LayerNo_t nLayer = 0; nLayer < cats.n_layers(); ++nLayer ) {
            std::cout << "  #" << (int) nLayer << ": "
                      << cats.n_points(nLayer) << std::endl;
        }

        cats.evaluate(filter, nMissingLayers);
        cats.collect(collector, minLength);
        cats.reset();
        std::cout << "event #" << nIt << " done." << std::endl;
    }

    for( auto p = layers; p->l; ++p ) {
        delete p->l;
    }

    return 0;
}

