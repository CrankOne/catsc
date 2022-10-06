/* Testing fixture for TrickTrack track finding lib
 *
 * Randomly distributes points on set of layers, fills the TT graph, evolves it
 * and provides track candidates. */

#include <TVector3.h>
#include <TRandom.h>

#include "catsc/cats.h"

#include <iostream>
#include <fstream>
#include <cstdio>

/// Base class generating random spatial point
class iLayerPointGenerator {
public:
    /// Samples point
    virtual TVector3 get( TRandom & ) const = 0;
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

/** Filter callback type */
int
by_angle_filter( cats_HitID_t, const cats_Float_t * a
               , cats_HitID_t, const cats_Float_t * b
               , cats_HitID_t, const cats_Float_t * c
               , void * threshold_ ) {
    float * threshold = reinterpret_cast<float *>(threshold_);
    cats_Float_t p[3][3] = {
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
    return angle > *threshold ? 1 : 0;
}

struct TrackCandsWriteStruct {
    std::unordered_map<cats_HitID_t, std::tuple<float, float, float>> hits;
    FILE * outfile;
};

static void
_write_cands( cats_HitID_t * tc
            , size_t tcLen
            , void * tcws_ ) {
    TrackCandsWriteStruct * tcws = reinterpret_cast<TrackCandsWriteStruct*>(tcws_);
    for(size_t tcN = 0; tcN < tcLen; ++tcN) {
        fprintf( tcws->outfile, "%e\t%e\t%e\n"
               , std::get<0>(tcws->hits[tc[tcN]])
               , std::get<1>(tcws->hits[tc[tcN]])
               , std::get<2>(tcws->hits[tc[tcN]])
               );
    }
    fputs("\n\n", tcws->outfile);
}

int
main(int argc, char * argv[]) {
    // Init generators
    struct {
        iLayerPointGenerator * l;
        size_t nSamples;
    } layers[] = {
        { new GausPlaneGenerator( {0, 0, 0}, {0, 0, 0}
                              ,  0.23, 1.2
                              , -0.12, 1.43 )
        , 10 },
        { new GausPlaneGenerator( {-0.17, 0.23, 5}, {M_PI/6, M_PI/12, 0}
                              , -0.12, 3.2
                              ,  1.2,  0.75 )
        , 5 },
        { new GausPlaneGenerator( {0.7, 0.25, 12}, {M_PI/8, -M_PI/8, M_PI/4}
                              ,  0.07, 1.5
                              ,  0.2,  1.75 )
        , 18 },
        #if 1
        { new GausPlaneGenerator( {-0.02, 0.035, 17}, {-M_PI/9, -M_PI/67, M_PI/3}
                              ,  2.23, 1.5
                              ,  0.2,  0.26 )
        , 30 },
        #endif
        { nullptr, 0 }
    };
    // Generate points
    TrackCandsWriteStruct tcws;
    TRandom gr;
    cats_Layers * acc = cats_layers_create(4);
    {
        cats_LayerNo_t nLayer = 0;
        std::ofstream pointsOfstream("/tmp/points.dat");
        for( auto p = layers; p->l; ++p, ++nLayer ) {
            for( size_t nPoint = 0; nPoint < p->nSamples; ++nPoint ) {
                TVector3 point = p->l->get(gr);
                cats_HitID_t hitID = nLayer*1000 + nPoint;
                pointsOfstream << point(0) << "\t" << point(1) << "\t" << point(2)
                              << std::endl;
                tcws.hits.emplace( hitID
                                 , std::tuple<float, float, float>{point(0), point(1), point(2)});
                cats_layer_add_point( acc, nLayer
                                    , point(0), point(1), point(2)
                                    , hitID );
            }
        }
    }

    cats_CellsPool * cats = cats_cells_pool_create(4);
    float angleThreshold = .95*M_PI;
    cats_evolve(acc, cats, by_angle_filter, &angleThreshold, 1);

    FILE * outf = fopen("/tmp/cells-dump.json", "w");
    cats_dump_json(acc, outf);
    fclose(outf);

    tcws.outfile = fopen("/tmp/tracklets.dat", "w");
    cats_for_each_track_candidate(acc, 3, _write_cands, &tcws);
    fclose(tcws.outfile);

    cats_cells_pool_delete(cats);
    cats_layers_delete(acc);

    return 0;
}

