#ifndef VESSEL_SURFACE_RAY_H
#define VESSEL_SURFACE_RAY_H

#include "glucoseGrid.h"
#include <algorithm>

namespace vessel {
// Half-open projected-triangle coverage: adjacent triangles share the exact
// same edge predicate with opposite signs. An edge belongs to just one of
// them, so no distance-based merging of separate surface crossings is needed.
inline bool surfaceRayIntersection(Vec a, Vec b, Vec c, int axis,
                                   long double u, long double v, double& hit) {
    const int U=(axis+1)%3,V=(axis+2)%3;
    const auto edge=[&](const Vec& p,const Vec& q,long double x,long double y) {
        return (static_cast<long double>(p[U])-x)*(static_cast<long double>(q[V])-y)
             - (static_cast<long double>(p[V])-y)*(static_cast<long double>(q[U])-x);
    };
    const long double area=edge(a,b,c[U],c[V]);
    if(area==0)return false;
    if(area<0)std::swap(b,c);
    const auto accepts=[&](long double e,const Vec& p,const Vec& q) {
        const long double du=static_cast<long double>(q[U])-p[U];
        const long double dv=static_cast<long double>(q[V])-p[V];
        return e>0 || (e==0 && (dv>0 || (dv==0 && du<0)));
    };
    const long double wa=edge(b,c,u,v),wb=edge(c,a,u,v),wc=edge(a,b,u,v);
    if(!accepts(wa,b,c)||!accepts(wb,c,a)||!accepts(wc,a,b))return false;
    // Sum the point-relative edge values for consistent interpolation even
    // when the projection is very slender.
    hit=static_cast<double>((wa*a[axis]+wb*b[axis]+wc*c[axis])/(wa+wb+wc));
    return true;
}
}
#endif
