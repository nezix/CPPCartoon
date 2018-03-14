#include <vector>
#include <math.h>

#include "cpdb/cpdb.h"
#include "cartoon.h"
#include "PeptidePlane.hpp"
#include "Spline.hpp"

using namespace std;

inline v3 lerp(v3 v0, v3 v1, float t){
    return v0 * (1-t) + v1 * t;
}

v3 *ellipseProfile(int n, float w, float h) {
    v3 *result = new v3[n];
    for(int i=0;i<n;i++){
        float t = (float)i / (float)n;
        float a = t * 2.0f * M_PI + M_PI/4.0f;
        float x = cosf(a) * w / 2.0f;
        float y = sinf(a) * h / 2.0f;
        result[i].x = x;
        result[i].y = y;
        result[i].z = 0.0f;
    }
    return result;
}


v3 *rectangleProfile(int n, float w, float h) {
    v3 *result = new v3[n];
    float hw = w / 2.0f;
    float hh = h / 2.0f;
    v3 segments[][2] = 
    {
    
        {v3(hw, hh, 0.0f),
        v3(-hw, hh, 0.0f)},

        {v3(-hw, hh, 0.0f),
        v3(-hw, -hh, 0.0f)},

        {v3(-hw, -hh, 0.0f),
        v3(hw, -hh, 0.0f)},

        {v3(hw, -hh, 0.0f),
        v3(hw, hh, 0.0f)}
    };

    int m = n / 4;
    int cpt=0;
    for(int s=0;s<4;s++){
        for(int i=0;i<m;i++){
            float t = (float)i / (float)m;
            v3 p = lerp(segments[s][0], segments[s][1], t);
            result[cpt++] = p;
        }
    }
    return result;
}




v3 *roundedRectangleProfile(int n, float w, float h) {
    v3 *result = new v3[n];
    float r = h / 2.0f;
    float hw = w/2.0f - r;
    float hh = h / 2.0f;

    v3 segments[][2] = 
    {
    
        {v3(hw, hh, 0.0f),
        v3(-hw, hh, 0.0f)},

        {v3(-hw, 0.0f, 0.0f),
        v3(0.0f, 0.0f, 0.0f)},

        {v3(-hw, -hh, 0.0f),
        v3(hw, -hh, 0.0f)},

        {v3(hw, 0.0f, 0.0f),
        v3(0.0f, 0.0f, 0.0f)}
    };

    int m = n / 4;
    int cpt=0;

    for(int s=0;s<4;s++){
        for(int i=0;i<m;i++){
            float t = (float)i / (float)m;
            v3 p;
            if(s == 0 || s == 2){
                p = lerp(segments[s][0], segments[s][1], t);
            }
            else if(s==1){
                float a = M_PI/2.0f + M_PI*t;
                float x = cosf(a) * r;
                float y = sinf(a) * r;
                p = segments[s][0] + v3(x, y, 0.0f);
            }
            else if(s==3){
                float a = 3*M_PI/2.0f + M_PI*t;
                float x = cosf(a) * r;
                float y = sinf(a) * r;
                p = segments[s][0] + v3(x, y, 0.0f);
            }
            result[cpt++] = p;
        }
    }
    return result;
}

void scaleProfile(v3 *&p, float s, int lenP) {
    for(int i=0;i<lenP;i++){
        p[i] = p[i] * s;
    }
}

void translateProfile(v3 *&p, float dx, float dy, int lenP) {
    for(int i=0;i<lenP;i++){
        p[i] = p[i] + v3(dx, dy, 0.0f);
    }
}

void segmentProfiles(PeptidePlane *pp1, PeptidePlane *pp2, int n, v3 *&p1, v3 *&p2){
    char type0 = pp1->Residue1->ss;
    char type1, type2;
    Transition(pp1, &type1, &type2);

    float offset1 = ribbonOffset;
    float offset2 = ribbonOffset;

    if (pp1->Flipped) {
        offset1 = -offset1;
    }
    if (pp2->Flipped) {
        offset2 = -offset2;
    }

    switch (type1) {
        case HELIX:
            if (type0 == STRAND) {
                p1 = roundedRectangleProfile(n, 0.0f, 0.0f);
            }
            else {
                p1 = roundedRectangleProfile(n, ribbonWidth, ribbonHeight);
            }
            translateProfile(p1, 0.0f, offset1, n);
            break;
        case STRAND:
            if (type2 == STRAND) {
                p1 = rectangleProfile(n, arrowWidth, arrowHeight);
            } else {
                p1 = rectangleProfile(n, arrowHeadWidth, arrowHeight);
            }
            break;
        default:
            if (type0 == STRAND) {
                p1 = ellipseProfile(n, 0.0f, 0.0f);
            } else {
                p1 = ellipseProfile(n, tubeSize, tubeSize);
            }
            break;
    }
    switch(type2) {
        case HELIX:
            p2 = roundedRectangleProfile(n, ribbonWidth, ribbonHeight);
            translateProfile(p2, 0.0f, offset2, n);
            break;
        case STRAND:
            p2 = rectangleProfile(n, arrowWidth, arrowHeight);
            break;
        default:
            p2 = ellipseProfile(n, tubeSize, tubeSize);
            break;
    }
    if(type1 == STRAND && type2 != STRAND) {
       p2 = rectangleProfile(n, 0.0f, arrowHeight);
    }
}

void segmentColors(PeptidePlane *pp, v3 *c1, v3 *c2) {
    // const minTemp = 10
    // const maxTemp = 50
    // f1 := pp.Residue2.Atoms["CA"].TempFactor
    // f2 := pp.Residue3.Atoms["CA"].TempFactor
    // t1 := fauxgl.Clamp((f1-minTemp)/(maxTemp-minTemp), 0, 1)
    // t2 := fauxgl.Clamp((f2-minTemp)/(maxTemp-minTemp), 0, 1)
    // c1 = fauxgl.MakeColor(Viridis.Color(t1))
    // c2 = fauxgl.MakeColor(Viridis.Color(t2))
    // return
    char type1, type2;
    Transition(pp, &type1, &type2);
    switch (type1) {
        case HELIX:
            *c1 = v3(1.0f, 0.71f, 0.2f);
            break;
        case STRAND:
            *c1 = v3(0.96f, 0.45f, 0.21f);
            break;
        default:
            *c1 = v3(0.02f, 0.47f, 0.47f);
            break;
    }
    switch (type2) {
        case HELIX:
            *c2 = v3(1.0f, 0.71f, 0.2f);
            break;
        case STRAND:
            *c2 = v3(0.96f, 0.45f, 0.21f);
            break;
        default:
            *c2 = v3(0.02f, 0.47f, 0.47f);
            break;
    }
    if (type1 == STRAND) {
        *c2 = *c1;
    }
}

void triangulateQuad( vector<int> &triangles,   vector<v3> &vertices,  vector<v3> &colors,
                     v3 p1, v3 p2, v3 p3, v3 p4, v3 c1, v3 c2, v3 c3, v3 c4){
    
    vertices.push_back(p1);
    int idp1 = vertices.size()-1;
    vertices.push_back(p2);
    vertices.push_back(p3);
    vertices.push_back(p4);

    colors.push_back(c1);
    colors.push_back(c2);
    colors.push_back(c3);
    colors.push_back(c4);

    //Add 2 triangles
    triangles.push_back(idp1);
    triangles.push_back(idp1+1);
    triangles.push_back(idp1+2);

    triangles.push_back(idp1);
    triangles.push_back(idp1+2);
    triangles.push_back(idp1+3);
}

void createSegmentMesh(Mesh &mesh, int curi, int n, PeptidePlane *pp1, PeptidePlane *pp2, PeptidePlane *pp3, PeptidePlane *pp4) {

    char type0 = pp2->Residue1->ss;
    char type1, type2;
    Transition(pp2, &type1, &type2);
    v3 c1, c2;
    segmentColors(pp2, &c1, &c2);
    v3 *profile1, *profile2;
    int lenProf1 = profileDetail;
    int lenProf2 = profileDetail;
    segmentProfiles(pp2, pp3, profileDetail, profile1, profile2);

    auto easeFunc = &Linear;
    if (!(type1 == STRAND && type2 != STRAND)) {

        easeFunc = &InOutQuad;
    }
    if (type0 == STRAND && type1 != STRAND) {
        easeFunc = &OutCirc;
    }
    // if type1 != STRAND && type2 == STRAND {
    //  easeFunc = &ease.InOutSquare
    // }
    if (curi == 0) {
        profile1 = ellipseProfile(profileDetail, 0.0f, 0.0f);
        easeFunc = &OutCirc;
    } else if (curi == n-1) {
        profile2 = ellipseProfile(profileDetail, 0.0f, 0.0f);
        easeFunc = &InCirc;
    }

    v3 **splines1 = new v3*[lenProf1];
    v3 **splines2 = new v3*[lenProf2];
    for(int i=0;i<lenProf1;i++){
        v3 p1 = profile1[i];
        v3 p2 = profile2[i];
        splines1[i] = splineForPlanes(pp1, pp2, pp3, pp4, splineSteps, p1.x, p1.y);
        splines2[i] = splineForPlanes(pp1, pp2, pp3, pp4, splineSteps, p2.x, p2.y);
    }


    for(int i = 0; i < splineSteps; i++){
        float t0 = easeFunc(((float)i) / splineSteps);
        float t1 = easeFunc(((float)i+1) / splineSteps);
        if(i == 0 && type1 == STRAND && type2 != STRAND) {
            v3 p00 = splines1[0][i];
            v3 p10 = splines1[profileDetail/4][i];
            v3 p11 = splines1[2*profileDetail/4][i];
            v3 p01 = splines1[3*profileDetail/4][i];
            triangulateQuad(mesh.triangles, mesh.vertices, mesh.colors, p00, p01, p11, p10, c1, c1, c1, c1);
        }
        for(int j = 0; j < profileDetail ; j++){
        // for(int j = 0; j < profileDetail && j < lenProf1 && j < lenProf2; j++){
            v3 p100 = splines1[j][i];
            v3 p101 = splines1[j][i+1];
            v3 p110 = splines1[(j+1)%profileDetail][i];
            v3 p111 = splines1[(j+1)%profileDetail][i+1];
            v3 p200 = splines2[j][i];
            v3 p201 = splines2[j][i+1];
            v3 p210 = splines2[(j+1)%profileDetail][i];
            v3 p211 = splines2[(j+1)%profileDetail][i+1];
            v3 p00 = lerp(p100, p200, t0);
            v3 p01 = lerp(p101, p201, t1);
            v3 p10 = lerp(p110, p210, t0);
            v3 p11 = lerp(p111, p211, t1);
            v3 c00 = lerp(c1, c2, t0);
            v3 c01 = lerp(c1, c2, t1);
            v3 c10 = lerp(c1, c2, t0);
            v3 c11 = lerp(c1, c2, t1);
            triangulateQuad(mesh.triangles, mesh.vertices, mesh.colors, p10, p11, p01, p00, c10, c11, c01, c00);

        }
    }
    // mesh.triangles = triangles;
    // mesh.colors = colors;
    // mesh.vertices = vertices;

}


Mesh createChainMesh(chain *C) {

    Mesh mesh;
    PeptidePlane *planes = new PeptidePlane[C->size-2];
    v3 previous;
    int nbPlanes = 0;
    for (int i = 0; i < C->size-2; i++) {

        residue *r1 = &C->residues[i];
        residue *r2 = &C->residues[i+1];
        residue *r3 = &C->residues[i+2];
        PeptidePlane plane = NewPeptidePlane(r1, r2, r3);
        if(plane.Residue1 != NULL){
            // TODO: better handling missing required atoms
            planes[nbPlanes++] = plane;
        }
    }
    for (int i = 0; i < nbPlanes; i++) {
        PeptidePlane *p = &planes[i];
        if(i > 0 && p->Side.dotProduct(previous) < 0.0f ){
            Flip(p);
        }
        previous = p->Side;
    }
    
    int n = nbPlanes - 3;
    for(int i = 0; i < n; i++) {
        // TODO: handle ends better
        PeptidePlane pp1 = planes[i];
        PeptidePlane pp2 = planes[i+1];
        PeptidePlane pp3 = planes[i+2];
        PeptidePlane pp4 = planes[i+3];

        createSegmentMesh(mesh, i, n, &pp1, &pp2, &pp3, &pp4);
    }
    return mesh;
}