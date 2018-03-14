
const int splineSteps = 32/2;
const int profileDetail = 16/2;

const float ribbonWidth = 2.0f;
const float ribbonHeight = 0.125f;
const float ribbonOffset = 1.5f;
const float arrowHeadWidth = 3.0f;
const float arrowWidth = 2.0f;
const float arrowHeight = 0.5f;
const float tubeSize = 0.75f;


#ifndef M_PI
   #define M_PI 3.14159265358979323846
#endif



struct Mesh
{
    std::vector<v3> vertices;
    std::vector<v3> colors;
    std::vector<int> triangles;
};

Mesh createChainMesh(chain *C);