
#include <fstream>
#include <vector>

#include "cpdb/cpdb.h"
#include "cartoon.h"

#include <omp.h>

#include <chrono>
#include <iomanip>

using namespace std;

void writeToObj(string fileName, vector<Mesh> &meshes) {

    FILE *fich = NULL;

    cerr << "Writting to " << fileName << endl;

    fich = fopen(fileName.c_str(), "w");
    if (fich == NULL)
    {
        printf("Error opening file  '%s'!\n", fileName.c_str());
        exit(1);
    }


    for (int i = 0; i < meshes.size(); i++) {
        vector<v3> verts = meshes[i].vertices;
        vector<v3> colors = meshes[i].colors;

        for (int j = 0; j < verts.size(); j++) {
            // myfile << "v "<<meshes[i].vertices[j].x<<" "<<meshes[i].vertices[j].y<<" "<<meshes[i].vertices[j].z<<endl;
            fprintf(fich, "v %.3f %.3f %.3f %.3f %.3f %.3f\n", verts[j].x, verts[j].y, verts[j].z, colors[j].x, colors[j].y, colors[j].z);
        }
    }
    int cpt = 1;
    for (int i = 0; i < meshes.size(); i++) {
        vector<int> tri = meshes[i].triangles;
        for (int j = 0; j < tri.size(); j += 3) {
            // myfile << "f "<<cpt+tri[j]<<" "<<cpt+tri[j+1]<<" "<<cpt+tri[j+2]<<endl;
            fprintf(fich, "f %d %d %d\n", cpt + tri[j], cpt + tri[j + 1], cpt + tri[j + 2] );
        }
        cpt += meshes[i].vertices.size();
    }

    fclose(fich);

    cerr << "Wrote " << meshes.size() << " meshes to " << fileName << endl;
}

//Assuming there are always 2 atoms per residue
vector<Mesh> computeCartoonMesh(int nbChain, int *nbResPerChain, float *CA_OPositions, char *ssTypePerRes) {
    vector<Mesh> meshes(nbChain);
    if (nbChain <= 0) {
        return meshes;
    }

    #pragma omp parallel for num_threads(nbChain)
    for (int chainId = 0; chainId < nbChain; chainId++) {
        Mesh m = createChainMesh(chainId, nbResPerChain, CA_OPositions, ssTypePerRes);
        meshes[omp_get_thread_num()] = m;
    }

    // return &meshes[0];
    return meshes;
}

int main(int argc, char const *argv[]) {

    if (argc < 2) {
        std::cout << "Usage : " << argv[0] << " file.pdb " << endl;
        exit(-1);
    }

    pdb *P;
    P = initPDB();

    parsePDB((char *)argv[1], P, (char *)"");

    auto start = std::chrono::high_resolution_clock::now();

#if 0
    // vector<Mesh> meshes;
    vector<Mesh> meshes(P->size);
    #pragma omp parallel for num_threads(P->size)
    for (int chainId = 0; chainId < P->size; chainId++) {
        chain C = P->chains[chainId];

        Mesh m = createChainMesh(C);
        // meshes.push_back(m);
        meshes[omp_get_thread_num()] = m;
    }


#else
    vector<int> nbResPerChain(P->size);
    vector<float> CAOPos;
    vector<char> allSS;
    for (int chainId = 0; chainId < P->size; chainId++) {
        chain *C = &P->chains[chainId];
        nbResPerChain[chainId] = C->size;
        for (int r = 0; r < C->size; r++) {
            residue *R = &C->residues[r];
            atom *CA = getAtom(*R, (char *)"CA");
            atom *O = getAtom(*R, (char *)"O");
            char ss = R->ss;
            CAOPos.push_back(CA->coor.x);
            CAOPos.push_back(CA->coor.y);
            CAOPos.push_back(CA->coor.z);
            CAOPos.push_back(O->coor.x);
            CAOPos.push_back(O->coor.y);
            CAOPos.push_back(O->coor.z);
            allSS.push_back(ss);
        }
    }


    vector<Mesh> meshes = computeCartoonMesh(P->size, &nbResPerChain[0], &CAOPos[0], &allSS[0]);


#endif
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( stop - start ).count();

    cerr << duration << endl;
    writeToObj("output.obj", meshes);

    return 0;
}