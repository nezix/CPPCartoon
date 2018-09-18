
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





//Assumes having 2 positions per residue, CA and O
vector<Mesh> createMeshes(v3 *positions, int *ids, char *types, int nbChain, int *residuesPerChain) {

    int *idResInPos = new int[nbChain];

    int id = 0;
    for (int i = 0; i < nbChain; i++) {
        idResInPos[i] = id;
        id += residuesPerChain[i];
    }
    vector<Mesh> meshes;
    meshes.resize(nbChain);

    // #pragma omp parallel for num_threads(nbChain)
    for (int chainId = 0; chainId < nbChain; chainId++) {
        Mesh m = createChainMesh2(chainId, positions, ids, types, nbChain, residuesPerChain, idResInPos);
        meshes[chainId] = m;
    }
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

    vector<v3> v;
    vector<int> is;
    vector<char> tps;
    vector<int> resPC;

    for (int chainId = 0; chainId < P->size; chainId++) {
        chain C = P->chains[chainId];
        resPC.push_back(C.size);
        for (int i = 0; i < C.size; i++) {
            const residue &R = C.residues[i];
            atom *CA = getAtom(R, (char *)"CA");
            atom *O = getAtom(R, (char *)"O");
            v.push_back(CA->coor);
            v.push_back(O->coor);
            is.push_back(R.id);
            tps.push_back(R.ss);
        }
    }

    v3 *positions = &v[0];;
    int *ids = &is[0];
    char *types = &tps[0];
    int nbChain = P->size;
    int *residuesPerChain = &resPC[0];


    auto start = std::chrono::high_resolution_clock::now();

    vector<Mesh> meshes = createMeshes(positions, ids, types, nbChain, residuesPerChain);

    // // vector<Mesh> meshes;
    // vector<Mesh> meshes(P->size);
    // #pragma omp parallel for num_threads(P->size)
    // for (int chainId = 0; chainId < P->size; chainId++) {
    //     chain C = P->chains[chainId];

    //     Mesh m = createChainMesh(C);
    //     // meshes.push_back(m);
    //     meshes[omp_get_thread_num()] = m;
    // }

    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( stop - start ).count();

    cerr << duration << endl;
    writeToObj("output.obj", meshes);

    return 0;
}

