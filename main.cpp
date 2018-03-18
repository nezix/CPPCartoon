
#include <fstream>
#include <vector>

#include "cpdb/cpdb.h"
#include "cartoon.h"

#include <chrono>
#include <iomanip> 

using namespace std;

void writeToObj(string fileName, vector<Mesh> &meshes){

    FILE *fich = NULL;

    cerr << "Writting to "<<fileName<<endl;

    fich = fopen(fileName.c_str(),"w");
    if (fich == NULL)
    {
        printf("Error opening file  '%s'!\n",fileName.c_str());
        exit(1);
    }


    for(int i=0;i<meshes.size();i++){
        vector<v3> verts = meshes[i].vertices;
        for(int j=0;j<verts.size();j++){
            // myfile << "v "<<meshes[i].vertices[j].x<<" "<<meshes[i].vertices[j].y<<" "<<meshes[i].vertices[j].z<<endl;
            fprintf(fich, "v %.3f %.3f %.3f\n",verts[j].x,verts[j].y,verts[j].z);
        }
    }
    int cpt=1;
    for(int i=0;i<meshes.size();i++){
        vector<int> tri = meshes[i].triangles;
        for(int j=0;j<tri.size();j+=3){
            // myfile << "f "<<cpt+tri[j]<<" "<<cpt+tri[j+1]<<" "<<cpt+tri[j+2]<<endl;
            fprintf(fich, "f %d %d %d\n",cpt+tri[j],cpt+tri[j+1],cpt+tri[j+2] );
        }
        cpt+=meshes[i].vertices.size();
    }

    fclose(fich);



    // ofstream myfile;
    // cerr << "Writting to "<<fileName<<endl;
    // myfile.open (fileName);

    // for(int i=0;i<meshes.size();i++){
    //     for(int j=0;j<meshes[i].vertices.size();j++){
    //         myfile << "v "<<meshes[i].vertices[j].x<<" "<<meshes[i].vertices[j].y<<" "<<meshes[i].vertices[j].z<<endl;
    //     }
    // }
    // int cpt=1;
    // for(int i=0;i<meshes.size();i++){
    //     for(int j=0;j<meshes[i].triangles.size();j+=3){
    //         myfile << "f "<<cpt+meshes[i].triangles[j]<<" "<<cpt+meshes[i].triangles[j+1]<<" "<<cpt+meshes[i].triangles[j+2]<<endl;
    //     }
    //     cpt+=meshes[i].vertices.size();
    // }

    // myfile.close();
    cerr << "Wrote "<<meshes.size()<<" meshes to "<<fileName<<endl;
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
    vector<Mesh> meshes;
    chain C;

    for (int chainId = 0; chainId < P->size; chainId++) {
        C = P->chains[chainId];
        // for(int r=0;r<C->size;r++){
        //     if(C->residues[r].ss == STRAND)
        //         // cerr << "'"<<(int)C->residues[r].ss<<"'"<<endl;
        //         cerr << "'"<<C->id<<" / "<<C->residues[r].id<<endl;
        // }
        Mesh m = createChainMesh(C);
        meshes.push_back(m);
    }
    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( stop - start ).count();

    cerr << duration<<endl;
    writeToObj("output.obj",meshes);


    return 0;
}