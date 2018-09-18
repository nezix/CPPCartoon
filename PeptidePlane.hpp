
struct simpleRes
{
    int id = -1;
    char type;
    v3 posCA;
    v3 posO;
};

struct PeptidePlane2 {
    simpleRes Residue1;
    simpleRes Residue2;
    simpleRes Residue3;
    v3 Position;
    v3 Normal;
    v3 Forward;
    v3 Side;
    bool Flipped;
};



void Transition(const PeptidePlane2 &pp, char &type1, char &type2) {

    char t1 = pp.Residue1.type;
    char t2 = pp.Residue2.type;
    char t3 = pp.Residue3.type;
    type1 = t2;
    type2 = t2;
    if (t2 > t1 && t2 == t3){
        type1 = t1;
    }
    if (t2 > t3 && t1 == t2){
        type2 = t3;
    }
}

void Flip(PeptidePlane2 &pp) {
    pp.Side = pp.Side * -1;
    pp.Normal = pp.Normal * -1;
    pp.Flipped = !pp.Flipped;
}


PeptidePlane2 NewPeptidePlane(int chainId, v3 positions[], int ids[], char types[], int idResInPos[], int id, int id1, int id2){
    PeptidePlane2 newPP;

    if(id < 0 || id1 < 0 || id2 < 0){
        std::cout << "Something went wrong"<<std::endl;
        return newPP;
    }
    //CA of residue id
    v3 CA1 = positions[(idResInPos[chainId] + id)*2];
    //O of residue id
    v3 O1 = positions[(idResInPos[chainId] + id)*2+1];
    //CA of residue id1
    v3 CA2  = positions[(idResInPos[chainId] + id1)*2];

    simpleRes r1;
    r1.id = ids[idResInPos[chainId] + id];
    r1.type = types[idResInPos[chainId] + id];
    r1.posCA = CA1;
    r1.posO = O1;

    simpleRes r2;
    r2.id = ids[idResInPos[chainId] + id1];
    r2.type = types[idResInPos[chainId] + id1];
    r2.posCA = CA2;
    r2.posO = positions[(idResInPos[chainId] + id1)*2+1];

    simpleRes r3;
    r3.id = ids[idResInPos[chainId] + id2];
    r3.type = types[idResInPos[chainId] + id2];
    r3.posCA = positions[(idResInPos[chainId] + id1)*2];
    r3.posO = positions[(idResInPos[chainId] + id2)*2+1];

    v3 a = (CA2 - CA1).normalized();
    v3 b = (O1 - CA1).normalized();
    v3 c = v3::crossProduct(a, b).normalized();
    v3 d = v3::crossProduct(c, a).normalized();
    v3 p = (CA1 + CA2)/ 2.0f;
    newPP.Residue1 = r1;
    newPP.Residue2 = r2;
    newPP.Residue3 = r3;
    newPP.Position = p;
    newPP.Normal = c;
    newPP.Forward = a;
    newPP.Side = d;
    newPP.Flipped = false;

    return newPP;
}

