
struct PeptidePlane {
    residue *Residue1;
    residue *Residue2;
    residue *Residue3;
    v3 Position;
    v3 Normal;
    v3 Forward;
    v3 Side;
    bool Flipped;
};

PeptidePlane NewPeptidePlane(residue *r1, residue *r2, residue *r3){
    PeptidePlane newPP;

    atom *CA1 = getAtom(r1, (char *)"CA");
    atom *O1 = getAtom(r1, (char *)"O");
    atom *CA2  = getAtom(r2, (char *)"CA");

    if(CA1 == NULL || O1 == NULL || CA2 == NULL)
        return newPP;

    v3 ca1 = CA1->coor;
    v3 o1 = O1->coor;
    v3 ca2 = CA2->coor;

    v3 a = (ca2 - ca1).normalized();
    v3 b = (o1 - ca1).normalized();
    v3 c = v3::crossProduct(a, b).normalized();
    v3 d = v3::crossProduct(c, a).normalized();
    v3 p = (ca1 + ca2)/ 2.0f;
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

void Transition(PeptidePlane *pp, char *type1, char *type2) {

    char t1 = pp->Residue1->ss;
    char t2 = pp->Residue2->ss;
    char t3 = pp->Residue3->ss;
    *type1 = t2;
    *type2 = t2;
    if (t2 > t1 && t2 == t3){
        *type1 = t1;
    }
    if (t2 > t3 && t1 == t2){
        *type2 = t3;
    }
}

void Flip(PeptidePlane *pp) {
    pp->Side = pp->Side * -1;
    pp->Normal = pp->Normal * -1;
    pp->Flipped = !pp->Flipped;
}
