vector vl(int samples; float _maxdist, scatter, absorb)
{

//    float maxdist = 15;
//    float scatter = 0.1;
//    float absorb = 0;
//    int samples = 16;

    float maxdist = _maxdist;
    float ext = scatter + absorb;
    vector nI = normalize(I);

    maxdist = min(maxdist, length(I));

    vector tmp = 0;

    for (int i=0; i<samples; i++)
    {   
        float r = rand(i+SID);
        float dist = maxdist * (i+r)/samples;
        vector pos = nI * dist;
        illuminance(pos, {0,0,0})
        {
            shadow(Cl, pos, L);
            tmp += Cl*scatter*exp((-dist-length(L)) * ext)*0.5;
        }
    }


    return tmp*maxdist/samples;

}


