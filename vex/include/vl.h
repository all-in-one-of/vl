vector vl_old(int samples; float _maxdist, scatter, absorb)
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
        float pdf = 1/maxdist;
        illuminance(pos, {0,0,0})
        {
            shadow(Cl, pos, L);
            tmp += Cl*scatter*exp((-dist-length(L)) * ext)*0.5/pdf;
        }
    }


    return tmp/samples;

}





vector vl(int samples; float _maxdist, scatter, absorb)
{

    float maxdist = _maxdist;
    float ext = scatter + absorb;
    vector nI = normalize(I);
    maxdist = min(maxdist, length(I));
    vector tmp = 0;

    int lights[] = getlights();
    
    for (int i=0; i<len(lights); i++)
    {
    
        // getlight pos
        vector posL = ptransform( getlightname( lights[i] ) , "space:camera" , {0,0,0});
        
        //get delta and D
        float delta = dot( posL, nI);
        float D = length( posL-nI*delta);
        
        float thetaA = atan( 0 - delta, D );
        float thetaB = atan( maxdist-delta, D );
        
        for (int sample=0; sample<samples; sample++)
        {   
            float r = (rand(i+SID) + sample) /samples;

            float tt = D*tan( (1-r)*thetaA + r*thetaB );
            
            //float dist = maxdist * (i+r)/samples;
            float pdf = D/((thetaB-thetaA)*(D*D + tt*tt));
            
            
            float dist = (delta + tt);
            vector pos = nI * dist;

            illuminance(pos, {0,0,0}, "lightmask", getlightname( lights[i] ) )
            {
                shadow(Cl, pos, L);
                tmp += Cl*scatter*exp((-dist-length(L)) * ext)*0.5/pdf;
            }
            //tmp = pdf;
        }
    
    }

    return tmp/samples;
}

