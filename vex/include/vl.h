vector vl_old(int samples; float _maxdist, scatter, absorb)
{

    float maxdist = _maxdist;
    float ext = scatter + absorb;
    vector nI = normalize(I);

    maxdist = min(maxdist, length(I));

    vector tmp = 0;

    // loop over samples
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


float misweights(float pdf, otherpdf)
{
    return pdf/(pdf + otherpdf);
    
    //return pdf*pdf/(pdf*pdf + otherpdf*otherpdf);
}


vector vl(int samples; float _maxdist, scatter, absorb; vector perlight[])
{

    float maxdist = _maxdist;
    float ext = scatter + absorb;
    vector nI = normalize(I);
    maxdist = min(maxdist, length(I));
    float pdf_u = 1/maxdist;
    vector tmp = 0; //accumulate all distributions here

    int lights[] = getlights();

    // loop over all light
    for (int i=0; i<len(lights); i++)
    {
        vector tmpl = 0;// accumulate distribution of current light
        // getlight pos
        vector posL = ptransform( getlightname( lights[i] ) , "space:camera" , {0,0,0});
        
        // get delta and D
        float delta = dot( posL, nI);
        float D = length( posL-nI*delta);
        
        // get angles
        float thetaA = atan( 0 - delta, D );
        float thetaB = atan( maxdist-delta, D );
        
        // loop over samples
        for (int sample=0; sample<samples; sample++)
        {   
        
            //EQUIANGULAR
            // compute stratified sample and pdf for it
            float r = (rand(i+sample+SID) + sample) /samples;
            float tt = D*tan( (1-r)*thetaA + r*thetaB );
            float pdf_e = D/((thetaB-thetaA)*(D*D + tt*tt));
            float dist_e = (delta + tt);
            vector pos_e = nI * dist_e;
            vector tmp_e = 0;

            // illuminance loop only for current light
            illuminance(pos_e, {0,0,0}, "lightmask", getlightname( lights[i] ) )
            {
                shadow(Cl, pos_e, L);
                tmp_e += Cl*scatter*exp((-dist_e-length(L)) * ext)*0.5/pdf_e;
            }
            tmpl +=tmp_e * misweights(pdf_e, pdf_u);
            
            // UNIFORM
            vector tmp_u = 0;
            float dist_u = maxdist * r;
            vector pos_u = nI * dist_u;
            illuminance(pos_u, {0,0,0}, "lightmask", getlightname( lights[i] ) )
            {
                shadow(Cl, pos_u, L);
                tmp_u += Cl*scatter*exp((-dist_u-length(L)) * ext)*0.5/pdf_u;
            }
            tt = delta - dist_u;
            pdf_e = D/((thetaB-thetaA)*(D*D + tt*tt));
            tmpl +=tmp_u * misweights(pdf_u, pdf_e);
            
            
        }
        
        tmpl /= samples;
        perlight[i] = tmpl;// export perlight values in array
        tmp += tmpl;
    }
    
    return tmp;
}
