// Gillespie Algorithm, linear aggregation(primary nucleation, add-on and fragmentation ) with diffusion in 2d grid. A single compartment could hold many particles and species. Diffusion here means when one of the particle across from one compartment to another.
// Ace Shen, 01/05/2018

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <random>
#include <cstdlib>
#include <pthread.h>
using namespace std;

double totalpropensity( vector< vector< vector<double> > > a, vector< vector<double> > m, int lmax, void* args);
double lpropensity( vector< vector<double> > a0, vector< vector<double> > a, vector< vector<double> > m, int i, void* args);
double diffconst(int i, void* args);
int write2file(string fname, vector<double> data);
void* single_traj(void* args);
typedef struct outputData {
    int maxl;
}outputData;

// Structure for input parameters
typedef struct Parameters {
    int isample, nstep, nmonomer, primnuc, xcompartment, ycompartment, randini, Ndatapoint, p0, n2;
    double kn, kp, kf, km, k2, Dscale;
    float xlength, ylength, Nt;
    string datafilepath;
    
    Parameters( int isampleIn, int nstepIn, int nmonomerIn, int primnucIn, int xcompartmentIn, int ycompartmentIn, int randiniIn, float NtIn, int NdatapointIn, double knIn, double kpIn, double kfIn, double kmIn, double k2In, int n2In, float xlengthIn, float ylengthIn, double DscaleIn, int p0In, string datafilepathIn) {
        
        this->isample = isampleIn;
        this->nstep = nstepIn;
        this->nmonomer = nmonomerIn;
        this->primnuc = primnucIn;
        this->xcompartment = xcompartmentIn;
        this->ycompartment = ycompartmentIn;
        this->randini = randiniIn;
        this->Nt = NtIn;
        this->Ndatapoint = NdatapointIn;
        this->kn = knIn;
        this->kp = kpIn;
        this->kf = kfIn;
        this->km = kmIn;
        this->k2 = k2In;
        this->n2 = n2In;
        this->xlength = xlengthIn;
        this->ylength = ylengthIn;
        this->Dscale = DscaleIn;
        this->p0 = p0In;
        this->datafilepath = datafilepathIn;
    }
}Parameters;

int main(){
    cout.precision(10);
    //================== Read initial parameters from parameters.ini =======================
    ifstream paraini("parameters.ini");
    string name;
    double para_ini[19+1];
    int ipara = 1;
    while( paraini >> name >> para_ini[ipara] ){
        ipara++;
    }
    int nsample = para_ini[1];
    int nstep = para_ini[2];
    int nmonomer = para_ini[3];
//    float xlength = para_ini[4]; // in micron
//    float ylength = para_ini[5];
//    int xcompartment = para_ini[6];
//    int ycompartment = para_ini[7];
    float xlength = 1; // in micron
    float ylength = 1;
    int xcompartment = 1;
    int ycompartment = 1;
    int primnuc = para_ini[4];
    double kn = para_ini[5]; // [ 1/s * [dm^2/mol]^(nc-1) ]
    double kp = para_ini[6]; // [ 1/s * dm^2/mol ] (This is ismilar to 1/s * l / mol, l = dm^3)
    double kf = para_ini[7]; // [ 1/s ]
    double km = para_ini[8];
    double k2 = para_ini[9];
    int n2 = para_ini[10];
//    int randini = para_ini[15]; // randini = 0 => monomers are at the center; randini = 1 => randomly distributed;
    int randini = 0;
    float Nt = para_ini[11];
    int Ndatapoint = para_ini[12];
//    double Dscale = para_ini[18];
    double Dscale = 0;
    int p0 = para_ini[13];
    
    string datafilepath;
    datafilepath = name;
//    // create as many subdirectories under data/sample/ as nsample
//    for (int i=1; i<=nsample; i++)
//    {
//        string dir_isample = "mkdir -p " + datafilepath;
//        dir_isample += to_string(i);
//        system(dir_isample.c_str());
//    }
    
    //
    //kn = kn / pow(6.e13* xlength * ylength, primnuc-1); // convert unit to [1/s]
    //kp = kp / (6.e13 * xlength * ylength); // convert unit to [1/s]
    // km = .........
    
    // pack the input parameters as params for threads
    vector<Parameters*> params(nsample, nullptr);
    // Start pthread thing
    pthread_t p_threads[nsample]; // nsample threads...
    // ===== sample -> thread
    for(int isample=1; isample < nsample + 1; isample++){
        params[isample] = new Parameters(isample, nstep, nmonomer, primnuc, xcompartment, ycompartment, randini, Nt, Ndatapoint, kn, kp, kf, km, k2, n2, xlength, ylength, Dscale, p0, datafilepath);
        pthread_create(&p_threads[isample], NULL, single_traj, (void*)params[isample]);
    }// end of sample;
    
    ofstream maxlfile;
    maxlfile.open(datafilepath+"maxl.txt");
    outputData* m = NULL;
    for(int isample=1; isample < nsample +1; isample++){
        pthread_join(p_threads[isample], (void**)&m);
        maxlfile << m->maxl << "\t";
    }
    maxlfile.close();
    
    vector<double> time(Ndatapoint+1,0);
    for(int i=0;i<Ndatapoint+1;i++){
        time.at(i) = i * Nt/Ndatapoint;
    }
    ofstream tfile;
    tfile.open(datafilepath + "time.dat");
    for(int i=0;i<Ndatapoint+1;i++){
        tfile << time.at(i) << " ";
    }
    tfile.close();
    
    ofstream parafile;
    parafile.open(datafilepath+"para.txt");
    parafile << "nsample = " << nsample << endl;
    parafile << "nstep = " << nstep << endl;
    parafile << "nmonomer = " << nmonomer << endl;
//    parafile << "x(y)l = " << xlength << endl;
//    parafile << "x(y)compartment = " << xcompartment << endl;
    parafile << "nc = " << primnuc << endl;
//    parafile << "kn = " << kn * pow(6.e13* xlength * ylength, primnuc-1) << endl;
    parafile << "kn [1/s] = " << kn << endl;
//    parafile << "Input kp [ dm^2/(mol*s) ] = " << kp * (6.e13 * xlength * ylength) << endl;
    parafile << "kp [1/s] = " << kp  << endl;
    parafile << "kf = " << kf << endl;
    parafile << "km = " << km << endl;
    parafile << "k2 = " << k2 << endl;
    parafile << "n2 = " << n2 << endl;
//    parafile << "randini = " << randini << endl;
    parafile << "Nt = " << Nt << endl;
    parafile << "Ndatapoint = " << Ndatapoint << endl;
//    parafile << "Dscale = " << Dscale << endl;
    parafile << "p0 = " << p0 << endl;
    parafile.close();
    
    return(0);
}


double totalpropensity( vector< vector< vector<double> > > a, vector< vector<double> > m, int lmax, void* args){
    Parameters* params = (Parameters*) args;
    double t = 0;
    for(int j=0; j<=lmax; j++){
         t = t + lpropensity( a.at(0), a.at(j), m, j, params );
    }
    return t;
}

double lpropensity( vector< vector<double> > a0, vector< vector<double> > a, vector< vector<double> > m, int i, void* args){
    Parameters* params = (Parameters*) args;
    int xc = params->xcompartment;
    int yc = params->ycompartment;
    int nc = params->primnuc;
    int n2 = params->n2;
    int nmonomer = params->nmonomer;
    
    double dd;
    dd = diffconst(i, params);
    
    // ip = ( ipx+, ipx-, ipy+, ipy-, (ipprimarynuc or ipaddition), ipfragmentation, ipdissociation )
    vector<double> ip(7,0);
    //cout << ip.at(0) << endl;
    
    double a0sum = 0, asum = 0, a0asum = 0, a0powprimsum = 0, a0powsecMsum = 0; //  sum(a0^primnuc), sum(  ), M = total monomers in agg.
    for( int x=0; x< xc; x++ ){
        for(int y=0; y< yc; y++){
            asum = asum + a.at(x).at(y);
            a0asum = a0asum + a.at(x).at(y) * a0.at(x).at(y); // a0 * a
            if(a0.at(x).at(y) >= nc){
                a0powprimsum = a0powprimsum + pow(a0.at(x).at(y), nc);}
                //cout << pow(a0.at(x).at(y), nc)<< endl;
            if(a0.at(x).at(y) >= n2){
                a0powsecMsum = a0powsecMsum + pow(a0.at(x).at(y), n2)*m.at(x).at(y);}
        }
    }
    
    ip.at(0) = dd * ( asum - accumulate( a.at(xc-1).begin(),a.at(xc-1).end() ,0) );
    ip.at(1) = dd * ( asum - accumulate( a.at(0).begin(),a.at(0).end() ,0) );
    
    double asumyend = 0, asumybegin = 0;
    for( int x=0; x< xc; x++ ){
        asumyend = asumyend + a.at(x).at(yc-1);
        asumybegin = asumybegin + a.at(x).at(0);
    }
    ip.at(2) = dd * ( asum - asumyend );
    ip.at(3) = dd * ( asum - asumybegin );
    
    // ipprim (ipsec), ipplusone, ipfrag, ipminusone
    // Following equation 3 in doi: 10.1119/1.4870004
    if(i == 0){
        ip.at(4) = params->kn * a0powprimsum + params->k2 * a0powsecMsum;
        ip.at(5) = 0;
        ip.at(6) = 0;
    }
    else if(i == 1){
        ip.at(4) = 2. * params->kp * a0asum;
        ip.at(5) = (i + nc -1 -1) * params->kf * asum;;
        ip.at(6) = 0;
    }
    else{
        ip.at(4) = 2. * params->kp * a0asum;
        ip.at(5) = (i + nc -1 -1) * params->kf * asum; // The length of an aggregate is l = (i+nc-1). There are (l-1) sites for the fragmentation to occur
        ip.at(6) = 2. * params->km * asum;
     }
    //cout << i << "\t"<< ip.at(0) << "\t" << ip.at(4) << "\t" << ip.at(6) << "\t" << endl;
    return accumulate( ip.begin(), ip.end(), 0.);
}

double diffconst(int i, void* args){
    Parameters* params = (Parameters*) args;
    
    // For calculating diffusion coefficients
    float hx = params->xlength/params->xcompartment;
    float hy = params->ylength/params->ycompartment;
    double Dscale = params->Dscale;
    if(hx!=hy){
        cout << "Warning: hx and hy are not equal! Stop." << endl;
        return(0);
    }
    // The diffusion coefficient for monomer in here is about 27 micron^2/s
    // Reference: PRL 112, 098101 (2014), diffusion coefficient for a rigid rod
    // D = k_B*T/(6*pi*eta*R), R = length/(2*log(2*length/diameter))
    // Assume T = 310 K,  eta = 0.5 mPa s, so k_B * T / (6*pi*eta) = 0.454 micron^3 s^-1
    double l = 0; // the length of an aggregate
    double d = 0, R = 0; // The rate constant for diffustion in master equation. d = D/h^2, where D is the original diffusion coefficient.
    // Assume diameter = 10 nm = 0.01 micron
    // double diameter = 0.01;
    // Assume diameter = 1 nm = 0.001 micron
    double diameter = 0.001;
    // the length of an aggregate is (number * diameter)
    if(i==0){
        l = (i + 1)*diameter;
    }
    else{
        l = (i + params->primnuc - 1)*diameter;
    }
    R = l / (2*log10(2*l/diameter));
    d = 0.454 / (R * hx * hx) ;
    d = d * Dscale; // eta->100*eta; to slow down diffusion
    //cout << d * params->hx * params->hx << endl;
    return d;
}

int write2file(string fname, vector< vector<double> > data){
    ofstream file;
    file.open(fname, std::ofstream::app);
    for(int i=0; i<data.size(); i++){
        for(int j=0; j<data.at(i).size(); j++){
            file << data.at(i).at(j) << " ";
        }
        file << endl;
    }
    file.close();
    return(0);
}


void* single_traj(void* args) {
    
    Parameters* params = (Parameters*) args;
    int isample = params->isample;
    int nstep = params->nstep;
    int nmonomer = params->nmonomer;
    float xlength = params->xlength; // in micron
    float ylength = params->ylength;
    int xcompartment = params->xcompartment;
    int ycompartment = params->ycompartment;
    int primnuc = params->primnuc;
    double kn = params->kn; // [1/s * 1/M^(nc-1)]
    double kp = params->kp; // [1/M * 1/s] or [1/s * liter/mol]
    double kf = params->kf; // [1/s]
    double km = params->km;
    double k2 = params->k2;
    int n2 = params->n2;
    int randini = params->randini;
    float Nt = params->Nt;
    double Ndatapoint = params->Ndatapoint;
    int p0 = params->p0;
    string datafilepath = params->datafilepath;
    
    random_device rd;     // only used once to initialise (seed) engine
    mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    
    uniform_real_distribution<double> random (0, 1.);
    uniform_real_distribution<double> random2 (0, 1.);
    uniform_int_distribution<int> xran (0,xcompartment-1);
    uniform_int_distribution<int> yran (0,ycompartment-1);
    
    int nfiles = 0;
    //=============== initialize vectors for monomer and nucleation seed ===============
    // adataxy.at(x).at(y)
    vector<double> adatay(ycompartment, 0); // initialze ycompartment data with initial value 0;
    vector< vector<double> > adataxy(xcompartment, adatay); // initialize xcompartment data, adatay for each x;
    // adata2d.at(length).at(x).at(y)
    vector< vector< vector<double> > > adata2d(2, adataxy); // initialize two states, # of monomer and seeds, for each xycompartmnet;
    switch(randini){
    // == Initialize nmonomers at the center of the box ==
        case 0: {
                adata2d.at(0).at((xcompartment+1)/2 -1).at((ycompartment+1)/2 -1) = nmonomer;
                break;
        }
    // == Initialize nmonomers randomly distributed ==
        case 1: {
                for(int imono = 1; imono <= nmonomer; imono++){
                    int xloc=0, yloc=0;
                    xloc = xran(rng);
                    yloc = yran(rng);
                    adata2d.at(0).at(xloc).at(yloc) = adata2d.at(0).at(xloc).at(yloc) + 1;
                }
                break;
        }
    }
    // initialize nucleation seeds adata2d.at(1).at....
    adata2d.at(1).at((xcompartment+1)/2 -1).at((ycompartment+1)/2 -1) = p0;
    //    adata2d.at(1).at(5).at(5) = 10;
    //    adata2d.at(1).at(4).at(5) = 10;
    //    adata2d.at(1).at(5).at(4) = 10;
    //adata2d.at(1).at(4).at(4) = 10;
    //================================================================================
    int imax = 1;
    double t = 0; // Time
    double tpropensity = 0.;
    int count = 0;
    vector<double> time(Ndatapoint+1,0);
    for(int i=0;i<Ndatapoint+1;i++){
        time.at(i) = i * Nt/Ndatapoint;
    }
    
    // ===== data holder =======
    vector<double> iagg(Ndatapoint+1, 0);
    vector< vector<double> > outd(imax+1, iagg);
    
    
    // ====== write initial adata2d into data holder ==========
    for(int i = 0; i < imax+1; i++){
        for(int x=0; x < adata2d.at(i).size(); x++){
            for(int y=0; y < adata2d.at(i).at(x).size(); y++){
                outd.at(i).at(count) = adata2d.at(i).at(x).at(y);
            }
        }
    }
    count = count + 1;
    // ================================================
    // The main for loop. Each cycle there will be one reaction, diffusion, aggregation or fragmentation, happening in one of the compartment.
    for(int istep=0; istep < nstep +1; istep++ ){
        // Total number of monomers in aggs. in each compartment. used in the calculations of secondary nuc.
        vector<double> mdatay(ycompartment, 0); // initialze ycompartment data with initial value 0;
        vector< vector<double> > mdata2d(xcompartment, mdatay); // initialize xcompartment data, adatay for
        // Update mdata2d
        for( int x=0; x< xcompartment; x++ ){
            for(int y=0; y< ycompartment; y++){
                for(int il =1; il <= imax; il++){
                    mdata2d.at(x).at(y) = adata2d.at(il).at(x).at(y);
                }
                //cout << mdata2d.at(x).at(y) << endl;
            }
        }
        // update data holder
        if(outd.size() < adata2d.size()){
            int oldsize = outd.size();
            outd.resize(adata2d.size());
            for(int is=oldsize; is<outd.size(); is++){
                outd.at(is) = iagg;
            }
        }
        
        double r1 = 0, r2 = 0;
        r1 = random(rng) + 1e-30; // random indicator that will lead us to the compartment which will take action;
        r2 = random2(rng) + 1e-30;
        
        tpropensity = totalpropensity(adata2d, mdata2d, imax, params);
        if( tpropensity > 0 ){
            t = t + log(1/r1) / tpropensity;
            //cout << "t = " << t << endl;
        }
        else if( tpropensity == 0 ){
            // If no reaction is going to happen, write the current data repeatly until Nt and then stop. This will guarantee the same output data size.
            // cout << "At t = " << t << " and istep = " << istep << ", total propensity = 0." << endl;
            // ======== write adata2d into data holder =========================================
            for (int ic = count; ic <= Ndatapoint; ic++){
                for(int i = 0; i < imax+1; i++){
                    for(int x=0; x < adata2d.at(i).size(); x++){
                        for(int y=0; y < adata2d.at(i).at(x).size(); y++){
                            outd.at(i).at(ic) = adata2d.at(i).at(x).at(y);
                        }
                    }
                }
            }
            //  << "Propensity == 0; Simulation finished!" << endl;
            break;
        }
        else{
            cout << "Negative total propensity: " << tpropensity << endl;
            break;
        }
        
        // Check if the time sequence for storing data is suitable for the system
        if( t > time.at(Ndatapoint)*10 ){
            // cout << " Current delta_t = " << t << endl;
            // cout << " This system time interval is 10 times larger than the time interval for collecting data. You should change the Nt. " << endl;
            // ======== write adata2d into data files =========================================
            for (int ic = count; ic <= Ndatapoint; ic++){
                for(int i = 0; i < imax+1; i++){
                    for(int x=0; x < adata2d.at(i).size(); x++){
                        for(int y=0; y < adata2d.at(i).at(x).size(); y++){
                            outd.at(i).at(ic) = adata2d.at(i).at(x).at(y);
                        }
                    }
                }
            }
            break;
        }
        // ======== write adata2d into data files =========================================
        while( t >= time.at(count-1) and t >= time.at(count) ){
            for(int i = 0; i < imax+1; i++){
                for(int x=0; x < adata2d.at(i).size(); x++){
                    for(int y=0; y < adata2d.at(i).at(x).size(); y++){
                        outd.at(i).at(count) = adata2d.at(i).at(x).at(y);
                    }
                }
            }
            count = count + 1;
            if ( count == Ndatapoint+1 ){
                break; }
        }
        if ( count == Ndatapoint+1){
            // cout << "Simulation finished! istep = " << istep << endl;
            break; }
        // =================================================================================
        int i = -1; // the i-the index for adata2d.at(i). The length of an aggregate is i + nc - 1.
        double ipartialsum = 0, cpartialsum = 0;
        // cout << " Random = " << r2 * tpropensity << endl;
        // The big picture here is that we need to find out the l-aggregate is going to happen sth, and we start from determining which l-aggregate, then for the moving-to-the-right reaction, we locate the the compartment that is going to loose one l-aggregate to the one right next to it, then for the moving-to-the-left reaction......
        // First, we find out the l-aggregate
        while( ipartialsum <= r2 * tpropensity and i < imax ){
            i = i + 1;
            ipartialsum = ipartialsum + lpropensity(adata2d.at(0), adata2d.at(i), mdata2d, i, params);
            
        }
        // Now our random indicator stops at the l-aggregate, actually, lpartialsum is now larger than the random indicator (we just hopped over it), so we need to go back one l-step and move with a smaller step toward the indicator
        cpartialsum = ipartialsum - lpropensity(adata2d.at(0), adata2d.at(i), mdata2d, i, params);
        // For moving-to-the-right reaction
        int cxy = -1;
        while( cpartialsum <= tpropensity * r2 and cxy < xcompartment*ycompartment - ycompartment - 1 ){
            cxy = cxy + 1;
            cpartialsum = cpartialsum + adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment) * diffconst(i, params);
        }
        // cout << "cx+ " << cpartialsum << endl;
        // Now (cx = cxy/ycompartment, cy = cxy%ycompartment) is the compartment that one l-aggregate is moving to its right
        if( r2 * tpropensity < cpartialsum){
            adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment) - 1;
            adata2d.at(i).at(cxy/ycompartment + 1).at(cxy%ycompartment) = adata2d.at(i).at(cxy/ycompartment + 1).at(cxy%ycompartment) + 1;
//            if(adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment) < 0){
//                cout << "cx+" << endl;
//            }
        }
        else{
            // For moving-to-the-left reaction
            cxy = -1+ycompartment;
            while( cpartialsum <= tpropensity * r2 and cxy < xcompartment*ycompartment - 1 ){
                cxy = cxy + 1;
                cpartialsum = cpartialsum + adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment) * diffconst(i, params);
            }
            // cout << "cx- " << cpartialsum << endl;
            if( r2 * tpropensity < cpartialsum){
                adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment) - 1;
                adata2d.at(i).at(cxy/ycompartment - 1).at(cxy%ycompartment) = adata2d.at(i).at(cxy/ycompartment - 1).at(cxy%ycompartment) + 1;
//                if(adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment) < 0){
//                    cout << "cx-" << endl;
//                }
            }
            else{
                // For cy+ reaction
                int cxy = -1;
                while( cpartialsum <= tpropensity * r2 and cxy < xcompartment*ycompartment - xcompartment - 1 ){
                    cxy = cxy + 1;
                    cpartialsum = cpartialsum + adata2d.at(i).at(cxy%xcompartment).at(cxy/xcompartment) * diffconst(i, params);
                }
                //cout << "cy+ " << cpartialsum << endl;
                // Now (cx = cxy%ycompartment, cy = cxy/ycompartment) is the compartment that one l-aggregate is moving y+
                if( r2 * tpropensity < cpartialsum){
                    adata2d.at(i).at(cxy%xcompartment).at(cxy/xcompartment) = adata2d.at(i).at(cxy%xcompartment).at(cxy/xcompartment) - 1;
                    adata2d.at(i).at(cxy%xcompartment).at(cxy/xcompartment + 1) = adata2d.at(i).at(cxy%xcompartment).at(cxy/xcompartment + 1) + 1;
//                    if(adata2d.at(i).at(cxy%xcompartment).at(cxy/xcompartment) < 0){
//                        cout << "cy+" << endl;
//                    }
                }
                else{
                    // For cy- reaction
                    cxy = -1 + xcompartment;
                    while( cpartialsum <= tpropensity * r2 and cxy < xcompartment*ycompartment - 1 ){
                        cxy = cxy + 1;
                        cpartialsum = cpartialsum + adata2d.at(i).at(cxy%xcompartment).at(cxy/xcompartment) * diffconst(i, params);
                    }
                    //cout << "cy- " << cpartialsum << endl;
                    if( r2 * tpropensity < cpartialsum){
                        adata2d.at(i).at(cxy%xcompartment).at(cxy/xcompartment) = adata2d.at(i).at(cxy%xcompartment).at(cxy/xcompartment) - 1;
                        adata2d.at(i).at(cxy%xcompartment).at(cxy/xcompartment - 1) = adata2d.at(i).at(cxy%xcompartment).at(cxy/xcompartment - 1) + 1;
//                        if(adata2d.at(i).at(cxy%xcompartment).at(cxy/xcompartment) < 0){
//                            cout << "cy-" << endl;
//                        }
                    }
                    else{
                        // For add one process & primary nucleation
                        cxy = -1;
                        while( cpartialsum <= r2 * tpropensity and cxy < xcompartment*ycompartment - 1  ){
                            cxy = cxy + 1;
                            // Skip those compartment that has monomers less than primnuc
                            if(i==0){
                                if(adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) >= primnuc){
                                    cpartialsum = cpartialsum + pow(adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) , primnuc) * kn;
                                }
                                else{
                                    cpartialsum = cpartialsum;
                                }
                            }
                            else{
                                cpartialsum = cpartialsum + 2. * adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment) * adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) * kp;
                            }
                        }
                        //cout << i << " imax: " << imax << " add on " << cpartialsum << endl;
                        if( r2 * tpropensity < cpartialsum ){
                            if(i == 0){ // Primary nucleation
                                adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) - primnuc;
//                                if(adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) < 0){
//                                    cout << cxy/ycompartment << "\t" << cxy%ycompartment<< "\t" << "i = 0, addon" << endl;
//                                    return(0);
//                                }
                            }
                            else if(i < imax){ // Elongation
                                adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment) - 1;
                                adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) - 1;
//                                if(adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) < 0){
//                                    cout << "i =/= 0, addon" << endl;
//                                }
                            }
                            else{ 
                                imax = imax + 1;
                                // Resize the adata2d vector to include the i+1 aggregate
                                adata2d.resize(i + 2);
                                adata2d.at(i+1) = adataxy;
                                adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment) - 1;
                                adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) - 1;
//                                if(adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) < 0){
//                                    cout << "i = imax, addon" << endl;
//                                }
                            }
                            //
                            adata2d.at(i+1).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(i+1).at(cxy/ycompartment).at(cxy%ycompartment) + 1;
                        }
                        else{
                            // For fragmentation
                            cxy = -1;
                            while( cpartialsum <= tpropensity * r2 and cxy < xcompartment*ycompartment - 1 ){
                                cxy = cxy + 1;
                                cpartialsum = cpartialsum + (i+primnuc-2) * kf * adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment);
                            }
                            if( r2 * tpropensity < cpartialsum ){
//                                if(i == 0){cout << i << " Warning: Fragmentation with monomer or nuc " << endl;}
                                // cout << "frag " << cpartialsum << endl;
                                // cxy looses one l-aggregate
                                //cout << i+primnuc-1 << " -agg is being frag!!!! " <<endl;
                                adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment) - 1;
//                                if(adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment) < 0){
//                                    cout << "frag" << endl;
//                                }
                                uniform_int_distribution<int> random3(1,i+primnuc-2); // generate a random number between 1 and l-1 = i+nc-1-1
                                int j = random3(rng); // j is the location that the fragmentation occurs
				                //cout << "FFFFFFFFragment  " <<endl;
                                //cout << "i = " << i << "----, j = " << j << endl;
                                // if the breage happen at a location that results in a shorter-than-nc aggregate, then this smaller aggregate is assumed to dissolve into monomers immediately.
                                // l = i+nc-1; l1 = j; l2 = l-j;
                                if( j < primnuc and i+primnuc-1-j < primnuc ){//cout << "aa create total " << i+primnuc-1 << "monomers!"<<endl;
                                    adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) + i+primnuc-1;}
                                else if( j < primnuc and i+primnuc-1-j >= primnuc ){//cout << "bb create " << j << "monomers and one" << i-j + primnuc -1 << "-agg" <<endl;
                                    adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) + j;
                                    adata2d.at(i-j).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(i-j).at(cxy/ycompartment).at(cxy%ycompartment) + 1;}
                                else if( j >= primnuc and i+primnuc-1-j < primnuc ){//cout << "cc create " << i+primnuc-1-j << "monomers and one" << j << "-agg" <<endl;
                                    adata2d.at(j-primnuc+1).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(j-primnuc+1).at(cxy/ycompartment).at(cxy%ycompartment) + 1;
                                    adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) + i+primnuc-1-j;}
                                else{//cout << "dd create one " << j << "-agg and one" << i-j + primnuc -1 << "-agg" <<endl;
                                    adata2d.at(j-primnuc+1).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(j-primnuc+1).at(cxy/ycompartment).at(cxy%ycompartment) + 1;
                                    adata2d.at(i-j).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(i-j).at(cxy/ycompartment).at(cxy%ycompartment) + 1;}
                            }
                            else{
                                // For dissociation
                                // In case tpropensity = 0, there should be no reaction
                                cxy = -1;
                                while( cpartialsum <= tpropensity * r2 and cxy < xcompartment*ycompartment - 1 ){
                                    cxy = cxy + 1;
                                    cpartialsum = cpartialsum + 2. * km * adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment);
                                }
                                
                                if( r2 * tpropensity < cpartialsum ){
//                                    if(i == 0 or i == 1){cout << i << " Warning: Dissociation with monomer or nuc " << endl;}
//                                    if(km == 0){ cout << " warning: km=0 but dissociation is happening!"<<endl; }
                                    adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(i).at(cxy/ycompartment).at(cxy%ycompartment) - 1;
                                    adata2d.at(i-1).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(i-1).at(cxy/ycompartment).at(cxy%ycompartment) + 1;
                                    adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) + 1;
                                }
                                else{
                                    // For secondary nucleation
                                    cxy = -1;
                                    while( cpartialsum <= tpropensity * r2 and cxy < xcompartment*ycompartment - 1 ){
                                        cxy = cxy + 1;
                                        if(adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) >= n2){
                                            cpartialsum = cpartialsum + k2 * pow(adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) , n2)* mdata2d.at(cxy/ycompartment).at(cxy%ycompartment);
                                        }
                                    }
                                    if( r2 * tpropensity < cpartialsum ){
//                                        if(i != 0){cout << i << " Warning: Secondary nucleation" << endl;}
//                                        if(k2 == 0){ cout << " warning: k2=0 but sec nuc is happening!"<<endl; }
                                        adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) - n2;
//                                        if(adata2d.at(0).at(cxy/ycompartment).at(cxy%ycompartment) < 0){
//                                            cout << cxy/ycompartment << "\t" << cxy%ycompartment<< "\t" << "i = 0, sec nuc" << endl;
//                                            return(0);
//                                        }
                                        if(imax >= n2+1-primnuc){
                                            adata2d.at(n2+1-primnuc).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(n2+1-primnuc).at(cxy/ycompartment).at(cxy%ycompartment) + 1;
                                        }
                                        else{
                                            int imax_temp = 0;
                                            imax_temp = imax;
                                            imax = n2+1-primnuc;
                                            // Resize the adata2d vector to include the i+1 aggregate
                                            adata2d.resize(imax + 1);
                                            for(int jj=imax_temp+1; jj<=imax; jj++ ){
                                                adata2d.at(jj) = adataxy;
                                            }
                                            adata2d.at(imax).at(cxy/ycompartment).at(cxy%ycompartment) = adata2d.at(imax).at(cxy/ycompartment).at(cxy%ycompartment) + 1;
                                        }
                                    }
//                                    else{cout<<" Total Propensities do not match! " << endl;}
                                }
                            }
                        }
                    }
                }
            }
        }
//        if(istep == nstep){cout << "Finished with nstep!" << endl;}
    }// end of for loop
    cout << "Max length: " << primnuc << " + " <<  adata2d.size() - 2 << endl;
   // maxlforsample[isample-1] = int(adata2d.size() - 1);
    int* ml = new int;
    *ml = int(adata2d.size() - 1);
    //cout << " Final: " << adata2d.size() << " , " << *ml << ", " << imax << endl;
    
    string fname = datafilepath + "adata_" + to_string(isample) + ".dat";
    write2file(fname, outd);
    
    adata2d.clear();
    outd.clear();
    pthread_exit( (void *)ml );
}


