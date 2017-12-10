#include <iostream>
#include <vector>
#include <cstdio>
#include "utility.h"
#include "operator.h"
#include <random>
#include <chrono>


using namespace std;

long acceptance = 0;
int nBeads = 10;
int nParticles = 1;
long nSteps = 500000; //number of MC steps
double simulation(int nBeads);
long genmoves = 0;



/**
 *
 *
 * @param p : probability - function
 * @param init : initial values of input parameters
 * @param nSteps : number of marcov chain monte carlo steps
 * @return : vector of size nSteps with values
 */
vector<vector<vector<double> > > MCMC(double (*prob)(vector<vector<double> >,vector<vector<double> >, int, int),  double alpha, long nSteps) //propagator for the whole system between time t and t+1
{
    std::mt19937_64 rng;
    // initialize the random number generator with time-dependent seed
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(ss);
    // initialize a uniform distribution between 0 and 1
    std::uniform_real_distribution<double> unif(0, 1);
    // ready to generate random numbers

    vector<vector<vector<double> > > R;
    R.assign(nParticles,vector< vector<double> > (nBeads,vector<double>(3,0)));
    double oxd=1;
    //R=pos_init(R,0,0,0,1,nBeads);//R, x, y, z, particle index, nBeads
    //R=pos_init(R,oxd,0,0,2,nBeads);
    R=pos_init(R,oxd/2,0,0,0,nBeads);
    double plast = 0;
    double pnext = 0;
    int k = 0;
    int p = 0;
    vector<vector<vector<double> > > MChain;
    MChain.assign(nSteps, vector< vector<double> >(nParticles,vector<double>(3,0)));
    vector<vector<double> > next_part(nBeads,vector<double>(3,0)) ; // vector of 3d beads for particle 0
    vector<vector<double> > next_beads(nParticles,vector<double>(3,0)) ; // vector of 3d coordinates of a bead for all particles
    vector<vector<double> > last_part(nBeads,vector<double>(3,0));
    vector<vector<double> > last_beads(nParticles,vector<double>(3,0));

    int l = 0;
    while (l<nSteps)
    {
        pnext = 0;
        plast = 0;
        if (l>nSteps)
        {
            k = static_cast<int>((nBeads + 1) * unif(rng));   //select a random bead
        } else {
            k = static_cast<int>((nBeads) * unif(rng));
        }
        last_part=vec_ext(R,0,0,nBeads);
        last_beads=vec_ext(R,k,1,nParticles);
        next_part=last_part;
        next_beads=last_beads;

        double dx=0.01*(-1+2*unif(rng));
        //double dy=(-1+2*unif(rng));
        //double dz=0.01*(-1+2*unif(rng));
        if (k == nBeads)
        {
            next_part.at(p).at(0) += alpha * dx;//divide displacement by alpha when
            //next_part.at(p).at(1) += alpha * dy;// moving all beads at once
            //next_part.at(p).at(2) += alpha * dz;
            for (int j = 0; j < nBeads; ++j) {

                next_beads.at(j).at(0) += alpha * dx;
                //next_beads.at(j).at(1) += alpha * dy;
                //next_beads.at(j).at(2) += alpha * dz;
            }
            pnext += prob(next_part, next_beads, k, 0);
            plast += prob(last_part, last_beads,k,0);
            plast = plast / nBeads;
            pnext = pnext / nBeads;

        }

        else
        {


            next_part.at(k).at(0) += dx;
            //next_part.at(k).at(1) += dy;
            //next_part.at(k).at(2) += dz;
            next_beads.at(p).at(0) += dx;
            //next_beads.at(p).at(1) += dy;
            //next_beads.at(p).at(2) += dz;

            pnext = prob(next_part, next_beads,k,0);
            plast = prob(last_part, last_beads,k,0);
            plast = plast / nBeads;
            pnext = pnext / nBeads;
            printf("pnext %f  plast %f \n",pnext, plast);
        }
        if (pnext > plast || unif(rng)<pnext/plast)
        {
            acceptance++;
            last_beads = next_beads;
            last_part = next_part;
            R.at(p).at(k).at(0)+=dx;
            //R.at(p).at(k).at(1)+=dy;
            //R.at(p).at(k).at(2)+=dz;
        }
        if (k == nBeads)
        {
            for (int j = 0; j < nBeads; ++j)
            {
                MChain.at(l)= last_beads;
                for (int i = 0; i < nParticles; ++i) {
                    double d=distance(MChain.at(l).at(i));
                    //printf("%i %i %f %f %f %f \n",l,i,MChain.at(i).at(l).at(0),MChain.at(i).at(l).at(1),MChain.at(i).at(l).at(2),d);
                }
                l=l+1;
            }
            genmoves +=1;
        }
        else
        {
            MChain.at(l)=last_beads;
            for (int i = 0; i < nParticles; ++i) {

                double d = distance(MChain.at(l).at(i));
                //printf("%i %i  %f %f %f %f \n",l,i,MChain.at(l).at(i).at(0),MChain.at(l).at(i).at(1),MChain.at(l).at(i).at(2),d);
            }
            l=l+1;
        }
    }
    return MChain;
}

double simulation(int nBeads) {


    //srand( (unsigned int)( time(nullptr) ) );
    ofstream myfile;
    //char fn [100];
    //snprintf (fn, sizeof fn, "simulation_nBeads_%02d.txt", nBeads);
    myfile.open("run1.txt");



    //set up number of steps, number of beads, initial values, delta

    double alpha = 0.1;


    vector<vector<vector<double> > > MCserie;
    MCserie.assign(nSteps, vector< vector<double> >(nParticles,vector<double>(3,0)));

    //calculate markov chain
    MCserie = MCMC(&prob,alpha,nSteps);

    //save markov chain into file




    //beginning of post analysis

    double avg;
    avg = 0;
    double avgsq = 0;
    double sigma = 0;
    double sum = 0;
    for(int i=0;i<nSteps;i++) {
        //    for (int j = 0; j < nParticles; ++j) {
        //vector<double> ox1(3,2);
        //vector<double> ox2(3,0);
        avg = runavg(i, distance(MCserie.at(i).at(0)), avg);
        avgsq = runavg(i, distance(MCserie.at(i).at(0)) * distance(MCserie.at(i).at(0)), avgsq);
        sigma = avgsq - (avg * avg);
        sum = sum + MCserie.at(i).at(0).at(0);
        if (i % (100 * nBeads) == 0) {
            myfile << i << " " << MCserie.at(i).at(0).at(0) << " " << MCserie.at(i).at(0).at(1) << " "
                   << MCserie.at(i).at(0).at(2) << " " << sum << " " << avg << " " << sigma << endl;
            //printf("%i %f %f \n",i,avg,sigma);
            //fprintf("%d \n",genmoves);
            //    }
        }
    }
    myfile.close();


    //end of post analysis



    //print out averages
    //   printf("  N   = %d\n",nBeads);
    //   printf(" <X>  = %f\n",average(MCserie));
    //   printf("<X*X> = %f\n",averageOfSquares(MCserie));


    //   return avg;
    return 0;

}






//////////////////////////////////////////////////////////////////////////////////////////////////////////////





int main(){


    /*
     for (int i = 0; i < 20 ; ++i)
     {
         double j = double(i);
         j = 10*pow(2,j);
         string filenumber = to_string(j);
         string fn = "simulation_nBeads" + filenumber + ".txt";              //ereasing previous files (do not work)
         const char* filename = fn.c_str();
         remove(filename);
     }
      */
    double threshold = 0.05 ;
    double avg1 = simulation(nBeads);
    double avg2 = 0;
    double diff = 1;
//    while (diff>threshold){
//        nSteps=nSteps*2;
//        //avg2=simulation(nBeads);
//        //diff = abs(avg1-avg2);
//        avg1=avg2;
//        nBeads=2*nBeads;
//        diff=0;
//    }

    /*
vector<vector<int> > heatmap;
    heatmap = DataAnalysis2D(heatmap);
    vector<double> Free_E_distrib;
    Free_E_distrib2D = Free_E2D(heatmap);
*/
vector<int> bargraph;
 bargraph = DataAnalysis1D(bargraph);
    for (i=0; i<division; i++) {
    printf("bargraph: %i \n", bargraph.at(i));
    }
    printf("bargraph size%i \n", bargraph.size());

      Free_E_distrib = Free_E(bargraph);
    printf("acceptance = %i",acceptance);
    return 0;
}
