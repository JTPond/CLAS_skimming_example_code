/*
    As simple as possible CLAS merging example:
        gamma P -> N* omega -> n pi+ pi+pi-(pi0)
        gamma P -> Delta* omega -> n pi+ pi+pi-(pi0)
    Author: Josh Pond
    Email: jp4cj@virginia.edu, jpond@jlab.org
*/
//C++
#include <iostream>
#include <cmath>
#include <vector>

//Import CLAS6 headers
#include <Vec.h>
#include <lorentz.h>
#include <pputil.h>
#include <clasEvent.h>

//Import root headers
#include "TString.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TBranch.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TH2D.h"
#include "TH1F.h"


using namespace std;


int main(int __argc,char *__argv[]){
// Allocate hold variables
    Float_t run,evt,mm2Pi0,gammaE,mpX,mpY,mpZ,mpR,mpT,tranmp,S,T,U,pipPhi,pipTheta,pip2Phi,pip2Theta,pimPhi,pimTheta,neuPhi,neuTheta,mOmega,mNstar,mm2neuX,mm2INminusNstar,pipr,pip2r,pimr,neur,dTpip,dTpip2,dTpim,dTneu,neuX,neuY,neuZ,neuE,pipX,pipY,pipZ,pipE,pip2X,pip2Y,pip2Z,pip2E,pimX,pimY,pimZ,pimE,evtX,evtY,evtZ,cTOmega;

// Allocate new matraces
    Float_t pip_mat[5][5];
    Float_t pim_mat[5][5];
    Float_t pip2_mat[5][5];
    Float_t neu_mat[5][5];

// Allocate hold matreces
    Float_t pip_mat_n[5][5];
    Float_t pim_mat_n[5][5];
    Float_t pip2_mat_n[5][5];
    Float_t neu_mat_n[5][5];

// Read input file names from the command line
    Int_t Argc = __argc;
    char **Input = __argv;
    std::vector<TString> w;
    w.assign(__argv, __argv + __argc);
    TApplication* theApp = new TApplication("App", &__argc, __argv);
   
    extern int optind;

// Open new file 
    TFile *outfile = new TFile((char *) "g12_N_pmN_merged.root","recreate");

// Get new TNtuple and new mat branches
    const int nData = 53;
    Float_t nt_value[nData];
    TNtuple *pp = new TNtuple("nStar","nS -> n* omega -> PPi- Pi+Pi-Pi0","run:evt:mm2Pi0:gammaE:mpX:mpY:mpZ:mpR:mpT:tranmp:S:T:U:pipPhi:pipTheta:pip2Phi:pip2Theta:pimPhi:pimTheta:neuPhi:neuTheta:mOmega:mNstar:mm2neuX:mm2INminusNstar:pipr:pip2r:pimr:neur:dTpip:dTpip2:dTpim:dTneu:neuX:neuY:neuZ:neuE:pipX:pipY:pipZ:pipE:pip2X:pip2Y:pip2Z:pip2E:pimX:pimY:pimZ:pimE:evtX:evtY:evtZ:cTOmega");
    pp->Branch("pip_mat",pip_mat,"pip_mat[5][5]/F");
    pp->Branch("pip2_mat",pip2_mat,"pip2_mat[5][5]/F");
    pp->Branch("pim_mat",pim_mat,"pim_mat[5][5]/F");
    pp->Branch("neu_mat",neu_mat,"neu_mat[5][5]/F");


    for(int n_arg = optind; n_arg < Argc; n_arg++){ //Loop over input files
        TString input =w[n_arg];
        TFile inFile(input); // open the input file               

        if(TTree *tree = (TTree*)inFile.Get("nStar")){ //Check for the right input tree name

        // Get all branch addresses (i.e associate each branch with a variable memory address)
            tree->SetBranchAddress("run",&run);
            tree->SetBranchAddress("evt",&evt);
            tree->SetBranchAddress("mm2Pi0",&mm2Pi0);
            tree->SetBranchAddress("gammaE",&gammaE);
            tree->SetBranchAddress("mpX",&mpX);
            tree->SetBranchAddress("mpY",&mpY);
            tree->SetBranchAddress("mpZ",&mpZ);
            tree->SetBranchAddress("mpR",&mpR);
            tree->SetBranchAddress("mpT",&mpT);
            tree->SetBranchAddress("tranmp",&tranmp);
            tree->SetBranchAddress("S",&S);
            tree->SetBranchAddress("T",&T);
            tree->SetBranchAddress("U",&U);
            tree->SetBranchAddress("pipPhi",&pipPhi);
            tree->SetBranchAddress("pipTheta",&pipTheta);
            tree->SetBranchAddress("pip2Phi",&pip2Phi);
            tree->SetBranchAddress("pip2Theta",&pip2Theta);
            tree->SetBranchAddress("pimPhi",&pimPhi);
            tree->SetBranchAddress("pimTheta",&pimTheta);
            tree->SetBranchAddress("neuPhi",&neuPhi);
            tree->SetBranchAddress("neuTheta",&neuTheta);
            tree->SetBranchAddress("mOmega",&mOmega);
            tree->SetBranchAddress("mNstar",&mNstar);
            tree->SetBranchAddress("mm2neuX",&mm2neuX);
            tree->SetBranchAddress("mm2INminusNstar",&mm2INminusNstar);
            tree->SetBranchAddress("pipr",&pipr);
            tree->SetBranchAddress("pip2r",&pip2r);
            tree->SetBranchAddress("pimr",&pimr);
            tree->SetBranchAddress("neur",&neur);
            tree->SetBranchAddress("dTpip",&dTpip);
            tree->SetBranchAddress("dTpip2",&dTpip2);
            tree->SetBranchAddress("dTpim",&dTpim);
            tree->SetBranchAddress("dTneu",&dTneu);
            tree->SetBranchAddress("neuX",&neuX);
            tree->SetBranchAddress("neuY",&neuY);
            tree->SetBranchAddress("neuZ",&neuZ);
            tree->SetBranchAddress("neuE",&neuE);
            tree->SetBranchAddress("pipX",&pipX);
            tree->SetBranchAddress("pipY",&pipY);
            tree->SetBranchAddress("pipZ",&pipZ);
            tree->SetBranchAddress("pipE",&pipE);
            tree->SetBranchAddress("pip2X",&pip2X);
            tree->SetBranchAddress("pip2Y",&pip2Y);
            tree->SetBranchAddress("pip2Z",&pip2Z);
            tree->SetBranchAddress("pip2E",&pip2E);
            tree->SetBranchAddress("pimX",&pimX);
            tree->SetBranchAddress("pimY",&pimY);
            tree->SetBranchAddress("pimZ",&pimZ);
            tree->SetBranchAddress("pimE",&pimE);
            tree->SetBranchAddress("evtX",&evtX);
            tree->SetBranchAddress("evtY",&evtY);
            tree->SetBranchAddress("evtZ",&evtZ);
            tree->SetBranchAddress("missingOmega",&cTOmega);
            tree->SetBranchAddress("pip_mat_n",&pip_mat_n); 
            tree->SetBranchAddress("pim_mat_n",&pim_mat_n);
            tree->SetBranchAddress("pip2_mat_n",&pip2_mat_n);
            tree->SetBranchAddress("neu_mat_n",&neu_mat_n);
            Int_t entries= tree->GetEntries(); 
            for (int i=0; i<entries; i++){ // Loop over entries in file
                tree->GetEntry(i);// Criticle!!! 
                // Take value out of the hold variable and put it in the new array
		            nt_value[0] = run;
                    nt_value[1] = evt;
                    nt_value[2] = mm2Pi0;
                    nt_value[3] = gammaE;
                    nt_value[4] = mpX;
                    nt_value[5] = mpY;
                    nt_value[6] = mpZ;
                    nt_value[7] = mpR;
                    nt_value[8] = mpT;
                    nt_value[9] = tranmp;
                    nt_value[10] = S;
                    nt_value[11] = T;
                    nt_value[12] = U;
                    nt_value[13] = pipPhi;
                    nt_value[14] = pipTheta;
                    nt_value[15] = pip2Phi;
                    nt_value[16] = pip2Theta;
                    nt_value[17] = pimPhi;
                    nt_value[18] = pimTheta;
                    nt_value[19] = neuPhi;
                    nt_value[20] = neuTheta;
                    nt_value[21] = mOmega;
                    nt_value[22] = mNstar;
                    nt_value[23] = mm2neuX;
                    nt_value[24] = mm2INminusNstar;
                    nt_value[25] = pipr;
                    nt_value[26] = pip2r;
                    nt_value[27] = pimr;
                    nt_value[28] = neur;
                    nt_value[29] = dTpip;
                    nt_value[30] = dTpip2;
                    nt_value[31] = dTpim;
                    nt_value[32] = dTneu;
                    nt_value[33] = neuX;
                    nt_value[34] = neuY;
                    nt_value[35] = neuZ;
                    nt_value[36] = neuE;
                    nt_value[37] = pipX;
                    nt_value[38] = pipY;
                    nt_value[39] = pipZ;
                    nt_value[40] = pipE;
                    nt_value[41] = pip2X;
                    nt_value[42] = pip2Y;
                    nt_value[43] = pip2Z;
                    nt_value[44] = pip2E;
                    nt_value[45] = pimX;
                    nt_value[46] = pimY;
                    nt_value[47] = pimZ;
                    nt_value[48] = pimE;
                    nt_value[49] = evtX;
                    nt_value[50] = evtY;
                    nt_value[51] = evtZ;
                    nt_value[52] = cTOmega;
                    for(int j = 0; j<5; j++){
                        for(int k = 0; k<5; k++){
                            pip_mat[j][k] = pip_mat_n[j][k];
                            pim_mat[j][k] = pim_mat_n[j][k];
                            pip2_mat[j][k] = pip2_mat_n[j][k];
                            neu_mat[j][k] = neu_mat_n[j][k]; 
                        }
                    }
                    outfile->cd();
                    pp->Fill(nt_value); // Fill the TNtuple
            }//for entry
        }//if correct tree
    }//for file
    pp->Write(); // Write the tuple to the file
    outfile->Close();

    return 0;
   
}//main

