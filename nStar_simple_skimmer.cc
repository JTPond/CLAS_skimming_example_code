/*
    As simple as possible CLAS skimming example:
        gamma P -> N* omega -> n pi+ pi+pi-(pi0)
        gamma P -> Delta* omega -> n pi+ pi+pi-(pi0)
    Author: Josh Pond, Will Phelps, and jb. 
    Email: jp4cj@virginia.edu, jpond@jlab.org
*/
#include <iostream>
#include <cmath>

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

//Import corrections header, modify for specific needs
#include "/home/ptgroup/g12_corrections/All_Corrections/g12_corrections.hpp"

//Global definitions
TFile *f;
TNtuple *pp; //This example saves to a TNtuple, not a TTree. This can be modified
const int nData = 52;
Float_t nt_value[nData]; // The array that will fill the Ntuple
Float_t pip_mat[5][5];  // The arrays that will hold the TBER matrices. 
Float_t pip2_mat[5][5];
Float_t pim_mat[5][5];
Float_t neu_mat[5][5];

clas::g12::MomentumCorrection pcor; //Create a momentum correction object

#define RAD_TO_DEGREE 57.295820614 // 180.0/pi

//The function that returns a momentum corrected four vector, and was written by W. Phelps
fourVec pcorFourVec_new( fourVec particle, int part_id ){
    float proton_mass = 0.938272; //GeV
    float pion_mass = 0.139570; //GeV
    float kaon_mass = 0.493667; //GeV
    //Assign mass based on pid.
    double PID_mass =-1000;
    PID_mass = (part_id == 14                )?proton_mass:PID_mass;
    PID_mass = (part_id == 11|| part_id == 12)?kaon_mass  :PID_mass;
    PID_mass = (part_id == 8 || part_id == 9 )?pion_mass  :PID_mass;

        float id = part_id;
        float p =  particle.r(); // total momentum
        float px =  particle.x();
        float py =  particle.y();

        float pz = sqrt(p*p - px*px - py*py);
        float phi = particle.V().phi();
        float theta = acos(pz/p);

        /// Momentum correction ////////////////////////////////////////
        float new_p = p + pcor.pcor(phi,id);
        /// ////////////////////////////////////////////////////////////

        float new_px = (new_p / p) * px;
        float new_py = (new_p / p) * py;
        float new_pz = (new_p / p) * pz;

    double pE = sqrt( PID_mass*PID_mass + new_p*new_p );
    fourVec correctedFourVec;
    correctedFourVec.set(pE, new_px, new_py, new_pz);

    return correctedFourVec;
}

//returns sector (int) of the CLAS drift chambers with inut of phi (float)
int getSector(float phi){
    int sector = 0;

    if(std::abs(phi) <= 30.) sector = 1;
    else if(phi > 0.){
        if(phi <= 90.) sector = 2;
        else if(phi <= 150) sector = 3;
        else sector = 4;
    }
    else {
    // phi < 0
        if(std::abs(phi) <= 90.) sector = 6;
        else if(std::abs(phi) <= 150.) sector = 5;
        else sector = 4;
    }
    return sector;
}

// This function returns the TAGR bank, and was written by "jb", 
tagr_t *my_get_photon_tagr(clasTAGR_t *TAGR,clasBID_t *TBID, threeVec myVert){
    /* This routine only works for time-based tracks! */
    float best_diff=ST_TAG_COINCIDENCE_WINDOW;
    float tprop=0.0;
    tagr_t *tagr = NULL;
    clasTBTR_t *TBTR=(clasTBTR_t *)getBank(&bcs_,"TBTR");
    float g11targz=-10; //this was for g11
    float g12targz=-90; //jb v11
    int i, j;
    /* Exit from function if missing the requisite banks... */
    if(!TAGR || !TBTR || !TBID) return(NULL);
    for (i=0;i<TBID->bank.nrow;i++){
        int trk_ind=TBID->bid[i].track;
        if(trk_ind){
            tprop=(myVert.z() - g12targz)/LIGHT_SPEED;  //jb v13
            if (TBID->bid[i].st.stat){
                float mytime=-1000; //jb v14
                float myenergy=-1000; //jb v14
                for(j=0;j<TAGR->bank.nrow;j++){
                    float diff=fabs(TBID->bid[i].st.vtime-(TAGR->tagr[j].tpho+tprop));
                    if (( diff<ST_TAG_COINCIDENCE_WINDOW&&diff<best_diff || (abs(TAGR->tagr[j].tpho - mytime)<0.011&& (TAGR->tagr[j].erg) > myenergy)  )  &&  (TAGR->tagr[j].stat==7 || TAGR->tagr[j].stat==15)){  //jb v14
                        best_diff=diff;
                        tagr=&(TAGR->tagr[j]);
                        mytime=tagr->tpho; //jb v14
                        myenergy=tagr->erg; //jb v14
                    }
                }
            }
        }
    }
    return(tagr);
}

//This is the main Event analyzer
void ProcessEvent(clasEvent &evt){
	fourVec pip, pip2, pim, neu; //These are CLAS 4 momenta / 4 vectors for the analyzed particles. They behave in the expected way.
	clasParticle pipx, pip2x, pimx, neux; //These are CLAS particle objects, they have many useful members.
	fourVec beam,beam1,target;
	fourVec X;
	double eBeamLo = 4.4;
	double eBeamHi = 7.0;
    float tpho; //Photon energy
    target = evt.target().get4P();

    clasTAGR_t *TAGR = (clasTAGR_t *)getBank(&bcs_,"TAGR"); //These are the BOS banks 
    clasBID_t *TBID  =  (clasBID_t *)getBank(&bcs_,"TBID");
    clasTBER_t *TBER = (clasTBER_t *)getBank(&bcs_,"TBER");

    threeVec evtVert=evt.V(); //event vertex 
    tagr_t *tagr = my_get_photon_tagr(TAGR, TBID, evtVert);
    double beamp = -1000.0;
    if(tagr){
        beamp = tagr->erg;
        tpho=tagr->tpho;
    }
    beam1.set(beamp,threeVec(0.0,0.0,beamp)); //set the initial tagged photon

    if(evt.N(PiPlus)==2 && evt.N(PiMinus)==1 && evt.N(Neutron)==1){ //Part Bank PID cut

        evt.eLoss((char*)"g11a",1); //Event energy loss correction
        double ebeamcor = clas::g12::corrected_beam_energy(evt.run(),beam1.t()); //Beam energy correction

        beam.set(ebeamcor,threeVec(0.0,0.0,ebeamcor)); //set corrected photon

        pipx = evt.cp(PiPlus, 1); //Get the 4 momenta of the particles
        pip2x = evt.cp(PiPlus, 2);
        pimx = evt.cp(PiMinus,1);
        neux = evt.cp(Neutron,1);

        pip = pcorFourVec_new(pipx.p(), 8); // momentum corrections (takes a 4 momenta and the particle id number)
        pip2 = pcorFourVec_new(pip2x.p(), 8);
        pim = pcorFourVec_new(pimx.p(), 9);
        neu = neux.p(); //no momentum corrections on neutral particles (See presentation)

        //PION SORTING, i.e. make a decision as to which pion came from the omega, and which from the baryon

        double delOmega1 = std::abs(.782 - sqrt((beam + target - pip - neu).lenSq()));
        double delOmega2 = std::abs(.782 - sqrt((beam + target - pip2 - neu).lenSq()));

        if(delOmega1 > delOmega2){
            clasParticle p_temp = pip2x;
            pip2x = pipx;
            pipx = p_temp;
            pip = pcorFourVec_new(pipx.p(), 8);
            pip2 = pcorFourVec_new(pip2x.p(), 8);
        }

        X = beam + target - pip - pip2 - pim - neu; //get missing 4 momenta

        //Time of flight panel knock out cuts (see presentation)
        Int_t TOFKOpip = clas::g12::pass_g12_TOF_knockout(pip.theta()*RAD_TO_DEGREE,pip.phi()*RAD_TO_DEGREE,getSector(pip.phi()*RAD_TO_DEGREE));
        Int_t TOFKOpip2 = clas::g12::pass_g12_TOF_knockout(pip2.theta()*RAD_TO_DEGREE,pip2.phi()*RAD_TO_DEGREE,getSector(pip2.phi()*RAD_TO_DEGREE));
        Int_t TOFKOpim = clas::g12::pass_g12_TOF_knockout(pim.theta()*RAD_TO_DEGREE,pim.phi()*RAD_TO_DEGREE,getSector(pim.phi()*RAD_TO_DEGREE));

        //Fiducial cuts (see presentation)
        Float_t FCNpip = clas::g12::g12_PosParticle_fiducial_cuts(pip.r(),pip.theta()*RAD_TO_DEGREE, pip.phi()*RAD_TO_DEGREE,(char*)"nominal");
        Float_t FCNpip2 = clas::g12::g12_PosParticle_fiducial_cuts(pip2.r(),pip2.theta()*RAD_TO_DEGREE, pip2.phi()*RAD_TO_DEGREE,(char*)"nominal");
        Float_t FCNpim = clas::g12::g12_NegParticle_fiducial_cuts(pim.r(),pim.theta()*RAD_TO_DEGREE, pim.phi()*RAD_TO_DEGREE,(char*)"nominal");

        // vertex components
        Float_t evtX = evt.V().x();
        Float_t evtY = evt.V().y();
        Float_t evtZ = evt.V().z();

        // get vertex time of each particle
        Float_t pip_ttag = pipx.stVtime();
        Float_t pip2_ttag = pip2x.stVtime();
        Float_t pim_ttag = pimx.stVtime();

        // get event vertex time
        Float_t evt_stt = evt.vtime();

        Float_t pipbeta = pipx.beta(); //get particle beta (v/c) calculated two ways
        Float_t pip2beta = pip2x.beta(); // beta: calculated from measured momentum
	    Float_t pimbeta = pimx.beta();
	    Float_t neubeta = neux.beta();
	    Float_t pipBeta = pipx.Beta(); // Beta: calculated from timing 
	    Float_t pip2Beta = pip2x.Beta();
	    Float_t pimBeta = pimx.Beta();
	    Float_t neuBeta = neux.Beta();


        //get a delta T, or the difference between the measured and calculated time of flight for each particle
        Float_t M = 0.139569, c=29.9792458;
	    Float_t dTpip = ((pipx.scTOF()) - (pipx.scPathLen()/c)*sqrt(1.0 + pow((M/pip.r()),2)));
	    Float_t dTpip2 = ((pip2x.scTOF()) - (pip2x.scPathLen()/c)*sqrt(1.0 + pow((M/pip2.r()),2)));
	    Float_t dTpim = ((pimx.scTOF()) - (pimx.scPathLen()/c)*sqrt(1.0 + pow((M/pim.r()),2)));
	    Float_t dTneu = ((neux.scTOF()) - (neux.scPathLen()/c)*sqrt(1.0 + pow((.939565/neu.r()),2)));

        // impose cuts. If an event does not meet the requirements in this if statement, then it is not stored. 
        if (TOFKOpip && TOFKOpip2 && TOFKOpim && FCNpip != 0.0 && FCNpip2 != 0.0 && FCNpim != 0.0
            && (evtX*evtX + evtY*evtY) < 4.0 && evtZ > -110.0 && evtZ < -70.0
            && abs(pip_ttag - evt_stt) < 1.002 && abs(pip2_ttag - evt_stt) < 1.002 && abs(pim_ttag - evt_stt) < 1.002
            && (pipbeta-pipBeta) > -.02 && (pipbeta-pipBeta) < .02 && (pip2beta-pip2Beta) > -.02 && (pip2beta-pip2Beta) < .02 && (pimbeta-pimBeta) > -.02 && (pimbeta-pimBeta) < .02
            && (dTpip*dTpip + dTpim*dTpim) < 4.0 ){

	        //Missing Mass Squared / Mass
            nt_value[0] = evt.run();
            nt_value[1] = evt.event();
            nt_value[2] = (beam + target - pip - pip2 - pim - neu).lenSq(); //mm2
            nt_value[3] = ebeamcor;
            //Missing Momentum
            nt_value[4] = X.x();
            nt_value[5] = X.y();
            nt_value[6] = X.z();
            nt_value[7] = X.t();
            nt_value[8] = X.r();
            nt_value[9] = sqrt(X.x()*X.x() + X.y()*X.y());
            //Mandelstrum variables 
            nt_value[10] = (beam + target).lenSq(); //S
            nt_value[11] = (target - neu - pip).lenSq(); //T
            nt_value[12] = (beam - neu - pip).lenSq();//U
            //Angular components 
            nt_value[13] = pip.phi();
            nt_value[14] = pip.theta();
            nt_value[15] = pip2.phi();
            nt_value[16] = pip2.theta();
            nt_value[17] = pim.phi();
            nt_value[18] = pim.theta();
            nt_value[19] = neu.phi();
            nt_value[20] = neu.theta();
            // Invariant mass squared 
            nt_value[21] = (pip2 + pim + X).lenSq(); //omega
            nt_value[22] = (pip + neu).lenSq(); //N*
            nt_value[23] = (beam + target - pip - pip2 - pim).lenSq(); //mm neu ^2
            nt_value[24] = (beam + target - neu - pip).lenSq(); // missing omega 
            //total momentum
	        nt_value[25] = pip.r();
	        nt_value[26] = pip2.r();
	        nt_value[27] = pim.r();
	        nt_value[28] = neu.r();
            // delta T
            nt_value[29] = dTpip;
            nt_value[30] = dTpip2;
            nt_value[31] = dTpim;
            nt_value[32] = dTneu;
            // 4 momenta components
            nt_value[33] = neu.x();
            nt_value[34] = neu.y();
            nt_value[35] = neu.z();
            nt_value[36] = neu.t();
            nt_value[37] = pip.x();
            nt_value[38] = pip.y();
            nt_value[39] = pip.z();
            nt_value[40] = pip.t();
            nt_value[41] = pip2.x();
            nt_value[42] = pip2.y();
            nt_value[43] = pip2.z();
            nt_value[44] = pip2.t();
            nt_value[45] = pim.x();
            nt_value[46] = pim.y();
            nt_value[47] = pim.z();
            nt_value[48] = pim.t();
            // vertex components
            nt_value[49] = evt.V().x();
            nt_value[50] = evt.V().y();
            nt_value[51] = evt.V().z();
            // TBER matrices, used later
            matrix<double> pipErrMatrix = pipx.TBERmatrix();
            matrix<double> pip2ErrMatrix = pip2x.TBERmatrix();
            matrix<double> pimErrMatrix = pimx.TBERmatrix();
            matrix<double> neuErrMatrix = neux.TBERmatrix();
            for(int j = 0; j<5; j++){
                for(int k = 0; k<5; k++){
                    pip_mat[j][k] = pipErrMatrix.el(j,k);
                    pip2_mat[j][k] = pip2ErrMatrix.el(j,k);
                    pim_mat[j][k] = pimErrMatrix.el(j,k);
                    neu_mat[j][k] = neuErrMatrix.el(j,k);
                }
            }
            f->cd(); //access the Tfile
            pp->Fill(nt_value);// fill the TNtuple

	    }//cuts
    }//initial PID
}

int main(int argc,char **argv){
    
    // creates static output file named outPutFile.root, will be renamed to a unique name when Auger moves it to disk (see presentation)
    f = new TFile((char *) "outPutFile.root","recreate");

    // create the TNtuple object, (name, description, identifier for every branch)
    pp = new TNtuple("nStar","nS -> n* omega -> nPi+Pi+Pi-Pi0","run:evt:mm2Pi0:gammaE:mpX:mpY:mpZ:mpR:mpT:tranmp:S:T:U:pipPhi:pipTheta:pip2Phi:pip2Theta:pimPhi:pimTheta:neuPhi:neuTheta:mOmega:mNstar:mm2neu:mm2INminusNstar:pipr:pip2r:pimr:neur:dTpip:dTpip2:dTpim:dTneu:neuX:neuY:neuZ:neuE:pipX:pipY:pipZ:pipE:pip2X:pip2Y:pip2Z:pip2E:pimX:pimY:pimZ:pimE:evtX:evtY:evtZ");

    const int narray = 5;
    //add TBER mat branches
    pp->Branch("pip_mat",pip_mat,"pip_mat[5][5]/F");
    pp->Branch("pip2_mat",pip2_mat,"pip2_mat[5][5]/F");
    pp->Branch("pim_mat",pim_mat,"pim_mat[5][5]/F");
    pp->Branch("neu_mat",neu_mat,"neu_mat[5][5]/F");

    //initiate the BOS reading
    initbos();

    for (int i = 1; i < argc; ++i) { //loop over args (input files)
        char *argptr = argv[i];
        if (*argptr != '-') { // check to make sure it's not a flag
            clasEvent event(argptr,&bcs_,1,0); // get first event
            if (event.status()) {
                while (event.read(1) > 0) {//loop over events 
                    clasHEAD_t *HEAD;
                    if (event.type() == 1) {
                        if (HEAD = (clasHEAD_t *)getBank(&bcs_, "HEAD")) { // check for complete HEAD bank
                            ProcessEvent(event); // run analyzer
                        }
                        event.clean();
                    }                    
                }
            }else return 1;
        }
    }
    pp->Write(); //write full TNtuple to file
    f->Close(); // close the file
    return 0;
}
