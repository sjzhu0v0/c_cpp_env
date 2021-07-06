#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TObject.h"
#include "TVector3.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TRandom3.h"  

#include "Pythia8/Pythia.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
using namespace Pythia8;

void fillParticle(int id, double px, double py, double pz,
  Event& event, ParticleData& pdt, Rndm& rndm, bool atRest = false,
  bool hasLifetime = false) {

  // Reset event record to allow for new event.
  event.reset();

  // Select particle mass; where relevant according to Breit-Wigner.
  double mm = pdt.mSel(id);

  double ee;
  // Special case when particle is supposed to be at rest.
  if (atRest) {
    ee = mm;
  }
  ee=sqrt(mm*mm+px*px+py*py+pz*pz);
  // Angles as input or uniform in solid angle.
  

  // Store the particle in the event record.
  int iNew = event.append( id, 1, 0, 0, px,
    py, pz, ee, mm);

  // Generate lifetime, to give decay away from primary vertex.
  if (hasLifetime) event[iNew].tau( event[iNew].tau0() * rndm.exp() );

}

int main()
{
    TFile *f1= new TFile("bmeson.root");
    TFile *f = new TFile("MC.root", "recreate");

    int id = 0;
    double px = 0.;
    double py = 0.;
    double pz = 0.;
    double e = 0.;
    double vx = 0.;
    double vy = 0.;
    double vz = 0.;
    int trackMC = 0;
    
    const string xmlDB = "/home/ephy/work/sjzhu/PYTHIA8303/share/Pythia8/xmldoc";
    Pythia pythia(xmlDB.c_str());
    Event &event = pythia.event;
    Settings& settings = pythia.settings;
    ParticleData& pdt = pythia.particleData;
    
    pythia.readFile("EP_NCR.cmnd");//common switchs
    pythia.readString("ProcessLevel:all = off");/*key switchs in examples/main21.cc*/
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("521:onMode = off");
    pythia.readString("521:onIfMatch = -421 211");
    pythia.readString("421:onMode = off");
    pythia.readString("421:onIfMatch = -321 211");
    pythia.init();

    TTree *t_MC = new TTree("t_MC", "Ep_out");
    TTree *t_track=new TTree("Track","Track");
    TTree *data = (TTree*) f1->Get("t1");
    
    t_MC->Branch("id", &id, "id/I");
    t_MC->Branch("px", &px, "px/D");
    t_MC->Branch("py", &py, "py/D");
    t_MC->Branch("pz", &pz, "pz/D");
    t_MC->Branch("e", &e, "e/D");
    t_MC->Branch("vx", &vx, "vx/D");
    t_MC->Branch("vy", &vy, "vy/D");
    t_MC->Branch("vz", &vz, "vz/D");
    t_track->Branch("trackMC",&trackMC,"trackMC/I");

    data->SetBranchAddress("px",&px);
    data->SetBranchAddress("py",&py);
    data->SetBranchAddress("pz",&pz);
    int num_data=data->GetEntries();

    int  maxNumberOfEvents = settings.mode("Main:numberOfEvents");
    int count = 0;int totalrun=0;bool exist=false;

    TRandom3 rndgen;

    while(count < maxNumberOfEvents)
    {   
        totalrun++;
        data->GetEntry(rndgen.Uniform(0,num_data));
        fillParticle(521,px,py,pz,event,pdt,pythia.rndm,false,true);
        if (pythia.next())
        {
            trackMC = 0;
            count++;
            exist=false;
            for (int i0 = 0; i0 < event.size(); i0++)
            {
                if (event[i0].status() > 0)
                {
                    id = event[i0].id();
                    px = event[i0].px();
                    py = event[i0].py();
                    pz = event[i0].pz();
                    e  = event[i0].e();
                    vx = event[i0].xProd();
                    vy = event[i0].yProd();
                    vz = event[i0].zProd();
                    trackMC++;t_MC->Fill();
                }
            }
            //cout << t_MC->GetEntries() << endl;
            t_track->Fill();
        }
    }

    f->cd();
    t_MC->Write();
    t_track->Write();
    f->Close();
    std::cout << count << "     "<< totalrun<< std::endl;
    std::cout << "success";

    return 0;
}
