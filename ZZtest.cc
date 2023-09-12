// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {
  /// @brief Normalised cross sections for 5.02 TeV for ZZ diboson production
  class ZZtest : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ZZtest);


    /// Book histograms and initialise projections before the run
    void init() {

      // final state of all stable particles
      Cut particle_cut = (Cuts::abseta < 2.4) and (Cuts::pT > 0.*MeV);
      FinalState fs(particle_cut);

      // select charged leptons
      ChargedLeptons charged_leptons(fs);

      // select final state photons for dressed lepton clustering
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      // select final state prompt charged leptons
      PromptFinalState prompt_leptons(charged_leptons);
      prompt_leptons.acceptMuonDecays(true);
      prompt_leptons.acceptTauDecays(true);

    

      // Dressed leptons from selected prompt charged leptons and photons:
      		// For electrons |eta|< 2.5 and pT > 8 Gev    &   for muons |eta|< 2.4 and pT > 8 Gev  
      Cut lepton_cut = ((Cuts::abseta < 2.4) and 
                        (Cuts::pT > 8.*GeV) and 
                        (((Cuts::abspid == PID::ELECTRON) and (Cuts::abseta < 2.5)) or 
                        (Cuts::abspid == PID::MUON)));

      DressedLeptons dressed_leptons(
        prompt_leptons, 0.1,
        lepton_cut, true
      );
      declare(dressed_leptons, "DressedLeptons");

   

      // book histograms
      book(_hist_ZZ, "d01-x01-y01");
    }


    /// @brief Perform the per-event analysis
    void analyze(const Event& event) {
      vector<DressedLepton> dressedLeptons = apply<DressedLeptons>(
        event,
        "DressedLeptons"
      ).dressedLeptons();

      // Require at least two dressed leptons
      	// Two leptons and two neutrinos or four leptons
      if (dressedLeptons.size() < 2) vetoEvent;

      //  Require a lepton-antilepton pair with opposite sign and same flavour
     
      int c=0;
      for (int iL = 0; iL < int(dressedLeptons.size()); iL++) {
        for (int jL = iL + 1; jL < int(dressedLeptons.size()); jL++) {
          if (abs(dressedLeptons.at(iL).pid() == dressedLeptons.at(jL).pid())) {
          	c += 1;
         }
        }
      } 
     if (dressedLeptons.size() == 4 and (c!=2 or c!=6)) vetoEvent;
     if (dressedLeptons.size() == 2 and c!=1) vetoEvent;
          
      // Require that the leading lepton has at least 20 GeV of pT 
      if (dressedLeptons.at(0).pt() <= 20.*GeV) vetoEvent;
      
      // Require that all identified lepton pairs have at least 4 GeV of invariant mass
      for (int iL = 0; iL < int(dressedLeptons.size()); iL++) {
        for (int jL = iL + 1; jL < int(dressedLeptons.size()); jL++) {
        if (abs(dressedLeptons.at(iL).pid() + dressedLeptons.at(jL).pid()) == 0) {
          if ((dressedLeptons.at(iL).momentum() + 
          dressedLeptons.at(jL).momentum()).mass() <= 4.*GeV) vetoEvent;
        }
        }
      }

      // Require that all identified lepton pairs have at least an invariant mass of:
        // 60 < mll' < 120.
      for (int iL = 0; iL < int(dressedLeptons.size()); iL++) {
        for (int jL = iL + 1; jL < int(dressedLeptons.size()); jL++) {
        if (abs(dressedLeptons.at(iL).pid() + dressedLeptons.at(jL).pid()) == 0) {
          if (((dressedLeptons.at(iL).momentum() + 
          dressedLeptons.at(jL).momentum()).mass() <= 60.*GeV) or ((dressedLeptons.at(iL).momentum() 
          + dressedLeptons.at(jL).momentum()).mass() >= 120.*GeV)) vetoEvent;
        }
        }
      }


      // fill the histograms
      _hist_ZZ->fill(5020);
    }

    /// @brief Normalise histograms after the run
    void finalize() {
      normalize(_hist_ZZ);
    }

  private:
    // Declaration of histograms
    Histo1DPtr _hist_ZZ;
  };

  RIVET_DECLARE_PLUGIN(ZZtest);
}
