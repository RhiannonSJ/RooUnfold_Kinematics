#ifndef UNFOLDING_FUNCTIONS_CPP
#define UNFOLDING_FUNCTIONS_CPP

#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/unfolding_functions.h"
using namespace xsec;

// Load the list of events from the gst tree and push onto a vector of events
void LoadEventList( TTree *gst, std::vector< Event > &event_list ){

    // Get the branches from the tree
    //  neu
    //  cc
    //  nc
    //  nf
    //  pdgf
    //  fspl
    //  cthl
    //  cthf
    //  ef
    //  el
    //  qel
    //  mec
    //  res
    //  dis
    //  coh
    //  npip
    //  npim
    //  npi0

    TBranch *b_neu  = gst->GetBranch("neu");
    TBranch *b_cc   = gst->GetBranch("cc");
    TBranch *b_nc   = gst->GetBranch("nc");
    TBranch *b_nf   = gst->GetBranch("nf");
    TBranch *b_pdgf = gst->GetBranch("pdgf");
    TBranch *b_fspl = gst->GetBranch("fspl");
    TBranch *b_cthl = gst->GetBranch("cthl");
    TBranch *b_cthf = gst->GetBranch("cthf");
    TBranch *b_ef   = gst->GetBranch("Ef");
    TBranch *b_el   = gst->GetBranch("El");
    TBranch *b_qel  = gst->GetBranch("qel");
    TBranch *b_mec  = gst->GetBranch("mec");
    TBranch *b_res  = gst->GetBranch("res");
    TBranch *b_dis  = gst->GetBranch("dis");
    TBranch *b_coh  = gst->GetBranch("coh");
    TBranch *b_npip = gst->GetBranch("nfpip");
    TBranch *b_npim = gst->GetBranch("nfpim");
    TBranch *b_npi0 = gst->GetBranch("nfpi0");

    int n_entries = gst->GetEntries();

    // Loop over event list in tree and fill all components of the event needed
    for( int i = 0; i < n_entries; ++i ){
        
        // Get entry in tree
        gst->GetEntry(i);

        // Get the integers from branches
        int neu     = b_neu->GetLeaf("neu")->GetValue();
        bool cc     = b_cc->GetLeaf("cc")->GetValue();
        int nf      = b_nf->GetLeaf("nf")->GetValue();
        int fspl    = b_fspl->GetLeaf("fspl")->GetValue();
        double cthl = b_cthl->GetLeaf("cthl")->GetValue();
        double el   = b_el->GetLeaf("El")->GetValue();
        bool qel    = b_qel->GetLeaf("qel")->GetValue();
        bool mec    = b_mec->GetLeaf("mec")->GetValue();
        bool res    = b_res->GetLeaf("res")->GetValue();
        bool dis    = b_dis->GetLeaf("dis")->GetValue();
        bool coh    = b_coh->GetLeaf("coh")->GetValue();
        int npip    = b_npip->GetLeaf("nfpip")->GetValue();
        int npim    = b_npim->GetLeaf("nfpim")->GetValue();
        int npi0    = b_npi0->GetLeaf("nfpi0")->GetValue();

        // Find the physical process of the event and enumerate according to 
        //  QE           : 0
        //  MEC          : 1
        //  RES          : 2
        //  DIS          : 3
        //  COH          : 4
        //  Non RES, 1pi : 5
        //  Other        : 6

        int physical_process; 
        
        if ( qel ){
            physical_process = 0;    
        }
        else if ( mec ){
            physical_process = 1;    
        }
        else if ( res ){
            physical_process = 2;    
        }
        else if ( dis ){
            physical_process = 3; 
        }
        else if ( coh ){
            physical_process = 4; 
        }
        else if ( !res && ( npip + npim + npi0 == 1 ) ){
            physical_process = 5; 
        }
        else{
            physical_process = 6;
        }

        // Initialise new Event object with cc, the neutrino flavour, physical process 
        Event event( cc, neu, physical_process ); 
        
        double mass;
        // Add the final state primary lepton to the event
        if( fspl == 13 || fspl == -13 ){
            mass = 0.105658;
            event.Add( Particle( fspl, el - mass, cthl ) );
        }
        else if( fspl == 11 || fspl == -11 ){ 
             mass = 0.000511;
            event.Add( Particle( fspl, el - mass, cthl ) );
        }

        // Loop over final state particles and add the protons and pions to the event
        for( int j = 0; j < nf; ++j ){

            // Get the energy and angle doubles and pdgf
            int pdgf    = b_pdgf->GetLeaf("pdgf")->GetValue(j);
            double ef   = b_ef->GetLeaf("Ef")->GetValue(j); 
            double cthf = b_cthf->GetLeaf("cthf")->GetValue(j);

            // For each particle type, get the mass and add the particle to the event
            if( pdgf == 2212 ){
                mass = 0.938272;
                event.Add( Particle( pdgf, ef - mass, cthf ) ); 
            }
            else if( pdgf == 111 ){
                mass = 0.134976;
                event.Add( Particle( pdgf, ef - mass, cthf ) ); 
            } 
            else if( pdgf == 211 || pdgf == -211 ){
                mass = 0.139570;
                event.Add( Particle( pdgf, ef - mass, cthf ) ); 
            }
        }
        // Add current event to event list
        event_list.push_back(event);
    }
}

void GetResponse( const std::vector< Event >       & event_list,
                  const Interaction                & interaction,
                  const std::vector< Interaction > & background,
                  RooUnfoldResponse                & response ) {
    
    for( unsigned int i = 0; i < event_list.size(); ++i ) {
        
        // Temporary event for ease
        Event ev = event_list[i];

        if( ev.CheckIfTrue( interaction ) ) {
       
            // Get the primary PDG code
            Particle primary = ev.GetMostEnergeticParticleByPDG( interaction.GetPrimaryPDG() );
                    
            if( ev.CheckIfReconstructed( interaction ) ){
                response.Fill( primary.GetCosSmeared(), primary.GetTSmeared(), primary.GetCos(), primary.GetT() );
            }
            else{
                response.Miss( primary.GetCos(), primary.GetT() ); 
            }
        }
        else{
            if( ev.CheckIfReconstructed( interaction ) ){
                Particle reco_primary = ev.GetMostEnergeticParticleByPDG( interaction.GetPrimaryPDG() );
                response.Fake( reco_primary.GetCosSmeared(), reco_primary.GetTSmeared() );
            }
            for( unsigned int j = 0; j < background.size(); ++j ){
                if( ev.CheckIfReconstructed( background[j] ) ){ 

                    Particle bg_primary = ev.GetMostEnergeticParticleByPDG( background[j].GetPrimaryPDG() );

                    // Generate a random probability to see if the backgrounds should be counted
                    ROOT::Math::Random<ROOT::Math::GSLRngMT> *_rand = new ROOT::Math::Random<ROOT::Math::GSLRngMT>;
                    _rand->SetSeed( time( NULL ) );
                    double prob = _rand->Uniform();
   
                    // If the probability given is less than the random probability generated
                    //  count the background
                    if( prob <= background[j].GetProbability() ){
                        response.Fake( bg_primary.GetCosSmeared(), bg_primary.GetTSmeared() );
                    }
                }
            }
        }
    }
} 

void GetTrueRecoHists( const std::vector< Event >       & event_list,
                       const Interaction                & interaction,
                       const std::vector< Interaction > & background,
                       TH2D *true_hist,
                       TH2D *reco_hist ) {
    
    for( unsigned  int i = 0; i < event_list.size(); ++i ) {
        
        // Temporary event for ease
        Event ev = event_list[i];

        if( ev.CheckIfTrue( interaction ) ) {
       
            // Get the primary PDG code
            Particle primary = ev.GetMostEnergeticParticleByPDG( interaction.GetPrimaryPDG() );
            
            true_hist->Fill( primary.GetCos(), primary.GetT() );
                    
            if( ev.CheckIfReconstructed( interaction ) ){
                reco_hist->Fill( primary.GetCosSmeared(), primary.GetTSmeared() );
            }
        }
        else{
            if( ev.CheckIfReconstructed( interaction ) ){
                Particle reco_primary = ev.GetMostEnergeticParticleByPDG( interaction.GetPrimaryPDG() );
                reco_hist->Fill( reco_primary.GetCosSmeared(), reco_primary.GetTSmeared() );
            }
            for( unsigned int j = 0; j < background.size(); ++j ){
                if( ev.CheckIfReconstructed( background[j] ) ){ 
                    
                    Particle bg_primary = ev.GetMostEnergeticParticleByPDG( background[j].GetPrimaryPDG() );
                    
                    // Generate a random probability to see if the backgrounds should be counted
                    ROOT::Math::Random<ROOT::Math::GSLRngMT> *_rand = new ROOT::Math::Random<ROOT::Math::GSLRngMT>;
                    _rand->SetSeed( time( NULL ) );
                    double prob = _rand->Uniform();
    
                    // If the probability given is less than the random probability generated
                    //  count the background
                    if( prob <= background[j].GetProbability() ){
                        reco_hist->Fill( bg_primary.GetCosSmeared(), bg_primary.GetTSmeared() );
                    }
                }
            }
        }
    }
} 

void GetRecoEventList( const std::vector< Event >       & event_list,
                       const Interaction                & interaction,
                       const std::vector< Interaction > & background,
                       std::vector< Particle >          & primary_list,
                       std::vector< Event >             & reco_event_list ) {
    
    for( unsigned  int i = 0; i < event_list.size(); ++i ) {
        
        // Temporary event for ease
        Event ev = event_list[i];

        if( ev.CheckIfTrue( interaction ) ) {
       
            if( ev.CheckIfReconstructed( interaction ) ){
                
                primary_list.push_back(  ev.GetMostEnergeticParticleByPDG( interaction.GetPrimaryPDG() ) );
                
                reco_event_list.push_back( ev );
            }
        }
        else{
            if( ev.CheckIfReconstructed( interaction ) ){
             
                primary_list.push_back(  ev.GetMostEnergeticParticleByPDG( interaction.GetPrimaryPDG() ) );
                
                reco_event_list.push_back( ev );
            }
            for( unsigned int j = 0; j < background.size(); ++j ){
                if( ev.CheckIfReconstructed( background[j] ) ){ 
                    
                    // Generate a random probability to see if the backgrounds should be counted
                    ROOT::Math::Random<ROOT::Math::GSLRngMT> *_rand = new ROOT::Math::Random<ROOT::Math::GSLRngMT>;
                    _rand->SetSeed( time( NULL ) );
                    double prob = _rand->Uniform();
    
                    // If the probability given is less than the random probability generated
                    //  count the background
                    if( prob <= background[j].GetProbability() ){
                        
                        primary_list.push_back(  ev.GetMostEnergeticParticleByPDG( background[j].GetPrimaryPDG() ) );
             
                        reco_event_list.push_back( ev );
                    }
                }
            }
        }
    }
} 
void Set2DHistInfo( TH2D *hist, const char x_axis[1024], const char y_axis[1024], const char title[1024], const char draw_opt[1024] ){


    gStyle->SetPalette(55);
    gStyle->SetNumberContours(250);

    hist->SetStats(kFALSE);
    hist->GetXaxis()->SetTitle( x_axis );
    hist->GetYaxis()->SetTitle( y_axis );
    hist->SetTitle( title );
    hist->Draw( draw_opt );

}
#endif
