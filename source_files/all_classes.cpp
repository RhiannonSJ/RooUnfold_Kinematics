#ifndef ALL_CLASSES_CPP
#define ALL_CLASSES_CPP

#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/all_class_headers.h"

namespace xsec{

    // ====================================================================================
    //                                 PARTICLE CLASS
    // ====================================================================================

    // Constructor
    //      
    //      PDG code
    //      Kinetic energy
    //      Opening angle

    Particle::Particle( int PDG, double T, double cos_theta ){
        
        // Define the member variables
        m_PDG              = PDG;
        m_T                = T;
        m_cos              = cos_theta;
        m_is_reconstructed = this->CheckIfReconstructed();
        this->Smear();
    }

    // ====================================================================================
    //                                     Getters
    // ====================================================================================
    
    //  PDG
    int Particle::GetPDG() const{
        return m_PDG;
    }
            
    //  True KE
    double Particle::GetT() const{
        return m_T;
    }

    //  True cos_theta
    double Particle::GetCos() const{
        return m_cos;
    }

    //  Smeared KE
    double Particle::GetTSmeared() const{
        return m_smeared_T;
    }

    //  True cos_theta
    double Particle::GetCosSmeared() const{
        return m_smeared_cos;
    }

    //  Reconstructed?
    bool Particle::GetIfReconstructed() const{
        return m_is_reconstructed;
    }

    // ====================================================================================
    //                                    Operators
    // ====================================================================================
    
    std::ostream& operator<<( std::ostream& os, const Particle& P ){
    
        os << std::setprecision(3) << std::setw(4)
           << "--------------------------------" << std::endl
           << " PDG           : " << P.GetPDG() << std::endl
           << "--------------------------------" << std::endl
           << " True " << std::endl
           << "            KE : " << P.GetT() << std::endl
           << "    cos(theta) : " << P.GetCos() << std::endl
           << "--------------------------------" << std::endl
           << " Smeared " << std::endl
           << "            KE : " << P.GetTSmeared() << std::endl
           << "    cos(theta) : " << P.GetCosSmeared() << std::endl
           << "--------------------------------" << std::endl
           << " Reconstructed?  " << P.GetIfReconstructed() << std::endl
           << "--------------------------------" << std::endl;
           
        return os;
    
    }
    
    // ====================================================================================
    //                              Function definitions
    // ====================================================================================
    
    //  Check if the particle is reconstructed based on cuts and efficiency
    bool Particle::CheckIfReconstructed(){

        if( m_PDG != 13 ){
            return  ( m_T >= 0.05 );
        }

        //CCQE efficiency as function of momentum from Fig10.b of uboone DocDB 7561
        double effmom[] = {0.79,0.95,0.975,0.98,0.98,0.99,0.99,0.995,1.0,1.0};
    
        //CCQE efficiency as function of theta from Fig10.c of uboone DocDB 7561
        double efftheta[] = {1,1,1,1,1,1,0.995,0.985,1,0.985,0.995,1,1,1,1,0.995,0.995,0.995,0.99,0.98,0.97,0.97,0.975,0.98,0.975,0.97,0.95,0.945,0.94,0.92,0.87,0.91};
    
        //Find position of muon momentum in momentum efficiency array
        int posmom = 9; 
    
        for( int i = 0; i < 10; i++ ){
    
            if( m_T < 0.1 * ( i + 1 ) ) {
                posmom = i; 
                break;
            }
        }
    
        //Find position of muon angle in angle efficiency array
        int postheta = 0; 
    
        for( int j = 0; j < 32; j++ ){
    
            if( TMath::ACos( m_cos ) < 0.09817 * ( j + 1 ) ) {
                postheta = j; 
                break;
            }
        }
    
        //Don't have correlations between angle and momentum efficiencies, take the smallest value to be conservative
    
        double efficiency = 1; 
    
        if( effmom[ posmom ] < efftheta[ postheta ] ){
            efficiency = effmom[ posmom ];
        }
        else{
            efficiency = efftheta[ postheta ];
        }

        //Generate random number between 0 and 1 
        // Initialise random number generator
        ROOT::Math::Random<ROOT::Math::GSLRngMT> *_rand = new ROOT::Math::Random<ROOT::Math::GSLRngMT>;
        _rand->SetSeed( time( NULL ) );
       
        double prob = _rand->Uniform();
        
        delete _rand;
    
        return ( prob <= efficiency );

    } // CheckReco

    // Containment
    bool Particle::IsContained(){
        
        // Values from http://pdg.lbl.gov/2012/AtomicNuclearProperties/MUON_ELOSS_TABLES/muonloss_289.png
        double momentum[] = {0.047,0.056,0.068,0.085,0.100,0.153,0.176,0.222,0.287,0.392,0.495,0.900,1.101,1.502,2.103}; //GeV
        double range[] = {0.007,0.013,0.024,0.047,0.076,0.221,0.304,0.482,0.761,1.236,1.708,3.534,4.415,6.126,8.610}; //m
        
        // Find position of muon momentum in momentum array
        int pos = 15;
    
        for( int i = 0; i < 15; i++ ){
    
            if( m_T < momentum[ i ] ){
                pos = i; 
                break;
            }
        }
    
        //Use linear interpolation to calculate expected range in LAr
        double rmu = 0;
    
        if( pos == 15 ){
            rmu = 9;
        }
    
        else if( pos == 0 ){
            rmu = 0.007;
        }
    
        else{
            rmu = range[ pos - 1 ] + ( m_T - momentum[ pos - 1 ] ) * ( range[ pos ] - range[ pos - 1 ] ) / ( momentum[ pos ] - momentum[ pos - 1 ] );
        }
    
        // Initialise random number generator
        ROOT::Math::Random<ROOT::Math::GSLRngMT> *_rand = new ROOT::Math::Random<ROOT::Math::GSLRngMT>;
        _rand->SetSeed( time( NULL ) );

        //Generate a random position inside detector
        double zstart = 5.312 * _rand->Uniform();
        double ystart = 4.294 * _rand->Uniform();
        double xstart = 4.312 * _rand->Uniform();
    
        double phi    = 6.283 * _rand->Uniform();
        double theta  = TMath::ACos( m_cos );
    
        delete _rand;
        
        //Calculate end position
        double zend = zstart + rmu * m_cos;
        double yend = ystart + rmu * TMath::Sin( theta ) * TMath::Sin( phi );
        double xend = xstart + rmu * TMath::Sin( theta ) * TMath::Cos( phi );
    
        if( xend >= 0 && xend <= 4.312 && yend >= 0 && yend <= 4.294 && zend >= 0 && zend <= 5.312 ){
            return 1;
        }
        else{
            return 0;
        }
        
    } // Contained
    
    //  Smearing        
    double Particle::RangeSmear(){
    
        //For contained muons use range based bias and resolution
        //Values from Fig 12 of https://arxiv.org.png/1703.06187.png
        double bias[]       = {-0.0035,-0.0059,-0.0047,-0.0059,-0.0035,-0.0029,-0.0076,-0.0059,0.0006};
        double resolution[] = {0.017,0.021,0.023,0.026,0.025,0.030,0.030,0.040,0.032};

        // Set the mass based on PDG code
        double m;

        if ( m_PDG == 13 || m_PDG == -13 ){
            m = 0.105658;
        }
        else if ( m_PDG == 11 || m_PDG == -11 ){
            m = 0.000511;
        }
        else if( m_PDG == 2212 ){
            m = 0.938272;
        }
        else if( m_PDG == 111 ){
            m = 0.134976;
        } 
        else if( m_PDG == 211 || m_PDG == -211 ){
            m = 0.139570;
        }
        else{
            m = 1;
            std::cout << " Not defined PDG of either muon, proton or pion " << std::endl;
            std::cout << "      PDG : " << m_PDG << std::endl;
            std::cin.get();
        }
        
        double e = m_T + m;

        int pos = 8;

        for( int i = 0; i < 9; i++ ){

            if( e < ( 0.33 + 0.186 * ( i + 1 ) ) ) { 
                pos = i; 
                break;
            }
        }

        // Initialise random number generator
        ROOT::Math::GSLRngMT *_rand = new ROOT::Math::GSLRngMT;
        _rand->Initialize();
        _rand->SetSeed( time( NULL ) );

        double e_smear = e + _rand->Gaussian( resolution[ pos ] * e ) + bias[ pos ] * e;
        double T_smear = e_smear - m;

        delete _rand;
    
        if(T_smear<0) {
            T_smear = 0;
        }
    
        return T_smear;
    } // Range

    double Particle::McsSmear(){

        //For exiting muons use multiple coulomb scattering bias and resolution
        //Values from Fig 5 of https://arxiv.org.png/1703.06187.png
        double bias[]       = {0.0273,0.0409,0.0352,0.0250,0.0227,0.0068,0.0364,0.0273,0.0227};//0.0409
        double resolution[] = {0.127,0.145,0.143,0.141,0.164,0.177,0.250,0.266,0.341};//0.145
    
        int pos = 8;
    
        for( int i = 0; i < 9; i++ ){
    
            if( m_T < ( 0.34 + 0.41 * ( i + 1 ) ) ) {
                pos = i; 
                break;
            }
    
        }
   
        // Initialise random number generator
        ROOT::Math::GSLRngMT *_rand = new ROOT::Math::GSLRngMT;
        _rand->Initialize();
        _rand->SetSeed( time( NULL ) );
    
        double var = _rand->Gaussian( resolution [ pos ] * m_T );
        double T_smear = m_T + var + bias[ pos ] * m_T;

        delete _rand;
    
        if( T_smear < 0 ){
            T_smear = 0; 
        }
            
        return T_smear;
    } // MCS

    void Particle::Smear(){
       
        // Find out if the event was contained
        // If contained: use Range smear
        //               smear cos as normal
        // If exiting  : use MCS smear
        //               smear cos as normal
                    
        // -------------------------------------------------------
        //                   Kinetic energy
        // -------------------------------------------------------
        // Calculate the mean and sigma for the LogNormal function
        //     zeta  = TMath::Log( m * ( 1 / sqrt( 1 + ( var / pow( m, 2 ) ) ) ) );
        //     sigma = sqrt( log( 1 + ( var / pow( m, 2 ) ) ) );
        //     m     = expectation value = Tl
        //     var   = variance = s.d.^2 = ( Tl * 0.1 ) ^ 2
        
        if ( this->IsContained() ){
            m_smeared_T = this->RangeSmear();  
        }
        else{
            m_smeared_T = this->McsSmear();
        }
        // -------------------------------------------------------
        //                  Cos theta
        // -------------------------------------------------------                
        // Calculate the mean and sigma for the LogNormal function
        //      theta = acos(costtheta)
        //      var   = 5 degrees
                 
        ROOT::Math::GSLRngMT *_rand = new ROOT::Math::GSLRngMT;
        _rand->Initialize();
        _rand->SetSeed( time( NULL ) );
    
        double sd_thetamu    = TMath::Pi() / 36; // 5 degrees
        double gaus_theta    = TMath::ACos( m_cos ) + _rand->Gaussian( sd_thetamu );
        m_smeared_cos        = TMath::Cos( gaus_theta );
        
        delete _rand;

    } // Smear
        
    // ====================================================================================
    //                                 EVENT CLASS
    // ====================================================================================

    // Constructor
    Event::Event( std::vector< Particle > parts, bool is_cc, int nu_flav, int physical_proc ){
        m_particle_vect = parts;
        m_is_cc         = is_cc;
        m_neutrino_flav = nu_flav;
        m_physical_proc = physical_proc;
    }

    Event::Event( bool is_cc, int nu_flav, int physical_proc ){ 
        m_is_cc         = is_cc;
        m_neutrino_flav = nu_flav;
        m_physical_proc = physical_proc;
    }

    // ====================================================================================
    //                                     Getters
    // ====================================================================================
    
    // Getters
    //  List of particles
    std::vector< Particle > Event::GetListOfParticles() const{
        return m_particle_vect;
    }
    
    // Get a particle in the list
    Particle Event::GetParticle( int i ) const{
        return m_particle_vect[i];
    }
            
    // Get a particle in the list
    Particle Event::GetMostEnergeticParticleByPDG( int pdg ) const{
        
        // Vector for particles with correct pdg
        std::vector< Particle > correct_pdg;
        std::vector< double > energy;

        // Loop over particles in the event 
        for(  unsigned int i = 0; i < m_particle_vect.size(); ++i ){
        
            // Check if the particle has the correct pdg
            if( m_particle_vect[i].GetPDG() == pdg ){
                correct_pdg.push_back( m_particle_vect[i] );
                energy.push_back( m_particle_vect[i].GetT() );
            }
        }

        // Return the particle with the maximum energy
        int position = std::distance( energy.begin(), std::max_element( energy.begin(), energy.end() ) ); 
        return correct_pdg[ position ];
    }

    // Get particles by the pdg code
    std::vector< Particle > Event::GetParticlesByPDG( int pdg ) const{
            
        std::vector< Particle > particles;

        // Loop over vector of particles and push back
        // any with the chosen pdg code
        for(  unsigned int i = 0; i < m_particle_vect.size(); ++i ){
            
            // Temporary particle
            Particle p = m_particle_vect[i];
            if( p.GetPDG() == pdg ){
                particles.push_back( p );
            }   
        }
    }

    // Size of vector
    int Event::GetLength() const{
        return m_particle_vect.size();
    }

    // CC or NC
    bool Event::GetIsCC() const{
        return m_is_cc;
    }
    
    // Incoming neutrino flavour 
    int Event::GetNeutrinoFlavour() const{
     
        // Possible flavours:
        //  nu_mu    :  14
        //  nu_e     :  12
        //  nu_mu_bar: -14
        //  nu_e_bar : -12
        return m_neutrino_flav;
    }

    // Physical process 
    int Event::GetPhysicalProcess() const{

        // Possible processes:
        //  QE         : 0
        //  MEC        : 1
        //  RES        : 2
        //  DIS        : 3
        //  COH        : 4
        //  No RES, 1pi: 5
        //  Other      : 6
        return m_physical_proc;
    }


    // ====================================================================================
    //                                    Operators
    // ====================================================================================
    
    // Overload []
    Particle Event::operator[]( unsigned int i ) const{

        if ( i >= m_particle_vect.size() ){
            std::cout << "Error: trying to access an element which is out of bounds " << std::endl;
            exit(1);
        }

        return m_particle_vect[i];

    }
    // Overload ostream
      
    std::ostream& operator<<( std::ostream& os, const Event& ev ){
     
        os << "--------------------------------" << std::endl
           << ( ev.GetIsCC() ? " CC " : " NC " ) << std::endl
           << " Physical process : " << ev.GetPhysicalProcess() << std::endl;

        for( int i = 0; i < ev.GetLength(); ++i ){
            os << ev[i];
        }
     
        return os;
    }

    // ====================================================================================
    //                              Function definitions
    // ====================================================================================
    
    // Add particle
    void Event::Add( Particle P ){
        m_particle_vect.push_back( P );
    }
    // Checking functions
    // Does the chosen signal topology match the true
    // topology of the event
    bool Event::CheckIfTrue( Interaction interaction ){ 
    
        // If at the truth level the topology should be cc/nc and isn't, return false
        if( interaction.GetIsCC() != m_is_cc ) return false; 

        top_map topology = interaction.GetMap();

        // Loop over the map
        for( top_map::iterator it = topology.begin(); it != topology.end(); ++it ){

            // Define temporary variables for the current map element
            std::vector< int > PDG_codes = it->first; 
            int n_total                  = it->second;

            // Count the number of particles in the current event with the same PDG codes 
            // as given by the chosen topology
            int counter = 0;
           
            // Loop over particles in current event
            for(  unsigned int i = 0; i < m_particle_vect.size(); ++i ){
            
                // Count the number of particles in event which match the topology
                if( std::find( PDG_codes.begin(), PDG_codes.end(), m_particle_vect[i].GetPDG() ) != PDG_codes.end() ) ++counter;
            }
            if( n_total != -1 ){
                if( counter != n_total ) return false;
            }
            else{
                if( counter == 0 ) return false;
            }
        }
        
        return true; 
    
    }
            
    // Does the chosen signal or background topology match
    // the reconstructed topology of the event
    bool Event::CheckIfReconstructed( Interaction interaction ){

        top_map topology = interaction.GetMap();

        // Loop over the map
        for( top_map::iterator it = topology.begin(); it != topology.end(); ++it ){

            // Define temporary variables for the current map element
            std::vector< int > PDG_codes = it->first; 
            int n_total                  = it->second;

            // Count the number of particles in the current event with the same PDG codes 
            // as given by the chosen topology
            int counter = 0;
            
            // Loop over particles in current event
            for(  unsigned int i = 0; i < m_particle_vect.size(); ++i ){

                // Check if the particle is visible in the detector before adding to the counter
                if( m_particle_vect[i].GetIfReconstructed() ){
                 
                    // Count the number of particles in event which match the topology
                    if( std::find( PDG_codes.begin(), PDG_codes.end(), m_particle_vect[i].GetPDG() ) != PDG_codes.end() ) ++counter;
                }
            }
            if( n_total != -1 ){
                if( counter != n_total ) return false;
            }
            else{
                if( counter == 0 ) return false;
            }
        } 
        
        return true; 
    }
    
    // ====================================================================================
    //                               INTERACTION CLASS
    // ====================================================================================

    // Constructor
    Interaction::Interaction( top_map topology, int primary_pdg, bool is_cc, double probability ){
        m_topology    = topology;
        m_primary_pdg = primary_pdg;
        m_is_cc       = is_cc;
        m_probability = probability;

        bool primary_found = false; 

        for( top_map::iterator it = topology.begin(); it != topology.end(); ++it ){

            // Define temporary variables for the current map element
            std::vector< int > PDG_codes = it->first; 
            int n_total                  = it->second;

            for ( unsigned int i = 0; i < PDG_codes.size(); ++i ){

                if( std::find( PDG_codes.begin(), PDG_codes.end(), primary_pdg ) != PDG_codes.end() ){
                    
                    primary_found = true;
                    
                    if( PDG_codes.size() == 1 ){
                        if( n_total < 1 && n_total != -1 ){
                            std::cerr << " Must have at least one of the primary particle in the chosen signal " <<  std::endl;
                            exit(1);
                        }   
                    }
                    else{
                        std::cerr << " Must only contain 1 PDG option " << std::endl;
                        exit(1);
                    }
                }
            }
        }
        if( !primary_found ){
            std::cerr << " Must at least contain one somewhere " << std::endl;
            exit(1);
         
        }
    }

    // ====================================================================================
    //                                     Getters
    // ====================================================================================
    
    // Getters
    // Get the map
    top_map Interaction::GetMap() const{
        return m_topology;
    }

    // Get the primary PDG code
    int Interaction::GetPrimaryPDG() const{
        return m_primary_pdg;
    }

    // Get CC
    bool Interaction::GetIsCC() const{
        return m_is_cc;
    }

    // Get probability
    double Interaction::GetProbability() const{
        return m_probability;
    }

} // namespace : xsec
    
#endif
