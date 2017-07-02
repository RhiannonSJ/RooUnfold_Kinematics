/* Particle class
 *
 *  Hold the information for a single particle
 *      - PDG code
 *      - Kinetic energy
 *      - Opening angle
 *      - Smeared kinetic energy
 *      - Smeared opening angle
 *      - Is Visible?
 *      - Is Reconstructed?
 *
 * =============================================
 *
 *  Author: Rhiannon Jones
 *  Date  : July 2017
 *
 * =============================================
 */


#ifndef PARTICLE_CLASS_CPP
#define PARTICLE_CLASS_CPP

// Header file for all std::library and ROOT headers needed
#include "/hepstore/rjones/Exercises/Unfolding_General/headers/particle_class.h"

namespace xsec{


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

        if ( m_PDG == 13 ){
            m = 0.105658;
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
    
} // namespace : xsec
    
    
#endif // PARTICLE_CLASS_H
