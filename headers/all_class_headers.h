#ifndef ALL_CLASS_HEADERS_H
#define ALL_CLASS_HEADERS_H

#include "/hepstore/rjones/Exercises/RooUnfold_Kinematics/headers/roo_unfold.h"

namespace xsec{

    // Typedef for the map
    typedef std::map< std::vector< int >, int > top_map; 

    class Particle {
    
        public: 

            // Constructor
            //      
            //      PDG code
            //      Kinetic energy
            //      Opening angle
            //      Smeared kinetic energy
            //      Smeared opening angle
            //      Is visible
            //      Is reconstructed
            //
            Particle( int PDG, double T, double cos_theta );

            // Destructor
            ~Particle(){}

            // Getters
            //  PDG
            int GetPDG() const;
            
            //  True KE
            double GetT() const;

            //  True cos_theta
            double GetCos() const;

            //  Smeared KE
            double GetTSmeared() const;

            //  Smeared cos_theta
            double GetCosSmeared() const;
            
            //  Reconstructed?
            bool GetIfReconstructed() const;
        
            friend std::ostream& operator<<( std::ostream& , const Particle& );
        
        private:

            // Member variables
            int    m_PDG;
            double m_T;
            double m_cos;
            double m_smeared_T;
            double m_smeared_cos;
            bool   m_is_reconstructed;

            // Functions
            //  Check if the particle is reconstructed based on cuts and efficiency
            bool CheckIfReconstructed();

            // Check if the particle is contained within the TPC
            bool IsContained();
            
            //  Smearing
            double RangeSmear();
            double McsSmear();
            void Smear();


    }; // class: Particle

    class Interaction{
    
        public : 

            // Constructor
            Interaction( top_map topology, int primary_pdg, bool is_cc, double probability = 1. );

            // Getters
            // Get the map
            top_map GetMap() const;

            // Get the primary PDG code
            int GetPrimaryPDG() const;

            // Get CC
            bool GetIsCC() const;

            // Get probability
            double GetProbability() const;

        private :

            // Member variables
            top_map m_topology;
            double m_probability;
            int m_primary_pdg;
            bool m_is_cc;

    }; // Interaction class
    
    class Event{
    
        public: 

            // Constructor
            Event( std::vector< Particle > parts, bool is_cc, int nu_flav, int physical_proc );

            // Default constructor
            //  nu_flav:       muon
            //  physical_proc: QE
            Event( bool is_cc = true, int nu_flav = 14, int physical_proc = 0 );

            // Destructor
            ~Event(){}

            // Getters
            //  List of particles
            std::vector< Particle > GetListOfParticles() const;
           
            // Get a particle in the list
            Particle GetParticle( int i ) const;
            
            // Get a particle in the list
            Particle GetMostEnergeticParticleByPDG( int pdg ) const;
    
            // Get particles by the pdg code
            std::vector< Particle > GetParticlesByPDG( int pdg ) const;

            //  Length of list
            int GetLength() const;
            
            // Get Is CC
            bool GetIsCC() const;
            
            // Get physical process
            int GetPhysicalProcess() const;
            
            // Get physical process
            int GetNeutrinoFlavour() const;
            
            // Functions
            //  Add particle
            void Add( Particle P );

            // Checking functions
            // Does the chosen signal topology match the true
            // topology of the event
            bool CheckIfTrue( Interaction interaction );
            
            // Does the chosen signal or background topology match
            // the reconstructed topology of the event
            bool CheckIfReconstructed( Interaction interaction );

            // Overload []
            Particle operator[]( unsigned int i ) const; 

            // Overload ostream
            friend std::ostream& operator<<( std::ostream& , const Event& );

        private:

            // Member variables
            std::vector< Particle > m_particle_vect;
            bool                    m_is_cc;
            int                     m_physical_proc;
            int                     m_neutrino_flav;
            
    }; // class: Event

} // namespace : xsec

#endif
