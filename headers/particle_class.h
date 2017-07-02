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


#ifndef PARTICLE_CLASS_H
#define PARTICLE_CLASS_H

// Header file for all std::library and ROOT headers needed
#include "/hepstore/rjones/Exercises/Unfolding_General/headers/roo_unfold.h"

namespace xsec{

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

} // namespace : xsec


#endif // PARTICLE_CLASS_H
