/* Event class
 *
 *  Hold all particles in an event
 *
 * =============================================
 *
 *  Author: Rhiannon Jones
 *  Date  : July 2017
 *
 * =============================================
 */


#ifndef EVENT_CLASS_H
#define EVENT_CLASS_H

// Header file for all std::library and ROOT headers needed
#include "/hepstore/rjones/Exercises/Unfolding_General/headers/roo_unfold.h"
#include "/hepstore/rjones/Exercises/Unfolding_General/headers/particle_class.h"

namespace xsec{

    class Event{
    
        public: 

            // Constructor
            //      vector< Particle > 
            Event( std::vector< Particle > parts );

            // Destructor
            ~Event(){}

            // Getters
            //  List of particles
            std::vector< Particle > GetListOfParticles() const;
            
            //  Length of list
            int GetLength() const;
            
            // Overload []
            Particle operator[]( int i ) const; 

            // Overload ostream
            friend std::ostream& operator<<( std::ostream& , const Event& );
            
            // Functions
            //  Add particle
            void Add( Particle &P );


        private:

            // Member variables
            std::vector< Particle > event;
            
    }; // class: Event

} // namespace : xsec


#endif // EVENT_CLASS_H
