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


#ifndef EVENT_CLASS_CPP
#define EVENT_CLASS_CPP

// Header file for all std::library and ROOT headers needed
#include "/hepstore/rjones/Exercises/Unfolding_General/headers/event_class.h"
#include "/hepstore/rjones/Exercises/Unfolding_General/headers/particle_class.h"

namespace xsec{

        // Constructor
        Event::Event( std::vector< Particle > parts ){
            event = parts;
        }

        // Getters
        //  List of particles
        std::vector< Particle > Event::GetListOfParticles() const{
            return event;
        }

        // Size of vector
        int Event::GetLength() const{
            return event.size();
        }

        // Overload []
        Particle Event::operator[]( int i ) const{

            if ( i < 0 || i >= event.size() ){
                std::cout << "Error: trying to access an element which is out of bounds " << std::endl;
                exit(1);
            }

            return event[i];

        }

        // Overload ostream
        
        std::ostream& operator<<( std::ostream& os, const Event& ev ){
     
            for( int i = 0; i < ev.GetLength(); ++i ){
                os << ev[i] << std::endl; 
            }
            return os;

        }

        // Functions
        // Add particle
        void Event::Add( Particle &P ){
            event.push_back( P );
        }

} // namespace : xsec


#endif // EVENT_CLASS_H
