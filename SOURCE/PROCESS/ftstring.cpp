#include "ftstring.h"

template<typename T>
string to_string( const T & Value )
{
    // utiliser un flux de sortie pour créer la chaîne
    ostringstream oss;
    // écrire la valeur dans le flux
    oss << Value;
    // renvoyer une string
    return oss.str();
}

template <class T>
void from_string(T& t, 
                 const string& s, 
                 ios_base& (*f)(ios_base&))
{
  istringstream iss(s);
  iss >> f >> t;
}
