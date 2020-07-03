#ifndef __FTSTRING__
#define __FTSTRING__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

template<typename T>
string to_string( const T & Value );

template <class T>
void from_string(T& t, 
                 const string& s, 
                 ios_base& (*f)(ios_base&));
                 
               
#endif
