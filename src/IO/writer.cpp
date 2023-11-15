/*
 * writer.cpp
 *
 *  Created on: May 30, 2023
 *      Author: pemueller
 */

#include "writer.h"
#include <iostream>
#include <cmath>
#include "../testing/pgTests.h"

string getNumberAsString( complex<double> input )
{
	string real = to_string(input.real());
	string imag = to_string(input.imag());

    // Remove trailing zeroes
	real = real.substr(0, real.find_last_not_of('0')+1);
    // If the decimal point is now the last character, remove that as well
    if(real.find('.') == real.size()-1)
    {
    	real = real.substr(0, real.size()-1);
    }

    // Remove trailing zeroes
	imag = imag.substr(0, imag.find_last_not_of('0')+1);
    // If the decimal point is now the last character, remove that as well
    if(imag.find('.') == imag.size()-1)
    {
    	imag = imag.substr(0, imag.size()-1);
    }

    string result = "";

    if ( abs(input.real()) > 0.00001 )
    {
    	result += real;
    	if ( input.imag() > 0.00001 )
    		result += '+';
    }

	if ( abs(input.imag()) > 0.00001 )
	{
		if ( imag != "1" && imag != "-1")
			result += imag;
		if ( imag == "-1" )
			result += '-';
		result += 'i';
	}

//	cout << "\"" << result << "\"" << endl;
//	cout << result.size() << endl;

    // if result is empty, it shall be zero
	if (result.size() == 0)
		return "0";

    return result;
}

string getIrrepVectorAsString( VectorXcd irrep )
{
	string result="";

	for ( size_t i = 0; i < irrep.size(); ++i )
	{
			result += getNumberAsString(irrep(i)) + ' ';
	}
	result = result.substr(0,result.size()-1);// remove last space character

//	cout << irrep.transpose() << " -> " << result << endl;

	return result;
}


