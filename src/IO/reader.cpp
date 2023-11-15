/*
 * reader.cpp
 *
 *  Created on: Feb 7, 2023
 *      Author: pemueller
 */

#include "reader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <boost/make_shared.hpp>


Molecule readMolecule( string filename )
{
#ifdef DEBUG
	cout << __FUNCTION__ << endl;
#endif
	ifstream infile(filename);
	string line, inString;
	stringstream ss;
	vector<AtomPtr> atoms;

	getline( infile, line );

	ss << line;
	ss >> inString;

	const int nAtoms = stoi(inString);

	getline( infile, line );

	for ( int index = 0; index < nAtoms; ++index )
	{
		ss.clear();
		ss.str("");

		getline( infile, line );

		ss << line;

		string elementSymbol;

		ss >> elementSymbol;

		Vector3d coordinates;
		for ( size_t j = 0; j < 3; ++j )
		{
			ss >> coordinates(j);
		}

		atoms.push_back(boost::make_shared<Atom>(index, elementSymbol, coordinates));
	}

	return Molecule(0, atoms);
}
