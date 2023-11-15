//============================================================================
// Name        : symmetrizer.cpp
// Author      : me
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>

#include <boost/make_shared.hpp>

#include "characterTable/symmetryElements.h"
#include "common/molecule.h"
#include "IO/reader.h"
#include "characterTable/pointGroup.h"
#include "common/constants.h"
#include "IO/writer.h"

#include "testing/pgTests.h"

using namespace std;

void testSymmetryOperations()
{
if (debugLevel > 3) {
	cout << __FUNCTION__ << endl;
}
	Vector3d coordinates (2.0,0.0,2.0);

	cout << "original point: " << coordinates.transpose() << endl;

	SymmetryOperation inversion( "inversion", Vector3d::Zero() );

	Vector3d inversionResult = inversion.applySymmetryOperation(coordinates);

	cout << "result of the inversion: " << inversionResult.transpose() << endl;

	SymmetryOperation rotation( "rotation", Vector3d::Zero(), Vector3d(1,0,0), 3 );
	Vector3d rotationResult1 = rotation.applySymmetryOperation(coordinates);
	Vector3d rotationResult2 = rotation.applySymmetryOperation(rotationResult1);
	Vector3d rotationResult3 = rotation.applySymmetryOperation(rotationResult2);
	cout << "result of the first rotation: " << rotationResult1.transpose() << endl;
	cout << "result of the second rotation: " << rotationResult2.transpose() << endl;
	cout << "result of the third rotation: " << rotationResult3.transpose() << endl;

	SymmetryOperation reflection( "reflection", Vector3d::Zero(), Vector3d(1,0,0));
	Vector3d reflectionResult1 = reflection.applySymmetryOperation(coordinates);
	Vector3d reflectionResult2 = reflection.applySymmetryOperation(reflectionResult1);
	cout << "result of the first reflection: " << reflectionResult1.transpose() << endl;
	cout << "result of the second reflection: " << reflectionResult2.transpose() << endl;

	SymmetryOperation rotationReflection( "rotationReflection", Vector3d::Zero(), Vector3d(1,0,0), 3 );
	Vector3d rotationReflection1 = rotationReflection.applySymmetryOperation(coordinates);
	Vector3d rotationReflection2 = rotationReflection.applySymmetryOperation(rotationReflection1);
	Vector3d rotationReflection3 = rotationReflection.applySymmetryOperation(rotationReflection2);
	Vector3d rotationReflection4 = rotationReflection.applySymmetryOperation(rotationReflection3);
	Vector3d rotationReflection5 = rotationReflection.applySymmetryOperation(rotationReflection4);
	Vector3d rotationReflection6 = rotationReflection.applySymmetryOperation(rotationReflection5);
	cout << "result of the first rotation reflection: " << rotationReflection1.transpose() << endl;
	cout << "result of the second rotation reflection: " << rotationReflection2.transpose() << endl;
	cout << "result of the third rotation reflection: " << rotationReflection3.transpose() << endl;
	cout << "result of the fourth rotation reflection: " << rotationReflection4.transpose() << endl;
	cout << "result of the fifth rotation reflection: " << rotationReflection5.transpose() << endl;
	cout << "result of the sixth rotation reflection: " << rotationReflection6.transpose() << endl;

}

void testDistance( MolPtr molecule )
{
	vector<SymmOpPtr> sos = molecule->getPointGroup()->getSymmetryOperations();

	vector<AtomPtr> atoms = molecule->getAtoms();

	for ( auto so : sos )
	{
		if ( !(so->getType() == "rotation" || so->getType() == "reflection") )
			continue;

		cout << "Distance between symmetry operation " << so->getSymbol() << " at " << so->getAxis().transpose() << endl;
		for ( auto atom: atoms )
		{
			cout << "atom " << atom->getElementSymbolForPrinting() << " " << molecule->getPointGroup()->getDistanceFromSO(atom, so) << endl;
		}
	}

	for ( auto so : sos )
	{
		if ( !(so->getType() == "rotation" || so->getType() == "reflection") )
			continue;

		cout << "Weigths of atoms for symmetry operation " << so->getSymbol() << " at " << so->getAxis().transpose() << endl;

		for ( auto atom: atoms )
		{
			cout << "atom " << atom->getElementSymbolForPrinting() << " " << molecule->getPointGroup()->weightWithAtom( molecule->getPointGroup()->getDistanceFromSO(atom, so) , atom ) << endl;
		}
	}

}

void writePointGroup( MolPtr molecule, string filenamePrev )
{
if (debugLevel > 3) {
	cout << __FUNCTION__ << endl;
}
	string filename = filenamePrev.substr(0, filenamePrev.length() - 4) + ".txt";

	cout << "writing " << filename << endl;
	ofstream outfile(filename);

	outfile << "Molecule "  << molecule->getName() << " contains " << molecule->getAtoms().size() << " atoms" << endl;
	outfile << "Found point group " << molecule->getPointGroup()->getName() << endl;
	for ( auto so : molecule->getPointGroup()->getSymmetryOperations() )
	{
		outfile << "Contains " << so->getSymbol() << /*so.getType() << " of order " << so.getOrderOfRotation() <<*/ " around " << so->getCenter().transpose() << " " << so->getAxis().transpose() << " with error " << so->getError() << endl;
	}

	for ( int i = 0; i < molecule->getPointGroup()->getSymmetryOperations().size(); ++i )
	{
		for ( int j = 0; j < molecule->getPointGroup()->getSymmetryOperations().size(); ++j )
			outfile << molecule->getPointGroup()->getSymmetryElementMatrix()[i][j] << " ";
		outfile << endl;
	}

	for ( int j = 0; j < molecule->getPointGroup()->getCharacterTable()->getClasses().size(); ++j )
	{
		vector<SymmOpPtr> c = molecule->getPointGroup()->getCharacterTable()->getClasses()[j];

		outfile << "Class " << j + 1;

		for ( int i = 0; i < c.size(); ++i )
		{
			outfile << " " << c[i]->getSymbol();
		}

		outfile << endl;
	}

	outfile << molecule->getPointGroup()->getMainRotationAxes().size() << " main rotation axes in total" << endl;
	for ( auto mra : molecule->getPointGroup()->getMainRotationAxes() )
	{
		outfile << "Main rotation axis " << mra->getSymbol() << " around " << mra->getAxis().transpose() << endl;
	}

	outfile << endl << "character table" << endl;
	for ( auto so : molecule->getPointGroup()->getCharacterTable()->getClasses() )
	{
		if ( so.size() > 1 )
			outfile << so.size();
		outfile << so[0]->getSymbol() << " ";
	}
	outfile << endl;

	for ( size_t i = 0; i <  molecule->getPointGroup()->getCharacterTable()->getIrreps().size(); ++i )
	{
		auto irrep = molecule->getPointGroup()->getCharacterTable()->getIrreps()[i];
		outfile << irrep->getSymbol() << " ";
		for ( auto characters : irrep->getCharacters() )
		{
			outfile << getIrrepVectorAsString(characters) << endl;
		}
		outfile << endl;
	}

	outfile.close();
}

void setDebugLevel( int argc, char** argv )
{
	int pos = -1;
	for ( int i = 1; i < argc; ++i )
	{

		string levelAsString = argv[i];
		if ( levelAsString == "-d" )
		{
			pos = i;
			break;
		}
	}

	try {
		debugLevel = stoi(argv[pos+1]);

		cout << "Debugging verbosity set to " << debugLevel << endl;
	}
	catch (...)
	{ debugLevel = 0; };
}

void setPointGroup( int argc, char** argv )
{
	int pos = -1;
	for ( int i = 1; i < argc; ++i )
	{
		string levelAsString = argv[i];
		if ( levelAsString == "-pg" )
		{
			pos = i;
			break;
		}
	}

	if (pos == -1)
		return;

	try {
		suggestedPointGroup = argv[pos+1];

		if ( !isPointGroupValid(suggestedPointGroup) )
		{
			cout << "Manually selected point group " << suggestedPointGroup << " is invalid!" << endl;
			cout << "Falling back to automatic determination" << endl;
			suggestedPointGroup = "";
		}
		else
		{
			cout << "Manually selected point group " << suggestedPointGroup << " of type " << suggestedPointGroupType << endl;
			cout << "I will do my best" << endl;
		}
	}
	catch (...)
	{ suggestedPointGroup = ""; };
}

int main( int argc, char** argv )
{
if (debugLevel > 3) {
	cout << __FUNCTION__ << endl;
}
	// read molecule
	const string filename = argv[1];


	setDebugLevel(argc, argv);

	cout << endl << "Reading " << filename << endl;

	MolPtr mol = boost::make_shared<Molecule>(readMolecule(filename));

	setPointGroup( argc, argv );

//	cout << "Molecule "  << mol->getName() << " contains " << mol->getAtoms().size() << " atoms" << endl;

	mol->findPointGroup();

	//testDistance(mol);

	writePointGroup(mol, filename);
//	cout << "Found point group " << mol->getPointGroup().getName() << " for molecule " << mol->getName() << endl;

/*	for ( auto so : mol->getPointGroup().getSymmetryOperations() )
		cout << "Contains " << so.getType() << " of order " << so.getOrderOfRotation() << " around " << so.getCenter().transpose() << " " << so.getAxis().transpose() << endl;*/

	return 0;
}
