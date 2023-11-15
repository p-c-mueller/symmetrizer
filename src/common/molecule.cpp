/*
 * molecule.cpp
 *
 *  Created on: Feb 6, 2023
 *      Author: pemueller
 */

#include "molecule.h"
#include "math.h"
#include "constants.h"
#include <iostream>
#include <Eigen/Geometry>
#include <map>
#include <algorithm>
#include <boost/make_shared.hpp>
#include <boost/algorithm/string.hpp>

Molecule::Molecule(){
	_moleculeInit();
}

void Molecule::_moleculeInit()
{
	_index = -1;
	_atoms.clear();
	_coordinatesReal.setConstant(-1.0);
	_transformMatrix = MatrixXd::Zero(0,0);
}

Molecule::Molecule( int index, vector<AtomPtr> atoms)
{
	_index = index;
	_atoms = atoms;
	_coordinatesReal = Vector3d::Zero();	// Treat the center of the atoms as position of the molecule.
	for ( size_t i = 0; i < _atoms.size(); ++i )
	{
		_coordinatesReal += _atoms[i]->getCoordinatesReal()/_atoms.size();
	}
	_pointGroup = boost::make_shared<PointGroup>(_atoms);
}
Molecule::~Molecule(){}

void Molecule::findPointGroup()
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	if ( suggestedPointGroup != "" )
	{
		try{
			findManuallyDefinedPointGroup();
			cleanUpPointGroup();

			string name = _pointGroup->getName();

			boost::algorithm::to_lower(name);

			if ( name == suggestedPointGroup ) // sucessful
			{
				cout << "Manually selected point group " << suggestedPointGroup << " found" << endl;
				return;
			}
			else
			{
				cout << "Manually selected point group " << suggestedPointGroup << " not found" << endl;
				cout << "Based on your input, I found " << name << endl;
				cout << "Please check if your choice is reasonable or try adjusting the symmetry precision" << endl;
				cout << "Falling back to automatic determination" << endl;
				_pointGroup->clear();
			}
		}
		catch(...)
		{
			cout << "Crashed during search for " << suggestedPointGroup << endl;
			cout << "Falling back to automatic determination" << endl;
		}
	}

	findSymmetryElements();

	cleanUpPointGroup();
}

void Molecule::findManuallyDefinedPointGroup()
{
	string pg = suggestedPointGroupType;
	size_t order = intFromString(suggestedPointGroup);
	if ( pg == "td" ) order = 3;
	if ( pg == "oh" ) order = 4;
	if ( pg == "ih" ) order = 5;

	if (debugLevel > 0) cout << "checking for manually defined point group " << suggestedPointGroup << " of type " << suggestedPointGroupType << " and order " << order << endl;
	// C1
	checkIdentity();

	// has inversion? -> Ci; Cnh, n even;
	if ( pg == "ci" ||
			( pg == "cnh" && order % 2 == 0) ||
			pg == "sn" ||
			( pg == "dnh" && order % 2 == 0 ) ||
			( pg == "dnd" && order % 2 == 1) ||
			pg == "oh" || pg == "ih"
	)
		checkInversion();

	// has Cn?
	if ( !( pg == "c1" || pg == "ci" || pg == "cs" ) )
		findRotationAxes(order);

	// has n C2?
	if ( pg == "dn" || pg == "dnh" || pg == "dnd" || pg == "td" || pg == "oh" || pg == "ih" )
		searchC2Axes();

	// has Sn?
	if ( pg == "sn" || pg == "cnh" || pg == "dnh" || pg == "dnd" || pg == "td" || pg == "oh" || pg == "ih" )
		findRotationReflection(order);

	// has sigma_h?
	if ( pg == "cnh" || pg == "dnh" || pg == "oh")
		searchSigmaH();

	// has sigma_v or sigma_d?
	if ( pg == "cnv" || pg == "dnd" || pg == "dnh" || pg == "td" || pg == "ih" )
		searchSigmaV();

	if ( isLinear() )
		searchSigmaLinear();

	// has sigma of other type?
	if ( pg == "cs" || pg == "ih" )
		searchSigma();
}

void Molecule::findSymmetryElements()
{
	if ( isLinear() )
	{
		findSymmetryElementsLinear();
		return;
	}

	checkIdentity();
	checkInversion();
	findRotationAxes();
	searchC2Axes();

	findRotationReflection();

	findReflection();
}

void Molecule::findSymmetryElementsLinear()
{
	checkIdentity();
	checkInversion();
	findRotationAxesLinear();
	findRotationReflectionLinear();
	searchSigmaLinear();
	searchC2AxesLinear();
}

void Molecule::cleanUpPointGroup()
{
	_pointGroup->setSymmetryOperationClasses();

	_pointGroup->findNameOfPointGroup();
	_pointGroup->generateCharacterTable();

	// set primes according to coordinate system
	_pointGroup->setPrimes();

	_pointGroup->setSymmetryOperationClasses();
	if ( isLinear() )
		_pointGroup->getCharacterTable()->combineC2AndSigmaV();
	_pointGroup->getCharacterTable()->sortClasses();

	//_pointGroup->generateCharacterTable();
	_pointGroup->generateIrreps();
	_pointGroup->getCharacterTable()->generateMullikenSymbols();
	_pointGroup->getCharacterTable()->sortIrreps();
}

// for debugging only. remove later
void Molecule::printAllSymmetryOperations()
{
	for ( auto so : _pointGroup->getSymmetryOperations() )
	{
		cout << "         found " << so->getSymbol() << endl;
	}
}

PGPtr Molecule::getPointGroup()
{
	return _pointGroup;
}

Vector3d Molecule::getCenter()
{
	return _coordinatesReal;
}

void Molecule::setTransformMatrix(MatrixXd transformMatrix)
{
	_transformMatrix = transformMatrix;
}

vector<AtomPtr> Molecule::getAtoms()
{
	return _atoms;
}

int Molecule::getIndex()
{
	return _index;
}

AtomPtr Molecule::getAtom( int index )
{
	return _atoms[index];
}

AtomPtr Molecule::getAtom( size_t index )
{
	return _atoms[index];
}

size_t Molecule::getNumberOfAtoms()
{
	return _atoms.size();
}

string Molecule::getName()
{
	string result;

	vector<pair<string, int> > elements;

	for ( size_t i = 0; i < _atoms.size(); ++i )
	{
		addAtomToName(i, elements);
	}


	for ( size_t i = 0; i < elements.size(); ++i )
	{
		result += elements[i].first;
		if ( elements[i].second != 1 )
			result += to_string(elements[i].second);
	}

	return result;
}

void Molecule::addAtomToName(int i, vector<pair<string,int> >& elements)
{
	for ( size_t j = 0; j < elements.size(); ++j )
	{
		if ( elements[j].first == _atoms[i]->getElementSymbolForPrinting() )
		{
			++elements[j].second;
			return;
		}
	}

	elements.push_back(make_pair(_atoms[i]->getElementSymbolForPrinting(), 1));
}

double Molecule::distanceFrom(MolPtr other) // todo: include translation
{
	size_t thisAtoms = this->getNumberOfAtoms();
	size_t otherAtoms = other->getNumberOfAtoms();

	MatrixXd results = MatrixXd::Zero( thisAtoms, otherAtoms );

	for ( size_t i = 0; i < thisAtoms; ++i )
	{
		AtomPtr thisAtom = _atoms[i];
		for ( size_t j = 0; j < otherAtoms; ++j )
		{
			AtomPtr otherAtom = other->getAtom(j);
			results(i,j) = thisAtom->distanceFrom(otherAtom);
		}
	}

	return results.minCoeff();
}

bool Molecule::operator==(const Molecule& other) const {
	if ( _atoms.size() != other._atoms.size() )
		return false;

	VectorXd bools = VectorXd::Ones(_atoms.size());

	for ( size_t a = 0; a < _atoms.size(); ++a )
	{
		AtomPtr atomA = _atoms[a];
		for ( size_t b = 0; b < _atoms.size(); ++b )
		{
			AtomPtr atomB = other._atoms[b];

			if ( (*atomA) == (*atomB) )
			{
				bools( atomA->getIndex() ) = 0;
			}
		}
		
		if ( bools( atomA->getIndex() ) != 0 )
			return false;
	}

	if ( bools.sum() == 0 )
		return true;
	else
		return false;
}

bool Molecule::operator !=(const Molecule& other) const {
	return !(*this == other);
}

bool Molecule::isLinear()
{
	if ( _atoms.size() == 1 ) // todo : prevent single atoms from becoming molecules
		return false;

	if ( _atoms.size() == 2 )
		return true;

	Vector3d ref = _atoms[0]->getCoordinatesReal() - _atoms[1]->getCoordinatesReal();
	ref.normalize();

	for ( size_t i = 2; i < _atoms.size(); ++i )
	{
		Vector3d v2 = _atoms[0]->getCoordinatesReal() - _atoms[i]->getCoordinatesReal();
		v2.normalize();

		if ( abs(ref.dot(v2)) < 0.999 )
			return false;
	}

	return true;
}

bool Molecule::isPlanar()
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	if ( _atoms.size() < 3 )
		return false;

	Vector3d normalVector = Vector3d::Zero();

	int j = 2;
	while (true )
	{
		normalVector = (_atoms[0]->getCoordinatesReal() - _atoms[1]->getCoordinatesReal()).cross(_atoms[0]->getCoordinatesReal() - _atoms[j]->getCoordinatesReal());
		normalVector.normalize();
		++j;
		if (normalVector.norm() > 0 || j >= _atoms.size())
			break;
	}

	if ( normalVector.norm() < symmetryPrecision ) // linear molecule
		return false;


	for ( int i = 0; i < _atoms.size(); ++i )
	{
		Vector3d point = (_atoms[i]->getCoordinatesReal() - _coordinatesReal);
		point.normalize();

		if ( point.dot(normalVector) > symmetryPrecision )
			return false;
	}

	return true;
}

double Molecule::getErrorOfSymmetryOperation( SymmetryOperation so )
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	vector<AtomPtr> newAtoms;

	for ( auto a : _atoms )
	{
		Vector3d coordinatesNew = so.applySymmetryOperation(a->getCoordinatesReal());

		AtomPtr newAtom = boost::make_shared<Atom>( a->getIndex(), a->getElementSymbol(), coordinatesNew );
		newAtoms.push_back(newAtom);
	}

	MolPtr other = boost::make_shared<Molecule>( _index, newAtoms );
	double result = 0;
	{
		for ( size_t a = 0; a < _atoms.size(); ++a )
		{
			AtomPtr atomA = _atoms[a];
			VectorXd dist = VectorXd::Zero( _atoms.size() );
			for ( size_t b = 0; b < _atoms.size(); ++b )
			{
				AtomPtr atomB = other->getAtoms()[b];
				if ( atomA->getElementSymbolForPrinting() == atomB->getElementSymbolForPrinting() )
					dist(atomB->getIndex()) = atomA->distanceFrom(atomB);
				else
					dist(atomB->getIndex()) = 9999;
			}
				result += dist.minCoeff();
		}
	}

	return result / (double) _atoms.size();
}

void Molecule::checkIdentity()
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	if ( debugLevel > 0 ) cout << "Checking identity" << endl;
	SymmOpPtr so = boost::make_shared<SymmetryOperation>("identity", _coordinatesReal, Vector3d::Zero());

	double error = getErrorOfSymmetryOperation(*so);

	if ( error < symmetryPrecision )
	{
		so->setError(error);
		_pointGroup->addSymmetryOperation(so);
	}
}

void Molecule::checkInversion()
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	if ( debugLevel > 0 ) cout << "Checking inversion" << endl;
	SymmOpPtr so = boost::make_shared<SymmetryOperation>("inversion", _coordinatesReal);
	double error = getErrorOfSymmetryOperation(*so);

	if ( error < symmetryPrecision )
	{
		so->setError(error);
		_pointGroup->addSymmetryOperation(so);
	}
}


// returns the maximum number of atoms with the same element
size_t Molecule::getMaximumPossibleRotationAxis( size_t order )
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}

	if (  _pointGroup->getMainRotationAxes().size() != 0 )
	{
		if ( debugLevel > 0 ) cout << "max rotation axis is " << _pointGroup->getMainRotationAxes()[0]->getOrderOfRotation() << " according to mra" << endl;
		return _pointGroup->getMainRotationAxes()[0]->getOrderOfRotation();
	}

	if (order != 0)
		return order;

	map<string, size_t> elementTable;

	for ( auto atom : _atoms )
	{
		++elementTable[atom->getElementSymbol()];
	}

	return max_element( elementTable.begin(), elementTable.end(),[] (const std::pair<string,size_t>& a, const std::pair<string,size_t>& b)->bool{ return a.second < b.second; } )->second;
}


void Molecule::findRotationAxes( size_t order )
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	if ( debugLevel > 0 ) cout << "Looking for rotation axes with order " << order << endl;
	searchMainRotationAxis(order);
	checkForMainRotationAxis();
	_pointGroup->determineMainRotationAxis();
}

void Molecule::findRotationAxesLinear( size_t order )
{
	if ( debugLevel > 0 ) cout << "searching linear Cn" << endl;
	_pointGroup->setLinear();

	findRotationAxes(infinityRotationSteps);

/*	Vector3d axis = _atoms[0]->getCoordinatesReal() - _atoms[1]->getCoordinatesReal();
	axis.normalize();

	int n = infinityRotationSteps;

	for ( int i = n; i > 2; --i )
	{
		SymmOpPtr so = boost::make_shared<SymmetryOperation>( "rotation", _coordinatesReal, axis, i);
		double error = getErrorOfSymmetryOperation(*so);

		if ( error < symmetryPrecision )
		{
			so->setError(error);
			_pointGroup->addSymmetryOperation(so);
		}
	}*/

	//_pointGroup->determineMainRotationAxis();

	if ( debugLevel > 0 ) cout << "Done with linear Cn" << endl;
}

void Molecule::searchMainRotationAxis( size_t order )
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	// determine the maximum possible order of rotation
	const size_t maxOrder = getMaximumPossibleRotationAxis(order);
	if ( debugLevel > 0 ) cout << "Checking for mra with max " << maxOrder << endl;

	if ( isPlanar() )
	{
		for ( size_t i = maxOrder; i > 1; --i )
		{
			Vector3d normalVector = Vector3d::Zero();

			int j = 2;
			while (true )
			{
				normalVector = (_atoms[0]->getCoordinatesReal() - _atoms[1]->getCoordinatesReal()).cross(_atoms[0]->getCoordinatesReal() - _atoms[j]->getCoordinatesReal());
				normalVector.normalize();
				++j;
				if (normalVector.norm() > 0 || j >= _atoms.size())
					break;
			}

			SymmOpPtr so = boost::make_shared<SymmetryOperation>( "rotation", _coordinatesReal, normalVector, i );

			double error = getErrorOfSymmetryOperation(*so);

			if ( error < symmetryPrecision )
			{
				so->setError(error);
				_pointGroup->addSymmetryOperation(so);
				return;
			}
		}
	}

	for ( size_t i = maxOrder; i > 1; --i )
	{
		if ( maxOrder % i != 0 )
			continue;

		if (debugLevel > 0) cout << i << " of " << maxOrder << endl;
		for ( size_t a = 0; a < _atoms.size(); ++a )
		{
			for ( size_t b = a; b < _atoms.size(); ++b )
			{
				Vector3d rotationAxis = (_atoms[a]->getCoordinatesReal() + _atoms[b]->getCoordinatesReal()) / 2 - _coordinatesReal;

				rotationAxis.normalize();

				if ( rotationAxis.norm() > 1e-12 )
				{
					SymmOpPtr so = boost::make_shared<SymmetryOperation>( "rotation", _coordinatesReal, rotationAxis, i );

					double error = getErrorOfSymmetryOperation(*so);

					if ( error < symmetryPrecision )
					{
						so->setError(error);
						_pointGroup->addSymmetryOperation(so);
					}
				}
			}
		}

		for ( size_t a = 0; a < _atoms.size(); ++a )
		{
			for ( size_t b = a + 1; b < _atoms.size(); ++b )
			{
				for ( size_t c = b + 1; c < _atoms.size(); ++c )
				{
					Vector3d rotationAxis = (_atoms[b]->getCoordinatesReal() - _atoms[a]->getCoordinatesReal()).cross(_atoms[c]->getCoordinatesReal() - _atoms[a]->getCoordinatesReal());
					rotationAxis.normalize();

					if ( rotationAxis.norm() > 1e-12 )
					{
						SymmOpPtr so = boost::make_shared<SymmetryOperation>( "rotation", _coordinatesReal, rotationAxis, i );
						double error = getErrorOfSymmetryOperation(*so);

						if ( error < symmetryPrecision )
						{
							so->setError(error);
							_pointGroup->addSymmetryOperation(so);
						}
					}
				}
			}
		}
	}
}

void Molecule::searchC2Axes()
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	if ( debugLevel > 0 ) cout << "Checking for C2" << endl;
	const int order = 2;

	for ( size_t a = 0; a < _atoms.size(); ++a )
	{
		for ( size_t b = a; b < _atoms.size(); ++b )
		{
			Vector3d rotationAxis = (_atoms[a]->getCoordinatesReal() + _atoms[b]->getCoordinatesReal()) / 2 - _coordinatesReal;

			if ( rotationAxis.norm() > 1e-12 )
			{
				SymmOpPtr so = boost::make_shared<SymmetryOperation>( "rotation", _coordinatesReal, rotationAxis, order );

				double error = getErrorOfSymmetryOperation(*so);

				if ( error < symmetryPrecision )
				{
					so->setError(error);
					_pointGroup->addSymmetryOperation(so);
				}
			}
		}
	}

	for ( size_t a = 0; a < _atoms.size(); ++a )
	{
		for ( size_t b = a + 1; b < _atoms.size(); ++b )
		{
			for ( size_t c = b + 1; c < _atoms.size(); ++c )
			{
				Vector3d rotationAxis = (_atoms[b]->getCoordinatesReal() - _atoms[a]->getCoordinatesReal()).cross(_atoms[c]->getCoordinatesReal() - _atoms[a]->getCoordinatesReal());

				if ( rotationAxis.norm() > 1e-12 )
				{
					SymmOpPtr so = boost::make_shared<SymmetryOperation>( "rotation", _coordinatesReal, rotationAxis, order );
					double error = getErrorOfSymmetryOperation(*so);

					if ( error < symmetryPrecision )
					{
						so->setError(error);
						_pointGroup->addSymmetryOperation(so);
					}
				}
			}
		}
	}
}

void Molecule::searchC2AxesLinear()
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	if ( debugLevel > 0 ) cout << "Checking for C2 in linear" << endl;

	for ( auto it : _pointGroup->getSymmetryOperations() )
	{
		if (it->getType() != "reflection")
			continue;

		Vector3d v = it->getAxis();

		if ( debugLevel > 0 ) cout << "vector " << v.transpose() << endl;

		SymmOpPtr so = boost::make_shared<SymmetryOperation>( "rotation", _coordinatesReal, v, 2 );

		double error = getErrorOfSymmetryOperation(*so);

		if ( error < symmetryPrecision )
		{
			if ( debugLevel > 0 ) cout << "adding C2" << endl;
			so->setError(error);
			_pointGroup->addSymmetryOperation(so);
		}
	}
}

void Molecule::checkForMainRotationAxis()
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	if ( _pointGroup->hasMultipleC2())
	{
		for ( auto so1 : _pointGroup->getSymmetryOperations() )
		{
			for ( auto so2 : _pointGroup->getSymmetryOperations() )
			{
				Vector3d v1, v2;

				if ( (so1->getType() == "rotation") && ( so1->getOrderOfRotation() == 2 ) && (so2->getType() == "rotation") && ( so2->getOrderOfRotation() == 2 ) && ( so1->getAxis() != so2->getAxis() ) )
				{
					checkForMainRotationAxis( Vector3d(so1->getAxis().cross(so2->getAxis())) );

					return;
				}
			}
		}
	}
	else
	{
		for ( auto so1 : _pointGroup->getSymmetryOperations() )
		{
			if ( (so1->getType() == "rotation") && ( so1->getOrderOfRotation() == 2 ))
			{
				checkForMainRotationAxis(so1->getAxis());

				return;
			}
		}
	}
}

void Molecule::checkForMainRotationAxis( Vector3d input )
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	// determine the maximum possible order of rotation
	const size_t maxOrder = getMaximumPossibleRotationAxis();

	for ( size_t order = maxOrder; order > 1; --order )
	{
		if ( input.norm() > 1e-12 )
		{
			SymmOpPtr so = boost::make_shared<SymmetryOperation>( "rotation", _coordinatesReal, input, order );
			double error = getErrorOfSymmetryOperation(*so);

			if ( error < symmetryPrecision )
			{
				so->setError(error);
				_pointGroup->addSymmetryOperation(so);
			}
		}
	}
}

void Molecule::findReflection()
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	if ( debugLevel > 0 ) cout << "Searching reflections" << endl;
	if ( _pointGroup->getOrderOfMainRotationAxis() < 2 )
	{
			searchSigmaH();
			searchSigmaV();
	}
	else
	{
		searchSigma();
	}
}

void Molecule::searchSigmaH()
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	if ( debugLevel > 0 ) cout << "Checking for sigma h" << endl;
	// take main rotation axis and see if there is a reflection plane with the same normal vector as the rotation axis
	for ( auto rot : _pointGroup->getSymmetryOperations() )
	{
		if ( rot->getType() == "rotation" )
		{
			if ( rot->getOrderOfRotation() == _pointGroup->getOrderOfMainRotationAxis() )
			{
				SymmOpPtr so = boost::make_shared<SymmetryOperation>( "reflection", _coordinatesReal, rot->getAxis() );

				double error = getErrorOfSymmetryOperation(*so);

				if ( error < symmetryPrecision )
				{
					so->setError(error);
					_pointGroup->addSymmetryOperation(so);
				}
			}
		}
	}
}

void Molecule::searchSigmaV()
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	if ( debugLevel > 0 ) cout << "searching sigma v" << endl;
	// take main rotation axis and search reflection planes with a normal vector perpendicular to the axis
	for ( auto mra : _pointGroup->getMainRotationAxes() )
		searchSigmaVAround(mra->getAxis());
}

void Molecule::searchSigmaVAround( Vector3d axis )
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	for ( auto atom1 : _atoms )
	{
		Vector3d normalVector = (atom1->getCoordinatesReal() - _coordinatesReal).cross(axis);
		if (normalVector.norm() < 1e-3)
			continue;

		SymmOpPtr so = boost::make_shared<SymmetryOperation>( "reflection", _coordinatesReal, normalVector );

		double error = getErrorOfSymmetryOperation(*so);

		if ( error < symmetryPrecision )
		{
			so->setError(error);
			_pointGroup->addSymmetryOperation(so);
		}
	}

	for ( int i = 0; i < _atoms.size(); ++i )
	{
		AtomPtr atom1 = _atoms[i];
		for ( int j = i; j < _atoms.size(); ++j )
		{
			AtomPtr atom2 = _atoms[j];

			Vector3d normalVector = (atom1->getCoordinatesReal() - atom2->getCoordinatesReal());//.cross(axis);
			if (normalVector.norm() < 1e-3)
				continue;

			SymmOpPtr so = boost::make_shared<SymmetryOperation>( "reflection", _coordinatesReal, normalVector );

			double error = getErrorOfSymmetryOperation(*so);

			if ( error < symmetryPrecision )
			{
				so->setError(error);
				_pointGroup->addSymmetryOperation(so);
			}
		}
	}
}

void Molecule::searchSigma()
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	if (debugLevel > 0) cout << "searching for sigmas the old way" << endl;

	for ( int i = 0; i < _atoms.size(); ++i )
	{
		AtomPtr atom1 = _atoms[i];
		for ( int j = i; j < _atoms.size(); ++j )
		{
			AtomPtr atom2 = _atoms[j];
			Vector3d normalVector = (atom1->getCoordinatesReal() - atom2->getCoordinatesReal());
			if (normalVector.norm() < 1e-3)
				continue;

			SymmOpPtr so = boost::make_shared<SymmetryOperation>( "reflection", _coordinatesReal, normalVector );

			double error = getErrorOfSymmetryOperation(*so);

			if ( error < symmetryPrecision )
			{
				so->setError(error);
				_pointGroup->addSymmetryOperation(so);
			}
		}
	}

	for ( int i = 0; i < _atoms.size(); ++i )
	{
		AtomPtr atom1 = _atoms[i];
		for ( int j = i; j < _atoms.size(); ++j )
		{
			AtomPtr atom2 = _atoms[j];
			for ( int k = j;  k < _atoms.size(); ++k )
			{
				AtomPtr atom3 = _atoms[k];

				Vector3d normalVector = (atom1->getCoordinatesReal() - atom2->getCoordinatesReal()).cross(atom1->getCoordinatesReal() - atom3->getCoordinatesReal());

				if ( debugLevel > 0 ) cout << "checking for sigma with " << normalVector.transpose() << endl;

				if (normalVector.norm() < 1e-3)
					continue;

				SymmOpPtr so = boost::make_shared<SymmetryOperation>( "reflection", _coordinatesReal, normalVector );

				double error = getErrorOfSymmetryOperation(*so);

				if ( error < symmetryPrecision )
				{
					if ( debugLevel > 0 ) cout << "added sigma with " << normalVector.transpose() << endl;

					so->setError(error);
					_pointGroup->addSymmetryOperation(so);
				}
			}
		}
	}
}

void Molecule::searchSigmaLinear()
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}

	if ( debugLevel > 0 ) cout << "Checking reflections in linear" << endl;

	// sigma_h
	for ( auto rot : _pointGroup->getSymmetryOperations() )
	{
		if ( rot->getType() == "rotation" )
		{
			SymmOpPtr so = boost::make_shared<SymmetryOperation>( "reflection", _coordinatesReal, rot->getAxis() );

			double error = getErrorOfSymmetryOperation(*so);

			if ( error < symmetryPrecision )
			{
				so->setError(error);
				_pointGroup->addSymmetryOperation(so);
			}
		}
	}

	// then sigma_v; as all of them are equal, we only need one of them
	Vector3d normalVector = _atoms[0]->getCoordinatesReal() - _atoms[1]->getCoordinatesReal();
	normalVector.normalize();

	Vector3d v2 = Vector3d::Zero();
	Vector3d v1 = Vector3d::Zero();
	int i = 0;

	while ( true )
	{
		v2 = Vector3d{rand() % 100, rand() % 100, rand() % 100};
		v2.normalize();
		v1 = v2.cross(normalVector);
		if ( v1.dot(normalVector) < symmetryPrecision || i == 2 )
			break;
		++i;
	}

	SymmetryOperation E ("identity");

	int orderForRotation = infinityRotationSteps;
	if ( orderForRotation % 2 == 0 )
		orderForRotation *= 2;

	SymmetryOperation Cn ("rotation",_pointGroup->getMainRotationAxis()->getCenter(),_pointGroup->getMainRotationAxis()->getAxis(), orderForRotation);

	for ( size_t i = 0; i < _pointGroup->getMainRotationAxis()->getOrderOfRotation(); ++i )
	{
		SymmOpPtr so = boost::make_shared<SymmetryOperation>( "reflection", _coordinatesReal, E.getTransformMatrix() * v1 );

		if ( debugLevel > 0 ) cout << "trying to add " << so->getAxis().transpose() << endl;

		double error = getErrorOfSymmetryOperation(*so);

		if ( error < symmetryPrecision )
		{
			if ( debugLevel > 0 ) cout << "added" << endl;
			so->setError(error);
			_pointGroup->addSymmetryOperation(so);
		}

		E = E * Cn;
	}

	if ( debugLevel > 0 ) cout << "sigma_v check done" << endl;
}

void Molecule::findRotationReflection( size_t order )
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	if ( debugLevel > 0 ) cout << "Checking Sn" << endl;

	size_t snOrder = getMaximumPossibleRotationAxis(order);
	if ( snOrder % 2 == 1 )
	{
		if ( debugLevel > 0 ) cout << "main rotation axis has an odd order " << snOrder;
		snOrder = 2 * getMaximumPossibleRotationAxis(order);
		if ( debugLevel > 0 ) cout << " therefore, I will use with order " << snOrder << " for the Sn" << endl;
	}

	if ( debugLevel > 0 ) cout << "with order " << snOrder << endl;

	// determine the maximum possible order of rotation
	const size_t maxOrder = snOrder;

	for ( size_t i = maxOrder; i > 1; --i )
	{
		if ( maxOrder % i != 0 )
			continue;

		if ( debugLevel > 0 ) cout << "considering S" << i << endl;

		for ( auto cn : _pointGroup->getSymmetryOperations() )
		{
			if ( cn->getType() != "rotation" )
				continue;

			if ( debugLevel > 0 ) cout << "around " << cn->getAxis().transpose() << endl;
			SymmOpPtr so = boost::make_shared<SymmetryOperation>( "rotationReflection", _coordinatesReal, cn->getAxis(), i );
			double error = getErrorOfSymmetryOperation(*so);

			if ( error < symmetryPrecision )
			{
				if ( debugLevel > 0 ) cout << "and found it." << endl;
				so->setError(error);
				_pointGroup->addSymmetryOperation(so);
				if ( isLinear() )
					return;
			}
			else
				if (debugLevel > 0) cout << "and rejected it with error " << error << endl;
		}
	}
}

void Molecule::findRotationReflectionLinear()
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	if ( debugLevel > 0 ) cout << "Checking Sn linear" << endl;

	findRotationReflection(infinityRotationSteps);
}

void Molecule::checkForMainRotationReflectionAxis()
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	if( isLinear() )
		return;

	if ( _pointGroup->hasMultipleC2())
	{
		for ( auto so1 : _pointGroup->getSymmetryOperations() )
		{
			for ( auto so2 : _pointGroup->getSymmetryOperations() )
			{
				Vector3d v1, v2;

				if ( (so1->getType() == "rotation") && ( so1->getOrderOfRotation() == 2 ) && (so2->getType() == "rotation") && ( so2->getOrderOfRotation() == 2 ) && ( so1->getAxis() != so2->getAxis() ) )
				{
					checkForMainRotationReflectionAxis( Vector3d(so1->getAxis().cross(so2->getAxis())) );

					return;
				}
			}
		}
	}
}

void Molecule::checkForMainRotationReflectionAxis( Vector3d input )
{
if( debugLevel > 3 ) {
	cout << "Molecule::" << __FUNCTION__ << endl;
}
	// determine the maximum possible order of rotation
	const size_t maxOrder = getMaximumPossibleRotationAxis();

	for ( size_t order = maxOrder; order > 1; --order )
	{
		if ( input.norm() > 1e-12 )
		{
			SymmOpPtr so = boost::make_shared<SymmetryOperation>( "rotationReflection", _coordinatesReal, input, order );
			double error = getErrorOfSymmetryOperation(*so);

			if ( error < symmetryPrecision )
			{
				so->setError(error);
				_pointGroup->addSymmetryOperation(so);
			}
		}
	}
}
