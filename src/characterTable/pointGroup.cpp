/*
 * pointGroup.cpp
 *
 *  Created on: Feb 6, 2023
 *      Author: pemueller
 */
#include "pointGroup.h"
#include <math.h>
#include <numeric>
#include <iostream>
#include "../common/constants.h"
#include <boost/make_shared.hpp>
#include <Eigen/Dense>

PointGroup::PointGroup( vector<AtomPtr> atoms )
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	_symmetryOperations.clear();
	_name = "";
	_mainRotationAxis.resize(0);
	_atoms = atoms;
}

PointGroup::~PointGroup(){
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
}

void PointGroup::clear()
{
	_symmetryOperations.clear();
	_mainRotationAxis.clear();
	_name ="";
	_type = "";
	_characterTable->clear();
}

void PointGroup::setLinear()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	_isLinear = true;
}

void PointGroup::generateCharacterTable()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}

	if ( _mainRotationAxis.size() == 0 )
		_mainRotationAxis.push_back( boost::make_shared<SymmetryOperation>("identity") );

	_characterTable = boost::make_shared<CharacterTable>(_type, _mainRotationAxis[0]);

	vector<SymmOpPtr> soVector = _symmetryOperations;

	while ( soVector.size() > 0 )
	{
		SymmOpPtr so1 = soVector[0];

		vector<SymmOpPtr> newClass;

		if ( so1->getType() == "rotation" && so1->getOrderOfRotation() == 0 )
		{
			newClass.push_back(so1);
			soVector.erase(soVector.begin());
			_characterTable->addClass(newClass);
			continue;
		}

		for ( auto so2 : _symmetryOperations )
		{
			MatrixXd similarityTransform = so2->getTransformMatrix() * so1->getTransformMatrix() * so2->getTransformMatrix().inverse();

			for ( int i = 0; i < soVector.size(); ++i )
			{
				if (( similarityTransform - soVector[i]->getTransformMatrix()).cwiseAbs().sum() / (double) _symmetryOperations.size() < symmetryPrecision)
				{
					newClass.push_back(soVector[i]);
					soVector.erase(soVector.begin() + i);
				}
			}
		}

		_characterTable->addClass(newClass);
	}
}

CTPtr PointGroup::getCharacterTable()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	return _characterTable;
}

vector<SymmOpPtr> PointGroup::getSymmetryOperations()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	return _symmetryOperations;
}

// not needed probably
MatrixXs PointGroup::getSymmetryElementMatrix()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	MatrixXs result(boost::extents[_symmetryOperations.size()][_symmetryOperations.size()]);

	for ( size_t i = 0; i < _symmetryOperations.size(); ++i )
	{
		for ( size_t j = i; j < _symmetryOperations.size(); ++j )
		{
			string symbol = "NULL";

			MatrixXd similarityTransform = _symmetryOperations[i]->getTransformMatrix() * _symmetryOperations[j]->getTransformMatrix();// *_symmetryOperations[i]->getTransformMatrix().inverse();

			for ( size_t k = 0; k < _symmetryOperations.size(); ++k )
			{

				if ( (similarityTransform - _symmetryOperations[k]->getTransformMatrix()).cwiseAbs().sum() / (double)_symmetryOperations.size() / (double)_symmetryOperations.size() < symmetryPrecision )
				{
					symbol = _symmetryOperations[k]->getSymbol() + to_string(k+1);
				}
			}

			result[i][j] = symbol;
			result[j][i] = symbol;
		}

	}

	return result;
}

void PointGroup::determineMainRotationAxis()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	int order = 0;

	// We look for the highest Cn here. If point group is S2n, it will be updated later.
	for ( auto so : _symmetryOperations )
	{
		if ( so->getType() == "rotation" && so->getOrderOfRotation() == 0 ) // in case we have a C_infty (treated as C_0), break here
			break;

		if ( so->getType() == "rotation" && so->getOrderOfRotation() > order )
		{
			order = so->getOrderOfRotation();
		}
	}

	for ( auto so : _symmetryOperations )
	{
		if (  so->getOrderOfRotation() == order && !isAxisInMainRotationAxis(so) )
		{
			_mainRotationAxis.push_back(so);
		}
	}

	if ( debugLevel > 0 ) cout << "found " << _mainRotationAxis.size() << " main rotation axes with order " << _mainRotationAxis[0]->getOrderOfRotation() << endl;
}

bool PointGroup::isAxisInMainRotationAxis( SymmOpPtr axis )
{
	for ( auto so : _mainRotationAxis )
	{
		if ( abs(axis->getAxis().dot(so->getAxis())) > 0.9999  )
			return true;
	}

	return false;
}

void PointGroup::setSymmetryOperationClasses()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	if ( debugLevel > 0 ) cout << "setting symmetry operation classes for " << _type << endl;
	for ( auto so : _symmetryOperations )// todo create fallback for D2h
	{
		if ( so->getType() == "rotation" )
		{
			string soClass = "C"+to_string(so->getOrderOfRotation())+"^"+to_string(so->getNumberOfApplications());
			for ( int i = 0; i < so->getNumberOfPrimes(); ++i )
				soClass+='\'';

			so->setSymbol(soClass);
		}
		if ( so->getType() == "rotationReflection" )
		{
			string numberOfApplicationsAsString = to_string(so->getNumberOfApplications());
			if ( so->getNumberOfApplications() < 0 && so->getOrderOfRotation() % 2 == 0)
			{
				if (debugLevel > 1)
					cout << "Number of applications is below zero: " << so->getNumberOfApplications() << endl <<
					"Fixing it to " << so->getNumberOfApplications() + so->getOrderOfRotation() << endl;
				numberOfApplicationsAsString = to_string( so->getNumberOfApplications() + so->getOrderOfRotation() );
			}

			string soClass = "S"+to_string(so->getOrderOfRotation())+"^"+numberOfApplicationsAsString;
			so->setSymbol(soClass);
		}
		if ( so->getType() == "identity" )
		{
			string soClass = "E";
			so->setSymbol(soClass);
		}
		if ( so->getType() == "inversion")
		{
			string soClass = "i";
			so->setSymbol(soClass);
		}
		if ( so->getType() == "reflection" )
		{
			string soClass = "sigma";

			if ( isSigmaV(so) )
				soClass = "sigma_v";
			if ( isSigmaD(so) )
				soClass = "sigma_d";
			if ( isSigmaH(so) )
				soClass = "sigma_h";

			for ( int i = 0; i < so->getNumberOfPrimes(); ++i )
				soClass+='\'';

			so->setSymbol(soClass);
		}
	}
}

void PointGroup::setPrimes()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	// this is only important if the same kind of symmetry operation appears in multiple classes.
	// happens only in Cnv, Dn, Dnh and Dnd with n even
	if ( (_type == "Cnv"/* || _type == "Cuv"*/) && _mainRotationAxis[0]->getOrderOfRotation() % 2 == 0 )
		setPrimesCnv();
	if ( _type == "Dn" && _mainRotationAxis[0]->getOrderOfRotation() % 2 == 0 )
		setPrimesDn();
	if ( (_type == "Dnh" /*|| _type == "Duh"*/) && _mainRotationAxis[0]->getOrderOfRotation() % 2 == 0 )
		setPrimesDndh();
	if ( _type == "Dnd" && _mainRotationAxis[0]->getOrderOfRotation() % 2 == 0 && _mainRotationAxis[0]->getOrderOfRotation() != 2 ) // exclude D2d
		setPrimesDndh();
}


double PointGroup::getDistanceFromSO(AtomPtr atom, SymmOpPtr so)
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	if ( so->getType() == "reflection" )
		return ( abs((atom->getCoordinatesReal() - so->getCenter()).dot(so->getAxis())) );

	if ( so->getType() == "rotation" )
		return ( ( so->getAxis().cross( so->getCenter() - atom->getCoordinatesReal() ) ).norm() );
}

VectorXd PointGroup::getDistanceFromSO( SymmOpPtr so )
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	VectorXd result = VectorXd::Zero( _atoms.size() );

	if ( so->getType() == "reflection" )
		for ( size_t i = 0; i < _atoms.size(); ++i )
			result(i) = abs(_atoms[i]->getCoordinatesReal().dot(so->getAxis()) + so->getAxis().dot(so->getCenter()));

	if ( so->getType() == "rotation" )
		for ( size_t i = 0; i < _atoms.size(); ++i )
			result (i) = ( so->getAxis().cross( so->getCenter() - _atoms[i]->getCoordinatesReal() ) ).norm();

	return result;
}

double PointGroup::weightWithAtom(double distance, AtomPtr atom)
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	double result = 1;

	result = result + distance;

	return result / (double) atom->getAtomicNumber();
}

VectorXd PointGroup::weightWithAtom(VectorXd distances)
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	VectorXd result = VectorXd::Ones(distances.size());

	result = result + distances;

	for ( int i = 0; i < distances.size(); ++i )
	{
		result (i) = result (i) / (double) _atoms[i]->getAtomicNumber();
	}

	return result;
}



bool PointGroup::isGreater(vector<SymmOpPtr> v1, vector<SymmOpPtr> v2)
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	if ( v1.size() != v2.size() )
		return v1.size() < v2.size();

	if ( v1[0]->getAxis().dot(_mainRotationAxis[0]->getAxis()) > 0.999  )
		return true;
	if ( v2[0]->getAxis().dot(_mainRotationAxis[0]->getAxis()) > 0.999  )
		return false;

	return weightWithAtom(getDistanceFromSO( v1[0] )).minCoeff() >  weightWithAtom(getDistanceFromSO( v2[0] )).minCoeff();
}

string PointGroup::getName()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	return _name;
}

string PointGroup::getType()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	return _type;
}



void PointGroup::addSymmetryOperation(SymmOpPtr so)
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	if ( !isSymmetryOperationPresent(so) )
	{
		if ( ((so->getType() == "rotation") || ( so->getType() == "rotationReflection" )) )
		{
			const int order = so->getOrderOfRotation();
			for ( int i = 1; i <= so->getOrderOfRotation(); ++i )
			{
				int d = gcd(order, i);
				if ( so->getType() == "rotationReflection" && ( order % 2 == 0 && i % 2 == 0) )
					continue;
				if ( order/d == 1 )
					continue;

				SymmOpPtr so2 = boost::make_shared<SymmetryOperation> ( so->getType(), so->getCenter(), so->getAxis(), order/d, i/d );

				string sym;
				if ( so->getType() == "rotationReflection" )
					sym = "S";
				else
					sym = "C";
				if ( debugLevel > 0)
					cout << "trying to add a " << sym << "_" << order/d << "^" << i/d << endl;
				if ( debugLevel > 1)
					cout << "transform matrix: " << endl << so2->getTransformMatrix() << endl;

				if ( !isSymmetryOperationPresent(so2) )
				{
					_symmetryOperations.push_back(so2);
					if ( debugLevel > 1)
						cout << "added a " << sym << "_" << order/d << "^" << i/d << endl << "transform matrix " << endl << so2->getTransformMatrix() << endl;
				}
				else
					if ( debugLevel > 0 )
						cout << "rejected a " << sym << "_" << order/d << "^" << i/d << endl;
				so2 = boost::make_shared<SymmetryOperation> ( so->getType(), so->getCenter(), so->getAxis(), order/d, -i/d );
				if ( !isSymmetryOperationPresent(so2) )
				{
					_symmetryOperations.push_back(so2);
					if ( debugLevel > 0 )
						cout << "added a " << sym << "_" << order/d << "^-" << i/d << endl;
					if ( debugLevel > 1 )
						cout << "transform matrix " << endl << so2->getTransformMatrix() << endl;
				}
				else
					if ( debugLevel > 0 )
						cout << "rejected a " << sym << "_" << order/d << "^-" << i/d << endl;
			}
		}
		else
		{
			_symmetryOperations.push_back(so);
		}
	}
}

bool PointGroup::isThereARotationOfHigherOrder(SymmOpPtr input)
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	for ( auto so : _symmetryOperations )
	{
		if( *so > *input )
			return true;
	}
	return false;
}

bool PointGroup::isSymmetryOperationPresent( SymmOpPtr input )
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	for ( auto so : _symmetryOperations )
	{
		if( debugLevel > 1 )
			cout << endl << input->getSymbol() << " - " << so->getSymbol() << endl << *input - *so << endl;

		if ( *input == *so )
		{
			if( debugLevel > 1 )
				cout << "are equal" << endl;

			return true;
		}
		else
		{
			if( debugLevel > 1 )
				cout << "are not equal" << endl;
		}
	}

	return false;
}

bool PointGroup::isSigmaV( SymmOpPtr input )
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	if ( _type == "Cs" )
		return false;

	for ( auto mra : _mainRotationAxis )
	{
		if ( abs(input->getAxis().dot( mra->getAxis()) ) < symmetryPrecision )
			return true;
	}

	return false;
}

bool PointGroup::isSigmaH( SymmOpPtr input )
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}

	if ( _type == "Cs" )
		return true;

	for ( auto so : _mainRotationAxis )
	{
		if ( abs(input->getAxis().dot( so->getAxis() )) > 1 - symmetryPrecision )
			return true;
	}

	return false;
}

bool PointGroup::isSigmaD( SymmOpPtr input )	// fixme: apply the correct requirements for sigma_d
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}

	if ( _type == "Cs" )
		return false;

	for ( auto so : _symmetryOperations )
	{
		if ( !( so->getType() == "rotation" && so->getOrderOfRotation() == 2 ) )
			continue;
		for ( auto so2 : _symmetryOperations )
		{
			if ( !( so2->getType() == "rotation" && so2->getOrderOfRotation() == 2 ) || so == so2)
				continue;

			Vector3d vectorSum = (so->getAxis() + so2->getAxis()) / (so->getAxis() + so2->getAxis()).norm();
			if ( ( vectorSum.dot(input->getAxis()) )  < symmetryPrecision )
			{
					return true;
			}
		}
	}

	return false;
}

int PointGroup::getOrderOfMainRotationAxis()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	if ( _mainRotationAxis.size() == 0 )
	{
/*		cout << "no main rotation axis found" << endl;
		for ( auto so : _symmetryOperations )
		{
			if ( so->getType() == "rotation" )
				cout << "C" << so->getOrderOfRotation() << endl;
		}
*/
		return -1;
	}

	return _mainRotationAxis[0]->getOrderOfRotation();
}

SymmOpPtr PointGroup::getMainRotationAxis()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	return _mainRotationAxis[0];
}

vector<SymmOpPtr> PointGroup::getMainRotationAxes()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	return _mainRotationAxis;
}

size_t PointGroup::getPositionOfSymmetryOperation(SymmOpPtr so)
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	vector<vector<SymmOpPtr>> classes = _characterTable->getClasses();
	for ( size_t i = 0; i < classes.size(); ++i )
	{
		for ( auto so2 : classes[i] )
		if ( *so2 == *so )
		{
			if ( debugLevel > 0 )
				cout << so->getSymbol() << " == " << so2->getSymbol() << endl;

			return i;
		}
	}

	cout << _type << endl;
	cout << "ERROR: Symmetry operation " << so->getSymbol() << " not found!" << endl;
	cout << "ERROR: Axis is " << so->getAxis().transpose() << endl;
	cout << "ERROR: Transform matrix is" << endl << so->getTransformMatrix() << endl;
/*
	cout << "ERROR: Candidates are:" << endl;

	for ( size_t i = 0; i < classes.size(); ++i )
	{
		for ( auto so2 : classes[i] )
		{
			cout << "ERROR: " << so2->getSymbol() << " around " << so2->getAxis().transpose() << endl;
			cout << "ERROR: with transform matrix" << endl << so2->getTransformMatrix() << endl;
		}

	}
*/
	return 0;
}
