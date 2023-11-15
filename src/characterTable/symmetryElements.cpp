/*
 * symmetryElements.cpp
 *
 *  Created on: Feb 6, 2023
 *      Author: pemueller
 */

#include "symmetryElements.h"
#include "../common/constants.h"
#include <math.h>
#include <iostream>
#include <Eigen/Dense>

SymmetryOperation::SymmetryOperation()
{
	_type = "";
	_center = Vector3d::Zero();
	_axis = Vector3d::Zero();
	_orderOfRotation = -1;
	_error = 0;
	_numberOfApplications = 0;
	_transformMatrix = Matrix3d::Zero();
	_primes = 0;
}

SymmetryOperation::SymmetryOperation( string type, Vector3d center, Vector3d axis, int order, int numberOfApplications, int primes )
{
if (debugLevel > 3) {
	cout << "SymmetryOperation::" << __FUNCTION__ << endl;
}
	SymmetryOperation();
	_type = type;
	_center = (center);
	_error = 0.0;

	if ( _type == "identity" || numberOfApplications == 0 ) // todo: split into multiple functions
	{
		_primes = 0;
		_axis = Vector3d::Zero();
		_orderOfRotation = 0;
		_symbol = "E";
		_numberOfApplications = numberOfApplications;
		_transformMatrix = Vector3d::Ones().asDiagonal();
	}
	else if ( type == "inversion" )
	{
		_primes = 0;
		_axis = Vector3d::Zero();
		_orderOfRotation = 0;
		_symbol = "i";
		_numberOfApplications = numberOfApplications;
		_transformMatrix = -1 * Vector3d::Ones().asDiagonal();
	}
	else if (type == "rotation")
	{
		_orderOfRotation = order;
		_axis = (axis) / (axis).norm();
		string numberOfApplicationsAsString = to_string(numberOfApplications);

		if ( numberOfApplications < 0)
		{
			if (debugLevel > 1)
				cout << "Number of applications is below zero: " << numberOfApplications << endl <<
				"Fixing it to " << numberOfApplications + order << endl;
			numberOfApplicationsAsString = to_string( numberOfApplications + order );
		}

		_symbol = "C_"+to_string(order)+"^"+numberOfApplicationsAsString;
		_numberOfApplications = numberOfApplications;
		_primes = primes;
	    double x = _axis[0];
		double y = _axis[1];
		double z = _axis[2];

		double cos_a = cos(M_PI * _numberOfApplications * 2 / (double)_orderOfRotation);
		double sin_a = sin(M_PI * _numberOfApplications * 2 / (double)_orderOfRotation);

		_transformMatrix(0,0) = cos_a + x * x * (1 - cos_a);
		_transformMatrix(0,1) = -1.0 * z * sin_a + x * y * (1 - cos_a);
		_transformMatrix(0,2) = y * sin_a + x * z * (1 - cos_a);
		_transformMatrix(1,0) = z * sin_a + x * y * (1 - cos_a);
		_transformMatrix(1,1) = cos_a + y * y * (1 - cos_a);
		_transformMatrix(1,2) = -1.0 * x * sin_a + y * z * (1 - cos_a);
		_transformMatrix(2,0) = -1.0 * y * sin_a + z * x * (1 - cos_a);
		_transformMatrix(2,1) = x * sin_a + z * y * (1 - cos_a);
		_transformMatrix(2,2) = cos_a + z * z * (1 - cos_a);
	}
	else if ( type == "reflection" )
	{
		_primes = 0;
		_orderOfRotation = 0;
		_axis = (axis) / (axis).norm();
		_symbol = "sigma";
		_numberOfApplications = numberOfApplications;
		Matrix3d NN = 2 * _axis * _axis.transpose();
		Matrix3d ones = Vector3d::Ones().asDiagonal();
		_transformMatrix = ones - NN;
	}
	else if ( type == "rotationReflection" )
	{
		if ( debugLevel > 1 )
			cout << "adding a " << type << endl;
		_primes = 0;
		_orderOfRotation = order;
		_axis = (axis) / (axis).norm();
		string numberOfApplicationsAsString = to_string(numberOfApplications);
		if ( numberOfApplications < 0 && order % 2 == 0)
		{
			if (debugLevel > 1)
				cout << "Number of applications is below zero: " << numberOfApplications << endl <<
				"Fixing it to " << numberOfApplications + order << endl;
			numberOfApplicationsAsString = to_string( numberOfApplications + order );
		}

		_symbol = "S_"+to_string(order)+"^"+numberOfApplicationsAsString;
		_numberOfApplications = numberOfApplications;

		double x = _axis[0];
		double y = _axis[1];
		double z = _axis[2];

		double cos_a = cos(M_PI * _numberOfApplications * 2 / (double)_orderOfRotation);
		double sin_a = sin(M_PI * _numberOfApplications * 2 / (double)_orderOfRotation);
		Matrix3d rotationMatrix;
		rotationMatrix(0,0) = cos_a + x * x * (1 - cos_a);
		rotationMatrix(0,1) = -1.0 * z * sin_a + x * y * (1 - cos_a);
		rotationMatrix(0,2) = y * sin_a + x * z * (1 - cos_a);
		rotationMatrix(1,0) = z * sin_a + x * y * (1 - cos_a);
		rotationMatrix(1,1) = cos_a + y * y * (1 - cos_a);
		rotationMatrix(1,2) = -1.0 * x * sin_a + y * z * (1 - cos_a);
		rotationMatrix(2,0) = -1.0 * y * sin_a + z * x * (1 - cos_a);
		rotationMatrix(2,1) = x * sin_a + z * y * (1 - cos_a);
		rotationMatrix(2,2) = cos_a + z * z * (1 - cos_a);

		Matrix3d NN = 2 * _axis * _axis.transpose();
		Matrix3d ones = Vector3d::Ones().asDiagonal();
		for ( size_t i = 0; i < abs<int>(_numberOfApplications); ++i )
			rotationMatrix = (ones - NN) * rotationMatrix;
		_transformMatrix = rotationMatrix;
	}
	else
		throw( "Unknown symmetry operation "+type+". This is a bug. Please contact the author!" );
}

SymmetryOperation::SymmetryOperation( string type, SymmetryOperation other, int order, int numberOfApplications )
{
if (debugLevel > 3) {
	cout << "SymmetryOperation::" << __FUNCTION__ << endl;
}
	SymmetryOperation();
	_type = type;
	_center = other._center;
	_axis = other._axis;
	_numberOfApplications = other._numberOfApplications;
	_primes = other._primes;

	_error = 0.0;

	if ( _type == "identity" || numberOfApplications == 0 )
	{
		_symbol = "E";
		_numberOfApplications = numberOfApplications;
		_transformMatrix = Vector3d::Ones().asDiagonal();
	}
	else if ( type == "inversion" )
	{
		_orderOfRotation = 0;
		_symbol = "i";
		_numberOfApplications = numberOfApplications;
		_transformMatrix = -1 * Vector3d::Ones().asDiagonal();
	}
	else if (type == "rotation")
	{
		string numberOfApplicationsAsString = to_string(numberOfApplications);

		if ( numberOfApplications < 0)
		{
			if (debugLevel > 1)
				cout << "Number of applications is below zero: " << numberOfApplications << endl <<
				"Fixing it to " << numberOfApplications + order << endl;
			numberOfApplicationsAsString = to_string( numberOfApplications + order );
		}

		_symbol = "C_"+to_string(order)+"^"+numberOfApplicationsAsString;
	    double x = _axis[0];
		double y = _axis[1];
		double z = _axis[2];

		double cos_a = cos(M_PI * _numberOfApplications * 2 / (double)_orderOfRotation);
		double sin_a = sin(M_PI * _numberOfApplications * 2 / (double)_orderOfRotation);

		_transformMatrix(0,0) = cos_a + x * x * (1 - cos_a);
		_transformMatrix(0,1) = -1.0 * z * sin_a + x * y * (1 - cos_a);
		_transformMatrix(0,2) = y * sin_a + x * z * (1 - cos_a);
		_transformMatrix(1,0) = z * sin_a + x * y * (1 - cos_a);
		_transformMatrix(1,1) = cos_a + y * y * (1 - cos_a);
		_transformMatrix(1,2) = -1.0 * x * sin_a + y * z * (1 - cos_a);
		_transformMatrix(2,0) = -1.0 * y * sin_a + z * x * (1 - cos_a);
		_transformMatrix(2,1) = x * sin_a + z * y * (1 - cos_a);
		_transformMatrix(2,2) = cos_a + z * z * (1 - cos_a);
	}
	else if ( type == "reflection" )
	{
		_symbol = "sigma";
		_numberOfApplications = numberOfApplications;
		Matrix3d NN = 2 * _axis * _axis.transpose();
		Matrix3d ones = Vector3d::Ones().asDiagonal();
		_transformMatrix = ones - NN;
	}
	else if ( type == "rotationReflection" )
	{
		if ( debugLevel > 1 )
			cout << "adding a " << type << endl;
		string numberOfApplicationsAsString = to_string(numberOfApplications);
		if ( numberOfApplications < 0 && order % 2 == 0)
		{
			if (debugLevel > 1)
				cout << "Number of applications is below zero: " << numberOfApplications << endl <<
				"Fixing it to " << numberOfApplications + order << endl;
			numberOfApplicationsAsString = to_string( numberOfApplications + order );
		}

		_symbol = "S_"+to_string(order)+"^"+numberOfApplicationsAsString;

		double x = _axis[0];
		double y = _axis[1];
		double z = _axis[2];

		double cos_a = cos(M_PI * _numberOfApplications * 2 / (double)_orderOfRotation);
		double sin_a = sin(M_PI * _numberOfApplications * 2 / (double)_orderOfRotation);
		Matrix3d rotationMatrix;
		rotationMatrix(0,0) = cos_a + x * x * (1 - cos_a);
		rotationMatrix(0,1) = -1.0 * z * sin_a + x * y * (1 - cos_a);
		rotationMatrix(0,2) = y * sin_a + x * z * (1 - cos_a);
		rotationMatrix(1,0) = z * sin_a + x * y * (1 - cos_a);
		rotationMatrix(1,1) = cos_a + y * y * (1 - cos_a);
		rotationMatrix(1,2) = -1.0 * x * sin_a + y * z * (1 - cos_a);
		rotationMatrix(2,0) = -1.0 * y * sin_a + z * x * (1 - cos_a);
		rotationMatrix(2,1) = x * sin_a + z * y * (1 - cos_a);
		rotationMatrix(2,2) = cos_a + z * z * (1 - cos_a);

		Matrix3d NN = 2 * _axis * _axis.transpose();
		Matrix3d ones = Vector3d::Ones().asDiagonal();
		for ( size_t i = 0; i < abs<int>(_numberOfApplications); ++i )
			rotationMatrix = (ones - NN) * rotationMatrix;
		_transformMatrix = rotationMatrix;
	}
	else
		throw( "Unknown symmetry operation "+type+". This is a bug. Please contact the author!" );
}

SymmetryOperation::~SymmetryOperation(){}

string SymmetryOperation::getType()
{
if (debugLevel > 3) {
	cout << "SymmetryOperation::" << __FUNCTION__ << endl;
}
	return _type;
}

string SymmetryOperation::getSymbol()
{
if (debugLevel > 3) {
	cout << "SymmetryOperation::" << __FUNCTION__ << endl;
}
	return _symbol;
}

void SymmetryOperation::setSymbol(string symbol)
{
if (debugLevel > 3) {
	cout << "SymmetryOperation::" << __FUNCTION__ << endl;
}
	_symbol = symbol;
}

int SymmetryOperation::getOrderOfRotation()
{
if (debugLevel > 3) {
	cout << "SymmetryOperation::" << __FUNCTION__ << endl;
}
	return _orderOfRotation;
}

Vector3d SymmetryOperation::getAxis()
{
if (debugLevel > 3) {
	cout << "SymmetryOperation::" << __FUNCTION__ << endl;
}
	return _axis;
}

Vector3d SymmetryOperation::getCenter()
{
if (debugLevel > 3) {
	cout << "SymmetryOperation::" << __FUNCTION__ << endl;
}
	return _center;
}

int SymmetryOperation::getNumberOfApplications()
{
if (debugLevel > 3) {
	cout << "SymmetryOperation::" << __FUNCTION__ << endl;
}
	return _numberOfApplications;
}

Matrix3d SymmetryOperation::getTransformMatrix()
{
if (debugLevel > 3) {
	cout << "SymmetryOperation::" << __FUNCTION__ << endl;
}
	return _transformMatrix;
}

void SymmetryOperation::setPrimes( int primes )
{
if (debugLevel > 3) {
	cout << "SymmetryOperation::" << __FUNCTION__ << endl;
}
	_primes = primes;
}

int SymmetryOperation::getNumberOfPrimes()
{
if (debugLevel > 3) {
	cout << "SymmetryOperation::" << __FUNCTION__ << endl;
}
	return _primes;
}

Vector3d SymmetryOperation::applySymmetryOperation( Vector3d input)
{
if (debugLevel > 3) {
	cout << "SymmetryOperation::" << __FUNCTION__ << endl;
}
	return _transformMatrix * (input - _center) + _center;
}

bool SymmetryOperation::operator==(const SymmetryOperation other) const
{
if (debugLevel > 3) {
	cout << "SymmetryOperation::" << __FUNCTION__ << endl;
}

	return ((_transformMatrix - other._transformMatrix).cwiseAbs().maxCoeff() < 0.001 );
}

bool SymmetryOperation::operator>(const SymmetryOperation other) const // probably not needed anymore
{
	if (_type != other._type)
		return false;

	if ((_orderOfRotation % other._orderOfRotation) != 0)
		return false;

	 if (( (acos(_axis.dot(other._axis) ) / M_PI * 180) > 5) && ( (acos(_axis.dot(other._axis) ) / M_PI * 180) < 175 ))
		 return false;

	return true;
}

Matrix3d SymmetryOperation::operator-( const SymmetryOperation other ) const
{
	return _transformMatrix - other._transformMatrix;
}

SymmetryOperation SymmetryOperation::operator*(const SymmetryOperation other)
{
	SymmetryOperation result = other;

	result._transformMatrix = _transformMatrix * other._transformMatrix;

	result._symbol = _symbol+"*"+other._symbol;

	return result;
}

double SymmetryOperation::getError()
{
if (debugLevel > 3) {
	cout << "SymmetryOperation::" << __FUNCTION__ << endl;
}
	return _error;
}

void SymmetryOperation::setError( double error )
{
	_error = error;
}
