/*
 * decisionTree.cpp
 *
 *  Created on: Jun 19, 2023
 *      Author: pemueller
 */

#include "pointGroup.h"
#include <iostream>

#include "../testing/pgTests.h"

void PointGroup::findNameOfPointGroup()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	if ( isLinear() && suggestedPointGroup == "" )
	{
		if ( hasInversion() )
		{
			_name = "Duh";
			_type = "Duh";
		}
		else
		{
			_name = "Cuv";
			_type = "Cuv";
		}
	}
	else
	{
		if ( hasTwoOrMoreCn() )
		{
			if ( hasInversion() )
			{
				if ( hasC5()  )
				{
					_name = "Ih";
					_type = "Ih";
				}
				else
				{
					_name = "Oh";
					_type = "Oh";
				}
			}
			else
			{
				_name = "Td";
				_type = "Td";
			}
		}
		else
		{
			if ( hasCn() )
			{
				if ( isDihedral() )
				{
					if ( hasSigmaHDihedral() )
					{
						_name = "D"+to_string(_mainRotationAxis[0]->getOrderOfRotation())+"h";
						_type = "Dnh";
					}
					else
					{
						if ( hasNSigmaD() )
						{
							_name = "D"+to_string(_mainRotationAxis[0]->getOrderOfRotation())+"d";
							_type = "Dnd";
							if( _mainRotationAxis[0]->getOrderOfRotation() % 2 == 0 )
								setS2nAsMainRotationAxis();
						}
						else
						{
							_name = "D"+to_string(_mainRotationAxis[0]->getOrderOfRotation());
							_type = "Dn";
						}
					}
				}
				else
				{
					if ( hasSigmaH() )
					{
						_name = "C"+to_string(_mainRotationAxis[0]->getOrderOfRotation())+"h";
						_type = "Cnh";
					}
					else
					{
						if ( hasNSigmaV() )
						{
							_name = "C"+to_string(_mainRotationAxis[0]->getOrderOfRotation())+"v";
							_type = "Cnv";
						}
						else
						{
							if ( hasS2N() )
							{
								_name = "S"+to_string(2 * _mainRotationAxis[0]->getOrderOfRotation())+"";
								_type = "Sn";
								setS2nAsMainRotationAxis();
							}
							else
							{
								_name = "C"+to_string(_mainRotationAxis[0]->getOrderOfRotation());
								_type = "Cn";
							}
						}
					}
				}
			}
			else
			{
				if ( hasSigma() )
				{
					_name = "Cs";
					_type = "Cs";
				}
				else
				{
					if ( hasInversion() )
					{
						_name = "Ci";
						_type = "Ci";
					}
					else
					{
						_name = "C1";
						_type = "C1";
					}

				}
			}
		}
	}
}

bool PointGroup::isLinear()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}

	return _isLinear;
}

bool PointGroup::hasInversion()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	for ( auto so : _symmetryOperations )
	{
		if ( ( so->getType() == "inversion") )
			return true;
	}

	return false;
}

bool PointGroup::hasTwoOrMoreCn()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	if ( _mainRotationAxis.size() > 1 && _mainRotationAxis[0]->getOrderOfRotation() > 2 )
			return true;

	return false;
}

bool PointGroup::hasC5()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	for ( auto so : _symmetryOperations )
	{
		if ( (so->getType() == "rotation") && ( so->getOrderOfRotation() == 5 ))
			return true;
	}

	return false;
}

bool PointGroup::hasCn()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	for ( auto so : _symmetryOperations )
	{
		if ( (so->getType() == "rotation") )
			return true;
	}

	return false;
}

bool PointGroup::hasMultipleC2()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	int counter = 0;

	for ( auto so : _symmetryOperations )
	{
		if ( (so->getType() == "rotation") && ( so->getOrderOfRotation() == 2 ))
		{
			if ( counter == 0 )
				++counter;
			else
				return true;
		}
	}

	return false;
}

bool PointGroup::hasSigma()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	for ( auto so : _symmetryOperations )
	{
		if ( (so->getType() == "reflection") )
			return true;
	}

	return false;
}

bool PointGroup::isDihedral()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	int it = 0;
	for ( auto so : _symmetryOperations )
	{
		if ( so-> getOrderOfRotation() != 2 )
			continue;

		Vector3d v1 = so->getAxis();
		Vector3d v2 = _mainRotationAxis[0]->getAxis();

		if ( abs( v1.dot(v2) < 0.0001) )
			++it;
	}


	if ( it == _mainRotationAxis[0]->getOrderOfRotation()  )
		return true;

	if ( debugLevel > 0)
		cout << "Point group is not dihedral. Found " << it << " reflection planes but expected " << _mainRotationAxis[0]->getOrderOfRotation() << endl;

	return false;
}

bool PointGroup::hasSigmaH()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	for ( auto so : _symmetryOperations )
	{
		Vector3d v1 = so->getAxis();
		Vector3d v2 = _mainRotationAxis[0]->getAxis();
		if ( so->getType() == "reflection"
				&& ( abs(v1.dot(v2)) > 0.9999 ))
			return true;
	}

	return false;
}

bool PointGroup::hasSigmaHDihedral()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	int it = 0;

	for ( auto so : _symmetryOperations )
	{
		if ( so->getType() == "reflection" )
			++it;
	}

	return (it == _mainRotationAxis[0]->getOrderOfRotation() + 1);
}

bool PointGroup::hasNSigmaD()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	int it = 0;
	for ( auto so : _symmetryOperations )
	{

		Vector3d v1 = so->getAxis();
		Vector3d v2 = _mainRotationAxis[0]->getAxis();
		if (
				so->getType() == "reflection"
				&& (( (acos( (v1.dot(v2) / ( v1.norm() * v2.norm() )) ) / M_PI * 180) > 85)
					|| ( (acos( (v1.dot(v2) / ( v1.norm() * v2.norm() )) ) / M_PI * 180) < 95 ))
			)
			++it;
	}

	if ( it == _mainRotationAxis[0]->getOrderOfRotation() )
		return true;
	else
		return false;
}

bool PointGroup::hasNSigmaV()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	int it = 0;
	for ( auto so : _symmetryOperations )
	{
		Vector3d v1 = so->getAxis();
		Vector3d v2 = _mainRotationAxis[0]->getAxis();
		if ( so->getType() == "reflection" && ( abs(v1.dot(v2)) < 0.0001 ) )
			++it;
	}

	return ( it == _mainRotationAxis[0]->getOrderOfRotation() );
}

bool PointGroup::hasS2N()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	for ( auto it : _symmetryOperations )
	{
		if ( it->getOrderOfRotation() == 2* _mainRotationAxis[0]->getOrderOfRotation() && it->getType() == "rotationReflection")
			return true;
	}

	return false;
}


