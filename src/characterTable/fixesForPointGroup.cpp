/*
 * fixesForPointGroup.cpp
 *
 *  Created on: Jun 19, 2023
 *      Author: pemueller
 */


#include "pointGroup.h"

#include "../testing/pgTests.h"

void PointGroup::setPrimes(vector<SymmOpPtr> soVector, int primes)
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	for ( SymmOpPtr so : soVector )
	{
		so->setPrimes(primes);
	}
}

void PointGroup::setPrimesCnv()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
if ( debugLevel > 0 ) cout << "Setting primes for " << _type << endl;
	vector<vector<SymmOpPtr>> so;

	vector<VectorXd> atomDistances;
	// todo exclude mra
	// todo make break points if differentiation of classes is possible
	// todo update class labels according to importance/weight

	for ( auto soClass : _characterTable->getClasses() )
	{
		if ( soClass[0]->getType() != "reflection" ) // only sigma_v/sigma_d can occur sigma_v = '; sigma_d = ''
			continue;

		so.push_back(soClass);
	}

	if ( isGreater(so[0], so[1]) )
	{
		setPrimes(so[0], 2);
		setPrimes(so[1], 1);
	}
	else
	{
		setPrimes(so[0], 1);
		setPrimes(so[1], 2);
	}
}

void PointGroup::setPrimesDn()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	if ( debugLevel > 0 ) cout << "Setting primes for " << _type << endl;

	vector<vector<SymmOpPtr>> so;

	vector<VectorXd> atomDistances;
	// todo exclude mra
	// todo make break points if differentiation of classes is possible
	// todo update class labels according to importance/weight

	for ( auto soClass : _characterTable->getClasses() )
	{
		if ( soClass[0]->getType() != "rotation" || soClass[0]->getOrderOfRotation() != 2 ) // only C2 axes
			continue;

		so.push_back(soClass);
	}

	VectorXi primes = VectorXi::Zero(so.size());

	if (debugLevel > 0)
		cout << "found " << primes.size() << " C2 axes" << endl;

	for ( size_t i = 0; i < so.size(); ++i )
	{
		for ( size_t j = i; j < so.size(); ++j )
		{
			if ( isGreater(so[i], so[j]) )
				++primes(i);
			else
				++primes(j);
		}
	}

	primes = primes - VectorXi::Ones(primes.size());

	if ( debugLevel > 0 )
		cout << "no of primes are " << primes.transpose() << endl;

	for ( size_t i = 0; i < so.size(); ++i )
		setPrimes(so[i], primes(i));

	if ( debugLevel > 0 )
		for ( size_t i = 0; i < so.size(); ++i )
			cout << "primes of " << i + 1 << ": " << so[i][0]->getNumberOfPrimes() << endl;

	if ( _mainRotationAxis[0]->getOrderOfRotation() == 2 )
		fixD2();
}

void PointGroup::setPrimesDndh()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}

 	if ( debugLevel > 0 ) cout << "Setting primes for " << _type << endl;

	vector<vector<SymmOpPtr>> so;

	vector<VectorXd> atomDistances;
	// todo exclude mra
	// todo make break points if differentiation of classes is possible
	// todo update class labels according to importance/weight

	for ( auto soClass : _characterTable->getClasses() ) // C2 axes and sigma_d/sigma_v; just take one of them and apply the same to the other one.
	{
		if ( soClass[0]->getType() != "reflection" ) // only sigma_v/sigma_d can occur sigma_v = '; sigma_d = ''
			continue;

		so.push_back(soClass);
	}

	VectorXi primes = VectorXi::Zero(so.size());

/*	if ( _type == "Dnh" )
		primes = -1 * primes; // the sigma_h is supposed to have 0 primes.
	else
		primes.setZero();*/

	for ( size_t i = 0; i < so.size(); ++i )
	{
		for ( size_t j = i; j < so.size(); ++j )
		{
			if ( isGreater(so[i], so[j]) )
				++primes(i);
			else
				++primes(j);
		}
	}

	for ( int i = 0; i < primes.size(); ++i )
		if ( primes(i) > 2 )
			primes(i) = 0;


	for ( size_t i = 0; i < so.size(); ++i )
		setPrimes(so[i], primes(i));

	if ( debugLevel > 0 ) {cout << "Primes are " << primes.transpose() << endl;
		for ( auto it : so )
			cout << it[0]->getNumberOfPrimes() << " ";
		cout << endl;
	}

	// set the primes for the C2 axes according to the sigmas.

	for ( vector<SymmOpPtr> soClass : _characterTable->getClasses() )	// soClass contains the C2 axes
	{
		if ( soClass[0]->getType() != "rotation" || soClass[0]->getAxis() == _mainRotationAxis[0]->getAxis() )
			continue;

		for ( vector<SymmOpPtr> soClass2 : so )	// soClass2 contains the sigma
		{
			for ( SymmOpPtr c2 : soClass )
			{
				for ( SymmOpPtr sigma : soClass2 )
				{
					if ( abs(c2->getAxis().dot(sigma->getAxis())) > 0.99 )
					{
						c2->setPrimes(sigma->getNumberOfPrimes());
					}
				}
			}
		}
	}

	if ( _type == "Dnd" && _mainRotationAxis[0]->getOrderOfRotation() == 2 )
		fixD2d();

	if ( _type == "Dnh" && _mainRotationAxis[0]->getOrderOfRotation() == 2 )
		fixD2h();
}

void PointGroup::fixD2()
{
	// In the D2 point group, there are three types of C2 axes. The one with no primes shall be the main rotation axis.
	_mainRotationAxis.clear();

	for ( auto it : _symmetryOperations )
	{
		if ( it->getOrderOfRotation() == 2 && it->getNumberOfPrimes() == 0 )
			_mainRotationAxis.push_back(it);
	}
}

void PointGroup::fixD2d()
{
	// In the D2d point group, there are two types of C2 axes. One class contains the main rotation axis, the other one contains two axes.
	for ( size_t i = 0; i < _mainRotationAxis.size(); ++i )
	{
		for ( size_t j = 0; j < _characterTable->getClasses().size(); ++j)
		{
			if ( _mainRotationAxis[i] == _characterTable->getClasses()[j][0] && _characterTable->getClasses()[j].size() == 1 )
			{
				SymmOpPtr mra = _mainRotationAxis[i];
				_mainRotationAxis.clear();
				_mainRotationAxis.push_back(mra);
			}
		}
	}
}

void PointGroup::fixD2h()
{
	// The D2h point group contains three C2 axes, each in their own class.
	// The main rotation axis shall be perpenticular to the sigma_h

	vector<vector<SymmOpPtr>> so;

	for ( auto soClass : _characterTable->getClasses() )
	{
		if ( soClass[0]->getType() != "reflection" ) // only sigma_v/sigma_d can occur sigma_v = '; sigma_d = ''
			continue;

		so.push_back(soClass);
	}

	//cout << "Fixing Dh2 having " << so.size() << " different C2 axes" << endl;

	if ( isGreater(so[0], so[1]) )
	{
		setPrimes(so[0], 2);
		setPrimes(so[1], 1);
	}
	else
	{
		setPrimes(so[0], 1);
		setPrimes(so[1], 2);
	}

}

void PointGroup::setS2nAsMainRotationAxis()
{
	_mainRotationAxis.clear();

	for ( auto it : _symmetryOperations )
	{
		if ( it->getType() == "rotationReflection" && it->getNumberOfApplications() == 1  )
		{
			_mainRotationAxis.push_back(it);
			return;
		}
	}
}
