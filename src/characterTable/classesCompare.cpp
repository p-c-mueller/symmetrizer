/*
 * classesCompare.cpp
 *
 *  Created on: Jun 19, 2023
 *      Author: pemueller
 */

#include "characterTable.h"
#include "../common/constants.h"

#include <boost/make_shared.hpp>
#include <iostream>

bool CharacterTable::lt(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2)
{
	if ( debugLevel > 0 ) cout << "lt" << _type << endl;
	if ( _type == "Duh" )
		return ltDuh( class1, class2 );

	if ( _type == "Cuv" )
		return ltCuv( class1, class2 );

	if ( _type == "Ih" )
		return ltIh( class1, class2 );

	if ( _type == "Oh" )
		return ltOh( class1, class2 );

	if ( _type == "Td" )
		return ltTd( class1, class2 );

	if ( _type == "Dnh" )
		return ltDnh( class1, class2 );

	if ( _type == "Dnd" )
		return ltDnd( class1, class2 );

	if ( _type == "Dn" )
		return ltDn( class1, class2 );

	if ( _type == "Cnh" )
		return ltCnh( class1, class2 );

	if ( _type == "Cnv" )
		return ltCnv( class1, class2 );

	if ( _type == "Cn" )
		return ltCn( class1, class2 );

	if ( _type == "Sn" )
		return ltSn( class1, class2 );

	if ( _type == "Cs" )
		return ltCs( class1, class2 );

	if ( _type == "Ci" )
		return ltCi( class1, class2 );

	if ( _type == "C1" )
		return true;

//	throw (
	cout << "ERROR: Missing comparison operator for point group type " << _type << endl;
//	);
}

bool CharacterTable::ltDuh( vector<SymmOpPtr> class1, vector<SymmOpPtr> class2 )
{
	int a = -1;
	int b = -1;

	if ( class1[0]->getSymbol() == "E" ) a = 0;
	if ( class2[0]->getSymbol() == "E" ) b = 0;

	size_t i = 0;
	for ( int j = infinityRotationSteps; j > 0; --j)
	{
		if ( class1[0]->getOrderOfRotation() == j && class1[0]->getAxis().dot(_mainRotationAxis->getAxis()) > 0.999) a = i;
		if ( class2[0]->getOrderOfRotation() == j && class2[0]->getAxis().dot(_mainRotationAxis->getAxis()) > 0.999) b = i;
		++i;
	}

	if ( class2[0]->getSymbol() == "sigma_h" && b == -1) b = i;
	++i;


	if ( class1[0]->getType() == "reflection" && a == -1 ) a = i;
	if ( class2[0]->getType() == "reflection" && b == -1) b = i;
	++i;

	if ( class1[0]->getType() == "inversion" ) a = i;
	if ( class2[0]->getType() == "inversion" ) b = i;
	++i;

	// Sn
	for ( int j = 2 * infinityRotationSteps; j > 0; --j )
	{
		if ( class1[0]->getOrderOfRotation() == j && class1[0]->getType() == "rotationReflection" ) a = i;
		if ( class2[0]->getOrderOfRotation() == j && class2[0]->getType() == "rotationReflection" ) b = i;

		++i;
	}

	if ( class1[0]->getType() == "rotation" && a == -1 ) a = i;
	if ( class2[0]->getType() == "rotation" && b == -1 ) b = i;
	++i;

	if ( debugLevel > 0 )
		cout << class1[0]->getSymbol() << ": " << a << " " << class2[0]->getSymbol() << ": " << b << endl;

	return a < b;
}

bool CharacterTable::ltCuv(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2)
{
	return ltDuh(class1, class2);

	/*// todo
	if ( class1[0]->getSymbol() == "E" )
		return true;
	if ( class2[0]->getSymbol() == "E" )
		return false;

	return false;*/
}

bool CharacterTable::ltIh(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2)
{
	// todo
	if ( class1[0]->getSymbol() == "E" )
		return true;
	if ( class2[0]->getSymbol() == "E" )
		return false;

	return false;
}

//#include <boost/algorithm/>

bool CharacterTable::ltOh(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2)
{
	int a = -1;
	int b = -1;

	if ( class1[0]->getSymbol() == "E" ) a = 0;
	if ( class2[0]->getSymbol() == "E" ) b = 0;

	if ( class1[0]->getSymbol() == "C3^1" ) a = 1;
	if ( class2[0]->getSymbol() == "C3^1" ) b = 1;

	if ( class1[0]->getSymbol() == "C2^1" && class1.size() == 6 ) a = 2;
	if ( class2[0]->getSymbol() == "C2^1" && class2.size() == 6 ) b = 2;

	if ( class1[0]->getSymbol() == "C4^1" ) a = 3;
	if ( class2[0]->getSymbol() == "C4^1" ) b = 3;

	if ( class1[0]->getSymbol() == "C2^1" && class1.size() == 3 ) a = 4;
	if ( class2[0]->getSymbol() == "C2^1" && class2.size() == 3 ) b = 4;

	if ( class1[0]->getSymbol() == "i" ) a = 5;
	if ( class2[0]->getSymbol() == "i" ) b = 5;

	if ( class1[0]->getSymbol() == "S4^1" ) a = 6;
	if ( class2[0]->getSymbol() == "S4^1" ) b = 6;

	if ( class1[0]->getSymbol() == "S6^1" ) a = 7;
	if ( class2[0]->getSymbol() == "S6^1" ) b = 7;

	if ( class1[0]->getSymbol() == "sigma_h" ) a = 8;
	if ( class2[0]->getSymbol() == "sigma_h" ) b = 8;

	if ( class1[0]->getSymbol() == "sigma_d" ) a = 9;
	if ( class2[0]->getSymbol() == "sigma_d" ) b = 9;

	if ( debugLevel > 0 )
		cout << class1[0]->getSymbol() << ": " << a << " " << class2[0]->getSymbol() << ": " << b << endl;


	return a < b;
}

bool CharacterTable::ltTd(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2)
{
	int a = -1;
	int b = -1;

	if ( class1[0]->getSymbol() == "E" ) a = 0;
	if ( class2[0]->getSymbol() == "E" ) b = 0;

	if ( class1[0]->getSymbol() == "C3^1" ) a = 1;
	if ( class2[0]->getSymbol() == "C3^1" ) b = 1;

	if ( class1[0]->getSymbol() == "C2^1" ) a = 2;
	if ( class2[0]->getSymbol() == "C2^1" ) b = 2;

	if ( class1[0]->getSymbol() == "S4^1" ) a = 3;
	if ( class2[0]->getSymbol() == "S4^1" ) b = 3;

	if ( class1[0]->getSymbol() == "sigma_d" ) a = 4;
	if ( class2[0]->getSymbol() == "sigma_d" ) b = 4;

	return a < b;
}

bool CharacterTable::ltDnh(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2)
{
	int a = -1;
	int b = -1;

	// n odd
	// first Cn; then Cn * sigma_h
	if ( _mainRotationAxis->getOrderOfRotation() % 2 == 1 )
	{
		SymmetryOperation Cn ("identity");

		int i = 0;
		while ( i < _mainRotationAxis->getOrderOfRotation() )
		{
			if ( *class1[0] == Cn ) a = i;
			if ( *class2[0] == Cn ) b = i;

			Cn = Cn * *_mainRotationAxis;

			++i;
		}

		if ( class1[0]->getSymbol() == "C2^1" && class1.size() == _mainRotationAxis->getOrderOfRotation() ) a = i;
		if ( class2[0]->getSymbol() == "C2^1" && class2.size() == _mainRotationAxis->getOrderOfRotation() ) b = i;
		++i;

		SymmetryOperation sigma_h ("reflection", _mainRotationAxis->getCenter(), _mainRotationAxis->getAxis());
		Cn = SymmetryOperation("identity");

		for ( int j = 0; j < _mainRotationAxis->getOrderOfRotation(); ++j )
		{
			if ( *class1[0] == (Cn * sigma_h )) a = i;
			if ( *class2[0] == (Cn * sigma_h )) b = i;

			Cn = Cn * *_mainRotationAxis;

			++i;
		}

		if ( class1[0]->getType() == "reflection" && class1.size() == _mainRotationAxis->getOrderOfRotation() ) a = i;
		if ( class2[0]->getType() == "reflection" && class2.size() == _mainRotationAxis->getOrderOfRotation() ) b = i;
		++i;
	}
	// n even
	// first Cn; then Cn * i
	else
	{
		SymmetryOperation Cn ("identity");
		int i = 0;
		while ( i < _mainRotationAxis->getOrderOfRotation() )
		{
			if ( *class1[0] == Cn && a == -1 ) a = i;
			if ( *class2[0] == Cn && b == -1 ) b = i;

			Cn = Cn * *_mainRotationAxis;

			++i;
		}

		if ( class1[0]->getSymbol() == "C2^1" && class1[0]->getNumberOfPrimes() == 1 ) a = i;
		if ( class2[0]->getSymbol() == "C2^1" && class2[0]->getNumberOfPrimes() == 1 ) b = i;
		++i;

		if ( class1[0]->getSymbol() == "C2^1" && class1[0]->getNumberOfPrimes() == 2 ) a = i;
		if ( class2[0]->getSymbol() == "C2^1" && class2[0]->getNumberOfPrimes() == 2 ) b = i;
		++i;

		// Cn * i
		Cn = SymmetryOperation("identity");
		for ( int j = 0; j < _mainRotationAxis->getOrderOfRotation(); ++j )
		{
			if ( *class1[0] == (Cn * SymmetryOperation("inversion") ) && a == -1) a = i;
			if ( *class2[0] == (Cn * SymmetryOperation("inversion") ) && b == -1) b = i;

			Cn = Cn * *_mainRotationAxis;

			++i;
		}

		if ( class1[0]->getType() == "reflection" && class1[0]->getNumberOfPrimes() == 1 ) a = i;
		if ( class2[0]->getType() == "reflection" && class2[0]->getNumberOfPrimes() == 1 ) b = i;
		++i;

		if ( class1[0]->getType() == "reflection" && class1[0]->getNumberOfPrimes() == 2 ) a = i;
		if ( class2[0]->getType() == "reflection" && class2[0]->getNumberOfPrimes() == 2 ) b = i;
		++i;
	}

	if ( debugLevel > 0 ) cout << "Class 1: " << class1[0]->getSymbol() << "_" << class1[0]->getNumberOfPrimes() << " == " << a << endl <<
								"Class 2: " << class2[0]->getSymbol() << "_" <<class2[0]->getNumberOfPrimes() << " == " << b << endl;

	return a < b;
}

bool CharacterTable::ltDnd(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2)
{
	int a = -1;
	int b = -1;

	// n odd
	if ( _mainRotationAxis->getOrderOfRotation() % 2 == 1 )
	{
		// Cn -> Cn^k; nC2
		SymmetryOperation Cn ("identity");

		int i = 0;
		while ( i < _mainRotationAxis->getOrderOfRotation() )
		{
			if ( *class1[0] == Cn ) a = i;
			if ( *class2[0] == Cn ) b = i;

			Cn = Cn * *_mainRotationAxis;

			++i;
		}

		if ( class1[0]->getSymbol() == "C2^1" ) a = i;
		if ( class2[0]->getSymbol() == "C2^1" ) b = i;
		++i;

		// Cn * i
		Cn = SymmetryOperation("identity");
		for ( int j = 0; j < _mainRotationAxis->getOrderOfRotation(); ++j )
		{
			if ( *class1[0] == (Cn * SymmetryOperation("inversion") )) a = i;
			if ( *class2[0] == (Cn* SymmetryOperation("inversion") )) b = i;

			Cn = Cn * *_mainRotationAxis;

			++i;
		}
		if ( class1[0]->getSymbol() == "sigma_d" ) a = i;
		if ( class2[0]->getSymbol() == "sigma_d" ) b = i;
	}
	else // n even
	{
		// Sn^k; C2'; sigma
		SymmetryOperation Cn("identity");
		int i = 0;
		while( i < _mainRotationAxis->getOrderOfRotation() )
		{
			if ( *class1[0] == Cn ) a = i;
			if ( *class2[0] == Cn ) b = i;

			Cn = Cn * *_mainRotationAxis;
			++i;
		}

		if ( class1[0]->getSymbol() == "C2^1" && class1.size() == _mainRotationAxis->getOrderOfRotation()/2 ) a = i;
		if ( class2[0]->getSymbol() == "C2^1" && class2.size() == _mainRotationAxis->getOrderOfRotation()/2 ) b = i;
		++i;

		if ( class1[0]->getType() == "reflection" ) a = i;
		if ( class2[0]->getType() == "reflection" ) b = i;
	}

	if ( debugLevel > 0 )
		cout << class1[0]->getSymbol() << ": " << a << " " << class2[0]->getSymbol() << ": " << b << endl;

	return a < b;
}

bool CharacterTable::ltDn(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2)
{
	// fixme D2 not working
	// C_n^k with k = 0 < n; then C_2
	int a = -1;
	int b = -1;
	size_t i = 0;

	SymmetryOperation Cn ("identity");

	while ( i < Cn.getOrderOfRotation())
	{
		if ( isIn(boost::make_shared<SymmetryOperation>(Cn),class1))
			a = i;
		if ( isIn(boost::make_shared<SymmetryOperation>(Cn),class2))
			b = i;
		Cn = Cn * *_mainRotationAxis;

		++i;
	}

	if ( class1[0]->getType()== "rotation" && class1[0]->getAxis() != _mainRotationAxis->getAxis() && class1[0]->getNumberOfPrimes() == 1 )
		a = i;
	if ( class2[0]->getType()== "rotation" && class1[0]->getAxis() != _mainRotationAxis->getAxis() && class2[0]->getNumberOfPrimes() == 1 )
		b = i;

	++i;

	if ( class1[0]->getType()== "rotation" && class1[0]->getAxis() != _mainRotationAxis->getAxis() && class1[0]->getNumberOfPrimes() == 2 )
		a = i;
	if ( class2[0]->getType()== "rotation" && class1[0]->getAxis() != _mainRotationAxis->getAxis() && class2[0]->getNumberOfPrimes() == 2 )
		b = i;
	return a < b;
}

bool CharacterTable::ltCnh(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2)
{
	int a = -1;
	int b = -1;

	// n odd
	// C_n^k with k = 0 < n; then sigma_h then S_n^k
	if ( _mainRotationAxis->getOrderOfRotation() % 2 == 1 )
	{
		SymmetryOperation Cn ("identity");
		size_t i = 0;
		while ( i < _mainRotationAxis->getOrderOfRotation())
		{
			if ( isIn(boost::make_shared<SymmetryOperation>(Cn),class1))
				a = i;
			if ( isIn(boost::make_shared<SymmetryOperation>(Cn),class2))
				b = i;
			Cn = Cn * *_mainRotationAxis;

			++i;
		}

		if ( class1[0]->getType() == "reflection" )
			a = i;
		if ( class2[0]->getType() == "reflection" )
			b = i;
		++i;

		const SymmetryOperation Sn ( "rotationReflection", _mainRotationAxis->getCenter(), _mainRotationAxis->getAxis(), _mainRotationAxis->getOrderOfRotation() );

		SymmetryOperation so = Sn;
		size_t j = 0;
		while ( j < 2 * Cn.getOrderOfRotation())
		{
			if ( isIn(boost::make_shared<SymmetryOperation>(so),class1) && a == -1)
				a = i;
			if ( isIn(boost::make_shared<SymmetryOperation>(so),class2) && b == -1)
				b = i;
			so = so * Sn;
			++j;
			++i;
		}
	}
	// n even
	// C_n^k with k = 0 < n; then i * C_n^k with k = 0 < n
	else
	{
		SymmetryOperation Cn ("identity");

		if ( debugLevel > 0)
			cout << "ltCnh even " << Cn.getOrderOfRotation() << endl;

		size_t i = 0;
		while ( i < _mainRotationAxis->getOrderOfRotation()  )
		{
			if ( debugLevel > 0)
				cout << i << endl;

			if ( isIn(boost::make_shared<SymmetryOperation>(Cn),class1))
				a = i;
			if ( isIn(boost::make_shared<SymmetryOperation>(Cn),class2))
				b = i;
			Cn = Cn * *_mainRotationAxis;
			++i;
		}

		Cn = SymmetryOperation("identity");
		SymmetryOperation inversion("inversion");
		i = 0;
		while ( i < _mainRotationAxis->getOrderOfRotation() )
		{
			if ( debugLevel > 0)
				cout << i + _mainRotationAxis->getOrderOfRotation() << endl;

			SymmetryOperation CnTimesI = Cn *  inversion;
			if ( isIn(boost::make_shared<SymmetryOperation>(CnTimesI),class1))
				a = i + _mainRotationAxis->getOrderOfRotation();
			if ( isIn(boost::make_shared<SymmetryOperation>(CnTimesI),class2))
				b = i + _mainRotationAxis->getOrderOfRotation();
			Cn = Cn * *_mainRotationAxis;
			++i;
		}
	}
	if ( debugLevel > 1)
		cout << a << " " << b << endl;

	return a < b;
}

bool CharacterTable::ltCnv(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2)
{
	// C_n^k with k = 0 < n/2; then sigma_v, then sigma_d
	int a = -1;
	int b = -1;

	SymmetryOperation Cn ("identity");
	size_t i = 0;
	while ( i < Cn.getOrderOfRotation() )
	{
		if (isIn(boost::make_shared<SymmetryOperation>(Cn),class1))
			a = i;
		if (isIn(boost::make_shared<SymmetryOperation>(Cn),class2))
			b = i;
		Cn = Cn * *_mainRotationAxis;
		 ++i;
	}

	// fixme: primes do not seem to be read in correctly
	if ( class1[0]->getType()== "reflection" && class1[0]->getNumberOfPrimes() == 1 )
		a = i;
	if ( class2[0]->getType()== "reflection" && class2[0]->getNumberOfPrimes() == 1 )
		b = i;

	++i;

	if ( class1[0]->getType()== "reflection" && class1[0]->getNumberOfPrimes() == 2 )
		a = i;
	if ( class2[0]->getType()== "reflection" && class2[0]->getNumberOfPrimes() == 2 )
		b = i;

	return a < b;
}

bool CharacterTable::ltCn(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2)
{
	// C_n^k with k = 0 < n
	// E is always first.
	if ( class1[0]->getSymbol() == "E" )
		return true;
	if ( class2[0]->getSymbol() == "E" )
		return false;

	return ((double)class1[0]->getOrderOfRotation()/(double)class1[0]->getNumberOfApplications() > (double)class2[0]->getOrderOfRotation()/(double)class2[0]->getNumberOfApplications());
}

bool CharacterTable::ltSn(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2)
{
	// S_n^k with k = 0 < n
	int a = -1;
	int b = -1;

	SymmetryOperation Sn("identity");

	if ( debugLevel > 1 )
		cout << endl << _mainRotationAxis->getSymbol() << endl;

	for ( size_t i = 0; i < _mainRotationAxis->getOrderOfRotation(); ++i )
	{
		if (debugLevel > 1)
			cout << "run " << i+1 << " of " << _mainRotationAxis->getOrderOfRotation() << endl;

		if ( isIn(boost::make_shared<SymmetryOperation>(Sn),class1))
			a = i;
		if ( isIn(boost::make_shared<SymmetryOperation>(Sn),class2))
			b = i;
		Sn = Sn * *_mainRotationAxis;
	}

	if ( debugLevel > 1)
		cout << "a=" << a << " b= " << b << endl;

	return a < b;
}

bool CharacterTable::ltCs(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2)
{
	// C_n^k with k = 0 < n
	// E is always first.
	if ( class1[0]->getSymbol() == "E" )
		return true;
	return false;
}

bool CharacterTable::ltCi(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2)
{
	// E is always first.
	if ( class1[0]->getSymbol() == "E" )
		return true;

	return false;
}
