/*
 * characterTable.cpp
 *
 *  Created on: Feb 23, 2023
 *      Author: pemueller
 */

#include "characterTable.h"
#include "../common/constants.h"
#include <boost/make_shared.hpp>
#include <iostream>
#include <algorithm>
#include "../IO/writer.h"
#include "../IO/reader.h"

int numberOfRuns = 1;

CharacterTable::CharacterTable(string type, SymmOpPtr mra)
{
	_classes.clear();
	_irreps.clear();
	_order = 0;
	_type = type;
	_mainRotationAxis = mra;
	_irreps.clear();
}

CharacterTable::~CharacterTable(){}

void CharacterTable::clear()
{
	_classes.clear();
	_irreps.clear();
	_order = 0;
	_type = "";
	_mainRotationAxis = boost::make_shared<SymmetryOperation>("identity");
}

void CharacterTable::addClass(vector<SymmOpPtr> input)
{
if (debugLevel > 3) {
	cout << "CharacterTable::" << __FUNCTION__ << endl;
	cout << numberOfRuns << " runs" << endl;
	++numberOfRuns;
}
	_classes.push_back(input);
}

vector<vector<SymmOpPtr>> CharacterTable::getClasses()
{
if (debugLevel > 3) {
	cout << "CharacterTable::" << __FUNCTION__ << endl;
}
	return _classes;
}

void CharacterTable::combineC2AndSigmaV()
{
	if ( debugLevel > 0 ) cout << "combining C2 and reflections for linear molecule" << endl;

	for ( size_t j = 0; j < _classes.size(); ++j )
	{
		if ( (_classes[j][0]->getType() == "rotation" || _classes[j][0]->getType() == "reflection" )&& _classes[j][0]->getAxis().dot(_mainRotationAxis->getAxis()) < 0.001 )
		{
			for ( size_t i = j+1; i < _classes.size(); ++i )
			{
				if ( _classes[i][0]->getType() == _classes[j][0]->getType() && _classes[i][0]->getAxis().dot(_mainRotationAxis->getAxis()) < 0.001 )
				{
					_classes[j].insert(_classes[j].end(),_classes[i].begin(),_classes[i].end());
					_classes.erase(_classes.begin()+i);
				}
			}
		}
	}
}

void CharacterTable::addIrrep( VectorXcd irrep )
{
if (debugLevel > 3) {
	cout << "CharacterTable::" << __FUNCTION__ << endl;
}
	if ( _type == "Duh" || _type == "Cuv" ) // fallback for linear molecules. The below function does not work there.
	{
		_irreps.push_back(boost::make_shared<Irrep>(irrep));
		return;
	}


	if ( !irrepAlreadyExists(irrep) )
	{
		// if the irrep is degenerate, add it to the existing one, if we have any

		for ( size_t i = 0; i < _irreps.size(); ++i )
		{
			auto irrepVector = _irreps[i];
			for ( auto ir : irrepVector->getCharacters() )
			{
				if ( !irrepsOrthogonal(ir, irrep) )
				{
					_irreps[i]->add(irrep);
					return;
				}
			}
		}

		// add a new irrep, if
		_irreps.push_back(boost::make_shared<Irrep>(irrep));
	}
}

vector<boost::shared_ptr<Irrep>> CharacterTable::getIrreps()
{
if (debugLevel > 3) {
	cout << "CharacterTable::" << __FUNCTION__ << endl;
}
	return _irreps;
}

bool CharacterTable::irrepAlreadyExists( VectorXcd irrep )
{
if (debugLevel > 3) {
	cout << "CharacterTable::" << __FUNCTION__ << endl;
	cout << "irrep = " << irrep.transpose() << endl;
}
	for ( auto irrepVector : _irreps )
	{
		for ( auto ir : irrepVector->getCharacters() )
		{
			if ( (ir  - irrep).cwiseAbs().maxCoeff() < symmetryPrecision )
				return true;
		}
	}

	return false;
}

bool CharacterTable::irrepsOrthogonal(VectorXcd irrep1, VectorXcd irrep2)
{
if (debugLevel > 3) {
	cout << "CharacterTable::" << __FUNCTION__ << endl;
}
	VectorXcd v1 = irrep1;
	VectorXcd v2 = irrep2;

	for ( int i = 0; i < _classes.size(); ++i )
	{
if (debugLevel > 3) {
		cout << (double)_classes[i].size() << " symmetry operations in class " << i+1 << endl;
}
		v1(i) = v1(i) * (double)_classes[i].size();
	}

	if ( debugLevel > 0 )
	{
		cout << "irreps " << v1.transpose() << " and " << v2.transpose() << " are ";

		if ( abs(v1.conjugate().dot(v2)) > symmetryPrecision )
			cout << "not ";
		cout << "orthogonal" << endl;
	}

	return ( abs(v1.conjugate().dot(v2)) < symmetryPrecision );
}

void CharacterTable::sortClasses( )
{
	int iter = 0;
	vector<vector<SymmOpPtr>> classes;

	if ( debugLevel > 0 ) cout << "sorting classes of " << _type << endl;

	do {
		++iter;
		classes = _classes;

		if ( debugLevel > 0 ) cout << "iteration " << iter << endl <<  to_string(classes != _classes) << endl;


		for ( size_t i = 0; i < _classes.size() - 1; ++i )
		{
			if ( lt(_classes[i+1], _classes[i] ) )
			{
				if ( debugLevel > 0)
					cout << "swapping " << _classes[i][0]->getSymbol() << " and " << _classes[i+1][0]->getSymbol() << endl;
				iter_swap(_classes.begin()+i+1, _classes.begin()+i);
			}
		}
		if ( debugLevel > 0 ) cout << "iteration " << iter << " done" << endl;

	}
	while ( classes != _classes );
	if ( debugLevel > 0 ) {cout << "Classes are sorted properly. Exiting now." << endl;
		cout << "classes are ";
		for ( auto it : _classes )
			cout << it.size() << it[0]->getSymbol() << " ";
		cout << endl;
	}
}

void CharacterTable::sortIrreps( )
{
	int iter = 0;
	vector<boost::shared_ptr<Irrep>> irreps;

	if ( debugLevel > 0 )
		cout << "sorting irreps of " << _type << endl;

	do {
		++iter;
		if ( debugLevel > 0 )
			cout << "iteration " << iter << endl;

		irreps = _irreps;

		for ( size_t i = 0; i < _irreps.size() - 1; ++i )
		{
			if ( _irreps[i+1]->lt(*_irreps[i] ) )
			{
				if ( debugLevel > 0)
					cout << "swapping " << _irreps[i]->getSymbol() << " and " << _irreps[i+1]->getSymbol() << endl;
				iter_swap(_irreps.begin()+i+1, _irreps.begin()+i);
			}
		}
	}
	while ( irreps != _irreps );
}

bool CharacterTable::allSymbolsDifferent()
{
	if ( debugLevel > 0 )
		cout << "checking if Mulliken symbols are unique" << endl;

	// check if A1 has been set
	if ( _irreps[0]->getNumber() == "1" && _irreps[1]->getNumber() == "" )
		return false;


	for ( size_t i = 0; i < _irreps.size(); ++i )
		for ( size_t j = i+1; j < _irreps.size(); ++j )
			if ( _irreps[i]->getSymbol() == _irreps[j]->getSymbol() )
			{
				if (debugLevel>0)
					cout << _irreps[i]->getSymbol() << " equals " << _irreps[j]->getSymbol() << endl;
				return false;
			}

	if ( debugLevel > 0 ) cout << "They are!" << endl;
	return true;
}

void CharacterTable::generateMullikenSymbols()
{
	if ( debugLevel > 0 )
		cout << "generating Mulliken symbols" << endl;

	// Degeneracy
	for ( size_t i = 0; i < _irreps.size(); ++i )
	{
		if ( _irreps[i]->getCharacters(0)(0).real() == 1 && _irreps[i]->getCharacters().size() == 1 )
		{
			if ( _mainRotationAxis->getOrderOfRotation() != 0 && _irreps[i]->getCharacters(0)(1).real() == -1 )
				_irreps[i]->setDegeneracy("B");
			else
				_irreps[i]->setDegeneracy("A");

		}

		if ( (_irreps[i]->getCharacters(0)(0).real() == 2 && _irreps[i]->getCharacters().size() == 1) || (_irreps[i]->getCharacters(0)(0).real() == 1 && _irreps[i]->getCharacters().size() == 2) )
			_irreps[i]->setDegeneracy("E");

		if ( (_irreps[i]->getCharacters(0)(0).real() == 3 && _irreps[i]->getCharacters().size() == 1) || (_irreps[i]->getCharacters(0)(0).real() == 1 && _irreps[i]->getCharacters().size() == 3) )
			_irreps[i]->setDegeneracy("T");

		if ( (_irreps[i]->getCharacters(0)(0).real() == 4 && _irreps[i]->getCharacters().size() == 1) || (_irreps[i]->getCharacters(0)(0).real() == 1 && _irreps[i]->getCharacters().size() == 4) )
			_irreps[i]->setDegeneracy("G");

		if ( (_irreps[i]->getCharacters(0)(0).real() == 5 && _irreps[i]->getCharacters().size() == 1) || (_irreps[i]->getCharacters(0)(0).real() == 1 && _irreps[i]->getCharacters().size() == 5) )
			_irreps[i]->setDegeneracy("H");
		if ( _irreps[i]->getDegeneracy() == "" ) _irreps[i]->setDegeneracy("X");

	}

	if ( debugLevel > 0 )
	{
		cout << "step one" << endl << "Mulliken symbols are" << endl;
		for ( auto it : _irreps )
			cout << it->getSymbol() << endl;
	}


	if (allSymbolsDifferent())
		return;

	// gerade / ungerade
	for ( size_t i = 0; i < _classes.size(); ++i )
	{
		if ( _classes[i][0]->getType() == "inversion" )
			posOfI = i;
	}

	if ( posOfI > -1 )
		for ( size_t i = 0; i < _irreps.size(); ++i )
		{
				if ( _irreps[i]->getCharacters(0)(posOfI).real() > 0 )
					_irreps[i]->setGerade("g");
				else
					_irreps[i]->setGerade("u");
		}
	else
		if (debugLevel > 0) cout << "no inversion found" << endl;



	if ( debugLevel > 0 ){cout << "step two" << endl << "Mulliken symbols are" << endl;for ( auto it : _irreps )cout << it->getSymbol() << endl;}

	if (allSymbolsDifferent())
		return;

	// primes
	for ( size_t i = 0; i < _classes.size(); ++i )
	{
		if ( _classes[i][0]->getSymbol() == "sigma_h" )
			posOfSigma = i;
	}

	// exception for cs
	if ( _type == "Cs" )
	{
		for ( size_t i = 0; i < _classes.size(); ++i )
		{
			if ( _classes[i][0]->getSymbol() == "sigma" )
				posOfSigma = i;
		}
	}

	if ( posOfSigma > -1 && posOfI == -1 )
		for ( size_t i = 0; i < _irreps.size(); ++i )
		{
				if ( _irreps[i]->getCharacters(0)(posOfSigma).real() > 0 )
					_irreps[i]->setPrimes("'");
				else
					_irreps[i]->setPrimes("''");
		}
	else
		if (debugLevel > 0) cout << "no sigma_h found" << endl;

	if (allSymbolsDifferent())
		return;

	// numbers
	if ( debugLevel > 0 ){cout << "step three" << endl << "Mulliken symbols are" << endl;for ( auto it : _irreps )cout << it->getSymbol() << endl;}

	pair<int,int> numbersA{1,1};
	pair<int,int> numbersB{1,1};
	pair<int,int> numbersE{1,1};
	pair<int,int> numbersT{1,1};
	pair<int,int> numbersG{1,1};
	pair<int,int> numbersH{1,1};

	// set the fully symmetric irrep to A1
	for ( size_t i = 0; i < _irreps.size(); ++i )
	{
		if ( _irreps[i]->getCharacters(0).real().minCoeff() > 0 )
		{
			if ( debugLevel > 0 ) cout << "setting irrep " << i + 1 << " to A1" << endl;
			_irreps[i]->setNumber(to_string(numbersA.first));
			++numbersA.first;
		}
	}

	int runs = 1;
	while ( !allSymbolsDifferent() )
	{
		if ( debugLevel > 0 )
			cout << "run " << runs << endl;
		if ( debugLevel > 0 ){cout << "Mulliken symbols are" << endl;for ( auto it : _irreps )cout << it->getSymbol() << endl;}

		vector<boost::shared_ptr<Irrep>> symbols = _irreps;
		vector<vector<boost::shared_ptr<Irrep>>> equals;

		start:
		if (debugLevel > 0)
		{
			cout << "equals:" << endl;
			for ( auto it : equals )
			{
				for ( auto it2 : it )
					cout << it2->getSymbol() << " ";
				cout << endl;
			}
		}


		while ( symbols.size() > 0 )
		{
			if ( debugLevel > 0 ) cout << "symbols left " << symbols.size() << endl;

			for ( size_t j = 0; j < equals.size(); ++j )
			{
				if ( symbols[0]->getSymbol() == equals[j][0]->getSymbol() )
				{
					equals[j].push_back(symbols[0]);
					symbols.erase(symbols.begin());
					goto start;
				}
			}
			equals.push_back(vector<boost::shared_ptr<Irrep>>{symbols[0]});
			symbols.erase(symbols.begin());
			goto start;
		}

		for ( size_t i = 0; i < equals.size(); ++i )
		{
			if ( debugLevel > 0 ) cout << "Symbol " << i + 1 << " occurs " << equals[i].size() << " times" << endl;

			if ( equals[i].size() == 1 && equals[i][0]->getDegeneracy() != "A")
			{
				if ( debugLevel > 0 ) cout << "therefore, I stop here." << endl;
				continue;
			}

			if ( debugLevel > 0 ) cout << "therefore, I continue." << endl;

			for ( size_t j = 0; j < equals[i].size(); ++j )
			{
				if ( equals[i][j]->getDegeneracy() == "A" && equals[i][j]->getNumber() == "" )
				{
					if ( ( equals[i][j]->getGerade() == "" && equals[i][j]->getPrimes() == "") || equals[i][j]->getGerade() == "g" || equals[i][j]->getPrimes() == "'" )
					{
						if ( debugLevel > 0 ) cout << "setting " << equals[i][j]->getSymbol() << " number " << i + 1 << " to A" << numbersA.first << endl;
						equals[i][j]->setNumber(to_string(numbersA.first));
						++numbersA.first;
					}
					else
					{
						if ( debugLevel > 0 ) cout << "setting " << equals[i][j]->getSymbol() << " number " << i + 1 << " to A" << numbersA.second << endl;
						equals[i][j]->setNumber(to_string(numbersA.second));
						++numbersA.second;
					}
				}
				if ( equals[i][j]->getDegeneracy() == "B" && equals[i][j]->getNumber() == "" )
				{
					if ( ( equals[i][j]->getGerade() == "" && equals[i][j]->getPrimes() == "") || equals[i][j]->getGerade() == "g" || equals[i][j]->getPrimes() == "'" )
					{
						equals[i][j]->setNumber(to_string(numbersB.first));
						if ( debugLevel > 0 ) cout << "setting " << equals[i][j]->getSymbol() << " number " << i + 1 << " to B" << numbersB.first << endl;
						++numbersB.first;
					}
					else
					{
						if ( debugLevel > 0 ) cout << "setting " << equals[i][j]->getSymbol() << " number " << i + 1 << " to B" << numbersB.second << endl;
						equals[i][j]->setNumber(to_string(numbersB.second));
						++numbersB.second;
					}
				}
				if ( equals[i][j]->getDegeneracy() == "E" && equals[i][j]->getNumber() == "" )
				{
					if ( ( equals[i][j]->getGerade() == "" && equals[i][j]->getPrimes() == "") || equals[i][j]->getGerade() == "g" || equals[i][j]->getPrimes() == "'" )
					{
						if ( debugLevel > 0 ) cout << "setting " << equals[i][j]->getSymbol() << " number " << i + 1 << " to E" << numbersE.first << endl;
						equals[i][j]->setNumber(to_string(numbersE.first));
						++numbersE.first;
					}
					else
					{
						if ( debugLevel > 0 ) cout << "setting " << equals[i][j]->getSymbol() << " number " << i + 1 << " to E" << numbersE.second << endl;
						equals[i][j]->setNumber(to_string(numbersE.second));
						++numbersE.second;
					}
				}
				if ( equals[i][j]->getDegeneracy() == "T" && equals[i][j]->getNumber() == "" )
				{
					if ( ( equals[i][j]->getGerade() == "" && equals[i][j]->getPrimes() == "") || equals[i][j]->getGerade() == "g" || equals[i][j]->getPrimes() == "'" )
					{
						if ( debugLevel > 0 ) cout << "setting " << equals[i][j]->getSymbol() << " number " << i + 1 << " to T" << numbersT.first << endl;
						equals[i][j]->setNumber(to_string(numbersT.first));
						++numbersT.first;
					}
					else
					{
						if ( debugLevel > 0 ) cout << "setting " << equals[i][j]->getSymbol() << " number " << i + 1 << " to T" << numbersT.second << endl;
						equals[i][j]->setNumber(to_string(numbersT.second));
						++numbersT.second;
					}
				}
				if ( equals[i][j]->getDegeneracy() == "G" && equals[i][j]->getNumber() == "" )
				{
					if ( ( equals[i][j]->getGerade() == "" && equals[i][j]->getPrimes() == "") || equals[i][j]->getGerade() == "g" || equals[i][j]->getPrimes() == "'" )
					{
						if ( debugLevel > 0 ) cout << "setting " << equals[i][j]->getSymbol() << " number " << i + 1 << " to G" << numbersG.first << endl;
						equals[i][j]->setNumber(to_string(numbersG.first));
						++numbersG.first;
					}
					else
					{
						if ( debugLevel > 0 ) cout << "setting " << equals[i][j]->getSymbol() << " number " << i + 1 << " to G" << numbersG.second << endl;
						equals[i][j]->setNumber(to_string(numbersG.second));
						++numbersG.second;
					}
				}
				if ( equals[i][j]->getDegeneracy() == "H" && equals[i][j]->getNumber() == "" )
				{
					if ( ( equals[i][j]->getGerade() == "" && equals[i][j]->getPrimes() == "") || equals[i][j]->getGerade() == "g" || equals[i][j]->getPrimes() == "'" )
					{
						if ( debugLevel > 0 ) cout << "setting " << equals[i][j]->getSymbol() << " number " << i + 1 << " to H" << numbersH.first << endl;
						equals[i][j]->setNumber(to_string(numbersH.first));
						++numbersH.first;
					}
					else
					{
						if ( debugLevel > 0 ) cout << "setting " << equals[i][j]->getSymbol() << " number " << i + 1 << " to H" << numbersH.second << endl;
						equals[i][j]->setNumber(to_string(numbersH.second));
						++numbersH.second;
					}
				}
			}
		}

		++runs;
		if ( runs > 1000 )
			break;
	}
}
