/*
 * atom.cpp
 *
 *  Created on: Feb 6, 2023
 *      Author: pemueller
 */

#include "elements.h"
#include "atom.h"
#include "constants.h"
#include <iostream>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/case_conv.hpp>

Atom::Atom() {
	_initAtom();
}

void Atom::_initAtom() {
	_index = -1;
	_elementSymbol.clear();
	_coordinatesReal.setConstant(-1.0);
	_atomicNumber = 0;
}

Atom::Atom(int index, const string& elementSymbol, const Vector3d& coordinatesReal) {
	_initAtom();
	_index = index;
	_elementSymbol = elementSymbol;
	boost::trim(_elementSymbol);
	boost::to_lower(_elementSymbol);
	_atomicNumber = _convertElementSymbolToAtomicNumber(_elementSymbol);
	_coordinatesReal = coordinatesReal;
}

Atom::~Atom() {
}

const Vector3d& Atom::getCoordinatesReal() const {
	return _coordinatesReal;
}

int Atom::getIndex() const {
	return _index;
}

int Atom::getAtomicNumber() const {
	return _atomicNumber;
}

const string& Atom::getElementSymbol() const {
	return _elementSymbol;
}

string Atom::getElementSymbolForPrinting() const {
	string elementSymbol = _elementSymbol;
	elementSymbol[0] = toupper(elementSymbol[0]);
	return elementSymbol;
}

int Atom::_convertElementSymbolToAtomicNumber(const string &symbol) {
	return elementSymbols.right.at(symbol);
}

double Atom::distanceFrom(AtomPtr other) {

		return (_coordinatesReal - other->_coordinatesReal).norm();
}

bool Atom::operator ==(const Atom& other) const {
	return ((_elementSymbol == other._elementSymbol) && ((_coordinatesReal - other._coordinatesReal ).norm() < symmetryPrecision ));
}

bool Atom::operator !=(const Atom& other) const {
	return !(*this == other);
}


