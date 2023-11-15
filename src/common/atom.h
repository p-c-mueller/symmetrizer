/*
 * atom.h
 *
 *  Created on: Feb 6, 2023
 *      Author: pemueller
 */

#ifndef ATOM_H_
#define ATOM_H_

#include <Eigen/Core>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <boost/serialization/vector.hpp>

using namespace Eigen;
using namespace std;

class Atom {
public:
//	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Atom();
	Atom(int index, const string& elementSymbol, const Vector3d& coordinatesFrac);
	virtual ~Atom();

	const Vector3d& getCoordinatesReal() const;

	int getIndex() const;

	int getAtomicNumber() const;

	string getElementSymbolForPrinting() const;

	const string& getElementSymbol() const;

	double distanceFrom(boost::shared_ptr<Atom> other);

	bool operator==(const Atom& other) const;
	bool operator!=(const Atom& other) const;

private:
	int _index;
	string _elementSymbol;
	Vector3d _coordinatesReal;
	int _atomicNumber;
	int _convertElementSymbolToAtomicNumber(const string &symbol);
	void _initAtom();
};

typedef boost::shared_ptr<Atom> AtomPtr;

#endif /* ATOM_H_ */
