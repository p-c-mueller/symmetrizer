/*
 * symmetryElements.h
 *
 *  Created on: Feb 6, 2023
 *      Author: pemueller
 */

#ifndef SYMMETRYELEMENTS_H_
#define SYMMETRYELEMENTS_H_


#include <Eigen/Core>
#include <string>
#include <boost/shared_ptr.hpp>
#include "../testing/pgTests.h"

using namespace std;
using namespace Eigen;

class SymmetryOperation
{
public:
	SymmetryOperation();
	SymmetryOperation( string type, Vector3d center = Vector3d::Zero(), Vector3d rotationAxis = Vector3d::Zero(), int order = -1, int numberOfApplications = 1, int primes = 0 );
	SymmetryOperation( string type, SymmetryOperation other, int order = -1 , int numberOfApplications = 1);
	virtual ~SymmetryOperation();

	string getType();
	void setSymbol(string symbol);
	string getSymbol();
	int getOrderOfRotation();
	Vector3d getAxis();
	Vector3d getCenter();
	int getNumberOfApplications();
	Matrix3d getTransformMatrix();
	int getNumberOfPrimes();
	void setPrimes( int primes );

	Vector3d applySymmetryOperation( Vector3d input );

	bool operator==(const SymmetryOperation other) const;
	bool operator>(const SymmetryOperation other) const;
	bool operator<(const SymmetryOperation other) const;
	Matrix3d operator-( const SymmetryOperation other ) const;

	SymmetryOperation operator*(const SymmetryOperation other);

	// todo: error as mean displacement of atoms
	double getError();
	void setError(double error);

private:
	Vector3d applySymmetryOperation( Vector3d input, string type );
	Vector3d applyIdentity(Vector3d input);
	Vector3d applyInversion(Vector3d input);
	Vector3d applyRotation(Vector3d input);
	Vector3d applyReflection(Vector3d input);
	Vector3d applyRotationReflection(Vector3d input);

	string _symbol;
	string _type;
	int _orderOfRotation;
	int _primes;
	int _numberOfApplications;
	Vector3d _center;
	Vector3d _axis;
	double _error;
	Matrix3d _transformMatrix;
};

typedef boost::shared_ptr<SymmetryOperation> SymmOpPtr;


#endif /* SYMMETRYELEMENTS_H_ */
