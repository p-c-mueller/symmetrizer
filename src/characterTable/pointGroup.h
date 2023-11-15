/*
 * pointGroup.h
 *
 *  Created on: Feb 6, 2023
 *      Author: pemueller
 */

#ifndef POINTGROUP_H_
#define POINTGROUP_H_
#include "symmetryElements.h"
#include "characterTable.h"
#include "../common/atom.h"
#include <vector>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include "../testing/pgTests.h"

typedef boost::multi_array<string,2> MatrixXs;

class PointGroup
{
public:
	PointGroup( vector<AtomPtr> atoms );
	~PointGroup();
	void clear();

	void findNameOfPointGroup();
	void addSymmetryOperation( SymmOpPtr so );
	bool hasMultipleC2();
	vector<SymmOpPtr> getSymmetryOperations();
	void determineMainRotationAxis();
	void setSymmetryOperationClasses();
	string getName();
	string getType();
	SymmOpPtr getMainRotationAxis();
	vector<SymmOpPtr> getMainRotationAxes();
	int getOrderOfMainRotationAxis();
	MatrixXs getSymmetryElementMatrix();
	void generateCharacterTable();
	CTPtr getCharacterTable();
	bool isSymmetryOperationPresent( SymmOpPtr so );
	void setLinear();
	void setPrimes();
	VectorXd getDistanceFromSO( SymmOpPtr so );
	double getDistanceFromSO(AtomPtr atom, SymmOpPtr so);
	VectorXd weightWithAtom(VectorXd distances);
	double weightWithAtom(double distances, AtomPtr atoms);
	void generateIrreps();

private:
	// Defined in generateIrreps.cpp
	void generateIrrepsC1();
	void generateIrrepsCs();
	void generateIrrepsCi();
	void generateIrrepsCn();
	void generateIrrepsCnh();
	void generateIrrepsCnv();
	void generateIrrepsCuv();
	void generateIrrepsSn();
	void generateIrrepsDn();
	void generateIrrepsDnh();
	void generateIrrepsDnd();
	void generateIrrepsDuh();
	void generateIrrepsTd();
	void generateIrrepsOh();
	void generateIrrepsIh();

	// Defined in fixesForPointGroup.cpp
	void setPrimesCnv();
	void setPrimesDn();
	void setPrimesDndh();
	void setPrimes(vector<SymmOpPtr> soVector, int primes);
	void fixD2();
	void fixD2d();
	void fixD2h();
	void setS2nAsMainRotationAxis();

	size_t getPositionOfSymmetryOperation(SymmOpPtr so);

	bool isThereARotationOfHigherOrder(SymmOpPtr input);
	bool isAxisInMainRotationAxis( SymmOpPtr axis );

	// decision tree for determination of point group
	bool isLinear();
	bool hasInversion();
	bool hasTwoOrMoreCn();
	bool hasC5();
	bool hasCn();
	bool hasSigma();
	bool isDihedral();
	bool hasSigmaH();
	bool hasSigmaHDihedral();
	bool hasNSigmaD();
	bool hasNSigmaV();
	bool hasS2N();

	bool isSigmaV( SymmOpPtr so );
	bool isSigmaH( SymmOpPtr so );
	bool isSigmaD( SymmOpPtr so );

	bool _isLinear = false;
	vector<SymmOpPtr> _symmetryOperations;
	vector<SymmOpPtr> _mainRotationAxis;
	string _name;
	string _type;
	CTPtr _characterTable;
	vector<AtomPtr> _atoms;
	bool isGreater(vector<SymmOpPtr> v1, vector<SymmOpPtr> v2);
};

typedef boost::shared_ptr<PointGroup> PGPtr;

#endif /* POINTGROUP_H_ */
