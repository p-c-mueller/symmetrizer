/*
 * molecule.h
 *
 *  Created on: Feb 6, 2023
 *      Author: pemueller
 */

#ifndef MOLECULE_H_
#define MOLECULE_H_


#include "atom.h"
#include <Eigen/Core>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <boost/serialization/vector.hpp>
#include "../characterTable/symmetryElements.h"
#include "../characterTable/pointGroup.h"
#include "../testing/pgTests.h"

using namespace Eigen;
using namespace std;

class BungeSTO;
typedef boost::shared_ptr<BungeSTO> BungeSTOPtr;

class UnitCell;
typedef boost::shared_ptr<UnitCell> UnitCellPtr;

class Molecule : public Atom {
public:
	Molecule();
	Molecule( int index, vector<AtomPtr> atoms);
	~Molecule();

	void setTransformMatrix( MatrixXd transformMatrix );
	void writeMolecularOrbitals( size_t spin );

	double getNumberOfElectrons();

	vector<AtomPtr> getAtoms();

	int getIndex();
	AtomPtr getAtom( int index );
	AtomPtr getAtom( size_t index );
	vector<BungeSTO*> getLocalBasisFunctions();
	size_t getNumberOfAtoms();
	size_t getNumberOfLocalBasisFunctions();
	size_t getIndexOfLocalBasisFunction(size_t mu, VectorXi lookupTable);
	MatrixXd getCoordinateSystem();
	string getName();
	void addAtomToName(int i, vector<pair<string,int> >& elements);

	double distanceFrom(boost::shared_ptr<Molecule> other);

	Vector3d getCenter();

	void findPointGroup();
	PGPtr getPointGroup();
	double getErrorOfSymmetryOperation( SymmetryOperation so );

	bool operator==(const Molecule &other) const;
	bool operator!=(const Molecule& other) const;

	// debugging stuff
	void printAllSymmetryOperations();


private:
	void _moleculeInit();

	void findManuallyDefinedPointGroup();
	void findSymmetryElements();
	void findSymmetryElementsLinear();
	void cleanUpPointGroup();

	void checkIdentity();
	void checkInversion();
	void findRotationAxes(size_t order = 0);
	void findRotationAxesLinear( size_t order = 0 );
	void searchMainRotationAxis(size_t order = 0);
	void searchC2Axes();
	void searchC2AxesLinear();
	void checkForMainRotationAxis();
	void checkForMainRotationAxis( Vector3d input );
	void checkForMainRotationReflectionAxis();
	void checkForMainRotationReflectionAxis(Vector3d input);
	void findReflection();
	void searchSigmaH();
	void searchSigmaV();
	void searchSigma();
	void searchSigmaLinear();
	void searchSigmaVAround( Vector3d axis );
	void findRotationReflection(size_t order = 0);
	void findRotationReflectionLinear();
	size_t getMaximumPossibleRotationAxis( size_t order = 0 );

	MatrixXd _transformMatrix;
	PGPtr _pointGroup;

	bool isLinear();
	bool isPlanar();
	int _index;
	vector<AtomPtr> _atoms;
	Vector3d _coordinatesReal;
};

typedef boost::shared_ptr<Molecule> MolPtr;




#endif /* MOLECULE_H_ */
