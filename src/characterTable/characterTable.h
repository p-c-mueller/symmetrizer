/*
 * characterTable.h
 *
 *  Created on: Feb 23, 2023
 *      Author: pemueller
 */

#ifndef CHARACTERTABLE_H_
#define CHARACTERTABLE_H_

#include <Eigen/Core>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "symmetryElements.h"
#include "../testing/pgTests.h"
#include "irrep.h"

using namespace Eigen;
using namespace std;

class CharacterTable
{
public:
	CharacterTable(string type, SymmOpPtr mra);
	~CharacterTable();
	void clear();

	void sortClasses();
	void sortIrreps( );
	void generateIrreps();
	void generateMullikenSymbols();
	void combineC2AndSigmaV();

	void addClass(vector<SymmOpPtr> input);
	vector<vector<SymmOpPtr>> getClasses();
	void addIrrep( VectorXcd irrep );
	vector<boost::shared_ptr<Irrep>> getIrreps();

private:
	bool allSymbolsDifferent();

	bool irrepAlreadyExists( VectorXcd irrep );
	bool irrepsOrthogonal(VectorXcd irrep1, VectorXcd irrep2);

	bool lt(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2);
	bool ltDuh( vector<SymmOpPtr> class1, vector<SymmOpPtr> class2 );
	bool ltCuv(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2);
	bool ltIh(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2);
	bool ltOh(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2);
	bool ltTd(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2);
	bool ltDnh(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2);
	bool ltDnd(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2);
	bool ltDn(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2);
	bool ltCnh(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2);
	bool ltCnv(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2);
	bool ltCn(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2);
	bool ltSn(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2);
	bool ltCs(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2);
	bool ltCi(vector<SymmOpPtr> class1, vector<SymmOpPtr> class2);

	bool isIn( SymmOpPtr so, vector<SymmOpPtr> soV )
	{
		for ( auto it : soV )
		{
			if ( (it->getTransformMatrix() - so->getTransformMatrix()).cwiseAbs().maxCoeff() < 0.001 )
				return true;
		}
		return false;
	}

	int posOfSigma = -1;
	int posOfI = -1;

	vector<vector<SymmOpPtr>> _classes;
	vector<boost::shared_ptr<Irrep>> _irreps;
	int _order;
	string _type;
	SymmOpPtr _mainRotationAxis;
};

typedef boost::shared_ptr<CharacterTable> CTPtr;

#endif /* CHARACTERTABLE_H_ */
