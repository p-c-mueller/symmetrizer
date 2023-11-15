/*
 * generateIrreps.cpp
 *
 *  Created on: Jun 19, 2023
 *      Author: pemueller
 */

#include "pointGroup.h"
#include "../common/constants.h"
#include <iostream>
#include "../testing/pgTests.h"
#include "../IO/writer.h"

void PointGroup::generateIrreps()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
	cout << "generating irreps for " << _type << endl;
}
	if ( _type == "C1" )
		generateIrrepsC1();

	if ( _type == "Cs" )
		generateIrrepsCs();

	if ( _type == "Ci" )
		generateIrrepsCi();

	if ( _type == "Cn" )
		generateIrrepsCn();

	if ( _type == "Cnh" )
		generateIrrepsCnh();

	if ( _type == "Cnv" )
		generateIrrepsCnv();

	if ( _type == "Cuv" )
		generateIrrepsCuv();

	if ( _type == "Sn" )
		generateIrrepsSn();

	if ( _type == "Dn" )
		generateIrrepsDn();

	if ( _type == "Dnh" )
		generateIrrepsDnh();

	if ( _type == "Dnd" )
		generateIrrepsDnd();

	if ( _type == "Duh" )
		generateIrrepsDuh();

	if ( _type == "Td" )
		generateIrrepsTd();

	if ( _type == "Oh" )
		generateIrrepsOh();

	if ( _type == "Ih" )
		generateIrrepsIh();
}

void PointGroup::generateIrrepsC1()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	VectorXcd irrep = VectorXcd::Zero(1);
	irrep << 1;
	_characterTable->addIrrep(irrep);
}

void PointGroup::generateIrrepsCs()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	VectorXcd irrep = VectorXcd::Zero(2);
	irrep << 1,1;
	_characterTable->addIrrep(irrep);
	irrep.setZero();
	irrep << 1,-1;
	_characterTable->addIrrep(irrep);
}

void PointGroup::generateIrrepsCi()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	VectorXcd irrep = VectorXcd::Zero(2);
	irrep << 1,1;
	_characterTable->addIrrep(irrep);
	irrep.setZero();
	irrep << 1,-1;
	_characterTable->addIrrep(irrep);
}

void PointGroup::generateIrrepsCn()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	const int order = _characterTable->getClasses().size();
	VectorXcd irrep = VectorXcd::Zero(order);
	auto mra = _mainRotationAxis[0];
	const int n = mra->getOrderOfRotation();
	for ( int mu = -n; mu <= n; ++mu )
	{
		for ( int k = 0; k < n; ++k )
		{
			SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
			irrep(getPositionOfSymmetryOperation(so)) = exp( - TWOPITimesI * (double)k * (double)mu / (double)n );
		}

		_characterTable->addIrrep(irrep);
		irrep.setZero();
	}
}

void PointGroup::generateIrrepsCnh()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	const int order = _characterTable->getClasses().size();
	VectorXcd irrep = VectorXcd::Zero(order);
	auto mra = _mainRotationAxis[0];
	const int n = mra->getOrderOfRotation();

	if (n % 2 == 1 )
	{
		for ( int mu = -n; mu <= n; ++mu )
		{
			for ( int k = 0; k < 2*n; ++k )
			{
				SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotationReflection", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
				complex<double> omega = exp(- PiTimesI / (double)n);
				irrep(getPositionOfSymmetryOperation(so)) = pow(pow( omega, (n + 2) * mu ), k);
			}

			_characterTable->addIrrep(irrep);
			irrep.setZero();
		}
	}
	else
	{
		for ( int i = 0; i < 2; ++i )
		{
			for ( int mu = -n; mu <= n; ++mu )
			{
				for ( int k = 0; k < n; ++k )
				{
					SymmetryOperation rotation("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
					SymmetryOperation identity("identity", mra->getCenter());
					SymmOpPtr so = boost::make_shared<SymmetryOperation>(rotation * identity);
					irrep(getPositionOfSymmetryOperation(so)) = exp( - TWOPITimesI * (double)k * (double)mu / (double)n );
				}

				for ( int k = 0; k < n; ++k )
				{
					SymmetryOperation rotation("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
					SymmetryOperation reflection("reflection", mra->getCenter(), mra->getAxis());
					SymmOpPtr so = boost::make_shared<SymmetryOperation>(rotation * reflection);
					irrep(getPositionOfSymmetryOperation(so)) = exp( - TWOPITimesI * (double)k * (double)mu / (double)n ) * pow(-1, i);
				}
				_characterTable->addIrrep(irrep);
				irrep.setZero();
			}
		}
	}
}

void PointGroup::generateIrrepsCnv()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	const int order = _characterTable->getClasses().size();
	VectorXcd irrep = VectorXcd::Zero(order);
	auto mra = _mainRotationAxis[0];
	const int n = mra->getOrderOfRotation();
	if ( n % 2 == 0 )
	{
		// non-degenerate irreps
		for ( int rho = 1; rho > -2; rho -= 2 )
		{
			for ( int nu = 1; nu > -2; nu -= 2 )
			{
				for ( int k = 0; k < n; ++k )
				{
					SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
					irrep(getPositionOfSymmetryOperation(so)) = pow(rho, k);
				}
				for ( auto so : _symmetryOperations )
				{
					if ( so->getType() == "reflection" )
					{
if (debugLevel > 3) {
						cout <<  rho << "^" << so->getNumberOfPrimes() + 1 << " = " << pow(rho, so->getNumberOfPrimes() + 1) << endl;
}
						irrep(getPositionOfSymmetryOperation(so)) = nu * pow(rho, so->getNumberOfPrimes() + 1);
					}
				}
				_characterTable->addIrrep(irrep);
				irrep.setZero();
			}
		}

		// degenerate irreps
		for ( int mu = 1; mu < n/2; ++mu )
		{
			for ( int k = 0; k < n; ++k )
			{
				SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
				irrep(getPositionOfSymmetryOperation(so)) = 2 * cos( TWOPI * (double) k * (double) mu / (double) n );
			}
			for ( auto so : _symmetryOperations )
			{
				if ( so->getType() == "reflection" )
				{
					irrep(getPositionOfSymmetryOperation(so)) = 0;
				}
			}

			if (debugLevel > 0)
				cout << irrep.transpose() << endl;

			_characterTable->addIrrep(irrep);
			irrep.setZero();
		}
	}
	else
	{
		for ( int nu = 1; nu > -2; nu-=2 )
		{
			for ( int k = 0; k < n; ++k )
			{
				SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
				irrep(getPositionOfSymmetryOperation(so)) = 1;
			}
			for ( auto so : _symmetryOperations )
			{
				if ( so->getType() == "reflection" )
				{
					irrep(getPositionOfSymmetryOperation(so)) = nu;
				}
			}
			_characterTable->addIrrep(irrep);
			irrep.setZero();
		}
		// degenerate irreps
		for ( int mu = 1; mu < n; ++mu )
		{
			for ( int k = 0; k < n; ++k )
			{
				SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
				irrep(getPositionOfSymmetryOperation(so)) = 2 * cos( TWOPI * (double) k * (double) mu / (double) n );
			}
			for ( auto so : _symmetryOperations )
			{
				if ( so->getType() == "reflection" )
				{
					irrep(getPositionOfSymmetryOperation(so)) = 0;
				}
			}
			_characterTable->addIrrep(irrep);
			irrep.setZero();
		}
	}
}

void PointGroup::generateIrrepsCuv()
{
	generateIrrepsCnv();
}

void PointGroup::generateIrrepsSn()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	const int order = _characterTable->getClasses().size();
	VectorXcd irrep = VectorXcd::Zero(order);
	auto mra = _mainRotationAxis[0];
	const int n = mra->getOrderOfRotation();
	if ( (n/2) % 2 == 0 )
	{
		for ( int mu = -n; mu <= n; ++mu )
		{
			for ( int k = 0; k < 2*n; ++k )
			{
				SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotationReflection", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
				complex<double> omega = exp(- PiTimesI / (double)(n/2));
				irrep(getPositionOfSymmetryOperation(so)) = pow(pow( omega, ((n/2) + 1) * mu ), k);
			}

			_characterTable->addIrrep(irrep);
			irrep.setZero();
		}
	}
	else
	{
		for ( int i = 0; i < 2; ++i )
		{
			for ( int mu = -n; mu <= n; ++mu )
			{
				for ( int k = 0; k < n; ++k )
				{
					SymmetryOperation rotation("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation()/2, k);
					SymmetryOperation identity("identity", mra->getCenter());
					SymmOpPtr so = boost::make_shared<SymmetryOperation>(rotation * identity);
					irrep(getPositionOfSymmetryOperation(so)) = exp( - TWOPITimesI * (double)k * (double)mu / (double)(n/2) );
				}

				for ( int k = 0; k < n; ++k )
				{
					SymmetryOperation rotation("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation()/2, k);
					SymmetryOperation reflection("inversion", mra->getCenter(), mra->getAxis());
					SymmOpPtr so = boost::make_shared<SymmetryOperation>(rotation * reflection);
					irrep(getPositionOfSymmetryOperation(so)) = exp( - TWOPITimesI * (double)k * (double)mu / (double)(n/2) ) * pow(-1, i);
				}
				_characterTable->addIrrep(irrep);
				irrep.setZero();
			}
		}
	}
}

void PointGroup::generateIrrepsDn()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	const int order = _characterTable->getClasses().size();
	VectorXcd irrep = VectorXcd::Zero(order);
	auto mra = _mainRotationAxis[0];
	const int n = mra->getOrderOfRotation();
	if ( n % 2 == 0 )
	{
		// non-degenerate irreps
		for ( int rho = 1; rho > -2; rho -= 2 )
		{
			for ( int nu = 1; nu > -2; nu -= 2 )
			{
				for ( int k = 0; k < n; ++k )
				{
					SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
					irrep(getPositionOfSymmetryOperation(so)) = pow(rho, k);
				}
				for ( auto so : _symmetryOperations )
				{
					if ( so->getOrderOfRotation() == 2 && so->getAxis().dot(mra->getAxis()) < symmetryPrecision )
					{
						irrep(getPositionOfSymmetryOperation(so)) = nu * pow(rho, so->getNumberOfPrimes() + 1);
					}
				}
				_characterTable->addIrrep(irrep);
				irrep.setZero();
			}
		}

		// degenerate irreps
		for ( int mu = 1; mu < n - 1; ++mu )
		{
			for ( int k = 0; k < n; ++k )
			{
				SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
				irrep(getPositionOfSymmetryOperation(so)) = 2 * cos( TWOPI * (double) k * (double) mu / (double) n );
			}
			for ( auto so : _symmetryOperations )
			{
				if ( so->getOrderOfRotation() == 2 && so->getAxis().dot(mra->getAxis()) < symmetryPrecision )
				{
					irrep(getPositionOfSymmetryOperation(so)) = 0;
				}
			}

			_characterTable->addIrrep(irrep);
			irrep.setZero();
		}

	}
	else
	{
		// non-degenerate irreps
		for ( int nu = 1; nu > -2; nu-=2 )
		{
			for ( int k = 0; k < n; ++k )
			{
				SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
				irrep(getPositionOfSymmetryOperation(so)) = 1;
			}
			for ( auto so : _symmetryOperations )
			{
				if ( so->getOrderOfRotation() == 2 && so->getAxis().dot(mra->getAxis()) < symmetryPrecision )
				{
					irrep(getPositionOfSymmetryOperation(so)) = nu;
				}
			}

			_characterTable->addIrrep(irrep);
			irrep.setZero();
		}
		// degenerate irreps
		for ( int mu = 1; mu < n; ++mu )
		{
			for ( int k = 0; k < n; ++k )
			{
				SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
				irrep(getPositionOfSymmetryOperation(so)) = 2 * cos( TWOPI * (double) k * (double) mu / (double) n );
			}
			for ( auto so : _symmetryOperations )
			{
				if ( so->getOrderOfRotation() == 2 && so->getAxis().dot(mra->getAxis()) < symmetryPrecision )
				{
					irrep(getPositionOfSymmetryOperation(so)) = 0;
				}
			}

			_characterTable->addIrrep(irrep);
			irrep.setZero();
		}
	}
}

void PointGroup::generateIrrepsDnh()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	const int order = _characterTable->getClasses().size();
	VectorXcd irrep = VectorXcd::Zero(order);
	auto mra = _mainRotationAxis[0];
	const int n = mra->getOrderOfRotation();
	SymmOpPtr inversion = boost::make_shared<SymmetryOperation >("inversion", mra->getCenter());

	if ( n % 2 == 0 )
	{
		for ( int i = 1; i > -2; i-=2 )
		{
			// Dn * E
			// non-degenerate irreps
			for ( int rho = 1; rho > -2; rho -= 2 )
			{
				for ( int nu = 1; nu > -2; nu -= 2 )
				{
					for ( int k = 0; k < n; ++k )
					{
						SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
						irrep(getPositionOfSymmetryOperation(so)) = pow(rho, k);
					}
					for ( auto so : _symmetryOperations )
					{
						if ( so->getOrderOfRotation() == 2 && so->getAxis().dot(mra->getAxis()) < symmetryPrecision )
						{
							irrep(getPositionOfSymmetryOperation(so)) = nu * pow(rho, so->getNumberOfPrimes() + 1);
						}
					}

					// Dn * i
					for ( int k = 0; k < n; ++k )
					{
						SymmetryOperation rotation ("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
						SymmOpPtr so = boost::make_shared<SymmetryOperation>(*inversion * rotation);
						irrep(getPositionOfSymmetryOperation(so)) = pow(rho, k) * i;
					}
					for ( auto so : _symmetryOperations )
					{
						if ( so->getOrderOfRotation() == 2 && so->getAxis().dot(mra->getAxis()) < symmetryPrecision )
						{
							SymmOpPtr so2 = boost::make_shared<SymmetryOperation>(*so * *inversion);
							irrep(getPositionOfSymmetryOperation(so2)) = nu * pow(rho, so->getNumberOfPrimes() + 1) * i;
						}
					}
					_characterTable->addIrrep(irrep);
					irrep.setZero();
				}
			}

			// degenerate irreps
			for ( int mu = 1; mu < n/2; ++mu )
			{
				for ( int i = 1; i > -2; i-=2 )
				{
					// Dn * E
					for ( int k = 0; k < n; ++k )
					{
						SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
						irrep(getPositionOfSymmetryOperation(so)) = 2 * cos( TWOPI * (double) k * (double) mu / (double) n );
					}

					// Dn * i
					for ( int k = 0; k < n; ++k )
					{
						SymmetryOperation rotation ("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
						SymmOpPtr so = boost::make_shared<SymmetryOperation>(rotation * *inversion);
						irrep(getPositionOfSymmetryOperation(so)) = 2 * cos( TWOPI * (double) k * (double) mu / (double) n ) * i;
					}

					_characterTable->addIrrep(irrep);
					irrep.setZero();
				}
			}
		}
	}
	else
	{
		// non-degenerate irreps
		for ( int nu = 1; nu > -2; nu -= 2 )
		{
			// Sn
			for ( int k = 0; k <= n; ++k )
			{
				SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotationReflection", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
				irrep(getPositionOfSymmetryOperation(so)) = 1;
			}

			// Cn
			for ( int k = 0; k < n; ++k )
			{
				SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
				irrep(getPositionOfSymmetryOperation(so)) = 1;
			}

			// C2
			for ( auto c2 : _symmetryOperations )
			{
				if ( c2->getOrderOfRotation() == 2 && c2->getAxis().dot(mra->getAxis()) < symmetryPrecision )
				{
					irrep(getPositionOfSymmetryOperation(c2)) = nu;
				}
			}

			// sigma_v
			for ( auto sigma : _symmetryOperations )
			{
				if ( sigma->getType() == "reflection" && sigma->getAxis().dot(mra->getAxis()) < symmetryPrecision )
				{
					irrep(getPositionOfSymmetryOperation(sigma)) = nu;
				}
			}
			_characterTable->addIrrep(irrep);
			irrep.setZero();
		}

		// degenerate irreps
		for ( int mu = 1; mu < n; ++mu )
		{
			// Sn
			for ( int k = 0; k <= n; ++k )
			{
				SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotationReflection", mra->getCenter(), mra->getAxis(), n, k);
				irrep(getPositionOfSymmetryOperation(so)) = 2 * cos( M_PI * (2*(double)k + (double)n)* (double)mu / (double)n );
			}

			// Cn
			for ( int k = 0; k < n; ++k )
			{
				SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotation", mra->getCenter(), mra->getAxis(), n, k);
				irrep(getPositionOfSymmetryOperation(so)) = 2 * cos( TWOPI * (double)k * (double)mu / (double)n );
			}

			_characterTable->addIrrep(irrep);
			irrep.setZero();
		}
	}
}

void PointGroup::generateIrrepsDuh()
{
	if( debugLevel > 0 )
		cout << "generating irreps for " << _type << endl;

	const int order = _characterTable->getClasses().size();
	VectorXcd irrep = VectorXcd::Zero(order);
	auto mra = _mainRotationAxis[0];
	const int n = mra->getOrderOfRotation();
	SymmOpPtr inversion = boost::make_shared<SymmetryOperation >("inversion", mra->getCenter());
	SymmOpPtr sigma;

	for ( auto so : _symmetryOperations )
	{
		if ( so->getType() == "reflection" && so->getAxis().dot(mra->getAxis()) < 0.001) // don't want to have the sigma_h here
		{
			sigma = so;
			break;
		}
	}

	for ( int mu = 0; mu < 2; ++mu )
	{
		for ( int rho = 0; rho < 2; ++rho )
		{
			if ( debugLevel > 0 ) cout << getIrrepVectorAsString(irrep) << endl;
			// non-degenerate irreps = A1g, A2g, A1u, A2u
			// E
			irrep(getPositionOfSymmetryOperation(boost::make_shared<SymmetryOperation>("identity"))) = 1;

			if ( debugLevel > 0 ) cout << getIrrepVectorAsString(irrep) << endl;

			// Cn

			for ( size_t i = n; i > 1; --i )
			{
				if ( i % n != 0 )
					continue;

				SymmOpPtr Cn = boost::make_shared<SymmetryOperation>( "rotation", Vector3d::Zero(), mra->getAxis(), i );
				irrep(getPositionOfSymmetryOperation(Cn)) = 1;
			}

			if ( debugLevel > 0 ) cout << getIrrepVectorAsString(irrep) << endl;

			// sigma_v
			irrep(getPositionOfSymmetryOperation(sigma)) = pow(-1, mu);
			if ( debugLevel > 0 ) cout << getIrrepVectorAsString(irrep) << endl;

			// i
			irrep(getPositionOfSymmetryOperation(boost::make_shared<SymmetryOperation>("inversion"))) = pow(-1, rho);
			if ( debugLevel > 0 ) cout << getIrrepVectorAsString(irrep) << endl;

			// Sn
			for ( size_t i = n; i > 1; --i )
			{
				if ( n % i != 0 )
					continue;

				SymmOpPtr Cn = boost::make_shared<SymmetryOperation>( "rotationReflection", Vector3d::Zero(), mra->getAxis(), i );
				irrep(getPositionOfSymmetryOperation(Cn)) = pow(-1, rho);
			}
			if ( debugLevel > 0 ) cout << getIrrepVectorAsString(irrep) << endl;

			// C2
			irrep(getPositionOfSymmetryOperation(boost::make_shared<SymmetryOperation>("rotation",*sigma,2))) = pow(-1, mu) * pow(-1, rho);
			if ( debugLevel > 0 ) cout << getIrrepVectorAsString(irrep) << endl;

			_characterTable->addIrrep(irrep);
			irrep.setZero();
		}
	}

	// degenerate irreps = Eng, Enu
	for ( int step = 0; step < infinityRotationSteps; ++step )
	{
		for ( int mu = 0; mu < 2; ++mu )
		{
			for ( int rho = 0; rho < 2; ++rho )
			{
				// E
				irrep(getPositionOfSymmetryOperation(boost::make_shared<SymmetryOperation>("identity"))) = 2;

				// Cn

				// sigma_v

				// i

				// Sn

				// C2
			}
		}

		_characterTable->addIrrep(irrep);
		irrep.setZero();
	}
}

void PointGroup::generateIrrepsDnd()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	const int order = _characterTable->getClasses().size();
	VectorXcd irrep = VectorXcd::Zero(order);
	auto mra = _mainRotationAxis[0];
	const int n = mra->getOrderOfRotation();

	if( debugLevel > 0 )
		cout << "generating irreps for " << _type << " n = " << n << endl;


	SymmOpPtr inversion = boost::make_shared<SymmetryOperation>("inversion", mra->getCenter());
	if ( n % 2 == 0 )
	{
		const int nOverTwo = n/2;
		// non-degenerate irreps
		for ( int rho = 1; rho > -2; rho-=2 )
		{
			for ( int nu = 1; nu > -2; nu-=2 )
			{
				// S_n
				for ( int k = 0; k < 2*n; ++k )
				{
					SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotationReflection", mra->getCenter(), mra->getAxis(), n, k);
					irrep(getPositionOfSymmetryOperation(so)) = pow(rho, k);
				}

				// C_2
				for ( auto c2 : _symmetryOperations )
				{
					if ( c2->getOrderOfRotation() == 2 && c2->getAxis().dot(mra->getAxis()) < symmetryPrecision )
					{
						irrep(getPositionOfSymmetryOperation(c2)) = nu;
					}
				}

				// sigma
				for ( auto sigma : _symmetryOperations )
				{
					if ( sigma->getType() == "reflection" && sigma->getAxis().dot(mra->getAxis()) < symmetryPrecision )
					{
						irrep(getPositionOfSymmetryOperation(sigma)) = nu * rho;
					}
				}
				_characterTable->addIrrep(irrep);
				irrep.setZero();
			}
		}

		// degenerate irreps
		for ( int mu = 1; mu < n/2; ++mu )
		{
			// S_n
			for ( int k = 0; k < n; ++k )
			{
				SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotationReflection", mra->getCenter(), mra->getAxis(), n, k);
				irrep(getPositionOfSymmetryOperation(so)) = 2 * cos( (double)(nOverTwo+1) * (double)k * (double)mu * M_PI / (double)nOverTwo );
			}

			// C_2
			for ( auto c2 : _symmetryOperations )
			{
				if ( c2->getOrderOfRotation() == 2 && c2->getAxis().dot(mra->getAxis()) < symmetryPrecision )
				{
					irrep(getPositionOfSymmetryOperation(c2)) = 0;
				}
			}

			// sigma
			for ( auto sigma : _symmetryOperations )
			{
				if ( sigma->getType() == "reflection" && sigma->getAxis().dot(mra->getAxis()) < symmetryPrecision )
				{
					irrep(getPositionOfSymmetryOperation(sigma)) = 0;
				}
			}
			_characterTable->addIrrep(irrep);
			irrep.setZero();
		}
	}
	else
	{

		for ( int i = 1; i > -2; i-=2 )
		{
			// non-degenerate irreps
			for ( int nu = 1; nu > -2; nu -= 2 )
			{
				// Dn * E
				for ( int k = 0; k < n; ++k )
				{
					SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
					irrep(getPositionOfSymmetryOperation(so)) = 1;
				}
				for ( auto so : _symmetryOperations )
				{
					if ( so->getOrderOfRotation() == 2 && so->getAxis().dot(mra->getAxis()) < symmetryPrecision )
					{
						irrep(getPositionOfSymmetryOperation(so)) = nu;
					}
				}

				// Dn * i
				for ( int k = 0; k < n; ++k )
				{
					SymmetryOperation rotation ("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
					SymmOpPtr so = boost::make_shared<SymmetryOperation>(*inversion * rotation);
					//if ( abs(irrep(getPositionOfSymmetryOperation(so))) == complexZero )
						irrep(getPositionOfSymmetryOperation(so)) = i;
				}
				for ( auto so : _symmetryOperations )
				{
					if ( so->getOrderOfRotation() == 2 && so->getAxis().dot(mra->getAxis()) < symmetryPrecision )
					{
						SymmOpPtr so2 = boost::make_shared<SymmetryOperation>(*so * *inversion);
						//if ( irrep(getPositionOfSymmetryOperation(so)) == complexZero )
							irrep(getPositionOfSymmetryOperation(so2)) = nu * i;
					}
				}

				_characterTable->addIrrep(irrep);
				irrep.setZero();
			}

			// degenerate irreps
			for ( int mu = 1; mu < n; ++mu )
			{
				for ( int i = 1; i > -2; i-=2 )
				{
					// Dn * E
					for ( int k = 0; k < n; ++k )
					{
						SymmOpPtr so = boost::make_shared<SymmetryOperation>("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
						irrep(getPositionOfSymmetryOperation(so)) = 2 * cos( TWOPI * (double) k * (double) mu / (double) n );
					}
					for ( auto so : _symmetryOperations )
					{
						if ( so->getOrderOfRotation() == 2 && so->getAxis().dot(mra->getAxis()) < symmetryPrecision )
						{
							irrep(getPositionOfSymmetryOperation(so)) = 0;
						}
					}

					// Dn * i
					for ( int k = 0; k < n; ++k )
					{
						SymmetryOperation rotation ("rotation", mra->getCenter(), mra->getAxis(), mra->getOrderOfRotation(), k);
						SymmOpPtr so = boost::make_shared<SymmetryOperation>(rotation * *inversion);
						irrep(getPositionOfSymmetryOperation(so)) = 2 * cos( TWOPI * (double) k * (double) mu / (double) n ) * i;
					}
					for ( auto so : _symmetryOperations )
					{
						if ( so->getOrderOfRotation() == 2 && so->getAxis().dot(mra->getAxis()) < symmetryPrecision )
						{
							SymmOpPtr so2 = boost::make_shared<SymmetryOperation>(*so * *inversion);
							irrep(getPositionOfSymmetryOperation(so2)) = 0;
						}
					}

					_characterTable->addIrrep(irrep);
					irrep.setZero();
				}
			}
		}
	}
}

void PointGroup::generateIrrepsTd()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	// non-degenerate irreps
	VectorXcd v = VectorXcd::Zero(5);
	v << 1,1,1,1,1;
	_characterTable->addIrrep(v);
	v.setZero();
	v << 1,1,1,-1,-1;
	_characterTable->addIrrep(v);
	// double degenerate irreps
	v.setZero();
	v << 2,-1,2,0,0;
	_characterTable->addIrrep(v);

	// triple degenerate irreps
	v.setZero();
	v << 3,0,-1,1,-1;
	_characterTable->addIrrep(v);
	v.setZero();
	v << 3,0,-1,-1,1;
	_characterTable->addIrrep(v);
}

void PointGroup::generateIrrepsOh()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	VectorXcd v = VectorXcd::Zero(10);

	// non-degenerate irreps
	v << 1,1,1,1,1,1,1,1,1,1;
	_characterTable->addIrrep(v);
	v.setZero();
	v << 1 ,1,-1,-1,1,1,-1,1,1,-1;
	_characterTable->addIrrep(v);
	v.setZero();
	v << 1 ,1 ,1 ,1 ,1 ,-1 ,-1 ,-1 ,-1 ,-1;
	_characterTable->addIrrep(v);
	v.setZero();
	v << 1 ,1 ,-1 ,-1 ,1 ,-1 ,1 ,-1 ,-1 ,1;
	_characterTable->addIrrep(v);

	// double degenerate irreps
	v.setZero();
	v << 2 ,-1 ,0 ,0 ,2 ,2 ,0 ,-1 ,2 ,0;
	_characterTable->addIrrep(v);
	v.setZero();
	v << 2 ,-1 ,0 ,0 ,2 ,-2 ,0 ,1 ,-2 ,0;
	_characterTable->addIrrep(v);

	// triple degenerate irreps
	v.setZero();
	v << 3 ,0 ,-1 ,1 ,-1 ,3 ,1 ,0 ,-1 ,-1 ;
	_characterTable->addIrrep(v);
	v.setZero();
	v << 3 ,0 ,1 ,-1 ,-1 ,3 ,-1 ,0 ,-1 ,1 ;
	_characterTable->addIrrep(v);
	v.setZero();
	v << 3 ,0 ,-1 ,1 ,-1 ,-3 ,-1 ,0 ,1 ,1;
	_characterTable->addIrrep(v);
	v.setZero();
	v << 3 ,0 ,1 ,-1 ,-1 ,-3 ,1 ,0 ,1 ,-1;
	_characterTable->addIrrep(v);
}

void PointGroup::generateIrrepsIh()
{
if (debugLevel > 3) {
	cout << "PointGroup::" << __FUNCTION__ << endl;
}
	// todo
	const double theta = M_PI / 5;
	VectorXcd v = VectorXcd::Zero(10);

	// non-degenerate irreps
	v << 1,1,1,1,1,1,1,1,1,1;
	_characterTable->addIrrep(v);
	v.setZero();
	v << 1,-1,1,1,1,1,-1,-1,-1,-1;
	_characterTable->addIrrep(v);


	// triple degenerate irreps
	v.setZero();
	v << 3,3,2*cos(theta),2*cos(3*theta),0,-1 ,2*cos(3*theta),2*cos(theta),0,-1 ;
	_characterTable->addIrrep(v);

	v.setZero();
	v << 3,3,2*cos(3*theta),2*cos(theta),0,-1 ,2*cos(theta),2*cos(3*theta),0,-1;
	_characterTable->addIrrep(v);

	v.setZero();
	v << 3,-3,2*cos(theta),2*cos(3*theta),0,-1 ,-2*cos(3*theta),-2*cos(theta),0,1;
	_characterTable->addIrrep(v);

	v.setZero();
	v << 3,-3,2*cos(3*theta),2*cos(theta),0,-1 ,-2*cos(theta),-2*cos(3*theta),0,1;
	_characterTable->addIrrep(v);


	// fourfold degenerate irreps
	v.setZero();
	v << 4,4,-1,-1,1,0,-1,-1,1,0;
	_characterTable->addIrrep(v);

	v.setZero();
	v << 4,-4,-1,-1,1,0,1,1,-1,0;
	_characterTable->addIrrep(v);


	// fivefold degenerate irreps
	v.setZero();
	v << 5,5,0,0,-1,1,0,0,-1,1;
	_characterTable->addIrrep(v);

	v.setZero();
	v << 5,-5,0,0,-1,1,0,0,1,-1;
	_characterTable->addIrrep(v);
}
