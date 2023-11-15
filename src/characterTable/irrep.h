/*
 * irrep.h
 *
 *  Created on: Jun 20, 2023
 *      Author: pemueller
 */

#ifndef CHARACTERTABLE_IRREP_H_
#define CHARACTERTABLE_IRREP_H_

#include <string>
#include <vector>
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

class Irrep{
public:
	class MullikenSymbol
	{
	public:
		MullikenSymbol(){
			_degeneracy = "";
			_gerade = "";
			_number = "";
			_primes = "";
		}
		~MullikenSymbol(){}

		string getSymbol(){return _degeneracy + _number + _gerade + _primes;}

		string getDegeneracy(){return _degeneracy;}
		string getNumber(){ return _number;}
		string getGerade() { return _gerade;}
		string getPrimes() { return _primes; }

		void setDegeneracy( string s ) {_degeneracy = s;}
		void setGerade( string s ) {_gerade = s;}
		void setNumber( string s ) {_number = s;}
		void setPrimes( string s ) { _primes = s; }


	private:
		string _degeneracy;
		string _gerade;
		string _number;
		string _primes;
	};

	Irrep(){
		_characters.clear();
		_symbol = MullikenSymbol();
	}
	Irrep(VectorXcd v)
	{
		Irrep();
		_characters.push_back(v);
	}
	~Irrep(){}

	void setDegeneracy( string s ) { _symbol.setDegeneracy(s);}
	void setGerade( string s ) {_symbol.setGerade(s);}
	void setNumber( string s ) {_symbol.setNumber(s);}
	void setPrimes( string s ) { _symbol.setPrimes(s); }

	string getDegeneracy(){ return _symbol.getDegeneracy();}
	string getNumber(){ return  _symbol.getNumber();}
	string getGerade() { return  _symbol.getGerade();}
	string getPrimes() { return  _symbol.getPrimes(); }

	void add(VectorXcd v){_characters.push_back(v);}
	vector<VectorXcd> getCharacters(){return _characters;}
	VectorXcd getCharacters(int i){return _characters[i];}
	string getSymbol(){ return _symbol.getSymbol(); }

	bool lt ( Irrep other )
	{
		if ( debugLevel > 0 )
			cout << "comparing " << _symbol.getSymbol() << " and " << other._symbol.getSymbol() << endl;

		int a = -1;
		int b = -1;
		if ( _symbol.getGerade() == "g" )
			++a;
		if (other._symbol.getGerade() == "g" )
			++b;
		if ( _symbol.getGerade() == "u" )
			a += 2;
		if ( other._symbol.getGerade() == "u" )
			b += 2;

		if ( a != b)
			return a < b;

		if ( _symbol.getPrimes() == "\'" )
			++a;
		if (other._symbol.getPrimes() == "\'" )
			++b;
		if ( _symbol.getPrimes() == "\'\'" )
			a += 2;
		if ( other._symbol.getPrimes() == "\'\'" )
			b += 2;

		if ( a != b)
			return a < b;

		if ( _symbol.getDegeneracy() == "A" )
			++a;
		if ( other._symbol.getDegeneracy() == "A" )
			++b;
		if ( _symbol.getDegeneracy() == "B" )
			a+=2;
		if ( other._symbol.getDegeneracy() == "B" )
			b+=2;
		if ( _symbol.getDegeneracy() == "E" )
			a+=3;
		if ( other._symbol.getDegeneracy() == "E" )
			b+=3;
		if ( _symbol.getDegeneracy() == "T" )
			a+=4;
		if ( other._symbol.getDegeneracy() == "T" )
			b+=4;
		if ( _symbol.getDegeneracy() == "G" )
			a+=5;
		if ( other._symbol.getDegeneracy() == "G" )
			b+=5;
		if ( _symbol.getDegeneracy() == "H" )
			a+=6;
		if ( other._symbol.getDegeneracy() == "H" )
			b+=6;

		if ( a != b)
			return a < b;

		a += stoi(_symbol.getNumber());
		b += stoi(other._symbol.getNumber());

		return a < b;
	}

private:
	vector<VectorXcd> _characters;
	MullikenSymbol _symbol;
};


#endif /* CHARACTERTABLE_IRREP_H_ */
