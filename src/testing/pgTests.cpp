/*
 * pgTests.cpp
 *
 *  Created on: Jun 19, 2023
 *      Author: pemueller
 */

#include <boost/algorithm/string.hpp>

#include "pgTests.h"
#include "../common/constants.h"

int debugLevel;

string suggestedPointGroup;
string suggestedPointGroupType;


bool isCharacterTableCorrect()
{
	return true;
}

bool isPointGroupValid(string pg)
{
	 boost::algorithm::to_lower(pg);

	 // Check for polyhedral point groups first
	if ( pg == "td" || pg == "oh" || pg == "ih" )
	{
		suggestedPointGroupType = pg;
		return true;
	}

	// C1, Ci, Cs
	if ( pg == "c1" || pg == "ci" || pg == "cs" )
	{
		suggestedPointGroupType = pg;
		return true;
	}

	// Cn, Sn, Dn
	if ( (pg[0] == 'c' || pg[0] == 's' || pg[0] == 'd' ) && regex_match(pg.substr(1, pg.size()), regex("[0-9]+")) )
	{
		suggestedPointGroupType = pg.substr(0,1) + 'n';
		return true;
	}

	// Cnh, Dnh
	if ( (pg[0] == 'c' || pg[0] == 'd' ) && regex_match(pg.substr(1, pg.size()-2), regex("[0-9]+")) && pg.back() == 'h')
	{
		suggestedPointGroupType = pg.substr(0,1) + "nh";
		return true;
	}

	// Cnv
	if ( (pg[0] == 'c' ) && regex_match(pg.substr(1, pg.size()-2), regex("[0-9]+")) && pg.back() == 'v' )
	{
		suggestedPointGroupType = "cnv";
		return true;
	}

	// Dnd
	if ( ( pg[0] == 'd' ) && regex_match(pg.substr(1, pg.size()-2), regex("[0-9]+")) && pg.back() == 'd' )
	{
		suggestedPointGroupType = "dnd";
		return true;
	}

	if ( pg == "cuv" || pg == "duh" )
	{
		suggestedPointGroupType = pg.substr(0,1) + to_string(infinityRotationSteps) + pg.back();
		return true;
	}

	return false;
}

int intFromString(string in)
{
	string result;
	for ( size_t i = 0; i < in.size(); ++i )
		if (regex_match(in.substr(i,i), regex("[0-9]+")))
			result += in.substr(i,i);

	if ( result != "" )
		return stoi(result);
	else
		return 0;
}
