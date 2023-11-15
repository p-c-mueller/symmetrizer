/*
 * pgTests.h
 *
 *  Created on: Jun 19, 2023
 *      Author: pemueller
 */

#ifndef TESTING_PGTESTS_H_
#define TESTING_PGTESTS_H_

#include <iostream>
#include <string>
#include <regex>

using namespace std;

extern int debugLevel;
extern string suggestedPointGroup;
extern string suggestedPointGroupType;

extern int intFromString( string in );

bool isCharacterTableCorrect();

bool isPointGroupValid(string pg);

#endif /* TESTING_PGTESTS_H_ */
