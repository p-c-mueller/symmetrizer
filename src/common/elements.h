/*
 * elements.h
 *
 *  Created on: Feb 6, 2023
 *      Author: pemueller
 */

#ifndef ELEMENTS_H_
#define ELEMENTS_H_

#include <string>
#include <boost/assign.hpp>
#include <boost/bimap.hpp>

using namespace std;

const boost::bimap<int, string> elementSymbols = boost::assign::list_of<boost::bimap<int,string>::relation>
	(0,"x")
	(1,"h") (2,"he") (3,"li") (4,"be") (5,"b") (6,"c") (7,"n") (8,"o")
	(9,"f") (10,"ne") (11,"na") (12,"mg") (13,"al") (14,"si") (15,"p") (16,"s")
	(17,"cl") (18,"ar") (19,"k") (20,"ca") (21,"sc") (22,"ti") (23,"v") (24,"cr")
	(25,"mn") (26,"fe") (27,"co") (28,"ni") (29,"cu") (30,"zn") (31,"ga") (32,"ge")
	(33,"as") (34,"se") (35,"br") (36,"kr") (37,"rb") (38,"sr") (39,"y") (40,"zr")
	(41,"nb") (42,"mo") (43,"tc") (44,"ru") (45,"rh") (46,"pd") (47,"ag") (48,"cd")
	(49,"in") (50,"sn") (51,"sb") (52,"te") (53,"i") (54,"xe") (55,"cs") (56,"ba")
	(57,"la") (58,"ce") (59,"pr") (60,"nd") (61,"pm") (62,"sm") (63,"eu") (64,"gd")
	(65,"tb") (66,"dy") (67,"ho") (68,"er") (69,"tm") (70,"yb") (71,"lu") (72,"hf")
	(73,"ta") (74,"w") (75,"re") (76,"os") (77,"ir") (78,"pt") (79,"au") (80,"hg")
	(81,"tl") (82,"pb") (83,"bi") (84,"po") (85,"at") (86,"rn") (87,"fr") (88,"ra")
	(89,"ac") (90,"th") (91,"pa") (92,"u") (93,"np") (94,"pu") (95,"am") (96,"cm")
	(97,"bk") (98,"cf") (99,"es") (100,"fm") (101,"md") (102,"no") (103,"lr")
	(104,"rf") (105,"db") (106,"sg") (107,"bh") (108,"hs") (109,"mt") (110,"ds")
	(111,"rg") (112,"cn") (113,"nh") (114,"fl") (115,"mc");

const map<string, string> electronConfigs = boost::assign::map_list_of
		("h","1s1") ("he","1s2") ("li","1s2 2s1") ("be","1s2 2s2")
		("b","1s2 2s2 2p1") ("c","1s2 2s2 2p2") ("n","1s2 2s2 2p3") ("o","1s2 2s2 2p4")
		("f","1s2 2s2 2p5") ("ne","1s2 2s2 2p6") ("na","[ne] 3s1") ("mg","[ne] 3s2")
		("al","[ne] 3s2 3p1") ("si","[ne] 3s2 3p2") ("p","[ne] 3s2 3p3") ("s","[ne] 3s2 3p4")
		("cl","[ne] 3s2 3p5") ("ar","[ne] 3s2 3p6") ("k","[ar] 4s1") ("ca","[ar] 4s2")
		("sc","[ar] 3d1 4s2") ("ti","[ar] 3d2 4s2") ("v","[ar] 3d3 4s2") ("cr","[ar] 3d5 4s1")
		("mn","[ar] 3d5 4s2") ("fe","[ar] 3d6 4s2") ("co","[ar] 3d7 4s2") ("ni","[ar] 3d8 4s2")
		("cu","[ar] 3d10 4s1") ("zn","[ar] 3d10 4s2") ("ga","[ar] 3d10 4s2 4p1") ("ge","[ar] 3d10 4s2 4p2")
		("as","[ar] 3d10 4s2 4p3") ("se","[ar] 3d10 4s2 4p4") ("br","[ar] 3d10 4s2 4p5") ("kr","[ar] 3d10 4s2 4p6")
		("rb","[kr] 5s1") ("sr","[kr] 5s2") ("y","[kr] 4d1 5s2") ("zr","[kr] 4d2 5s2")
		("nb","[kr] 4d4 5s1") ("mo","[kr] 4d5 5s1") ("tc","[kr] 4d5 5s2") ("ru","[kr] 4d7 5s1")
		("rh","[kr] 4d8 5s1") ("pd","[kr] 4d10") ("ag","[kr] 4d10 5s1") ("cd","[kr] 4d10 5s2")
		("in","[kr] 4d10 5s2 5p1") ("sn","[kr] 4d10 5s2 5p2") ("sb","[kr] 4d10 5s2 5p3") ("te","[kr] 4d10 5s2 5p4")
		("i","[kr] 4d10 5s2 5p5") ("xe","[kr] 4d10 5s2 5p6") ("cs","[xe] 6s1") ("ba","[xe] 6s2")
		("la","[xe] 5d1 6s2") ("ce","[xe] 4f1 5d1 6s2") ("pr","[xe] 4f3 6s2") ("nd","[xe] 4f4 6s2")
		("pm","[xe] 4f5 6s2") ("sm","[xe] 4f6 6s2") ("eu","[xe] 4f7 6s2") ("gd","[xe] 4f7 5d1 6s2")
		("tb","[xe] 4f9 6s2") ("dy","[xe] 4f10 6s2") ("ho","[xe] 4f11 6s2") ("er","[xe] 4f12 6s2")
		("tm","[xe] 4f13 6s2") ("yb","[xe] 4f14 6s2") ("lu","[xe] 4f14 5d1 6s2") ("hf","[xe] 4f14 5d2 6s2")
		("ta","[xe] 4f14 5d3 6s2") ("w","[xe] 4f14 5d4 6s2") ("re","[xe] 4f14 5d5 6s2") ("os","[xe] 4f14 5d6 6s2")
		("ir","[xe] 4f14 5d7 6s2") ("pt","[xe] 4f14 5d9 6s1") ("au","[xe] 4f14 5d10 6s1") ("hg","[xe] 4f14 5d10 6s2")
		("tl","[xe] 4f14 5d10 6s2 6p1") ("pb","[xe] 4f14 5d10 6s2 6p2") ("bi","[xe] 4f14 5d10 6s2 6p3") ("po","[xe] 4f14 5d10 6s2 6p4")
		("at","[xe] 4f14 5d10 6s2 6p5") ("rn","[xe] 4f14 5d10 6s2 6p6") ("fr","[rn] 7s1") ("ra","[rn] 7s2")
		("ac","[rn] 6d1 7s2") ("th","[rn] 6d2 7s2") ("pa","[rn] 5f2 6d1 7s2") ("u","[rn] 5f3 6d1 7s2")
		("np","[rn] 5f4 6d1 7s2") ("pu","[rn] 5f6 7s2") ("am","[rn] 5f7 7s2") ("cm","[rn] 5f7 6d1 7s2")
		("bk","[rn] 5f9 7s2") ("cf","[rn] 5f10 7s2") ("es","[rn] 5f11 7s2") ("fm","[rn] 5f12 7s2")
		("md","[rn] 5f13 7s2") ("no","[rn] 5f14 7s2") ("lr","[rn] 5f14 6d1 7s2") ("rf","[rn] 5f14 6d2 7s2")
		("db","[rn] 5f14 6d3 7s2") ("sg","[rn] 5f14 6d4 7s2") ("bh","[rn] 5f14 6d5 7s2") ("hs","[rn] 5f14 6d6 7s2")
		("mt","[rn] 5f14 6d7 7s2") ("ds","[rn] 5f14 6d8 7s2") ("rg","[rn] 5f14 6d9 7s2") ("cn","[rn] 5f14 6d10 7s2")
		("nh","[rn] 5f14 6d10 7s2 7p1") ("fl","[rn] 5f14 6d10 7s2 7p2") ("mc","[rn] 5f14 6d10 7s2 7p3");




#endif /* ELEMENTS_H_ */
