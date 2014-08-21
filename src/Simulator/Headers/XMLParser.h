#pragma once
#include <string>
#include "cinder/Xml.h"
#include "cinder/Vector.h"
#include "Rod.h"

using namespace ci;
using namespace std;

class XMLParser {
public:
	XMLParser();
	XMLParser(string p);
	vector<Rod> parse();
	void write();

	string path;

private:
	Rod parseRod(XmlTree & r);
	Vec3f parseVec(string s);
};
