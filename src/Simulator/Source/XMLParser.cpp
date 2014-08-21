#include "XMLParser.h"
#include "cinder/app/AppNative.h"
#include "Rod.h"

using namespace ci;
using namespace std;

// Constructs a new XMLParser with an empty path.
XMLParser::XMLParser() {
	path = "";
}

// Constructs a new XMLParser that parses the file given by string p.
XMLParser::XMLParser(string p) {
	path = p;
}

// Parses the XML file given in the path.
vector<Rod> XMLParser::parse() {
	app::console() << "Printing to console." << endl;
	XmlTree doc(loadFile(path));
	XmlTree rodList = doc.getChild("rodList");
	vector<Rod> rodVector;
	for (XmlTree::Iter child = rodList.begin(); child != rodList.end(); ++child) {
		rodVector.push_back(parseRod(*child));
	}
	return rodVector;
}

// Parses the rod given by the XmlTree r.
Rod XMLParser::parseRod(XmlTree & r) {
	//Rod rod;
	size_t numCtrlPts = stoi(r.getChild("numCtrlPts").getValue());
	size_t numCoordinates = 3 * numCtrlPts;
	vector<Vec3f> positions;
	VecXe rodPos(numCoordinates);
	XmlTree posList = r.getChild("pos");
	size_t k = 0;
	for (XmlTree::Iter vec = posList.begin(); vec != posList.end(); ++vec) {
		Vec3f coords = parseVec(vec->getValue());
		positions.push_back(parseVec(vec->getValue()));
		rodPos(k) = coords[0];
		rodPos(k+1) = coords[1];
		rodPos(k+2) = coords[2];
		k += 3;
	}
	//rod.pos = positions;

	//vector<Vec3f> velocities;
	k = 0;
	VecXe rodVel(numCoordinates);
	if (r.hasChild("vel")) {
		XmlTree velList = r.getChild("vel");
		for (XmlTree::Iter vec = velList.begin(); vec != velList.end(); ++vec) {
			Vec3f coords = parseVec(vec->getValue());
			//velocities.push_back(parseVec(vec->getValue()));
			rodVel(k) = coords[0];
			rodVel(k + 1) = coords[1];
			rodVel(k + 2) = coords[2];
			k += 3;
		}
	} else {
		for (int k = 0; k != numCoordinates; k++) {
			//velocities.push_back(Vec3f::zero());
			rodVel(k) = 0.0f;
		}
	}
	//rod.vel = velocities;

	XmlTree ref = r.getChild("ref");
	//rod.ref = parseVec(ref.getChild("vec3f").getValue());
	Vec3f refVec = parseVec(ref.getChild("vec3f").getValue());
	Vec3e u0 = Vec3e(refVec[0], refVec[1], refVec[2]);
	
	Rod rod = Rod(rodPos, u0);
	//rod.console();
	return rod;
}

// Returns a vector constructed from the coordinates of s.
Vec3f XMLParser::parseVec(string s) {
	string delimiter = ", ";
	size_t pos = 0;
	string token;
	vector<float> coordinates;
	while ((pos = s.find(delimiter)) != string::npos) {
		token = s.substr(0, pos);
		float fl = stof(token);
		coordinates.push_back(fl);
		s.erase(0, pos + delimiter.length());
	}
	coordinates.push_back(stof(s));
	Vec3f vec = Vec3f(coordinates[0], coordinates[1], coordinates[2]);
	//Vec3f vec = Vec3f(0.0f, 0.0f, 0.0f);
	return vec;
}

// Outputs a string in XML format.
void XMLParser::write() {
	XmlTree rodList("rodList", "");
	
	XmlTree rod("rod", "");
	rod.push_back(XmlTree("numCtrlPts", "3"));
	XmlTree pos("pos", "");
	pos.push_back(XmlTree("vec3f", "0.0f, 0.0f, 0.0f"));
	pos.push_back(XmlTree("vec3f", "1.0f, 2.0f, 3.0f"));
	pos.push_back(XmlTree("vec3f", "5.0f, 5.0f, 5.0f"));
	rod.push_back(pos);
	XmlTree vel("vel", "");
	vel.push_back(XmlTree("vec3f", "0.0f, 0.0f, 0.0f"));
	vel.push_back(XmlTree("vec3f", "1.0f, 1.0f, 1.0f"));
	vel.push_back(XmlTree("vec3f", "-1.0f, -1.0f, -1.0f"));
	rod.push_back(vel);
	XmlTree cross("cross", "");
	cross.push_back(XmlTree("cricle", "5.0f"));
	rod.push_back(cross);
	XmlTree ref("ref", "");
	ref.push_back(XmlTree("vec3f", "1.0f, 1.0f, 1.0f"));
	rod.push_back(ref);
	rodList.push_back(rod);

	rod = XmlTree("rod", "");
	rod.push_back(XmlTree("numCtrlPts", "2"));
	pos = XmlTree("pos", "");
	pos.push_back(XmlTree("vec3f", "-15.0f, 6.72f, 25.0f"));
	pos.push_back(XmlTree("vec3f", "1.3f, -2.9f, 3.72586f"));
	rod.push_back(pos);
	vel = XmlTree("vel", "");
	vel.push_back(XmlTree("vec3f", "0.0f, 0.0f, 0.0f"));
	vel.push_back(XmlTree("vec3f", "1.0f, 1.0f, 1.0f"));
	rod.push_back(vel);
	cross = XmlTree("cross", "");
	cross.push_back(XmlTree("square", "9.0f"));
	rod.push_back(cross);
	ref = XmlTree("ref", "");
	ref.push_back(XmlTree("vec3f", "1.0f, 1.0f, 1.0f"));
	rod.push_back(ref);
	rodList.push_back(rod);

	rodList.write(writeFile("C:/Users/Travis/Desktop/example.txt"));
}
