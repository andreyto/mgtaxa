//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the MGTAXA package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>

typedef int int32;

using namespace std;


int main(int argc, char** argv) {
	if( argc < 3 || argv[1] != std::string("--out-path") ) {
		std::cerr << "Usage: " << argv[0] << " --out-path path\n";
		return 1;
	}
	std::string outPath = argv[2];
	std::string outPathVal = outPath + ".val";
	std::string outPathHdr = outPath + ".hdr";
	const int fileBufSize = 1000000;
	char outBuff[fileBufSize];
	char inpBuff[fileBufSize];
	FILE *outval = fopen(outPathVal.c_str(),"w");
	setbuffer(outval,outBuff,fileBufSize);
	setbuffer(stdin,inpBuff,fileBufSize);
	int nRows = 0, nCols = 0, nElem = 0;
	int rowElem = 2080 + 5;
	for(float x = 0; scanf("%f",&x) == 1; nCols++, nElem++) {
		fwrite(&x,sizeof(x),1,outval);
		nRows = nElem / rowElem;
		if( nElem % (rowElem * 10000) == 0 ) {
			std::cout << nRows << " records, " << nElem << " elements\n";
		}
		if( nRows >= 50000 ) {
			break;
		}
	}

	fclose(outval);
	
	std::ofstream outhdr(outPathHdr.c_str());
	outhdr << 
	"nRows=" << nRows << '\n' <<
	"nCols=" << nCols << '\n';
	outhdr.close();

	return 0;
}

