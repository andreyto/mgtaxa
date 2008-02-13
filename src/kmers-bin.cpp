#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>

typedef int int32;

template<typename T>
void
writeRaw(std::ostream& out,const T& x, int nitems=1) {
	out.write(reinterpret_cast<const char*>(&x),sizeof(x)*nitems);
}

int main(int argc, char** argv) {
	if( argc < 3 || argv[1] != std::string("--out-path") ) {
		std::cerr << "Usage: " << argv[0] << " --out-path path\n";
		return 1;
	}
	std::string outPath = argv[2];
	std::string outPathVal = outPath + ".val";
	std::string outPathHdr = outPath + ".hdr";
	std::ofstream outval(outPathVal.c_str()); 
	int nRows = 0, nCols = 0, nElem = 0;
	char lineBuff[10000];
	lineBuff[10000-1]=0;

	for(std::string s(lineBuff); std::getline(std::cin,s); nRows++) {

		std::istringstream istr(s);
		int32 taxid = 0;
		istr >> taxid;
		writeRaw(outval,taxid);
		for(int i = 1, dummy = 0; i < 5; i++) {
			istr >> dummy;
		}
		nCols = 0;
		for(float x = 0; istr >> x; nCols++, nElem++) {
			//std::cout << x << ' ';
			writeRaw(outval,x);
		}
		if( nRows > 0 && nRows % 10000 == 0 ) {
			std::cout << nRows << " records, " << nElem << " elements\n";
		}
		if( nRows >= 50000 ) {
			break;
		}
	}

	outval.close();

	std::ofstream outhdr(outPathHdr.c_str());
	outhdr << 
	"nRows=" << nRows << '\n' <<
	"nCols=" << nCols << '\n';
	outhdr.close();

	return 0;
}

