# include <cstdlib>
# include <string>

using namespace std;

int main(){

	////////////////////////////////////////////////////////////////////////////////////////
	//
	//  Get data from SCAN DFT
	//
	////////////////////////////////////////////////////////////////////////////////////////

	string make_dat= string(" > SCAN_DAT");
	const char* md = make_dat.c_str();
	system(md);
	string get_dat= string("grep \"# of k-points:\" PROCAR | awk '{print $4}' >> SCAN_DAT;") +
	string("grep \"# of k-points:\" PROCAR | awk '{print $8}' >> SCAN_DAT;") +
	string("grep \"# of k-points:\" PROCAR | awk '{print $12}' >> SCAN_DAT;") +
	string("grep \"weight\" PROCAR | awk '{print $4,$5,$6}' >> SCAN_DAT;") + 
	string("grep \"occ.\" PROCAR | awk '{print $5}' >> SCAN_DAT;") +
        string("grep -A 4 \"ion\" PROCAR | grep -v \"tot\" | grep -v \"ion\" | awk 'NF==10{print}{}' >> SCAN_DAT;");
	const char* gd = get_dat.c_str();
	system(gd);

	
	return 0;
}
