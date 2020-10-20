# include <cstdlib>
# include <string>

using namespace std;

int main(){

	////////////////////////////////////////////////////////////////////////////////////////
	//
	//  Get data from SCAN DFT
	//
	////////////////////////////////////////////////////////////////////////////////////////

	string make_dat= string(" > DAT");
	const char* md = make_dat.c_str();
	system(md);
	string get_dat= string("grep \"# of k-points:\" PROCAR | awk '{print $4}' >> DAT;") +
	string("grep \"# of k-points:\" PROCAR | awk '{print $8}' >> DAT;") +
	string("grep \"# of k-points:\" PROCAR | awk '{print $12}' >> DAT;") +
	string("grep \"weight\" PROCAR | awk '{print $4,$5,$6}' >> DAT;") + 
	string("grep \"occ.\" PROCAR | awk '{print $5}' >> DAT;") +
	string("grep \"tot\" PROCAR | grep -v ion | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10}' >> DAT;");
	const char* gd = get_dat.c_str();
	system(gd);

	
	return 0;
}
