// usage probename file, FASTA, output to COut 
#include <iostream>
#include <fstream>
#include <limits>  
#include <iterator>  
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>

using namespace std;

int main(int argc, char **argv)
{
    string transcript;    // string to store transcript proporties
    ifstream FASTAin;
    ifstream probeIn(argv[1], ios::binary);

    //reading probeIn into vector ProbeNames 
    if (!probeIn.is_open())
    {
        cout << "probein not open \n";
        return 0;
    }

    vector<string> ProbeNames;
    string probeName;
    while (probeIn.good())
    {
        getline(probeIn, probeName);
        ProbeNames.push_back(probeName.c_str());
    }

    // sort the vector with probenames for quicker finding;
    sort (ProbeNames.begin(), ProbeNames.end());

    //copy(ProbeNames.begin(), ProbeNames.end(), ostream_iterator<string>(cout, "\n"));
    // read fasta file
    FASTAin.open(argv[2]);
    if (!FASTAin.is_open())
    {
        cout << "probein not open \n";
        return 0;
    }

    int count = 0;
    int found = 0;
    string line;
    while (FASTAin.good())
    {
        getline(FASTAin, line);
  
            count++;
            if(line[0] == '>' && find_if(ProbeNames.begin(), ProbeNames.end(), [&](const string &str){
            	return (!str.empty() && line.find(string("HJ-ProbeName \"")+str) != string::npos);
            }) != ProbeNames.end())
            {
        	    found++; 
                // ignore fasta sequence & putback > if the inputstream is still good
                FASTAin.ignore(std::numeric_limits<std::streamsize>::max(), '>');
                if(FASTAin.good()) FASTAin.putback('>');
            }
            else
				 cout << line << '\n';
    }
    //cout << "Found " << found << "/" << count << " number of transcripts\n";

    probeIn.close();
    FASTAin.close();
    return 0;
}
