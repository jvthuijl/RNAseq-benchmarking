#include <iostream>
#include <fstream>  
#include <string>
#include <vector>
#include <string.h>
#include <stdlib.h>

using namespace std;
char bases[] = "ACGT";    // Bases
int mutation_rate = 2;    // 2 of 1000

std::string remove_extension(const std::string& filename) {
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot); 
}

int main(int argc, char **argv)
{
    srand(42);
    string transcript;    // string to store transcript proporties
    ofstream FASTAout (remove_extension(argv[1]) + "_mutated.fasta"); // output file
    ifstream FASTAin;
    FASTAin.open(argv[1]);
    if (FASTAin.is_open())
    {
        int pos = 0;          
        char c;       
        while (FASTAin.good()) {
            FASTAin.get(c);
            switch (c)
            {
                case '>':
                    FASTAin.unget();
                    getline (FASTAin, transcript);
                    if (pos > 0) FASTAout << '\n';
                    FASTAout << transcript << '\n';
                    pos = 0;
                    break;
                case '\n': // skip char if newline
                    break;
                default:
                    if (rand() % 1000 < mutation_rate)
                    {
                        char * basePointer = strchr(bases, c);
                        int basePos = basePointer - bases + 1;
                        cout << c; // old nucleotide
                        c = bases[(basePos + (rand() % 2 )) % 4]; // new 
                        cout << '\t' << c << '\t' <<  transcript <<  '\t' << pos << '\n';
                    }   

                    FASTAout << c;

                    if ((pos + 1)  % 60 == 0 && pos)
                        FASTAout << '\n';
                    pos++;
                    break;
            }     
        }

    }
    FASTAin.close();
    return 0;
}
