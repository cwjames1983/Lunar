/** 
 * Function that takes an input track and splits it into as many equal-length subtracks as you want
 * In this version of the code, I am modifying the 4m long track called "trivial_tracks4.dat"
 */

#include <iostream>
#include <fstream>
#include <sstream>

int main(void) {
    int max_sub_tracks = 8; //Set this to the number of equal-length sub tracks you want
    int c=3e8;
    std::ifstream in("trivial_tracks4.dat");

    char* base = "subtrack_";
    char filename[256];

    int temp, i, j;
    double x[9];

    in >> temp;

    for(i=0;i<9;i++) {
        in>>x[i];
    }
    for(i=0;i<max_sub_tracks;i++) {
        sprintf(filename, "SubTracks/%s%d.in", base, i);
        std::ofstream out(filename);
        out<<"1\n";
        for(j=0;j<9;j++) {
            if(j==0) {
                out<<(x[7]-x[3])/max_sub_tracks * i /c<<" ";
            }
            else if(j==3) {
                out<<(x[7]-x[3])/max_sub_tracks * i<<" ";
            }
            else if(j==4||j==8) {
                out<<x[j]<<" ";
            }
            else if(j==7) {
                out<<(x[7]-x[3])/max_sub_tracks * (i+1)<<" ";
            }
            else {
                out<<"0 ";
            }
        }
    }

    return 0;
}
