#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include "evomodel.h"
#include <stdlib.h>

using namespace std;

int main(int argc, char* argv[4]) {

  // Arguments
  //  if(argc != 8){cerr << "Usage: <treefile> <bgA> <bgC> <bgG> <bgT> <column> <conpat>" << endl;exit(1);}
  if(argc != 6){cerr << "Usage: <treefile> <bgA> <bgC> <bgG> <bgT>" << endl; exit(1);}

  vector <double> bg;
  char* treefile = argv[1];
  bg.push_back(atof(argv[2]));
  bg.push_back(atof(argv[3]));
  bg.push_back(atof(argv[4]));
  bg.push_back(atof(argv[5]));
  string column = string("a");
  string conpat = string(argv[7]);
  string tree;

  // Reading tree
  ifstream in(treefile);
  assert(in);
  while (!in.eof()){
    in >> tree;
  }

  // creating evomolde class
  evomodel myevo(tree);

  // getting the score
  double score = myevo.get_owscore(column,conpat,bg);
  cerr << score << endl;
}
