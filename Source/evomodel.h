#ifndef _evomodel_h_included_
#define _evomodel_h_included_
#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <stdlib.h>
using namespace std;


// Class evomodel
// --------------------------------------------
class evomodel {

 public:
  ~evomodel();
  evomodel(string);
  double get_owscore(string, string, vector <double>);

 private:

  // Structures and global variables
  struct node {
    int base;
    int fgmodel;
    double q;
    string name;
    node* father;
    vector <node*> children;
  };

  string tree;
  node* root;
  vector <node*> nodes;
  vector <node*> innodes;
  vector <node*> exnodes;
  vector <string> bases;
  vector <string> basesgap;
  vector <string> mapbase;
  vector <double> gammavec;
  vector <double> gammasum;
  vector <vector <double> > gammabase;
  int numnodes ;
  int numinnodes;
  int numexnodes;

  // Functions
  void insertnode(node*);
  void readtree(string, node*);
  vector <string> allcombases(vector <node*>, vector <string>, int, vector <string>);
  vector <double> getcoeff(vector <double>);
  vector <double> backmodel(vector <int>&, vector <int>&, vector <vector <double>*>&, vector <double>&);
  double foremodel(vector <int>&, vector <int>&, vector <vector <double>*>&, vector <double>&, vector <double>&);
  int getconpat (node*);
};

#endif
