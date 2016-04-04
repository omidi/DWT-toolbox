#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include "evomodel.h"
#include <stdlib.h>

using namespace std;


// Constructor
evomodel::evomodel(string tree) {

  cerr << "Initialize EvoModel ... ";

  // Initialitation
  numnodes = 0;
  numinnodes = 0;
  numexnodes = 0;
  root = new node;
  root -> q = 0;
  root -> father = NULL;
  root -> fgmodel = -2;
  root -> name = "root";
  nodes.push_back(root);
  innodes.push_back(root);
  numnodes++;
  numinnodes++;
  bases.push_back("0");
  bases.push_back("1");
  bases.push_back("2");
  bases.push_back("3");
  basesgap.push_back("0");
  basesgap.push_back("1");
  basesgap.push_back("2");
  basesgap.push_back("3");
  basesgap.push_back("4");
  mapbase.push_back("A");
  mapbase.push_back("C");
  mapbase.push_back("G");
  mapbase.push_back("T");
  mapbase.push_back("-");

  // Reading tree
  readtree(tree,root);

  cerr << "done." << endl;
}


// Destructor
evomodel::~evomodel()
{
  delete root;
  cerr << "EvoModel has just been freed." << endl;
}


// Reading treefile
void evomodel::insertnode(node* father) {
  node* newnode = new node;
  newnode -> father = father;
  newnode -> fgmodel = -2;
  nodes.push_back(newnode);
  numnodes++;
}
void evomodel::readtree(string tree, node* nodetmp) {

  string name;
  node* newnode;
  int n = 0;
  int length = tree.length();

  // Constructing tree
  while (n < length) {

    // Creating a new node
    if (tree[n] == '(') {
      insertnode(nodetmp);
      newnode = nodes.back();
      (nodetmp -> children).push_back(newnode);
      nodetmp = newnode;
    }

    // Jumping up from current node to its father
    else if (tree[n] == ')') {
      nodetmp = nodetmp -> father;
    }

    // Jumping up to the father and creating a new node
    else if (tree[n] == ',') {
      nodetmp = nodetmp -> father;
      insertnode(nodetmp);
      newnode = nodes.back();
      (nodetmp -> children).push_back(newnode);
      nodetmp = newnode;
    }

    // Setting distance and the name
    else if (tree[n] == ':') {

      // Checking if the node is internal or external
      if (tree[n-1] == ')') {
	nodetmp -> name = "ancester";
	innodes.push_back(nodetmp);
	numinnodes++;
	name = "";
      }
      else {
	nodetmp -> name = name;
	exnodes.push_back(nodetmp);
	numexnodes++;
	name = "";
      }

      // Getting the distance
      string dis = "";
      while (tree[n+1] != ',' && tree[n+1] != ')') {
	n++;
	dis.push_back(tree[n]);
      }
      nodetmp -> q = exp( -atof(dis.c_str()) );
    }
    else {
      name.push_back(tree[n]);
    }
    n++;
  }
}

// Creating all possible ancesters
vector <string> evomodel::allcombases(vector <node*> setnodes, vector <string> inbases, int numanc, vector <string> setbases) {
  int n = 0;
  int i = 0;
  int len = inbases[0].length();
  int size = inbases.size();
  int numbases = setbases.size();
  vector <string> newset;

  if (len != numanc) {

    if (setnodes[len] -> fgmodel == -1) {
      for (n = 0; n < size; n++) {
	newset.push_back(inbases[n]+"4");
      }
    }
    else {
      for (n = 0; n < size; n++) {
	for (i = 0; i < numbases; i++) {
	  newset.push_back(inbases[n]+setbases[i]);
	}
      }
    }
    inbases = allcombases(setnodes,newset,numanc,setbases);
  }
  return(inbases);
}

// Computing coefficients
vector <double> evomodel::getcoeff(vector <double> setq) {
  int n = 0;
  int k = 0;
  int numq = setq.size();
  vector <double> coeffold;
  vector <double> coeff;

  coeffold.push_back(1);
  for (n = 1; n <= numq; n++) {
    coeffold.push_back(0);
  }


  for (n = 0; n < numq; n++) {
    coeff.clear();
    for (k = 0; k <= n+1; k++) {

      if (k == 0) {
	coeff.push_back(setq[n]*coeffold[k]);
      }
      else {
	coeff.push_back(setq[n]*coeffold[k]+coeffold[k-1]);
      }
    }
    coeffold = coeff;
  }
  return(coeff);
}

// Computing back model given a set of ancs.
vector <double> evomodel::backmodel(vector <int> &numeq, vector <int> &numdif, vector <vector <double>*> &setq, vector <double> &bg) {

  node* father;
  vector <double> vec;
  vector <double> probs;
  double prob = 0;
  double probset = 1;
  double q = 0;
  int fgmodel = 0;
  int n = 0;
  int fbase;
  int base;
  double tmp = 0;;

  base = nodes[0] -> base;
  prob = bg[base];

  for (n = 1; n < numnodes; n++) {
    father = nodes[n] -> father;
    fgmodel = nodes[n] -> fgmodel;
    fbase = father -> base;
    base = nodes[n] -> base;
    q = nodes[n] -> q;

    if (fgmodel != -1) {

      if (fbase == base) {
	tmp = (q+((1-q)*(bg[base])));
	prob = prob * (q+((1-q)*(bg[base])));
	if (fgmodel == 1) {
	  setq[base] -> push_back(q/(1-q));
	  numeq[base]++;
	}
	else {
	  probset = probset * (q+((1-q)*(bg[base])));
	}
      }
      else {
	tmp = ((1-q)*(bg[base]));
	prob = prob * ((1-q)*(bg[base]));
	if (fgmodel == 1) {
	  numdif[base]++;
	}
	else {
	  probset = probset * ((1-q)*(bg[base]));
	}
      }
    }
  }
  probs.push_back(prob);
  probs.push_back(probset);
  return(probs);
}

// Computing fore model
double evomodel::foremodel(vector <int> &numeq, vector <int> &numdif, vector <vector <double>*> &setq, vector <double> &bg, vector <double> &pseudo) {

  int n,base;
  int numq = 0;
  int rootbase = (nodes[0] -> base);
  double prob = 0;
  double product = 1;
  int lambdasum = 0;
  vector <int> num(4,0);
  vector <int> lambda;
  vector <vector <double> > coeff;
  vector <double> coeffbase;

  // Computing coefficients
  for (n = 0; n < 4; n++) {
    numq = setq[n] -> size();
    coeffbase.clear();

    if (numq > 0) {coeffbase = getcoeff(*setq[n]);}
    else          {coeffbase.push_back(1);}

    coeff.push_back(coeffbase);
  }

  // Computing prob
  for (num[0] = 0; num[0] <= numeq[0]; num[0]++) {
    for (num[1] = 0; num[1] <= numeq[1]; num[1]++) {
      for (num[2] = 0; num[2] <= numeq[2]; num[2]++) {
	for (num[3] = 0; num[3] <= numeq[3]; num[3]++) {

	  // Initialitation
	  lambda.clear();
	  lambdasum = 0;
	  product = 1;

	  // Computing counts
	  for (base = 0; base < 4; base++) {
	    if (rootbase == base) {
	      lambda.push_back(num[base] + numdif[base] + 1);
	    }
	    else {
	      lambda.push_back(num[base] + numdif[base]);
	    }

	    lambdasum = lambdasum + lambda[base];
	    product = product * (coeff[base][num[base]] * (gammabase[base][lambda[base]]/gammabase[base][0]));
	  }

	  // Summing probs
	  product = product * (gammasum[0]/gammasum[lambdasum]);
	  prob = prob + product;
	}
      }
    }
  }
  return(prob);
}

// Felling the tree according to the conservation pattern
int evomodel::getconpat (node* nodetmp) {
  vector <node*> children = nodetmp -> children;
  int num = children.size();
  int fgmodel = 0;
  int fgchild = 0;
  int n = 0;

  for (n = 0; n < num; n++) {
    fgmodel = nodetmp -> fgmodel;
    fgchild = children[n] -> fgmodel;
    if (fgchild == -2)     {fgchild = getconpat(children[n]);}
    if (fgchild > fgmodel) {nodetmp -> fgmodel = fgchild;}
  }

  return(nodetmp -> fgmodel);
}

// Getting otherw score (MAIN)
double evomodel::get_owscore(string column, string conpat, vector <double> bg) {

  // Checking if the input is correct
  string conall;
  conall = column[0];
  double test0 = bg[0]+bg[1]+bg[2]+bg[3];
  //  if ((conall != "a" && column.length() != abs(numexnodes)) || conpat.length() != abs(numexnodes)) {cerr << "Discrepances in the number of species" << endl;exit(1);}
  //if (test0 != 1.0)                                                                                {cerr << "Bg frequencies do not sum up to 1" << endl;exit(1);}

  string tmp = "";
  for(int l=0;l<abs(numexnodes);++l) {

    tmp.append("1");
  }

  conpat = tmp;

  // Declaration
  double proback = 0;
  double profore = 0;
  double pseudosum = 0;
  double K = 1;
  int numinbases = 0;
  int numexbases = 0;
  int n,i,k,j;
  string base;
  string cpat;
  string outseq;
  vector <double> probs;
  vector <int> numeq;
  vector <int> numdif;
  vector <string> inbases;
  vector <string> exbases;
  vector <double> pseudo;
  vector <vector <double>*> setq;
  gammavec.clear();
  gammasum.clear();
  gammabase.clear();
  cout.precision(10);

  // Printing out the names of the orgaisms
  if (conall == "a") {
    for (n = 0; n < numexnodes; n++) {
      cout << ">" << exnodes[n] -> name << endl;
    }
  }

  // Setting prior parameters
  vector<double> *test[4];
  for (n = 0; n < 4; n++) {
    test[n] = new vector<double>;
    setq.push_back(test[n]);
    pseudo.push_back(bg[n]);
    pseudosum = pseudosum + pseudo[n];
  }

  // Creating vectors of gammas -> This should be in the constructor
  for (i = 0; i < 4; i++) {
    gammavec.clear();

    for (n = 0; n <= 2*numexnodes+2; n++) {
      gammasum.push_back(exp(lgamma(pseudosum + n)));
      gammavec.push_back(exp(lgamma(pseudo[i] + n)));
    }
    gammabase.push_back(gammavec);
  }

  // Creating all possible ancesters
  if (conall == "a") {
    exbases = allcombases(exnodes,basesgap,numexnodes,basesgap);
    numexbases = exbases.size();
    numinbases = inbases.size();
  }
  else {
    exbases.push_back(column);
    numexbases = exbases.size();
  }

  // Computing score for all columns
  for (j = 0; j < numexbases; j++) {

    // Initialitation
    proback = 0;
    profore = 0;
    K = 1;

    // Filling leaves
    for (n = 0; n < numexnodes; n++) {
      base = exbases[j].substr(n,1);
      exnodes[n] -> base = atoi(base.c_str());

      if (base != "4") {
	cpat = conpat.substr(n,1);
	exnodes[n] -> fgmodel = atoi(cpat.c_str());
      }
      else {
	exnodes[n] -> fgmodel = -1;
      }
    }
    for (n = 0; n < numinnodes; n++) {
      innodes[n] -> fgmodel = -2;
    }
    getconpat(root);

    // Running if the column has at least one base
    if (innodes[0] -> fgmodel != -1) {

      // Creating all possible ancesters
      inbases = allcombases(innodes,bases,numinnodes,bases);
      numinbases = inbases.size();


      // Computing prefactor
      for (n = 1; n < numnodes; n++) {
	if (nodes[n] -> fgmodel == 1) {
	  K = K * (1 - (nodes[n] -> q));
	}
      }

      // Summing over all ancesters
      for (i = 0; i < numinbases; i++) {

	// Initialitation
	numeq.clear();
	numdif.clear();
	for (k = 0; k < 4; k++) {
	  setq[k] -> clear();
	  numeq.push_back(0);
	  numdif.push_back(0);
	}

	// Filling ancesters
	for (n = 0; n < numinnodes; n++) {
	  base = inbases[i].substr(n,1);
	  innodes[n] -> base = atoi(base.c_str());
	}

	// Computing back model
	probs = backmodel(numeq,numdif,setq,bg);
	proback = proback + probs[0];

	// Computing fore model
	profore = profore + K*probs[1]*foremodel(numeq,numdif,setq,bg,pseudo);
      }
    }

    // If the column is all gaps
    else {
      proback = 1;
      profore = 1;
    }

    // Printing out scores
    if (conall == "a") {
      outseq.clear();
      for (n = 0; n < numexnodes; n++) {
	base = exbases[j].substr(n,1);
	outseq.append(mapbase[atoi(base.c_str())]);
      }
      cout << scientific << outseq << "\t" << proback << "\t" << profore << "\t" << profore/proback << endl;
    }
  }

  if (conall == "a") {return(0);}
  else               {return(profore);}
}

