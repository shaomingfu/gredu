#ifndef __EXTREMITY_H__
#define __EXTREMITY_H__

#include "gene.h"

#include <string> 

using namespace std;

class extremity
{
public:
	gene * g;			// pointer to the gene
	bool b;				// true for head, false for tail

public:
	extremity();
	extremity(gene * _g, bool _b);
	extremity(const extremity & e);
	
public:
	bool operator==(const extremity & e) const;
	bool operator!=(const extremity & e) const;
	bool forward() const;
	bool complement(const extremity & e) const;
	int label() const;
	string label_str() const;
	string label_tex() const;
	int weak_compare(const extremity & e) const;
};

#endif
