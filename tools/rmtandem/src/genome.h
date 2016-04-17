#ifndef __GENOME_H__
#define __GENOME_H__

#include "genome_base.h"
#include "chrm.h"
#include "operation.h"

#include <vector>
#include <map>
#include <set>
#include <string>

using namespace std;

class genome : public genome_base
{
public:
	int load(const string & file);
	int write(const string & file, const string & prefix);

	int remove_tandem();
	bool remove_tandem_no_sign(chrm * ch, int p);
	bool remove_tandem_linear(chrm * ch, int p);
	bool remove_tandem_linear_triple(chrm * ch);
};

#endif
