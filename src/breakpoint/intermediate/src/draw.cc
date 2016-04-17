#include "draw.h"

int draw_header(ofstream & fout)
{
	fout<<"\\documentclass{llncs}\n";
	fout<<"\\usepackage{tikz}	\n";
	fout<<"\\usetikzlibrary{calc}\n";
	fout<<"\\usetikzlibrary{shapes.geometric}\n";
	fout<<"\\usetikzlibrary{fit}\n";
	fout<<"\\usetikzlibrary{decorations}\n";
	fout<<"\\usepgflibrary{decorations.shapes}\n";
	fout<<"\\usepgflibrary{decorations.pathmorphing}\n";
	fout<<"\\begin{document}\n";
	fout<<"{\\begin{tikzpicture}[myrectangle/.style={draw, rectangle, minimum size=1.1em, inner sep = 0.5mm}, >=stealth]\n";
	fout<<"\\def\\cola{red}\n";
	fout<<"\\def\\colb{blue}\n";
	fout<<"\\def\\colc{green}\n";
	fout<<"\\def\\cold{gray}\n";
	fout<<"\\def\\cole{purple}\n";
	fout<<"\\def\\colf{brown}\n";
	fout<<"\\def\\colx{black}\n";
	return 0;
}

int draw_footer(ofstream & fout)
{
	fout<<"\\end{tikzpicture}}\n";
	fout<<"\\end{document}\n";
	return 0;
}
